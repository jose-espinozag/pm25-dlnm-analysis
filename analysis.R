###############################################################
# Core analytical workflow for distributed lag non-linear models (DLNM)
# Functional form selection → Lag-response selection → Primary model
# → Effect modification → Stratified analyses
###############################################################
library(dlnm)
library(splines)
library(survival)
library(data.table)
library(lmtest)

setDT(bd_integral)
###############################################################
# 0. Global settings
###############################################################
max_lag_final <- 60
lag_knots_final <- c(15, 30)
###############################################################
# 1. Log-transformed PM2.5 variable and lag structure
###############################################################
bd_integral[, PM_RF_log := log(PM_RF)]

for(i in 1:max_lag_final){
  
  bd_integral[, paste0("PM_RF_log_Lag", i) :=
                log(get(paste0("PM_RF_Lag", i)))]
}
###############################################################
# 2. Cross-basis constructor
###############################################################
create_cb <- function(data,
                      varname,
                      df_var,
                      max_lag,
                      lag_knots){
  
  lag_matrix <- data[, c(varname,
                         paste0(varname, "_Lag", 1:max_lag)),
                     with = FALSE]
  
  if(df_var == 1){
    
    crossbasis(
      lag_matrix,
      lag = max_lag,
      argvar = list(fun = "lin"),
      arglag = list(fun = "ns",
                    knots = lag_knots)
    )
    
  } else {
    
    crossbasis(
      lag_matrix,
      lag = max_lag,
      argvar = list(fun = "ns",
                    df = df_var),
      arglag = list(fun = "ns",
                    knots = lag_knots)
    )
  }
}
###############################################################
# 3. Lag–response sensitivity analyses
#    (Lag 30 vs 45 vs 60)
###############################################################
run_lag_model <- function(max_lag,
                          lag_knots,
                          linear_lag = FALSE){
  
  if(linear_lag){
    
    cbPM <- crossbasis(
      bd_integral[, c("PM_RF",
                      paste0("PM_RF_Lag", 1:max_lag)),
                  with = FALSE],
      lag = max_lag,
      argvar = list(fun = "ns", df = 3),
      arglag = list(fun = "lin")
    )
    
    cbt2m <- crossbasis(
      bd_integral[, c("t2m",
                      paste0("t2m_Lag", 1:max_lag)),
                  with = FALSE],
      lag = max_lag,
      argvar = list(fun = "ns", df = 3),
      arglag = list(fun = "lin")
    )
    
    cbHR <- crossbasis(
      bd_integral[, c("HR",
                      paste0("HR_Lag", 1:max_lag)),
                  with = FALSE],
      lag = max_lag,
      argvar = list(fun = "ns", df = 3),
      arglag = list(fun = "lin")
    )
    
  } else {
    
    cbPM <- create_cb(bd_integral, "PM_RF", 3, max_lag, lag_knots)
    cbt2m <- create_cb(bd_integral, "t2m", 3, max_lag, lag_knots)
    cbHR <- create_cb(bd_integral, "HR", 3, max_lag, lag_knots)
  }
  
  clogit(
    CASO ~ cbPM + cbt2m + cbHR +
      holiday + strata(N),
    data = bd_integral
  )
}

model_lag30 <- run_lag_model(max_lag = 30, lag_knots = c(10, 20))
model_lag45 <- run_lag_model(max_lag = 45, lag_knots = c(15, 30))
model_lag60 <- run_lag_model(max_lag = 60, lag_knots = c(15, 30))

model_lag60_linear <- run_lag_model(max_lag = 60, lag_knots = c(15, 30), linear_lag = TRUE)
###############################################################
# 4. Final lag specification
###############################################################
max_lag <- 60
lag_knots <- c(15, 30)
###############################################################
# 5. Exposure–response functional form assessment
#    Linear vs log-linear vs spline
###############################################################
cbPM_lin <- create_cb(bd_integral, "PM_RF", 1, max_lag, lag_knots)
cbPM_log <- create_cb(bd_integral, "PM_RF_log", 1, max_lag, lag_knots)
cbPM_ns3 <- create_cb(bd_integral, "PM_RF", 3, max_lag, lag_knots)
cbt2m <- create_cb(bd_integral, "t2m", 3, max_lag, lag_knots)
cbHR <- create_cb(bd_integral, "HR", 3, max_lag, lag_knots)

model_lin <- clogit(CASO ~ cbPM_lin + cbt2m + cbHR + holiday + strata(N), data = bd_integral)

model_log <- clogit(CASO ~ cbPM_log + cbt2m + cbHR + holiday + strata(N), data = bd_integral)

model_ns3 <- clogit(CASO ~ cbPM_ns3 + cbt2m + cbHR + holiday + strata(N), data = bd_integral)
###############################################################
# 6. Formal assessment of non-linearity
###############################################################
lr_linearity <- lrtest(
  model_lin,
  model_ns3
)

p_linearity <- lr_linearity$`Pr(>Chisq)`[2]
###############################################################
# 7. Spline flexibility assessment
###############################################################
df_values <- 2:5

spline_models <- list()
spline_results <- data.table()

for(i in seq_along(df_values)){
  
  df_tmp <- df_values[i]
  
  cbPM_tmp <- create_cb(
    bd_integral,
    "PM_RF",
    df_tmp,
    max_lag,
    lag_knots
  )
  
  model_tmp <- clogit(
    CASO ~ cbPM_tmp + cbt2m + cbHR +
      holiday + strata(N),
    data = bd_integral
  )
  
  spline_models[[paste0("df", df_tmp)]] <- model_tmp
  
  spline_results <- rbind(
    spline_results,
    data.table(
      df_PM = df_tmp,
      k = length(coef(model_tmp)),
      Minus2LogLik = -2 * as.numeric(logLik(model_tmp)),
      AIC = AIC(model_tmp)
    )
  )
}
###############################################################
# 8. Sequential likelihood ratio tests
###############################################################
lrt_results <- data.table(
  Comparison = character(),
  Delta_k = numeric(),
  LRT = numeric(),
  p_value = numeric()
)

for(i in 2:length(df_values)){
  
  model_prev <- spline_models[[i - 1]]
  model_curr <- spline_models[[i]]
  
  lrt_tmp <- lrtest(model_prev, model_curr)
  
  lrt_results <- rbind(
    lrt_results,
    data.table(
      Comparison = paste0(
        "df ",
        df_values[i - 1],
        " vs df ",
        df_values[i]
      ),
      Delta_k = length(coef(model_curr)) -
        length(coef(model_prev)),
      LRT = lrt_tmp$Chisq[2],
      p_value = lrt_tmp$`Pr(>Chisq)`[2]
    )
  )
}
###############################################################
# 9. Final PM2.5 specification
###############################################################
best_df_PM <- spline_results[
  which.min(AIC),
  df_PM
]

cbPM_final <- create_cb(
  bd_integral,
  "PM_RF",
  best_df_PM,
  max_lag,
  lag_knots
)

argvar_PM <- attr(cbPM_final, "argvar")
arglag_PM <- attr(cbPM_final, "arglag")

argvar_t2m <- attr(cbt2m, "argvar")
arglag_t2m <- attr(cbt2m, "arglag")

argvar_HR <- attr(cbHR, "argvar")
arglag_HR <- attr(cbHR, "arglag")

final_model_main <- clogit(
  CASO ~ cbPM_final + cbt2m + cbHR +
    holiday + strata(N),
  data = bd_integral
)
###############################################################
# 10. Effect modification across residential proximity thresholds
###############################################################
proximity_thresholds <- c(100, 200, 300, 400, 500)

interaction_results <- data.table()
interaction_models <- list()

for(th in proximity_thresholds){
  
  proximity_var <- paste0("proximity_", th)
  
  model_no_interaction <- clogit(
    CASO ~ cbPM_final + cbt2m + cbHR +
      holiday + strata(N),
    data = bd_integral
  )
  
  model_interaction <- clogit(
    as.formula(
      paste0(
        "CASO ~ cbPM_final * ",
        proximity_var,
        " + cbt2m + cbHR + holiday + strata(N)"
      )
    ),
    data = bd_integral
  )
  
  lr_interaction <- lrtest(
    model_no_interaction,
    model_interaction
  )
  
  interaction_results <- rbind(
    interaction_results,
    data.table(
      Threshold_m = th,
      LRT = lr_interaction$Chisq[2],
      p_value = lr_interaction$`Pr(>Chisq)`[2]
    )
  )
  
  interaction_models[[proximity_var]] <- model_interaction
}
###############################################################
# 11. Stratified analyses
###############################################################
run_dlnm_strat <- function(data_subset){
  
  PM_lag <- data_subset[, c(
    "PM_RF",
    paste0("PM_RF_Lag", 1:max_lag)),
    with = FALSE]
  
  t2m_lag <- data_subset[, c(
    "t2m",
    paste0("t2m_Lag", 1:max_lag)),
    with = FALSE]
  
  HR_lag <- data_subset[, c(
    "HR",
    paste0("HR_Lag", 1:max_lag)),
    with = FALSE]
  
  cbPM <- crossbasis(PM_lag, lag = max_lag, argvar = argvar_PM, 
arglag = arglag_PM)
  
  cbt2m <- crossbasis(t2m_lag, lag = max_lag, argvar = argvar_t2m, arglag = arglag_t2m)
  
  cbHR <- crossbasis(HR_lag, lag = max_lag, argvar = argvar_HR, 
arglag = arglag_HR)
  
  clogit(
    CASO ~ cbPM + cbt2m + cbHR +
      holiday + strata(N),
    data = data_subset
  )
}

run_stratified <- function(var){
  
  split_data <- split(
    bd_integral,
    bd_integral[[var]]
  )
  
  lapply(
    split_data,
    function(d)
      run_dlnm_strat(as.data.table(d))
  )
}

models_sex <- run_stratified("SEXO")
models_age <- run_stratified("EDAD_cat")
models_season <- run_stratified("SEASON_med")
models_nse <- run_stratified("NSE_med")

models_dist_400 <- run_stratified("proximity_400")
###############################################################
# 12. Double-stratified analyses
###############################################################
run_double_strat <- function(var1, var2){
  
  bd_integral[, group :=
                interaction(
                  get(var1),
                  get(var2),
                  drop = TRUE
                )]
  
  split_data <- split(
    bd_integral,
    bd_integral$group
  )
  
  models <- lapply(
    split_data,
    function(d)
      run_dlnm_strat(as.data.table(d))
  )
  
  bd_integral[, group := NULL]
  
  models
}

models_sex_dist <- run_double_strat("SEXO", "proximity_400")
models_age_dist <- run_double_strat("EDAD_cat", "proximity_400")
models_season_dist <- run_double_strat("SEASON_med", "proximity_400")
models_nse_dist <- run_double_strat("NSE_med", "proximity_400")
###############################################################
# End of analytical workflow
###############################################################
