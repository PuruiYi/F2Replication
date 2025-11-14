# DATA TO CONSTRUCT RESPONSE FUNCTIONS IN FIGURE 2

# Load required packages
library(haven)
library(lfe)
library(marginaleffects)
library(broom)
library(dplyr)

# ---------- ESTIMATE GLOBAL RESPONSE WITH BASELINE REGRESSION SPECIFICATION
########################################################
#  CONSTRUCTION OF FIGURE 2 Panel A
########################################################

# Load data
data <- read_dta("data/input/GrowthClimateDataset.dta")
# Generate temp variable
data$temp <- data$UDel_temp_popweight
# Create squared precipitation term
data$UDel_precip_popweight_2 <- data$UDel_precip_popweight^2

# Create country-specific time trends if _yi_ and _y2_ variables don't exist
# Check if these variables exist in the dataset
yi_vars <- grep("^_yi_", names(data), value = TRUE)
y2_vars <- grep("^_y2_", names(data), value = TRUE)

if (length(yi_vars) == 0) {
  # Create time variables if they don't exist
  if (!"time" %in% names(data)) {
    data$time <- data$year - 1985
    data$time2 <- data$time^2
  }
}

# Generate second degree polynomial model in the form of: y = ax^2 + bx + c 
if (length(yi_vars) > 0 & length(y2_vars) > 0) {
  # If _yi_ and _y2_ variables exist, include them with backticks
  yi_vars_quoted <- paste0("`", yi_vars, "`")
  y2_vars_quoted <- paste0("`", y2_vars, "`")
  
  formula_str <- paste0("growthWDI ~ temp + I(temp^2) + UDel_precip_popweight + ",
                        "UDel_precip_popweight_2 + ",
                        paste(yi_vars_quoted, collapse = " + "), " + ",
                        paste(y2_vars_quoted, collapse = " + "),
                        " | iso_id + year | 0 | iso_id")
  model1 <- felm(as.formula(formula_str), data = data)
} else {
  # Otherwise create country-specific time trends on the fly
  model1 <- felm(growthWDI ~ temp + I(temp^2) + UDel_precip_popweight + 
                   UDel_precip_popweight_2 | 
                   iso_id + year + iso_id:time + iso_id:time2 | 0 | iso_id, 
                 data = data)
}

# Extract coefficients for temp and temp^2
coef_temp <- coef(model1)["temp"]
coef_temp2 <- coef(model1)["I(temp^2)"]

# Calculate optimal temperature
optimal_temp <- -coef_temp / (2 * coef_temp2)   # x-value of vertex: -b/2a 
print(paste("Optimal temperature:", optimal_temp))

# ---------- marginal analysis
# Generate margins/predictions at different temperature values
temp_seq <- seq(-5, 35, by = 1)

# Create prediction data - use mean values for other variables
pred_data <- data.frame(
  temp = temp_seq,
  UDel_precip_popweight = mean(data$UDel_precip_popweight, na.rm = TRUE)
)
pred_data$UDel_precip_popweight_2 <- pred_data$UDel_precip_popweight^2

# Get coefficients from model
all_coefs <- coef(model1)

# Identify the climate variable coefficients (assume no fixed effects in estimation)
coef_names <- c("temp", "I(temp^2)", "UDel_precip_popweight", "UDel_precip_popweight_2")
coef_names_available <- coef_names[coef_names %in% names(all_coefs)]

if (length(coef_names_available) == 0) {
  stop("Climate coefficients not found in model")
}

# Get the relevant coefficients
coefs <- all_coefs[coef_names_available]

# Create design matrix matching the available coefficients
X_pred <- model.matrix(~ temp + I(temp^2) + UDel_precip_popweight + 
                         UDel_precip_popweight_2, data = pred_data)

# Reorder to match coefficient order
common_names <- intersect(colnames(X_pred), names(coefs))
X_pred <- X_pred[, common_names, drop = FALSE]
coefs <- coefs[common_names]
cat(" coefs :", coefs, "\n")

# Calculate predictions
pred_data$estimate <- as.vector(X_pred %*% coefs)

# Calculate standard errors
vcov_mat <- vcov(model1)
vcov_subset <- vcov_mat[common_names, common_names, drop = FALSE]
se_vec <- sqrt(diag(X_pred %*% vcov_subset %*% t(X_pred)))
pred_data$stderr <- se_vec

# 90% confidence intervals (z = 1.645)
pred_data$min90 <- pred_data$estimate - 1.645 * pred_data$stderr
pred_data$max90 <- pred_data$estimate + 1.645 * pred_data$stderr

# ---------- Format output for Figure 2, Panel A
global_response <- pred_data %>%
  # Change the name of the header
  mutate(x = temp,
         dof = model1$df.residual,  # degrees of freedom
         parm = paste0("_at#", row_number())  # parameter to comply with STATA format
         ) %>%  
  select(estimate, stderr, dof, min90, max90, parm, x)

# Write out results
write.csv(global_response, "data/output/estimatedGlobalResponse.csv", 
          row.names = FALSE)

# Write out main dataset subset
main_dataset <- data %>%
  select(UDel_temp_popweight, Pop, TotGDP, growthWDI, 
         GDPpctile_WDIppp, continent, iso, countryname, year)

write.csv(main_dataset, "data/output/mainDataset.csv", row.names = FALSE)

# Write out coefficients
coefficients_df <- data.frame(
  temp = coef_temp,
  temp_squared = coef_temp2
)
write.csv(coefficients_df, "data/output/estimatedCoefficients.csv", 
          row.names = FALSE)


# -----------------  DATA FOR FIGURE 2, PANELS B, D, E  -----------------------

outcome_vars <- c("growthWDI", "AgrGDPgrowthCap", "NonAgrGDPgrowthCap")
heterogeneity_results <- list()

for (var in outcome_vars) {
  
  # Reload data
  data <- read_dta("data/input/GrowthClimateDataset.dta")
  
  # Create time variables
  data$time <- data$year - 1985
  data$time2 <- data$time^2
  
  # Generate temperature variable
  data$temp <- data$UDel_temp_popweight
  
  # Create precipitation squared
  data$UDel_precip_popweight_2 <- data$UDel_precip_popweight^2
  
  # Create poor country indicator (bottom 50% of GDP per capita)
  data$poorWDIppp <- as.numeric(data$GDPpctile_WDIppp < 50)
  data$poorWDIppp[is.na(data$GDPpctile_WDIppp)] <- NA
  data$interact <- data$poorWDIppp
  
  # Check if _yi_ and _y2_ variables exist
  yi_vars <- grep("^_yi_", names(data), value = TRUE)
  y2_vars <- grep("^_y2_", names(data), value = TRUE)
  
  # Run regression with interactions
  if (length(yi_vars) > 0 & length(y2_vars) > 0) {
    # If trend variables exist in dataset, use backticks
    yi_vars_quoted <- paste0("`", yi_vars, "`")
    y2_vars_quoted <- paste0("`", y2_vars, "`")
    
    formula_str <- paste0(var, " ~ temp + I(temp^2) + ",
                          "UDel_precip_popweight + UDel_precip_popweight_2 + ",
                          paste(yi_vars_quoted, collapse = " + "), " + ",
                          paste(y2_vars_quoted, collapse = " + "),
                          " | iso_id + year | 0 | iso_id")
    
  } else {
    # Create trends on the fly
    formula_str <- paste0(var, " ~ temp + I(temp^2) + ",
                          "UDel_precip_popweight + UDel_precip_popweight_2 | ",
                          "iso_id + year + iso_id:time + iso_id:time2 | 0 | iso_id")
  }
  
  model <- felm(as.formula(formula_str), data = data)
  
  # Generate margins at different temperatures for each group
  temp_seq <- seq(0, 30, by = 1)
  
  # Create prediction data for both groups
  pred_data_list <- lapply(c(0, 1), function(interact_val) {
    pred_df <- data.frame(
      temp = temp_seq,
      interact = interact_val,
      UDel_precip_popweight = mean(data$UDel_precip_popweight, na.rm = TRUE)
    )
    pred_df$UDel_precip_popweight_2 <- pred_df$UDel_precip_popweight^2
    
    
    
    # Create design matrix
    X_pred <- model.matrix(~ interact * (temp + I(temp^2) + 
                                           UDel_precip_popweight + UDel_precip_popweight_2), 
                           data = pred_df)
    
    # Get coefficients that match the design matrix
    all_coefs <- coef(model)
    common_names <- intersect(colnames(X_pred), names(all_coefs))
    #cat("  Common names:", common_names, "\n")
    
    if (length(common_names) == 0) {
      warning(paste("No matching coefficients found for", var))
      return(NULL)
    }
    
    # Subset to matching coefficients
    X_pred <- X_pred[, common_names, drop = FALSE]
    cat("  X pred:", X_pred, "\n")
    coefs <- all_coefs[common_names]
    cat(" coefs :", coefs, "\n")
    cat(" all coefs :", all_coefs, "\n")
    
    # Predictions
    pred_df$estimate <- as.vector(X_pred %*% coefs)
    print(pred_df$estimate)
    
    # Standard errors
    vcov_mat <- vcov(model)
    vcov_subset <- vcov_mat[common_names, common_names, drop = FALSE]
    se_vec <- sqrt(diag(X_pred %*% vcov_subset %*% t(X_pred)))
    pred_df$stderr <- se_vec
    
    # 90% CI
    pred_df$min90 <- pred_df$estimate - 1.645 * pred_df$stderr
    pred_df$max90 <- pred_df$estimate + 1.645 * pred_df$stderr
    
    return(pred_df)
  })
  
  # Combine results
  result_df <- bind_rows(pred_data_list) %>%
    # Change the name of the header
    mutate(x = temp,
           model = var,
           dof = model1$df.residual,  # degrees of freedom
           parm = paste0("_at#", row_number())  # parameter to comply with STATA format
    ) %>%  
    select(x, interact, estimate, stderr, min90, max90, model)
  
  heterogeneity_results[[var]] <- result_df
}

# Combine all results
effect_heterogeneity <- bind_rows(heterogeneity_results)
write.csv(effect_heterogeneity, "data/output/EffectHeterogeneity.csv", 
          row.names = FALSE)


# -----------------  DATA FOR FIGURE 2, PANEL C  -----------------------

# Reload data
data <- read_dta("data/input/GrowthClimateDataset.dta")

# Create time variables
data$time <- data$year - 1985
data$time2 <- data$time^2

# Generate temperature variable
data$temp <- data$UDel_temp_popweight

# Create precipitation squared
data$UDel_precip_popweight_2 <- data$UDel_precip_popweight^2

# Create early period indicator
data$early <- as.numeric(data$year < 1990)
data$interact <- data$early

# Check if _yi_ and _y2_ variables exist
yi_vars <- grep("^_yi_", names(data), value = TRUE)
y2_vars <- grep("^_y2_", names(data), value = TRUE)

# Run regression with time period interactions
if (length(yi_vars) > 0 & length(y2_vars) > 0) {
  # If trend variables exist in dataset, use backticks
  yi_vars_quoted <- paste0("`", yi_vars, "`")
  y2_vars_quoted <- paste0("`", y2_vars, "`")
  
  formula_str <- paste0("growthWDI ~ interact * (temp + I(temp^2) + ",
                        "UDel_precip_popweight + UDel_precip_popweight_2) + ",
                        paste(yi_vars_quoted, collapse = " + "), " + ",
                        paste(y2_vars_quoted, collapse = " + "),
                        " | iso_id + year | 0 | iso_id")
} else {
  # Create trends on the fly
  formula_str <- paste0("growthWDI ~ interact * (temp + I(temp^2) + ",
                        "UDel_precip_popweight + UDel_precip_popweight_2) | ",
                        "iso_id + year + iso_id:time + iso_id:time2 | 0 | iso_id")
}

model <- felm(as.formula(formula_str), data = data)

# Generate margins
temp_seq <- seq(0, 30, by = 1)

# Create prediction data for both time periods
pred_data_list <- lapply(c(0, 1), function(interact_val) {
  pred_df <- data.frame(
    temp = temp_seq,
    interact = interact_val,
    UDel_precip_popweight = mean(data$UDel_precip_popweight, na.rm = TRUE)
  )
  pred_df$UDel_precip_popweight_2 <- pred_df$UDel_precip_popweight^2
  
  # Create design matrix
  X_pred <- model.matrix(~ interact * (temp + I(temp^2) + 
                                         UDel_precip_popweight + UDel_precip_popweight_2), 
                         data = pred_df)
  
  # Get coefficients that match the design matrix
  all_coefs <- coef(model)
  common_names <- intersect(colnames(X_pred), names(all_coefs))
  
  if (length(common_names) == 0) {
    warning("No matching coefficients found for time heterogeneity")
    return(NULL)
  }
  
  # Subset to matching coefficients
  X_pred <- X_pred[, common_names, drop = FALSE]
  coefs <- all_coefs[common_names]
  
  # Predictions
  pred_df$estimate <- as.vector(X_pred %*% coefs)
  
  # Standard errors
  vcov_mat <- vcov(model)
  vcov_subset <- vcov_mat[common_names, common_names, drop = FALSE]
  se_vec <- sqrt(diag(X_pred %*% vcov_subset %*% t(X_pred)))
  pred_df$stderr <- se_vec
  
  # 90% CI
  pred_df$min90 <- pred_df$estimate - 1.645 * pred_df$stderr
  pred_df$max90 <- pred_df$estimate + 1.645 * pred_df$stderr
  
  return(pred_df)
})

# Format results
time_heterogeneity <- bind_rows(pred_data_list) %>%
  mutate(x = temp) %>%
  select(x, interact, estimate, stderr, min90, max90)

write.csv(time_heterogeneity, "data/output/EffectHeterogeneityOverTime.csv", 
          row.names = FALSE)

print("Analysis complete. All outputs saved to data/output/")