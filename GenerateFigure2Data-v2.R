# DATA TO CONSTRUCT RESPONSE FUNCTIONS IN FIGURE 2

library(haven)
library(lfe) # for fixed effect 
#library(marginaleffects)
library(broom)
library(dplyr)

# ---------- ESTIMATE GLOBAL RESPONSE WITH BASELINE REGRESSION SPECIFICATION
########################################################
#  CONSTRUCTION OF FIGURE 2 Panel A
########################################################

# Load data
data <- read_dta("data/input/GrowthClimateDataset.dta")
# Generate variables
data$temp <- data$UDel_temp_popweight
data$temp_sq <- data$UDel_temp_popweight^2
data$precip <- data$UDel_precip_popweight
data$precip_sq <- data$UDel_precip_popweight^2

# Generate fixed effects non-linear model in the form of: y = ax^2 + bx + c 
# Dependent Variable: growthWDI (economic growth rate)
# Independent Variables: temp (Temperature)
#                        temp_sq (Temperature squared)
#                        precip (Precipitation level)
#                        precip_sq (Precipitation level squared)
#                        i.year (Year fixed effects, controls for global trends)
#                        `_yi_` (Country-specific time trends - linear)
#                        `_y2_` (Country-specific time trends - quadratic)
#                        i.iso_id (Country fixed effects)
#                        cluster(iso_id) (Standard errors clustered by country)

# Patterns that match all variables that start with _yi_ and _y2_
yi_vars <- grep("^_yi_", names(data), value = TRUE)
y2_vars <- grep("^_y2_", names(data), value = TRUE)
# include _yi_ and _y2_ variables with back ticks
yi_vars_quoted <- paste0("`", yi_vars, "`")
y2_vars_quoted <- paste0("`", y2_vars, "`")

# formula string with predictors and fixed effect variables 
formula_str <- paste0(
  "growthWDI ~ temp + temp_sq + precip + precip_sq + ",
  paste(yi_vars_quoted, collapse = " + "), " + ",
  paste(y2_vars_quoted, collapse = " + "),
  " | year + iso_id | 0 | iso_id"
)
# Fixed effect linear model method with above independent variables
model <- felm(as.formula(formula_str), data = data)
summary(model)

# Extract model coefficients
coefs <- coef(model)
coef_temp <- coef(model)["temp"]
coef_temp_sq <- coef(model)["temp_sq"]

# Calculate optimal temperature
optimal_temp <- -coef_temp / (2 * coef_temp_sq)   # x-value of vertex: -b/2a 
print(paste("Optimal temperature:", optimal_temp))


# ---------- marginal analysis
temp_seq <- seq(-5, 35, by = 1)
# Generate predictions
# Create prediction data - use mean values for precipitation
pred_data <- data.frame(
  temp = temp_seq,
  temp_sq = temp_seq^2
  #precip = mean(data$precip, na.rm = TRUE)
)
#pred_data$precip_sq <- pred_data$precip^2

# Get coefficients from model
all_coefs <- coef(model)
coef_names <- c("temp", "temp_sq", "precip", "precip_sq")
# Create design matrix matching the available coefficients
#X_pred <- model.matrix(~ temp + temp_sq + precip + precip_sq 
#                     ,data = pred_data)
X_pred <- model.matrix(~ temp + temp_sq, data = pred_data)


common_names <- intersect(colnames(X_pred), names(coefs))
X_pred <- X_pred[, common_names, drop = FALSE]
coefs <- coefs[common_names]
cat(" coefs :", coefs, "\n")

# Calculate predictions
pred_data$estimate <- as.vector(X_pred %*% coefs)

# Calculate standard errors
vcov_mat <- vcov(model)
vcov_subset <- vcov_mat[common_names, common_names, drop = FALSE]
print(vcov_subset)
pred_data$stderr <- sqrt(diag(X_pred %*% vcov_subset %*% t(X_pred)))

# 90% confidence intervals (z = 1.645)
z_crit <- qnorm(0.90)
pred_data$min90 <- pred_data$estimate - z_crit  * pred_data$stderr
pred_data$max90 <- pred_data$estimate + z_crit  * pred_data$stderr

global_response <- pred_data %>%
  # Change the name of the header
  mutate(dof = model$df.residual,  # degrees of freedom
         parm = paste0("_at#", row_number())  # parameter to comply with STATA format
  ) %>%  
  select(temp, estimate, stderr, min90, max90, dof, parm)
write.csv(global_response, "data/output/estimatedGlobalResponse.csv", row.names = FALSE)

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

# Save all coefs in a csv file
write.csv( tidy(model) , "data/output/allCoefficients.csv" )

write.csv(coefficients_df, "data/output/estimatedCoefficients.csv", 
          row.names = FALSE)


