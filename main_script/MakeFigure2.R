# SCRIPT TO MAKE FIGURE 2
# This file calls data created in GenerateFigure2Data.R

require(dplyr)

# For any spatial operations that maptools used to provide, use:
# require(sf)         # Modern replacement for spatial operations
# require(terra)      # For raster operations

"%&%" <- function(x, y) paste(x, y, sep = "")  # define a function for easy string pasting

# Create output directory if needed
dir.create("figures/MainFigs_Input", recursive = TRUE, showWarnings = FALSE)
pdf(file = "figures/MainFigs_Input/Figure2.pdf", width = 10, height = 5.5, useDingbats = FALSE)

mat <- matrix(c(1, 1, 2, 3, 1, 1, 4, 5), nrow = 2, byrow = TRUE)
layout(mat)
ax <- 1.5  # scaling for axes
par(mar = c(4, 4, 2, 1))

########################################################
#  Panel A
########################################################

resp <- read.csv("data/output/estimatedGlobalResponse.csv")
dta <- read.csv("data/output/mainDataset.csv")
smpl <- is.na(dta$growthWDI) == FALSE & is.na(dta$UDel_temp_popweight) == FALSE   # main estimation sample
coef <- read.csv("data/output/estimatedCoefficients.csv")

# center response at optimum
x <- resp$x
mx <- max(resp$estimate, na.rm = TRUE)
est <- resp$estimate - mx
min90 <- resp$min90 - mx
max90 <- resp$max90 - mx

# Remove any NA values to ensure consistent lengths
valid_idx <- !is.na(x) & !is.na(est) & !is.na(min90) & !is.na(max90)
x <- x[valid_idx]
est <- est[valid_idx]
min90 <- min90[valid_idx]
max90 <- max90[valid_idx]

ctys <- c('USA', 'CHN', 'DEU', 'JPN', 'IND', 'NGA', 'IDN', 'BRA', 'FRA', 'GBR')
ctt <- c('US', 'CHN', "GER", "JPN", 'IND', 'NGA', 'INDO', 'BRA', 'FRA', 'UK')

# initialize plot
plot(1, xlim = c(-2, 30), ylim = c(-0.4, 0.1), type = "n", las = 1, cex.axis = 1.3,
     xlab = "Annual average temperature (°C)", ylab = expression(Delta*" in ln(GDP/cap)"),
     main = "Panel A: Global Temperature vs. GDP/cap Growth")

# plot CI and main effect
polygon(c(x, rev(x)), c(min90, rev(max90)), col = "lightblue", border = NA)
lines(x, est, lwd = 2)

# add vertical average temperature lines for selected countries
for (j in 1:length(ctys)) {
  tt <- mean(dta$UDel_temp_popweight[dta$iso == ctys[j]], na.rm = TRUE)
  segments(tt, -0.23, tt, 0.15, lwd = 0.5)
}

# Add histograms at bottom
# first calculate percent of population and global gdp produced at each temperature bin, for our estimation sample
bins <- seq(-7, 30, 0.5)
histtemp <- dta$UDel_temp_popweight[smpl]
histpop <- dta$Pop[smpl]
histgdp <- dta$TotGDP[smpl]
pop <- gdp <- c()
for (j in 1:(length(bins) - 1)) {
  lo <- bins[j]
  hi <- bins[j + 1]
  pop <- c(pop, sum(histpop[histtemp >= lo & histtemp < hi]))
  gdp <- c(gdp, sum(histgdp[histtemp >= lo & histtemp < hi]))
}
pop <- pop / sum(pop)
gdp <- gdp / sum(gdp)

# parameters that set where histograms go
dis <- 0.055
base <- -0.3

# now make histograms
# temperature
zz <- hist(histtemp, plot = FALSE, breaks = bins)
cts <- zz$counts / max(zz$counts) * 0.05  # sets the height of the tallest bar to 0.05
rect(bins, base, bins + 0.5, base + cts, col = "red")
text(-1, base + 0.025, "Global distribution of temperature observations", pos = 4, cex = 1)
# population
cts <- pop / max(pop) * 0.05
rect(bins, base - dis * (1), bins + 0.5, base - dis * (1) + cts, col = "grey")
text(-1, base - dis * (1) + 0.025, "Global distribution of population", pos = 4, cex = 1)
# GDP
cts <- gdp / max(gdp) * 0.05
rect(bins, base - dis * (2), bins + 0.5, base - dis * (2) + cts, col = "black")
text(-1, base - dis * (2) + 0.025, "Global distribution of GDP", pos = 4, cex = 1)

########################################################
#  Panel B - Rich vs Poor Countries
########################################################
resp <- read.csv("data/output/EffectHeterogeneity.csv")
poor <- dta$GDPpctile_WDIppp < 50
rich <- dta$GDPpctile_WDIppp >= 50

resp <- resp[resp$x >= 5, ]  # dropping estimates below 5C, since so little poor country exposure down there
mods <- unique(as.character(resp$model))

m <- "growthWDI"
plot(1, xlim = c(5, 30), ylim = c(-0.35, 0.1), type = "n", las = 1, 
     cex.axis = 1.3, cex.lab = 1.3, ylab = "Change in growth rate", 
     xlab = "Annual average temperature (°C)",
     main = "Rich vs Poor Countries")

smp <- resp$model == m & resp$interact == 1  # poor countries
xx <- resp$x[smp]
mx <- max(resp$estimate[smp], na.rm = TRUE)
est <- resp$estimate[smp] - mx
min90 <- resp$min90[smp] - mx
max90 <- resp$max90[smp] - mx

# Remove NA values
valid_idx <- !is.na(xx) & !is.na(est) & !is.na(min90) & !is.na(max90)
xx <- xx[valid_idx]
est <- est[valid_idx]
min90 <- min90[valid_idx]
max90 <- max90[valid_idx]

polygon(c(xx, rev(xx)), c(min90, rev(max90)), col = "lightblue", border = NA)
lines(xx, est, lwd = 2, col = "steelblue3")

# now add rich countries
smp <- resp$model == m & resp$interact == 0  # rich countries
xx <- resp$x[smp]
mx <- max(resp$estimate[smp], na.rm = TRUE)
est <- resp$estimate[smp] - mx

# Remove NA values
valid_idx <- !is.na(xx) & !is.na(est)
xx <- xx[valid_idx]
est <- est[valid_idx]

lines(xx, est, lwd = 2, col = "red")

# now add histograms of temperature exposures at the base
bins <- seq(-7, 30, 0.5)
poortemp <- dta$UDel_temp_popweight[smpl == TRUE & poor == TRUE]
richtemp <- dta$UDel_temp_popweight[smpl == TRUE & rich == TRUE]
base <- -0.3
zz <- hist(richtemp, plot = FALSE, breaks = bins)
cts <- zz$counts / max(zz$counts) * 0.05  # sets the height of the tallest bar to 0.05
rect(bins, base, bins + 0.5, base + cts, border = "red", col = NA)
text(4, base + 0.055, "rich country temperatures", pos = 4, cex = 0.8)
base <- -0.35
zz1 <- hist(poortemp, plot = FALSE, breaks = bins)
cts <- zz1$counts / max(zz1$counts) * 0.05
rect(bins, base, bins + 0.5, base + cts, col = "lightblue")
text(4, base + 0.03, "poor country temperatures", pos = 4, cex = 0.8)



########################################################
#  Panel C - Early vs Late Period
########################################################
resp <- read.csv("data/output/EffectHeterogeneityOverTime.csv")
early <- dta$year < 1990

smp <- resp$interact == 1  # early period
xx <- resp$x[smp]
mx <- max(resp$estimate[smp], na.rm = TRUE)
est <- resp$estimate[smp] - mx	
min90 <- resp$min90[smp] - mx
max90 <- resp$max90[smp] - mx

# Remove NA values
valid_idx <- !is.na(xx) & !is.na(est) & !is.na(min90) & !is.na(max90)
xx <- xx[valid_idx]
est <- est[valid_idx]
min90 <- min90[valid_idx]
max90 <- max90[valid_idx]

plot(1, xlim = c(5, 30), ylim = c(-0.35, 0.1), type = "n", las = 1, 
     cex.axis = 1.3, cex.lab = 1.3, ylab = "Change in growth rate", 
     xlab = "Annual average temperature (°C)",
     main = "Early vs Late Period")

polygon(c(xx, rev(xx)), c(min90, rev(max90)), col = "lightblue", border = NA)
lines(xx, est, lwd = 2, col = "steelblue3")

# now add point estimate for later period
smp <- resp$interact == 0  # later period
xx <- resp$x[smp]
mx <- max(resp$estimate[smp], na.rm = TRUE)
est <- resp$estimate[smp] - mx

# Remove NA values
valid_idx <- !is.na(xx) & !is.na(est)
xx <- xx[valid_idx]
est <- est[valid_idx]

lines(xx, est, lwd = 2, col = "red")

# now add histograms of temperature exposures at the base
bins <- seq(-7, 30, 0.5)
earlytemp <- dta$UDel_temp_popweight[smpl == TRUE & early == TRUE]
latetemp <- dta$UDel_temp_popweight[smpl == TRUE & early == FALSE]
base <- -0.3
zz <- hist(earlytemp, plot = FALSE, breaks = bins)
cts <- zz$counts / max(zz$counts) * 0.05  # sets the height of the tallest bar to 0.05
rect(bins, base, bins + 0.5, base + cts, border = "red", col = NA)
text(4, base + 0.025, "1960-1989 temperatures", pos = 4, cex = 0.8)
base <- -0.35
zz1 <- hist(latetemp, plot = FALSE, breaks = bins)
cts <- zz1$counts / max(zz1$counts) * 0.05
rect(bins, base, bins + 0.5, base + cts, col = "lightblue")
text(4, base + 0.025, "1990-2010 temperatures", pos = 4, cex = 0.8)


########################################################
#  Panels D, E - Agricultural vs Non-Agricultural GDP
########################################################

resp <- read.csv("data/output/EffectHeterogeneity.csv")
poor <- dta$GDPpctile_WDIppp < 50
rich <- dta$GDPpctile_WDIppp >= 50
resp <- resp[resp$x >= 5, ]  # dropping estimates below 5C, because so little poor country exposure there
mods <- unique(as.character(resp$model))
toplot <- c("AgrGDPgrowthCap", "NonAgrGDPgrowthCap")
panel_titles <- c("Agricultural GDP", "Non-Agricultural GDP")

for (i in 1:length(toplot)) {
  m <- toplot[i]
  plot(1, xlim = c(5, 30), ylim = c(-0.35, 0.1), type = "n", las = 1, 
       cex.axis = 1.3, cex.lab = 1.3, ylab = "Change in growth rate", 
       xlab = "Annual average temperature (°C)",
       main = panel_titles[i])
  
  smp <- resp$model == m & resp$interact == 1  # poor countries
  xx <- resp$x[smp]
  mx <- max(resp$estimate[smp], na.rm = TRUE)
  est <- resp$estimate[smp] - mx  
  min90 <- resp$min90[smp] - mx
  max90 <- resp$max90[smp] - mx
  
  # Remove NA values
  valid_idx <- !is.na(xx) & !is.na(est) & !is.na(min90) & !is.na(max90)
  xx <- xx[valid_idx]
  est <- est[valid_idx]
  min90 <- min90[valid_idx]
  max90 <- max90[valid_idx]
  
  polygon(c(xx, rev(xx)), c(min90, rev(max90)), col = "lightblue", border = NA)
  lines(xx, est, lwd = 2, col = "steelblue3")
  
  # now add rich countries
  smp <- resp$model == m & resp$interact == 0  # rich countries
  xx <- resp$x[smp]
  mx <- max(resp$estimate[smp], na.rm = TRUE)
  est <- resp$estimate[smp] - mx
  
  # Remove NA values
  valid_idx <- !is.na(xx) & !is.na(est)
  xx <- xx[valid_idx]
  est <- est[valid_idx]
  
  lines(xx, est, lwd = 2, col = "red")
  
  legend("bottomleft", legend = c("Poor countries", "Rich countries"), 
         col = c("steelblue3", "red"), lwd = 2, bty = "n")
}

dev.off()

cat("Figure 2 saved to figures/MainFigs_Input/Figure2.pdf\n")