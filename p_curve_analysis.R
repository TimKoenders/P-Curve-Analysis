# Introduction ------------------------------------------------------------
# This is the R Code behind the p-curve app 4.052
# Written by Uri Simonsohn (urisohn@gmail.com) / Adapted by Tim Koenders (tim.koenders@wu.ac.at)
# To run it you need to store the p-values in a column vector in a csv
# Then, you just run the function on that file. The app generates various text files with the results in table format and saves the figures 
# as .PNG files on the same folder as where you put the .txt file. 


# Clean_up and packages ------------------------------------------------------------------
rm(list = ls())
gc()

pacman::p_load(
  tidyverse,
  dplyr,
  stringr,
  poibin,
  scales
)


# Set your working directory ----------------------------------------------
setwd("C:\\Users\\koend\\OneDrive\\Bureaublad\\WU 2023-2024\\Institute_Cognition\\P-curve analysis\\Routput")


# Create data for p-curve analysis ------------------------------------------------------------
# load the file with the p_values from python
df <- read_csv("p_values.csv")
# double check p-values who are over multiple lines

# Clean p-values
d <- df %>%
  separate(p_value, into = c("NA", "p_value"), sep = "([<=>])\\s*") %>%
  mutate(p_value = str_remove(p_value, ".")) %>%
  mutate(p_value = ifelse(p_value == "000", "001", p_value)) %>%
  mutate(p_value = as.numeric(p_value) / 1000)

# Create a dataframe of p-values
p_values_vector <- d$p_value
p_values_df <- data.frame(p_value = p_values_vector)

# Generate the output filename 
output_filename <- "p_curve_analysis"


# P_curve analysis --------------------------------------------------------
# Function to perform p-curve analysis with p-values as input
pcurveAnalysis <- function(p_values_vector, filek= output_filename) {
  ksig <- length(p_values_vector)  # Total number of p-values
  khalf <- 0  # Number of p-values in the range [0.025, 0.05] (modify this if needed)
  
  # Define functions for calculations
  prop33 <- function(p) {
    return(p_values_vector <= p)
  }
  
  # Calculate power and confidence interval
  hat <- mean(p_values_vector <= 0.05)
  power.ci.lb <- hat - 1.96 * sqrt(hat * (1 - hat) / ksig)
  power.ci.ub <- hat + 1.96 * sqrt(hat * (1 - hat) / ksig)
  
  # Save power estimate and confidence intervals to a text file
  power_results <- c(power.ci.lb, hat, power.ci.ub)
  write(power_results, paste("POWERHAT_", filek, ".txt", sep = ""), sep = "\n")
  
  # Calculate expected p-curve for 33% power
  gcdf1 <- prop33(0.01)
  gcdf2 <- prop33(0.02)
  gcdf3 <- prop33(0.03)
  gcdf4 <- prop33(0.04)
  green1 <- mean(gcdf1, na.rm = TRUE) * 3
  green2 <- mean(gcdf2 - gcdf1, na.rm = TRUE) * 3
  green3 <- mean(gcdf3 - gcdf2, na.rm = TRUE) * 3
  green4 <- mean(gcdf4 - gcdf3, na.rm = TRUE) * 3
  green5 <- mean(1 / 3 - gcdf4, na.rm = TRUE) * 3
  green <- 100 * c(green1, green2, green3, green4, green5)
  
  # Put p-values into bins (0.01 to 0.05)
  ps <- ceiling(p_values_vector * 100) / 100
  blue <- c()
  for (i in c(0.01, 0.02, 0.03, 0.04, 0.05)) {
    blue <- c(blue, sum(ps == i, na.rm = TRUE) / ksig * 100)
  }
  
  # Red line
  red <- c(20, 20, 20, 20, 20)
  
  # Create the plot
  png(filename = paste(filek, ".png", sep = ""), width = 2600, height = 2400, res = 400)
  x <- c(0.01, 0.02, 0.03, 0.04, 0.05)
  par(mar = c(6, 5.5, 1.5, 3))
  moveup <- max(max(blue[2:5]) - 66, 0)
  ylim <- c(0, 105 + moveup)
  legend.top <- 100 + moveup
  plot(x, blue, type = 'l', col = 'dodgerblue2', main = "",
       lwd = 2, xlab = "", ylab = "", xaxt = "n", yaxt = "n", xlim = c(0.01, 0.051),
       ylim = ylim, bty = 'L', las = 1, axes = F)
  x_ <- c(".01", ".02", ".03", ".04", ".05")
  axis(1, at = x, labels = x_)
  y_ <- c("0%", "25%", "50%", "75%", "100%")
  y <- c(0, 25, 50, 75, 100)
  axis(2, at = y, labels = y_, las = 1, cex.axis = 1.2)
  mtext("Percentage of test results", font = 2, side = 2, line = 3.85, cex = 1.25)
  mtext("p            ", font = 4, side = 1, line = 2.3, cex = 1.25)
  points(x, blue, type = "p", pch = 20, bg = "dodgerblue2", col = "dodgerblue2")
  text(x + 0.00075, blue + 3.5, scales::percent(round(blue) / 100), col = 'black', cex = 0.75)
  lines(x, red, type = 'l', col = 'firebrick2', lwd = 1.5, lty = 3)
  lines(x, green, type = 'l', col = 'springgreen4', lwd = 1.5, lty = 5)
  tab1 = 0.017
  tab2 = tab1 + 0.0015
  gap1 = 9
  gap2 = 4
  font.col = 'gray44'
  text.blue = paste0("Power estimate: ", percent(hat), ", CI(",
                     percent(power.ci.lb), ",",
                     percent(power.ci.ub), ")")
  text(tab1, legend.top, adj = 0, cex = 0.85, bquote("Observed " * italic(p) * "-curve"))
  text(tab2, legend.top - gap2, adj = 0, cex = 0.68, text.blue, col = font.col)
  text.red = bquote("Null of no effect: " * italic(p))
  text(tab1, legend.top - gap1, adj = 0, cex = 0.85, text.red)
  text(tab2, legend.top - gap1 - gap2, adj = 0, cex = 0.68, text.red, col = font.col)
  text.green = bquote("Null of 33% power: " * italic(p))
  text(tab1, legend.top - 2 * gap1, adj = 0, cex = 0.85, text.green)
  text(tab2, legend.top - 2 * gap1 - gap2, adj = 0, cex = 0.68, text.green, col = font.col)
  segments(x0 = tab1 - 0.005, x1 = tab1 - 0.001, y0 = legend.top, y1 = legend.top, col = 'dodgerblue2', lty = 1, lwd = 1.5)
  segments(x0 = tab1 - 0.005, x1 = tab1 - 0.001, y0 = legend.top - gap1, y1 = legend.top - gap1, col = 'firebrick2', lty = 3, lwd = 1.5)
  segments(x0 = tab1 - 0.005, x1 = tab1 - 0.001, y0 = legend.top - 2 * gap1, y1 = legend.top - 2 * gap1, col = 'springgreen4', lty = 2, lwd = 1.5)
  rect(tab1 - 0.0065, legend.top - 2 * gap1 - gap2 - 3, tab1 + 0.032, legend.top + 3, border = 'gray85')
  msgx = bquote("Note: The observed " * italic(p) * "-curve includes " * ksig *
                  " statistically significant (" * italic(p) * " < .05) results.")
  mtext(msgx, side = 1, line = 4, cex = 0.65, adj = 0)
  kns = ksig
  ns_msg = bquote("There were no non-significant results entered.")
  if (kns == 1) ns_msg = bquote("There was one additional result entered but excluded from " * italic(p) * "-curve because it was " * italic(p) * " > .05.")
  if (kns > 1)  ns_msg = bquote("There were " * kns * " additional results entered but excluded from " * italic(p) * "-curve because they were " * italic(p) * " > .05.")
  mtext(ns_msg, side = 1, line = 4.75, cex = 0.65, adj = 0)
  dev.off()
  
  # Save p-value calculations to text files
  table_calc = data.frame(p_value = p_values_vector, observed_p_curve = p_values_vector <= 0.05)
  write.table(table_calc, sep = "\t", row.names = FALSE, file = paste("Calculations_", filek, ".txt", sep = ""))
  
  # Save results behind p-curve figure
  headers2 = c("p-value", "Observed (blue)", "Power 33% (Green)", "Null of no effect (Red)")
  table_figure = data.frame(x, blue, green, red)
  write.table(table_figure, sep = "\t", row.names = FALSE, file = paste("FigNumbers_", filek, ".txt", sep = ""))
  
  # Cumulative p-curves
  dropk = function(pp, k, droplow) {
    pp = pp[!is.na(pp)]
    n = length(pp)
    pp = sort(pp)
    if (k == 0) ppk = pp
    if (droplow == 1 & k > 0) {
      ppk = pp[(1 + k):n]
      ppmin = min(pp[k], k / (n + 1))
      ppk = (ppk - ppmin) / (1 - ppmin)
    }
    if (droplow == 0 & k > 0) {
      ppk = pp[1:(n - k)]
      ppmax = max(pp[n - k + 1], (n - k) / (n + 1))
      ppk = ppk / ppmax
    }
    ppk = pmax(ppk, 0.00001)
    ppk = pmin(ppk, 0.99999)
    Z = sum(qnorm(ppk)) / sqrt(n - k)
    return(pnorm(Z))
  }
  
  droplow.r = droplow.33 = drophigh.r = drophigh.33 = c()
  
  for (i in 0:(round(ksig / 2) - 1)) {
    droplow.r = c(droplow.r, dropk(pp = p_values_vector, k = i, droplow = 1))
    drophigh.r = c(drophigh.r, dropk(pp = p_values_vector, k = i, droplow = 0))
    droplow.33 = c(droplow.33, dropk(pp = p_values_vector, k = i, droplow = 1))
    drophigh.33 = c(drophigh.33, dropk(pp = p_values_vector, k = i, droplow = 0))
  }
  
  plotdrop = function(var, col) {
    k = length(var)
    plot(0:(k - 1), var, xlab = "", ylab = "", type = "b", yaxt = "n", xaxt = "n", main = "",
         cex.main = 1.15, ylim = c(0, 1), col = col)
    points(0, var[1], pch = 19, cex = 1.6)
    abline(h = 0.05, col = "red")
    axis(2, c(0.05, 2:9 / 10), labels = c('.05', '.2', '.3', '.4', '.5', '6', '7', '.8', '.9'), las = 1, cex.axis = 1.5)
    axis(1, c(0:(k - 1)), las = 1, cex.axis = 1.4)
  }
  
  png(filename = paste(filek, "_cumulative.png", sep = ""), width = 4000, height = 4000, res = 400)
  
  par(mfrow = c(3, 2), mar = c(4, 3, 0, 2), mgp = c(2.5, 1, 0), oma = c(5, 14, 5, 1))
  
  plotdrop(droplow.r, col = "dodgerblue2")
  mtext(side = 2, line = 4.5, bquote(italic(P) * "-value of overall test"), font = 2, cex = 1.25)
  mtext(side = 2, line = 3, "(Stouffer's method)", font = 3, cex = 1)
  mtext("Right skew\n\n", line = 7, side = 2, cex = 1.2, las = 1.25, col = "dodgerblue2")
  mtext(bquote("(Full " * italic(p) * "-curve)"), line = 7, side = 2, cex = 1, las = 1, col = "dodgerblue2")
  
  plotdrop(drophigh.r, col = "firebrick2")
  mtext(side = 2, line = 4.5, bquote(italic(P) * "-value of overall test"), font = 2, cex = 1.25)
  mtext(side = 2, line = 3, "(Stouffer's method)", font = 3, cex = 1)
  mtext("Left skew\n\n", line = 7, side = 2, cex = 1.2, las = 1.25, col = "firebrick2")
  mtext(bquote("(Full " * italic(p) * "-curve)"), line = 7, side = 2, cex = 1, las = 1, col = "firebrick2")
  
  plotdrop(droplow.33, col = "dodgerblue2")
  mtext(side = 2, line = 4.5, bquote(italic(P) * "-value of overall test"), font = 2, cex = 1.25)
  mtext(side = 2, line = 3, "(Stouffer's method)", font = 3, cex = 1)
  mtext("Right skew\n\n", line = 7, side = 2, cex = 1.2, las = 1.25, col = "dodgerblue2")
  mtext(bquote("(33% power " * italic(p) * "-curve)"), line = 7, side = 2, cex = 1, las = 1, col = "dodgerblue2")
  
  plotdrop(drophigh.33, col = "firebrick2")
  mtext(side = 2, line = 4.5, bquote(italic(P) * "-value of overall test"), font = 2, cex = 1.25)
  mtext(side = 2, line = 3, "(Stouffer's method)", font = 3, cex = 1)
  mtext("Left skew\n\n", line = 7, side = 2, cex = 1.2, las = 1.25, col = "firebrick2")
  mtext(bquote("(33% power " * italic(p) * "-curve)"), line = 7, side = 2, cex = 1, las = 1, col = "firebrick2")
  
  dev.off()
}
pcurveAnalysis(p_values_vector, filek= output_filename)
