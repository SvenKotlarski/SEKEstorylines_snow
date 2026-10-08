#!/usr/bin/env Rscript

results <- read.delim("../RESULTS/snowfall_sorted.txt")
top5 <- head(results, 5)
reds <- c("#67000D", "#A50F15", "#CB181D", "#EF3B2C", "#FB6A4A")

png("../RESULTS/snowfall_sorted.png", width = 1800, height = 1400, res = 180)
layout(1:2, heights = c(4, 1.5))
par(mar = c(4, 4.5, 3, 1))
plot(
  range(results$snowfall), c(0, 1),
  type = "n",
  xaxt = "n",
  yaxt = "n",
  xlab = "Mean seasonal snowfall (mm/day)",
  ylab = "",
  main = "Mean NOV-APR snowfall in CH above 1000 m (GWL3)"
)
axis(1, at = seq(0.4, 4.4, by = 0.2))
segments(results$snowfall, 0, results$snowfall, 1, col = "grey60")
segments(top5$snowfall, 0, top5$snowfall, 1, col = reds, lwd = 2)

# Legend table below the plot
par(mar = c(0.5, 1, 0, 1))
plot.new()
plot.window(xlim = c(0, 1), ylim = c(0, 6))
x <- c(0.08, 0.18, 0.40, 0.50)
text(x, 5.65, c("Rank", "Snowfall (mm/day)", "Year", "Model"), adj = 0, font = 2)
y <- 5:1 - 0.35
points(rep(0.045, 5), y, pch = 16, col = reds)
text(x[1], y, 1:5, adj = 0)
text(x[2], y, top5$snowfall, adj = 0)
text(x[3], y, top5$year, adj = 0)
text(x[4], y, top5$model, adj = 0)
dev.off()
