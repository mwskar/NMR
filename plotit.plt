set title "Fit Results"
set xlabel "x"
set ylabel "y"
plot "analysis.txt" using 1:2 with lines title "Analysis"
replot 1650 with lines title "Baseline"
set output "plot.ps"
set terminal postscript enhanced color landscape
replot
set output "plot.png"
set terminal png
replot

