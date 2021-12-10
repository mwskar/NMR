set title "Filter Results"
set xlabel "x"
set ylabel "y"
plot "ordered.dat" using 1:2 with lines title "Ordered"
replot "filt.dat" using 1:2 with lines title "Filter"
set output "plot.ps"
set terminal postscript enhanced color landscape
replot
set output "plot.png"
set terminal png
replot

