#Jack Messerli-Wallace Project 4

set st da line
set key top left
set xzeroaxis
set mxtics 2
#set xtics font ", 8"
set yzeroaxis
set mytics 2

#1 is time, 2 is phi, 3 is omega
#set title ''
#set xlabel 'time [s]'
#set ylabel 'phi [rad]'
#plot 'data.out' u 1:2 title 'phi vs. time' w p ps 0.3
#pause -1
set title 'Phi vs Omega'
set size square
set xlabel 'phi [rad]'
set ylabel 'omega [rad/s]'
set xrange [:]
set yrange [:]
plot 'data.out' u 2:3 w p ps 0.1 notitle
#pause -1 #pause is not needed on windows gnuplot