# Jack Messerli-Wallace Project 1
# Some gnuplot installation go with points as the default.  
set st da points
# One can change the font, but let's take the defaults for now.
set xlabel 'x position (steps)'
set ylabel 'y position (steps)'
# Let's cchange the default position of the legend.
set key top left
# The next statement plots the zero axis.
set xzeroaxis
set mxtics 2
set yzeroaxis
set mytics 2
#set x and y range
set xrange [-450:450]
set yrange [-450:450]

set title 'Milk after 100 Steps'
plot 'step.100' u 2:3 t '100 Steps' w p pt 7 ps 0.2
# The next statement is needed to keep your beautiful picture on the screen until you type something else.
pause -1
# plot next
set title 'Milk after 200 Steps'
plot 'step.200' u 2:3 t '200 Steps' w p pt 7 ps 0.2
pause -1
#
set title 'Milk after 500 Steps'
plot 'step.500' u 2:3 t '500 Steps' w p pt 7 ps 0.2
pause -1
#
set title 'Milk after 1000 Steps'
plot 'step.1000' u 2:3 t '1000 Steps' w p pt 7 ps 0.2
pause -1
#
set title 'Milk after 5000 Steps'
plot 'step.5000' u 2:3 t '5000 Steps' w p pt 7 ps 0.2
pause -1
#
set title 'Milk after 10000 Steps'
plot 'step.10000' u 2:3 t '10000 Steps' w p pt 7 ps 0.2
pause -1