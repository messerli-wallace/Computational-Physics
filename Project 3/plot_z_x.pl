#Jack Messerli-Wallace Project 2

#square plot
#set size square 
set size ratio 0.7
#lines
set st da line
# One can change the font, but let's take the defaults for now.
set xlabel 'x position (m)'
set ylabel 'z position (m)'
# Let's cchange the default position of the legend.
set key top right
# The next statement plots the zero axis.
set xzeroaxis
set mxtics 2
set yzeroaxis
set mytics 5
#set x and y range
#set xrange [-450:450]
#set yrange [0.9:5] 
set xrange [:]
set yrange [-0.1:]


set title 'Trajectory of the Baseball (x vs. z)'
plot 'baseball_3.666' u 1:3 title 'Launched with a spin of 168.5 rad/s'
# The next statement is needed to keep your beautiful picture on the screen until you type something else.
pause -1