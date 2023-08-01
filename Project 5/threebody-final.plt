set title 'The Effect of Jupiter on Earth'
set xlabel "x(AU))"
set ylabel "y(AU)"
set size square
set xzeroaxis
set yzeroaxis
set xrange [-6:6]
set yrange [-6:6]
plot "threebody-final.out" using 2:3 title 'Earth orbit' with l lc rgb "#FF0000", \
     "threebody-final.out" using 4:5 title 'Jupiter orbit' with l lc rgb "#0000FF", \
     
pause -1
 
