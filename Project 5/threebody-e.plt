set title 'The Effect of Jupiter (1000x) on Earth'
set xlabel "x(AU)"
set ylabel "y(AU)"
set size square
set xzeroaxis
set yzeroaxis
set xrange [-6:6]
set yrange [-6:6]
plot "threebody-final.out.1000" using 2:3 title 'Earth orbit' with l lc rgb "#FF0000", \
     "threebody-final.out.1000" using 4:5 title 'Jupiter orbit' with l lc rgb "#0000FF", \
     "00_point" u 1:2 t 'Center (Sun)' w p lt 1 pt 7 lc rgb "#000000"
pause -1
 
