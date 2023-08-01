set title 'Orbit of The Three Planets, Jupiter 1000x'
set xlabel "x(AU)"
set ylabel "y(AU)"
set size square
set xzeroaxis
set yzeroaxis
#set xrange [-6:6]
#set yrange [-6:6]
set term png
set output 'output.png'
plot "threebody.out" using 2:3 title 'Earth orbit' with l lc rgb "#FF0000", \
     "threebody.out" using 4:5 title 'Jupiter orbit' with l lc rgb "#0000FF", \
     "threebody.out" using 6:7 title 'Sun orbit' with l lc rgb "#000000"
#set xrange [-0.006:0.006]
#set yrange [-0.006:0.006]
#set title 'Orbit of the Sun'
#set output 'output_sun.png'
#plot "threebody.out" using 6:7 title 'Sun orbit' with l lc rgb "#000000"
