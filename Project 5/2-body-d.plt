set size square
set xr [-1.6:1.6]; set yr [-1.6:1.6]
set xzeroaxis; set yzeroaxis
set xlabel 'x [AU]'
set ylabel 'y [AU]'
set title "Orbit of Earth"
plot 'earth.out' u 2:3 ps 0.05, \
'earth.out.b.1.9' u 2:3 ps 0.05, \
'earth.out.b.2.1' u 2:3 ps 0.05
pause -1