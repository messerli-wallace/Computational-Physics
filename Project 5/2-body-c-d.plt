set size square
set xr [-0.55:0.55]; set yr [-0.55:0.55]
set xzeroaxis; set yzeroaxis
set xlabel 'x [AU]'
set ylabel 'y [AU]'
set title "Orbit of Earth"
plot 'earth.out' u 2:3 ps 0.1
pause -1