set xlabel 'r'
# set ylabel 'Energy'

set title 'Harmonic Oscillator Solution'
set key left bottom

set xrange [85:110]
set yrange [-1:1]


plot 'wfn.out' u 1:3 w lines title 'n=1',\
'wfn.out' u 1:5 w lines title 'n=3'
pause -1
plot 'wfn.out' u 1:8 w lines title 'n=6',\
'wfn.out' u 1:12 w lines title 'n=10'
pause -1