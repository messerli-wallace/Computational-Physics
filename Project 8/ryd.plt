set xlabel 'r'
# set ylabel 'Energy'

set title 'Rydberg Solution'
set key right bottom

# set xrange [85:110]
# set yrange [-1:1]


# plot 'wfn.out' u 1:4 w lines title 'n=2',\
# 'wfn.out' u 1:10 w lines title 'n=8'
# pause -1
# plot 'wfn.out' u 1:8 w lines title 'n=6',\
# 'wfn.out' u 1:12 w lines title 'n=10'
# pause -1

plot 'wfn.out' u 1:15 w lines title 'n=13'