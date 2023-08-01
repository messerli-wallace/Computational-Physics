set xlabel 'r'

set title 'Coulomb Potential Solution'
set key right bottom

plot 'wfn.out' u 1:3 w lines title 'l=0, n=1',\
'wfn.out' u 1:5 w lines title 'l=0, n=2'
# 'wfn.out' u 1:8 w lines title 'l=0, n=5'
pause -1
plot 'wfn.out' u 1:15 w lines title 'l=1, n=5',\
'wfn.out' u 1:23 w lines title 'l=2, n=4'
pause -1