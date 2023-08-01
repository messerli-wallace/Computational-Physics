set grid
#set hidden3d
set xlabel 'x [cm]'
set ylabel 'y [cm]'
set zlabel 'potential [V]' offset -3,-3,0

#splot 'fort.1' u 1:2:3 with lines lt rgb 'black'
splot 'fort.1' u 1:2:3 with lines palette