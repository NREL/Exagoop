set terminal pngcairo dashed size 800,500
set output "fig1.png"
set termoption font "Helvetica,18"
set key spacing 1.2
set lmargin 9
set rmargin 5
set bmargin 3.5
set xlabel "X Coordinate (m)"
set ylabel "Y Coordinate (m)"

set yrange [0.45:0.54]
#set xtics (0,2e-5,5e-5,7e-5,0.0001)
#set format x "%5.2e"

plot 'rawdeflection.dat' u 1:2 w p lw 1 ps 2 pt 6 dt 1 title "Simulated raw",\
'meandeflection.dat' u 1:2 w l lw 3 dt 1 title "Simulated mean",\
 'meandeflection.dat' u 1:(0.5066-0.0377*((0.7-$1)**4+4*0.7**3*$1-0.7**4)) w l lw 4 dt 1 title "Euler Beam Theory"
