set terminal pdf size 11cm,11cm
set output "Ulog.pdf"

set xlabel "y^+" font ",16"
set ylabel "U^+" font ",16"
set tics font ",14"
set xrange [0.3:180]
set yrange [0:22]
set size square
set key top left font ",16"
set logscale x
set ytics nomirror
#set y2tics 2
#set my2tics 2

#set bmargin 3

k=0.40
B=5.6

Re_tau=180
visc=0.002
TWO_H=128.0
utau=2*Re_tau*visc/TWO_H
yplus=visc/utau
nx=384

dx=TWO_H/nx

print dx, yplus
print utau

plot \
x with lines dashtype 3 lw 1 lc rgb "black" notitle ,\
1/k*log(x)+B with lines dashtype 3 lw 1 lc rgb "black" notitle ,\
"JAXA1.dat" u 2:3 with lines lw 3 lc rgb "red", \
"VELMEAN.dat" using (($1-0.5)*dx/yplus):($2/utau) every 2 with points ps 0.5 pt 5 lw 3 dt 1 lc rgb "forest-green" title "DUGKS, Re_{τ}=180"
