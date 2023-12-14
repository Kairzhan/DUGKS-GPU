set terminal pdf size 11cm,11cm
set output "Ustress.pdf"

set xlabel "y/2H" font ",16"
set ylabel "-<u^'_x u^'_y>/u^*^2" font ",16"
set xrange [0:0.5]
set yrange [0:1]
set size square
set key top right font ",16"
set tics font ",14"

Re_tau=180
visc=0.002
TWO_H=128.0
utau=2*Re_tau*visc/TWO_H
yplus=visc/utau
nx=384

dx=TWO_H/nx

plot \
1-2*x with lines lw 2 lc rgb "black" notitle,\
"JAXA2.dat" using ($2/360):4 with lines lw 2 lc rgb "red" title "JAXA, CHAN180", \
"VELMEAN.dat"   using (($1-0.5)*dx/TWO_H):(-$3/utau**2) every 4 with points ps 0.6 pt 5 lw 3 lc rgb "forest-green" title  "DUGKS, Re_{τ}=180"
