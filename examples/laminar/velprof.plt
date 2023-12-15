set terminal pdf size 11cm,11cm
set output "U.pdf"

set xlabel "y" font ",16"
set ylabel "U" font ",16"
set tics font ",14"
set size square

set key top left font ",16"

u0=0.1
visc=0.166
L=64
F=8.10546875e-06
nx=64
dx=L/nx

set xrange [0:1]

plot \
40*x*(1-x)*u0 with lines lw 2.5 lc rgb "black" title "Analytic" ,\
"VELMEAN.dat" using (($1-0.5)*dx/L):($2/u0) every 1 with points ps 0.5 pt 5 lw 3 dt 1 lc rgb "forest-green" title "DUGKS"
