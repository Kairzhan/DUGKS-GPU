set terminal pdf size 11cm,11cm 
set output "Urms.pdf"

set xlabel "y/2H" font ",16"
set ylabel "u^+_{rms}, v^+_{rms}, w^+_{rms}" font ",16"
set tics font ",14"
set ytics 0.25
set xrange [0:0.5]
set yrange [0:3]
set size square
set key top right font ",16"

#set tmargin 0
set bmargin 3

Re_tau=180
visc=0.002
TWO_H=128.0
utau=2*Re_tau*visc/TWO_H
yplus=visc/utau
nx=384

dx=TWO_H/nx

H=180

plot \
"JAXA1.dat" using ($2/H/2):(sqrt($4)) with lines lw 2 lc rgb "red" title "JAXA CHAN180", \
"JAXA1.dat" using ($2/H/2):(sqrt($5)) with lines lw 2 lc rgb "red" notitle, \
"JAXA2.dat" using ($2/H/2):(sqrt($3)) with lines lw 2 lc rgb "red" notitle, \
"VELMEAN.dat" using (($1-0.5)*dx/TWO_H):(sqrt($4)/utau) every 4 with points ps 0.6 pt 5 lw 3 dt 6 lc rgb "forest-green" title "DUGKS,  Re_{τ}=180, u^+_{rms}", \
"VELMEAN.dat" using (($1-0.5)*dx/TWO_H):(sqrt($5)/utau) every 4 with points ps 0.6 pt 5 lw 3 dt 7 lc rgb "forest-green" title "v^+_{rms}", \
"VELMEAN.dat" using (($1-0.5)*dx/TWO_H):(sqrt($6)/utau) every 4 with points ps 0.6 pt 5 lw 3 dt 5 lc rgb "forest-green" title "w^+_{rms}"
