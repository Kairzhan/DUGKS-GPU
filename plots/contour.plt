# Terminal type defines output format
# other options include 'svg', 'png', 'ps' etc.
#set terminal pdf size 18cm, 9cm
#set output "Contour.pdf"
set terminal png size 3149, 500 font "Arial, 16pt" #18cm, 9cm
set output "Contour.png"

unset key
set cblabel "U_x" 

set xlabel "X"
set ylabel "Y"

# Other color styles are available as well
#load 'jet.pal'
load 'moreland.pal'
#load 'spectral.pal'

set xrange [0:804]
set yrange [0:128]

#set title "Streamwise velocity contours"

# Plot data from the file 'c.dat', 1st and 3rd columns contain X, Y data,
# while 5th column contains required values, note, that 5th column is 
# scaled and such airthmetics will require braces and dollar sign ($)
# before column number.
#plot 'UVW.dat' using 3:1:($7) with image title "A"
plot 's2Y00005000.tec' skip 3 using 3:1:($7) with image title "A"
