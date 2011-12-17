#!/usr/bin/gnuplot

print "Rudimentary Gnuplot script for visualizing"
print "particles from a sequence of ascii data"
print "files with the name 'pxxeNN.xy',"
print "where NN represents a 2-digit-number,"
print "starting from 00"
print "The script simply loops over all files and"
print "draws a 2D scatterplot. The particles are"
print "colored according to their charge."



# construct timestep counter
if (!exists("timestep")) timestep = 0
#if (!exists("timestep")) timestep = 1

print sprintf("=========== Timestep %2.2d ==========", timestep)

# filename constructor
filename(itime) = sprintf("pxxe%2.2d.xy", itime)
#filenameX(itime) = sprintf("exsi%2.2d.xy", itime)
#filenameZ(itime) = sprintf("ezsi%2.2d.xy", itime)


# define function to generate rgb color value form input ranges 0.0..1.0
rgb(r,g,b) = int(255*r)*65536 + int(255*g)*256 + int(255*b)

# color mapping function for charges: negative: blue, positive: red, neutral: grey
chargecolor(q) = q<0 ? rgb(0.,0.,1.) : ( q>0 ? rgb(1.,0.,0.) : rgb(.3,.3,.3) )

# set aspect ratio to 1:1:1
set size 1,1

actfilename = filename(timestep)
#actfilenameX = filenameX(timestep)
#actfilenameZ = filenameZ(timestep)

# set title
titletext = sprintf("Showing file %s", actfilename)
set title titletext

# 2D-color-plot: x, y, q
set xrange[0:500]
set yrange[-5:25]
#set zrange[-1:1]
#set style line 1 lw 3 lc rgbcolor "red"
set style line 2 lw 3 lc rgbcolor "blue"
set xtics("50" 50,"60" 60, "440" 440,"450" 450)
set ytics("-5" -5,"0" 0,"25" 25)
set grid
plot actfilename using 1:2 with dots ls 2 
#plot actfilename using 1:2 with points pt .7 ps .1 lc rgb variable notitle
#plot actfilenameX using 1:2 with lines ls 1 variable notitle
#replot actfilenameZ using 1:2 with lines ls 2 variable notitle

timestep = timestep + 1;

# check whether new file exists. if not: restart from the beginning
actfilename = filename(timestep)
#actfilenameX = filenameX(timestep)
#actfilenameZ = filenameZ(timestep)
fileexists = system(sprintf('if [ -e "%s" ]; then echo "1"; else echo "0"; fi', actfilename))
#fileexists = system(sprintf('if [ -e "%s" ]; then echo "1"; else echo "0"; fi', actfilenameX))
if (fileexists == 0) timestep = 0
#if (fileexists == 0) timestep = 1

pause 0.6
reread
