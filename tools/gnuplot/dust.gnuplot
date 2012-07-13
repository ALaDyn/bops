#!/usr/bin/gnuplot

print "Rudimentary Gnuplot script for visualizing"
print "particles from a sequence offlat ascii data"
print "files with the name 'particles_xxxxxx.dat',"
print "where xxxxxx represents a six-digit-number,"
print "starting from 000000."
print "The script simply loops over all files and"
print "draws a 3D scatterplot. The particles are"
print "colored according to their charge."


# columns in outputfiles:
# /-----------------------------------------------\
# |  1|2|3      |4|5|6|7 |8 |9 |10|11|12 |13|14|15|
# |cpu|p|pelabel|x|y|z|vx|vy|vz|q |m |pot|ex|ey|ez|
# \-----------------------------------------------/

# construct timestep counter
if (!exists("timestep")) timestep = 0

print sprintf("=========== Timestep %6.6d ==========", timestep)

# filename constructor
filename(itime) = sprintf("particles_%6.6d.dat", itime)

# define function to generate rgb color value form input ranges 0.0..1.0
rgb(r,g,b) = int(255*r)*65536 + int(255*g)*256 + int(255*b)

# color mapping function for charges: negative: blue, positive: red, neutral: grey
chargecolor(q) = q<0 ? rgb(0.,0.,1.) : ( q>0 ? rgb(1.,0.,0.) : rgb(.3,.3,.3) )

# set aspect ratio to 1:1:1
set size 1,1

actfilename = filename(timestep)

# set title
titletext = sprintf("Showing file %s", actfilename)
set title titletext

# 3D-color-plot: x, y, z, q
set xrange[-1.:3.5]
set yrange[-.2:1.2]
set zrange[-.3:.3]
splot actfilename using 4:5:6:(chargecolor($10)) with points pt 7 ps .5 lc rgb variable notitle

timestep = timestep + 5;

# check wether new file exists. if not: restart from the beginning
actfilename = filename(timestep)
fileexists = system(sprintf('if [ -e "%s" ]; then echo "1"; else echo "0"; fi', actfilename))
if (fileexists == 0) timestep = 0

pause 0.3
reread
