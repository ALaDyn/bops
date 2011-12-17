set terminal postscript enhanced colour
set output 'multiplot.ps'
set multiplot
#set data style lines
set size .35, .25
set origin 0.,.7
set yrange[0:2]
plot "ninc00.xy","ninc03.xy" w l #densita' ioni
set origin 0.35,.75
set yrange[0:2]
plot "nenc00.xy", "nenc03.xy" w l  #densita' elettroni
set origin 0.0,0.5
set autoscale
plot "exsi01.xy","exsi02.xy"   #campo elettrostatico
set origin 0.35,0.5
set autoscale
plot "ezsi01.xy"    #campo elettromagnetico pol. s
set origin 0.7,0.5
set autoscale
plot "bysi01.xy"    #B_y
set origin 0.0,0.25
set autoscale
plot "fuep02.xy"    
set origin 0.35,0.25
set autoscale
#plot "fuip02.xy"
set origin 0.7,0.25
set autoscale
plot "uinc00.xy","ubac00.xy"
set origin 0,0
plot "pxxe02.xy" with dots
set origin 0.35,0
#plot "pxxi02.xy" with dots
unset multiplot

