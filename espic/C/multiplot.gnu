 
set multiplot;
#set data style lines
set size .4, .4; 
set origin 0.,.5;  set yrange[0:2];  plot "density1.data" w l, "density2.data" w l, "density3.data" w l, "density4.data" w l; 
set origin 0.0,0.0;  set autoscale; plot  "field1.data" w l, "field2.data" w l, "field3.data" w l,"field4.data" w l;  


set nomultiplot

