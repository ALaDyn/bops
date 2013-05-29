 
set multiplot;
set size .35, .25; 
set origin 0.,.75;  set yrange[0:2];  plot "density1.data" w l; 
set origin 0.35,.75;  set yrange[0:2];  plot "density2.data" w l; 
set origin 0.7,.75;  set yrange[0:2];  plot "density3.data" w l; 
set origin 0.0,0.5;  set autoscale; plot "field1.data" w l;  
set origin 0.35,0.5;  set autoscale; plot "field2.data" w l;  
set origin 0.7,0.5;  set autoscale; plot "field3.data" w l;  
set origin 0.0,0.25;  set autoscale; plot "laser1.data" w l;  
set origin 0.35,0.25;  set autoscale; plot "laser2.data" w l;  
set origin 0.7,0.25;  set autoscale; plot "laser3.data" w l;  
set origin 0,0; plot "phase1.data" with dots;
set origin 0.35,0; plot "phase2.data" with dots;
set origin 0.7,0; plot "phase3.data" with dots;

set nomultiplot

