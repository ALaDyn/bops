 
set multiplot;
set data style lines
set size .35, .25; 
set origin 0.,.75;  set yrange[0:2];  plot "density1.dat"; 
set origin 0.35,.75;  set yrange[0:2];  plot "density2.dat"; 
set origin 0.7,.75;  set yrange[0:2];  plot "density3.dat"; 
set origin 0.0,0.5;  set autoscale; plot "field1.dat";  
set origin 0.35,0.5;  set autoscale; plot "field2.dat";  
set origin 0.7,0.5;  set autoscale; plot "field3.dat";  
set origin 0.0,0.25;  set autoscale; plot "laser1.dat";  
set origin 0.35,0.25;  set autoscale; plot "laser2.dat";  
set origin 0.7,0.25;  set autoscale; plot "laser3.dat";  
set origin 0,0; plot "phase1.dat" with dots;
set origin 0.35,0; plot "phase2.dat" with dots;
set origin 0.7,0; plot "phase3.dat" with dots;

set nomultiplot

