#!/bin/sh

#Calculate the Zb using a center wavenumber of 5.336

Proj ESA.xyz 1 125 20 mag.dat
Minc mag.dat -2875.06055/2875.06055/-3429.12256/3633.48022 2.8 2.8 mag_test.grd
Cwt mag_test.grd 5.336 4.0 10 spectra_cwt_5.336.dat 
Lsfit spectra_cwt_5.336.dat Zt 3 43 32 37 Zt_5.336_dat Zt_5.336_uncertainty.dat
Lsfit spectra_cwt_5.336.dat Z0 3 43 22 31 Z0_5.336_dat  Z0_5.336_uncertainty.dat
Depth Z0_5.336.dat Zt_5.336.dat Z0_5.336_uncertainty.dat Zt_5.336_uncertainty.dat Zb_5.336.dat Zb_5.336_uncertainty.dat
Proj Zb_5.336.dat -1 125 20 Z_5.336.dat

xyzname=Z_5.336.dat
grdname=Z_5.336.grd
cptname=Z_cwt_average.cpt
filename=Z_5.336.ps
toponame=topo.grd

R=100/150/-10/50


gmt surface $xyzname -G$grdname -R$R -I0.05/0.05
gmt surface topo.xyz -G$toponame -R$R -I0.05/0.05
gmt grdgradient $toponame -Ne0.7 -A90 -G$toponame.light
gmt grdimage -R$R -Ba20WSEN -JU14c $grdname -C$cptname -E300 -K -V -I$toponame.light >$filename
gmt psmask NaN.xyz -R$R -JU14c -I0.05/0.05 -Ggrey -K -O -V >> $filename
gmt psmask -C -K -O >>$filename
gmt pscoast -R$R  -JU14c -B10/10 -O -K -Dh -Ggrey -A20000 -I1 -W1/0.01i -N3/1p >>$filename
gmt psscale -D3i/-0.3i/4i/0.1ih -O -E  -I -C$cptname -Ba5/:"(km)": >>$filename

