#!/bin/sh

#Calculate the Zb using a window size of 302.4km

Proj ESA.xyz 1 125 20 mag.dat
Minc mag.dat -2875.06055/2875.06055/-3429.12256/3633.48022 2.8 2.8 mag_test.grd
Wft mag_test.grd 116 36 0.006 0.15 spectra_fft_300.dat 
Lsfit spectra_fft_300.dat Zt 3 27 7 14 Zt_300.dat Zt_300_uncertainty.dat
Lsfit spectra_fft_300.dat Z0 3 27 1 6 Z0_300.dat  Z0_300_uncertainty.dat
Depth Z0_300.dat Zt_300.dat Z0_300_uncertainty.dat Zt_300_uncertainty.dat Zb_300.dat Zb_300_uncertainty.dat
Proj Zb_300.dat -1 125 20 Z_300.dat

xyzname=Z_300.dat
grdname=Z_300.grd
cptname=Z_fft_average.cpt
filename=Z_300.ps
toponame=topo.grd

R=100/150/-10/50


gmt surface $xyzname -G$grdname -R100/150/-10/50 -I0.05/0.05
gmt surface topo.xyz -G$toponame -R100/150/-10/50 -I0.05/0.05
gmt grdgradient $toponame -Ne0.7 -A90 -G$toponame.light
gmt grdimage -R$R -Ba10WSEN -JU14c $grdname -C$cptname -E300 -K -V -I$toponame.light >$filename
gmt psmask NaN.xyz -R$R -JU14c -I0.05/0.05 -Ggrey -K -O -V >> $filename
gmt psmask -C -K -O >>$filename
gmt pscoast -R$R  -JU14c -B10/10 -O -K -Dh -Ggrey -A30000 -I1 -W1/0.01i -N3/1p >>$filename
gmt psscale -D3i/-0.3i/4i/0.1ih -O -E  -I -C$cptname -Ba5/:"(km)": >>$filename

