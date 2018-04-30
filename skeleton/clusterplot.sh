#!/bin/bash

rm -f .gmt*
gmtset BASEMAP_TYPE = plain
gmtset HEADER_OFFSET           = -0.1c
gmtset ANNOT_OFFSET_PRIMARY    = 0.05
gmtset LABEL_OFFSET            = 0.05c
gmtset ANNOT_FONT_SIZE = 12p
gmtset LABEL_FONT_SIZE = 12p
gmtset HEADER_FONT_SIZE = 18p
gmtset PLOT_DEGREE_FORMAT = +DF
gmtset PAPER_MEDIA         = a4

ofile="event.orig"
psfile=arraycluster.ps
conffile="event.stations"
centroidfile="event.centroid"

R=g
lon=`awk -F','  '{print $2}' $ofile`
lat=`awk -F','  '{print $1}' $ofile`
J=E$lon/$lat/180/18
echo $lon,$lat
echo $J

grdmath -Rg -I1 $lon $lat SDIST 111.13 MUL = dist.nc
pscoast -R$R -J$J  -W1/0 -Dc -V -S128/144/184 -G192/192/160 -U  -K > $psfile

awk '{print $3,$2}' $centroidfile | psxy -R$R -J$J -St0.2 -W1p -Gred -O -K >> $psfile

while read counter clat clon; do
	(echo $lon $lat; echo $clon $clat) | psxy -R$R -J$J -O -K -Wthickest/red >> $psfile
done < $centroidfile


while read name cla clo rank; do

if [ ${rank:-0} == "0" ];then
	XY=208/32/144
elif    [ ${rank:-0} == "1" ]; then   
	XY=255/0/0
elif [ ${rank:-0} == "2" ]; then
	XY=0/255/0
elif [ ${rank:-0} == "3" ];then      
	XY=0/0/255
elif [ ${rank:-0} == "4" ];then 
	XY=255/255/0
elif [ ${rank:-0} == "5" ];then
	XY=0/255/255
elif [ ${rank:-0} == "6" ];then
	XY=255/0/255
elif [ ${rank:-0} == "7" ];then
	XY=238/0/167
elif [ ${rank:-0} == "8" ];then
	XY=106/90/205
elif [ ${rank:-0} == "9" ];then
	XY=255/165/0
elif [ ${rank:-0} == "10" ];then
	XY=205/92/92
elif [ ${rank:-0} == "11" ];then
	XY=255/255/255	
else
	XY=30/100/20
fi

echo $clo $cla | psxy -R$R -J$J -Sd0.3 -G${XY}  -W1/0/0/0  -K -O >> $psfile
done < $conffile



while read rank cla clo; do

if [ ${rank:-0} == "0" ];then
	XY=208/32/144
elif    [ ${rank:-0} == "1" ]; then   
	XY=255/0/0
elif [ ${rank:-0} == "2" ]; then
	XY=0/255/0
elif [ ${rank:-0} == "3" ];then      
	XY=0/0/255
elif [ ${rank:-0} == "4" ];then 
	XY=255/255/0
elif [ ${rank:-0} == "5" ];then
	XY=0/255/255
elif [ ${rank:-0} == "6" ];then
	XY=255/0/255
elif [ ${rank:-0} == "7" ];then
	XY=238/48/167
elif [ ${rank:-0} == "8" ];then
	XY=106/90/205
elif [ ${rank:-0} == "9" ];then
	XY=255/165/0
elif [ ${rank:-0} == "10" ];then
	XY=205/92/92
elif [ ${rank:-0} == "11" ];then
	XY=255/255/255	
else
	XY=30/100/20
fi

echo $clo $cla | psxy -R$R -J$J -Sa0.5 -G${XY}  -W1/0/0/0  -K -O >> $psfile
done < ${centroidfile}



echo $lon $lat | psxy -R -J -O -K -Sa0.2i -Gyellow -Wthin >> $psfile

grdcontour dist.nc -A3000+v+ukm+kwhite -Glz-/z+ -S8 -C1000 -O  -J -Wathin,white -Wcthinnest,white,- >> $psfile

rm -f dist.nc

gv $psfile &

cp $psfile ../..
