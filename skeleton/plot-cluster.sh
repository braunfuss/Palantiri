#!/bin/bash

rm -f .gmt*
gmtset BASEMAP_TYPE = plain
gmtset HEADER_OFFSET           = -0.1c
gmtset ANNOT_OFFSET_PRIMARY    = 0.05
gmtset LABEL_OFFSET            = 0.05c
gmtset ANNOT_FONT_SIZE = 14p
gmtset LABEL_FONT_SIZE = 10p
gmtset HEADER_FONT_SIZE = 18p
gmtset PLOT_DEGREE_FORMAT = +DF
gmtset PAPER_MEDIA         = a4


input="plot_cluster.dat"

#Variablen
PS=cluster.ps

#haiti
WESN=-180/180/-85/85
SCALING=0.05c

#-B1p/100
gmtset BASEMAP_TYPE = fancy

psbasemap -R$WESN -Jm$SCALING -Lf287/18/48.3/100 -B45:"":/45:""::."Cluster Overview":WSEn -V -K > $PS
pscoast -R -Jm -B -W1.0p -G237/221/168 -S169/227/249 -Dc -N1 -Ia -O -V -K >> $PS

stafile="eq.xy"
staname="stanme.xyn"
rm $staname
rm $stafile

while read line; do

cluster=`echo $line |  awk -F' ' '{print $1}'`
lat=`echo $line |  awk -F' ' '{print $3}'`
lon=`echo $line |  awk -F' ' '{print $4}'`
#mag=`echo $line | awk -F, '{print $2}' | sed s/'"'/''/g| awk -F' ' '{print $1}'`

echo $lon","$lat","$cluster
echo $lon","$lat,","$cluster  >>$stafile

done < $input


#pstext $staname -R$WESN -P -Jm$SCALING -G255/255/255 -S4/0/0/0 -O -K  >> $PS
awk '{print $1,$2}' /home/sysop/.seiscomp3/platebound.gmt | psxy -R$WESN -Jm$SCALING -W2p/100/100/150 -A -M  -K -O >> $PS


if    [ ${cluster} == "1" ]; then   
	XY=255/0/0
	awk -F, '{print $1, $2} ' $stafile | psxy -R$WESN -Jm$SCALING -Sa0.75 -G255/0/0 -G${XY} -W2 -K -O >> $PS
elif [ ${cluster} == "2" ]; then
	XY=0/255/0
	awk -F, '{print $1, $2} ' $stafile | psxy -R$WESN -Jm$SCALING -Sa0.75 -G255/0/0 -G${XY} -W2 -K -O >> $PS
elif [ ${cluster} == "3" ];then      
	XY=0/0/255
	awk -F, '{print $1, $2} ' $stafile | psxy -R$WESN -Jm$SCALING -Sa0.75 -G255/0/0 -G${XY} -W2 -K -O >> $PS
elif [ ${cluster} == "4" ];then 
	XY=255/255/0
	awk -F, '{print $1, $2} ' $stafile | psxy -R$WESN -Jm$SCALING -Sa0.75 -G255/0/0 -G${XY} -W2 -K -O >> $PS
elif [ ${cluster} == "5" ];then
	XY=0/255/255
	awk -F, '{print $1, $2} ' $stafile | psxy -R$WESN -Jm$SCALING -Sa0.75 -G255/0/0 -G${XY} -W2 -K -O >> $PS
elif [ ${cluster} == "6" ];then
	XY=255/0/255
	awk -F, '{print $1, $2} ' $stafile | psxy -R$WESN -Jm$SCALING -Sa0.75 -G255/0/0 -G${XY} -W2 -K -O >> $PS
elif [ ${cluster} == "7" ];then
	XY=238/48/167
	awk -F, '{print $1, $2} ' $stafile | psxy -R$WESN -Jm$SCALING -Sa0.75 -G255/0/0 -G${XY} -W2 -K -O >> $PS
elif [ ${cluster} == "8" ];then
	XY=106/90/205
	awk -F, '{print $1, $2} ' $stafile | psxy -R$WESN -Jm$SCALING -Sa0.75 -G255/0/0 -G${XY} -W2 -K -O >> $PS
elif [ ${cluster} == "9" ];then
	XY=255/165/0
	awk -F, '{print $1, $2} ' $stafile | psxy -R$WESN -Jm$SCALING -Sa0.75 -G255/0/0 -G${XY} -W2 -K -O >> $PS
elif [ ${cluster} == "10" ];then
	XY=205/92/92
	awk -F, '{print $1, $2} ' $stafile | psxy -R$WESN -Jm$SCALING -Sa0.75 -G255/0/0 -G${XY} -W2 -K -O >> $PS
elif [ ${cluster} == "11" ];then
	XY=208/32/144
	awk -F, '{print $1, $2} ' $stafile | psxy -R$WESN -Jm$SCALING -Sa0.75 -G255/0/0 -G${XY} -W2 -K -O >> $PS
elif [ ${cluster} == "12" ];then
	XY=255/255/255
	awk -F, '{print $1, $2} ' $stafile | psxy -R$WESN -Jm$SCALING -Sa0.75 -G255/0/0 -G${XY} -W2 -K -O >> $PS	
else
	XY=30/100/20
	awk -F, '{print $1, $2} ' $stafile | psxy -R$WESN -Jm$SCALING -Sa0.75 -G255/0/0 -G${XY} -W2 -K -O >> $PS
fi


#awk 'substr($1,1,1)!="#"&&substr($1,1,1)!="%" {print $4,$3}' ./$i | psxy -R$R -J$J -St0.2 -G${XY}  -W1/0/0/0  -K -O >> $psfile

#awk -F, '{print $1, $2} ' $stafile | psxy -R$WESN -Jm$SCALING -Sa0.75 -G255/0/0 -G${XY} -W2 -K -O >> $PS

#psbasemap -R$WESN -Jm$SCALING -Lf286/17.2/48.3/100 -B45:"":/45:""::."Haiti 12.01.2010":WSEn -V -K -O >> $PS
lon1=-72.55
lat1=18.37
#echo $lon1 $lat1 | psxy -R$WESNc -Jm$SCALING -Sa0.75 -W2p,0 -Gred -N -O >> $PS


gv $PS




