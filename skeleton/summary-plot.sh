#!/bin/bash

#set verbose
export TZ=UTC

advance=200

conffile=stations_0.cfg

duration=`awk '$1=="%duration:" {print $2}' $conffile`				 #duration
step=`awk '$1=="%step:" {print $2}' $conffile`   					 # time step
maxstep=`awk '$1=="%step:" {print int('$duration'/$2)}' $conffile`            # Anzahl der Zeitschritte
gridspacing=`awk '$1=="%gridspacing:" {print $2}' $conffile`                       # time step

# get start time for SH:
#t1=`awk '$1=="%source:" {print $5}' $conffile | sed -e's/_/ /g' | sed -e's/:/ /g' | sed -e's/-/ /g' |  awk '{printf("%04i %04i %02i %02i %02i %02i", $3,$2,$1,$4,$5,$6)}'`
#t=`awk  -v t="$t1" 'BEGIN {print mktime(t)}'`
#start=`awk 'BEGIN {print strftime("%d-%b-%Y_%H:%M:%S",'$t')}'`

stime=`awk '$1=="%source:"{print $5}' $conffile` 

# prepare migregion.dat (5x5 degree around hypocentre, gridinterval: 0.2 degree):
lat=`awk '$1=="%source:" {print $2}' $conffile`
lon=`awk '$1=="%source:" {print $3}' $conffile`
swlat=`echo $lat | awk '{if(($1<85)&&($1>-85)){print $1-5.} if($1<-85) {print -80.} if($1>85) {print 80.}}'`
swlon=`echo $lon | awk '{print $1-5.}'`
nelat=`echo $swlat | awk '{print $1+10.}'`
nelon=`echo $swlon | awk '{print $1+10.}'`

# grid interval:
I=`awk '$1=="%gridspacing:" {print $2}' $conffile`   					 # gridspacing
#I=0.2

directory=`pwd | awk -F '/' '{print $NF}'`

# check if computations are ready on nodes an start plotting
# prerequisites:
makecpt -Crainbow -T0/1.1/0.1 -Z > col.cpt
rm -f .gmt*
gmtset BASEMAP_TYPE = plain
gmtset HEADER_OFFSET           = -0.3c
gmtset ANNOT_OFFSET_PRIMARY    = 0.05
gmtset LABEL_OFFSET            = 0.05c
gmtset ANNOT_FONT_SIZE = 8p
gmtset LABEL_FONT_SIZE = 8p
gmtset PLOT_DEGREE_FORMAT +DF
gmtset PAPER_MEDIA         = a4
R=$swlon/$nelon/$swlat/$nelat
J=Q$lon/15
echo $J
echo $R
#echo $start
echo $stime

psfile=`pwd | awk -F '/' '{print $NF"-summary-test.ps"}'`

rm -f .gmt*
gmtset BASEMAP_TYPE = plain
gmtset HEADER_OFFSET           = -0.1c
gmtset ANNOT_OFFSET_PRIMARY    = 0.05
gmtset LABEL_OFFSET            = 0.05c
gmtset ANNOT_FONT_SIZE = 8p
gmtset LABEL_FONT_SIZE = 8p
gmtset HEADER_FONT_SIZE = 18p
gmtset PLOT_DEGREE_FORMAT = +DF
gmtset PAPER_MEDIA         = a4

r=`awk '{print $1}' sembmax.txt | sort -g | tail -1`
l=`awk '{print $4*1.1}' sembmax.txt | sort -g | tail -1`
Rs=0/$r/0/$l

pstext -JX1  -R0/1/0/1 -N -P -K <<END > $psfile
#7.75  26 12 0 0 MC `pwd`/$psfile
#7.75 25.5 12 0 0 MC program: `which $0`
#7.75 25. 12 0 0 MC Zeit: `date +%D_%H:%M:%S` on `hostname`
4.25 -1 12 0 0 BC  Time [s] from start at $stime
#4.25 -1 12 0 0 BC  Time [s] from start at $start
END

#  semblance maximum until that time:
            Js=X8.5/4
            dstart=`awk 'substr($1,1,1)!="#"{print $2}' duration.txt | head -1`	
            end=`awk 'substr($1,1,1)!="#"{print $2}' duration.txt | tail -1`		

            istart=`awk 'substr($1,1,1)!="#"{print $1}' duration.txt | head -1`
            iend=`awk 'substr($1,1,1)!="#"{print $1}' duration.txt | tail -1`
            duration=`awk 'substr($1,1,1)=="#"&&$3=="duration" {printf("%d", int($5+0.5))}' duration.txt`
            l=`awk '$1>='$dstart'&&$1<='$end' {printf("%d \n", int($8+0.5))}' sembmax.txt | sort  -g  | tail -1`

             echo  $dstart 0  $dstart 10  $end 10  $end 0 $duration $l | psxy -R$Rs -J$Js -W1p,0 -G255/193/193 -L -O -K -N -X0 >> $psfile		


            awk '{print $1, $4}' sembmax.txt |\
            psxy -R$Rs -J$Js -W1p,0  -Ba100f10:"":/a0.1f0.001:"Semblance"::."Semblance Maximum vs. Time":WSn -O -K -N >> $psfile

           awk 'substr($1,1,1)!="#" {print $2, $3}' duration.txt |   psxy -R$Rs -J$Js -W2p,red  -O -K -N >> $psfile						
   
# 3D- Entfernung:
Rs=`awk '{print $1,$4*1.1}' sembmax.txt  | minmax -I0.001`
awk '{print $1,$8}' sembmax.txt |\
psxy $Rs -J$Js -Sx0.2 -W1p,0 -K -O -Ba50f10:"":/a100f10:"Distance  [km]"::."":E  >> $psfile

# Richtung (Azimuth):
sr=`awk '{print $1}' sembmax.txt | sort -g  | tail -1`
sl=`awk '{print $1}' sembmax.txt | sort -rg | tail -1`
Rr=$sl/$sr/0/360

awk '{print $1,$7}' sembmax.txt |\
psxy -R$Rr -J$Js -Sc0.1 -W1p,0 -K -O >> $psfile
psbasemap -R$Rr -J$Js -Ba50f10:"":/a90f10:"Azimuth  [degree]"::."":E -K -O -X2 >> $psfile
 
# Polaransicht
R=g
lat=`awk '$1=="%source:" {print $2}' ./stations_0.cfg `
lon=`awk '$1=="%source:"{print $3}' ./stations_0.cfg `
J=E$lon/$lat/120/5.5

pscoast -O -K -R$R -J$J  -W1/0 -Dc -V -S128/144/184 -G192/192/160 -X10.5 -Y-1 >> $psfile

for i in *.dat ; do

rank=`awk 'NR==1{print $1}' $i ` 

if    [ ${rank:-0} == "1" ]; then   
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
	XY=208/32/144
elif [ ${rank:-0} == "12" ];then
	XY=255/255/255	
else
	XY=30/100/20
fi

awk 'substr($1,1,1)!="#"&&substr($1,1,1)!="%" {print $4,$3}' ./$i |\

psxy -R$R -J$J -St0.2 -G${XY}  -W1/0/0/0  -K -O >> $psfile

gmtset POLAR_CAP               = 90/90
psbasemap  -R0/360/0/90 -JE0/90/5.5 -K -O -Bg30:"":/g15:"":wsne >> $psfile

done


# area:
gmtset BASEMAP_TYPE = fancy
idx=`awk 'substr($1,1,1)!="#" {printf("%3.3d \n", $1)}' duration.txt`
name=`ls *000.grd | tail -1 | awk -F '_' '{print $1}'`
R=`grdinfo -I0.01 *_000.grd`
J=Q$lon/15
contour=0.9
makecpt -Cno_green -T0/$duration/5 -Z -N > col.cpt

psbasemap $R -J$J -Ba0.2g${gridspacing}:"":/a0.2g${gridspacing}:""::."Semblance - Time Map":WSen -O -K -Y7 -X-12.5 >> $psfile
awk '{print $1,$2}' platebound.gmt | psxy $R -J$J -W1p/0 -A -M  -K -O >> $psfile
sw=0

for i in $idx; do
      grdfile=${name}_${i}.grd
      ii=`echo $i | awk '{printf("%i",$1)}'`
      zeit=`awk 'substr($1,1,1)!="#"&&$1=='$ii' {print int($2-'$dstart')}' duration.txt`

# normalisation:
      smax=`grdinfo -M $grdfile | awk 'FNR==8{ print $12}'`
      grdmath $grdfile $smax DIV = tmp.grd
# create contousr:
      contfile=cont-$i
      grdcontour tmp.grd -J$J $R -C$contour -D$contfile > tmp
      contfiles=`ls $contfile"_$contour"_*.xyz`
#       echo $i $ii $contfile $zeit

# Darstellung Zeitschritte:
#   Farbe:
      W=`awk 'substr($1,1,1)!="#" {printf("%10.5f %3i %3i %3i \n", $1,$2,$3,$4)}' col.cpt |   awk '$1<='$zeit' {print ""$2"/"$3"/"$4}' | tail -1`
      for j in $contfiles; do
            awk '{print $1,$2}' $j |\
            psxy $R -J$J -W2p,$W -N -K -O -M >> $psfile
      done

# get area:
      if [ $sw -eq 0 ]; then
            sw=1
            grdmath tmp.grd $contour GE = area.grd
      else
            grdmath tmp.grd $contour GE area.grd ADD $contour GE = temp.grd
            mv temp.grd area.grd
      fi
      cp area.grd area-$i.grd


done

awk 'substr($1,1,1)!="#" {print $4,$5,$1 }' duration.txt |\
psxy $R -J$J -W2p,0 -Sx0.5 -K -O >> $psfile
flaeche=`ls area-*.grd | tail -1`
a0=`grdvolume -C0.9999 -Sk $flaeche | awk '{printf("%10i",$2+0.5)}'  | awk '{print $1}'`

pstext -JX1  -R0/1/0/1 -N -P -O -K  <<END >> $psfile
0.2 0.2 18 0 1 BL duration: $duration sec
0.2 1.8 18 0 1 BL rupture length: $l km
14.9 0.1 10 0 0 BR \251 University of Potsdam, `date +%Y-%b-%d_%H:%M:%S`
END

#pstext -JX1  -R0/1/0/1 -N -P -O -K  <<END >> $psfile
#0.2 0.2 18 0 1 BL duration: $duration sec
#0.2 1.0 18 0 1 BL rupture area: $a0 km\262
#0.2 1.8 18 0 1 BL rupture length: $l km
#14.9 0.1 10 0 0 BR \251 University of Potsdam, `date +%Y-%b-%d_%H:%M:%S`
#END

pscoast -J$J $R  -O -Dh -N1 -W1/0 -K -P -Ia>> $psfile
echo $lon $lat |\
psxy -R$Rc -J$J -Sa1 -W2p,0 -Gred -N -K -O >> $psfile

B=`echo $duration | awk '{printf("a%if%i",$1/5,$1/20)}'`
psscale -Ccol.cpt -D11.8/2.3/6/0.5h -B$B:"rupture time [s]":/:"": -K -P -O >> $psfile

#closing file
pstext -JX1  -R0/1/0/1 -N -O <<END>> $psfile
END

endname=$1_summary.jpg
convert $psfile $endname
echo $endname 

#rm -f $grdfile

exit
