#!/bin/bash

path=$PWD"/"$0
cpath=`echo $path | sed 's/plot-semb-lutz.sh//g'`
cd $cpath


#makecpt -Crainbow -T-0.0/1/0.1 -Z > cpt
rm -f .gmtdef*

gmtdefaults -D > .gmtdefaults4
gmtset  BASEMAP_TYPE = fancy
gmtset HEADER_FONT_SIZE        = 15p
gmtset HEADER_OFFSET           = 0.0
gmtset  LABEL_FONT_SIZE         = 12p
gmtset  ANNOT_FONT_SIZE    = 12c
gmtset LABEL_OFFSET    = 0.1
gmtset ANNOT_OFFSET_PRIMARY    = 0.05
gmtset PAPER_MEDIA     Custom_7ix6i

#maxsemb=`awk '{print $3}' sembmaxvalue.txt | sort -g | tail -1 | awk '{printf("%g", $1)}' `
maxcumul=`awk 'BEGIN {sum=0} {sum+=$3} END {printf("%g",sum)}' sembmaxvalue.txt `
maxcumulint=` echo $maxcumul | awk '{print 0.1*$1}' `
maxcumul50=` echo $maxcumul | awk '{print 0.5*$1}' `
maxcumulrange=` echo $maxcumul | awk '{print -0.1*$1"/"1.1*$1}' `
maxsemb=`awk 'BEGIN {m=-1} {if($3>m) {m=$3}} END {printf("%g",m)}' sembmaxvalue.txt `
maxsembint=` echo $maxsemb | awk '{print 0.1*$1}' `
maxsembrange=` echo $maxsemb | awk '{print -0.1*$1"/"1.1*$1}' `
maxcumulatick=`echo $maxcumul | awk '{printf("%.2f",$1)}' `
maxcumulftick=`echo $maxcumulatick | awk '{print $1/10}' `
maxsembatick=`echo $maxsemb | awk '{printf("%.3f",$1)}' `
maxsembftick=`echo $maxsembatick | awk '{print $1/10}' `

origtime=`grep source stations_0.cfg | awk '{print $5}' | awk -F'[:]' 'BEGIN {a["Jan"]=01;a["Feb"]=02;a["Mar"]=03;a["Apr"]=04;a["May"]=05;a["Jun"]=06;a["Jul"]=07;a["Aug"]=08;a["Sep"]=09;a["Oct"]=10;a["Nov"]=11;a["Dec"]=12} {print mktime(sprintf("%04d %02d %02d %02d %02d %02d",$3,a[$2],$1,$4,$5,$6))}' `
advance=200
step=`grep step stations_0.cfg | awk '{print $2}' `
lat=`grep source stations_0.cfg | awk '{print $2}' `
lon=`grep source stations_0.cfg | awk '{print $3}' `
depth=`grep source stations_0.cfg | awk '{print $4}' `

echo $maxcumul $maxsemb $maxsembrange $maxcumulrange

makecpt -Cseis -I -Z -T0/$maxcumul/$maxcumulint > cumul.cpt

c=1
for i in `ls *.ASC`; do

    name=`echo $i | awk -F '.' '{printf("semblance_%3.3i",$1)}'`

# name:
    psfile=$name.ps
    pngfile=$name.png
    grdfile=$name.grd
    echo $i $psfile

# Gebietsgrenzen:
    swlat=`awk '$2=="southwestlat:"{print $3}' $i`
    nelat=`awk '$2=="southwestlat:"{print $3+$5*($7-1)}' $i`
    swlon=`awk '$2=="southwestlon:"{print $3}' $i`
    nelon=`awk '$2=="southwestlon:"{print $3+$5*($7-1)}' $i`
    ilon=`awk '$2=="southwestlon:"{print $5}' $i`
    ilat=`awk '$2=="southwestlat:"{print $5}' $i`
    jlon=`awk 'BEGIN {print ('$swlon'-(-1)*'$nelon')/2}'`

# projektion und begrenzung:
    J=Q$jlon/6
    R=$swlon/$nelon/$swlat/$nelat
    I=$ilon/$ilat


# darstellung der xyz Werte:
#     psbasemap -O -K -R$R -J$J -B1::/1::WSne -X1 -Y2 >>$psfile
     psbasemap -K -P -R$R -J$J -B2::/2::WSne >$psfile


# dasselbe nocheinmal als grid:
    awk 'substr($1,1,1)!="#" {print $2,$1,$3}' $i | xyz2grd -R$R -I$I -G$grdfile
if [ $c == 1 ]; then
cp $grdfile cumul.grd
else
grdmath $grdfile cumul.grd ADD = cumul.grd
fi

timestr=`echo $c $start | awk '{print strftime("%d-%b-%Y_%H:%M:%S",'$origtime'+$1*'$step'-'$advance')}' `

grd2cpt -Cseis -I -Z $grdfile > cpt
    #grdview $grdfile -R$R -J$J -O -K -Ccpt -Ba1:"":/a1:"":WSen -Qs -Y5 >>$psfile
    grdview $grdfile -R$R -J$J -O -K -Ccpt -Ba2:"":/a2:"":WSen -Qs >>$psfile
#  grdcontour $grdfile -R$R -J$J -O -K -A0.1 -G2 >>$psfile
    pscoast -R$R -J$J -A10000 -W1/0 -N1 -Di -O -K >>$psfile
    awk '{print $1,$2}' /home/ehlert/.seiscomp3/platebound.gmt | psxy -R$R -J$J -W1p/0 -A -M  -K -O >> $psfile
# Maximum:
    grdinfo -M $grdfile | awk '$11=="z_max:" {print $16,$19}' |\
    psxy -R$R -J$J -O -K -Sc0.5 -W3p,0 >>$psfile
# plot semblance time trace
    awk 'NR<='$c' {print $2, $3}' sembmaxvalue.txt |\
    psxy -JX6/3 -R0/700/$maxsembrange -W3/255/0/0 \
         -Ba120f30:"Time [s]":/a${maxsembatick}f${maxsembftick}:"Semblance":WSne -Y8. -O -K >> $psfile
    pstext -JX -R0/1/0/1 -O -K << END >> $psfile
0.7 0.8 10 0 0 5 absolute
END

# Skala:
    #psscale -Ccpt -D12/2/4/0.5 -B0.1:"semblance":/:"":WSen -K -P -O >>$psfile
    grdview cumul.grd -R$R -J$J -O -K -X7.5 -Y-8. -Ccumul.cpt -Ba2:"":/a2:"":WSen -Qs >>$psfile
    grdcontour cumul.grd -R$R -J$J -O -K -C${maxcumulint} >>$psfile
    grdcontour cumul.grd -R$R -J$J -O -K -C${maxcumul50} -W5/0 >>$psfile
    pscoast -R$R -J$J -A10000 -W1/0 -N1 -Di -O -K >>$psfile
    awk '{print $1,$2}' /home/ehlert/.seiscomp3/platebound.gmt | psxy -R$R -J$J -W1p/0 -A -M  -K -O >> $psfile
# Maximum:
    grdinfo -M cumul.grd | awk '$11=="z_max:" {print $16,$19}' |\
    psxy -R$R -J$J -O -K -Sc0.5 -W3p,0 >>$psfile
# plot semblance time trace
    awk 'BEGIN {sum=0} {if(NR<='$c') {sum+=$3; print $2, sum}}' sembmaxvalue.txt |\
    psxy -JX6/3 -R0/700/${maxcumulrange} -W3/255/0/0 \
         -Ba120f30:"Time [s]":/a${maxcumulatick}f${maxcumulftick}WSne -Y8. -O -K >> $psfile
    pstext -JX -R0/1/0/1 -O -K -N << END >> $psfile
0.7 0.8 10 0 0 5 cumulative
-0.5 1.2 12 0 0 5 $timestr
END


# dummy closing file
    pstext -JX1 -R0/1/0/1 -N -O -P <<END>> $psfile
END
convert -trim -density 75x75 $psfile $pngfile

c=`echo $c | awk '{print $1+1}' `

#    rm -f $grdfile
done



ag="-anim.gif"
animation=$1$ag 

#convert -delay 20 semblance_???.png test.gif
convert -delay 20 semblance_???.png $animation




scp $animation web_tsunami@141.89.111.29:/Webpage2010/tl_files/arbeitsgruppen/ag_array_seismologie/ruptrack/img/new

echo $animation web_tsunami@141.89.111.29:/Webpage2010/tl_files/arbeitsgruppen/ag_array_seismologie/ruptrack/img/new


./plot-summary.sh $1
#gs semblance_???.ps



