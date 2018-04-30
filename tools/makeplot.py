from optparse import OptionParser
import config
import os
import fnmatch

parser = OptionParser(usage="%prog -f eventpath ")
parser.add_option("-f", "--evpath", type="string", dest="evpath", help="evpath")

(options, args) = parser.parse_args()

def writeplotfile(Origin,eventpath):
    
    p = os.path.join(eventpath,'metaplot.sh')
    fobj = open(p,'w')
    
    for i in os.listdir(eventpath):
        if fnmatch.fnmatch(i, '*.meta'):
            meta = i
    
    metaf = os.path.join(eventpath,meta)
    metap = os.path.join(eventpath,'metaplot.ps')
    part1='''
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

    psfile='''+metap+'''\n
    
    lat='''+Origin['lat']+'''\n
    lon='''+Origin['lon']+'''\n
    
    R=g
    J=E$lon/$lat/150/16

    pscoast -R$R -J$J  -W1/0 -Dc -V -S128/144/184 -G192/192/160 -K > $psfile

    awk '{print $6,$5}' '''+metaf+''' | psxy -R$R -J$J -St0.2 -W1p -Gred -O -K>> $psfile

    gmtset POLAR_CAP               = 90/90
    psbasemap  -R0/360/0/90 -JE0/90/16 -K -O -Bg30:"":/g15:"":wsne >> $psfile
    echo $lon $lat | psxy -R$Rc -J$J -Sa1 -W2p,0 -Gred -N -O >> $psfile

    gv $psfile &
    '''
    
    fobj.write(part1)
    fobj.close()
    
    cmd ="bash "+p
    os.system(cmd)
    

def init():
    C = config.Config(options.evpath)
    Origin = C.parseConfig('origin')
    
    return Origin

if __name__ == "__main__":
    
    ODict = init()
    writeplotfile(ODict,options.evpath)
    
    