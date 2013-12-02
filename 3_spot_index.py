#! /usr/bin/env python

## based off ts2_mod.py - modified to use Cataspot and new global parameters file.
# Mat Iadanza 2013-09-16

# predict spots in all images (including non major plain images) based on the unit cell vectors produced by findcell.py and refined using findparallel.py
# outputs batch#.csv that is used by revcoords_tf.py to produce an illustration of the spots.  Uses Brent's new (more accurate) matrix method
#
# Matt Iadanza 2013-07-08 

#### Imports:

import numpy
import math
import os
import json
import datetime

####### user input variables

lauezonethresh = 0.25

########### Variables

params = json.load(open('params.json')) 
imgsize = params["imgsize"]      
boxsize =  params["circrad"]
imgmaxres = params["imgmaxres"]
reslimit = params["reslimit"]
hvals = range(-params["hrange"], params["hrange"]+1)
kvals = range(-params["krange"], params["krange"]+1)
lvals = range(-params["lrange"], params["lrange"]+1)
a = numpy.array(params["aucvec"])
b = numpy.array(params["bucvec"])
c = numpy.array(params["cucvec"])
maxresrad = (imgmaxres/reslimit)*(0.5*imgsize)
circrad =  params["circrad"]   

### open logfile
logout = open("logfile.txt", "a")
now = datetime.datetime.now()
logout.write("\n3_spot_index \t%s\n" % now.strftime("%Y-%m-%d %H:%M"))
logout.write("vectors:\t%s\n\t\t%s\n\t\t%s\n" % (a,b,c))


################ build the list of all the matrix of hk coords

hklindices = []
for eachhval in kvals:
    for eachkval in hvals:
        for eachlval in lvals:
            hklindices.append((eachhval, eachkval, eachlval))


###############  Get the images to process from the list

with open('imagelist.txt') as pfile:
    imagestoprocess = pfile.read().splitlines()
print "processing "+str(len(imagestoprocess))+" images\n"


data = json.load(open("cataspot.json"))
for eachimage in imagestoprocess:
    spots2draw = []
## get the image specific parametrs

    imgsplit = eachimage.split(".")
    imagename = imgsplit[0]+".gif"
    output = file(imgsplit[0]+"_spots.csv", "w")
    
###### calculate the vectors for the reference spots and their dot product:

    xcenter = data["data"]["images"][eachimage]["beamcenter"][0][0]
    ycenter = data["data"]["images"][eachimage]["beamcenter"][0][1]
    theta =float(data["data"]["images"][eachimage]["tiltangle"])
    refpoints = data["data"]["images"][eachimage]["references"]
    bscentx = data["data"]["images"][eachimage]["beamstopcenter"][0][0]
    bscenty = data["data"]["images"][eachimage]["beamstopcenter"][0][1]
    bsrtopy = bscenty-100
    bsrboty = bscenty+100
    bsrad = params["beamstoprad"]
    imagename = imgsplit[0]+".gif"

    refvectors = []
    for i in refpoints:
        x = i[0] - xcenter
        y = (ycenter - i[1])*math.cos(theta*math.pi/180)
        z = (ycenter - i[1])*math.sin(theta*math.pi/180)
        refvectors.append(numpy.array([x,y,z]).reshape(3,1))

## make the unit matrix

    unitmatrix = numpy.array([[a[0],b[0],c[0]],[a[1],b[1],c[1]],[a[2],b[2],c[2]]])
    uinv =  numpy.linalg.inv(unitmatrix)
    r1 = numpy.dot(uinv, refvectors[0])
    r2 = numpy.dot(uinv, refvectors[1])
    
    r1round = numpy.array([int(round(r1.item((0,0)),0)),int(round(r1.item((1,0)),0)),int(round(r1.item((2,0)),0))])
    r2round = numpy.array([int(round(r2.item((0,0)),0)),int(round(r2.item((1,0)),0)),int(round(r2.item((2,0)),0))])

    checkplane = numpy.cross(r1round,r2round)
    print "%s\t" % eachimage,
    
## test every hkl index
    spots = []
    for i in hklindices:
        hkltest = numpy.array([i[0]*checkplane[0],i[1]*checkplane[1],i[2]*checkplane[2]])

        if  -lauezonethresh*numpy.linalg.norm(checkplane) < int(hkltest[0])+int(hkltest[1])+int(hkltest[2]) < numpy.linalg.norm(checkplane)*lauezonethresh:
            spots.append((i[0],i[1],i[2])) 
######## calculate the xyz coordinates of all of miller indices
    xydic = {}          # miller index returns point xy
    xyzvectors = {}        # miller index returns x,y,z vector

    for i in spots:
        x = (i[0]*a[0])+(i[1]*b[0])+(i[2]*c[0])
        y = (i[0]*a[1])+(i[1]*b[1])+(i[2]*c[1])
        z = (i[0]*a[2])+(i[1]*b[2])+(i[2]*c[2])
        xyzvectors[i] = numpy.array([x,y,z])
    
        xproj = x + xcenter
        yproj = ycenter - (y/(math.cos(theta*(math.pi/180))))
        xydic[i] = (round(xproj,1),round(yproj,1))
        if xydic[i][1] > 0 and xydic[i][0] > 0:
            if math.sqrt((math.ceil(xydic[i][0]) - xcenter)**2 + (math.ceil(xydic[i][1]) - ycenter)**2) < maxresrad:
                if math.ceil(xydic[i][0]) < bscentx+bsrad or (math.ceil(xydic[i][0]) > bscentx+bsrad and (bsrtopy-circrad > math.ceil(xydic[i][1]) or math.ceil(xydic[i][1]) > bsrboty+circrad)):
                    if (math.ceil(xydic[i][0]) - bscentx)**2 + (math.ceil(xydic[i][1])-bscenty)**2 > (bsrad**2)+circrad:
                        output.write("%s,%s,%s\t%s,%s\t%s\n" % (i[0],i[1],i[2],xydic[i][0],xydic[i][1],math.ceil(xydic[i][0])+(4096*math.ceil(xydic[i][1]))))
                        spots2draw.append((xydic[i][0],xydic[i][1]))
### draw the images:

    imagename = imgsplit[0]+".gif"     
    csvfile = imgsplit[0]+"_spots.csv"
    spotlist = []
    for i in spots2draw:
        spotlist.append("circle %s,%s %s,%s" % (i[0],i[1],float(i[0])+circrad, float(i[1])+circrad))
    
    print "\t%s spots" % len(spotlist)
    command  = ('convert -size 4096x4096 xc:Black -font Arial -pointsize 10 -stroke Red -strokewidth 2 -fill none -transparent black -draw "%s" +compress points_px.tif' % ' '.join(spotlist))
    os.system('%s' % command)
    os.system("convert %s %s -gravity center -compose over -composite %s" % (imagename, "points_px.tif",imgsplit[0]+"_spots.gif"))