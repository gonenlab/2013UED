#!/usr/bin/env python

###################


## refines unit cell vectors based on spot centering...
## matt iadanza 2013-09-26

## modifications of integrate2_batch.py so it can handle data from cataspot.
## Matt Iadanza 2013-09-16

# findsw intensity or spots specified in a csv file.  Version 2, trying to get more accurate intensities, wih improved
# background subtraction and some sor of thresholding

# Matt Iadanza 2013-07-22

################### Variables

sqthresh = 1.15     # spot quality threshold % of mean background
centererrorthresh = 3


###################imports
from PIL import Image
import os
import csv
import math
import numpy
import json
import datetime

## get global parameters from parfile.

data = json.load(open('cataspot.json'))
params = json.load(open('params.json'))
imgsize = params["imgsize"]      
boxsize =  params["circrad"]   
radius = range(1,boxsize+1)            
fullbox = 2*boxsize
output = open("refined.txt","w")


### open logfile
logout = open("logfile.txt", "a")
now = datetime.datetime.now()
logout.write("\n5_UCR_index \t%s\n" % now.strftime("%Y-%m-%d %H:%M"))
logout.write("Spot intensity threshold: %sx background\tCenter error threshold: %s px\n" % (sqthresh,centererrorthresh))


#### init and version check
os.system('clear')
vers = 1
catspotvers = data["metadata"]["file version"]
print "** Indexing for Unit Cell Vector Determination vers %s **" % round(vers,2)
if vers != catspotvers:
    print "datafile version mismatch %s/%s - may cause errors" % (round(vers,2),round(catspotvers,2))
else:
    print "version check -- passed"




## get list of images to process

with open('imagelist.txt') as pfile:
    imagestoprocess = pfile.read().splitlines()
print "processing "+str(len(imagestoprocess))+" images"

################### Start processing images

for eachimage in imagestoprocess:
    qualityspots = {}       # input x,y return intensity
    theta =float(data["data"]["images"][eachimage]["tiltangle"])
    
## get the image specific parametrs

    imgsplit = eachimage.split(".")
    spotthreshold = params["integrationthresh"]
    rawoutput = open(imgsplit[0]+"_rawint.txt", "w")
    
## open the image:

    pic = Image.open(eachimage)

## prepare to process each spot

    values = csv.reader(open(imgsplit[0]+'_spots_r.csv', 'rb'), delimiter='\t')
    indexdic = {}   # input xy return miller index
    coordsdic = {}  # input miller index return xy
    for row in values:
        indexdic[row[1]] = row[0]
        coordsdic[row[0]] = row[1]

## for each spot make a spot sub array 
    xy = []
    for eachspot in indexdic:    
        xy = eachspot.split(",")
        spotbox = (int(float(xy[0])-boxsize),int(float(xy[1])-boxsize),int(float(xy[0])+boxsize),int(float(xy[1])+boxsize))
        spotregion = pic.crop(spotbox)
        spot = numpy.array(spotregion.getdata()).reshape(spotregion.size[0], spotregion.size[1])
        rawoutput.write( "\nRaw spot %s\n" % indexdic[eachspot])
        rawoutput.write( '\n'.join('\t'.join(str(cell) for cell in row) for row in spot))
        
## calculate the background mask
        zeros = numpy.zeros((fullbox,fullbox), numpy.int)
        
        rawoutput.write( "\n-bgmatrix-\n")
        
        maskcoords = []
        maskrange = range(0,boxsize+1)
        center = [boxsize,boxsize]
        for i in maskrange:
            for n in maskrange:
                if math.sqrt((boxsize-i)**2 + (boxsize-n)**2) > boxsize:
                    maskcoords.append([i,n])
                    maskcoords.append([(2*boxsize-1)-i,n])
                    maskcoords.append([i,(2*boxsize-1)-n])
                    maskcoords.append([(2*boxsize-1)-i,(2*boxsize-1)-n])
        for i in maskcoords:
            zeros[i[0],i[1]] = 1
        bgmask = numpy.multiply(zeros,spot)
        rawoutput.write( '\n'.join('\t'.join(str(cell) for cell in row) for row in bgmask))
        bgmean = numpy.mean(bgmask[numpy.nonzero(bgmask)])
        rawoutput.write( "\nmean BG intensity: %s\n" % bgmean)

### calculate the spot mask and quality test spots
        zeros = numpy.zeros((fullbox,fullbox), numpy.int)
        
        rawoutput.write("\n-raw spotmatrix-\n")
        
        maskcoords = []
        maskrange = range(0,boxsize+1)
        center = [boxsize,boxsize]
        for i in maskrange:
            for n in maskrange:
                if math.sqrt((boxsize-i)**2 + (boxsize-n)**2) < boxsize:
                    maskcoords.append([i,n])
                    maskcoords.append([(2*boxsize-1)-i,n])
                    maskcoords.append([i,(2*boxsize-1)-n])
                    maskcoords.append([(2*boxsize-1)-i,(2*boxsize-1)-n])
        for i in maskcoords:
            zeros[i[0],i[1]] = 1
        spotmask = numpy.multiply(zeros,spot)
        rawoutput.write( '\n'.join('\t'.join(str(cell) for cell in row) for row in spotmask))

## calculate center score
        rows = spot.sum(axis = 0)   
        cols = spot.sum(axis = 1)
        displacement = (rows.argmax(axis = 0),cols.argmax(axis = 0))
        dispvalues = (displacement[0]-boxsize, displacement[1]-boxsize)
        dispscore = abs(dispvalues[0])+abs(dispvalues[1])
        
### background subtract
        bgsubmatrix = numpy.ndarray(shape=(2*boxsize,2*boxsize), dtype=float)
        bgsubmatrix.fill(bgmean)
        rawoutput.write( "\n-BG subtracted spotmatrix-\n")
        subtracted = numpy.subtract(spot,bgsubmatrix)
        subtractedmasked = numpy.multiply(subtracted,zeros)
        rawoutput.write( '\n'.join('\t'.join(str(cell) for cell in row) for row in subtractedmasked.round(0)))
        rawoutput.write("\nspot: %s\t i: %s quality: %s\n\n" % (indexdic[eachspot],numpy.sum(subtractedmasked),numpy.mean(spotmask[numpy.nonzero(spotmask)])/bgmean))
        if numpy.mean(spotmask[numpy.nonzero(spotmask)])/bgmean >= sqthresh and dispscore <= centererrorthresh:
            qualityspots[eachspot] = numpy.sum(subtractedmasked)
            hkl = indexdic[eachspot].split(",")
            print "index:%10s\t intensity: %14s \tquality: %14s center error: %s : %s" % (indexdic[eachspot],numpy.sum(subtractedmasked),numpy.mean(spotmask[numpy.nonzero(spotmask)])/bgmean,dispvalues,dispscore)

## print image summary, write data to output,and draw the picture 
    
    rawoutput.write( "- %s image summary -\n" % eachimage)
    rawoutput.write( "%s/%s spots kept - threshold: %s\n" % (len(qualityspots),len(indexdic),sqthresh))
    
    print ""   
    print "- %s image summary -" % eachimage
    print "%s/%s spots kept - threshold: %s" % (len(qualityspots),len(indexdic),sqthresh)

    drawlist = []
    for spot in qualityspots:
        output.write("%s\t%s\t%s\t%s\n" % (xy[0],xy[1],theta,eachimage))
        xy = []
        xy = spot.split(',')
        drawlist.append("circle %s,%s %s,%s" % (float(xy[0]),float(xy[1]),float(xy[0]),float(xy[1])+boxsize))
                



    command  = ('convert -size 4096x4096 xc:Black -stroke Yellow -strokewidth 2 -fill none -transparent black -draw "%s" points_px.gif' % ' '.join(drawlist))
    os.system('%s' % command)
    os.system("convert %s %s -gravity center -compose over -composite %s" % (imgsplit[0]+".gif", "points_px.gif",imgsplit[0]+"_ucr_spots.gif"))

## cleanup

    os.system('rm points_px.gif')
