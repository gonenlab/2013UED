#!/usr/bin/env python

###################

## modifications of integrate2_batch.py so it can handle data form cataspot.
## Matt Iadanza 2013-09-16

# findsw intensity or spots specified in a csv file.  Version 2, trying to get more accurate intensities, wih improved
# background subtraction and some sor of thresholding

# Matt Iadanza 2013-07-22

################### Variables

sqthresh = 1.00     # spot quality threshold % of mean background
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
output = open("combint.txt","w")

#### init and version check
os.system('clear')
vers = 1
catspotvers = data["metadata"]["file version"]
print "** Intensity measurement %s **" % round(vers,2)
if vers != catspotvers:
    print "datafile version mismatch %s/%s - may cause errors" % (round(vers,2),round(catspotvers,2))
else:
    print "version check -- passed"

### open logfile
logout = open("logfile.txt", "a")
now = datetime.datetime.now()
logout.write("\n7_measure intensities \t%s\n" % now.strftime("%Y-%m-%d %H:%M"))
logout.write("integration thresh:\t%s\ncircle radius:\t\t%s\nbeamstop radius:\t%s\n" % (params["integrationthresh"],boxsize,params["beamstoprad"]))

## get list of images to process

with open('imagelist.txt') as pfile:
    imagestoprocess = pfile.read().splitlines()
print "processing "+str(len(imagestoprocess))+" images"

################### Start processing images
batch = 0

for eachimage in imagestoprocess:
    imgroot = eachimage.split('.')
    batch = batch+1
    qualityspots = {}       # input x,y return intensity
    
## get the image specific parametrs

    imgsplit = eachimage.split(".")
    spotthreshold = params["integrationthresh"]
    rawoutput = open(imgsplit[0]+"_rawint.txt", "w")
    
## open the image:

    pic = Image.open(eachimage)

## prepare to process each spot

    values = csv.reader(open(imgsplit[0]+'_spots.csv', 'rb'), delimiter='\t')
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
        
### background subtract
        bgsubmatrix = numpy.ndarray(shape=(2*boxsize,2*boxsize), dtype=float)
        bgsubmatrix.fill(bgmean)
        rawoutput.write( "\n-BG subtracted spotmatrix-\n")
        subtracted = numpy.subtract(spot,bgsubmatrix)
        subtractedmasked = numpy.multiply(subtracted,zeros)
        rawoutput.write( '\n'.join('\t'.join(str(cell) for cell in row) for row in subtractedmasked.round(0)))
        rawoutput.write("\nspot: %s\t i: %s quality: %s\n\n" % (indexdic[eachspot],numpy.sum(subtractedmasked),numpy.mean(spotmask[numpy.nonzero(spotmask)])/bgmean))
        if numpy.mean(spotmask[numpy.nonzero(spotmask)])/bgmean >= sqthresh and numpy.sum(subtractedmasked) > 0:
            qualityspots[eachspot] = numpy.sum(subtractedmasked)
            hkl = indexdic[eachspot].split(",")
            output.write("%s\t%s\t%s\t-\t%s\t%s\n" % (hkl[0],hkl[1],hkl[2],batch,numpy.sum(subtractedmasked)))
            print "index:%10s\t intensity: %14s \tquality: %14s" % (indexdic[eachspot],numpy.sum(subtractedmasked),numpy.mean(spotmask[numpy.nonzero(spotmask)])/bgmean)

## print image summary, write data to output,and draw the picture 
    
    rawoutput.write( "- %s image summary -\n" % eachimage)
    rawoutput.write( "%s/%s spots kept - threshold: %s\n" % (len(qualityspots),len(indexdic),sqthresh))
    
    print ""   
    print "- %s image summary -" % eachimage
    print "%s/%s spots kept - threshold: %s" % (len(qualityspots),len(indexdic),sqthresh)

    for spot in qualityspots:
        xy = []
        drawlist = []
        for i in qualityspots:
            xy = i.split(",")
            drawlist.append("circle %s,%s %s,%s" % (float(xy[0]),float(xy[1]),float(xy[0]),float(xy[1])+boxsize))
    command  = ('convert -size 4096x4096 xc:Black -stroke Yellow -strokewidth 2 -fill none -transparent black -draw "%s" points_px.gif' % ' '.join(drawlist))
    os.system('%s' % command)
    os.system("convert %s %s -gravity center -compose over -composite %s" % (imgroot[0]+".gif", "points_px.gif",imgroot[0]+"_int.gif"))

## cleanup

    os.system('rm points_px.gif')
