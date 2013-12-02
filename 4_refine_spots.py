#! /usr/bin/env python

### use mass centering to find centers of spots

## Matt Iadanza 2013-07-23

#### imports:
import json
from PIL import Image
import os
import csv
import math
import numpy
import json
import datetime

### get yer data

with open('imagelist.txt') as pfile:
    imagestoprocess = pfile.read().splitlines()
print "processing "+str(len(imagestoprocess))+" images"

data = json.load(open('cataspot.json'))
params = json.load(open('params.json'))

boxsize = params["circrad"]
fullbox = 2*boxsize
drawlist = []
drawlist2 = []

### open logfile
logout = open("logfile.txt", "a")
now = datetime.datetime.now()
logout.write("\n4_refine_spots \t%s\n" % now.strftime("%Y-%m-%d %H:%M"))


for eachimage in imagestoprocess:
    print eachimage
    imgsplit = eachimage.split('.')
    spots = []
    split = []
    imageroot = imgsplit[0]
    pic = Image.open(eachimage)
    values = csv.reader(open(imageroot+'_spots.csv', 'rb'), delimiter='\t')
    output = open(imageroot+"_spots_r.csv", "w")
    rawout = open(imageroot+"_spotrefinement_raw.txt", "w")
    rawout.write("spot\t\toriginal xy\tnewxy\t\t\tdisplacement\n")
    for i in values:
        split = i[1].split(",")
        spots.append([float(split[0]),float(split[1]),i[0]])
    rootname = eachimage.split('.')

## for each spot make a spot sub array 
    for eachspot in spots:    
        spotbox = (int(float(eachspot[0])-boxsize),int(float(eachspot[1])-boxsize),int(float(eachspot[0])+boxsize),int(float(eachspot[1])+boxsize))
        spotregion = pic.crop(spotbox)
        spot = numpy.array(spotregion.getdata()).reshape(spotregion.size[0], spotregion.size[1])
   
## calculate the background mask
        zeros = numpy.zeros((fullbox,fullbox), numpy.int)
        
        
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
        bgmean = numpy.mean(bgmask[numpy.nonzero(bgmask)])

### calculate the spot mask 
        zeros = numpy.zeros((fullbox,fullbox), numpy.int)
        
        
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

## calculate center score
        rows = spot.sum(axis = 0)   
        cols = spot.sum(axis = 1)
        displacement = (rows.argmax(axis = 0),cols.argmax(axis = 0))
        dispvalues = (float(displacement[0]-boxsize), float(displacement[1]-boxsize))

## calculate the new xy
        newxy = []
        newxy = [eachspot[0]+dispvalues[0],eachspot[1]+dispvalues[1]]
        drawlist.append("circle %s,%s %s,%s" % (float(eachspot[0]),float(eachspot[1]),float(eachspot[0]),float(eachspot[1])+boxsize))
        drawlist2.append("circle %s,%s %s,%s" % (float(newxy[0]),float(newxy[1]),float(newxy[0]),float(newxy[1])+boxsize))
        
        rawout.write("%11s\t%s,%s\t%s  \t%s\n" % (eachspot[2],eachspot[0],eachspot[1],newxy,dispvalues))
        
        output.write("%s\t%s,%s\n" % (eachspot[2],newxy[0],newxy[1]))

#### draw the images
    command1  = ('convert -size 4096x4096 xc:Black -stroke Blue -strokewidth 1 -fill none -transparent black -draw "%s" points_old.gif' % ' '.join(drawlist))    
    command2 = ('convert -size 4096x4096 xc:Black -stroke Yellow -strokewidth 1 -fill none -transparent black -draw "%s" points_new.gif' % ' '.join(drawlist2))
    os.system('%s' % command1)
    os.system('%s' % command2)
    drawlist = []
    drawlist2 = []
    
    os.system("convert %s %s -gravity center -compose over -composite %s" % (imageroot+".gif", "points_old.gif",imageroot+"_refine_spots.gif"))
    os.system("convert %s %s -gravity center -compose over -composite %s" % (imageroot+"_refine_spots.gif", "points_new.gif",imageroot+"_refine_spots.gif"))

## cleanup

os.system('rm points_old.gif')
os.system('rm points_new.gif')