#! /usr/bin/env python
## imports

## Finds the average unit cell vector using find cell_params.json which is created by user or from imageJ files and i2cfiles.py

# Matt Iadanza 2013-07-05

import json
import os
import math
import numpy
from collections import Counter
import datetime

# get some variables:
data = json.load(open('cataspot.json'))
output = open("magfreqs.csv", "w")
drawplot = open("plot", "w")
params = json.load(open('params.json'))


### prepare the logfile:
logout = open("logfile.txt", "w")
logout.write("** MxED logfile **\nusing cataspot.json created: %s\nuser name: %s\n\nPrograms run and results\n------------------------\n" % (data["metadata"]["date"],data["metadata"]["username"]))
logout.write("Initial Parameters\n------------------\nimage size:\t%sx%s\nedge res:\t%s\nh,k,l range:\t%s,%s,%s\nres cutoff:\t%s\n\n" % (params["imgsize"],params["imgsize"],params["imgmaxres"],params["hrange"],params["krange"],params["lrange"],params["reslimit"]))
now = datetime.datetime.now()
logout.write("1_find_lengths \t%s\n" % now.strftime("%Y-%m-%d %H:%M"))

os.system('clear')

#### version check
vers = 1
catspotvers = data["metadata"]["file version"]
print "** Rough Unit Cell Determination vers %s **" % round(vers,2)
if vers != catspotvers:
    print "datafile version mismatch %s/%s - may cause errors" % (round(vers,2),round(catspotvers,2))
else:
    print "version check passed"
##### do it

print ""
print "make the spolist dictionarys"
spotlist= []   # input raw spot number returns numpy array with vector 
xydic = {}  # original x,y coords to image#,spot#
xyzdic = {} # calculated x,y,z coords to image#,spot#

with open('imagelist.txt') as f:
    images2process = f.read().splitlines()

spotrange = []
for eachimage in images2process:
    theta = float(data["data"]["images"][eachimage]["tiltangle"])
    for i in data["data"]["images"][eachimage]["spots"]:
        ox = i[0]
        oy = i[1]
        x = i[0] - data["data"]["images"][eachimage]["beamcenter"][0][0]
        y = -(oy - data["data"]["images"][eachimage]["beamcenter"][0][1])*math.cos(theta*math.pi/180.0)
        z = -(oy - data["data"]["images"][eachimage]["beamcenter"][0][1])*math.sin(theta*math.pi/180.0)
        spotrange.append(numpy.array([x,y,z]))

print "build spotlist: %s spots picked" % len(spotrange)

#### calculate the difference vectors for all combinations of spots

print"subtract every vector from every other"
diffvectors = []
for n in spotrange:
     for i in spotrange:
        diff = numpy.subtract(n,i)
        if diff[0] and diff[1] and diff[2] != 0:
            diffvectors.append(diff)

print "%s difference vectors found" % len(diffvectors)

### count the distribution of the difference vectors

count = Counter()
for i in diffvectors:
    mag = round(numpy.linalg.norm(i),1)
    if mag != 0:
        count[mag] += 1

ymax = float(max(count.values()))/2+5

for i in count:
    if 0 < i < 250:
        output.write("%s %s\n" % (i,count[i]/2))


########  make the plot
os.system("touch outgraph.eps")
print "drawing plot"
drawplot.write("set terminal postscript eps enhanced color font 'Verdana,10'\n")
drawplot.write("set yrange [0:%s]\n" % ymax)
drawplot.write('set title "Unit Cell Predictions" font "Arial,12"\n')
drawplot.write("set output 'outgraph.eps'\n")
drawplot.write("set style line 1 lc rgb '#0060ad' lt 1 lw 2 pt 7 ps 1.5\n")
drawplot.write("plot 'magfreqs.csv' with impulses\n")

os.system("gnuplot plot")
os.system("open outgraph.eps")
