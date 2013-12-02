#! /usr/bin/env python
## imports

## Finds the average unit cell vector using cataspot.json which is created by cataspot

# Matt Iadanza 20-07-05

import json
import os
import math
import numpy
import datetime


#### user input variables


ea = 55.0       # expected a unit cell length
athresh = 1  
eb = 55.0       # expected b unit cell length
bthresh = 1
ec = 112.0      # expected c dimension
cthresh = 1


#### get yer data:
data = json.load(open('cataspot.json'))
numberofimages = len(data["data"]["images"])
output = open("output_cellfind.txt", "w")
params = json.load(open('params.json'))
imgsize = params['imgsize']
maxres = params['imgmaxres']

### open logfile
logout = open("logfile.txt", "a")
now = datetime.datetime.now()
logout.write("\n2_calc_uc vectors \t%s\n" % now.strftime("%Y-%m-%d %H:%M"))
logout.write("vector lengths(a,b,c): %s,%s,%s    thresholds(a,b,c): %s,%s,%s\n" % (ea,eb,ec,athresh,bthresh,cthresh))

#### init and version check
os.system('clear')
vers = 1
catspotvers = data["metadata"]["file version"]
print "** Unit Cell Vector Determination vers %s **" % round(vers,2)
if vers != catspotvers:
    print "datafile version mismatch %s/%s - may cause errors" % (round(vers,2),round(catspotvers,2))
else:
    print "version check -- passed"
    
##### do it

print ""
print "make the spolist dictionarys"
with open('imagelist.txt') as f:
    images2process = f.read().splitlines()

### calculate 3-D coordinates for each spot

##------- ewald sphere correction (z dimension) function----------
def ewaldcorr(xdim,ydim):
    wavelength = 0.025
    oneoverlambda = 1/(wavelength * (1/(0.5 * imgsize * maxres)))
    a = oneoverlambda/math.sqrt(xdim**2+ydim**2+oneoverlambda**2)
    deltaz = (1-a)*oneoverlambda
    return deltaz
##-----------------------------------------------------------------

spotlist = []
for eachimage in images2process:
    theta = float(data["data"]["images"][eachimage]["tiltangle"])
    for i in data["data"]["images"][eachimage]["spots"]:
        ox = i[0]
        oy = i[1]
        x = i[0] - data["data"]["images"][eachimage]["beamcenter"][0][0]
        y = -(oy - data["data"]["images"][eachimage]["beamcenter"][0][1])*math.cos(theta*math.pi/180.0)
        z = -(oy - data["data"]["images"][eachimage]["beamcenter"][0][1])*math.sin(theta*math.pi/180.0) - ewaldcorr(x,y)
        spotlist.append(numpy.array([x,y,z]))

print "build spotlist: %s spots picked" % len(spotlist)
logout.write("%s spots\n" % len(spotlist))

#### calculate difference vectors and sort out those that could be unit cells 
print"subtract every vector from every other and determine the magnitude of results - keep possible unit cell vectors\n"

poucvs = []
for n in spotlist:
    for i in spotlist:
        diffvec = numpy.subtract(n,i)
        sub = numpy.linalg.norm(diffvec)
        if i[0] != n[0] and i[1] != n[1] and i[2] != n[2] and ((ea-athresh < sub < ea+athresh)or(eb-bthresh < sub < eb+bthresh)or(ec-cthresh < sub < ec+cthresh)):
            poucvs.append(diffvec)

print "%s possible unit cell vectors found" % len(poucvs)

### sort the unit cell vectors into parallel groups
### compare all vectors to xy normal test vectors 
### use results to sort into roughly (within threshold) parallel groups

print "sorting unit cell vectors into parallel groups"

testvector1 = numpy.array([0,0,10000])
testvector2= numpy.array([0,10000,0])
testvector3= numpy.array([10000,0,0])
testvector4 = numpy.array([5000,5000,5000])
oset1 =  []
oset2 = []
oset3 = []


##----------- function to calculate the angle betwen two vectors ----------

def calcang(a,b):
    v12 = numpy.dot(a,b)
    v1mag = numpy.linalg.norm(a)
    v2mag = numpy.linalg.norm(b)
    cosphi = abs((v12)/(v1mag*v2mag))
    return round((180/math.pi)*math.acos(round(cosphi,12)),0)
##-------------------------------------------------------------------------

##----------- function to calculate the angle betwen two vectors - nonabsoulte ----------

def calcangnonabs(a,b):
    v12 = numpy.dot(a,b)
    v1mag = numpy.linalg.norm(a)
    v2mag = numpy.linalg.norm(b)
    cosphi = ((v12)/(v1mag*v2mag))
    return round((180/math.pi)*math.acos(round(cosphi,12)),0)
##-------------------------------------------------------------------------



### initial grouping for parallel vectors########
group01 = []
group02 = []
group03 = []
group10 = []
group12 = []
group13 = []
group20 = []
group21 = []
group23 = []
group30 = []
group31 = []
group32 = []
##################################################

## identify which reference vectro each is closest to and furthest from (in terms of angles) use to classify into roughly parallel groups
orientations = {}
for i in poucvs:
    
    orientations[(i[0],i[1],i[2])] = [calcang(i,testvector1),calcang(i,testvector2),calcang(i,testvector3),calcang(i,testvector4)]

for i in orientations:
    scores = [orientations[i][0], orientations[i][1],orientations[i][2],orientations[i][3]]
    maxscore = scores.index(max(scores))
    minscore = scores.index(min(scores))
    i = numpy.array(i)
    if sum(numpy.cross(i,testvector1)) < 0:
            i = numpy.multiply(-1,i)
    globals()['group'+str(maxscore)+str(minscore)].append(i)




##### refine the groups - only keep vectors within 1 STD of eth mean for each of th four reference vectors
#### goodgroups
good01 = []
good02 = []
good03 = []
good10 = []
good12 = []
good13 = []
good20 = []
good21 = []
good23 = []
good30 = []
good31 = []
good32 = []

### ---------group calc function----------------------------------------------->>>>>
def groupcalc(x,y):
    
    ang1 = []
    ang2 = []
    ang3 = []
    ang4 = []
    
    
    
    
    for i in x:
            ang1.append(calcangnonabs(i,testvector1))
            ang2.append(calcangnonabs(i,testvector2))
            ang3.append(calcangnonabs(i,testvector3))
            ang4.append(calcangnonabs(i,testvector4))

        
        
        
    for i in x:
        score1 = int(abs(numpy.mean(ang1)-calcangnonabs(i,testvector1))/(numpy.std(ang1)+.00001))
        score2 = int(abs(numpy.mean(ang2)-calcangnonabs(i,testvector2))/(numpy.std(ang2)+.00001))
        score3 = int(abs(numpy.mean(ang3)-calcangnonabs(i,testvector3))/(numpy.std(ang3)+.00001))
        score4 = int(abs(numpy.mean(ang4)-calcangnonabs(i,testvector4))/(numpy.std(ang4)+.00001))
        

        if score1+score2+score3+score4 < 1:
            globals()['good'+str(y)].append(i)
  
###--------------------------------------------------------------------------------->>>>>


scaled01 = []
scaled02 = []
scaled03 = []
scaled10 = []
scaled12 = []
scaled13 = []
scaled20 = []
scaled21 = []
scaled23 = []
scaled30 = []
scaled31 = []
scaled32 = []






###### output everything to screen ####
## make raw output file later?

## rough groups
for i in ('01','02','03','10','12','13','20','21','23','30','31','32'):
    if len(globals()['group'+i]) > 0:
        groupcalc(globals()['group'+i],i)

####### scale the vectors


##----calc minimum funct---------------
def calcmin(x):
    allmags = []
    for vector in x:
        allmags.append(numpy.linalg.norm(vector))
    return min(allmags)
##--------------------------------------

for i in ('01','02','03','10','12','13','20','21','23','30','31','32'):
    for x in globals()['good'+i]:
        if numpy.linalg.norm(x) > 2*calcmin(globals()['good'+i]):
            globals()['scaled'+i].append(x/2)
        else:
            globals()['scaled'+i].append(x) 

## -- list printing function------
def printlist(x):
    print "-----"
    for i in x:
        print i
##-----------------------------









##--------- funct to calculate mean vectors -----------------
def vecmean(x):
    if len(x) > 0:
        iss = []
        js = []
        ks = []
        for i in x:
            iss.append(i[0])
            js.append(i[1])
            ks.append(i[2])
        mean = numpy.array([numpy.mean(iss),numpy.mean(js),numpy.mean(ks)])
        return [mean, numpy.linalg.norm(mean)]
    else:
        return ["no vector","n/a"]
##----------------------------------------------------------


candidates = []
## final mean vectors
print '\n------- Candidate vectors -------'
for i in ('01','02','03','10','12','13','20','21','23','30','31','32'):
    if len(globals()['good'+i]) > 2:
        print"%s\t%s  \tmag: %s  \tscore: %s" % (i,vecmean(globals()['scaled'+i])[0],vecmean(globals()['scaled'+i])[1], len(globals()['scaled'+i]))
        logout.write("%s\t%s  \tmag: %s  \tscore: %s\n" % (i,vecmean(globals()['scaled'+i])[0],vecmean(globals()['scaled'+i])[1], len(globals()['scaled'+i])))
        candidates.append([i,vecmean(globals()['scaled'+i])[0]])
print"\n------- Angles Between Candidates -------"
print "\t",
logout.write("\t")
for i in candidates:
    print str(i[0])+"\t",
    logout.write(str(i[0])+"\t")
print ""
logout.write("\n")
for i in candidates:
    anglist = []
    for n in candidates:
        anglist.append(str(calcangnonabs(i[1],n[1])))
    print "%s\t%s" % (i[0],'\t'.join(anglist))
    logout.write("%s\t%s\n" % (i[0],'\t'.join(anglist)))