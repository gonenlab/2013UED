#! /usr/bin/env python

## imports:

import numpy
import math

# Takes file of all combined intensities makes lists of intensities for corresponding symmetry mates in p422 symmetry.
# chooses which ones to sum based on a threshold (tfactor), sums them, and calculates Rmerge for the data set. version 2

# Matt Iadanza 2013-07-18

############
# some variables
###########

tfactor = 0.6                           # threshold for keeping spots.  MaxRmerge = 1-t 
hrange = range(0,31)                    #
krange = range(0,31)                    # the maximum possible values for h,k and l
lrange = range(0,31)                    #
rawoutput = open("rawoutput_merge2.txt", "w")
output = open("f_mergedint.txt", "w")
vers = 2

### open the combined intensities files 
with open('combint.txt') as pfile:
    lines = pfile.read().splitlines()


# make dictionary structure to contain the merged spots

spotint = {}    # input spot return all intensities for it

for h in hrange:
    for k in krange:
        for l in lrange:
            if h >= k:
                spotint[(h,k,l)] = []


## put intensity values into the spotint dict

for each in lines:
    data = each.split("\t")
    h,k,l = int(data[0]),int(data[1]),int(data[2])

# first deal with all nonzero miller indices
## hkl 'root' spots (+++) and their freidel pairs (---) -h -k -l

    if h > 0 and k > 0 and l >= 0 and  h >= k:
        spotint[(h,k,l)].append(float(data[5]))
    if h < 0 and k < 0 and l < 0 and  h < k:
        spotint[(-h,-k,-l)].append(float(data[5]))

## khl symmetry mates  (+++) and their freidel pairs (---) -k -h -l
    if h > 0 and k > 0 and l >= 0 and  h < k:
        spotint[(k,h,l)].append(float(data[5]))
    if h < 0 and k < 0 and l < 0 and  h >= k:
        spotint[(-k,-h,-l)].append(float(data[5]))

## h k -l (++-) symmetry mates and freidel pairs -h -k l (--+)
    if h > 0 and k > 0 and l < 0 and h >= k:
        spotint[(h,k,-l)].append(float(data[5]))
    if h < 0 and k < 0 and l >= 0 and h < k:
        spotint[(-h,-k,l)].append(float(data[5]))

## k h -l (++-) symmetry mates and freidel pairs  -k -h l  (--+) 
    if h > 0 and k > 0 and l < 0 and h < k:
        spotint[(k,h,-l)].append(float(data[5]))
    if h < 0 and k < 0 and l >= 0 and h >= k:
        spotint[(-k,-h,l)].append(float(data[5]))
        
## h -k l (+-+) symmetry maes and freidel pairs -h k -l (-+-)
    if h > 0 and k < 0 and l >= 0 and h >= -k:
        spotint[(h,-k,l)].append(float(data[5]))
    if h < 0 and k > 0 and l < 0 and -h >= k:
        spotint[(-h,k,-l)].append(float(data[5]))

## h -k -l (+--) symmetry mates and freidel pairs -h k l (-++)
    if h > 0 and k < 0 and l < 0 and h >= -k:
        spotint[(h,-k,-l)].append(float(data[5]))
    if h < 0 and k > 0 and l >= 0 and -h >= k:
        spotint[(-h,k,l)].append(float(data[5]))

## k -h -l (+--) symmetry mates and freidel pairs -k h l (-++)
    if h < 0 and k > 0 and l < 0 and -h < k:
        spotint[(k,-h,-l)].append(float(data[5]))
    if h > 0 and k < 0 and l >= 0 and h < -k:
        spotint[(-k,h,l)].append(float(data[5]))
        
## k -h l (-+-)symmetry mates and freidel pairs -k h -l
    if h < 0 and k > 0 and l >= 0 and -h < k:
        spotint[(k,-h,l)].append(float(data[5]))
    if h > 0 and k < 0 and l < 0 and h < -k:
        spotint[(-k,h,-l)].append(float(data[5]))

## next deal with the special-case zero miller indices
    if h == 0 and k != 0:
        if k < 0:
            k = k*-1
        if l < 0:
            l = l*-1
        spotint[(k,h,l)].append(float(data[5]))
    if k == 0 and h != 0:
        if h < 0:
            h = h*-1
        if l < 0:
            l = l*-1
        spotint[(h,k,l)].append(float(data[5]))
    if k == 0 and h == 0:
        if l < 0:
            l = l *-1
        spotint[(h,k,l)].append(float(data[5]))
        
## go over spotint dict and remove any entries with no intensity values

for i in spotint.keys():
    if len(spotint[i]) == 0:
        del(spotint[i])
    
count = 0    
for i in spotint:
    count = count+len(spotint[i])

## go over spotint dict and remove any intensities that are outside the threshold of the max value:

keepvals = {}            # input spot(root spot) and return all intensity values that meet the threshold
finalvalues = {}
vals = []

for i in spotint:
    maxval = max(spotint[i])
    keepvals[i] = []
    for n in spotint[i]:
        if n >= tfactor*float(maxval):
            keepvals[i].append(n)

### Calculate the statistcs

rmergerunning = []
rmergecount = 0

for i in keepvals:
    nums = []
    for n in keepvals[i]:
        nums.append(float(n))
    mean = numpy.mean(nums)
    rmergecount = rmergecount + len(nums)
    for n in nums:
            val = abs(n-mean)/mean
            rmergerunning.append(float(val))
    if len(keepvals[i]) == 1:
        intensity = float(nums[0])
        stdi = math.sqrt(intensity)
        sigi = stdi
        sigf = math.sqrt(stdi)
    if len(keepvals[i]) == 2:
        intensity = sum(nums)/len(nums)
        stdi = math.sqrt(intensity)
        sigi = stdi
        sigf = math.sqrt(stdi)
    if len(keepvals[i]) > 2:       
        intensity = sum(nums)/len(nums)
        stdi = numpy.std(nums)
        sigi = math.sqrt(intensity)
        sigf = math.sqrt(stdi)
    finalvalues[i] = (sigi,sigf,intensity,stdi)
    rawoutput.write("***** %s: %s\n--------\n%s\n\n" % (i, spotint[i],nums))
    output.write("%2i %2i %2i %12.4f %12.4f %12.4f %12.4f\n" % (i[0],i[1],i[2],finalvalues[i][0],finalvalues[i][1],finalvalues[i][2],finalvalues[i][3]))
                 


print "** thresholded merge p422 vers: %2.1f **" % vers
print "max allowable Rmerge: %s" % (1-tfactor)
print "%s intensities" % count
print "%s merged intensities" % len(finalvalues)
print "Rmerge: %s" % (sum(rmergerunning)/rmergecount)




###############################################################################
### Keep for debugging - which spots are going where???
#
### hkl 'root' spots (+++) and their freidel pairs (---) -h -k -l
#    if h >= 0 and k >= 0 and l >= 0 and  h >= k:
#        spotint[(h,k,l)].append((h,k,l))
#    if h <= 0 and k <= 0 and l <= 0 and  h <= k:
#        spotint[(-h,-k,-l)].append((h,k,l))
#
### khl symmetry mates  (+++) and their freidel pairs (---) -k -h -l
#    if h >= 0 and k >= 0 and l >= 0 and  h <= k:
#        spotint[(k,h,l)].append((h,k,l))
#    if h <= 0 and k <= 0 and l <= 0 and  h >= k:
#        spotint[(-k,-h,-l)].append((h,k,l))
#
### h k -l (++-) symmetry mates and freidel pairs -h -k l (--+)
#    if h >= 0 and k >= 0 and l <= 0 and h >= k:
#        spotint[(h,k,-l)].append((h,k,l))
#    if h <= 0 and k <= 0 and l >= 0 and h <= k:
#        spotint[(-h,-k,l)].append((h,k,l))
#
### k h -l (++-) symmetry mates and freidel pairs  -k -h l  (--+) 
#    if h >= 0 and k >= 0 and l <= 0 and h <= k:
#        spotint[(k,h,-l)].append((h,k,l))
#    if h <= 0 and k <= 0 and l >= 0 and h >= k:
#        spotint[(-k,-h,l)].append((h,k,l))
#        
### h -k l (+-+) symmetry maes and freidel pairs -h k -l (-+-)
#    if h >= 0 and k <= 0 and l >= 0 and h >= -k:
#        spotint[(h,-k,l)].append((h,k,l))
#    if h <= 0 and k >= 0 and l <= 0 and -h >= k:
#        spotint[(-h,k,-l)].append((h,k,l))
#
### h -k -l (+--) symmetry mates and freidel pairs -h k l (-++)
#    if h >= 0 and k <= 0 and l <= 0 and h >= -k:
#        spotint[(h,-k,-l)].append((h,k,l))
#    if h <= 0 and k >= 0 and l >= 0 and -h >= k:
#        spotint[(-h,k,l)].append((h,k,l))
#
### k -h -l (+--) symmetry mates and freidel pairs -k h l (-++)
#    if h <= 0 and k >= 0 and l <= 0 and -h <= k:
#        spotint[(k,-h,-l)].append((h,k,l))
#    if h >= 0 and k <= 0 and l >= 0 and h <= -k:
#        spotint[(-k,h,l)].append((h,k,l))
#        
### k -h l (-+-)symmetry mates and freidel pairs -k h -l
#    if h <= 0 and k >= 0 and l >= 0 and -h <= k:
#        spotint[(k,-h,l)].append((h,k,l))
#    if h >= 0 and k <= 0 and l <= 0 and h <= -k:
#        spotint[(-k,h,-l)].append((h,k,l))
