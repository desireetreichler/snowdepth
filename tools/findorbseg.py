# -*- coding: utf-8 -*-
"""
Created on Thu Mar 12 11:43:37 2020

@author: desireet
"""

import numpy as np

# find the ICESat-2 orbital segment of a coordinate (latitude)
# according to the table here: https://nsidc.org/data/ATL08/versions/2
def findorbseg(b, *args): #bottom, upper latitude / existing list of orbsegs to add on
	
	otherargs= list(*args)# take only the first other argument
	#print(args)
	#print(otherargs)
	#print(type(otherargs))
	if type(otherargs) == list or type(otherargs) == tuple:
		orbseg=otherargs
	else: orbseg=list([])
#	print(orbseg)
		# north
	if np.array(b<27) & np.array(b>=0):
		orbseg=checkappend(orbseg,[1,7])
	if np.array(b<59.5) & np.array(b>=27):
		orbseg=checkappend(orbseg,[2,6])
	if np.array(b<80) & np.array(b>=59.5):
		orbseg=checkappend(orbseg,[3,5])
	if np.array(b>=80):
		orbseg=checkappend(orbseg,[4])
		# south
	if np.array(b>-27) & np.array(b<=0):
		orbseg=checkappend(orbseg,[8,14])
	if np.array(b>-59.5) & np.array(b<=27):
		orbseg=checkappend(orbseg,[9,13])
	if np.array(b>-80) & np.array(b<=59.5):
		orbseg=checkappend(orbseg,[10,12])
	if np.array(b<=-80):
		orbseg=checkappend(orbseg,[11])
	return orbseg



def checkappend(orbseg,a):
	# check whether the number already is in orbseg and if not, add it
	for i in a:
		istrue = np.array(i) in np.array(orbseg)
		
		if not istrue:
			orbseg.append(i)
	return orbseg