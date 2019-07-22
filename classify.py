import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats
import six
from os.path import isfile
import h5py
np.set_printoptions(threshold=17000)
import requests
import math
import csv

#Load data corresponding to the age gradients in the 2 bins
grad1=np.loadtxt("age_grad5.csv")			#0 to 2Re
grad2=np.loadtxt("age_grad6.csv")			#2Re to 4Re

gal_type=[]									#Array that is assigned number based on how age grad changes between bins 
											#1=sign change, 2=same sign but change in gradient, 3=same gradient, 4=not enough stellar particles

length=grad1.shape[0]
print(length)

#Classify the galaxy types
for i in range(length):
	

	if (grad1[i]>0 and grad2[i]<0) or (grad1[i]<0 and grad2[i]>0):
		gal_type.append(1)
		
	elif grad2[i]>1.1*grad1[i] or grad2[i]<0.9*grad1[i]:				#Difference in grad defined as gradients being at least 10% different
		gal_type.append(2)
		
	elif (grad2[i]==0):
		gal_type.append(4)
	else:
		gal_type.append(3)
		
np.savetxt("gal_type.csv",gal_type,delimiter=',')
print("======================")
one=0
two=0
three=0
four=0
for i in range(len(gal_type)):
	if gal_type[i]==1:
		one=one+1
	elif gal_type[i]==2:
		two=two+1
	elif gal_type[i]==3:
		three=three+1
	elif gal_type[i]==4:
		four=four+1

print("1:",one)
print("2:",two)
print("3:",three)
print("4:",four)
