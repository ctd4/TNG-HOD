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


def best_fit_slope(xs,ys):
	m = (((np.mean(xs)*np.mean(ys)) - np.mean(xs*ys))/((np.mean(xs)**2) - np.mean(xs**2)))
	return m

baseUrl = 'http://www.tng-project.org/api/'
headers = {"api-key":"c9bf2fda1907189085c080d503cdb0d4"}

def get(path, params=None):
    # make HTTP GET request to path
    r = requests.get(path, params=params, headers=headers)

    # raise exception if response code is not HTTP SUCCESS (200)
    r.raise_for_status()

    if r.headers['content-type'] == 'application/json':
        return r.json() # parse json responses automatically

    if 'content-disposition' in r.headers:
        filename = r.headers['content-disposition'].split("filename=")[1]
        with open(filename, 'wb') as f:
            f.write(r.content)
        return filename # return the filename string

    return r

r = get(baseUrl)

#Functions
#===================================================================================
def gcPath(basePath, snapNum, chunkNum=0):
    """ Return absolute path to a group catalog HDF5 file (modify as needed). """
    gcPath = basePath + '/groups_%03d/' % snapNum
    filePath1 = gcPath + 'groups_%03d.%d.hdf5' % (snapNum, chunkNum)
    filePath2 = gcPath + 'fof_subhalo_tab_%03d.%d.hdf5' % (snapNum, chunkNum)

    if isfile(filePath1):
        return filePath1
    return filePath2


def offsetPath(basePath, snapNum):
    """ Return absolute path to a separate offset file (modify as needed). """
    offsetPath = basePath + '/../postprocessing/offsets/offsets_%03d.hdf5' % snapNum

    return offsetPath


#Groupcat functions
def loadObjects(basePath, snapNum, gName, nName, fields):
    """ Load either halo or subhalo information from the group catalog. """
    result = {}

    # make sure fields is not a single element
    if isinstance(fields, six.string_types):
        fields = [fields]

    # load header from first chunk
    with h5py.File(gcPath(basePath, snapNum), 'r') as f:

        header = dict(f['Header'].attrs.items())
        result['count'] = f['Header'].attrs['N' + nName + '_Total']

        if not result['count']:
            print('warning: zero groups, empty return (snap=' + str(snapNum) + ').')
            return result

        # if fields not specified, load everything
        if not fields:
            fields = list(f[gName].keys())

        for field in fields:
            # verify existence
            if field not in f[gName].keys():
                raise Exception("Group catalog does not have requested field [" + field + "]!")

            # replace local length with global
            shape = list(f[gName][field].shape)
            shape[0] = result['count']

            # allocate within return dict
            result[field] = np.zeros(shape, dtype=f[gName][field].dtype)

    # loop over chunks
    wOffset = 0

    for i in range(header['NumFiles']):
        f = h5py.File(gcPath(basePath, snapNum, i), 'r')

        if not f['Header'].attrs['N'+nName+'_ThisFile']:
            continue  # empty file chunk

        # loop over each requested field
        for field in fields:
            if field not in f[gName].keys():
                raise Exception("Group catalog does not have requested field [" + field + "]!")

            # shape and type
            shape = f[gName][field].shape

            # read data local to the current file
            if len(shape) == 1:
                result[field][wOffset:wOffset+shape[0]] = f[gName][field][0:shape[0]]
            else:
                result[field][wOffset:wOffset+shape[0], :] = f[gName][field][0:shape[0], :]

        wOffset += shape[0]
        f.close()

    # only a single field? then return the array instead of a single item dict
    if len(fields) == 1:
        return result[fields[0]]

    return result




def loadSubhalos(basePath, snapNum, fields=None):
    """ Load all subhalo information from the entire group catalog for one snapshot
       (optionally restrict to a subset given by fields). """

    return loadObjects(basePath, snapNum, "Subhalo", "subgroups", fields)


def loadHalos(basePath, snapNum, fields=None):
    """ Load all halo information from the entire group catalog for one snapshot
       (optionally restrict to a subset given by fields). """

    return loadObjects(basePath, snapNum, "Group", "groups", fields)

basePath="/home/callum/Documents/SummerProject/z=0"
snapNum=99

#Main Program
#===================================================================
Sub_effective_rad=np.array(loadSubhalos(basePath,snapNum,fields="SubhaloHalfmassRadType"))[:,4]
Subhalo_mass_star=np.array(loadSubhalos(basePath,snapNum,fields="SubhaloMassType"))[:,4]* 1e10 / 0.6774
SubID=np.arange(Subhalo_mass_star.shape[0])
mask=(Subhalo_mass_star>1e9)

Sub_effective_rad=Sub_effective_rad[mask]
SubID=SubID[mask]


for n in range(0,Sub_effective_rad.shape[0]):
	print(n)
	sub_prog_url = "http://www.tng-project.org/api/TNG100-1/snapshots/99/subhalos/"+str(SubID[n])+"/"
	sub_prog = get(sub_prog_url)
	#print(sub_prog['pos_x'], sub_prog['pos_y'],sub_prog['pos_z'])

	cutout_request = {'stars':'Coordinates,GFM_Metallicity'}				#Change GFM_Metallicity to GFM_StellarFormationTime for age calculations	
	cutout = get(sub_prog_url+"cutout.hdf5", cutout_request)




	#Extract fields from snapshot for that subhalo
	with h5py.File(cutout,'r') as f:
		x = f['PartType4']['Coordinates'][:,0] - sub_prog['pos_x']
		y = f['PartType4']['Coordinates'][:,1] - sub_prog['pos_y']
		z = f['PartType4']['Coordinates'][:,2] - sub_prog['pos_z']
		gasmet=f['PartType4']['GFM_Metallicity'][:]
	distances=(x**2+y**2+z**2)**0.5
	
	#gasmet=1-np.array(gasmet)												#used for age calculation by taking time when formed away from 1.
	
	
	#Define bins and calculate means
	bins=np.arange(0,1.5*Sub_effective_rad[n],(1.5*Sub_effective_rad[n])/10)
	metallicity_means=[]
	distance_bin=[]
	met_grad=[]
	for i in range(1,len(bins)):

		mask=(distances>bins[i-1]) & (distances<bins[i])




		if distances[mask].shape[0]==0:
			pass
		else:
			distance_bin.append(np.amax(distances[mask])-(1.5*Sub_effective_rad[n])/20)
			metallicity_means.append(np.mean(gasmet[mask]))

	
	with open("met_grad2.csv",'a') as f:
		writer = csv.writer(f)
		writer.writerow([str(np.polyfit(np.array(distance_bin),np.array(metallicity_means),1)[0])])
		



	




#Plot
#===========================
plt.figure(0)
plt.plot(distance_bin,metallicity_means,'-')
plt.xlabel("Subhalo Radii [ckpc/h]")
plt.ylabel("Mean Metallicity")
plt.title("Gas Metallicity Gradient")



