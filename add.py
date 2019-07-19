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

np.set_printoptions(threshold=170000000)
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


#================================================================
Subhalo_mass_star=np.array(loadSubhalos(basePath,snapNum,fields="SubhaloMassType"))[:,4]* 1e10 / 0.6774
SubhaloStarMetallicity=np.array(loadSubhalos(basePath,snapNum,fields="SubhaloStarMetallicity"))
ID=np.arange(Subhalo_mass_star.shape[0])
mask=(Subhalo_mass_star>1e9)
SubhaloStarMetallicity=SubhaloStarMetallicity[mask]
print("Parameter shape:",SubhaloStarMetallicity.shape)


SubhaloMetGrad=np.loadtxt("met_grad2.csv")
SubhaloAgeGrad=np.loadtxt("age_grad4.csv")
gal_type=np.loadtxt("gal_type.csv")
print("Met grad shape:",SubhaloMetGrad.shape)
print("Age grad shape:",SubhaloAgeGrad.shape)

met_grad=[None]*Subhalo_mass_star.shape[0]
print("met_grad shape:",len(met_grad))

j=0
for i in range(0,ID.shape[0]):
	#print(i)
	if mask[i]==False:
		met_grad[i]=0
	elif mask[i]==True:
		met_grad[i]=gal_type[j]			#Change variable depending on which data is being used.
		j=j+1

print("========================")
print(len(met_grad))
np.savetxt("gal_type_large.csv",met_grad,delimiter=',')



















