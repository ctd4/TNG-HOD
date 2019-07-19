import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
import six
from os.path import isfile
import h5py
np.set_printoptions(threshold=17000)

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


#HOD plotting function that includes centrals and satellites
def plotHOD(halo_mass_limit,stellar_mass_limit,No_bins,group_mass,ID,GroupFirstSub,SubhaloGrNr,Subhalo_mass_star):
	N=[]
	N_central=[]
	N_satellite=[]

	#Limit by Halo mass
	length_mask=(group_mass>halo_mass_limit)
	group_mass=group_mass[length_mask]
	ID=ID[length_mask]
	#GroupFirstSub=GroupFirstSub[length_mask]

	length=group_mass.shape[0]
	print("length:",length)



	#Count number of subhalos in each halo limiting by stellar mass
	for i in range(length):
		mask1=(SubhaloGrNr==ID[i])
		mask2=(Subhalo_mass_star>stellar_mass_limit)
		mask_star=mask1 & mask2



		#print(Subhalo_mass_star[mask_star].shape[0],"Halo number:",i)
		N.append(Subhalo_mass_star[mask_star].shape[0])
		if N[i]<1:
			N_central.append(0)
			N_satellite.append(0)
		else:
			N_central.append(1)
			N_satellite.append(N[i]-1)

	binnedN=stats.binned_statistic(np.log10(group_mass),N,statistic='mean',bins=No_bins)
	binN_edge=binnedN[1]

	binnedCentral=stats.binned_statistic(np.log10(group_mass),N_central,statistic='mean',bins=No_bins)
	binCen_edge=binnedCentral[1]

	binnedSatellite=stats.binned_statistic(np.log10(group_mass),N_satellite,statistic='mean',bins=No_bins)
	binSat_edge=binnedSatellite[1]

	plt.figure(0)
	plt.plot(binCen_edge[:-1]+(binCen_edge[2]-binCen_edge[1])/2,binnedCentral[0],'.',label="Central",color="darkorange")
	plt.plot(binSat_edge[:-1]+(binSat_edge[2]-binSat_edge[1])/2,binnedSatellite[0],'.',label="Satellite",color="dodgerblue")
	plt.plot(binN_edge[:-1]+(binN_edge[2]-binN_edge[1])/2,binnedN[0],'.',label="Total",color="g")
	#plt.xscale('log')
	plt.yscale('log')
	plt.xlabel('Log[Mass/[$M_\odot$]]')
	plt.ylabel('Mean N')
	plt.title('Halo Occupation Distribution')
	plt.legend()
	plt.show()
	return

#Function that plots HOD, splitting on a parameter and plotting satellites and centrals
def plotHODsplitSatellite(halo_mass_limit,stellar_mass_limit,No_bins,group_mass,ID,GroupFirstSub,SubhaloGrNr,Subhalo_mass_star,SplitParameter):
	N=[]
	N1=[]
	N2=[]
	ID1=[]
	ID2=[]
	
	N_central=[]
	N_satellite=[]
	N1_central=[]
	N1_satellite=[]
	N2_central=[]
	N2_satellite=[]

	i1=0
	i2=0


	#Limit by Halo mass
	length_mask=(group_mass>halo_mass_limit)
	group_mass=group_mass[length_mask]
	ID=ID[length_mask]

	length=group_mass.shape[0]
	print("length:",length)
	



	#Count number of subhalos in each halo
	for i in range(length):
		mask1=(SubhaloGrNr==ID[i])
		mask2=(Subhalo_mass_star>stellar_mass_limit)
		
		mask_total=mask1 & mask2



		#print(i)
		N.append(Subhalo_mass_star[mask_total].shape[0])
		if N[i]<1:
			N_central.append(0)
			N_satellite.append(0)
		else:
			N_central.append(1)
			N_satellite.append(N[i]-1)		
		
		#Two populations HOD is split into 
		if SplitParameter[GroupFirstSub[ID[i]]]>np.percentile(SplitParameter[mask2],75) and Subhalo_mass_star[GroupFirstSub[ID[i]]]>stellar_mass_limit:
			ID1.append(i)

			N1.append(Subhalo_mass_star[mask_total].shape[0])
			if N1[i1]<1:
				N1_central.append(0)
				N1_satellite.append(0)
			else:
				N1_central.append(1)
				N1_satellite.append(N1[i1]-1)
			i1=i1+1
			
						
		elif SplitParameter[GroupFirstSub[ID[i]]]<np.percentile(SplitParameter[mask2],25) and Subhalo_mass_star[GroupFirstSub[ID[i]]]>stellar_mass_limit:
			ID2.append(i)

			N2.append(Subhalo_mass_star[mask_total].shape[0])
			if N2[i2]<1:
				N2_central.append(0)
				N2_satellite.append(0)
			else:
				N2_central.append(1)
				N2_satellite.append(N2[i2]-1)
			i2=i2+1			


	print("Median mass:",np.median(Subhalo_mass_star[mask2]))
	print("Mean mass:",np.mean(Subhalo_mass_star[mask2]))
	print("=====================")
	print("N shape:",len(N))
	print("N1 shape:",len(N1))
	print("N2 shape:",len(N2))



	
	binnedN=stats.binned_statistic(np.log10(group_mass),N,statistic='mean',bins=No_bins)
	binN_edge=binnedN[1]
	binnedN1=stats.binned_statistic(np.log10(group_mass[ID1]),N1,statistic='mean',bins=No_bins)
	binN1_edge=binnedN1[1]
	binnedN2=stats.binned_statistic(np.log10(group_mass[ID2]),N2,statistic='mean',bins=No_bins)
	binN2_edge=binnedN2[1]
	
	binnedCentral=stats.binned_statistic(np.log10(group_mass),N_central,statistic='mean',bins=No_bins)
	binCen_edge=binnedCentral[1]

	binnedSatellite=stats.binned_statistic(np.log10(group_mass),N_satellite,statistic='mean',bins=No_bins)
	binSat_edge=binnedSatellite[1]
	
	binnedCentral1=stats.binned_statistic(np.log10(group_mass[ID1]),N1_central,statistic='mean',bins=No_bins)
	binCen_edge1=binnedCentral1[1]

	binnedSatellite1=stats.binned_statistic(np.log10(group_mass[ID1]),N1_satellite,statistic='mean',bins=No_bins)
	binSat_edge1=binnedSatellite1[1]
	
	binnedCentral2=stats.binned_statistic(np.log10(group_mass[ID2]),N2_central,statistic='mean',bins=No_bins)
	binCen_edge2=binnedCentral2[1]

	binnedSatellite2=stats.binned_statistic(np.log10(group_mass[ID2]),N2_satellite,statistic='mean',bins=No_bins)
	binSat_edge2=binnedSatellite2[1]

	plt.figure(0)
	plt.plot(binN_edge[:-1]+(binN_edge[2]-binN_edge[1])/2,binnedN[0],'-',label="Total",color="black")
	plt.plot(binN1_edge[:-1]+(binN1_edge[2]-binN1_edge[1])/2,binnedN1[0],'-',label="Largest 25%",color="g")
	plt.plot(binN2_edge[:-1]+(binN2_edge[2]-binN2_edge[1])/2,binnedN2[0],'-',label="Smallest 25%",color="r")
	
	plt.plot(binCen_edge[:-1]+(binCen_edge[2]-binCen_edge[1])/2,binnedCentral[0],':',label="Central",color="black")
	plt.plot(binSat_edge[:-1]+(binSat_edge[2]-binSat_edge[1])/2,binnedSatellite[0],'--',label="Satellite",color="black")
	
	plt.plot(binCen_edge1[:-1]+(binCen_edge1[2]-binCen_edge1[1])/2,binnedCentral1[0],':',color="g")
	plt.plot(binSat_edge1[:-1]+(binSat_edge1[2]-binSat_edge1[1])/2,binnedSatellite1[0],'--',color="g")
	
	plt.plot(binCen_edge2[:-1]+(binCen_edge2[2]-binCen_edge2[1])/2,binnedCentral2[0],':',color="r")
	plt.plot(binSat_edge2[:-1]+(binSat_edge2[2]-binSat_edge2[1])/2,binnedSatellite2[0],'--',color="r")
	
	plt.yscale('log')
	plt.xlabel('Log[Mass/[$M_\odot$]]')
	plt.ylabel('Mean N')
	plt.title('Halo Occupation Distribution')
	plt.legend()
	plt.show()
	return

 
#Extracting Datasets from TNG data
#=================================================================================
SubhaloPos=np.array(loadSubhalos(basePath,snapNum,fields="SubhaloPos"))
HaloPos=np.array(loadHalos(basePath,snapNum,fields="GroupPos"))


Subhalo_mass=np.array(loadSubhalos(basePath,snapNum,fields="SubhaloMass"))
Subhalo_mass_star=np.array(loadSubhalos(basePath,snapNum,fields="SubhaloMassType"))[:,4]* 1e10 / 0.6774
SubhaloGrNr=np.array(loadSubhalos(basePath,snapNum,fields="SubhaloGrNr"))
SubhaloSFR=np.array(loadSubhalos(basePath,snapNum,fields="SubhaloSFR"))
SubhaloGasMetallicity=np.array(loadSubhalos(basePath,snapNum,fields="SubhaloGasMetallicity"))
SubhaloSpin=np.array(loadSubhalos(basePath,snapNum,fields="SubhaloSpin"))[:,2]
SubhaloBHMdot=np.array(loadSubhalos(basePath,snapNum,fields="SubhaloBHMdot"))
SubhaloStarMetallicity=np.array(loadSubhalos(basePath,snapNum,fields="SubhaloStarMetallicity"))

group_mass=np.array(loadHalos(basePath,snapNum,fields="GroupMass"))* 1e10 / 0.6774
GroupFirstSub=np.array(loadHalos(basePath,snapNum,fields="GroupFirstSub"))
ID=np.arange(group_mass.shape[0])



print("Subhalo mass star:",Subhalo_mass_star)
print("Subhalo mass star shape:",Subhalo_mass_star.shape)






#Plots
#=================================================================================
	
#plotHOD(1e11,1e9,30,group_mass,ID,GroupFirstSub,SubhaloGrNr,Subhalo_mass_star)
plotHODsplitSatellite(1e11,1e9,30,group_mass,ID,GroupFirstSub,SubhaloGrNr,Subhalo_mass_star,SubhaloSFR)















