#%%
import pandas as pd
import numpy as np
import scipy as scipy
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.io as scipyio
import os
from scipy import stats
from sklearn.cluster import KMeans
import re



# %%

data = pd.read_csv("/Users/brandonho/Dropbox (Grant Brown's Lab)/Complete Data (Brandon Ho)/Year 7/wildtypeLocalizationScreenAnalysis/20210809_HUMMS_penetrance.csv")



# %%
mms = data.iloc[:,0:12]

mmsTemp = mms.copy()
mmsTemp['MMS0'] = mmsTemp['MMS0']/mmsTemp['Max_MMS']
mmsTemp['MMS30'] = mmsTemp['MMS30']/mmsTemp['Max_MMS']
mmsTemp['MMS60'] = mmsTemp['MMS60']/mmsTemp['Max_MMS']
mmsTemp['MMS90'] = mmsTemp['MMS90']/mmsTemp['Max_MMS']
mmsTemp['MMS120'] = mmsTemp['MMS120']/mmsTemp['Max_MMS']
mmsTemp['MMS180'] = mmsTemp['MMS180']/mmsTemp['Max_MMS']
mmsTemp['MMS240'] = mmsTemp['MMS240']/mmsTemp['Max_MMS']

mmsTemp['MeanKinetics'] = mmsTemp.iloc[:, 1:8].sum(axis=1)



# %%
sns.catplot(x='MMS_pentile', y='Max_MMS', data=mmsTemp, kind='violin')
sns.catplot(x='MMS_pentile', y='MeanKinetics', data=mmsTemp, kind='violin')





# %%
hu = data.iloc[:,12:23]

huTemp = hu.copy()
huTemp['HU0'] = huTemp['HU0']/huTemp['Max_HU']
huTemp['HU30'] = huTemp['HU30']/huTemp['Max_HU']
huTemp['HU60'] = huTemp['HU60']/huTemp['Max_HU']
huTemp['HU90'] = huTemp['HU90']/huTemp['Max_HU']
huTemp['HU120'] = huTemp['HU120']/huTemp['Max_HU']
huTemp['HU180'] = huTemp['HU180']/huTemp['Max_HU']
huTemp['HU240'] = huTemp['HU240']/huTemp['Max_HU']

huTemp['MeanKinetics'] = huTemp.iloc[:, 1:8].sum(axis=1)


# %%
sns.catplot(x='HU_pentiles', y='MeanKinetics', data=huTemp, kind='violin')
sns.catplot(x='HU_pentiles', y='Max_HU', data=huTemp, kind='violin')


# %%
sns.scatterplot(y='Max_MMS', x='MeanKinetics', data=mmsTemp)
sns.scatterplot(y='Max_HU', x='MeanKinetics', data=huTemp)





# %%
# estimate the derivative for each point up until the maximum?

# First, determine the list of all values
test = hu.iloc[312,0:7].astype('float64')

# Next determine the max value, and also the index of the max value
maxValue = max(test)
maxIndex = test.argmax()

# subset the vector to obtain all the values before and at the time
# where the maximum penetrance is exhibited. The maximum value is
# the second last entry
valuesBeforeMax = test[0:maxIndex+2]

# Determine all differences between elements in the values
val = [x/maxValue - valuesBeforeMax[i - 1]/maxValue for i, x in enumerate(valuesBeforeMax)][1:]

# Determine all differences between elements of the time vector
timepoints = [0, 30, 60, 90, 120, 180, 240]
timepointsBeforeMax = timepoints[0:maxIndex+2]
time = [x - timepointsBeforeMax[i - 1] for i, x in enumerate(timepointsBeforeMax)][1:]

res = [i / j for i, j in zip(val, time)]





# %%
# ask for AUC?

def AUC(x):
    # First, determine the list of all values
    xVals = x.astype('float64')

    # Next determine the max value, and normalize the vector
    maxValue = max(xVals)
    xVals = xVals/maxValue

    # Determine the minimum normalized value to subtract
    minValue = min(xVals)

    # Determine all differences between elements in the values
    val = [((j-minValue) + (xVals[i - 1]-minValue))/2 for i, j in enumerate(xVals)][1:]

    # Determine all differences between elements of the time vector
    timepoints = [0, 30, 60, 90, 120, 180, 240]
    time = [j - timepoints[i - 1] for i, j in enumerate(timepoints)][1:]

    # multiply each element between the two vectors, and sum them to get AUC
    res = [i * j for i, j in zip(val, time)]
    return(sum(res))




# %%
# time at which the slope is at it's maximum? Or what the maximum slope is even

# First, determine the list of all values
test = hu.iloc[312,0:7].astype('float64')

# Next determine the max value, and also the index of the max value
maxValue = max(test)
maxIndex = test.argmax()

# Determine all differences between elements in the values
val = [x/maxValue - test[i - 1]/maxValue for i, x in enumerate(test)][1:]

# Determine all differences between elements of the time vector
timepoints = [0, 30, 60, 90, 120, 180, 240]
time = [x - timepoints[i - 1] for i, x in enumerate(timepoints)][1:]

res = [i / j for i, j in zip(val, time)]



#%%

def maxRate(x, returnKey):
    # First, determine the list of all values
    xVals = x.astype('float64')

    # Next determine the max value, and also the index of the max value
    maxValue = max(xVals)

    # Determine all differences between elements in the values
    val = [x/maxValue - xVals[i - 1]/maxValue for i, x in enumerate(xVals)][1:]

    # Determine all differences between elements of the time vector
    timepoints = [0, 30, 60, 90, 120, 180, 240]
    time = [x - timepoints[i - 1] for i, x in enumerate(timepoints)][1:]

    res = [i / j for i, j in zip(val, time)]

    if returnKey == 'slope':
        return(max(res))
    else:
        return(timepoints[res.index(max(res))+1])
    




# %%
# apply the maxRate function to the entire table for HU
columns = ['Gene', 'HU0', 'HU30', 'HU60', 'HU90', 'HU120', 'HU180', 'HU240',
           'HU_pentiles']
hu = data[columns]

hu['maxRate'] = hu.iloc[:,1:8].apply(maxRate, returnKey='slope', axis=1)
hu['timeOfMax'] = hu.iloc[:,1:8].apply(maxRate, returnKey='time', axis=1)

os.chdir("/Users/brandonho/Desktop/")

sns.catplot(x='HU_pentiles', y='maxRate', data=hu, kind='violin')
plt.savefig("hu_rates.pdf")

sns.catplot(x='HU_pentiles', y='timeOfMax', data=hu, kind='violin')
plt.savefig("hu_timeTo.pdf")




# %%
# apply the maxRate function to the entire table for MMS
columns = ['Gene', 'MMS0', 'MMS30', 'MMS60', 'MMS90', 'MMS120', 'MMS180', 'MMS240',
           'MMS_pentile']
mms = data[columns]

mms['maxRate'] = mms.iloc[:,1:8].apply(maxRate, returnKey='slope', axis=1)
mms['timeOfMax'] = mms.iloc[:,1:8].apply(maxRate, returnKey='time', axis=1)

os.chdir("/Users/brandonho/Desktop/")

sns.catplot(x='MMS_pentile', y='maxRate', data=mms)
# plt.savefig("mms_rates.pdf")

sns.catplot(x='MMS_pentile', y='timeOfMax', data=mms)
# plt.savefig("mms_timeto.pdf")




# %%
# Making violin plots for figure 5a
os.chdir("/Users/brandonho/Desktop/")

data = pd.read_csv("/Users/brandonho/Dropbox (Grant Brown's Lab)/Complete Data (Brandon Ho)/Year 7/wildtypeLocalizationScreenAnalysis/20210809_HUMMS_penetrance.csv")
dataSubset = data[['Gene', 'Max_MMS', 'Max_HU']]

sns.catplot(y='Max_MMS', data=dataSubset, kind='violin')

plt.savefig('MMSvioline.pdf')



# %%
os.chdir("/Users/brandonho/Desktop/")

data = pd.read_csv("/Users/brandonho/Dropbox (Grant Brown's Lab)/Complete Data (Brandon Ho)/Year 7/wildtypeLocalizationScreenAnalysis/20210809_HUMMS_penetrance.csv")
dataSubset = data[['Gene', 'Max_MMS', 'Max_HU']]

sns.catplot(y='Max_HU', data=dataSubset, kind='violin')

plt.savefig('HUvioline.pdf')
# %%
