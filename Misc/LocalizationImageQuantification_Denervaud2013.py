#%%

# The following is to test how are pipeline works
# with images from other groups. First, we will
# test with Denervaud microfluidics screen

# IMPORT ALL PACKAGES REQUIRED
import pandas as pd
import numpy as np
import scipy as scipy
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.io as scipyio
from PIL import Image
import PIL
import glob
import os
import progressbar
import csv
import statistics
import time



# %%
# Set directory to location of fluorescence images (these are the folders
# that contain both dia and epi image headers)

image_path = "/Volumes/SeagatePortabl/DenervaudImages"
os.chdir(image_path)

# The following generates a dataframe that contains the image that needs
# to be opened and analyzed, the timepoint, the imageID (which includes
# well position, and the )

def grepChannel(x):
    # search string for dia or epi
    if x.find('dia') > 0:
        channel = 'dia'
    elif x.find('epi') > 0:
        channel = 'epi'
    return(channel)


def grepProtein(x):
    # with the appropriate channel, identify protein
    start = x.find('DenervaudImages') + len('DenervaudImages') + 1
    strSubset=x[start:len(x)]
    end = strSubset.find('/')
    return(strSubset[0:end])


def grepImageID(x):
    # with the appropriate channel, identify protein
    start = x.find('_pos') + len('_pos')
    end = x.find('.tif')
    return(x[start:end])

def grepTime(x):
    # with the appropriate channel, identify protein
    start = x.find('_t') + len('_t')
    end = x.find('.tif')
    return(x[start:end])



count = 0
imgIndex = []
for root, dirs, files in os.walk(os.getcwd()):
    for name in files:
        if name.endswith('.tif'):
            imgIndex.append({'Path': os.path.join(root, name)})
            count = count + 1
            print(count, end="\r")

imgIndex = pd.DataFrame(imgIndex)
imgIndex['Protein'] = pd.Series(imgIndex.iloc[:,0]).apply(grepProtein)
imgIndex['ImageID'] = pd.Series(imgIndex.iloc[:,0]).apply(grepImageID)
imgIndex['Channel'] = pd.Series(imgIndex.iloc[:,0]).apply(grepChannel)
imgIndex['Time'] = pd.Series(imgIndex.iloc[:,0]).apply(grepTime)




# %%
# And now we repeat for the segmentation files
# Set directory to location of segmentation images

segmentation_path = "/Volumes/SeagatePortabl/DenervaudSegmentation"
os.chdir(segmentation_path)

def grepSegPath(x):
    # with the appropriate channel, identify protein
    start = 0
    end = x.find('/Total')
    return(x[start:end])

count = 0
segIndex = []
for root, dirs, files in os.walk(os.getcwd()):
    for name in files:
        if name.endswith('.tiff'):
            segIndex.append({'File': os.path.join(root, name)})
            count = count + 1
            print(count, end="\r")

segIndex = pd.DataFrame(segIndex)
segIndex['ImageID'] = pd.Series(segIndex.iloc[:,0]).apply(grepImageID)
segIndex['Path'] = pd.Series(segIndex.iloc[:,0]).apply(grepSegPath)








# %%

threshold = 98

pixvalTP_DF = pd.DataFrame([])
count = 0

# We just need the fluorescence channels ("epi")
imgIndexEPI = imgIndex[imgIndex['Channel']=='epi']

start_index = int(len(imgIndexEPI)/2)
end_index = len(imgIndexEPI)

# FOR LOOP START #
############################################################################################
print('function started')
for i in np.arange(start_index, end_index, 1):

    os.chdir(image_path)
    img_pxMAT = np.asmatrix(Image.open(imgIndexEPI.iloc[i,0]))
    
    if np.amax(img_pxMAT) == 0:
        continue

    protein = imgIndexEPI.iloc[i,:]['Protein']
    imageID = imgIndexEPI.iloc[i,:]['ImageID']
    TimePoint = imgIndexEPI.iloc[i,:]['Time']
    channel = imgIndexEPI.iloc[i,:]['Channel']

    FilteredObjects_Path = segIndex[segIndex['ImageID'] == imageID].iloc[0,:]['Path']

    if len(FilteredObjects_Path) == 0:
        continue
    
    os.chdir(FilteredObjects_Path)
    
    FilteredObjects = []
    for root, dirs, files in os.walk(os.getcwd()):
        for j in files:
            if j.endswith('tiff'):
                FilteredObjects.append({'ObjectFile': j})
            
    FilteredObjects = pd.DataFrame(FilteredObjects)
    
    if len(FilteredObjects) == 0:
        continue

    objMeas_DF = pd.DataFrame([])

    for k in np.arange(0, len(FilteredObjects), 1):

        object_measurements = []
        object_pxMAT = np.asmatrix(Image.open(FilteredObjects.iloc[k,0]))
        # object_pxMAT.setflags(write=1)
        # object_pxMAT[object_pxMAT == 255] = 1
        object_pxMAT = object_pxMAT/255

        try:
            rawGFP_object = np.multiply(object_pxMAT, img_pxMAT)
        except:
            continue

        rawGFP_object = np.array(rawGFP_object, dtype = np.uint8)

        rawGFP_object = rawGFP_object[rawGFP_object > 0] 

        factor = np.median(rawGFP_object)
        normGFP_object = rawGFP_object/factor
        
        normGFP_object = np.percentile(normGFP_object, threshold)
        
        x95th = np.median(normGFP_object)
        CellSize = len(rawGFP_object)
                    
        object_measurements = {'StrainID' : imageID,
                               'Protein' : protein,
                               'TimePoint' : TimePoint,
                               'ObjectNumber' : FilteredObjects.iloc[k,0],
                               'x95thPercentile' : [x95th],
                               'CellSize' : [CellSize],
                               'CellMedian' : [factor],
                               'Channel' : channel}

        object_measurements = pd.DataFrame(object_measurements)

        
        objMeas_DF = pd.concat([objMeas_DF, object_measurements])

    pixvalTP_DF = pd.concat([pixvalTP_DF, objMeas_DF])
    count = count + 1
    print(count/end_index)
print('function ended')



#%%
desktopPath = "/Users/brandonho/Desktop/"
os.chdir(desktopPath)

pixvalTP_DF.to_csv("denervaud_part3.csv")


#%%
# Since the above data required separate parts since there are
# too many cells, we will load the data here, and concat the
# dataframes

data1 = pd.read_csv("/Users/brandonho/Dropbox (Grant Brown's Lab)/Complete Data (Brandon Ho)/Year 7/LiteratureScreenImages/Denervaud2013_FillHitList/20211213_denervaud_part1.csv")
data2 = pd.read_csv("/Users/brandonho/Dropbox (Grant Brown's Lab)/Complete Data (Brandon Ho)/Year 7/LiteratureScreenImages/Denervaud2013_FillHitList/20211213_denervaud_part2.csv")
data3 = pd.read_csv("/Users/brandonho/Dropbox (Grant Brown's Lab)/Complete Data (Brandon Ho)/Year 7/LiteratureScreenImages/Denervaud2013_FillHitList/20211213_denervaud_part3.csv")

data = pd.concat([data1, data2, data3])
data = data[['StrainID', 'Protein', 'TimePoint', 'x95thPercentile',
             'CellSize', 'CellMedian']]


#%%
# Will need to convert the timepoints into actual numbers

def convertTimeFrame(dataframe):
    allProteins = pd.DataFrame([])

    for i in np.unique(dataframe['Protein']):

        # Subset for the protein of interest
        proteinSubset = dataframe[dataframe['Protein']==i]

        # Then assign time frame number to the time written out as strings
        time = pd.DataFrame(np.unique(proteinSubset['TimePoint']))
        time.columns = ['TimePoint']
        time = time.sort_values(by='TimePoint')
        time['TimeFrame'] = np.arange(0, len(time), 1)

        # merge with the quantification dataframe, and return for function output
        proteinSubset = pd.merge(proteinSubset, time, on='TimePoint')
        allProteins = pd.concat([allProteins, proteinSubset])
    
    return(allProteins)

quantifiedCells = convertTimeFrame(data)


#%%
# Next, all of Denervaud's files are in the positions
# of their microfluidic chip, I'm assuming. I have generated
# a map which contains the ID of the locations and which these
# correspond to, protein wise

denervaudMap = pd.read_csv("/Users/brandonho/Dropbox (Grant Brown's Lab)/Complete Data (Brandon Ho)/Year 7/LiteratureScreenImages/Denervaud2013_FillHitList/DenervaudHits.csv")
denervaudMapSubset = denervaudMap[['Denervaud_HighMMS_Hits', 'FileName']]
denervaudMapSubset.columns = ['Gene', 'FileName']

finalDF = pd.merge(quantifiedCells, denervaudMapSubset,
                   left_on='Protein', right_on='FileName')




# %%
# Next, determine the threshold from untreated cells to define a
# cell that has a localization change event

def calcThreshold(protein,data):

    # First subset all the quantified cells for just the untreated timepoint
    untreatedSubset = data[data['TimeFrame']==0]

    # determine the statistical measures for protein specified
    proteinSubset = untreatedSubset[untreatedSubset['Gene']==protein]
    medianValue = np.mean(proteinSubset['x95thPercentile'])
    MAD = np.std(proteinSubset['x95thPercentile'])

    return(medianValue, MAD)

calcThreshold(protein='Lcd1', data=finalDF)



# %%
# Now, assign localization to each cell, whether there is a change or not

proteinIdentity = []
increasedCells = []
decreasedCells = []
localizationAssignment = []
timeFrame = []
totalCells = []

count = 0
for i in np.unique(finalDF['Gene']):

    # First, subset the dataframe for the protein of interest
    proteinSubset = finalDF[finalDF['Gene']==i]
    
    for j in np.unique(proteinSubset['TimeFrame']):
        count = count + 1

        # the for each timepoint, calculate number of cells that exceed threshold
        # for increases or decreases in localization
        timeSubset = proteinSubset[proteinSubset['TimeFrame']==j]

        # derive the statistical measures from previous function
        median = calcThreshold(protein=i, data=finalDF)[0]
        mad = calcThreshold(protein=i, data=finalDF)[1]
        increasedThreshold = median + 1.5*mad
        decreasedThreshold = median - 1.5*mad

        # now determine number of cells that exceed or decrease threshold
        increased = len(timeSubset[timeSubset['x95thPercentile']>increasedThreshold])
        decreased = len(timeSubset[timeSubset['x95thPercentile']<decreasedThreshold])
        total = len(timeSubset)

        proteinIdentity.append(i)
        timeFrame.append(j)
        increasedCells.append(increased)
        decreasedCells.append(decreased)
        totalCells.append(total)

        if increasedCells > decreasedCells:
            localizationAssignment.append('Increased')
        elif increasedCells < decreasedCells:
            localizationAssignment.append('Increased')
        else:
            localizationAssignment.append('NA')
    print(count)




# %%
desktopPath = "/Users/brandonho/Desktop"
os.chdir(desktopPath)

finalLocalization = {
    'Protein' : proteinIdentity,
    'Time' : timeFrame,
    'Increased' : increasedCells,
    'Decreased' : decreasedCells,
    'TotalCells' : totalCells,
    'LocalizationCall' : localizationAssignment
}

finalLocalization = pd.DataFrame(finalLocalization)
finalLocalization['IncreasedCellsPercentage']=finalLocalization['Increased']/finalLocalization['TotalCells']
finalLocalization['DecreasedCellsPercentage']=finalLocalization['Decreased']/finalLocalization['TotalCells']

# sns.lineplot(x='Time', y='IncreasedCellsPercentage', data=finalLocalization, hue='Protein')
# plt.savefig('DenervaudApplication.pdf')


# %%

DEN_matrix = finalLocalization.pivot(index='Protein', columns='Time',
                                     values='IncreasedCellsPercentage')







# %%

dataSubset = finalLocalization[finalLocalization['Protein']=='Dhh1']
sns.scatterplot(x='Time', y='IncreasedCellsPercentage', data=dataSubset)

plt.savefig("Dhh1.pdf")

# %%
