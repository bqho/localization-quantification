# %%
# Import modules required for analysis
import pandas as pd
import numpy as np
from PIL import Image
import PIL
import glob
import os


# Import custom modules
import AcquireMetadata.acquireMetadata as meta



#%%
# associate protein and treatment for each of the images
# using function from metadata acquisition module

imageIndex = meta.getImgIndex("/Volumes/Seagate Portabl/ChongAnalysis_OutputFiles")


# Filter out the segmented single cells from the original GFP full image

def imgOrSeg(x):
    if x.find("TotalCell") > 0:
        return("segmentation")
    else:
        return("image")

imageIndex['ImageType'] = pd.Series(imageIndex.iloc[:,0]).apply(imgOrSeg)

imgIndex = imageIndex[imageIndex['ImageType'] == 'image']
segIndex = imageIndex[imageIndex['ImageType'] == 'segmentation']






# %%
# Then run the main function of the localization analysis pipeline.
# Open image and segmentation, and determine the intensity distributions

threshold = 99

pixvalTP_DF = pd.DataFrame([])
start = 1000
end = len(imgIndex)

# count = 0


# FOR LOOP START #
############################################################################################
print('function started')
for i in np.arange(0, len(imgIndex), 1):

    img_pxMAT = np.asmatrix(Image.open(imgIndex.iloc[i,0]))
    
    if np.amax(img_pxMAT) == 0:
        continue

    Protein = imgIndex.iloc[i,:]['Protein']
    Treatment = imgIndex.iloc[i,:]['Treatment']

    # The segmentation images are in the same folder as the imgIndex
    segPath = imgIndex.iloc[i,:]['Folder']
    segImages = segIndex[segIndex['Folder'] == segPath]

    if len(segImages) == 0:
        continue

    objMeas_DF = pd.DataFrame([])

    for k in np.arange(0, len(segImages), 1):

        object_measurements = []
        
        object_pxMAT = np.asmatrix(Image.open(segImages.iloc[k,0]))
        object_pxMAT_modify = object_pxMAT.copy()
        object_pxMAT_modify.setflags(write=1)
        object_pxMAT_modify[object_pxMAT_modify == 255] = 1

        rawGFP_object = np.multiply(object_pxMAT_modify, img_pxMAT)
        rawGFP_object = np.array(rawGFP_object, dtype = np.uint16)

        rawGFP_object = rawGFP_object[rawGFP_object > 0] 

        if len(rawGFP_object) == 0:
            continue

        factor = np.median(rawGFP_object)
        normGFP_object = rawGFP_object/factor
        
        normGFP_object = np.percentile(normGFP_object, threshold)
        #normGFP_object = np.percentile(normGFP_object, threshold) - np.percentile(normGFP_object, (100-threshold))
        
        x95th = np.median(normGFP_object)
        CellSize = len(rawGFP_object)
                    
        object_measurements = {'Protein' : Protein,
                               'Treatment' : Treatment,
                               'ObjectNumber' : segImages.iloc[k,0],
                               'x95thPercentile' : [x95th],
                               'CellSize' : [CellSize],
                               'CellMedian' : [factor]}

        object_measurements = pd.DataFrame(object_measurements)
        object_measurements.columns = ['Protein', 'Treatment', 
                                       'ObjectNumber',
                                       '95thPercentile',
                                       'CellSize',
                                       'CellMedian']
        
        objMeas_DF = pd.concat([objMeas_DF, object_measurements])

    pixvalTP_DF = pd.concat([pixvalTP_DF, objMeas_DF])
    # count = count + 1
    # print(str(count) + " out of " + str(len(imgIndex)))
print('function ended')



desktop_path = "/Users/brandonho/Desktop/"
os.chdir(desktop_path)

pixvalTP_DF.to_csv('20211206_ChongHitsAnalyzed_ProperAnalysis_99Percentile_Part3.csv')


# %%
