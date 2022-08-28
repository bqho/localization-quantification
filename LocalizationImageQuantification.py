# %%
# Import modules required for analysis_________________________________________

# Import modules
import pandas as pd
import numpy as np
from PIL import Image
import PIL
import glob
import os

# Import custom modules
import acquire_metadata.acquireMetadata as meta



#%%
# Run functions to obtain metadata, associating protein identity, time of
# drug treatments, etc. for each image ________________________________________

# Apply function, enter path to images
image_path = "/Volumes/Seagate Portabl/ChongAnalysis_OutputFiles" #image path
imageIndex = meta.getImgIndex(image_path)


# Folder contains both the original raw images from the OPERA microscope,
# as well as the segmented images after running the CellProfiler pipeline.
# Following function designates if the image is raw image or segmented image
def imgOrSeg(x):
    if x.find("TotalCell") > 0:
        return("segmentation")
    else:
        return("image")

# Apply function, and generate dataframe with file meta data and image paths
imageIndex['ImageType'] = pd.Series(imageIndex.iloc[:,0]).apply(imgOrSeg)
imgIndex = imageIndex[imageIndex['ImageType'] == 'image']
segIndex = imageIndex[imageIndex['ImageType'] == 'segmentation']



# %%
# Run main function of the localization analysis pipeline.
# Opens image and segmentation, and determine the intensity distributions


def LOC_quant(imgIndex, segIndex, threshold):
    """Main localization quantification function. Takes the list of
    files in the image and segmentation dataframes generated above, and
    uses it as a guide to open the appropriate image and segmentation files,
    then calculating the pixel intensity distributions.

    Args:
        imgIndex (dataframe): Dataframe containing the file path for all image files
        segIndex (dataframe): Dataframe containing the file path for all segmentation
        files
        threshold (int): A value between 0-100 that represents the percentile
        of the intensity distribution to calculate per cell (the LOC score)
    """
    pixvalTP_DF = pd.DataFrame([])
    start = 0
    end = len(imgIndex)

    count = 0

    # A for loop to go through each image file, find the corresponding
    # segmentation file, and extract the single cell intensity distribution
    
    print('function started')

    for i in np.arange(0, len(imgIndex), 1):

        # Opens images and reads as numeric matrix
        image_matrix = np.asmatrix(Image.open(imgIndex.iloc[i,0]))
        
        if np.amax(image_matrix) == 0:
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
            
            # Pixels corresponding to identified cells by CellProfiler in the
            # segmentation image are given a value of 1, and outside of the cell
            # is given a value of 0
            cell_object = np.asmatrix(Image.open(segImages.iloc[k,0]))
            cell_object_matrix = cell_object.copy()
            cell_object_matrix.setflags(write=1)
            cell_object_matrix[cell_object_matrix == 255] = 1

            # Multipling the cell_object_matrix with the raw image_matrix results 
            # in all pixels outside of the cell being 0, and inside the cell unchanged
            raw_GFP_matrix = np.multiply(cell_object_matrix, image_matrix)
            raw_GFP_matrix = np.array(raw_GFP_matrix, dtype = np.uint16)

            # extract only the pixel intensities greater than 0, corresponding
            # to the intensities inside the cell
            raw_GFP_matrix = raw_GFP_matrix[raw_GFP_matrix > 0] 

            if len(raw_GFP_matrix) == 0:
                continue

            # Normalization: divide by meadian to account for abundance changes
            factor = np.median(raw_GFP_matrix)
            normalized_GFP_object = raw_GFP_matrix/factor
            
            # Most important line: determining the LOC score! This gets the value of the
            # designated percentile of the cell intensity distribution
            normalized_GFP_object = np.percentile(normalized_GFP_object, threshold)
            LOC_score = np.median(normalized_GFP_object)

            # Obtain the size of the cell, in # of pixels
            cell_size = len(raw_GFP_matrix)
                        
            object_measurements = {'Protein' : Protein,
                                   'Treatment' : Treatment,
                                   'ObjectNumber' : segImages.iloc[k,0],
                                   'x95thPercentile' : [LOC_score],
                                   'CellSize' : [cell_size],
                                   'CellMedian' : [factor]}

            object_measurements = pd.DataFrame(object_measurements)
            object_measurements.columns = ['Protein',
                                           'Treatment', 
                                           'ObjectNumber',
                                           '95thPercentile',
                                           'CellSize',
                                           'CellMedian']
            
            objMeas_DF = pd.concat([objMeas_DF, object_measurements])

        pixvalTP_DF = pd.concat([pixvalTP_DF, objMeas_DF])

        # Grey next section out if you do not want to keep track of the
        # progress - printing to console will significantly slow down
        # analysis if you have a large number of files

        # --- Optional section --- #
        count = count + 1
        print(str(count) + " out of " + str(len(imgIndex)))
        # --- Optional section --- #

    print('function ended')
    return(pixvalTP_DF)



#%%
# Save the resulting quantification in user-defined directory
directory_path = "/Users/brandonho/Desktop/" # write out full path of directory
os.chdir(directory_path)

# Name your file
localization_quantification = LOC_quant(imgIndex, segIndex, threshold = 99)
localization_quantification.to_csv('enter_file_name.csv')


# %%
