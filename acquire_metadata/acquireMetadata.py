#%%
# IMPORT ALL PACKAGES REQUIRED
import pandas as pd
import numpy as np
import os



#%%
# The following generates a dataframe that contains the image that needs
# to be opened and analyzed, the timepoint, the imageID (which includes
# well position, etc.) - this acts as a directory to reference and know
# which files to open and close

def stringInterest(x):
    start = x.find('CYCLoPs_')
    end = x.find('000')
    return(x[start:end])

def grepProtein(x):
    start = stringInterest(x).find('_') + 1
    end = stringInterest(x)[start:].find('_')
    return(stringInterest(x)[start:start+end])

def grepTreatment(x):
    start = stringInterest(x).find(grepProtein(x)) + len(grepProtein(x)) + 1
    end = stringInterest(x)[start:].find('_')
    return(stringInterest(x)[start:start+end])

def folder(x):
    return(x[0:x.rindex("/")+1])

# For the Chong et al data, an additonal step needs to be done. There are duplicated
# images, and so these will be removed from the analysis

def findSeries(x):
    if x.find(").flex") > 0:
        return("flag")
    else:
        return("pass")


# Now determine the image index
def getImgIndex(imageLocation):

    os.chdir(imageLocation)

    count = 0
    imgIndex = []
    for root, dirs, files in os.walk(os.getcwd()):
        for name in files:
            if name.endswith(".tiff"):
                imgIndex.append({'Path': os.path.join(root, name)})
                count = count + 1
                print(count, end="\r")

    imgIndex = pd.DataFrame(imgIndex)
    imgIndex['Folder'] = pd.Series(imgIndex.iloc[:,0]).apply(folder)
    imgIndex['Treatment'] = pd.Series(imgIndex.iloc[:,0]).apply(grepTreatment)
    imgIndex['Protein'] = pd.Series(imgIndex.iloc[:,0]).apply(grepProtein)
    imgIndex['Series'] = pd.Series(imgIndex.iloc[:,0]).apply(findSeries)

    imgIndex = imgIndex[imgIndex['Series'] == 'pass']

    return(imgIndex)