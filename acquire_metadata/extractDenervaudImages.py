#%% 
import os
import pandas as pd
import subprocess
import numpy as np



#%%
# The format  Denervaud images are downloaded from the server 
# is .tar.bz, and these will be extracted into the designated 
# folder containing all our segmentation images

# provide the path to the .tar.bz files
directory = "/Volumes/Seagate Portabl/DenervaudImages"
os.chdir(directory)

# Extract appropriate meta data, which will be used to create
# the new folder locations and also to keep track of which 
# file is being extracted
zipFiles = []
count = 0

for root, dirs, files in os.walk(os.getcwd()):
    for name in files:

        if name.endswith('.DS_Store'):
            continue

        else:    
            zipFiles.append({'Path': os.path.join(root, name)})
            count = count + 1
            print(count, end="\r")

zipFiles = pd.DataFrame(zipFiles)

# Functions to extract the name of the folders, and also the
# the position identity from Denervaud annotations
def grepName(x):
    start = x.find('/DenervaudImages/') + len('/DenervaudImages/')
    end = x.find('.tar')
    return(x[start:end])

def grepTarFile(x):
    start= x.find('/DenervaudImages/') + len('/DenervaudImages/')
    end= x.find('.tar.bzip') + len('.tar.bzip')
    return(x[start:end])

zipFiles['FileName'] = pd.Series(zipFiles.iloc[:,0]).apply(grepName)
zipFiles['TarName'] = pd.Series(zipFiles.iloc[:,0]).apply(grepTarFile)


# %%
# Now run file I/O, and create directories and extract zip files

baseDirectory = '/Volumes/SeagatePortabl/DenervaudImages/'
os.chdir(baseDirectory)

for i in np.arange(0, len(zipFiles), 1):
    newDirectory = "/Volumes/SeagatePortabl/DenervaudImages/"+zipFiles.iloc[i,1]
    os.mkdir(newDirectory)

    subprocess.run("tar -xf " + zipFiles.iloc[i,2] + " -C " + newDirectory, shell=True)