# %%

# IMPORT ALL PACKAGES REQUIRED
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


#%%
# Set the path, and function to load the data

def acqData(csvFile):
    data = pd.read_csv(csvFile, encoding = "latin")
    # data.dropna(subset = ["ORF (SGD)"], inplace=True)
    # data.columns = ['ORF', 'rad53d_log2', 'mec1d_log2', '322ArrayGene']
    # data.columns = inputColumns
    return(data)

data = acqData(csvFile="/Users/brandonho/Desktop/20211206_ChongHitsAnalyzed_ProperAnalysis_99Percentile_Part3.csv")


# Chong's treatment have labels in them, and we just want to know whether it
# was treated or not. The following determiens this

def getTreatment(x):
    return("".join(re.split("[^a-zA-Z]*", x)))

data['Treatment'] = data.iloc[0:]['Treatment'].apply(getTreatment)



#%%
# Now to quantify calls as increased or decreased on localization.
# This requires determining the thresholds for those calls using
# untreated cells for a given protein

def threshold(data):
    
    protein = []
    increasedThreshold = []
    decreasedThreshold = []
    # count = 0

    for i in np.unique(data['Protein']):
        # count = count + 1
        # print(str(count) + " out of " + str(len(np.unique(data['Protein']))))
        untreatedSubset = data[(data['Protein']==i) & (data['Treatment'] =="WT")]
        LOCvalues = untreatedSubset['95thPercentile']
        protein.append(i)
        increasedThreshold.append(np.median(LOCvalues) + 1.5*stats.median_abs_deviation(LOCvalues, scale='normal'))
        decreasedThreshold.append(np.median(LOCvalues) - 1.5*stats.median_abs_deviation(LOCvalues, scale='normal'))

    thresholdDF = {
        'Protein' : protein,
        'IncreasedThreshold' : increasedThreshold,
        'DecreasedThreshold' : decreasedThreshold
    }

    return(
        pd.DataFrame(thresholdDF)
    )

thresholdDF = threshold(data)




# %%
# Now compare with the thresholds, and designate cells as either increased
# or decreased in protein localization compared to the untreated

finalQuantification = pd.DataFrame([])
# count = 0

for i in np.unique(data['Protein']):

    protein = []
    increasedCells = []
    decreasedCells = []
    treatment = []
    totalCells = []

    proteinSubset = data[data['Protein']==i]

    for j in np.unique(proteinSubset['Treatment']):

        # Obtain the parameters - find the thresholds for given protein
        # to compare to
        treatmentSubset = proteinSubset[proteinSubset['Treatment'] == j]
        increasedThreshold = float(thresholdDF[thresholdDF['Protein'] == i]['IncreasedThreshold'])
        decreasedThreshold = float(thresholdDF[thresholdDF['Protein'] == i]['DecreasedThreshold'])

        # Determine the number of cells that exceed these defined
        # thresholdsf or increased and decreased localization
        increasedCells.append(len(treatmentSubset[treatmentSubset['95thPercentile'] > increasedThreshold]))
        decreasedCells.append(len(treatmentSubset[treatmentSubset['95thPercentile'] < decreasedThreshold]))
        treatment.append(j)
        protein.append(i)
        totalCells.append(len(treatmentSubset))

    # Assemble the above data into a dataframe, and merge with the final
    # quantification table so that we can group all proteins and treatments
    # together
    proteinCalls = {
        'Protein' : protein,
        'Treatment' : treatment,
        'IncreasedNumbers' : increasedCells,
        'DecreasedNumbers' : decreasedCells,
        'TotalNumbers' : totalCells
    }

    proteinCalls = pd.DataFrame(proteinCalls)
    finalQuantification = pd.concat([finalQuantification, proteinCalls])

    # Progress of the loop
    # count = count + 1
    # print(str(count) + " out of " + str(len(np.unique(data['Protein']))))



# %%
# Save file to target destination
desktopPath = "/Users/brandonho/Desktop/"
os.chdir(desktopPath)

finalQuantification.to_csv("localizationNumbers_Chong_99Percentile_1p5MADThreshold.csv")


# %%
# Finally, we can determine which ones reach our threshold
# for changes in localization. And that will be determined by 
# looking for protein localization changes greater than 1.5 fold change

protein = []
foldChange = []

finalQuantification['PercentUp'] = finalQuantification['IncreasedNumbers']/finalQuantification['TotalNumbers']

for i in np.unique(finalQuantification['Protein']):
    proteinSubset = finalQuantification[finalQuantification['Protein'] == i]
    HU = float(proteinSubset[proteinSubset['Treatment']=='HU']['PercentUp'])
    WT = float(proteinSubset[proteinSubset['Treatment']=='WT']['PercentUp'])
    foldChange.append(np.log2(HU/WT))
    protein.append(i)

finalCalls = {
    'Protein' : protein,
    'log2FoldChange' : foldChange
}

finalCalls = pd.DataFrame(finalCalls)
# finalCalls.to_csv("log2FoldChange_Chong_99Percentile_1p5MADThreshold.csv")



# %%
# For a better comparison to Chong, maybe it will be worth determining
# the actual changes in the number/percent of cells (difference) and
# not the fold change

protein = []
difference = []

finalQuantification['PercentUp'] = finalQuantification['IncreasedNumbers']/finalQuantification['TotalNumbers']

for i in np.unique(finalQuantification['Protein']):
    proteinSubset = finalQuantification[finalQuantification['Protein'] == i]
    HU = float(proteinSubset[proteinSubset['Treatment']=='HU']['PercentUp'])
    WT = float(proteinSubset[proteinSubset['Treatment']=='WT']['PercentUp'])
    difference.append(HU - WT)
    protein.append(i)

finalCalls = {
    'Protein' : protein,
    'Difference' : difference
}

finalCalls = pd.DataFrame(finalCalls)
finalCalls.to_csv("percentDifference_Chong_99Percentile_1p5MADThreshold.csv")













# %%
# Calculating the coefficient of variation and comparing with
# protein abundance in molecules per cell

data = pd.read_csv("/Users/brandonho/Dropbox (Grant Brown's Lab)/Complete Data (Brandon Ho)/Year 7/wildtypeLocalizationScreenAnalysis/20211212_WTHU_RawCalls.csv")

CVquantification = pd.DataFrame([])

count = 0

for i in np.unique(data['StrainID']):

    count = count + 1
    print(count)

    coefVar = []
    time = []
    strainID = []
    totalCells = []

    proteinSubset = data[data['StrainID']==i]

    for j in np.unique(proteinSubset['TimePoint']):

        # Calculate CV at each time point
        timeSubset = proteinSubset[proteinSubset['TimePoint'] == j]

        # Coefficien variation function
        cv = lambda x: np.std(x, ddof=1) / np.mean(x) * 100

        coefVar.append(cv(timeSubset['95thPercentile']))

        time.append(j)
        strainID.append(i)
        totalCells.append(len(timeSubset))

    # Assemble the above data into a dataframe, and merge with the final
    # quantification table so that we can group all proteins and treatments
    # together
    proteinCalls = {
        'Strain' : strainID,
        'TimePoint' : time,
        'CoefficientVariation' : coefVar,
        'TotalNumbers' : totalCells
    }

    proteinCalls = pd.DataFrame(proteinCalls)
    CVquantification = pd.concat([CVquantification, proteinCalls])



# %%
fileGuide = pd.read_csv("/Users/brandonho/Dropbox (Grant Brown's Lab)/Complete Data (Brandon Ho)/Year 6/Paper_Figures_Tables/CSVFilesForAnalysis/LocalizationQuantification/FileGuide.csv")
ORFs = pd.read_csv("/Users/brandonho/Dropbox (Grant Brown's Lab)/Complete Data (Brandon Ho)/Year 6/Paper_Figures_Tables/CSVFilesForAnalysis/LocalizationQuantification/Mec1 Data Complete.csv")

fileGuideSubset = fileGuide[['Strain', 'ORF']]
ORFsubset = ORFs[['ORF', 'Gene']]

strainIDToORF = pd.merge(fileGuideSubset, ORFsubset)

CVquantificationWGene = pd.merge(CVquantification, strainIDToORF)


CVfinal = CVquantificationWGene.pivot(index='Gene', columns='TimePoint', values='CoefficientVariation')



desktopPath = "/Users/brandonho/Desktop"
os.chdir(desktopPath)

CVfinal.to_csv("HU_CV.csv")










# %%
commonHits=pd.read_csv("/Users/brandonho/Dropbox (Grant Brown's Lab)/Complete Data (Brandon Ho)/Year 7/wildtypeLocalizationScreenAnalysis/HUandMMS_CommonHits.csv")
data=pd.read_csv("/Users/brandonho/Dropbox (Grant Brown's Lab)/Complete Data (Brandon Ho)/Year 7/wildtypeLocalizationScreenAnalysis/20210809_HUMMS_penetrance.csv")

targetColumns=['Gene',
               'MMS0', 'MMS30', 'MMS60', 'MMS90', 'MMS120', 'MMS180', 'MMS240',
               'HU0', 'HU30', 'HU60', 'HU90', 'HU120', 'HU180', 'HU240']

data=data[targetColumns]

data = data[data['Gene'].isin(commonHits['CommonHits'])]
dataBackup = data.copy()



# %%
# First, figure out what the average rate is for each protein

def avgRate(x):
    max = np.max(x)
    min = np.min(x)
    return((max-min)/240)

data['MMSavgRate']=data.iloc[:,1:8].apply(avgRate, axis=1)
data['HUavgRate']=data.iloc[:,8:15].apply(avgRate, axis=1)

data['MMSmax']=data.iloc[:,1:8].apply(np.max, axis=1)
data['HUmax']=data.iloc[:,8:15].apply(np.max, axis=1)




# %%
data['pMaxMMS_0']=data['MMS0']/data['MMSmax']
data['pMaxMMS_30']=data['MMS30']/data['MMSmax']
data['pMaxMMS_60']=data['MMS60']/data['MMSmax']
data['pMaxMMS_90']=data['MMS90']/data['MMSmax']
data['pMaxMMS_120']=data['MMS120']/data['MMSmax']
data['pMaxMMS_180']=data['MMS180']/data['MMSmax']
data['pMaxMMS_240']=data['MMS240']/data['MMSmax']

data['pMaxHU_0']=data['HU0']/data['HUmax']
data['pMaxHU_30']=data['HU30']/data['HUmax']
data['pMaxHU_60']=data['HU60']/data['HUmax']
data['pMaxHU_90']=data['HU90']/data['HUmax']
data['pMaxHU_120']=data['HU120']/data['HUmax']
data['pMaxHU_180']=data['HU180']/data['HUmax']
data['pMaxHU_240']=data['HU240']/data['HUmax']



# %%

MMS = data.iloc[:,19:25]
MMSmat = MMS.to_numpy()


HU = data.iloc[:,27:33]
HUmat = HU.to_numpy()

diff = pd.DataFrame(np.log2(MMSmat/HUmat))
diff.index=data['Gene']
diff[diff < -100] = 0
diff[diff > 100] = 0


# %%
# Now Kmeans this thing

kmeans = KMeans(
    init="random",
    n_clusters=4,
    n_init=10,
    max_iter=300,
    random_state=42
)

kmeans.fit(diff)
diff['Cluster'] = kmeans.labels_
diff = diff.sort_values(by='Cluster')

# sns.color_palette("icefire", as_cmap=True)
sns.heatmap(diff.iloc[:,0:6],
            # cmap=sns.diverging_palette(250, 30, l=65, center="dark", as_cmap=True),
            cmap='vlag',
            center=0,
            vmin=-3,
            vmax=3)

# os.chdir("/Users/brandonho/Desktop/")
# plt.savefig('heatmap.pdf')


# %%

sns.regplot(x='MMSavgRate', y='HUavgRate', data=data)
plt.savefig('correlationplot.pdf')



# %%

correlation = []

for i in np.unique(data['Gene']):

    dataSubset = data[data['Gene']==i]

    MMSvals = dataSubset.iloc[0,1:8]
    HUvals = dataSubset.iloc[0,8:15]

    correlation.append(stats.pearsonr(MMSvals, HUvals)[0])

correlation = pd.DataFrame(correlation)
correlation.columns = ['CorVal']



# %%
