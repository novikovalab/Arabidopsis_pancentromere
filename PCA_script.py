from Bio.Seq import Seq
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import pandas as pd

def simple_rev(s,n):
    '''
    Function to calculate kmer frequencies in repeat monomers
    n - kmer length
    s - sequence in Bio.Seq object
    '''
    dictionary = {}
    total = 0
    rev_s = Seq(s).reverse_complement()
    for i in range(len(s)-(n-1)): # (n-1) to get last element
        k = i+n
        if s[i:k] in dictionary:
            dictionary[s[i:k]] += 1
        else:
            dictionary.update({s[i:k]:1})
        total += 1 # doing it here to avoid sum(dictionary.values())
        if str(rev_s[i:k]) in dictionary:
            dictionary[str(rev_s[i:k])] += 1    
        else:
            dictionary.update({str(rev_s[i:k]):1})
        total += 1 # doing it here to avoid sum(dictionary.values())
    for key, value in dictionary.items():
        dictionary[key] = value/total
    return dictionary
    
# Reading the table, making the kmer matrix
montab = pd.read_csv("Pref_filt_repeats.bed", sep=" ", header=None)
montab["kmers"] = montab.iloc[:,4].apply(lambda x: simple_rev(x,5))
kmertab = montab["kmers"].apply(pd.Series)
kmertab = kmertab.fillna(0)

# Buiilding PCA
x = StandardScaler().fit_transform(kmertab)
pca = PCA(n_components=3)

principalComponents = pca.fit_transform(x)

principalDf2 = pd.DataFrame(data = principalComponents2, columns = ['principal component 1', 'principal component 2', 'principal component 3'])
pca.explained_variance_ratio_
fig = plt.figure(figsize = (8,8))

ax = fig.add_subplot(1,1,1) 
ax.scatter(principalDf2.loc[:, 'principal component 1']
               , principalDf2.loc[:, 'principal component 2'])
               
plt.savefig('The_pca.png')

# Consult the PCA to split monomers into repeat types, example
montab["typ"] = ""
montab[(principalDf["principal component 1"] < 11) & (principalDf["principal component 2"] < 9)].iloc[:,3].value_counts()
montab["typ"][(principalDf["principal component 1"] > 7) & (principalDf["principal component 2"] > 8)] = "AthCEN159"

# Project on TE8 PCA if not enough pAge2 repeats

kmertab = pd.concat([kmertab2, kmertab]).iloc[-(kmertab.shape[0]):,:]
kmertab = kmertab.fillna(0)
x = StandardScaler().fit_transform(kmertab)
test_pca = pca.transform(x)
fig = plt.figure(figsize = (8,8))
testDf = pd.DataFrame(data = test_pca
             , columns = ['principal component 1', 'principal component 2', 'principal component 3'])
ax = fig.add_subplot(1,1,1) 
ax.scatter(testDf.loc[:, 'principal component 1']
               , testDf.loc[:, 'principal component 2'])
               
# Making PCA with violinplots, example
sue_A = ["ASS3_AA", "AS150_AA", "AS530_AA", "N22666_AA"]
pallllete = {"lyrata": "#FDC067FF", "halleri": "#468892FF", "pedemontana": "#03A62CFF", "cebennensis": "#025940FF", "croatica": "#240E31FF", "arenosa": "#751C6DFF", "kamch_lyr": "#FD6F30FF", "kamch_hal": "#BADE86FF", "suecica_AA": "#CB6BCEFF", "suecica_AT": "#E995EBFF"}
fig = plt.figure(figsize = (16,16)) 
ax = fig.add_subplot(2,2,1)
sns.violinplot(x=principalDf.loc[:, 'principal component 1'], y=montab.loc[:, 'chrsm'], hue=montab.loc[:, 'species'], orient='h', palette=pallllete, rasterized=True)
plt.ylabel("chromosome")
plt.xlim([-15, 25])
ax = fig.add_subplot(2,2,3)
indicesToKeep =  (montab.loc[:,5].isin(sue_A))
ax.scatter(principalDf.loc[indicesToKeep, 'principal component 1'], principalDf.loc[indicesToKeep, 'principal component 2'], c = '#CB6BCEFF' , s = 50, rasterized=True)
indicesToKeep = (montab.loc[:,5]== "BaltAre")
ax.scatter(principalDf.loc[indicesToKeep, 'principal component 1'], principalDf.loc[indicesToKeep, 'principal component 2'], c = '#751C6DFF' , s = 50, rasterized=True)
indicesToKeep =  (montab.loc[:,5]=="Croatica")
ax.scatter(principalDf.loc[indicesToKeep, 'principal component 1'], principalDf.loc[indicesToKeep, 'principal component 2'], c = '#240E31FF' , s = 50, rasterized=True)
plt.ylabel("principal component 2")
plt.xlabel("principal component 1")
plt.xlim([-15, 25])
plt.ylim([-15, 27])
ax = fig.add_subplot(2,2,4)
sns.violinplot(y=principalDf.loc[:, 'principal component 2'], x=montab.loc[:, 'chrsm'], hue=montab.loc[:, 'species'], legend=False, palette=pallllete, rasterize=True)
plt.ylim([-15, 27])
plt.xlabel("chromosome")
