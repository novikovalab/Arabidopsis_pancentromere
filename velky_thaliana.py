import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
fig = plt.figure(figsize=(22,20))

n_samp=55
n = 0
te8_rep = pd.read_csv("mabs_assemblies/final_trash/KZ7742_typen_reverse.bed", sep="\t", header=None)
cen_pos = pd.read_csv('mabs_assemblies/v2_assemblies/expand_centromeres/KZ7742_cen_reverse.bed', sep="\t", header=None)
ns = pd.read_csv('mabs_assemblies/mapping/KZ7742_rev_Ns.bed', sep="\t", header=None)
for i in range(1,9):
    ax = fig.add_subplot(n_samp,8,n*8+i)
    sc = "scaffold_" + str(i)
    chr8_rep = te8_rep[(te8_rep.iloc[:,0] == sc)]
    plt.eventplot(chr8_rep[chr8_rep.iloc[:,3] == "Co_168"].iloc[:,1], linewidths=0.01, orientation='horizontal', colors="#BCBD22FF")
    plt.eventplot(chr8_rep[chr8_rep.iloc[:,3] == "Co_48"].iloc[:,1], linewidths=0.01, orientation='horizontal', colors="#1F83B4FF")
    ax.axes.xaxis.set_visible(False)
    ax.axes.yaxis.set_visible(False)
    ns_sc = ns[ns.iloc[:,0] == sc]
    plt.eventplot(ns_sc.iloc[:,1], orientation='horizontal', colors="black", lineoffsets=0.2, linewidths=0.5)
    plt.eventplot(ns_sc.iloc[:,1], orientation='horizontal', colors="black", lineoffsets=1.8, linewidths=0.5)
    middle = (cen_pos.iloc[i-1,2] + cen_pos.iloc[i-1,1])/2
    plt.xlim([middle - 5000000, middle + 5000000])
    plt.ylim([0.1, 1.9])

n=1


te8_rep = pd.read_csv("mabs_assemblies/final_trash/Col-0_typen.bed", sep="\t")
cen_pos = pd.read_csv('mabs_assemblies/v2_assemblies/Col-0_merged_centromeres.bed', sep="\t", header=None)
ns = pd.read_csv('mabs_assemblies/mapping/Col-0_Ns.bed', sep="\t", header=None)
for i,nc in zip([1,2,3,4,5], [1,3,5,6,7]):
    ax = fig.add_subplot(n_samp,8,n*8+nc)
    sc = "Chr" + str(i)
    chr8_rep = te8_rep[ (te8_rep.iloc[:,0].str.contains(sc))]
    plt.eventplot(chr8_rep[chr8_rep.iloc[:,3] == "AthCEN178"].iloc[:,1], linewidths=0.01, orientation='horizontal', colors="#12A2A8FF")
    plt.eventplot(chr8_rep[chr8_rep.iloc[:,3] == "AthCEN159"].iloc[:,1], linewidths=0.01, orientation='horizontal', colors="#C7519CFF")
    ax.axes.xaxis.set_visible(False)
    ax.axes.yaxis.set_visible(False)
    ns_sc = ns[ns.iloc[:,0].str.contains(sc)]
    plt.eventplot(ns_sc.iloc[:,1], orientation='horizontal', colors="black", lineoffsets=0.2, linewidths=0.5)
    plt.eventplot(ns_sc.iloc[:,1], orientation='horizontal', colors="black", lineoffsets=1.8, linewidths=0.5)
    middle = (cen_pos.iloc[i-1,2] + cen_pos.iloc[i-1,1])/2
    plt.xlim([middle - 5000000, middle + 5000000])
    plt.ylim([0.1, 1.9])
n=2


te8_rep = pd.read_csv("mabs_assemblies/final_trash/N22666_v1_AT_pca_typen.bed", sep="\t")
cen_pos = pd.read_csv('mabs_assemblies/v2_assemblies/expand_centromeres/N22666_v1_AT_centromere_exp.bed', sep="\t", header=None)
ns = pd.read_csv('mabs_assemblies/mapping/N22666_v1_AT_Ns.bed', sep="\t", header=None)
for i,nc in zip([1,2,3,4,5], [1,3,5,6,7]):
    ax = fig.add_subplot(n_samp,8,n*8+nc)
    sc = "Chr" + str(i)
    chr8_rep = te8_rep[ (te8_rep.iloc[:,0].str.contains(sc))]
    plt.eventplot(chr8_rep[chr8_rep.iloc[:,3] == "AthCEN178"].iloc[:,1], linewidths=0.01, orientation='horizontal', colors="#12A2A8FF")
    plt.eventplot(chr8_rep[chr8_rep.iloc[:,3] == "AthCEN159"].iloc[:,1], linewidths=0.01, orientation='horizontal', colors="#C7519CFF")
    ax.axes.xaxis.set_visible(False)
    ax.axes.yaxis.set_visible(False)
    ns_sc = ns[ns.iloc[:,0].str.contains(sc)]
    plt.eventplot(ns_sc.iloc[:,1], orientation='horizontal', colors="black", lineoffsets=0.2, linewidths=0.5)
    plt.eventplot(ns_sc.iloc[:,1], orientation='horizontal', colors="black", lineoffsets=1.8, linewidths=0.5)
    middle = (cen_pos.iloc[i-1,2] + cen_pos.iloc[i-1,1])/2
    plt.xlim([middle - 5000000, middle + 5000000])
    plt.ylim([0.1, 1.9])

n =5

te8_rep = pd.read_csv("mabs_assemblies/final_trash/AS150_v6_AT_fix_typen.bed", sep="\t")
cen_pos = pd.read_csv('mabs_assemblies/v2_assemblies/AS150_v6_AT_centromeres.bed', sep="\t", header=None)
ns = pd.read_csv('mabs_assemblies/mapping/AS150_v6_AT_Ns.bed', sep="\t", header=None)
for i,nc in zip([1,2,3,4,5], [1,3,5,6,7]):
    ax = fig.add_subplot(n_samp,8,n*8+nc)
    sc = "Chr" + str(i)
    chr8_rep = te8_rep[ (te8_rep.iloc[:,0].str.contains(sc))]
    plt.eventplot(chr8_rep[chr8_rep.iloc[:,3] == "AthCEN178"].iloc[:,1], linewidths=0.01, orientation='horizontal', colors="#12A2A8FF")
    plt.eventplot(chr8_rep[chr8_rep.iloc[:,3] == "AthCEN159"].iloc[:,1], linewidths=0.01, orientation='horizontal', colors="#C7519CFF")
    plt.eventplot(chr8_rep[chr8_rep.iloc[:,3] == "pAa"].iloc[:,1], orientation='horizontal', colors="#6F63BBFF")
    ax.axes.xaxis.set_visible(False)
    ax.axes.yaxis.set_visible(False)
    ns_sc = ns[ns.iloc[:,0].str.contains(sc)]
    plt.eventplot(ns_sc.iloc[:,1], orientation='horizontal', colors="black", lineoffsets=0.2, linewidths=0.5)
    plt.eventplot(ns_sc.iloc[:,1], orientation='horizontal', colors="black", lineoffsets=1.8, linewidths=0.5)
    middle = (cen_pos.iloc[i-1,2] + cen_pos.iloc[i-1,1])/2
    plt.xlim([middle - 5000000, middle + 5000000])
    plt.ylim([0.1, 1.9])

n=8


te8_rep = pd.read_csv("mabs_assemblies/final_trash/N22666_AA_typen_reclass.bed", sep="\t")
cen_pos = pd.read_csv('mabs_assemblies/v2_assemblies/expand_centromeres/N22666_v1_AA_centromere_exp.bed', sep="\t", header=None)
ns = pd.read_csv('mabs_assemblies/mapping/N22666_v1_AA_Ns.bed', sep="\t", header=None)
for i in range(1,9):
    ax = fig.add_subplot(n_samp,8,n*8+i)
    sc = "Chr" + str(i+5)
    chr8_rep = te8_rep[ (te8_rep.iloc[:,0].str.contains(sc))]
    plt.eventplot(chr8_rep[chr8_rep.iloc[:,3] == "pAge2"].iloc[:,1], linewidths=0.01, orientation='horizontal', colors="#FFAA0EFF")
    plt.eventplot(chr8_rep[chr8_rep.iloc[:,3] == "pAa"].iloc[:,1], linewidths=0.01, orientation='horizontal', colors="#6F63BBFF")
    plt.eventplot(chr8_rep[chr8_rep.iloc[:,3] == "pAge1"].iloc[:,1], linewidths=0.01, orientation='horizontal', colors="#2CA030FF")
    ns_sc = ns[ns.iloc[:,0] == sc]
    plt.eventplot(chr8_rep[chr8_rep.iloc[:,3] == "AthCEN159"].iloc[:,1], linewidths=0.01, orientation='horizontal', colors="#C7519CFF")
    plt.eventplot(chr8_rep[chr8_rep.iloc[:,3] == "AthCEN178"].iloc[:,1], linewidths=0.01, orientation='horizontal', colors="#12A2A8FF")
    ax.axes.xaxis.set_visible(False)
    ax.axes.yaxis.set_visible(False)
    plt.eventplot(ns_sc.iloc[:,1], orientation='horizontal', colors="black", lineoffsets=0.2, linewidths=0.5)
    plt.eventplot(ns_sc.iloc[:,1], orientation='horizontal', colors="black", lineoffsets=1.8, linewidths=0.5)
    middle = (cen_pos.iloc[i-1,2] + cen_pos.iloc[i-1,1])/2
    plt.xlim([middle - 5000000, middle + 5000000])
    plt.ylim([0.1, 1.9])
n=3


te8_rep = pd.read_csv("mabs_assemblies/final_trash/ASS3_v6_AT_typen.bed", sep="\t")
cen_pos = pd.read_csv('mabs_assemblies/v2_assemblies/ASS3_v6_AT_merged_centromeres.bed', sep="\t", header=None)
ns = pd.read_csv('mabs_assemblies/mapping/ASS3_v6_AT_Ns.bed', sep="\t", header=None)
for i,nc in zip([1,2,3,4,5], [1,3,5,6,7]):
    ax = fig.add_subplot(n_samp,8,n*8+nc)
    sc = "Chr" + str(i)
    chr8_rep = te8_rep[ (te8_rep.iloc[:,0].str.contains(sc))]
    plt.eventplot(chr8_rep[chr8_rep.iloc[:,3] == "AthCEN178"].iloc[:,1], linewidths=0.01, orientation='horizontal', colors="#12A2A8FF")
    plt.eventplot(chr8_rep[chr8_rep.iloc[:,3] == "AthCEN159"].iloc[:,1], linewidths=0.01, orientation='horizontal', colors="#C7519CFF")
    plt.eventplot(chr8_rep[chr8_rep.iloc[:,3] == "pAa"].iloc[:,1], orientation='horizontal', colors="#6F63BBFF", linewidths=0.01)
    ax.axes.xaxis.set_visible(False)
    ax.axes.yaxis.set_visible(False)
    ns_sc = ns[ns.iloc[:,0].str.contains(sc)]
    plt.eventplot(ns_sc.iloc[:,1], orientation='horizontal', colors="black", lineoffsets=0.2, linewidths=0.5)
    plt.eventplot(ns_sc.iloc[:,1], orientation='horizontal', colors="black", lineoffsets=1.8, linewidths=0.5)
    middle = (cen_pos.iloc[i-1,2] + cen_pos.iloc[i-1,1])/2
    plt.xlim([middle - 5000000, middle + 5000000])
    plt.ylim([0.1, 1.9])

    
n=4


te8_rep = pd.read_csv("mabs_assemblies/final_trash/AS530_v6_AT_typen_fix.bed", sep="\t")
cen_pos = pd.read_csv('mabs_assemblies/v2_assemblies/AS530_v6_AT_merged_centromeres.bed', sep="\t", header=None)
ns = pd.read_csv('mabs_assemblies/mapping/AS530_v6_AT_Ns.bed', sep="\t")
for i,nc in zip([1,2,3,4,5], [1,3,5,6,7]):
    ax = fig.add_subplot(n_samp,8,n*8+nc)
    sc = "Chr" + str(i)
    chr8_rep = te8_rep[ (te8_rep.iloc[:,0].str.contains(sc))]
    plt.eventplot(chr8_rep[chr8_rep.iloc[:,3] == "AthCEN178"].iloc[:,1], linewidths=0.01, orientation='horizontal', colors="#12A2A8FF")
    plt.eventplot(chr8_rep[chr8_rep.iloc[:,3] == "AthCEN159"].iloc[:,1], linewidths=0.01, orientation='horizontal', colors="#C7519CFF")
    plt.eventplot(chr8_rep[chr8_rep.iloc[:,3] == "pAa"].iloc[:,1], orientation='horizontal', colors="#6F63BBFF")
    ax.axes.xaxis.set_visible(False)
    ax.axes.yaxis.set_visible(False)
    ns_sc = ns[ns.iloc[:,0].str.contains(sc)]
    plt.eventplot(ns_sc.iloc[:,1], orientation='horizontal', colors="black", lineoffsets=0.2, linewidths=0.5)
    plt.eventplot(ns_sc.iloc[:,1], orientation='horizontal', colors="black", lineoffsets=1.8, linewidths=0.5)
    middle = (cen_pos.iloc[i-1,2] + cen_pos.iloc[i-1,1])/2
    plt.xlim([middle - 5000000, middle + 5000000])
    plt.ylim([0.1, 1.9])
n=28
te8_rep = pd.read_csv("mabs_assemblies/final_trash/dach2_v1_satellites_fixed_pAge2.bed", sep="\t", header=None)
cen_pos = pd.read_csv('mabs_assemblies/v2_assemblies/dach2_v1_merge_centromeres.bed', sep="\t", header=None)
ns = pd.read_csv('mabs_assemblies/mapping/dach2_v1_Ns.bed', sep="\t", header=None)
for i in range(1,9):
    ax = fig.add_subplot(n_samp,8,n*8+i)
    sc = "chr_" + str(i)
    chr8_rep = te8_rep[(te8_rep.iloc[:,0] == sc)]
    plt.eventplot(chr8_rep[chr8_rep.iloc[:,3] == "pAge1"].iloc[:,1], linewidths=0.01, orientation='horizontal', colors="#2CA030FF")
    plt.eventplot(chr8_rep[chr8_rep.iloc[:,3] == "pAa"].iloc[:,1], linewidths=0.01, orientation='horizontal', colors="#6F63BBFF")
    plt.eventplot(chr8_rep[chr8_rep.iloc[:,3] == "pAge2"].iloc[:,1], linewidths=0.01, orientation='horizontal', colors="#FFAA0EFF")
    plt.eventplot(chr8_rep[chr8_rep.iloc[:,3] == "AthCEN159"].iloc[:,1], linewidths=0.01, orientation='horizontal', colors="#C7519CFF")
    ax.axes.xaxis.set_visible(False)
    ax.axes.yaxis.set_visible(False)
    ns_sc = ns[ns.iloc[:,0] == sc]
    plt.eventplot(ns_sc.iloc[:,1], orientation='horizontal', colors="black", lineoffsets=0.2, linewidths=0.5)
    plt.eventplot(ns_sc.iloc[:,1], orientation='horizontal', colors="black", lineoffsets=1.8, linewidths=0.5)
    middle = (cen_pos.iloc[i-1,2] + cen_pos.iloc[i-1,1])/2
    plt.xlim([middle - 5000000, middle + 5000000])
    plt.ylim([0.1, 1.9])
n=29
te8_rep = pd.read_csv("mabs_assemblies/final_trash/AK1_v3_pca_typen.bed", sep="\t", header=None)
cen_pos = pd.read_csv('mabs_assemblies/v2_assemblies/expand_centromeres/AK1_v3_centr.bed', sep="\t", header=None)
ns = pd.read_csv('mabs_assemblies/mapping/AK1_v3_Ns.bed', sep="\t", header=None)
for i in range(1,9):
    ax = fig.add_subplot(n_samp,8,n*8+i)
    sc = "chr_" + str(i)
    chr8_rep = te8_rep[(te8_rep.iloc[:,0] == sc)]
    plt.eventplot(chr8_rep[chr8_rep.iloc[:,4] == "pAge1"].iloc[:,1], linewidths=0.01, orientation='horizontal', colors="#2CA030FF")
    plt.eventplot(chr8_rep[chr8_rep.iloc[:,4] == "pAa"].iloc[:,1], linewidths=0.01, orientation='horizontal', colors="#6F63BBFF")
    plt.eventplot(chr8_rep[chr8_rep.iloc[:,4] == "pAge2"].iloc[:,1], linewidths=0.01, orientation='horizontal', colors="#FFAA0EFF")
    plt.eventplot(chr8_rep[chr8_rep.iloc[:,4] == "AthCEN159"].iloc[:,1], linewidths=0.01, orientation='horizontal', colors="#C7519CFF")
    ax.axes.xaxis.set_visible(False)
    ax.axes.yaxis.set_visible(False)
    ns_sc = ns[ns.iloc[:,0] == sc]
    plt.eventplot(ns_sc.iloc[:,1], orientation='horizontal', colors="black", lineoffsets=0.2, linewidths=0.5)
    plt.eventplot(ns_sc.iloc[:,1], orientation='horizontal', colors="black", lineoffsets=1.8, linewidths=0.5)
    middle = (cen_pos.iloc[i-1,2] + cen_pos.iloc[i-1,1])/2
    plt.xlim([middle - 5000000, middle + 5000000])
    plt.ylim([0.1, 1.9])
n=31
te8_rep = pd.read_csv("mabs_assemblies/final_trash/PU6_pca_typen.bed", sep="\t", header=None)
cen_pos = pd.read_csv('mabs_assemblies/v2_assemblies/expand_centromeres/PU6_v3.bed', sep="\t", header=None)
for i in range(1,9):
    ax = fig.add_subplot(n_samp,8,n*8+i)
    sc = "scaffold_" + str(i)
    chr8_rep = te8_rep[ (te8_rep.iloc[:,0].str.contains(sc))]
    plt.eventplot(chr8_rep[chr8_rep.iloc[:,3] == "pAge1"].iloc[:,1], linewidths=0.01, orientation='horizontal', colors="#2CA030FF")
    plt.eventplot(chr8_rep[chr8_rep.iloc[:,3] == "pAge2"].iloc[:,1], linewidths=0.01, orientation='horizontal', colors="#FFAA0EFF")
    plt.eventplot(chr8_rep[chr8_rep.iloc[:,3] == "pAa"].iloc[:,1], linewidths=0.01, orientation='horizontal', colors="#6F63BBFF")
    ax.axes.xaxis.set_visible(False)
    ax.axes.yaxis.set_visible(False)
    middle = (cen_pos.iloc[i-1,2] + cen_pos.iloc[i-1,1])/2
    plt.xlim([middle - 5000000, middle + 5000000])
    plt.ylim([0.1, 1.9])

n=27
te8_rep = pd.read_csv("mabs_assemblies/final_trash/kron_v2_satellites_fixed_pAge2.bed", sep="\t", header=None)
cen_pos = pd.read_csv('mabs_assemblies/v2_assemblies/kron_v1_merge_centromeres.bed', sep="\t", header=None)
ns = pd.read_csv('mabs_assemblies/mapping/kron_v1_Ns.bed', sep="\t", header=None)
for i in range(1,9):
    ax = fig.add_subplot(n_samp,8,n*8+i)
    sc = "chr_" + str(i)
    chr8_rep = te8_rep[(te8_rep.iloc[:,0] == sc)]
    plt.eventplot(chr8_rep[chr8_rep.iloc[:,3] == "pAge1"].iloc[:,1], linewidths=0.01, orientation='horizontal', colors="#2CA030FF")
    plt.eventplot(chr8_rep[chr8_rep.iloc[:,3] == "pAa"].iloc[:,1], linewidths=0.01, orientation='horizontal', colors="#6F63BBFF")
    plt.eventplot(chr8_rep[chr8_rep.iloc[:,3] == "pAge2"].iloc[:,1], linewidths=0.01, orientation='horizontal', colors="#FFAA0EFF")
    ax.axes.xaxis.set_visible(False)
    ax.axes.yaxis.set_visible(False)
    ns_sc = ns[ns.iloc[:,0] == sc]
    plt.eventplot(ns_sc.iloc[:,1], orientation='horizontal', colors="black", lineoffsets=0.2, linewidths=0.5)
    plt.eventplot(ns_sc.iloc[:,1], orientation='horizontal', colors="black", lineoffsets=1.8, linewidths=0.5)
    middle = (cen_pos.iloc[i-1,2] + cen_pos.iloc[i-1,1])/2
    plt.xlim([middle - 5000000, middle + 5000000])
    plt.ylim([0.1, 1.9])
te8_rep = pd.read_csv("mabs_assemblies/final_trash/AK1_v3_pca_typen.bed", sep="\t", header=None)
cen_pos = pd.read_csv('mabs_assemblies/v2_assemblies/expand_centromeres/AK1_v3_centr.bed', sep="\t", header=None)
ns = pd.read_csv('mabs_assemblies/mapping/AK1_v3_Ns.bed', sep="\t", header=None)
n=18
for i in range(1,9):
    ax = fig.add_subplot(n_samp,8,n*8+i)
    sc = "chr_" + str(i+8)
    chr8_rep = te8_rep[(te8_rep.iloc[:,0] == sc)]
    plt.eventplot(chr8_rep[chr8_rep.iloc[:,4] == "pAge1"].iloc[:,1], linewidths=0.01, orientation='horizontal', colors="#2CA030FF")
    plt.eventplot(chr8_rep[chr8_rep.iloc[:,4] == "pAa"].iloc[:,1], linewidths=0.01, orientation='horizontal', colors="#6F63BBFF")
    plt.eventplot(chr8_rep[chr8_rep.iloc[:,4] == "pAge2"].iloc[:,1], linewidths=0.01, orientation='horizontal', colors="#FFAA0EFF")
    ax.axes.xaxis.set_visible(False)
    ax.axes.yaxis.set_visible(False)
    ns_sc = ns[ns.iloc[:,0] == sc]
    plt.eventplot(ns_sc.iloc[:,1], orientation='horizontal', colors="black", lineoffsets=0.2, linewidths=0.5)
    plt.eventplot(ns_sc.iloc[:,1], orientation='horizontal', colors="black", lineoffsets=1.8, linewidths=0.5)
    middle = (cen_pos.iloc[i-9,2] + cen_pos.iloc[i-9,1])/2
    plt.xlim([middle - 5000000, middle + 5000000])
    plt.ylim([0.1, 1.9])

n=21
te8_rep = pd.read_csv("mabs_assemblies/final_trash/dach2_v1_satellites_fixed_pAge2.bed", sep="\t", header=None)
cen_pos = pd.read_csv('mabs_assemblies/v2_assemblies/dach2_v1_merge_centromeres.bed', sep="\t", header=None)
ns = pd.read_csv('mabs_assemblies/mapping/dach2_v1_Ns.bed', sep="\t", header=None)
for i in range(1,9):
    ax = fig.add_subplot(n_samp,8,n*8+i)
    sc = "chr_" + str(i+8)
    chr8_rep = te8_rep[(te8_rep.iloc[:,0] == sc)]
    plt.eventplot(chr8_rep[chr8_rep.iloc[:,3] == "pAge1"].iloc[:,1], linewidths=0.01, orientation='horizontal', colors="#2CA030FF")
    plt.eventplot(chr8_rep[chr8_rep.iloc[:,3] == "pAa"].iloc[:,1], linewidths=0.01, orientation='horizontal', colors="#6F63BBFF")
    plt.eventplot(chr8_rep[chr8_rep.iloc[:,3] == "pAge2"].iloc[:,1], linewidths=0.01, orientation='horizontal', colors="#FFAA0EFF")
    ax.axes.xaxis.set_visible(False)
    ax.axes.yaxis.set_visible(False)
    ns_sc = ns[ns.iloc[:,0] == sc]
    plt.eventplot(ns_sc.iloc[:,1], orientation='horizontal', colors="black", lineoffsets=0.2, linewidths=0.5)
    plt.eventplot(ns_sc.iloc[:,1], orientation='horizontal', colors="black", lineoffsets=1.8, linewidths=0.5)
    middle = (cen_pos.iloc[i-9,2] + cen_pos.iloc[i-9,1])/2
    plt.xlim([middle - 5000000, middle + 5000000])
    plt.ylim([0.1, 1.9])
n=32
te8_rep = pd.read_csv("mabs_assemblies/final_trash/WS1_v3_from_TE8_typen.bed", sep="\t", header=None)
cen_pos = pd.read_csv('mabs_assemblies/v2_assemblies/expand_centromeres/WS1_v3.bed', sep="\t", header=None)
ns = pd.read_csv('mabs_assemblies/mapping/WS1_v3_Ns.bed', sep="\t", header=None)
for i in range(1,9):
    ax = fig.add_subplot(n_samp,8,n*8+i)
    sc = "scaffold_" + str(i)
    chr8_rep = te8_rep[ (te8_rep.iloc[:,0].str.contains(sc))]
    plt.eventplot(chr8_rep[chr8_rep.iloc[:,4] == "pAge1"].iloc[:,1], linewidths=0.01, orientation='horizontal', colors="#2CA030FF")
    plt.eventplot(chr8_rep[chr8_rep.iloc[:,4] == "pAge2"].iloc[:,1], linewidths=0.01, orientation='horizontal', colors="#FFAA0EFF")
    plt.eventplot(chr8_rep[chr8_rep.iloc[:,4] == "pAa"].iloc[:,1], linewidths=0.01, orientation='horizontal', colors="#6F63BBFF")
    ax.axes.xaxis.set_visible(False)
    ax.axes.yaxis.set_visible(False)
    ns_sc = ns[ns.iloc[:,0].str.contains(sc)]
    plt.eventplot(ns_sc.iloc[:,1], orientation='horizontal', colors="black", lineoffsets=0.2, linewidths=0.5)
    plt.eventplot(ns_sc.iloc[:,1], orientation='horizontal', colors="black", lineoffsets=1.8, linewidths=0.5)
    middle = (cen_pos.iloc[i-1,2] + cen_pos.iloc[i-1,1])/2
    plt.xlim([middle - 5000000, middle + 5000000])
    plt.ylim([0.1, 1.9])
    

n=20
te8_rep = pd.read_csv("mabs_assemblies/final_trash/kron_v2_satellites_fixed_pAge2.bed", sep="\t", header=None)
cen_pos = pd.read_csv('mabs_assemblies/v2_assemblies/kron_v1_merge_centromeres.bed', sep="\t", header=None)
ns = pd.read_csv('mabs_assemblies/mapping/kron_v1_Ns.bed', sep="\t", header=None)
for i in range(1,9):
    ax = fig.add_subplot(n_samp,8,n*8+i)
    sc = "chr_" + str(i+8)
    chr8_rep = te8_rep[(te8_rep.iloc[:,0] == sc)]
    plt.eventplot(chr8_rep[chr8_rep.iloc[:,3] == "pAge1"].iloc[:,1], linewidths=0.01, orientation='horizontal', colors="#2CA030FF")
    plt.eventplot(chr8_rep[chr8_rep.iloc[:,3] == "pAa"].iloc[:,1], linewidths=0.01, orientation='horizontal', colors="#6F63BBFF")
    plt.eventplot(chr8_rep[chr8_rep.iloc[:,3] == "pAge2"].iloc[:,1], linewidths=0.01, orientation='horizontal', colors="#FFAA0EFF")
    ax.axes.xaxis.set_visible(False)
    ax.axes.yaxis.set_visible(False)
    ns_sc = ns[ns.iloc[:,0] == sc]
    plt.eventplot(ns_sc.iloc[:,1], orientation='horizontal', colors="black", lineoffsets=0.2, linewidths=0.5)
    plt.eventplot(ns_sc.iloc[:,1], orientation='horizontal', colors="black", lineoffsets=1.8, linewidths=0.5)
    middle = (cen_pos.iloc[i-9,2] + cen_pos.iloc[i-9,1])/2
    plt.xlim([middle - 5000000, middle + 5000000])
    plt.ylim([0.1, 1.9])
n=19
te8_rep = pd.read_csv("mabs_assemblies/final_trash/ah_hkg_pca_typen.bed", sep="\t", header=None)
cen_pos = pd.read_csv('mabs_assemblies/v2_assemblies/ah_hkg_v1_centromeres.bed', sep="\t", header=None)
ns = pd.read_csv('mabs_assemblies/mapping/ah_hkg_Ns.bed', sep="\t", header=None)
for i in range(1,9):
    ax = fig.add_subplot(n_samp,8,n*8+i)
    sc = "chr" + str(i)
    chr8_rep = te8_rep[ (te8_rep.iloc[:,0].str.contains(sc))]
    
    plt.eventplot(chr8_rep[chr8_rep.iloc[:,5] == "pAge1"].iloc[:,1], linewidths=0.01, orientation='horizontal', colors="#2CA030FF")
    plt.eventplot(chr8_rep[chr8_rep.iloc[:,5] == "pAge2"].iloc[:,1], linewidths=0.01, orientation='horizontal', colors="#FFAA0EFF")
    plt.eventplot(chr8_rep[chr8_rep.iloc[:,5] == "pAa"].iloc[:,1], linewidths=0.01, orientation='horizontal', colors="#6F63BBFF")
    ns_sc = ns[ns.iloc[:,0] == sc]
    ax.axes.xaxis.set_visible(False)
    ax.axes.yaxis.set_visible(False)
    plt.eventplot(ns_sc.iloc[:,1], orientation='horizontal', colors="black", lineoffsets=0.2, linewidths=0.5)
    plt.eventplot(ns_sc.iloc[:,1], orientation='horizontal', colors="black", lineoffsets=1.8, linewidths=0.5)
    middle = (cen_pos.iloc[i-1,2] + cen_pos.iloc[i-1,1])/2
    plt.xlim([middle - 5000000, middle + 5000000])
    plt.ylim([0.1, 1.9])

n=22
te8_rep = pd.read_csv("mabs_assemblies/final_trash/LPT_v1_from_TE8_typen.bed", sep="\t", header=None)
cen_pos = pd.read_csv('mabs_assemblies/v2_assemblies/expand_centromeres/LPT_v1.bed', sep="\t", header=None)
ns = pd.read_csv('mabs_assemblies/mapping/LPT_v1_Ns.bed', sep="\t", header=None)
for i in range(1,9):
    ax = fig.add_subplot(n_samp,8,n*8+i)
    sc = "scaffold_" + str(i)
    chr8_rep = te8_rep[ (te8_rep.iloc[:,0].str.contains(sc))]
    plt.eventplot(chr8_rep[chr8_rep.iloc[:,4] == "pAa"].iloc[:,1], alpha=0.5, orientation='horizontal', colors="#6F63BBFF")
    plt.eventplot(chr8_rep[chr8_rep.iloc[:,4] == "pAge2"].iloc[:,1], alpha=0.5, orientation='horizontal', colors="#FFAA0EFF")
    
    plt.eventplot(chr8_rep[chr8_rep.iloc[:,4] == "pAge1"].iloc[:,1], alpha=0.5, orientation='horizontal', colors="#2CA030FF")
    ax.axes.xaxis.set_visible(False)
    ax.axes.yaxis.set_visible(False)
    ns_sc = ns[ns.iloc[:,0].str.contains(sc)]
    plt.eventplot(ns_sc.iloc[:,1], orientation='horizontal', colors="black", lineoffsets=0.2, linewidths=0.5)
    plt.eventplot(ns_sc.iloc[:,1], orientation='horizontal', colors="black", lineoffsets=1.8, linewidths=0.5)
    middle = (cen_pos.iloc[i-1,2] + cen_pos.iloc[i-1,1])/2
    plt.xlim([middle - 5000000, middle + 5000000])
    plt.ylim([0.1, 1.9])
    
n=23
te8_rep = pd.read_csv("mabs_assemblies/final_trash/PTP_v0_from_TE8_typen.bed", sep="\t", header=None)
cen_pos = pd.read_csv('mabs_assemblies/v2_assemblies/expand_centromeres/PTP_v0.bed', sep="\t", header=None)
ns = pd.read_csv('mabs_assemblies/mapping/PTP_v0_Ns.bed', sep="\t", header=None)
for i in range(1,9):
    ax = fig.add_subplot(n_samp,8,n*8+i)
    sc = "scaffold_" + str(i)
    chr8_rep = te8_rep[ (te8_rep.iloc[:,0].str.contains(sc))]
    
    plt.eventplot(chr8_rep[chr8_rep.iloc[:,4] == "pAge1"].iloc[:,1], alpha=0.5, orientation='horizontal', colors="#2CA030FF")
    plt.eventplot(chr8_rep[chr8_rep.iloc[:,4] == "pAge2"].iloc[:,1], alpha=0.5, orientation='horizontal', colors="#FFAA0EFF")
    plt.eventplot(chr8_rep[chr8_rep.iloc[:,4] == "pAa"].iloc[:,1], alpha=0.5, orientation='horizontal', colors="#6F63BBFF")
    ax.axes.xaxis.set_visible(False)
    ax.axes.yaxis.set_visible(False)
    ns_sc = ns[ns.iloc[:,0].str.contains(sc)]
    plt.eventplot(ns_sc.iloc[:,1], orientation='horizontal', colors="black", lineoffsets=0.2, linewidths=0.5)
    plt.eventplot(ns_sc.iloc[:,1], orientation='horizontal', colors="black", lineoffsets=1.8, linewidths=0.5)
    middle = (cen_pos.iloc[i-1,2] + cen_pos.iloc[i-1,1])/2
    plt.xlim([middle - 5000000, middle + 5000000])
    plt.ylim([0.1, 1.9])
n=15
te8_rep = pd.read_csv("mabs_assemblies/final_trash/Pais09_v1.2_typen.bed", sep="\t", header=None)
cen_pos = pd.read_csv('mabs_assemblies/v2_assemblies/Pais09_merged_centromeres.bed', sep="\t", header=None)
ns = pd.read_csv('mabs_assemblies/mapping/Pais09_Ns.bed', sep="\t", header=None)
for i in range(1,9):
    ax = fig.add_subplot(n_samp,8,n*8+i)
    sc = "Chr" + str(i)
    chr8_rep = te8_rep[ (te8_rep.iloc[:,0].str.contains(sc))]
    
    plt.eventplot(chr8_rep[chr8_rep.iloc[:,3] == "pAge1"].iloc[:,1], linewidths=0.01, orientation='horizontal', colors="#2CA030FF")
    plt.eventplot(chr8_rep[chr8_rep.iloc[:,3] == "pAge2"].iloc[:,1], linewidths=0.01, orientation='horizontal', colors="#FFAA0EFF")
    plt.eventplot(chr8_rep[chr8_rep.iloc[:,3] == "pAa"].iloc[:,1], linewidths=0.01, orientation='horizontal', colors="#6F63BBFF")
    ns_sc = ns[ns.iloc[:,0] == sc]
    ax.axes.xaxis.set_visible(False)
    ax.axes.yaxis.set_visible(False)
    plt.eventplot(ns_sc.iloc[:,1], orientation='horizontal', colors="black", lineoffsets=0.2, linewidths=0.5)
    plt.eventplot(ns_sc.iloc[:,1], orientation='horizontal', colors="black", lineoffsets=1.8, linewidths=0.5)
    middle = (cen_pos.iloc[i-1,2] + cen_pos.iloc[i-1,1])/2
    plt.xlim([middle - 5000000, middle + 5000000])
    plt.ylim([0.1, 1.9])
n=35
te8_rep = pd.read_csv("mabs_assemblies/final_trash/TE11_pca_typen.bed", sep="\t", header=None)
cen_pos = pd.read_csv('mabs_assemblies/v2_assemblies/expand_centromeres/TE11_v2.bed', sep="\t", header=None)
ns = pd.read_csv('mabs_assemblies/mapping/TE11_v2_Ns.bed', sep="\t", header=None)
for i in range(1,9):
    ax = fig.add_subplot(n_samp,8,n*8+i)
    sc = "scaffold_" + str(i)
    chr8_rep = te8_rep[ (te8_rep.iloc[:,0].str.contains(sc))]
    plt.eventplot(chr8_rep[chr8_rep.iloc[:,3] == "pAge1"].iloc[:,1], linewidths=0.01, orientation='horizontal', colors="#2CA030FF")
    plt.eventplot(chr8_rep[chr8_rep.iloc[:,3] == "pAge2"].iloc[:,1], linewidths=0.01, orientation='horizontal', colors="#FFAA0EFF")
    plt.eventplot(chr8_rep[chr8_rep.iloc[:,3] == "pAa"].iloc[:,1], linewidths=0.01, orientation='horizontal', colors="#6F63BBFF")
    ax.axes.xaxis.set_visible(False)
    ax.axes.yaxis.set_visible(False)
    ns_sc = ns[ns.iloc[:,0].str.contains(sc)]
    plt.eventplot(ns_sc.iloc[:,1], orientation='horizontal', colors="black", lineoffsets=0.2, linewidths=0.5)
    plt.eventplot(ns_sc.iloc[:,1], orientation='horizontal', colors="black", lineoffsets=1.8, linewidths=0.5)
    middle = (cen_pos.iloc[i-1,2] + cen_pos.iloc[i-1,1])/2
    plt.xlim([middle - 5000000, middle + 5000000])
    plt.ylim([0.1, 1.9])
    

n=34
te8_rep = pd.read_csv("mabs_assemblies/final_trash/TE4_pca_typen.bed", sep="\t", header=None)
cen_pos = pd.read_csv('mabs_assemblies/v2_assemblies/expand_centromeres/TE4_v1.bed', sep="\t", header=None)
ns = pd.read_csv('mabs_assemblies/mapping/TE4_v1_Ns.bed', sep="\t", header=None)
for i in range(1,9):
    ax = fig.add_subplot(n_samp,8,n*8+i)
    sc = "scaffold_" + str(i)
    chr8_rep = te8_rep[ (te8_rep.iloc[:,0].str.contains(sc))]
    plt.eventplot(chr8_rep[chr8_rep.iloc[:,3] == "pAge1"].iloc[:,1], linewidths=0.01, orientation='horizontal', colors="#2CA030FF")
    plt.eventplot(chr8_rep[chr8_rep.iloc[:,3] == "pAge2"].iloc[:,1], linewidths=0.01, orientation='horizontal', colors="#FFAA0EFF")
    plt.eventplot(chr8_rep[chr8_rep.iloc[:,3] == "pAa"].iloc[:,1], linewidths=0.01, orientation='horizontal', colors="#6F63BBFF")
    ax.axes.xaxis.set_visible(False)
    ax.axes.yaxis.set_visible(False)
    ns_sc = ns[ns.iloc[:,0].str.contains(sc)]
    plt.eventplot(ns_sc.iloc[:,1], orientation='horizontal', colors="black", lineoffsets=0.2, linewidths=0.5)
    plt.eventplot(ns_sc.iloc[:,1], orientation='horizontal', colors="black", lineoffsets=1.8, linewidths=0.5)
    middle = (cen_pos.iloc[i-1,2] + cen_pos.iloc[i-1,1])/2
    plt.xlim([middle - 5000000, middle + 5000000])
    plt.ylim([0.1, 1.9])
    
n=33
te8_rep = pd.read_csv("mabs_assemblies/final_trash/TE8_pca_typen.bed", sep="\t", header=None)
cen_pos = pd.read_csv('mabs_assemblies/v2_assemblies/expand_centromeres/TE8_v1.bed', sep="\t", header=None)
ns = pd.read_csv('mabs_assemblies/mapping/TE8_v1_Ns.bed', sep="\t", header=None)
for i in range(1,9):
    ax = fig.add_subplot(n_samp,8,n*8+i)
    sc = "scaffold_" + str(i)
    chr8_rep = te8_rep[ (te8_rep.iloc[:,0].str.contains(sc))]
    plt.eventplot(chr8_rep[chr8_rep.iloc[:,3] == "pAge1"].iloc[:,1], linewidths=0.01, orientation='horizontal', colors="#2CA030FF")
    plt.eventplot(chr8_rep[chr8_rep.iloc[:,3] == "pAge2"].iloc[:,1], linewidths=0.01, orientation='horizontal', colors="#FFAA0EFF")
    plt.eventplot(chr8_rep[chr8_rep.iloc[:,3] == "pAa"].iloc[:,1], linewidths=0.01, orientation='horizontal', colors="#6F63BBFF")
    ax.axes.xaxis.set_visible(False)
    ax.axes.yaxis.set_visible(False)
    ns_sc = ns[ns.iloc[:,0].str.contains(sc)]
    plt.eventplot(ns_sc.iloc[:,1], orientation='horizontal', colors="black", lineoffsets=0.2, linewidths=0.5)
    plt.eventplot(ns_sc.iloc[:,1], orientation='horizontal', colors="black", lineoffsets=1.8, linewidths=0.5)
    middle = (cen_pos.iloc[i-1,2] + cen_pos.iloc[i-1,1])/2
    plt.xlim([middle - 5000000, middle + 5000000])
    plt.ylim([0.1, 1.9])

n=16
te8_rep = pd.read_csv("mabs_assemblies/final_trash/Wall10_v1.2_typen.bed", sep="\t", header=None)
cen_pos = pd.read_csv('mabs_assemblies/v2_assemblies/Wall10_merged_centromeres.bed', sep="\t", header=None)
ns = pd.read_csv('mabs_assemblies/mapping/Wall10_Ns.bed', sep="\t", header=None)
for i in range(1,9):
    ax = fig.add_subplot(n_samp,8,n*8+i)
    sc = "Chr" + str(i)
    chr8_rep = te8_rep[ (te8_rep.iloc[:,0].str.contains(sc))]
    
    plt.eventplot(chr8_rep[chr8_rep.iloc[:,3] == "pAge1"].iloc[:,1], linewidths=0.01, orientation='horizontal', colors="#2CA030FF")
    plt.eventplot(chr8_rep[chr8_rep.iloc[:,3] == "pAge2"].iloc[:,1], linewidths=0.01, orientation='horizontal', colors="#FFAA0EFF")
    plt.eventplot(chr8_rep[chr8_rep.iloc[:,3] == "pAa"].iloc[:,1], linewidths=0.01, orientation='horizontal', colors="#6F63BBFF")
    ns_sc = ns[ns.iloc[:,0] == sc]
    ax.axes.xaxis.set_visible(False)
    ax.axes.yaxis.set_visible(False)
    plt.eventplot(ns_sc.iloc[:,1], orientation='horizontal', colors="black", lineoffsets=0.2, linewidths=0.5)
    plt.eventplot(ns_sc.iloc[:,1], orientation='horizontal', colors="black", lineoffsets=1.8, linewidths=0.5)
    middle = (cen_pos.iloc[i-1,2] + cen_pos.iloc[i-1,1])/2
    plt.xlim([middle - 5000000, middle + 5000000])
    plt.ylim([0.1, 1.9])
n=17
te8_rep = pd.read_csv("mabs_assemblies/final_trash/Ahal_1.3_pca_typen.bed", sep="\t", header=None)
cen_pos = pd.read_csv('mabs_assemblies/v2_assemblies/Ahal_1.3_rev_merge_centromeres.bed', sep="\t", header=None)
ns = pd.read_csv('mabs_assemblies/mapping/Ahal_1.3_rev_Ns.bed', sep="\t", header=None)
for i in range(1,9):
    ax = fig.add_subplot(n_samp,8,n*8+i)
    sc = "scaffold_" + str(i)
    chr8_rep = te8_rep[ (te8_rep.iloc[:,0].str.contains(sc))]
    
    plt.eventplot(chr8_rep[chr8_rep.iloc[:,3] == "pAge1"].iloc[:,1], linewidths=0.01, orientation='horizontal', colors="#2CA030FF")
    plt.eventplot(chr8_rep[chr8_rep.iloc[:,3] == "pAge2"].iloc[:,1], linewidths=0.01, orientation='horizontal', colors="#FFAA0EFF")
    plt.eventplot(chr8_rep[chr8_rep.iloc[:,3] == "pAa"].iloc[:,1], linewidths=0.01, orientation='horizontal', colors="#6F63BBFF")
    ns_sc = ns[ns.iloc[:,0] == sc]
    ax.axes.xaxis.set_visible(False)
    ax.axes.yaxis.set_visible(False)
    plt.eventplot(ns_sc.iloc[:,1], orientation='horizontal', colors="black", lineoffsets=0.2, linewidths=0.5)
    plt.eventplot(ns_sc.iloc[:,1], orientation='horizontal', colors="black", lineoffsets=1.8, linewidths=0.5)
    middle = (cen_pos.iloc[i-1,2] + cen_pos.iloc[i-1,1])/2
    plt.xlim([middle - 5000000, middle + 5000000])
    plt.ylim([0.1, 1.9])
n=10
te8_rep = pd.read_csv("mabs_assemblies/final_trash/AS150_v6_AA_pAa_typen.bed", sep="\t", header=None)
cen_pos = pd.read_csv('mabs_assemblies/v2_assemblies/AS150_v6_AA_merged_centromeres.bed', sep="\t", header=None)
ns = pd.read_csv('mabs_assemblies/mapping/AS150_v6_AA_Ns.bed', sep="\t", header=None)
for i in range(1,9):
    ax = fig.add_subplot(n_samp,8,n*8+i)
    sc = "Chr" + str(i+5)
    chr8_rep = te8_rep[ (te8_rep.iloc[:,0].str.contains(sc))]
    plt.eventplot(chr8_rep[chr8_rep.iloc[:,3] == "pAge2"].iloc[:,1], linewidths=0.01, orientation='horizontal', colors="#FFAA0EFF")
    plt.eventplot(chr8_rep[chr8_rep.iloc[:,3] == "pAa"].iloc[:,1], linewidths=0.01, orientation='horizontal', colors="#6F63BBFF")
    plt.eventplot(chr8_rep[chr8_rep.iloc[:,3] == "pAge1"].iloc[:,1], linewidths=0.01, orientation='horizontal', colors="#2CA030FF")
    ns_sc = ns[ns.iloc[:,0] == sc]
    plt.eventplot(chr8_rep[chr8_rep.iloc[:,3] == "AthCEN159"].iloc[:,1], linewidths=0.01, orientation='horizontal', colors="#C7519CFF")
    plt.eventplot(chr8_rep[chr8_rep.iloc[:,3] == "AthCEN178"].iloc[:,1], linewidths=0.01, orientation='horizontal', colors="#12A2A8FF")
    ax.axes.xaxis.set_visible(False)
    ax.axes.yaxis.set_visible(False)
    plt.eventplot(ns_sc.iloc[:,1], orientation='horizontal', colors="black", lineoffsets=0.2, linewidths=0.5)
    plt.eventplot(ns_sc.iloc[:,1], orientation='horizontal', colors="black", lineoffsets=1.8, linewidths=0.5)
    middle = (cen_pos.iloc[i-1,2] + cen_pos.iloc[i-1,1])/2
    plt.xlim([middle - 5000000, middle + 5000000])
    plt.ylim([0.1, 1.9])
n=11
te8_rep = pd.read_csv("mabs_assemblies/final_trash/AS530_v6_AA_pAa_typen.bed", sep="\t", header=None)
cen_pos = pd.read_csv('mabs_assemblies/v2_assemblies/AS530_v6_AA_merged_centromeres.bed', sep="\t", header=None)
ns = pd.read_csv('mabs_assemblies/mapping/AS530_v6_AA_Ns.bed', sep="\t", header=None)
for i in range(1,9):
    ax = fig.add_subplot(n_samp,8,n*8+i)
    sc = "Chr" + str(i+5)
    chr8_rep = te8_rep[ (te8_rep.iloc[:,0].str.contains(sc))]
    plt.eventplot(chr8_rep[chr8_rep.iloc[:,3] == "pAge2"].iloc[:,1], linewidths=0.01, orientation='horizontal', colors="#FFAA0EFF")
    plt.eventplot(chr8_rep[chr8_rep.iloc[:,3] == "pAa"].iloc[:,1], linewidths=0.01, orientation='horizontal', colors="#6F63BBFF")
    plt.eventplot(chr8_rep[chr8_rep.iloc[:,3] == "pAge1"].iloc[:,1], linewidths=0.01, orientation='horizontal', colors="#2CA030FF")
    ns_sc = ns[ns.iloc[:,0] == sc]
    plt.eventplot(chr8_rep[chr8_rep.iloc[:,3] == "AthCEN159"].iloc[:,1], linewidths=0.01, orientation='horizontal', colors="#C7519CFF")
    plt.eventplot(chr8_rep[chr8_rep.iloc[:,3] == "AthCEN178"].iloc[:,1], linewidths=0.01, orientation='horizontal', colors="#12A2A8FF")
    ax.axes.xaxis.set_visible(False)
    ax.axes.yaxis.set_visible(False)
    plt.eventplot(ns_sc.iloc[:,1], orientation='horizontal', colors="black", lineoffsets=0.2, linewidths=0.5)
    plt.eventplot(ns_sc.iloc[:,1], orientation='horizontal', colors="black", lineoffsets=1.8, linewidths=0.5)
    middle = (cen_pos.iloc[i-1,2] + cen_pos.iloc[i-1,1])/2
    plt.xlim([middle - 5000000, middle + 5000000])
    plt.ylim([0.1, 1.9])
n=9
te8_rep = pd.read_csv("mabs_assemblies/final_trash/ASS3_v6_AA_pAa_typen.bed", sep="\t", header=None)
cen_pos = pd.read_csv('mabs_assemblies/v2_assemblies/ASS3_v6_AA_merged_centromeres.bed', sep="\t", header=None)
ns = pd.read_csv('mabs_assemblies/mapping/ASS3_v6_AA_Ns.bed', sep="\t", header=None)
for i in range(1,9):
    ax = fig.add_subplot(n_samp,8,n*8+i)
    sc = "Chr" + str(i+5)
    chr8_rep = te8_rep[ (te8_rep.iloc[:,0].str.contains(sc))]
    
    
    plt.eventplot(chr8_rep[chr8_rep.iloc[:,3] == "pAge2"].iloc[:,1], linewidths=0.01, orientation='horizontal', colors="#FFAA0EFF")
    plt.eventplot(chr8_rep[chr8_rep.iloc[:,3] == "pAa"].iloc[:,1], linewidths=0.01, orientation='horizontal', colors="#6F63BBFF")
    plt.eventplot(chr8_rep[chr8_rep.iloc[:,3] == "pAge1"].iloc[:,1], linewidths=0.01, orientation='horizontal', colors="#2CA030FF")
    ns_sc = ns[ns.iloc[:,0] == sc]
    plt.eventplot(chr8_rep[chr8_rep.iloc[:,3] == "AthCEN159"].iloc[:,1], linewidths=0.01, orientation='horizontal', colors="#C7519CFF")
    plt.eventplot(chr8_rep[chr8_rep.iloc[:,3] == "AthCEN178"].iloc[:,1], linewidths=0.01, orientation='horizontal', colors="#12A2A8FF")
    ax.axes.xaxis.set_visible(False)
    ax.axes.yaxis.set_visible(False)
    plt.eventplot(ns_sc.iloc[:,1], orientation='horizontal', colors="black", lineoffsets=0.2, linewidths=0.5)
    plt.eventplot(ns_sc.iloc[:,1], orientation='horizontal', colors="black", lineoffsets=1.8, linewidths=0.5)
    middle = (cen_pos.iloc[i-1,2] + cen_pos.iloc[i-1,1])/2
    plt.xlim([middle - 5000000, middle + 5000000])
    plt.ylim([0.1, 1.9])
n=6
te8_rep = pd.read_csv("mabs_assemblies/final_trash/Croatica_from_TE8_typen.bed", sep="\t", header=None)
cen_pos = pd.read_csv('mabs_assemblies/v2_assemblies/Croatica_v1_merged_centromeres.bed', sep="\t", header=None)
ns = pd.read_csv('mabs_assemblies/mapping/croa_i4_Ns.bed', sep="\t", header=None)
for i in range(1,9):
    ax = fig.add_subplot(n_samp,8,n*8+i)
    sc = "scaffold_" + str(i)
    chr8_rep = te8_rep[ (te8_rep.iloc[:,0].str.contains(sc))]
    plt.eventplot(chr8_rep[chr8_rep.iloc[:,4] == "pAge2"].iloc[:,1], linewidths=0.01, orientation='horizontal', colors="#FFAA0EFF")
    plt.eventplot(chr8_rep[chr8_rep.iloc[:,4] == "pAa"].iloc[:,1], linewidths=0.01, orientation='horizontal', colors="#6F63BBFF")
    plt.eventplot(chr8_rep[chr8_rep.iloc[:,4] == "pAge1"].iloc[:,1], linewidths=0.01, orientation='horizontal', colors="#2CA030FF")
    ns_sc = ns[ns.iloc[:,0] == sc]
    plt.eventplot(chr8_rep[chr8_rep.iloc[:,4] == "AthCEN159"].iloc[:,1], linewidths=0.01, orientation='horizontal', colors="#C7519CFF")
    ax.axes.xaxis.set_visible(False)
    ax.axes.yaxis.set_visible(False)
    plt.eventplot(ns_sc.iloc[:,1], orientation='horizontal', colors="black", lineoffsets=0.2, linewidths=0.5)
    plt.eventplot(ns_sc.iloc[:,1], orientation='horizontal', colors="black", lineoffsets=1.8, linewidths=0.5)
    middle = (cen_pos.iloc[i-1,2] + cen_pos.iloc[i-1,1])/2
    plt.xlim([middle - 5000000, middle + 5000000])
    plt.ylim([0.1, 1.9])
n=7
te8_rep = pd.read_csv("mabs_assemblies/final_trash/BaltAre_pca_typen.bed", sep="\t", header=None)
cen_pos = pd.read_csv('mabs_assemblies/v2_assemblies/BaltAre_v1_merge_centromeres.bed', sep="\t", header=None)
ns = pd.read_csv('mabs_assemblies/mapping/BaltAre_v1_Ns.bed', sep="\t", header=None)
for i in range(1,9):
    ax = fig.add_subplot(n_samp,8,n*8+i)
    sc = "scaffold_" + str(i)
    chr8_rep = te8_rep[ (te8_rep.iloc[:,0].str.contains(sc))]
    plt.eventplot(chr8_rep[chr8_rep.iloc[:,3] == "pAge1"].iloc[:,1], linewidths=0.01, orientation='horizontal', colors="#2CA030FF")
    plt.eventplot(chr8_rep[chr8_rep.iloc[:,3] == "pAge2"].iloc[:,1], linewidths=0.01, orientation='horizontal', colors="#FFAA0EFF")
    plt.eventplot(chr8_rep[chr8_rep.iloc[:,3] == "pAa"].iloc[:,1], linewidths=0.01, orientation='horizontal', colors="#6F63BBFF")
    ns_sc = ns[ns.iloc[:,0] == sc]
    plt.eventplot(chr8_rep[chr8_rep.iloc[:,3] == "AthCEN159"].iloc[:,1], linewidths=0.01, orientation='horizontal', colors="#C7519CFF")
    ax.axes.xaxis.set_visible(False)
    ax.axes.yaxis.set_visible(False)
    plt.eventplot(ns_sc.iloc[:,1], orientation='horizontal', colors="black", lineoffsets=0.2, linewidths=0.5)
    plt.eventplot(ns_sc.iloc[:,1], orientation='horizontal', colors="black", lineoffsets=1.8, linewidths=0.5)
    middle = (cen_pos.iloc[i-1,2] + cen_pos.iloc[i-1,1])/2
    plt.xlim([middle - 5000000, middle + 5000000])
    plt.ylim([0.1, 1.9])
n=14
te8_rep = pd.read_csv("mabs_assemblies/final_trash/Ceb_c_from_TE8_typen.bed", sep="\t", header=None)
cen_pos = pd.read_csv('mabs_assemblies/v2_assemblies/expand_centromeres/Ceb_c_v1.bed', sep="\t", header=None)
ns = pd.read_csv('mabs_assemblies/mapping/Ceb_c_v1_Ns.bed', sep="\t", header=None)
for i in range(1,9):
    ax = fig.add_subplot(n_samp,8,n*8+i)
    sc = "scaffold_" + str(i)
    chr8_rep = te8_rep[ (te8_rep.iloc[:,0].str.contains(sc))]
    plt.eventplot(chr8_rep[chr8_rep.iloc[:,4] == "pAge1"].iloc[:,1], linewidths=0.01, orientation='horizontal', colors="#2CA030FF")
    plt.eventplot(chr8_rep[chr8_rep.iloc[:,4] == "pAa"].iloc[:,1], linewidths=0.01, orientation='horizontal', colors="#6F63BBFF")
    plt.eventplot(chr8_rep[chr8_rep.iloc[:,4] == "pAge2"].iloc[:,1], linewidths=0.01, orientation='horizontal', colors="#FFAA0EFF")
    plt.eventplot(chr8_rep[chr8_rep.iloc[:,4] == "AthCEN159"].iloc[:,1], linewidths=0.01, orientation='horizontal', colors="#C7519CFF")
    ax.axes.xaxis.set_visible(False)
    ax.axes.yaxis.set_visible(False)
    ns_sc = ns[ns.iloc[:,0].str.contains(sc)]
    plt.eventplot(ns_sc.iloc[:,1], orientation='horizontal', colors="black", lineoffsets=0.2, linewidths=0.5)
    plt.eventplot(ns_sc.iloc[:,1], orientation='horizontal', colors="black", lineoffsets=1.8, linewidths=0.5)
    middle = (cen_pos.iloc[i-1,2] + cen_pos.iloc[i-1,1])/2
    plt.xlim([middle - 5000000, middle + 5000000])
    plt.ylim([0.1, 1.9])

    
n=13


te8_rep = pd.read_csv("mabs_assemblies/final_trash/Ceb_d_from_TE8_typen.bed", sep="\t", header=None)
cen_pos = pd.read_csv('mabs_assemblies/v2_assemblies/expand_centromeres/Ceb_d_v2.bed', sep="\t", header=None)
ns = pd.read_csv('mabs_assemblies/mapping/Ceb_d_v2_Ns.bed', sep="\t", header=None)
for i in range(1,9):
    ax = fig.add_subplot(n_samp,8,n*8+i)
    sc = "scaffold_" + str(i)
    chr8_rep = te8_rep[ (te8_rep.iloc[:,0].str.contains(sc))]
    plt.eventplot(chr8_rep[chr8_rep.iloc[:,4] == "pAge1"].iloc[:,1], linewidths=0.01, orientation='horizontal', colors="#2CA030FF")
    plt.eventplot(chr8_rep[chr8_rep.iloc[:,4] == "pAa"].iloc[:,1], linewidths=0.01, orientation='horizontal', colors="#6F63BBFF")
    plt.eventplot(chr8_rep[chr8_rep.iloc[:,4] == "pAge2"].iloc[:,1], linewidths=0.01, orientation='horizontal', colors="#FFAA0EFF")
    plt.eventplot(chr8_rep[chr8_rep.iloc[:,4] == "AthCEN159"].iloc[:,1], linewidths=0.01, orientation='horizontal', colors="#C7519CFF")
    ax.axes.xaxis.set_visible(False)
    ax.axes.yaxis.set_visible(False)
    ns_sc = ns[ns.iloc[:,0].str.contains(sc)]
    plt.eventplot(ns_sc.iloc[:,1], orientation='horizontal', colors="black", lineoffsets=0.2, linewidths=0.5)
    plt.eventplot(ns_sc.iloc[:,1], orientation='horizontal', colors="black", lineoffsets=1.8, linewidths=0.5)
    middle = (cen_pos.iloc[i-1,2] + cen_pos.iloc[i-1,1])/2
    plt.xlim([middle - 5000000, middle + 5000000])
    plt.ylim([0.1, 1.9])

n=12


te8_rep = pd.read_csv("mabs_assemblies/final_trash/Ped_v2_TE8_typen.bed", sep="\t", header=None)
cen_pos = pd.read_csv('mabs_assemblies/v2_assemblies/expand_centromeres/Ped_v2.bed', sep="\t", header=None)
ns = pd.read_csv('mabs_assemblies/mapping/Ped_v2_Ns.bed', sep="\t", header=None)
for i in range(1,9):
    ax = fig.add_subplot(n_samp,8,n*8+i)
    sc = "scaffold_" + str(i)
    chr8_rep = te8_rep[ (te8_rep.iloc[:,0].str.contains(sc))]
    plt.eventplot(chr8_rep[chr8_rep.iloc[:,4] == "pAge1"].iloc[:,1], linewidths=0.01, orientation='horizontal', colors="#2CA030FF")
    plt.eventplot(chr8_rep[chr8_rep.iloc[:,4] == "pAa"].iloc[:,1], linewidths=0.01, orientation='horizontal', colors="#6F63BBFF")
    plt.eventplot(chr8_rep[chr8_rep.iloc[:,4] == "pAge2"].iloc[:,1], linewidths=0.01, orientation='horizontal', colors="#FFAA0EFF")
    plt.eventplot(chr8_rep[chr8_rep.iloc[:,4] == "AthCEN159"].iloc[:,1], linewidths=0.01, orientation='horizontal', colors="#C7519CFF")
    ax.axes.xaxis.set_visible(False)
    ax.axes.yaxis.set_visible(False)
    ns_sc = ns[ns.iloc[:,0].str.contains(sc)]
    plt.eventplot(ns_sc.iloc[:,1], orientation='horizontal', colors="black", lineoffsets=0.2, linewidths=0.5)
    plt.eventplot(ns_sc.iloc[:,1], orientation='horizontal', colors="black", lineoffsets=1.8, linewidths=0.5)
    middle = (cen_pos.iloc[i-1,2] + cen_pos.iloc[i-1,1])/2
    plt.xlim([middle - 5000000, middle + 5000000])
    plt.ylim([0.1, 1.9])

n=24
te8_rep = pd.read_csv("mabs_assemblies/final_trash/TSS_v1_from_TE8_typen.bed", sep="\t", header=None)
cen_pos = pd.read_csv('mabs_assemblies/v2_assemblies/expand_centromeres/TSS_v1.bed', sep="\t", header=None)
ns = pd.read_csv('mabs_assemblies/mapping/TSS_v1_Ns.bed', sep="\t", header=None)
for i in range(1,9):
    ax = fig.add_subplot(n_samp,8,n*8+i)
    sc = "scaffold_" + str(i)
    chr8_rep = te8_rep[ (te8_rep.iloc[:,0].str.contains(sc))]
    plt.eventplot(chr8_rep[chr8_rep.iloc[:,4] == "pAge2"].iloc[:,1], alpha=0.5, orientation='horizontal', colors="#FFAA0EFF")
    plt.eventplot(chr8_rep[chr8_rep.iloc[:,4] == "pAa"].iloc[:,1], alpha=0.5, orientation='horizontal', colors="#6F63BBFF")
    plt.eventplot(chr8_rep[chr8_rep.iloc[:,4] == "pAge1"].iloc[:,1], alpha=0.5, orientation='horizontal', colors="#2CA030FF")
    ax.axes.xaxis.set_visible(False)
    ax.axes.yaxis.set_visible(False)
    ns_sc = ns[ns.iloc[:,0].str.contains(sc)]
    plt.eventplot(ns_sc.iloc[:,1], orientation='horizontal', colors="black", lineoffsets=0.2, linewidths=0.5)
    plt.eventplot(ns_sc.iloc[:,1], orientation='horizontal', colors="black", lineoffsets=1.8, linewidths=0.5)
    middle = (cen_pos.iloc[i-1,2] + cen_pos.iloc[i-1,1])/2
    plt.xlim([middle - 5000000, middle + 5000000])
    plt.ylim([0.1, 1.9])

n=47


te8_rep = pd.read_csv("mabs_assemblies/final_trash/BOR_v1_from_TE8_typen.bed", sep="\t", header=None)
cen_pos = pd.read_csv('mabs_assemblies/v2_assemblies/expand_centromeres/BOR_v1.bed', sep="\t", header=None)
ns = pd.read_csv('mabs_assemblies/mapping/BOR_v1_Ns.bed', sep="\t", header=None)
for i in range(1,9):
    ax = fig.add_subplot(n_samp,8,n*8+i)
    sc = "scaffold_" + str(i)
    chr8_rep = te8_rep[ (te8_rep.iloc[:,0].str.contains(sc))]
    plt.eventplot(chr8_rep[chr8_rep.iloc[:,4] == "pAge1"].iloc[:,1], linewidths=0.01, orientation='horizontal', colors="#2CA030FF")
    plt.eventplot(chr8_rep[chr8_rep.iloc[:,4] == "pAge2"].iloc[:,1], linewidths=0.01, orientation='horizontal', colors="#FFAA0EFF")
    plt.eventplot(chr8_rep[chr8_rep.iloc[:,4] == "pAa"].iloc[:,1], linewidths=0.01, orientation='horizontal', colors="#6F63BBFF")
    plt.eventplot(chr8_rep[chr8_rep.iloc[:,4] == "AthCEN159"].iloc[:,1], linewidths=0.01, orientation='horizontal', colors="#C7519CFF")
    ax.axes.xaxis.set_visible(False)
    ax.axes.yaxis.set_visible(False)
    ns_sc = ns[ns.iloc[:,0].str.contains(sc)]
    plt.eventplot(ns_sc.iloc[:,1], orientation='horizontal', colors="black", lineoffsets=0.2, linewidths=0.5)
    plt.eventplot(ns_sc.iloc[:,1], orientation='horizontal', colors="black", lineoffsets=1.8, linewidths=0.5)
    middle = (cen_pos.iloc[i-1,2] + cen_pos.iloc[i-1,1])/2
    plt.xlim([middle - 5000000, middle + 5000000])
    plt.ylim([0.1, 1.9])

n=52
te8_rep = pd.read_csv("mabs_assemblies/final_trash/AL57_pca_typen_reclass.bed", sep="\t", header=None)
cen_pos = pd.read_csv('mabs_assemblies/v2_assemblies/AL57_v1_centromeres.bed', sep="\t", header=None)
ns = pd.read_csv('mabs_assemblies/mapping/AL57_v1_Ns.bed', sep="\t", header=None)
for i in range(1,9):
    ax = fig.add_subplot(n_samp,8,n*8+i)
    sc = "scaffold_" + str(i)
    chr8_rep = te8_rep[ (te8_rep.iloc[:,0].str.contains(sc))]
    plt.eventplot(chr8_rep[chr8_rep.iloc[:,3] == "pAge1"].iloc[:,1], linewidths=0.01, orientation='horizontal', colors="#2CA030FF")
    plt.eventplot(chr8_rep[chr8_rep.iloc[:,3] == "pAge2"].iloc[:,1], linewidths=0.01, orientation='horizontal', colors="#FFAA0EFF")
    plt.eventplot(chr8_rep[chr8_rep.iloc[:,3] == "pAa"].iloc[:,1], linewidths=0.01, orientation='horizontal', colors="#6F63BBFF")
    plt.eventplot(chr8_rep[chr8_rep.iloc[:,3] == "AthCEN159"].iloc[:,1], linewidths=0.01, orientation='horizontal', colors="#C7519CFF")
    ax.axes.xaxis.set_visible(False)
    ax.axes.yaxis.set_visible(False)
    ns_sc = ns[ns.iloc[:,0].str.contains(sc)]
    plt.eventplot(ns_sc.iloc[:,1], orientation='horizontal', colors="black", lineoffsets=0.2, linewidths=0.5)
    plt.eventplot(ns_sc.iloc[:,1], orientation='horizontal', colors="black", lineoffsets=1.8, linewidths=0.5)
    middle = (cen_pos.iloc[i-1,2] + cen_pos.iloc[i-1,1])/2
    plt.xlim([middle - 5000000, middle + 5000000])
    plt.ylim([0.1, 1.9])
n=51
te8_rep = pd.read_csv("mabs_assemblies/final_trash/AL08_v1_from_TE8_typen.bed", sep="\t", header=None)
cen_pos = pd.read_csv('mabs_assemblies/v2_assemblies/expand_centromeres/AL08_v2.bed', sep="\t", header=None)
ns = pd.read_csv('mabs_assemblies/mapping/AL08_v2_Ns.bed', sep="\t", header=None)
for i in range(1,9):
    ax = fig.add_subplot(n_samp,8,n*8+i)
    sc = "scaffold_" + str(i)
    chr8_rep = te8_rep[ (te8_rep.iloc[:,0].str.contains(sc))]
    plt.eventplot(chr8_rep[chr8_rep.iloc[:,4] == "pAge1"].iloc[:,1], linewidths=0.01, orientation='horizontal', colors="#2CA030FF")
    plt.eventplot(chr8_rep[chr8_rep.iloc[:,4] == "pAge2"].iloc[:,1], linewidths=0.01, orientation='horizontal', colors="#FFAA0EFF")
    plt.eventplot(chr8_rep[chr8_rep.iloc[:,4] == "pAa"].iloc[:,1], linewidths=0.01, orientation='horizontal', colors="#6F63BBFF")
    plt.eventplot(chr8_rep[chr8_rep.iloc[:,4] == "AthCEN159"].iloc[:,1], linewidths=0.01, orientation='horizontal', colors="#C7519CFF")
    ax.axes.xaxis.set_visible(False)
    ax.axes.yaxis.set_visible(False)
    ns_sc = ns[ns.iloc[:,0].str.contains(sc)]
    plt.eventplot(ns_sc.iloc[:,1], orientation='horizontal', colors="black", lineoffsets=0.2, linewidths=0.5)
    plt.eventplot(ns_sc.iloc[:,1], orientation='horizontal', colors="black", lineoffsets=1.8, linewidths=0.5)
    middle = (cen_pos.iloc[i-1,2] + cen_pos.iloc[i-1,1])/2
    plt.xlim([middle - 5000000, middle + 5000000])
    plt.ylim([0.1, 1.9])
    
n=53
te8_rep = pd.read_csv("mabs_assemblies/final_trash/AL27_v1_from_TE8_typen.bed", sep="\t", header=None)
cen_pos = pd.read_csv('mabs_assemblies/v2_assemblies/expand_centromeres/AL27_v1.bed', sep="\t", header=None)
#ns = pd.read_csv('mabs_assemblies/mapping/AL27_v1_Ns.bed', sep="\t", header=None)
for i in range(1,9):
    ax = fig.add_subplot(n_samp,8,n*8+i)
    sc = "scaffold_" + str(i)
    chr8_rep = te8_rep[ (te8_rep.iloc[:,0].str.contains(sc))]
    plt.eventplot(chr8_rep[chr8_rep.iloc[:,4] == "pAge1"].iloc[:,1], linewidths=0.01, orientation='horizontal', colors="#2CA030FF")
    plt.eventplot(chr8_rep[chr8_rep.iloc[:,4] == "pAge2"].iloc[:,1], linewidths=0.01, orientation='horizontal', colors="#FFAA0EFF")
    plt.eventplot(chr8_rep[chr8_rep.iloc[:,4] == "pAa"].iloc[:,1], linewidths=0.01, orientation='horizontal', colors="#6F63BBFF")
    plt.eventplot(chr8_rep[chr8_rep.iloc[:,4] == "AthCEN159"].iloc[:,1], linewidths=0.01, orientation='horizontal', colors="#C7519CFF")
    ax.axes.xaxis.set_visible(False)
    ax.axes.yaxis.set_visible(False)
    ns_sc = ns[ns.iloc[:,0].str.contains(sc)]
    plt.eventplot(ns_sc.iloc[:,1], orientation='horizontal', colors="black", lineoffsets=0.2, linewidths=0.5)
    plt.eventplot(ns_sc.iloc[:,1], orientation='horizontal', colors="black", lineoffsets=1.8, linewidths=0.5)
    middle = (cen_pos.iloc[i-1,2] + cen_pos.iloc[i-1,1])/2
    plt.xlim([middle - 5000000, middle + 5000000])
    plt.ylim([0.1, 1.9])
n=50
te8_rep = pd.read_csv("mabs_assemblies/final_trash/AL26_v1_from_TE8_typen.bed", sep="\t", header=None)
cen_pos = pd.read_csv('mabs_assemblies/v2_assemblies/AL26_v1_centromeres.bed', sep="\t", header=None)
ns = pd.read_csv('mabs_assemblies/mapping/AL26_v1_Ns.bed', sep="\t", header=None)
for i in range(1,9):
    ax = fig.add_subplot(n_samp,8,n*8+i)
    sc = "scaffold_" + str(i)
    chr8_rep = te8_rep[ (te8_rep.iloc[:,0].str.contains(sc))]
    plt.eventplot(chr8_rep[chr8_rep.iloc[:,4] == "pAge1"].iloc[:,1], linewidths=0.01, orientation='horizontal', colors="#2CA030FF")
    plt.eventplot(chr8_rep[chr8_rep.iloc[:,4] == "pAge2"].iloc[:,1], linewidths=0.01, orientation='horizontal', colors="#FFAA0EFF")
    plt.eventplot(chr8_rep[chr8_rep.iloc[:,4] == "pAa"].iloc[:,1], linewidths=0.01, orientation='horizontal', colors="#6F63BBFF")
    plt.eventplot(chr8_rep[chr8_rep.iloc[:,4] == "AthCEN159"].iloc[:,1], linewidths=0.01, orientation='horizontal', colors="#C7519CFF")
    ax.axes.xaxis.set_visible(False)
    ax.axes.yaxis.set_visible(False)
    ns_sc = ns[ns.iloc[:,0].str.contains(sc)]
    plt.eventplot(ns_sc.iloc[:,1], orientation='horizontal', colors="black", lineoffsets=0.2, linewidths=0.5)
    plt.eventplot(ns_sc.iloc[:,1], orientation='horizontal', colors="black", lineoffsets=1.8, linewidths=0.5)
    middle = (cen_pos.iloc[i-1,2] + cen_pos.iloc[i-1,1])/2
    plt.xlim([middle - 5000000, middle + 5000000])
    plt.ylim([0.1, 1.9])
n=49
te8_rep = pd.read_csv("mabs_assemblies/final_trash/AL86_v1_from_TE8_typen.bed", sep="\t", header=None)
cen_pos = pd.read_csv('mabs_assemblies/v2_assemblies/AL86_v1_centromeres.bed', sep="\t", header=None)
ns = pd.read_csv('mabs_assemblies/mapping/AL86_v1_Ns.bed', sep="\t", header=None)
for i in range(1,9):
    ax = fig.add_subplot(n_samp,8,n*8+i)
    sc = "scaffold_" + str(i)
    chr8_rep = te8_rep[ (te8_rep.iloc[:,0].str.contains(sc))]
    plt.eventplot(chr8_rep[chr8_rep.iloc[:,4] == "pAge1"].iloc[:,1], linewidths=0.01, orientation='horizontal', colors="#2CA030FF")
    plt.eventplot(chr8_rep[chr8_rep.iloc[:,4] == "pAge2"].iloc[:,1], linewidths=0.01, orientation='horizontal', colors="#FFAA0EFF")
    plt.eventplot(chr8_rep[chr8_rep.iloc[:,4] == "pAa"].iloc[:,1], linewidths=0.01, orientation='horizontal', colors="#6F63BBFF")
    plt.eventplot(chr8_rep[chr8_rep.iloc[:,4] == "AthCEN159"].iloc[:,1], linewidths=0.01, orientation='horizontal', colors="#C7519CFF")
    ax.axes.xaxis.set_visible(False)
    ax.axes.yaxis.set_visible(False)
    ns_sc = ns[ns.iloc[:,0].str.contains(sc)]
    plt.eventplot(ns_sc.iloc[:,1], orientation='horizontal', colors="black", lineoffsets=0.2, linewidths=0.5)
    plt.eventplot(ns_sc.iloc[:,1], orientation='horizontal', colors="black", lineoffsets=1.8, linewidths=0.5)
    middle = (cen_pos.iloc[i-1,2] + cen_pos.iloc[i-1,1])/2
    plt.xlim([middle - 5000000, middle + 5000000])
    plt.ylim([0.1, 1.9])
n=48
te8_rep = pd.read_csv("mabs_assemblies/final_trash/AL85_v1_from_TE8_typen.bed", sep="\t", header=None)
cen_pos = pd.read_csv('mabs_assemblies/v2_assemblies/AL85_v1_centromeres.bed', sep="\t", header=None)
ns = pd.read_csv('mabs_assemblies/mapping/AL85_v1_Ns.bed', sep="\t", header=None)
for i in range(1,9):
    ax = fig.add_subplot(n_samp,8,n*8+i)
    sc = "scaffold_" + str(i)
    chr8_rep = te8_rep[ (te8_rep.iloc[:,0].str.contains(sc))]
    plt.eventplot(chr8_rep[chr8_rep.iloc[:,4] == "pAge1"].iloc[:,1], linewidths=0.01, orientation='horizontal', colors="#2CA030FF")
    plt.eventplot(chr8_rep[chr8_rep.iloc[:,4] == "pAge2"].iloc[:,1], linewidths=0.01, orientation='horizontal', colors="#FFAA0EFF")
    plt.eventplot(chr8_rep[chr8_rep.iloc[:,4] == "pAa"].iloc[:,1], linewidths=0.01, orientation='horizontal', colors="#6F63BBFF")
    plt.eventplot(chr8_rep[chr8_rep.iloc[:,4] == "AthCEN159"].iloc[:,1], linewidths=0.01, orientation='horizontal', colors="#C7519CFF")
    ax.axes.xaxis.set_visible(False)
    ax.axes.yaxis.set_visible(False)
    ns_sc = ns[ns.iloc[:,0].str.contains(sc)]
    plt.eventplot(ns_sc.iloc[:,1], orientation='horizontal', colors="black", lineoffsets=0.2, linewidths=0.5)
    plt.eventplot(ns_sc.iloc[:,1], orientation='horizontal', colors="black", lineoffsets=1.8, linewidths=0.5)
    middle = (cen_pos.iloc[i-1,2] + cen_pos.iloc[i-1,1])/2
    plt.xlim([middle - 5000000, middle + 5000000])
    plt.ylim([0.1, 1.9])
n=46
te8_rep = pd.read_csv("mabs_assemblies/final_trash/Plech_pca_typen_reclass.bed", sep="\t", header=None)
cen_pos = pd.read_csv('mabs_assemblies/v2_assemblies/Plech_v1_centromeres.bed', sep="\t", header=None)
ns = pd.read_csv('mabs_assemblies/mapping/Plech_v1_Ns.bed', sep="\t", header=None)
for i in range(1,9):
    ax = fig.add_subplot(n_samp,8,n*8+i)
    sc = "scaffold_" + str(i)
    chr8_rep = te8_rep[ (te8_rep.iloc[:,0].str.contains(sc))]
    plt.eventplot(chr8_rep[chr8_rep.iloc[:,3] == "pAge1"].iloc[:,1], linewidths=0.01, orientation='horizontal', colors="#2CA030FF")
    plt.eventplot(chr8_rep[chr8_rep.iloc[:,3] == "pAge2"].iloc[:,1], linewidths=0.01, orientation='horizontal', colors="#FFAA0EFF")
    plt.eventplot(chr8_rep[chr8_rep.iloc[:,3] == "pAa"].iloc[:,1], linewidths=0.01, orientation='horizontal', colors="#6F63BBFF")
    plt.eventplot(chr8_rep[chr8_rep.iloc[:,3] == "AthCEN159"].iloc[:,1], linewidths=0.01, orientation='horizontal', colors="#C7519CFF")
    ax.axes.xaxis.set_visible(False)
    ax.axes.yaxis.set_visible(False)
    ns_sc = ns[ns.iloc[:,0].str.contains(sc)]
    plt.eventplot(ns_sc.iloc[:,1], orientation='horizontal', colors="black", lineoffsets=0.2, linewidths=0.5)
    plt.eventplot(ns_sc.iloc[:,1], orientation='horizontal', colors="black", lineoffsets=1.8, linewidths=0.5)
    middle = (cen_pos.iloc[i-1,2] + cen_pos.iloc[i-1,1])/2
    plt.xlim([middle - 5000000, middle + 5000000])
    plt.ylim([0.1, 1.9])

n=45
te8_rep = pd.read_csv("mabs_assemblies/final_trash/S06_pca_typen_reclass.bed", sep="\t", header=None)
cen_pos = pd.read_csv('mabs_assemblies/v2_assemblies/S06_v1_centromeres.bed', sep="\t", header=None)
ns = pd.read_csv('mabs_assemblies/mapping/S06_v1_Ns.bed', sep="\t", header=None)
for i in range(1,9):
    ax = fig.add_subplot(n_samp,8,n*8+i)
    sc = "scaffold_" + str(i)
    chr8_rep = te8_rep[ (te8_rep.iloc[:,0].str.contains(sc))]
    plt.eventplot(chr8_rep[chr8_rep.iloc[:,3] == "pAge1"].iloc[:,1], linewidths=0.01, orientation='horizontal', colors="#2CA030FF")
    plt.eventplot(chr8_rep[chr8_rep.iloc[:,3] == "pAge2"].iloc[:,1], linewidths=0.01, orientation='horizontal', colors="#FFAA0EFF")
    plt.eventplot(chr8_rep[chr8_rep.iloc[:,3] == "pAa"].iloc[:,1], linewidths=0.01, orientation='horizontal', colors="#6F63BBFF")
    plt.eventplot(chr8_rep[chr8_rep.iloc[:,3] == "AthCEN159"].iloc[:,1], linewidths=0.01, orientation='horizontal', colors="#C7519CFF")
    ax.axes.xaxis.set_visible(False)
    ax.axes.yaxis.set_visible(False)
    ns_sc = ns[ns.iloc[:,0].str.contains(sc)]
    plt.eventplot(ns_sc.iloc[:,1], orientation='horizontal', colors="black", lineoffsets=0.2, linewidths=0.5)
    plt.eventplot(ns_sc.iloc[:,1], orientation='horizontal', colors="black", lineoffsets=1.8, linewidths=0.5)
    middle = (cen_pos.iloc[i-1,2] + cen_pos.iloc[i-1,1])/2
    plt.xlim([middle - 5000000, middle + 5000000])
    plt.ylim([0.1, 1.9])
n=44
te8_rep = pd.read_csv("mabs_assemblies/final_trash/S18_pca_typen_reclass.bed", sep="\t", header=None)
cen_pos = pd.read_csv('mabs_assemblies/v2_assemblies/S18_v1_merge_centromeres.bed', sep="\t", header=None)
ns = pd.read_csv('mabs_assemblies/mapping/S18_v1_Ns.bed', sep="\t", header=None)
for i in range(1,9):
    ax = fig.add_subplot(n_samp,8,n*8+i)
    sc = "scaffold_" + str(i)
    chr8_rep = te8_rep[ (te8_rep.iloc[:,0].str.contains(sc))]
    plt.eventplot(chr8_rep[chr8_rep.iloc[:,3] == "pAge1"].iloc[:,1], linewidths=0.01, orientation='horizontal', colors="#2CA030FF")
    plt.eventplot(chr8_rep[chr8_rep.iloc[:,3] == "pAge2"].iloc[:,1], linewidths=0.01, orientation='horizontal', colors="#FFAA0EFF")
    plt.eventplot(chr8_rep[chr8_rep.iloc[:,3] == "pAa"].iloc[:,1], linewidths=0.01, orientation='horizontal', colors="#6F63BBFF")
    plt.eventplot(chr8_rep[chr8_rep.iloc[:,3] == "AthCEN159"].iloc[:,1], linewidths=0.01, orientation='horizontal', colors="#C7519CFF")
    ax.axes.xaxis.set_visible(False)
    ax.axes.yaxis.set_visible(False)
    ns_sc = ns[ns.iloc[:,0].str.contains(sc)]
    plt.eventplot(ns_sc.iloc[:,1], orientation='horizontal', colors="black", lineoffsets=0.2, linewidths=0.5)
    plt.eventplot(ns_sc.iloc[:,1], orientation='horizontal', colors="black", lineoffsets=1.8, linewidths=0.5)
    middle = (cen_pos.iloc[i-1,2] + cen_pos.iloc[i-1,1])/2
    plt.xlim([middle - 5000000, middle + 5000000])
    plt.ylim([0.1, 1.9])
n=42
te8_rep = pd.read_csv("mabs_assemblies/final_trash/N7_pca_typen_reclass.bed", sep="\t", header=None)
cen_pos = pd.read_csv('mabs_assemblies/v2_assemblies/N7_v1_centromeres.bed', sep="\t", header=None)
ns = pd.read_csv('mabs_assemblies/mapping/N7_v1_Ns.bed', sep="\t", header=None)
for i in range(1,9):
    ax = fig.add_subplot(n_samp,8,n*8+i)
    sc = "scaffold_" + str(i)
    chr8_rep = te8_rep[ (te8_rep.iloc[:,0].str.contains(sc))]
    plt.eventplot(chr8_rep[chr8_rep.iloc[:,3] == "pAge1"].iloc[:,1], linewidths=0.01, orientation='horizontal', colors="#2CA030FF")
    plt.eventplot(chr8_rep[chr8_rep.iloc[:,3] == "pAge2"].iloc[:,1], linewidths=0.01, orientation='horizontal', colors="#FFAA0EFF")
    plt.eventplot(chr8_rep[chr8_rep.iloc[:,3] == "pAa"].iloc[:,1], linewidths=0.01, orientation='horizontal', colors="#6F63BBFF")
    plt.eventplot(chr8_rep[chr8_rep.iloc[:,3] == "AthCEN159"].iloc[:,1], linewidths=0.01, orientation='horizontal', colors="#C7519CFF")
    ax.axes.xaxis.set_visible(False)
    ax.axes.yaxis.set_visible(False)
    ns_sc = ns[ns.iloc[:,0].str.contains(sc)]
    plt.eventplot(ns_sc.iloc[:,1], orientation='horizontal', colors="black", lineoffsets=0.2, linewidths=0.5)
    plt.eventplot(ns_sc.iloc[:,1], orientation='horizontal', colors="black", lineoffsets=1.8, linewidths=0.5)
    middle = (cen_pos.iloc[i-1,2] + cen_pos.iloc[i-1,1])/2
    plt.xlim([middle - 5000000, middle + 5000000])
    plt.ylim([0.1, 1.9])
n=43
te8_rep = pd.read_csv("mabs_assemblies/final_trash/I04_pca_typen_reclass.bed", sep="\t", header=None)
cen_pos = pd.read_csv('mabs_assemblies/v2_assemblies/I04_v1_centromeres.bed', sep="\t", header=None)
ns = pd.read_csv('mabs_assemblies/mapping/I04_v1_Ns.bed', sep="\t", header=None)
for i in range(1,9):
    ax = fig.add_subplot(n_samp,8,n*8+i)
    sc = "scaffold_" + str(i)
    chr8_rep = te8_rep[ (te8_rep.iloc[:,0].str.contains(sc))]
    plt.eventplot(chr8_rep[chr8_rep.iloc[:,3] == "pAge1"].iloc[:,1], linewidths=0.01, orientation='horizontal', colors="#2CA030FF")
    plt.eventplot(chr8_rep[chr8_rep.iloc[:,3] == "pAge2"].iloc[:,1], linewidths=0.01, orientation='horizontal', colors="#FFAA0EFF")
    plt.eventplot(chr8_rep[chr8_rep.iloc[:,3] == "pAa"].iloc[:,1], linewidths=0.01, orientation='horizontal', colors="#6F63BBFF")
    plt.eventplot(chr8_rep[chr8_rep.iloc[:,3] == "AthCEN159"].iloc[:,1], linewidths=0.01, orientation='horizontal', colors="#C7519CFF")
    ax.axes.xaxis.set_visible(False)
    ax.axes.yaxis.set_visible(False)
    ns_sc = ns[ns.iloc[:,0].str.contains(sc)]
    plt.eventplot(ns_sc.iloc[:,1], orientation='horizontal', colors="black", lineoffsets=0.2, linewidths=0.5)
    plt.eventplot(ns_sc.iloc[:,1], orientation='horizontal', colors="black", lineoffsets=1.8, linewidths=0.5)
    middle = (cen_pos.iloc[i-1,2] + cen_pos.iloc[i-1,1])/2
    plt.xlim([middle - 5000000, middle + 5000000])
    plt.ylim([0.1, 1.9])
n=38
te8_rep = pd.read_csv("mabs_assemblies/final_trash/BAM12.3_pca_typen.bed", sep="\t", header=None)
cen_pos = pd.read_csv('mabs_assemblies/v2_assemblies/expand_centromeres/BAM12_v2.bed', sep="\t", header=None)
ns = pd.read_csv('mabs_assemblies/mapping/BAM12_v2_Ns.bed', sep="\t", header=None)
for i in range(1,9):
    ax = fig.add_subplot(n_samp,8,n*8+i)
    sc = "scaffold_" + str(i)
    chr8_rep = te8_rep[ (te8_rep.iloc[:,0].str.contains(sc))]
    plt.eventplot(chr8_rep[chr8_rep.iloc[:,3] == "pAge1"].iloc[:,1], linewidths=0.01, orientation='horizontal', colors="#2CA030FF")
    plt.eventplot(chr8_rep[chr8_rep.iloc[:,3] == "pAa"].iloc[:,1], linewidths=0.01, orientation='horizontal', colors="#6F63BBFF")
    plt.eventplot(chr8_rep[chr8_rep.iloc[:,3] == "pAge2"].iloc[:,1], linewidths=0.01, orientation='horizontal', colors="#FFAA0EFF")
    ax.axes.xaxis.set_visible(False)
    ax.axes.yaxis.set_visible(False)
    ns_sc = ns[ns.iloc[:,0].str.contains(sc)]
    plt.eventplot(ns_sc.iloc[:,1], orientation='horizontal', colors="black", lineoffsets=0.2, linewidths=0.5)
    plt.eventplot(ns_sc.iloc[:,1], orientation='horizontal', colors="black", lineoffsets=1.8, linewidths=0.5)
    middle = (cen_pos.iloc[i-1,2] + cen_pos.iloc[i-1,1])/2
    plt.xlim([middle - 5000000, middle + 5000000])
    plt.ylim([0.1, 1.9])
n=41
te8_rep = pd.read_csv("mabs_assemblies/final_trash/NT8_pca_typen.bed", sep="\t", header=None)
cen_pos = pd.read_csv('mabs_assemblies/v2_assemblies/expand_centromeres/NT8_v2.bed', sep="\t", header=None)
ns = pd.read_csv('mabs_assemblies/mapping/NT8_v2_Ns.bed', sep="\t", header=None)
for i in range(1,9):
    ax = fig.add_subplot(n_samp,8,n*8+i)
    sc = "scaffold_" + str(i)
    chr8_rep = te8_rep[ (te8_rep.iloc[:,0].str.contains(sc))]
    plt.eventplot(chr8_rep[chr8_rep.iloc[:,3] == "pAge1"].iloc[:,1], linewidths=0.01, orientation='horizontal', colors="#2CA030FF")
    plt.eventplot(chr8_rep[chr8_rep.iloc[:,3] == "pAa"].iloc[:,1], linewidths=0.01, orientation='horizontal', colors="#6F63BBFF")
    plt.eventplot(chr8_rep[chr8_rep.iloc[:,3] == "pAge2"].iloc[:,1], linewidths=0.01, orientation='horizontal', colors="#FFAA0EFF")
    ax.axes.xaxis.set_visible(False)
    ax.axes.yaxis.set_visible(False)
    ns_sc = ns[ns.iloc[:,0].str.contains(sc)]
    plt.eventplot(ns_sc.iloc[:,1], orientation='horizontal', colors="black", lineoffsets=0.2, linewidths=0.5)
    plt.eventplot(ns_sc.iloc[:,1], orientation='horizontal', colors="black", lineoffsets=1.8, linewidths=0.5)
    middle = (cen_pos.iloc[i-1,2] + cen_pos.iloc[i-1,1])/2
    plt.xlim([middle - 5000000, middle + 5000000])
    plt.ylim([0.1, 1.9])
    
n=36
te8_rep = pd.read_csv("mabs_assemblies/final_trash/NT9_pca_typen.bed", sep="\t", header=None)
cen_pos = pd.read_csv('mabs_assemblies/v2_assemblies/expand_centromeres/NT9_v2.bed', sep="\t", header=None)
ns = pd.read_csv('mabs_assemblies/mapping/NT9_v2_Ns.bed', sep="\t", header=None)
for i in range(1,9):
    ax = fig.add_subplot(n_samp,8,n*8+i)
    sc = "scaffold_" + str(i)
    chr8_rep = te8_rep[ (te8_rep.iloc[:,0].str.contains(sc))]
    plt.eventplot(chr8_rep[chr8_rep.iloc[:,3] == "pAge1"].iloc[:,1], linewidths=0.01, orientation='horizontal', colors="#2CA030FF")
    plt.eventplot(chr8_rep[chr8_rep.iloc[:,3] == "pAa"].iloc[:,1], linewidths=0.01, orientation='horizontal', colors="#6F63BBFF")
    plt.eventplot(chr8_rep[chr8_rep.iloc[:,3] == "pAge2"].iloc[:,1], linewidths=0.01, orientation='horizontal', colors="#FFAA0EFF")
    ax.axes.xaxis.set_visible(False)
    ax.axes.yaxis.set_visible(False)
    ns_sc = ns[ns.iloc[:,0].str.contains(sc)]
    plt.eventplot(ns_sc.iloc[:,1], orientation='horizontal', colors="black", lineoffsets=0.2, linewidths=0.5)
    plt.eventplot(ns_sc.iloc[:,1], orientation='horizontal', colors="black", lineoffsets=1.8, linewidths=0.5)
    middle = (cen_pos.iloc[i-1,2] + cen_pos.iloc[i-1,1])/2
    plt.xlim([middle - 5000000, middle + 5000000])
    plt.ylim([0.1, 1.9])
    
n=37
te8_rep = pd.read_csv("mabs_assemblies/final_trash/al1_pca_typen.bed", sep="\t", header=None)
cen_pos = pd.read_csv('mabs_assemblies/v2_assemblies/expand_centromeres/al1_v2.1.bed', sep="\t", header=None)
ns = pd.read_csv('mabs_assemblies/mapping/al1_v2.1_Ns.bed', sep="\t", header=None)
for i in range(1,9):
    ax = fig.add_subplot(n_samp,8,n*8+i)
    sc = "scaffold_" + str(i)
    chr8_rep = te8_rep[ (te8_rep.iloc[:,0].str.contains(sc))]
    plt.eventplot(chr8_rep[chr8_rep.iloc[:,3] == "pAge1"].iloc[:,1], linewidths=0.01, orientation='horizontal', colors="#2CA030FF")
    plt.eventplot(chr8_rep[chr8_rep.iloc[:,3] == "pAa"].iloc[:,1], linewidths=0.01, orientation='horizontal', colors="#6F63BBFF")
    plt.eventplot(chr8_rep[chr8_rep.iloc[:,3] == "pAge2"].iloc[:,1], linewidths=0.01, orientation='horizontal', colors="#FFAA0EFF")
    ax.axes.xaxis.set_visible(False)
    ax.axes.yaxis.set_visible(False)
    ns_sc = ns[ns.iloc[:,0].str.contains(sc)]
    plt.eventplot(ns_sc.iloc[:,1], orientation='horizontal', colors="black", lineoffsets=0.2, linewidths=0.5)
    plt.eventplot(ns_sc.iloc[:,1], orientation='horizontal', colors="black", lineoffsets=1.8, linewidths=0.5)
    middle = (cen_pos.iloc[i-1,2] + cen_pos.iloc[i-1,1])/2
    plt.xlim([middle - 5000000, middle + 5000000])
    plt.ylim([0.1, 1.9])
n = 40

te8_rep = pd.read_csv("mabs_assemblies/final_trash/NT12_pca_typen.bed", sep="\t", header=None)
cen_pos = pd.read_csv('mabs_assemblies/v2_assemblies/expand_centromeres/NT12_v2.bed', sep="\t", header=None)
ns = pd.read_csv('mabs_assemblies/mapping/NT12_v2_Ns.bed', sep="\t", header=None)
for i in range(1,9):
    ax = fig.add_subplot(n_samp,8,n*8+i)
    sc = "scaffold_" + str(i)
    chr8_rep = te8_rep[ (te8_rep.iloc[:,0].str.contains(sc))]
    plt.eventplot(chr8_rep[chr8_rep.iloc[:,3] == "pAge1"].iloc[:,1], linewidths=0.01, orientation='horizontal', colors="#2CA030FF")
    plt.eventplot(chr8_rep[chr8_rep.iloc[:,3] == "pAa"].iloc[:,1], linewidths=0.01, orientation='horizontal', colors="#6F63BBFF")
    plt.eventplot(chr8_rep[chr8_rep.iloc[:,3] == "pAge2"].iloc[:,1], linewidths=0.2, orientation='horizontal', colors="#FFAA0EFF")
    ax.axes.xaxis.set_visible(False)
    ax.axes.yaxis.set_visible(False)
    ns_sc = ns[ns.iloc[:,0].str.contains(sc)]
    plt.eventplot(ns_sc.iloc[:,1], orientation='horizontal', colors="black", lineoffsets=0.2, linewidths=0.5)
    plt.eventplot(ns_sc.iloc[:,1], orientation='horizontal', colors="black", lineoffsets=1.8, linewidths=0.5)
    middle = (cen_pos.iloc[i-1,2] + cen_pos.iloc[i-1,1])/2
    plt.xlim([middle - 5000000, middle + 5000000])
    plt.ylim([0.1, 1.9])    
n=39
te8_rep = pd.read_csv("mabs_assemblies/final_trash/BAM12.1_v1_from_TE8_typen.bed", sep="\t", header=None)
cen_pos = pd.read_csv('mabs_assemblies/v2_assemblies/expand_centromeres/BAM12.1_v1.bed', sep="\t", header=None)
ns = pd.read_csv('mabs_assemblies/mapping/BAM12.1_v1_Ns.bed', sep="\t", header=None)
for i in range(1,9):
    ax = fig.add_subplot(n_samp,8,n*8+i)
    sc = "scaffold_" + str(i)
    chr8_rep = te8_rep[ (te8_rep.iloc[:,0].str.contains(sc))]
    plt.eventplot(chr8_rep[chr8_rep.iloc[:,4] == "pAge2"].iloc[:,1], linewidths=0.01, orientation='horizontal', colors="#FFAA0EFF")
    plt.eventplot(chr8_rep[chr8_rep.iloc[:,4] == "pAge1"].iloc[:,1], linewidths=0.01, orientation='horizontal', colors="#2CA030FF")
    plt.eventplot(chr8_rep[chr8_rep.iloc[:,4] == "pAa"].iloc[:,1], linewidths=0.01, orientation='horizontal', colors="#6F63BBFF")
    ax.axes.xaxis.set_visible(False)
    ax.axes.yaxis.set_visible(False)
    ns_sc = ns[ns.iloc[:,0].str.contains(sc)]
    plt.eventplot(ns_sc.iloc[:,1], orientation='horizontal', colors="black", lineoffsets=0.2, linewidths=0.5)
    plt.eventplot(ns_sc.iloc[:,1], orientation='horizontal', colors="black", lineoffsets=1.8, linewidths=0.5)
    middle = (cen_pos.iloc[i-1,2] + cen_pos.iloc[i-1,1])/2
    plt.xlim([middle - 5000000, middle + 5000000])
    plt.ylim([0.1, 1.9])

n=30
te8_rep = pd.read_csv("mabs_assemblies/final_trash/Kar_pca_typen.bed", sep="\t", header=None)
cen_pos = pd.read_csv('mabs_assemblies/v2_assemblies/Kar_v1_merge_centromeres.bed', sep="\t", header=None)
ns = pd.read_csv('mabs_assemblies/mapping/Kar_v1_Ns.bed', sep="\t", header=None)
for i in range(1,9):
    ax = fig.add_subplot(n_samp,8,n*8+i)
    sc = "scaffold_" + str(i)
    chr8_rep = te8_rep[ (te8_rep.iloc[:,0].str.contains(sc))]
    plt.eventplot(chr8_rep[chr8_rep.iloc[:,3] == "pAge1"].iloc[:,1], linewidths=0.01, orientation='horizontal', colors="#2CA030FF")
    plt.eventplot(chr8_rep[chr8_rep.iloc[:,3] == "pAge2"].iloc[:,1], linewidths=0.01, orientation='horizontal', colors="#FFAA0EFF")
    plt.eventplot(chr8_rep[chr8_rep.iloc[:,3] == "pAa"].iloc[:,1], linewidths=0.01, orientation='horizontal', colors="#6F63BBFF")
    ax.axes.xaxis.set_visible(False)
    ax.axes.yaxis.set_visible(False)
    ns_sc = ns[ns.iloc[:,0].str.contains(sc)]
    plt.eventplot(ns_sc.iloc[:,1], orientation='horizontal', colors="black", lineoffsets=0.2, linewidths=0.5)
    plt.eventplot(ns_sc.iloc[:,1], orientation='horizontal', colors="black", lineoffsets=1.8, linewidths=0.5)
    middle = (cen_pos.iloc[i-1,2] + cen_pos.iloc[i-1,1])/2
    plt.xlim([middle - 5000000, middle + 5000000])
    plt.ylim([0.1, 1.9])
n=25
te8_rep = pd.read_csv("mabs_assemblies/final_trash/MN47_v4_from_TE8_typen.bed", sep="\t", header=None)
cen_pos = pd.read_csv('mabs_assemblies/v2_assemblies/expand_centromeres/MN47_v4.bed', sep="\t", header=None)
ns = pd.read_csv('mabs_assemblies/mapping/MN47_v4_Ns.bed', sep="\t", header=None)
for i in range(1,9):
    ax = fig.add_subplot(n_samp,8,n*8+i)
    sc = "scaffold_" + str(i)
    chr8_rep = te8_rep[ (te8_rep.iloc[:,0].str.contains(sc))]
    plt.eventplot(chr8_rep[chr8_rep.iloc[:,4] == "pAge2"].iloc[:,1], alpha=0.5, orientation='horizontal', colors="#FFAA0EFF")
    plt.eventplot(chr8_rep[chr8_rep.iloc[:,4] == "pAa"].iloc[:,1], alpha=0.5, orientation='horizontal', colors="#6F63BBFF")
    plt.eventplot(chr8_rep[chr8_rep.iloc[:,4] == "pAge1"].iloc[:,1], alpha=0.5, orientation='horizontal', colors="#2CA030FF")
    ax.axes.xaxis.set_visible(False)
    ax.axes.yaxis.set_visible(False)
    ns_sc = ns[ns.iloc[:,0].str.contains(sc)]
    plt.eventplot(ns_sc.iloc[:,1], orientation='horizontal', colors="black", lineoffsets=0.2, linewidths=0.5)
    plt.eventplot(ns_sc.iloc[:,1], orientation='horizontal', colors="black", lineoffsets=1.8, linewidths=0.5)
    middle = (cen_pos.iloc[i-1,2] + cen_pos.iloc[i-1,1])/2
    plt.xlim([middle - 5000000, middle + 5000000])
    plt.ylim([0.1, 1.9])
    
n=26
te8_rep = pd.read_csv("mabs_assemblies/final_trash/NT1_v2_from_TE8_typen.bed", sep="\t", header=None)
cen_pos = pd.read_csv('mabs_assemblies/v2_assemblies/expand_centromeres/NT1_v2.bed', sep="\t", header=None)
ns = pd.read_csv('mabs_assemblies/mapping/NT1_v2_Ns.bed', sep="\t", header=None)
for i in range(1,9):
    ax = fig.add_subplot(n_samp,8,n*8+i)
    sc = "scaffold_" + str(i)
    chr8_rep = te8_rep[ (te8_rep.iloc[:,0].str.contains(sc))]
    plt.eventplot(chr8_rep[chr8_rep.iloc[:,4] == "pAge2"].iloc[:,1], alpha=0.5, orientation='horizontal', colors="#FFAA0EFF")
    plt.eventplot(chr8_rep[chr8_rep.iloc[:,4] == "pAge1"].iloc[:,1], alpha=0.5, orientation='horizontal', colors="#2CA030FF")
    plt.eventplot(chr8_rep[chr8_rep.iloc[:,4] == "pAa"].iloc[:,1], alpha=0.5, orientation='horizontal', colors="#6F63BBFF")
    ax.axes.xaxis.set_visible(False)
    ax.axes.yaxis.set_visible(False)
    ns_sc = ns[ns.iloc[:,0].str.contains(sc)]
    plt.eventplot(ns_sc.iloc[:,1], orientation='horizontal', colors="black", lineoffsets=0.2, linewidths=0.5)
    plt.eventplot(ns_sc.iloc[:,1], orientation='horizontal', colors="black", lineoffsets=1.8, linewidths=0.5)
    middle = (cen_pos.iloc[i-1,2] + cen_pos.iloc[i-1,1])/2
    plt.xlim([middle - 5000000, middle + 5000000])
    plt.ylim([0.1, 1.9])
n=54
te8_rep = pd.read_csv("mabs_assemblies/final_trash/AL56_v1_from_TE8_typen.bed", sep="\t", header=None)
cen_pos = pd.read_csv('mabs_assemblies/v2_assemblies/AL56_v1_centromeres.bed', sep="\t", header=None)
ns = pd.read_csv('mabs_assemblies/mapping/AL56_v1_Ns.bed', sep="\t", header=None)
for i in range(1,9):
    ax = fig.add_subplot(n_samp,8,n*8+i)
    sc = "scaffold_" + str(i)
    chr8_rep = te8_rep[ (te8_rep.iloc[:,0].str.contains(sc))]
    plt.eventplot(chr8_rep[chr8_rep.iloc[:,4] == "pAge1"].iloc[:,1], linewidths=0.01, orientation='horizontal', colors="#2CA030FF")
    plt.eventplot(chr8_rep[chr8_rep.iloc[:,4] == "pAge2"].iloc[:,1], linewidths=0.01, orientation='horizontal', colors="#FFAA0EFF")
    plt.eventplot(chr8_rep[chr8_rep.iloc[:,4] == "pAa"].iloc[:,1], linewidths=0.01, orientation='horizontal', colors="#6F63BBFF")
    plt.eventplot(chr8_rep[chr8_rep.iloc[:,4] == "AthCEN159"].iloc[:,1], linewidths=0.01, orientation='horizontal', colors="#C7519CFF")
    ax.axes.xaxis.set_visible(False)
    ax.axes.yaxis.set_visible(False)
    ns_sc = ns[ns.iloc[:,0].str.contains(sc)]
    plt.eventplot(ns_sc.iloc[:,1], orientation='horizontal', colors="black", lineoffsets=0.2, linewidths=0.5)
    plt.eventplot(ns_sc.iloc[:,1], orientation='horizontal', colors="black", lineoffsets=1.8, linewidths=0.5)
    middle = (cen_pos.iloc[i-1,2] + cen_pos.iloc[i-1,1])/2
    plt.xlim([middle - 5000000, middle + 5000000])
    plt.ylim([0.1, 1.9])
plt.savefig("Velky_plot_thal.png")
