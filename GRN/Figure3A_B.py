#!/usr/bin/env python
# coding: utf-8
# In[1]:
import pandas as pd
import sys,os
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.pyplot as plt
get_ipython().run_line_magic('config', "InlineBackend.figure_format = 'retina'")
# In[3]:
#
import scanpy as sc
sc.settings.set_figure_params(dpi=80, figsize=(4, 4),facecolor='white')
# In[4]:
sns.set_style(style='ticks')
# # Fig. 3A
# In[8]:
peak_df=pd.read_csv('../scATAC/SnapATAC2/final_filtered/peak_annotation.tsv',sep='\t', index_col=0, low_memory=False)
peak_df.head()
# In[9]:
peak_anno=list(peak_df['Detailed Annotation'])
for i, anno in enumerate(peak_anno):
    if anno=='Intergenic':
        continue
    if anno=='CpG':
        continue
    if 'CpG' in anno:
        peak_anno[i]='CpG'
    if '|LINE|' in anno:
        peak_anno[i]='LINE'
    if '|SINE|' in anno:
        peak_anno[i]='SINE'
    if '|LTR|' in anno:
        peak_anno[i]='LTR'
    if 'intron' in anno:
        peak_anno[i]='intron'
    if 'exon' in anno:
        peak_anno[i]='exon'
    if 'TTS' in anno:
        peak_anno[i]='TTS'
    if 'promoter-TSS' in anno:
        peak_anno[i]='promoter-TSS'
    if "3' UTR" in anno:
        peak_anno[i]='3’UTR'
    if "5' UTR" in anno:
        peak_anno[i]='5’UTR'
    if '|DNA|' in anno:
        peak_anno[i]='DNA repeats'
    if '|RNA|' in anno or '|tRNA|' in anno or '|rRNA|' in anno or '|snRNA|' in anno or '|scRNA|' in anno or '|srpRNA|' in anno:
        peak_anno[i]='RNA repeats'
    if '|Simple_repeat|' in anno or '|Low_complexity|' in anno or '|Satellite|' in anno or '|RC|' in anno:
        peak_anno[i]='Other repeats'

# In[11]:
from collections import Counter
Counter(peak_anno)
# In[ ]:
import matplotlib.pyplot as plt
# given dictionary
count_dict = {
    '5’UTR': 1047,
    'promoter-TSS': 23343,
    'Intergenic': 69911,
    '3’UTR': 4561,
    'Other repeats': 3909,
    'DNA repeats': 10042,
    'RNA repeats': 323,
    'LINE': 40426,
    'CpG': 7139,
    'SINE': 22301,
    'LTR': 33671,
    'Non-coding':331,
    'TTS': 5852,
    'Exon': 9729,
    'Intron': 105451
}
# Get the dictionary keys and values ​​as the pie chart labels and sizes, respectively.
labels = count_dict.keys()
sizes = count_dict.values()

# Define color
colors = plt.cm.tab20.colors  # Use the colors in the colormap

# Draw a pie chart图
plt.figure(figsize=(8, 8))  # Set canvas size
plt.pie(sizes, labels=labels, colors=colors, autopct='%1.1f%%', startangle=180, textprops={'fontsize': 10})

# Display in equal proportions, ensuring the pie chart is circular
plt.axis('equal')

# # Show legend
plt.legend(loc='upper left', bbox_to_anchor=(1, 1), fontsize=8)

# set tittle
plt.title('Distribution of Elements', fontsize=14)
#plt.savefig('../Figures/Figure2/PieChart_class_cCREs.pdf')
plt.show()
# # Fig. 3B
# ## Left panel
# In[ ]:
#
'../Natural_Cohort/downloaded/SCREEN/GRCh38-cCREs.bed' # this file is downloaded from SCREEN database
get_ipython().system('bedtools intersect -a ./peaks.bed -b ../Natural_Cohort/downloaded/SCREEN/GRCh38-cCREs.bed -wa -wb -f 0.20 > overlap_rDHSs_cCREs.bed')
# In[21]:
overlap_df=pd.read_csv('../scATAC/SnapATAC2/final_filtered/CRE_analysis/DHS/overlap_rDHSs_cCREs.bed', sep='\t', header=None)
overlap_df.head()
# In[22]:
# report unique cCREs in each dataset
A=overlap_df[0].astype(str)+ '-' + overlap_df[1].astype(str)+ '-' + overlap_df[2].astype(str) # our data
B=overlap_df[3].astype(str)+ '-' + overlap_df[4].astype(str)+ '-' + overlap_df[5].astype(str) # SCREEN database
len(A), len(set(A)), len(B), len(set(B))
# In[ ]:
import matplotlib.pyplot as plt
from matplotlib_venn import venn2
#
plt.figure(figsize=(6,6))

# Define the size of two sets
set1 = 338036
set2 = 1063878
venn_set_size = min(set1, set2) # keep two cycles in same size
overlap = 253511
set_colors = ('blue', 'green')
v = venn2(subsets=(set1, set2, overlap), 
          set_labels=('cCREs', 'rDHSs'),
          set_colors=set_colors)
v.get_patch_by_id('11').set_color(set_colors[0])
v.get_patch_by_id('11').set_alpha(0.4) 
#
#plt.savefig('../../../Figures/Figure2/Venn_DHSs_CREs_overlap.pdf')

# 显示图形
plt.show()
# ## Right panel
# In[14]:
overlap_df=pd.read_csv('../scATAC/SnapATAC2/final_filtered/CRE_analysis/DHS/overlap_rDHSs_cCREs.bed', sep='\t', header=None)
overlap_peaks=set(overlap_df[0].astype(str)+ ':' + overlap_df[1].astype(str)+ '-' + overlap_df[2].astype(str))
# In[17]:
peak_df=pd.read_excel('../scATAC/SnapATAC2/final_filtered/FinalAnnotationPeak.xlsx')
peak_df=peak_df.set_index('Peaks')
#
peak_df2=pd.DataFrame(peak_df.index)
peak_df2['celltype_nums']=list(peak_df.sum(axis=1)) # add number of celltypes in peak calling
peak_df2['DHS']=['ovlp' if i in overlap_peaks else 'non-ovlp' for i in peak_df.index] # add if overlap with rDHSs
peak_df2['celltype_nums_modify']=[i if i < 10 else 10 for i in peak_df2['celltype_nums']]
# In[18]:
count_df=pd.DataFrame(peak_df2.value_counts(["DHS", "celltype_nums_modify"])).reset_index().sort_values('celltype_nums_modify')
#
count_df1=count_df[count_df['DHS']=='ovlp']
count_df2=count_df[count_df['DHS']=='non-ovlp']
# In[ ]:
# plot PieChart of ovlp and non-ovlp
import matplotlib.pyplot as plt

# Get the dictionary keys and values ​​as the pie chart labels and sizes, respectively.
labels = count_df1['celltype_nums_modify'].values
sizes = count_df1['count'].values

# Define color
colors = plt.cm.tab20.colors
plt.figure(figsize=(6, 6))
plt.pie(sizes, labels=labels, colors=colors, autopct='%1.1f%%', startangle=90, textprops={'fontsize': 10})
plt.axis('equal')
plt.legend(loc='upper left', bbox_to_anchor=(1, 1), fontsize=8)
plt.title('Ovlp cCREs', fontsize=14)
#plt.savefig('../Figures/Figure2/PieChart_OvlpCREs_specificity.pdf')
plt.show()


# In[ ]:

import matplotlib.pyplot as plt
labels = count_df2['celltype_nums_modify'].values
sizes = count_df2['count'].values
colors = plt.cm.tab20.colors 
plt.figure(figsize=(6, 6))
plt.pie(sizes, labels=labels, colors=colors, autopct='%1.1f%%', startangle=90, textprops={'fontsize': 10})
plt.axis('equal')
plt.legend(loc='upper left', bbox_to_anchor=(1, 1), fontsize=8)
plt.title('Non-ovlp cCREs', fontsize=14)
#plt.savefig('../Figures/Figure2/PieChart_non-OvlpCREs_specificity.pdf')
plt.show()

