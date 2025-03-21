import pandas as pd
import numpy as np

import warnings
warnings.filterwarnings("ignore")

rna_df=pd.read_csv('table.4C-1.DEGs-newTPM_ER-HER2-.csv', sep='\t')
protein_df=pd.read_csv('table.5A.DEPs.ER-HER2-.txt', sep='\t')
p_protein_df=pd.read_csv('table.5B.DEPPs.ER-HER2-.txt', sep='\t')
me_df=pd.read_csv('table.3C.DMGs.ER-HER2-.txt', sep='\t')

gene_name_list=list(rna_df["GeneName"])

mtx=[]
for gene in gene_name_list:
    
    a=float(rna_df[rna_df["GeneName"]==gene]["log2FoldChange"])
    b=float(rna_df[rna_df["GeneName"]==gene]["prank"])

    try:
        c=float(me_df[me_df["GeneName"]==gene]["log2foldChange"])
        d=float(me_df[me_df["GeneName"]==gene]["pValues"])
    except:
        c='NA'
        d='NA'
    
    try:
        e=float(protein_df[protein_df["GeneName"]==gene]["log2FoldChange"])
        f=float(protein_df[protein_df["GeneName"]==gene]["Pvalue_Ranksum"])
    except:
        e='NA'
        f='NA'
        
    try:
        h=float(p_protein_df[p_protein_df["GeneName"]==gene]["log2FoldChange"])
        i=float(p_protein_df[p_protein_df["GeneName"]==gene]["Pvalue_Ranksum"])

    except:
        h='NA'
        i='NA'

    mtx.append([gene,d,c,b,a,f,e,i,h])

df=pd.DataFrame(mtx,columns=["gene","me_pval","me_log2fc","rna_pval","rna_log2fc","protein_pval","protein_log2fc","p_protein_pval","p_protein_log2fc"])


df.to_csv("data_ER-HER2-_1.csv", index=False, encoding="utf-8")

