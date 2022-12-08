#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  8 12:18:05 2022

@author: weiqiyao
"""
import PDBdata_Scraper as Scraper
import pandas as pd

df1 = Scraper.df_pro
df2 = Scraper.df_euk
df3 = Scraper.df_virus

df = pd.read_csv('PDB_TaxID_final.csv',index_col= False)
df = df.drop(columns='Unnamed: 0')

f = open("PDB_GC_dict.txt","w")

PDB_GC_dict = Scraper.Merge_pdb_GC(df, df1, df2, df3)

# write file
f.write( str(PDB_GC_dict) )

# close file
f.close()


