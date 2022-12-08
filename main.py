#!/usr/bin/env python
# coding: utf-8

import PDBdata_Scraper as Scraper
import pandas as pd

#df1 = Scraper.df_pro
#df2 = Scraper.df_euk
#df3 = Scraper.df_virus
pdbid_list = Scraper.idsdata_list

# Run demo to get the raw dict for pdbTM database
#pdbtm_dict = Scraper.Merge_pdb_GC(pdbid_list, df1, df2, df3)
#print(pdbtm_dict)
#test_list = pdbid_list[0:100]
#print(pdbid_list)


def get_pdb_tax_csv(pdblist):
    i=0
    pdb_taxid_list = []
    df = pd.DataFrame()
    for pdb in pdblist:
        i += 1
        pdb_taxid = Scraper.get_pdb_taxlist(pdb)
       
        pdb_taxid_list.append(pdb_taxid)
         
        print(i)
        print(pdb, pdb_taxid)
        
    df['PDB_ID'] = pdblist 
 
    df['TAX_ID'] = pdb_taxid_list
    df.to_csv('PDB_TaxID.csv')
    
    return 

get_pdb_tax_csv(pdbid_list)

