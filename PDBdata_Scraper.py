#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd
#import os
#import numpy as np
from Bio import SeqIO
import requests
#from io import StringIO
from Bio import SeqIO
import time
#from Bio import Entrez


#
# read pdbtm database and convert chain id into id for information of organism
ids_list = []
for record in SeqIO.parse("pdbtm_all.seq.fasta", "fasta"):
    id_chain = record.id
    ids = id_chain[0:4]
    ids_list.append(ids)
    idsdata_list = list(set(ids_list))
    #print(ids)
    #print(record.seq)


# some of the pdb data is not available on RCSB. Need to filter these
def check_pdb_access(pdb_id):
    #time.sleep(0.01)
    pdbquery = 'https://data.rcsb.org/graphql?query={entry(entry_id:"%s"){exptl {method}}}' %pdb_id 
    r = requests.get(pdbquery)
    info_dic = r.json()
    data = info_dic['data']
    if None is data['entry']:
        return 'PDB data unavailable'
    else:
        return 'available'
    
'''   
# function to scrape organism information and taxonomy id from RCSB
def get_pdb_info(pdb_id):
    # check availability of pdb data
    if check_pdb_access(pdb_id) == 'PDB data unavailable':
        return 'PDB data unavailable'
    else: 
        info_list = []
        taxid_list = []
        
        #specify query entry to extract taxonomy id and organism name based on the entry id
        pdbquery = 'https://data.rcsb.org/graphql?query={entry(entry_id:"%s"){polymer_entities {entity_src_gen{pdbx_gene_src_ncbi_taxonomy_id, pdbx_gene_src_scientific_name}}}}' %pdb_id
        #make the request
        r = requests.get(pdbquery)
        info_dic = r.json()
        #parse our needed info from the integrated dict 
        data = info_dic['data']
        entry = data['entry']
        # some of the pdb data has been removed from the RCSB. filter these
        
        pdb_info = entry['polymer_entities']
        #n = 0 
        for i in pdb_info:
            #n += 1
            #print (i)
            #print(n)
            if i['entity_src_gen'] == None:
                #print('skip')
                continue
            else:
                info_list.append(i['entity_src_gen'][0])
        
        #print(infolist)
            for info in info_list:
                #print(info)
                taxid = info['pdbx_gene_src_ncbi_taxonomy_id']
                taxid_list.append(taxid)
                #print(taxid)
        #remove duplicates
        taxid_list = list(set(taxid_list))
     
        return taxid_list
'''

'''
# function to scrape organism information and taxonomy id from RCSB
# if function 1 do not find the taxid. might be another format of API
def get_pdb_info2(pdb_id):
    if check_pdb_access(pdb_id) == 'PDB data unavailable':
        return 'PDB data unavailable'
    else:
        
        info_list = []
        taxid_list = []
        
        #specify query entry to extract taxonomy id and organism name based on the entry id
        pdbquery = 'https://data.rcsb.org/graphql?query={entry(entry_id:"%s")'        '{polymer_entities {pdbx_entity_src_syn '        ' {ncbi_taxonomy_id, organism_scientific}}}}' %pdb_id
        #make the request
        r = requests.get(pdbquery)
        info_dic = r.json()
        #parse our needed info from the integrated dict 
        data = info_dic['data']
        entry = data['entry']
        pdb_info = entry['polymer_entities']
        #n = 0 
        for i in pdb_info:
            #n += 1
            #print (i)
            #print(n)
            if i['pdbx_entity_src_syn'] == None:
                #print('skip')
                continue
            else:
                info_list.append(i['pdbx_entity_src_syn'][0])
        
        #print(infolist)
            for info in info_list:
                #print(info)
                taxid = info['ncbi_taxonomy_id']
                taxid_list.append(taxid)
                #print(taxid)
        #remove duplicates
        taxid_list = list(set(taxid_list))
        
        return taxid_list
'''

# merge all of the possible API
def get_pdb_info_total(pdb_id):

    entry = ['entity_src_gen','entity_src_nat','pdbx_entity_src_syn']
    taxonomy = {'entity_src_gen':'pdbx_gene_src_ncbi_taxonomy_id',
           'entity_src_nat':'pdbx_ncbi_taxonomy_id',
           'pdbx_entity_src_syn':'ncbi_taxonomy_id'}
    
    # check availability of pdb data
    if check_pdb_access(pdb_id) == 'PDB data unavailable':
        return 'PDB data unavailable'
    else: 
        info_list = []
        taxid_list = []
        for s2 in entry:
            #print(s2)
            s3 = taxonomy[s2]
            #print(s3)
            #time.sleep(0.01)
            pdbquery = 'https://data.rcsb.org/graphql?query={entry(entry_id:"%s"){polymer_entities {%s{%s}}}}' %(pdb_id,s2,s3)
            r = requests.get(pdbquery)
            info_dic = r.json()
            #parse our needed info from the integrated dict 
            pdb_info = info_dic['data']['entry']['polymer_entities']
            #print(pdb_id)
            for idx in pdb_info:
                if idx[s2] == None:
                    continue
                else:
                    info_list.append(idx[s2][0])
            if len(info_list) != 0:
                break
            
            
        for info in info_list:
            #print(info)
            taxid = info[s3]
            taxid_list.append(taxid)
            #print(taxid)
        #remove duplicates
        taxid_list = list(set(taxid_list))
     
        return taxid_list         
                    
                    



def get_pdb_taxlist(pdb_id):
    taxid_list = get_pdb_info_total(pdb_id)
    if taxid_list == 'PDB data unavailable':
        return 'PDB data unavailable'
    else:
        if len(taxid_list) == 0:
            return 'Not found taxid'
        else:
            return taxid_list
 


# import genome report from NCBI for GC content info
df_pro = pd.read_csv("prokaryotes.txt", sep="\t", low_memory=False)
df_euk= pd.read_csv("eukaryotes.txt", sep="\t", low_memory=False)
df_virus = pd.read_csv("viruses.txt", sep="\t",low_memory=False)


# simplify format of col value 
def col_str_converter(col_str):
    # covert format of taxid in the df
    t1 = col_str.replace('[','')
    t2 = t1.replace(']','')
    t3 = t2.replace("'",'')
    t4 = t3.split(',') #convert t4 as a list
    return t4


# search for genomic GC% in all of the available genome repo databases via TaxID
# taxid might be multiple, so the corresponding GC will be a list
def get_GC_info(taxid_col_str,df1,df2,df3):
    GC_taxlist = []
    taxid_list = col_str_converter(taxid_col_str)
    #print(taxid_list)
    if 'PDB data unavailable'in taxid_list:
        return 'PDB data unavailable'
    else:
        if ('Not found taxid' in taxid_list) or ('None' in taxid_list) :
            return 'taxid not availble'
        else: 
            for TaxID in taxid_list:
                int_taxid = int(float(TaxID))
                if int_taxid in df1['TaxID'].values:
                    GC_df = df1[df1['TaxID'] == int_taxid]['GC%'].values
                    GC = float(GC_df[0])
                    GC_taxlist.append(GC) 
                    #return float(GC)
                else:
                    if int_taxid in df2['TaxID'].values:
                        GC_df = df2[df2['TaxID'] == int_taxid]['GC%'].values
                        GC = float(GC_df[0])
                        GC_taxlist.append(GC) 
                        #return float(GC)
                    else:
                        if int_taxid in df3['TaxID'].values:
                            GC_df = df3[df3['TaxID'] == int_taxid]['GC%'].values
                            GC = float(GC_df[0])
                            GC_taxlist.append(GC) 
                            #return float(GC)
                        else:
                            GC = 'Not found'
                            GC_taxlist.append(GC) 
                            #return ('Not found')
            return GC_taxlist
            
# Function to scrape genomic GC content based on taxnomony
# Match TaxID with the genomic GC content of the organism and link to pdbid
def Merge_pdb_GC(pdb_tax_df, df1,df2,df3):
    GC_list_total = []
    #df = pdb_tax_df.copy()
    pdb_idlist = []
    n = 0
    for idx,row in pdb_tax_df.iterrows():
        n += 1
        taxid = row['TAX_ID']
        pdb_id = row['PDB_ID']
        GC = get_GC_info(taxid, df1, df2, df3)
        
        #taxid_list = get_pdb_taxlist(pdb_id)
        # check if we find out the taxid
    
        #taxid_list.append(taxid)
        GC_list_total.append(GC)
        pdb_idlist.append(pdb_id)
        print(n)
        #report failed search for GC
        if GC == 'taxid not availble':
            print(pdb_id,'FAILED taxid')
        if GC == 'PDB data unavailable':
            print(pdb_id,'FAILED PDB')
    
    pdb_GC_dic = {pdb_idlist[i]: GC_list_total[i] for i in range(len(pdb_idlist))}
    
    return pdb_GC_dic
    

