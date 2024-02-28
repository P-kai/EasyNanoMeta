# -*- coding: utf-8 -*-
"""
@author: kai
"""
import pandas as pd

import argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--i','-i', type=str, default='blast_result',help='Input data.')
    parser.add_argument('--data_size','-d',type=float,default='5',help="--data_size, -d")
    parser.add_argument('--title',type=str, default='Sample')
    parser.add_argument('--p','-p',type=str, default='result', help='The prefix of result.')
    parser.add_argument('--output','-o', type=str, default='./',help='Output direction.')
    opt = parser.parse_args()
    

i = opt.i
data_size = opt.data_size
title = opt.title

temp = pd.read_csv(i,sep='\t')

temp['name']=temp['#FILE']+temp['SEQUENCE']

temp.drop(columns=['#FILE','SEQUENCE','STRAND','COVERAGE_MAP','GAPS',
                 '%COVERAGE','%IDENTITY','DATABASE','ACCESSION','PRODUCT'],inplace=True)
order=['name','START','END','GENE','COVERAGE','RESISTANCE']
temp=temp[order]

temp['gene_len']=temp['END']-temp['START']+1

temp_backup=pd.DataFrame(temp)
temp_backup_1=temp_backup['COVERAGE'].str.split("/",expand=True)
temp_backup=temp_backup.join(temp_backup_1)
temp_backup=temp_backup.drop_duplicates("GENE")
temp_backup=pd.DataFrame(temp_backup)
temp_backup=temp_backup.rename(columns={0:'coverage_len',1:'Gene_len','gene_len':'single_gene_len'})
temp_backup.drop(columns=['name','START','END','coverage_len'],inplace=True)
temp_backup=temp_backup.set_index("GENE")

temp1=temp.groupby('GENE')['gene_len'].sum()
temp2=temp.groupby('GENE').size()
temp1=pd.DataFrame(temp1)
temp2=pd.DataFrame(temp2)

temp1=temp1.join(temp2,how='left',on='GENE')
temp1=temp1.rename(columns={0:'gene_num'})
temp_final=temp_backup.join(temp1,how='left')
temp_final['gene_copy_Gb']=temp_final['gene_len'].map(int) / temp_final['Gene_len'].map(int) / data_size
temp_final.drop(columns=['COVERAGE','single_gene_len'],inplace=True)
temp_final.index.name='Gene'
temp_final = temp_final.rename(columns={'RESISTANCE':'Resistance', 'gene_len': 'Sum_gene_len', 'gene_num': 'Gene_num', 
        'gene_copy_Gb': 'Gene_copy/Gb'})
new_columns = [title + '_' + column for column in temp_final.columns]
temp_final.columns = new_columns
pre = opt.p
out_dir = opt.output
temp_final.to_csv(f'{out_dir}/{pre}_abundance', sep='\t')
