# -*- coding: utf-8 -*-
"""
@author: kai
"""

import pandas as pd
import argparse
from itertools import combinations

import argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--i', type=str, default='blast_result')
    opt = parser.parse_args()
        

i = opt.i

df = pd.read_csv(i,sep='\t')

df.drop(columns=['#FILE', 'START', 'END', 'STRAND', 'COVERAGE','COVERAGE_MAP', 'GAPS',
                     '%COVERAGE', '%IDENTITY', 'DATABASE', 'ACCESSION', 'PRODUCT','RESISTANCE'], inplace=True)

df['Seq_gene'] = df['SEQUENCE'] + df['GENE']

df.drop_duplicates(subset='Seq_gene',inplace=True)
df = pd.DataFrame(df).set_index('SEQUENCE')
temp1 = df.groupby('SEQUENCE').size()
temp1 = pd.DataFrame(temp1,columns=['Seq_num'])

df = df.join(temp1, how='left')
df=df[df["Seq_num"]>1]

# 按SEQUENCE分组
grouped = df.groupby('SEQUENCE')

# 定义一个函数，对每个组应用combn函数
def combn_func(group):
    unique_genes = group['GENE'].unique()
    return list(combinations(unique_genes, 2))

# 使用group_map应用函数到每个分组并将结果转换为DataFrame
result = pd.DataFrame(grouped.apply(combn_func).explode().tolist(), columns=['GENE1', 'GENE2'])
# 输出结果

result['Genes'] = result['GENE1'] + "/" + result['GENE2']
result = result.groupby('Genes').size()
result = pd.DataFrame(result)
result = result.rename(columns={0: 'gene_num'})
result = result.reset_index()

result1 = result['Genes'].str.split("/", expand=True)

result1 = result1.rename(columns={0: 'GENE1', 1: 'GENE2' })

result1 = result.join(result1, how="left")

result1['Dump'] = result1.apply(lambda row: str(row['GENE1']) + "/" + str(row['GENE2']) if row['GENE1'] > row['GENE2'] else str(row['GENE2']) + "/" + str(row['GENE1']), axis=1)
result1 = result1.groupby('Dump')['gene_num'].sum()

result1.to_csv('./co-located_result', sep='\t')
