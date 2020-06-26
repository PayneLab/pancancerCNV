#!/usr/bin/env python
# coding: utf-8

# In[2]:


import cptac
import pandas as pd
import pyensembl
import scipy.stats as stats


# In[3]:


def run_fishers_test(df, cancer):
    results = pd.DataFrame(columns=['Mutation', f'{cancer}_odds', f'{cancer}_pvalue'])
    cols = list(df.columns)
    cols.remove('event')
    i = 0
    for col in cols:
        table = pd.crosstab(df[col], df['event'])
        oddsratio, pvalue = stats.fisher_exact(table)
        results.loc[i] = [col, oddsratio, pvalue]
        i += 1
    return results


# In[9]:


brca_mutations_qarm = pd.read_csv("brca_mutations_qarm.csv")
brca_mutations_parm = pd.read_csv("brca_mutations_parm.csv")
colon_mutations_qarm = pd.read_csv("colon_mutations_qarm.csv")
colon_mutations_parm = pd.read_csv("colon_mutations_parm.csv")
hnscc_mutations_qarm = pd.read_csv("hnscc_muations_qarm.csv")
hnscc_mutations_parm = pd.read_csv("hnscc_muations_parm.csv")
lscc_mutations_qarm = pd.read_csv("lscc_mutations_qarm.csv")
lscc_mutations_parm = pd.read_csv("lscc_mutations_parm.csv")
luad_mutations_qarm = pd.read_csv("luad_mutations_qarm.csv")
luad_mutations_parm = pd.read_csv("luad_mutations_parm.csv")
ovarian_mutations_qarm = pd.read_csv("ovarian_mutations_qarm.csv")
ovarian_mutations_parm = pd.read_csv("ovarian_mutations_parm.csv")


# In[7]:


start = time.time()
brca_fishers_q = run_fishers_test(brca_mutations_qarm, 'Brca')
brca_fishers_q.to_csv("brca_fishers_q.tsv", sep='\t')
print(time.time()  - start)


# In[ ]:


start = time.time()
brca_fishers_p = run_fishers_test(brca_mutations_parm, 'Brca')
brca_fishers_p.to_csv("luad_fishers_p.tsv", sep='\t')
print(time.time()  - start)


# In[ ]:


start = time.time()
colon_fishers_q = run_fishers_test(colon_mutations_qarm, 'Colon')
colon_fishers_q.to_csv("colon_fishers_q.tsv", sep='\t')
print(time.time()  - start)


# In[ ]:


start = time.time()
colon_fishers_p = run_fishers_test(colon_mutations_parm, 'Colon')
colon_fishers_p.to_csv("colon_fishers_p.tsv", sep='\t')
print(time.time()  - start)


# In[ ]:


start = time.time()
hnscc_fishers_q = run_fishers_test(hnscc_mutations_qarm, 'Hnscc')
hnscc_fishers_q.to_csv("hnscc_fishers_q.tsv", sep='\t')
print(time.time()  - start)


# In[ ]:


start = time.time()
hnscc_fishers_p = run_fishers_test(hnscc_mutations_parm, 'Hnscc')
hnscc_fishers_p.to_csv("hnscc_fishers_p.tsv", sep='\t')
print(time.time()  - start)


# In[ ]:


start = time.time()
lscc_fishers_q = run_fishers_test(lscc_mutations_qarm, 'Lscc')
lscc_fishers_q.to_csv("lscc_fishers_q.tsv", sep='\t')
print(time.time()  - start)


# In[ ]:


start = time.time()
lscc_fishers_p = run_fishers_test(lscc_mutations_parm, 'Lscc')
lscc_fishers_p.to_csv("lscc_fishers_p.tsv", sep='\t')
print(time.time()  - start)


# In[ ]:


start = time.time()
luad_fishers_q = run_fishers_test(luad_mutations_qarm, 'Luad')
luad_fishers_q.to_csv("luad_fishers_q.tsv", sep='\t')
print(time.time()  - start)


# In[ ]:


start = time.time()
luad_fishers_p = run_fishers_test(luad_mutations_parm, 'Luad')
luad_fishers_p.to_csv("luad_fishers_p.tsv", sep='\t')
print(time.time()  - start)


# In[ ]:


start = time.time()
ovarian_fishers_q = run_fishers_test(ovarian_mutations_qarm, 'Ovarian')
ovarian_fishers_q.to_csv("ovarian_fishers_q.tsv", sep='\t')
print(time.time()  - start)


# In[ ]:


start = time.time()
ovarian_fishers_p = run_fishers_test(ovarian_mutations_parm, 'Ovarian')
ovarian_fishers_p.to_csv("ovarian_fishers_p.tsv", sep='\t')
print(time.time()  - start)

