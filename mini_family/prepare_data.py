"""
    Add expression information to gRNAs and cassettes
"""

import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from Bio import SeqIO
from Bio import Align
from Bio.Align import substitution_matrices
import pandas as pd
import numpy as np
pd.set_option('display.max_rows', None)
pd.set_option('display.max_columns', None)
from sklearn.preprocessing import MinMaxScaler
import re
import os
import copy
import pickle
import gzip
import math
import networkx as nx
import matplotlib.pyplot as plt
from matplotlib import cm
import pygraphviz
import datetime
date = datetime.datetime.now().strftime("%Y-%m-%d")

from .common import *

def load_database(Tb_pickle_updated):
  tmp=pickle_load(Tb_pickle_updated)
  return(tmp[0],tmp[1],tmp[2],tmp[3])
#===================================================================================================================================

  
def get_gRNA_info (txt,adjust=0): #build dictionary from gRNA_with_expression files
    with open(txt) as handle:
        keys=[x for x in next(handle).strip('\n').split(' ') if x !='']
        gRNA_dict={f'neo_gRNA_{i}': {i:j for i,j in zip(keys,[x for x in line.strip('\n').split(' ') if x !=''][1:])}for i,line in enumerate(handle)}
        for key in gRNA_dict:
            init_site=int(gRNA_dict[key]['mRNA_end'])+int(gRNA_dict[key]['rel_pos'])
            if gRNA_dict[key]['mRNA_name'] not in ['CYB','MURF2']:
                gRNA_dict[key]['gene_mRNA_end']=init_site+adjust
            else:
                gRNA_dict[key]['gene_mRNA_end']=init_site
    return(gRNA_dict)  
#
def combine_alternatives(gRNA_dict):
    alt,paired={},[]
    other={g:gRNA_dict[g] for g in gRNA_dict if '_v' not in gRNA_dict[g]['mRNA_name']}
    v1={g:gRNA_dict[g] for g in gRNA_dict if '_v1' in gRNA_dict[g]['mRNA_name']}
    v2={g:gRNA_dict[g] for g in gRNA_dict if '_v2' in gRNA_dict[g]['mRNA_name']}
    print('number of gRNAs for v1 and v2 mRNAs',len(v1),len(v2))
    for g1 in v1:
        for g2 in v2:
            if v1[g1]['mO_name']==v2[g2]['mO_name'] and v1[g1]['cassette_label']==v2[g2]['cassette_label'] and v1[g1]['product']==v2[g2]['product']:
                dcs=int(v1[g1]['circle_start'])-int(v2[g2]['circle_start'])
                dce=int(v1[g1]['circle_end'])-int(v2[g2]['circle_end'])
                if dcs>2 or dce>2:
                    #this is just for manual checking
                    print('paired gRNA with >2bp difference on minicircle:',v1[g1]['name'],v2[g2]['name'])
                paired.append(g2)
    v2={k:v2[k] for k in v2 if k not in paired}
    print('number of gRNAs in v1 and only in v2',len(v1),len(v2))
    print('\ndistinct gRNAs for alternatively edited mRNAs:')
    for k in v2:
        print(v2[k]['name'])
    v1.update(v2)
    v1.update(other)
    return(v1)
#adjust the gRNA annotation so that the relative positions on protein sequence/mRNA sequence could match
#align mRNA sequence based on protein sequence?
def orf_translate(mRNA):
    max_lengths = [max([len(j) for j in re.split('\*', str(mRNA[orf:].translate(table=4)))]) for orf in range(3)]
    orf = max_lengths.index(max(max_lengths))
    translate=mRNA[orf:].translate(table=4)
    return(translate,orf)
#for aligning mRNAs by protein sequences
def adjust_mRNA_by_protein_seq(m,p,orf):
    #check if the number of bases is correct:
    if len(m)-len(p.replace('-',''))*3<=4:
        #print('numbers of neucldotide and amino acids agree')
        i,mrna=0,'-'*(3-orf)+m[:orf]
        for a in p.rstrip('-'):
            if a!='-':
                mrna+=m[i*3+orf:(i+1)*3+orf]
                i+=1
            else:
                mrna+='---'
        mrna+=m[(i+1)*3+orf:]
        mrna+=m.replace(mrna.replace('-',''),'')
        if m==mrna.replace('-',''):
            #print('add gaps correctly')
            return(mrna)
        else:
            print('mRNA adjustment error')
            return(mrna)
    else:
        print('numbers of neucldotide and amino acids do not agree')
#for aligned mRNAs
def count_gaps(m):
    i,gaps=0,{}
    for j,b in enumerate(m):
        if b!='-':
            gaps[i]=m[:j].count('-')
            i+=1
    return(gaps)
#
#alignment based on amino acid sequence conservation
def align_protein(m1,m2):
    p1,orf1=orf_translate(m1)
    p2,orf2=orf_translate(m2)
    aligner = Align.PairwiseAligner(scoring="blastp")
    alignments = aligner.align(p1, p2)
    a=alignments[0]
    align2=''.join([([s for s in l.split(' ') if s !='']+[''])[2] for l in a.format().strip('\n').split('\n') if 'query' in l])
    align1=''.join([([s for s in l.split(' ') if s !='']+[''])[2] for l in a.format().strip('\n').split('\n') if 'target' in l])
    #should use the aligned protein sequence
    m1=adjust_mRNA_by_protein_seq(m1,align1,orf1)
    m2=adjust_mRNA_by_protein_seq(m2,align2,orf2)
    gap1=count_gaps(m1)
    gap2=count_gaps(m2)
    return(align1,align2,m1,m2,gap1,gap2)
#
def align_protein_iter(smallu1,smallu2,s1,s2,aligntxt):
    smallu1=SeqIO.to_dict(SeqIO.parse(smallu1,'fasta'))
    smallu2=SeqIO.to_dict(SeqIO.parse(smallu2,'fasta'))
    gaps1,gaps2,align1,align2={},{},{},{}
    f=open(aligntxt,'w')
    for k in smallu1:
        if k in smallu2:
            print(f'{k} is found in both strains')
            m1=smallu1[k].seq
            m2=smallu2[k].seq
            p1,p2,m1,m2,gap1,gap2=align_protein(m1,m2)
            gaps1[k]=gap1
            gaps2[k]=gap2
            align1[k]=m1
            align2[k]=m2
            pairs=''.join(['|' if i==j and i!='-' else ' ' for i,j in zip(m1.upper(),m2.upper())])
            f.write(f">{k}\t{s1}\tvs\t{s2}\n-  {'  '.join(list(p1))}\n{m1}\n{pairs}\n{m2}\n-  {'  '.join(list(p2))}\n\n")
        else:
            gaps1[k]={i:0 for i in range(len(smallu1[k].seq))}
            align1[k]=str(smallu1[k].seq)
    f.close()
        #print(f'done{k}')
    for k in smallu2:
        if k not in smallu1:
            gaps2[k]={i:0 for i in range(len(smallu2[k].seq))}
            align2[k]=str(smallu2[k].seq)
    insertion1=read_editing_sites(align1)
    insertion2=read_editing_sites(align2)
    return(gaps1,gaps2,insertion1,insertion2)
#
def adjust_gRNA_by_alignment(gRNA_dict,gaps):
    import copy
    new=copy.deepcopy(gRNA_dict) #two dictionaries with the same values 
    for g in gRNA_dict:
        if new[g]['mRNA_name'] in gaps:
            s=int(gRNA_dict[g]['mRNA_start'])
            e=int(gRNA_dict[g]['mRNA_end'])
            ge=int(gRNA_dict[g]['gene_mRNA_end'])
            b=max(list(gaps[gRNA_dict[g]['mRNA_name']].keys()))
            sk=min([int(gRNA_dict[g]['mRNA_start']),b])
            ek=min([int(gRNA_dict[g]['mRNA_end']),b])
            gek=min([int(gRNA_dict[g]['gene_mRNA_end']),b])
            try:
                new[g]['mRNA_start']=s+gaps[new[g]['mRNA_name']][sk]
                new[g]['mRNA_end']=e+gaps[new[g]['mRNA_name']][ek]
                new[g]['gene_mRNA_end']=ge+gaps[new[g]['mRNA_name']][gek]
            except:
                print(new[g]['mRNA_name'],s,e)
    return(new)

##==============================================================================================================================

def main(config_file='config.yaml'):
    ############################################### FILES #########################################
  config = load_config(config_file)
  in_dir, work_dir = get_directories(config)[:2]
  unedited_in_file  = f"{in_dir}/{config['unedited mRNA fasta infile']}"
  edited_in_file    = f"{in_dir}/{config['edited mRNA fasta infile']}"
  gRNAs_q           = f"{in_dir}/{config['gRNAs for strain of interest']}" #q stands for query
  #data provided
  Tb_pickle_updated=f"{in_dir}/{config['dicts_pickle']}" 
  eatro_edited     = f"{in_dir}/{config['EATRO1125 edited mRNA fasta']}"
  tbg1_edited      = f"{in_dir}/{config['Tbg1 edited mRNA fasta']}"
  #outfiles
  initi_plot        = f"{work_dir}/{config['strain']}_initiation_site_starts.pdf"
  gRNA_dictq_adj    =f"{work_dir}/gRNA_dict.adjusted.pickle"
  
  #load all data available
  edited=SeqIO.to_dict(SeqIO.parse(edited_in_file,'fasta'))
  ###get gRNAs
  gRNA_dictq=get_gRNA_info (gRNAs_q,adjust=0)
  gRNA_dictq=combine_alternatives(gRNA_dictq)
  #print(gRNA_families)
  #if config['species']!='Trypanosoma brucei gambiense type 1':
  #  gaps1,gaps2,insertion1,insertion2=align_protein_iter(edited_in_file,tbg1_edited,config['strain'],'Tbg1',f"{work_dir}/{config['surfix']}_Tbg1_aligned.txt")
  #  print('align to Tbg1 edited mRNAs for adjustments')
  #else:
  #  gaps1,gaps2,insertion1,insertion2=align_protein_iter(edited_in_file,eatro_edited,config['strain'],'EATRI1125',f"{work_dir}/{config['surfix']}_EATRO1125_aligned.txt")
  #  print('align to EATRO1125 edited mRNAs for adjustments')
  gaps1,gaps2,insertion1,insertion2=align_protein_iter(edited_in_file,eatro_edited,config['surfix'],'EATRI1125',f"{work_dir}/{config['surfix']}_EATRO1125_aligned.txt")
  gRNA_dictq=adjust_gRNA_by_alignment(gRNA_dictq,gaps1)
  pickle_out(gRNA_dictq,gRNA_dictq_adj)