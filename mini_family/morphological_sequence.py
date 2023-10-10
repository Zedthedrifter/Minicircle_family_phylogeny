"""
    Add expression information to gRNAs and cassettes
"""

import matplotlib.pyplot as plt
import subprocess as sbp
from matplotlib.patches import Rectangle
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import Align
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

def load_database(min_fam_database):
  tmp=pickle_load(min_fam_database)
  return(tmp['mf_similarity'],tmp['mini_family'],tmp['morph_seq'])
#===================================================================================================================================
def compare_mf(mfqs,mfs): 
  matches=0
  for mfq in mfqs:
    hit=0
    for mf in mfs:
      if mfs[mf]['cassettes']==mfqs[mfq]['cassettes']:
        #print(mfs[mf]['cassettes'],mfqs[mfq]['cassettes'])
        hit+=1
        matches+=1
        break
    if hit==0:
      s1=max([sum([1 for i,j in zip(mfs[mf]['cassettes'].split(';'),mfqs[mfq]['cassettes'].split(';')) if i==j and i!='n'])
               /len([i for i in mfs[mf]['cassettes'].split(';') if i!='n']) #divided by the number of cannonical gRNAs in mf1
               for mf in mfs ]) #mf2: mf contained in a certain strain
      s2=max([sum([1 for i,j in zip(mfs[mf]['cassettes'].split(';'),mfqs[mfq]['cassettes'].split(';')) if i==j and i!='n'])
               /len([i for i in mfqs[mfq]['cassettes'].split(';') if i!='n']) #divided by the number of cannonical gRNAs in mf1
               for mf in mfs ]) #mf2: mf contained in a certain strain
      print(f"no exact match found for {mfqs[mfq]['cassettes']},highest similartiy score: {round(s1,3)}-{round(s2,3)}")
  print(f"total minicircle family: {len(mfqs)} \nmwith exact match:{matches} \nwithout exact match:{len(mfqs)-matches}")

#===================================================================================================================================
def by_cassette_position(mfqs,mfs): #must be present at the same cassette position
    similarity={}
    for mf in mfs: #mf1: all Tb mfs
        s=max([sum([1 for i,j in zip(mfs[mf]['cassettes'].split(';'),mfqs[mfq]['cassettes'].split(';')) if i==j and i!='n'])
               /len([i for i in mfs[mf]['cassettes'].split(';') if i!='n']) #divided by the number of cannonical gRNAs in mf1
               for mfq in mfqs ]) #mf2: mf contained in a certain strain
        similarity[mf]=round(s,3)
    return(similarity)

#===================================================================================================================================
def make_morph_seq(mfqs,mfs,similarity,config,rcutoff,note=''):
  seq=[]
  for mf in mfs:
    hit=0
    for mfq in mfqs:
      if mfs[mf]['cassettes']==mfqs[mfq]['cassettes']:
        hit+=1
        seq.append('P')
        break
    if hit==0 and similarity[mf]>=rcutoff:
      seq.append('R')
    elif hit==0 and similarity[mf]<rcutoff:
      seq.append('-')
  identity=f"{config['surfix']}_{config['strain']}.{note}_{config['continent']}_{config['country'].replace(' ','_')}_{config['Year of isolation']}"
  record=SeqRecord(Seq(''.join(seq)),id=identity,description=config['strain'],name=config['strain'])
  return record
      
#===================================================================================================================================
#===================================================================================================================================

def main(config_file='config.yaml'):
    ############################################### FILES #########################################
  config = load_config(config_file)
  in_dir, work_dir = get_directories(config)[:2]
  Tb_pickle_updated=f"{in_dir}/{config['dicts_pickle']}" 
  min_fam_database=f"{in_dir}/{config['minicircle family database']}" 
  gRNAs_q           = f"{work_dir}/{config['adjusted gRNAs']}" #q stands for query
  mini_famq = f"{work_dir}/{config['minicircle family pickle']}"
  rcutoff= config['cutoff for related minicircles']
  #outfiles
  morph_seq = f"{work_dir}/{config['morphological sequence']}" 
    
  #load data available
  mf_similarity,mfs,records=load_database(min_fam_database)
  records2=[]
  for record in records:
    seq=[]
    for i,j in zip (record.seq,mf_similarity[record.description].values()):
      #print(i,j)
      if i=='P':
        seq.append('P')
      elif i=='R' and j>=rcutoff:
        seq.append('R')
      else:
        seq.append('-')
    records2.append(SeqRecord(Seq(''.join(seq)),id=record.id, description=record.description,name=record.name))
  records=records2
  #use gap version instead
  #records=[SeqRecord(record.seq.replace('N','-'),id=record.id, description=record.description,name=record.name)for record in records] #change N to - only
  mfqs=pickle_load(mini_famq)
  compare_mf(mfqs,mfs)
  similarity=by_cassette_position(mfqs,mfs)
  record=make_morph_seq(mfqs,mfs,similarity,config,rcutoff)
  records.append(record)
  #test with reduced minicircle population
  def test_subset():
    mf50={k:mfqs[k] for i,k in enumerate(mfqs) if i<50}
    similarity=by_cassette_position(mf50,mfs)
    record=make_morph_seq(mf50,mfs,similarity,config,rcutoff,'50')
    records.append(record)
    #
    mf100={k:mfqs[k] for i,k in enumerate(mfqs) if i<100}
    similarity=by_cassette_position(mf100,mfs)
    record=make_morph_seq(mf100,mfs,similarity,config,rcutoff,'100')
    records.append(record)
    #
    mf150={k:mfqs[k] for i,k in enumerate(mfqs) if i<150}
    similarity=by_cassette_position(mf150,mfs)
    record=make_morph_seq(mf150,mfs,similarity,config,rcutoff,'150')
    records.append(record)
    #
    mf200={k:mfqs[k] for i,k in enumerate(mfqs) if i<200}
    similarity=by_cassette_position(mf200,mfs)
    record=make_morph_seq(mf200,mfs,similarity,config,rcutoff,'200')
    records.append(record)
    #
    mf250={k:mfqs[k] for i,k in enumerate(mfqs) if i<250}
    similarity=by_cassette_position(mf250,mfs)
    record=make_morph_seq(mf250,mfs,similarity,config,rcutoff,'250')
    records.append(record)
  test_subset()
  #
  SeqIO.write(records,morph_seq,'fasta')
  #use iqtree to make phylogeny
  if config['make phylogeny with default iqtree settings']:
    sbp.call(f"iqtree -s {morph_seq} -st MORPH -redo -bb 10000",shell=True) #can add -m option to specify the model