"""
    Add expression information to gRNAs and cassettes
"""

import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from Bio import SeqIO
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

def load_database(Tb_pickle_updated):
  tmp=pickle_load(Tb_pickle_updated)
  return(tmp[0],tmp[1],tmp[2],tmp[3])
#===================================================================================================================================

#get deletion and editing seuqences
def read_editing_sites(small_u):
    edits={}
    for k in small_u:
        sites=[0 if b=='u' else 1 for b in small_u[k]] #0 for insertions
        edits[k]=sites
    return(edits)
#deletion
def read_deletion_sites(deletion):
    delets={}
    for k in deletion:
        sites=[0 if b!='-' else 1 for b in deletion[k]] #0 for deletion sites
        delets[k]=sites
    return(delets)
#
#initiation site counts
def initiation_site_counts(gRNA_dict,small_u,ws=1):
    insertions=read_editing_sites(small_u)
    mrna0={k:[0]*(len(insertions[k])+20) for k in insertions} #extend the mRNA a bit longer for the initiation sites
    for k in gRNA_dict:
        init_site=gRNA_dict[k]['gene_mRNA_end']
        mrna0[gRNA_dict[k]['mRNA_name']][init_site]+=1
    counts={k:[sum(mrna0[k][i:i+ws]) for i in range(len(mrna0[k])-ws)] for k in mrna0}
    return(counts)
    
#define gRNA family boundaries
def gRNA_family_boundaries(counts):
    gRNA_families={}
    for k in counts:
        #find boundaries
        mrna=''.join([str(i) for i in counts[k]])
        starts=[match.start() for match in re.finditer('0{3,}',mrna)] #change from 3 to 1
        ends=[match.end() for match in re.finditer('0{3,}',mrna)]
        boundaries=[(i,j) for i,j in zip(ends[:-1],starts[1:])]
        gRNA_families[k]={f'{k}-{pair[0]}_{pair[1]}':{'bound':pair,'gRNA':[],'mRNA':k} for i,pair in enumerate(boundaries)}#
        gRNA_families={m:{k: v for k, v in sorted(gRNA_families[m].items(), key=lambda item: item[1]['bound'][0])} for m in gRNA_families}
        #print(gRNA_families)            
    return(gRNA_families)

#assign gRNA families
def assign_gRNA_families(gRNA_dict,gRNA_families,lext,rext): #gRNA_dict from query, gRNA_families from database
    import copy
    gRNA_familiesq,new,count=copy.deepcopy(gRNA_families),{k:{} for k in gRNA_families},0
    for k in gRNA_familiesq:
      for gf in gRNA_familiesq[k]:
        gRNA_familiesq[k][gf]['gRNA']=[]
    for k in gRNA_dict:
        gRNA_dict[k]['gRNA_family']=[]
        hit=0
        for f in gRNA_families[gRNA_dict[k]['mRNA_name']]:
            bound=gRNA_families[gRNA_dict[k]['mRNA_name']][f]['bound']
            if int(gRNA_dict[k]['gene_mRNA_end']) >= bound[0]-lext and int(gRNA_dict[k]['gene_mRNA_end']) <= bound[1]+rext:
                gRNA_dict[k]['gRNA_family'].append(f)
                gRNA_families[gRNA_dict[k]['mRNA_name']][f]['gRNA'].append(k) #add the gRNA from query dataset to the main dataset
                gRNA_familiesq[gRNA_dict[k]['mRNA_name']][f]['gRNA'].append(k)
                hit+=1
                break
        if hit==0:
            count+=1
            #print(f'{k},{gRNA_dict[k]} has not been represented in the dataset')
            mrna=gRNA_dict[k]['mRNA_name']
            newgf=f"{mrna}_{gRNA_dict[k]['gene_mRNA_end']}_{gRNA_dict[k]['gene_mRNA_end']+1}"
            print(newgf)
            #gRNA_families[mrna][newgf]={}
            tmp={}
            tmp['bound']=(gRNA_dict[k]['gene_mRNA_end'],gRNA_dict[k]['gene_mRNA_end']+1)
            tmp['gRNA']=[k]
            tmp['mRNA']=mrna
            tmp['minicircle_family']=[]
            gRNA_families[mrna][newgf]=gRNA_familiesq[mrna][newgf]=tmp
            gRNA_dict[k]['gRNA_family'].append(newgf) #also update the info to gRNA_dict so that the new gRNAs can be recognized later
            
    print(f"{count} gRNAs are not repreesnted by the database")        
    #get gRNA family dict that contains only the gf present in the query dataset
    for k in gRNA_familiesq:
      for gf in gRNA_familiesq[k]:
        if gRNA_familiesq[k][gf]['gRNA']!=[]:
          new[k][gf]=gRNA_familiesq[k][gf]
    return(gRNA_dict,gRNA_families,new)

#make sure all gRNAs are assigned
def check_gRNA_number(gRNA_dict,gRNA_families):
    in_fam=sum([len(gRNA_families[k][gf]['gRNA']) for k in gRNA_families for gf in gRNA_families[k]])
    print(f'number of gRNA in gRNA family dict: {in_fam} \nnumber of gRNA in gRNA dict:{len(gRNA_dict)}')

#visualize the boundaries
def get_gf_presence(combined,gRNA_fam,gRNA_dict):
  gf_presence={}
  mini_presence={mini: [idx for idx in combined[mini][~(combined[mini].isna())].index] for mini in combined.columns}
  for k in gRNA_fam:
    for gf in gRNA_fam[k]:
        total,neo=[],0
        for g in gRNA_fam[k][gf]['gRNA']:
            try:
                mini=gRNA_dict[g]['mO_name']
                total+=mini_presence.get(mini,[])
            except:
                neo+=1
        total=set(total)
        gf_presence[gf]=round((len(total)+neo)/(neo+len(combined)),2) #the ratio is used to color the patches for each gf
  return(gf_presence)

def plot_initiation_sites(counts,gRNA_fam,gRNA_dict,gf_presence,outfile,figw=100,figh=50):
    c,c_dict=0,{'I':'r','II':'skyblue','III':'orange','IV':'purple','V':'green','Orphan':'black'}
    print(c_dict)
    fig,axs = plt.subplots(len(counts),1,figsize=(figw,figh))
    for k in counts:
        ax=axs[c]
        if 'v' in k and 'v1' not in k:
          ax.set_title(f'Unique Initiation site starts on {k}',fontdict={'fontsize': 30, 'fontweight': 30})
        else: 
          ax.set_title(f'Initiation site starts on {k}',fontdict={'fontsize': 30, 'fontweight': 30})
        c+=1
        #plot init sites
        ax.plot(range(len(counts[k])),counts[k],label=k)
        ax.axhline(y=1, xmin=0, xmax=len(counts[k]),ls='--')
        ax.set_ylim(-2,max(counts[k])+1)
        #plot boundaries
        for gf in gRNA_fam[k]:
            gc={}
            ax.add_patch(Rectangle((gRNA_fam[k][gf]['bound'][0], 0), gRNA_fam[k][gf]['bound'][1]-gRNA_fam[k][gf]['bound'][0], max(counts[k])+1,
                                  facecolor = 'black',alpha=gf_presence[gf]*0.3)) #need to change patch color to show relative prevalence 
            ax.text(gRNA_fam[k][gf]['bound'][0],0,f"{gRNA_fam[k][gf]['bound'][0]}\nmf#:{len(gRNA_fam[k][gf].get('minicircle_family',['Maxi']))}")
            for g in gRNA_fam[k][gf]['gRNA']:
              if g in gRNA_dict:
                  if gRNA_dict[g]['gene_mRNA_end'] in gc:
                      gc[gRNA_dict[g]['gene_mRNA_end']].append(gRNA_dict[g]['cassette_label'])
                  else:
                      gc[gRNA_dict[g]['gene_mRNA_end']]=[]
                      gc[gRNA_dict[g]['gene_mRNA_end']].append(gRNA_dict[g]['cassette_label'])
            gc={k:gc[k] for k in gc} #plot the same position twice
            for init in gc:
                h=-0.4
                for cass in gc[init]:
                    if cass!='Maxi':
                        ax.scatter([init],[h],color=c_dict[cass])
                        h+=(-0.2)
    plt.savefig(outfile)
    plt.show()

#
#character minicircles by gRNA family
def assign_gRNA_fam_to_mini(gRNA_dict,gRNA_families): #use the new version
    mini_dict,cassettes,non_cannonical={},['I','II','III','IV','V','Orphan'],[]
    for g in gRNA_dict:
        if gRNA_dict[g]['mO_name'] not in mini_dict:
            mini_dict[gRNA_dict[g]['mO_name']]={c:{} for c in cassettes}
            mini_dict[gRNA_dict[g]['mO_name']][gRNA_dict[g]['cassette_label']]=gRNA_dict[g]
        else:
            mini_dict[gRNA_dict[g]['mO_name']][gRNA_dict[g]['cassette_label']]=gRNA_dict[g]
    mini_dict.pop('Maxicircle','None') #remove maxicircle      
    mini_df = pd.DataFrame(columns=cassettes,index=mini_dict.keys())
    for mini in mini_dict:
        for c in mini_dict[mini]:
            try:
                mini_df.loc[[mini],[c]]=mini_dict[mini][c]['gRNA_family'][0]
            except:
                non_cannonical.append(f"no cannonical gRNA found in cassette {c} in minicircle {mini}")
    mini_df=mini_df.sort_values(by=cassettes)
    return(mini_dict,mini_df)
#group minicircles by gRNA family
def make_minicircle_family(mini_dict):
    mini_families={}
    for mini in mini_dict:
        hit=0
        #create a unique descriptor for each minicircle
        try:
            cs=';'.join([mini_dict[mini][c].get('gRNA_family','n')[0] for c in mini_dict[mini]])
            if cs not in mini_families:
                mini_families[cs]=[]
                mini_families[cs].append(mini)
            else:
                mini_families[cs].append(mini)
        except:
            print(f'error for {mini} \n {mini_dict[mini]}')
    #sort by the number of encoded gRNAs
    tmp={k:k.count('n') for k in mini_families}
    tmp={k: v for k, v in sorted(tmp.items(), key=lambda item: item[1])}
    ordered=list(tmp.keys())
    mini_families={f"neo_mf_{i}":{'cassettes':k,'class':mini_families[k]} for i,k in enumerate(ordered)}
    return(mini_families)

##==============================================================================================================================

def main(config_file='config.yaml'):
    ############################################### FILES #########################################
  config = load_config(config_file)
  in_dir, work_dir = get_directories(config)[:2]
  Tb_pickle_updated=f"{in_dir}/{config['dicts_pickle']}" 
  unedited_in_file  = f"{in_dir}/{config['unedited mRNA fasta infile']}"
  edited_in_file    = f"{in_dir}/{config['edited mRNA fasta infile']}"
  gRNAs_q           = f"{work_dir}/{config['adjusted gRNAs']}" #q stands for query
  lext= config['extend boundary left']
  rext= config['extend boundary right']
  #outfiles
  initi_plot        = f"{work_dir}/{config['surfix']}_initiation_site_starts.pdf"
  mini_families_out = f"{work_dir}/{config['minicircle family pickle']}"
    
  #load all data available
  gRNA_dict,gRNA_families,mini_dict,renamed_mini_family=load_database(Tb_pickle_updated)
  edited=SeqIO.to_dict(SeqIO.parse(edited_in_file,'fasta'))
  metadf2,combined=process_df(f"{in_dir}/{config['minicircle copy number']}",f"{in_dir}/{config['metadata']}")
  ###get gRNAs
  gRNA_dictq=pickle_load(gRNAs_q)
  gRNA_dictq,gRNA_families,gRNA_familiesq=assign_gRNA_families(gRNA_dictq,gRNA_families,lext,rext) #using the gRNA families from database
  check_gRNA_number(gRNA_dictq,gRNA_familiesq)
  #visualize
  if config['visualize initiation site starting points']:
    #get starts of initiation sites and anchors
    counts=initiation_site_counts(gRNA_dictq,edited)
    gf_presence=get_gf_presence(combined,gRNA_families,gRNA_dict)
    plot_initiation_sites(counts,gRNA_families,gRNA_dictq,gf_presence,initi_plot,figw=100,figh=50)
  #
  mini_dictq,mini_dfq=assign_gRNA_fam_to_mini(gRNA_dictq,gRNA_familiesq)
  mini_families=make_minicircle_family(mini_dictq)
  pickle_out(mini_families,mini_families_out)
  #print(mini_families)