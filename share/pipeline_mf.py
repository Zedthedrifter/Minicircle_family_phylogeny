#!/usr/bin/env python3.8

import mini_family as mf
from os import system
import subprocess as sbp
import argparse

#add the whole package directory to pythonpath:export PYTHONPATH=${PYTHONPATH}:/home/zed/disk1/Tb_Fre_WGS/minicircle_family_pipeline_design/minicircle_family_composition
#add the pipeline directory to path:export PATH=${PATH}:/home/zed/disk1/Tb_Fre_WGS/minicircle_family_pipeline_design/minicircle_family_composition/mini_family
  
  
def main(config_file,work_dir):
  #mf.test()
 # mf.prepare_data(config_file)
 # mf.assign_minicircle_family(config_file)
  mf.morphological_sequence(config_file)


if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='run kDNA annotation')
  parser.add_argument('config_file', help='Name of the config file', metavar='config_file')
  parser.add_argument('work_dir', help='Name of the working directiory', metavar='work_dir')
  options = parser.parse_args()
  
  main(options.config_file, options.work_dir)
  

