"""
This script is meant to go through CCGD runs that have been backed up for a specific project. 
It will attempt to create symlinks for all BreaKmer output directories and files which have been 
generated for each sample in the project.  

Data is stored in two places (on ccgd-stor1):
  Backup:
  /data/backups/mnt/rcgroups/ccgd/Illumina/pipeline_riker/Project

  Active:
  /mnt/rcgroups/ccgd/Illumina/pipeline_riker/Project

Some runs don't have any BreaKmer output calls.
  
Folders are structured like:
[Tumor]-[Normal]/000N/analysis/Breakmer.0.0.6/000N/output/[Tumor]_[BreaKmer_0_0_N]_[type]_svs.out
"""

from glob import glob
import os
import sys
from os.path import join as jp

def setup_links(projectPath, targetPath, sampleFn, subProjectIds) :
    # Write a file containing the samples being analyzed.
    sampleF = open(jp(targetPath, sampleFn), 'w')
    all_samples = {}
    for spDirName in subProjectIds :
        spPath = jp(projectPath, spDirName)
        print spDirName
	for sampleDirName in os.listdir(spPath) :
	    print sampleDirName
	    ## Assume that we always want the output from the latest run
	    latestRun = sorted(os.listdir(jp(spPath, sampleDirName)))[-1]
	    breakmerPath = jp(spPath, sampleDirName, latestRun, 'analysis/Breakmer_0.0.6')
  	    if os.path.exists(breakmerPath) :
   	        latestAnalysis = sorted(os.listdir(breakmerPath))[-1]
	        breakmerPath = jp(breakmerPath, latestAnalysis, 'output')
   	        svsOutFs = glob(jp(breakmerPath, "*svs.out"))
                cmd = 'ln -s ' + breakmerPath + ' ' + jp(targetPath, sampleDirName + '_breakmer_out')
 		os.system(cmd)
                for svsF in svsOutFs :
 	          cmd = 'ln -s ' + svsF + ' ' + targetPath 
                  os.system(cmd)
            else :
              cmd = 'ln -s ' + jp(spPath, sampleDirName, latestRun, 'analysis') + ' ' + targetPath
              print breakmerPath, 'does not exist'
	    sampleF.write(spDirName + "\t" + sampleDirName + "\n")
    sampleF.close()

if __name__ == "__main__" :
    args = sys.argv[1:]
    projectId = args[0]
    subProjectId = args[1]
    print args
#    projectPath = args[0]
#    targetPath = args[1]
#    sampleFn = args[2]

    #Path to the project directory, which usually contains subproject directories.
    # projectPath = '/data/backups/mnt/rcgroups/ccgd/Illumina/pipeline_riker/Project/SGT0023_P50-NSCLC'
    # Destination path for the linked files.
    # targetPath = '/home/rpa4/SGT0023_P50-NSCLC_breakmer'

    projectPath = os.path.join('/ifs/rcgroups/ccgd/Illumina/pipeline_riker/Project', projectId)
    targetPath = os.path.join('/ifs/rcgroups/ccgd/bmw35/analysis/breakmer/breakmer_reports/', projectId + '_breakmer')
    sampleFn = os.path.join(targetPath, subProjectId + '_sample_list.txt')

    subProjectIds = [subProjectId]
    if subProjectId == 'all' :
      subProjectIds = os.listdir(projectPath)
    else:
      subProjectIds = subProjectId.split(',')
    setup_links(projectPath, targetPath, sampleFn, subProjectIds)

