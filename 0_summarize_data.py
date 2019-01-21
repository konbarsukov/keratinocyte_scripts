#!/usr/bin/python


'''
The MIT License (MIT)

Copyright (c) 2016 Charles Lin

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
'''


#Main method run script for mycn code

#See README for additional information on downloading and installing dependencies

#Requires linlab pipeline set of utilities
#

#Requires bamliquidator
#

#Requires

#==========================================================================
#=============================DEPENDENCIES=================================
#==========================================================================


import sys, os, string
# Get the script's full local path
whereAmI = os.path.dirname(os.path.realpath(__file__))
print(whereAmI)

# set up config_helper
from config.config_helper import Config
config = Config('./config.cfg')

pipeline_dir = config.pipeline_dir
sys.path.append(pipeline_dir)
print(pipeline_dir)

import pipeline_dfci
import utils
import string
import subprocess
import numpy
import os
import re
from collections import defaultdict
#==========================================================================
#============================PARAMETERS====================================
#==========================================================================



projectName = config.project_name
genome = config.genome
annotFile =  config.annotation_file

# project folders
projectFolder = config.project_folder

# standard folder names
gffFolder = config.gff_folder
macsFolder = config.macs_folder
macsEnrichedFolder = config.macs_enriched_folder
mappedEnrichedFolder = config.mapped_enriched_folder
mappedFolder = config.mapped_folder
wiggleFolder = config.wiggles_folder
metaFolder = config.meta_folder
metaRoseFolder = config.meta_rose_folder
fastaFolder = config.fasta_folder
bedFolder = config.beds_folder
figureCodeFolder = config.figure_code_folder
figuresFolder = config.figures_folder
geneListFolder = config.gene_list_folder
signalFolder = config.signal_tables_folder
tableFolder = config.tables_folder
genePlotFolder = config.gene_plot_folder

# mask Files
maskFile = config.mask_file

# genome firectory
genomeDirectory = config.genome_dir

# gtf file
gtfFile = config.gtf_file

# making folders
folderList = [gffFolder,macsFolder,macsEnrichedFolder,mappedEnrichedFolder,mappedFolder,wiggleFolder,metaFolder,metaRoseFolder,fastaFolder,figureCodeFolder,figuresFolder,geneListFolder,bedFolder,signalFolder,tableFolder, genePlotFolder]

for folder in folderList:
    pipeline_dfci.formatFolder(folder, True)



#==========================================================================
#============================LIST OF DATAFILES=============================
#==========================================================================

#this project will utilize multiple datatables
#data tables are organized largely by type/system
#some data tables overlap for ease of analysis

#ATAC-Seq
atac_dataFile = config.get_data_table('atac_table')

#ChIP-Seq
chip_dataFile = config.get_data_table('chip_table')

#RNA-Seq
rna_dataFile = config.get_data_table('rna_table')



#==========================================================================
#===========================MAIN METHOD====================================
#==========================================================================


def main():


    print('main analysis for MYCN project')

    print('changing directory to project folder')
    os.chdir(projectFolder)

    print('\n\n')
    print('#======================================================================')
    print('#======================I. CHECKING CHIP-SEQ DATA=======================')
    print('#======================================================================')
    print('\n\n')

    #This section sanity checks each data table and makes sure both bam and .bai files are accessible

    #for ChIP-Seq
    #edit all of the data files to absolute path the


    pipeline_dfci.summary(chip_dataFile)



    print('\n\n')
    print('#======================================================================')
    print('#======================II. CHECKING RNA-SEQ DATA=======================')
    print('#======================================================================')
    print('\n\n')

    #This section sanity checks each data table and makes sure both bam and .bai files are accessible

    #for RNA-Seq
    #edit all of the data files to absolute path the


    pipeline_dfci.summary(rna_dataFile)


    #if no processed expression present, runs cuffquant/cuffnorm/RNA-seq pipeline
    cufflinksFolder = utils.formatFolder('%scufflinks' % (projectFolder),True)
    analysis_name = 'NIBR_YvsO'
    young_rna_list = config.get_global_name('young_rna_list')
    old_rna_list = config.get_global_name('old_rna_list')
    print(old_rna_list)
    groupList = [young_rna_list, old_rna_list]
    bashFileName = '%s%s_rna_cufflinks.sh' % (cufflinksFolder,analysis_name)
    pipeline_dfci.makeCuffTable(rna_dataFile,analysis_name,gtfFile,cufflinksFolder,groupList,bashFileName)


    call_bashFileName = 'bash %s' % bashFileName
    proc = subprocess.Popen(call_bashFileName, shell=True)

    # wait for finishing cufflinks
    proc.wait()

    # if call_bashFileName returns 1 (fail), then exit with status 1
    if proc.returncode:
        print 'running %s failed' (call_bashFileName)
        sys.exit(1)



    print('\n\n')
    print('#======================================================================')
    print('#====================III. CHECKING ATAC-SEQ DATA=======================')
    print('#======================================================================')
    print('\n\n')

    #This section sanity checks each data table and makes sure both bam and .bai files are accessible

    #for RNA-Seq
    #edit all of the data files to absolute path the


    pipeline_dfci.summary(atac_dataFile)




#==========================================================================
#==================================THE END=================================
#==========================================================================


if __name__=="__main__":
    main()
