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


import sys, os
# Get the script's full local path
whereAmI = os.path.dirname(os.path.realpath(__file__))

# set up config_helper
from config.config_helper import Config
config = Config('./config.cfg')

pipeline_dir = config.pipeline_dir

sys.path.append(whereAmI)
sys.path.append(pipeline_dir)

import pipeline_dfci
import utils
import string
import numpy
import os
import re
from collections import defaultdict
#==========================================================================
#============================PARAMETERS====================================
#==========================================================================

# python path
py27_path = config.py27_path

projectName = config.project_name
genome = config.genome
annotFile =  config.annotation_file

#project folders
projectFolder = config.project_folder

#standard folder names
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

# mask Files
maskFile = config.mask_file

# genome firectory
genomeDirectory = config.genome_dir

#making folders
folderList = [gffFolder,macsFolder,macsEnrichedFolder,mappedEnrichedFolder,mappedFolder,wiggleFolder,metaFolder,metaRoseFolder,fastaFolder,figureCodeFolder,figuresFolder,geneListFolder,bedFolder,signalFolder,tableFolder]

for folder in folderList:
    pipeline_dfci.formatFolder(folder,True)



#==========================================================================
#============================LIST OF DATAFILES=============================
#==========================================================================

#this project will utilize multiple datatables
#data tables are organized largely by type/system
#some data tables overlap for ease of analysis




#ATAC-Seq
atac_dataFile = config.data_tables_dict['atac_table']

atac_dataFile_riesling = atac_dataFile.replace('.txt','_riesling.txt')



#==========================================================================
#===========================MAIN METHOD====================================
#==========================================================================


def main():


    print('main analysis for NIBR keratinocyte project atac-seq analysis')

    print('changing directory to project folder')
    os.chdir(projectFolder)

    print('\n\n')
    print('#======================================================================')
    print('#=======================I. CHECKING DATA TABLES========================')
    print('#======================================================================')
    print('\n\n')

    #This section sanity checks each data table and makes sure both bam and .bai files are accessible


    pipeline_dfci.summary(atac_dataFile)


    print('\n\n')
    print('#======================================================================')
    print('#=========================II. MAPPING FASTQS===========================')
    print('#======================================================================')
    print('\n\n')

    # #for atac need no unaligned and no discordant

    # #atac_params =  '--end-to-end --sensitive --no-unal --no-discordant --mm --met-stderr --time'
    # #pipeline_dfci.makeBowtieBashJobsSlurm(atac_dataFile,namesList = [],launch=True,overwrite=False,pCount=16,paramString=atac_params)




    print('\n\n')
    print('#======================================================================')
    print('#======================III. RUNNING RIESLING===========================')
    print('#======================================================================')
    print('\n\n')

    # #riesling is an ATAC-seq pipeline jointly developed by our lab and the Gordon lab at WUSTL

    # #it sanitizes the bams removing duplicate and mitochondrial reads

    # riesling_dir = utils.formatFolder('%sriesling/' % (projectFolder),True)
    # input_dir = utils.formatFolder('/storage/cylin/grail/bam/hg19/NIBR_YvsO/ATACseq_YvsO/')
    # output_dir = utils.formatFolder('/storage/cylin/grail/bam/hg19/NIBR_YvsO/ATACseq_YvsO/riesling/',True)

    # analysis_name = 'NIBR_ATAC_NEW'

    # riesling_bash_path = '%s%s_riesling.sh' % (riesling_dir,analysis_name)

    # riesling_bash = open(riesling_bash_path,'w')

    # riesling_bash.write('#!/usr/bin/bash\n\n')

    # #now write the sbatch headers
    # riesling_bash.write('#SBATCH -n 32\n#SBATCH --mem=512000\n')
    # riesling_bash.write('#SBATCH -o %s%s_reisling_slurm_%%j.out\n' % (riesling_dir,analysis_name))
    # riesling_bash.write('#SBATCH -e %s%s_reisling_slurm_%%j.err\n' % (riesling_dir,analysis_name))
    # riesling_bash.write('pwd; hostname; date\n\n')

    # riesling_bash.write('cd /storage/cylin/bin/riesling-pipeline/\n')
    # riesling_bash.write('%s 2-sanitize-bam.py -i %s -o %s -g %s -v\n' % (py27_path,input_dir, output_dir,genome))
    # riesling_bash.close()

    # print('writing riesling bam commands to %s' % (riesling_bash_path))



    print('\n\n')
    print('#======================================================================')
    print('#=====================IV. FORMATTING RIESLING BAMS=====================')
    print('#======================================================================')
    print('\n\n')

    # atac_table = utils.parseTable(atac_dataFile,'\t')
    # #now fix the path to the right bam
    # for i in range(1,len(atac_table)):
    #     atac_table[i][0] = atac_table[i][0] + 'riesling/'



    # atac_dataFile_riesling = atac_dataFile.replace('.txt','_riesling.txt')
    # new_table = utils.unParseTable(atac_table,atac_dataFile_riesling,'\t')

    # #now we need to index all of the bams
    # dataDict = pipeline_dfci.loadDataTable(atac_dataFile_riesling)

    # names_list = dataDict.keys()
    # bam_directory = dataDict[names_list[0]]['folder']
    # bam_file_list = ['%s%s' % (bam_directory, x) for x in os.listdir(bam_directory) if x.split('.')[-1] == 'bam']
    # print(bam_file_list)
    # for bam_path in bam_file_list:
    #     index_cmd = 'samtools index %s' % (bam_path)
    #     print(index_cmd)
    #     os.system(index_cmd)


    print('\n\n')
    print('#======================================================================')
    print('#========================V. RUNNING PEAK CALLING=======================')
    print('#======================================================================')
    print('\n\n')


    atac_dataFile_riesling = atac_dataFile.replace('.txt','_riesling.txt')
    #run_macs(atac_dataFile_riesling,False)




    print('\n\n')
    print('#======================================================================')
    print('#========================V. RUNNING CLUSTERING ========================')
    print('#======================================================================')
    print('\n\n')


    dataDict = pipeline_dfci.loadDataTable(atac_dataFile_riesling)
    atac_list = dataDict.keys()
    atac_list.sort()
    print(atac_list)


    analysis_name ='keratinocyte_atac'
    cluster_folder = utils.formatFolder('%sclustering' % (projectFolder),True)
    cluster_rose_folder = utils.formatFolder('%sclustering_rose' % (projectFolder),True)
    output_folder = utils.formatFolder('%s%s_clustering' % (cluster_folder,analysis_name),True)
    names_string = ','.join(atac_list)
    cluster_bash_path = '%s%s_clustering.sh' % (cluster_folder,analysis_name)
    cluster_bash = open(cluster_bash_path,'w')
    cluster_bash.write('#!/usr/bin/bash\n\n\n')
    cluster_bash.write('#SBATCH --mem=64000\n\n\n')

    cluster_cmd = '%s %sclusterEnhancer.py -d %s -i %s -r %s -o %s -e super -t 0 -n %s --mask %s' % (py27_path,pipeline_dir,atac_dataFile_riesling,names_string,cluster_rose_folder,output_folder,analysis_name,maskFile)
    cluster_bash.write(cluster_cmd+'\n\n')
    cluster_bash.close()

    # runs only if no output detected
    #run_bash(cluster_bash_path,output_path)




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~RUNNING MACS~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def run_macs(dataFile,useBackground=True):
    dataDict = pipeline_dfci.loadDataTable(dataFile)
    namesList = [name for name in dataDict.keys() if name.upper().count('WCE') ==0 and name.upper().count('INPUT') == 0]
    namesList.sort()
    print(namesList)
    pipeline_dfci.callMacs(dataFile,macsFolder,namesList,False,'1e-9',useBackground)
    os.chdir(projectFolder) # the silly call macs script has to change into the output dir
    #so this takes us back to the project folder

    #to check for completeness, we will try to find all of the peak files
    peak_calling_done = False
    while not peak_calling_done:
        dataDict = pipeline_dfci.loadDataTable(dataFile)
        namesList = [name for name in dataDict.keys() if name.upper().count('WCE') ==0 and name.upper().count('INPUT') == 0]
        for name in namesList:
            peak_path = '%s%s/%s_summits.bed' % (macsFolder,name,name)
            print('searching for %s' % (peak_path))
            if utils.checkOutput(peak_path,1,180):
                peak_calling_done =True
                print('found %s' % (peak_path))
                continue
            else:
                print('Error: peak calling timed out')
                sys.exit()

    #now format the macs output
    print('formatting macs output')
    dataDict = pipeline_dfci.loadDataTable(dataFile)
    namesList = [name for name in dataDict.keys() if name.upper().count('WCE') ==0 and name.upper().count('INPUT') == 0]
    pipeline_dfci.formatMacsOutput(dataFile,macsFolder,macsEnrichedFolder,wiggleFolder,'',useBackground)
    print('Finished running Macs 1.4.2')


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~RUN AND CHECK UTILITY~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#simple helper to run and check a bash script

def run_bash(bash_path,output_path,maxWait=30):
    '''
    runs a bash script and waits up to N minutes
    '''

    if not utils.checkOutput(output_path,0,0):

        print('running bash script %s' % (bash_path))
        os.system('bash %s' % (bash_path))
        if utils.checkOutput(output_path,1,30):
            print('run completed, output detected for %s at %s' % (bash_path,output_path))
    else:
        print('found prior output for %s at %s' % (bash_path,output_path))



#==========================================================================
#==================================THE END=================================
#==========================================================================


if __name__=="__main__":
    main()
