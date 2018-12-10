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
whereAmI = os.path.dirname(os.path.realpath(__file__)) +'/'
print(whereAmI)

pipeline_dir = '/'.join(whereAmI.split('/')[0:-2]) + '/pipeline/'
sys.path.append(pipeline_dir)
print(pipeline_dir)

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



projectName = 'NIBR_YvsO_cyl'
genome ='hg19'
annotFile = '%s/annotation/%s_refseq.ucsc' % (pipeline_dir,genome)

#project folders
projectFolder = '/storage/cylin/grail/projects/%s/' % (projectName) #PATH TO YOUR PROJECT FOLDER

#standard folder names
gffFolder ='%sgff/' % (projectFolder)
macsFolder = '%smacsFolder/' % (projectFolder)
macsEnrichedFolder = '%smacsEnriched/' % (projectFolder)
mappedEnrichedFolder = '%smappedEnriched/' % (projectFolder)
mappedFolder = '%smappedFolder/' % (projectFolder)
wiggleFolder = '%swiggles/' % (projectFolder)
metaFolder = '%smeta/' % (projectFolder)
metaRoseFolder = '%smeta_rose/' % (projectFolder)
fastaFolder = '%sfasta/' % (projectFolder)
bedFolder = '%sbed/' % (projectFolder)
figureCodeFolder = '%sfigureCode/' % (projectFolder)
figuresFolder = '%sfigures/' % (projectFolder)
geneListFolder = '%sgeneListFolder/' % (projectFolder)
bedFolder = '%sbeds/' % (projectFolder)
signalFolder = '%ssignalTables/' % (projectFolder)
tableFolder = '%stables/' % (projectFolder)

#mask Files
maskFile ='%smasks/hg19_encode_blacklist.bed' % (projectFolder)

#genomeDirectory
genomeDirectory = '%s/genomes/Homo_sapiens/UCSC/hg19/Sequence/Chromosomes/' % (projectFolder)
gtfFile = '%s/gtf/genes.gtf' % (projectFolder)

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
atac_dataFile = '%sdata_tables/NIBR_CHIP_TABLE.txt' % (projectFolder)

#ChIP-Seq
chip_dataFile = '%sdata_tables/NIBR_ATAC_TABLE.txt' % (projectFolder)

#RNA-Seq
rna_dataFile = '%sdata_tables/NIBR_RNA_TABLE.txt' % (projectFolder)



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
    groupList = [['Y_BC10_Y1','Y_BC11_Y2','Y_BC16_Y3'],['O_BC18_O1','O_BC25_O2','O_BC27_O3']]
    bashFileName = '%s%s_rna_cufflinks.sh' % (cufflinksFolder,analysis_name)
    pipeline_dfci.makeCuffTable(rna_dataFile,analysis_name,gtfFile,cufflinksFolder,groupList,bashFileName)


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
