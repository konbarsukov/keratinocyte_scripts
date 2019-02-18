#!/usr/bin/python


'''
The MIT License (MIT)

Copyright (c) 2017 Charles Lin

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


#Main method run script for NIBR YvsO project

#1_map_enhancers.py runs enhancer calling and finds BRD4 subpeaks

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

pipeline_dir = config.pipeline_folder

sys.path.append(whereAmI)
sys.path.append(pipeline_dir)

from helpers import mergeMacs14, call_crc, tf_edge_delta_out, tf_edge_delta_in, makeDegreeTable, get_tf_target_genes
import pipeline_dfci
import utils
import string
from subprocess import Popen
import numpy
from scipy import stats
import os
import re
from collections import defaultdict
#==========================================================================
#============================PARAMETERS====================================
#==========================================================================



projectName = config.project_name
genome = config.genome_name
analysis_name = config.get_analysis_name('crc_analysis')
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
genePlotFolder = config.gene_plot_folder

# mask Files
maskFile = config.mask_file

# genome firectory
genomeDirectory = config.genome_folder

#making folders
folderList = [gffFolder,macsFolder,macsEnrichedFolder,mappedEnrichedFolder,mappedFolder,wiggleFolder,metaFolder,metaRoseFolder,fastaFolder,figureCodeFolder,figuresFolder,geneListFolder,bedFolder,signalFolder,tableFolder,genePlotFolder]

for folder in folderList:
    pipeline_dfci.formatFolder(folder,True)

py27_path = config.py27_path
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
#===========================GLOBAL NAMES LISTS=============================
#==========================================================================

group1_k27ac_list = config.get_global_name('young_h3k27ac_list')
group2_k27ac_list = config.get_global_name('old_h3k27ac_list')
k27ac_list = group1_k27ac_list + group2_k27ac_list

group1_brd4_list = config.get_global_name('young_brd4_list')
group2_brd4_list = config.get_global_name('old_brd4_list')

group1_atac_list = config.get_global_name('young_atac_list')
group2_atac_list = config.get_global_name('old_atac_list')

chip_list_no_input = config.get_global_name('chip_list_no_input')


#for some reason needed to re-run macs for that brd4
#pipeline_dfci.callMacsSlurm(chip_dataFile,macsFolder,namesList = ['O_BRD4_1'],overwrite=False,pvalue='1e-9',useBackground =True,launch=False)

#==========================================================================
#===========================MAIN METHOD====================================
#==========================================================================


def main():


    print('main analysis for {} project'.format(projectName))

    print('changing directory to project folder')
    os.chdir(projectFolder)

    print('\n\n')
    print('#======================================================================')
    print('#========================I. CHECKING SEQ DATA==========================')
    print('#======================================================================')
    print('\n\n')

    #This section sanity checks each data table and makes sure both bam and .bai files are accessible

    #for ChIP-Seq
    #edit all of the data files to absolute path the

    #for ChIP-Seq
    pipeline_dfci.summary(chip_dataFile)

    #for RNA-Seq
    pipeline_dfci.summary(rna_dataFile)

    #for atac-seq
    pipeline_dfci.summary(atac_dataFile)




    print('\n\n')
    print('#======================================================================')
    print('#====================II. DETECT ENHANCER CLUSTERING====================')
    print('#======================================================================')
    print('\n\n')

    ##enhancer clustering is run by the 1_map_enhancers.py script

    #analysis_name ='keratinocyte_se'
    #cluster_folder = utils.formatFolder('%sclustering' % (projectFolder),True)
    #output_folder = utils.formatFolder('%s%s_clustering' % (cluster_folder,analysis_name),True)
    #output_path = '%s%s_%s_clusterTable.txt' % (output_folder,genome.upper(),analysis_name)

    #if utils.checkOutput(output_path,0.1,0.1):
    #    print('IDENTIFIED SE CLUSTERING AT %s' % (output_path))
    #else:
    #    print('ERROR: NO SE CLUSTERING AT %s' % (output_path))
    #    sys.exit()


    #analysis_name ='keratinocyte_all'
    #cluster_folder = utils.formatFolder('%sclustering' % (projectFolder),True)
    #output_folder = utils.formatFolder('%s%s_clustering' % (cluster_folder,analysis_name),True)
    #output_path = '%s%s_%s_clusterTable.txt' % (output_folder,genome.upper(),analysis_name)

    #if utils.checkOutput(output_path,0.1,0.1):
    #    print('IDENTIFIED ALL ENHANCER CLUSTERING AT %s' % (output_path))
    #else:
    #    print('ERROR: NO ALL ENHANCER CLUSTERING AT %s' % (output_path))
    #    sys.exit()



    print('\n\n')
    print('#======================================================================')
    print('#===================III. MERGING ATAC PEAKS BY Y AND O=================')
    print('#======================================================================')
    print('\n\n')

    #for now assumes pre-run macs14 output from riesling (3_process_atac.py)


    #now write out some paths
    group1_bed_path = '%s%s_group1_atac_-0_+0.bed' % (bedFolder, genome)
    group2_bed_path = '%s%s_group2_atac_-0_+0.bed' % (bedFolder, genome)

    #combined
    combined_bed_path = '%s%s_combined_atac_-0_+0.bed' % (bedFolder, genome)

    #look for these to already exist
    if utils.checkOutput(group1_bed_path,0.1,0.1) and utils.checkOutput(group2_bed_path,0.1,0.1) and utils.checkOutput(combined_bed_path,0.1,0.1):
        print('IDENTIFIED ATAC BEDS FOR GROUP1, GROUP2, COMBINED:')
        print(group1_bed_path,group2_bed_path,combined_bed_path)

    else:
        print('MAKING ATAC PEAKS')

        #for group1
        print('making group1 atac collection')
        group1_collection = mergeMacs14(macsEnrichedFolder,group1_atac_list)

        #for group2
        print('making group2 atac collection')
        group2_collection = mergeMacs14(macsEnrichedFolder,group2_atac_list)

        #making beds
        group1_bed = utils.locusCollectionToBed(group1_collection)
        group2_bed = utils.locusCollectionToBed(group2_collection)

        #making combined collection
        group1_loci = group1_collection.getLoci()
        group2_loci = group2_collection.getLoci()

        combined_loci = group1_loci + group2_loci
        combined_collection = utils.LocusCollection(combined_loci)
        stitched_collection = combined_collection.stitchCollection()
        combined_bed= utils.locusCollectionToBed(stitched_collection)

        #writing to disk
        utils.unParseTable(group1_bed,group1_bed_path,'\t')
        utils.unParseTable(group2_bed,group2_bed_path,'\t')
        utils.unParseTable(combined_bed,combined_bed_path,'\t')



    print('\n\n')
    print('#======================================================================')
    print('#===============IV. DETECTING ACTIVE GENES IN THE SYSTEM===============')
    print('#======================================================================')
    print('\n\n')

    # assumes active genes defined in 1_map_enhancers.py
    # k27ac macs enriched peak present at the promoter (+/- 1kb tss) AND by expression
    # fpkm > 10 to define an active gene

    activity_path = '%s%s_KERATINOCYTE_ACTIVE.txt' % (geneListFolder, genome)

    if utils.checkOutput(activity_path,0.1,0.1):
        print('IDENTIFIED KERATINOCYTE ACTIVE GENE LIST AT %s' % (activity_path))
    else:
        print('ERROR: CANNOT FIND KERATINOCYTE ACTIVE GENE LIST AT %s' % (activity_path))
        sys.exit()



    print('\n\n')
    print('#======================================================================')
    print('#========================V. RUNNING CIRCUITRY==========================')
    print('#======================================================================')
    print('\n\n')

    #running circuitry on the consensus system
    #creates a sbatch bash script
    crc_folder = '%scrc_atac/' % (projectFolder)

    # #for young
    # analysis_name = 'keratinocyte_young'
    # enhancer_path = '%smeta_rose/young_h3k27ac/young_h3k27ac_SuperEnhancers_ENHANCER_TO_GENE.txt' % (projectFolder)
    # subpeak_path = '%sbeds/HG19_young_atac_-0_+0.bed' % (projectFolder)
    # activity_path = '%sgeneListFolder/HG19_KERATINOCYTE_ACTIVE.txt' % (projectFolder)

    # call_crc(analysis_name,enhancer_path,subpeak_path,activity_path,crc_folder)


    # #for old
    # analysis_name = 'keratinocyte_old'
    # enhancer_path = '%smeta_rose/old_h3k27ac/old_h3k27ac_SuperEnhancers_ENHANCER_TO_GENE.txt' % (projectFolder)
    # subpeak_path = '%sbeds/HG19_old_atac_-0_+0.bed' % (projectFolder)
    # activity_path = '%sgeneListFolder/HG19_KERATINOCYTE_ACTIVE.txt' % (projectFolder)

    # call_crc(analysis_name,enhancer_path,subpeak_path,activity_path,crc_folder)

    # #for combined se
    # analysis_name = 'keratinocyte_combined'
    # enhancer_path = '%sclustering/keratinocyte_se_clustering/HG19_keratinocyte_se_clusterTable_ENHANCER_TO_GENE.txt' % (projectFolder)
    # subpeak_path = '%sbeds/HG19_combined_atac_-0_+0.bed' % (projectFolder)
    # activity_path = '%sgeneListFolder/HG19_KERATINOCYTE_ACTIVE.txt' % (projectFolder)

    # call_crc(analysis_name,enhancer_path,subpeak_path,activity_path,crc_folder)


    # #for combined all enhancers with clustering
    # analysis_name = 'keratinocyte_combined_all'
    # enhancer_path = '%sclustering/keratinocyte_all_clustering/HG19_keratinocyte_all_clusterTable_ENHANCER_TO_GENE.txt' % (projectFolder)
    # subpeak_path = '%sbeds/HG19_combined_atac_-0_+0.bed' % (projectFolder)
    # activity_path = '%sgeneListFolder/HG19_KERATINOCYTE_ACTIVE.txt' % (projectFolder)

    # call_crc(analysis_name,enhancer_path,subpeak_path,activity_path,crc_folder)


    #for combined  enhancers
    # analysis_name = 'keratinocyte_combined'
    enhancers_analysis = config.get_analysis_name('map_enhancers_analysis')
    enhancer_path = '%smeta_rose/%s/%s_AllEnhancers_ENHANCER_TO_GENE.txt' % (projectFolder, enhancers_analysis, enhancers_analysis)
    subpeak_path = '%sbeds/%s_combined_atac_-0_+0.bed' % (projectFolder, genome)
    activity_path = '%sgeneListFolder/%s_KERATINOCYTE_ACTIVE.txt' % (projectFolder, genome)

    call_crc(analysis_name,enhancer_path,subpeak_path,activity_path,crc_folder)
    call_crc_script = 'bash %scrc_atac/%s_crc.sh' % (projectFolder,analysis_name)
    proc = Popen(call_crc_script, shell=True)

    # wait for finishing crc
    proc.wait()

    # if crc script returned 1 (failed), then exit with status 1
    if proc.returncode:
        print 'running %s failed' % (call_crc_script)
        sys.exit(1)

    print('\n\n')
    print('#======================================================================')
    print('#========================VI. DELTA OUT BY EDGE=========================')
    print('#======================================================================')
    print('\n\n')


    #calculating delta out degree by brd4 change at edges. only take edges in the top 50%
    #at least 1 dataset
    # crc_folder = '%scrc_atac/%s' % (projectFolder,analysis_name)
    # analysis_name = '%s' (analysis_name)
    # tf_edge_brd4_delta(crc_folder,chip_dataFile,analysis_name,group1_brd4_list,group2_brd4_list)


    crc_folder = '%scrc_atac/%s' % (projectFolder, analysis_name)
    # analysis_name = '%s' % (analysis_name)
    tf_edge_delta_out(crc_folder,chip_dataFile,analysis_name,group1_brd4_list,group2_brd4_list)

    print('\n\n')
    print('#======================================================================')
    print('#========================VII. DELTA IN BY BRD4=========================')
    print('#======================================================================')
    print('\n\n')


    crc_folder = '%scrc_atac/%s' % (projectFolder, analysis_name)
    # analysis_name = '%s' % (analysis_name)
    tf_edge_delta_in(crc_folder,chip_dataFile,analysis_name,group1_brd4_list,group2_brd4_list)


    print('\n\n')
    print('#======================================================================')
    print('#========================VIII. CONNECTIVITY OF STRING TFS==============')
    print('#======================================================================')
    print('\n\n')


    #get the string TFs
    edge_signal_path = '%scrc_atac/%s/%s_EDGE_TABLE_signal_filtered.txt' % (projectFolder, analysis_name, analysis_name)
    string_clustering_path = config.string_clustering_file
    degree_table_path = '%scrc_atac/%s/%s_DEGREE_TABLE_STRING_TF_signal_filtered.txt' % (projectFolder, analysis_name, analysis_name)
    degree_table_path = makeDegreeTable(edge_signal_path,string_clustering_path,normalize=True,output=degree_table_path)


    print('\n\n')
    print('#======================================================================')
    print('#========================IX. TF COR HEATMAPS===========================')
    print('#======================================================================')
    print('\n\n')

    crc_folder = '%scrc_atac/%s/' % (projectFolder, analysis_name)
    # analysis_name = '%s' % (analysis_name)
    degree_path = '%s%s_DEGREE_TABLE.txt' % (crc_folder,analysis_name)

    degree_table = utils.parseTable(degree_path,'\t')


    max_total_degree = max([int(line[-1]) for line in degree_table[1:]])

    print(max_total_degree)

    cutoff= 0.75
    top_tfs = [line[0] for line in degree_table[1:] if float(line[-1])/max_total_degree > cutoff]
    top_tfs.sort()
    #print(top_tfs)
    print('identified %s top tfs with total degree > %s' % (len(top_tfs),cutoff))

    if not os.path.exists('%smotif_beds/tables' % (crc_folder)):
        os.mkdir('%smotif_beds/tables' % (crc_folder))


    top_tf_table = [[tf] for tf in top_tfs]
    top_tf_path ='%smotif_beds/tables/top_tf_list.txt' % (crc_folder)
    utils.unParseTable(top_tf_table,top_tf_path,'\t')
    print(top_tf_path)

    #now do an additional filtering with string
    string_interaction_path = config.string_interaction_file   #'/grail/string/string_interactions.tsv'
    string_table = utils.parseTable(string_interaction_path,'\t')
    string_tfs = utils.uniquify([line[0] for line in string_table[1:]] + [line[1] for line in string_table[1:]])
    print('identified %s string tfs with interactions' % (len(string_tfs)))

    string_tf_table = [[tf] for tf in string_tfs]
    string_tf_path ='%smotif_beds/tables/string_tf_list.txt' % (crc_folder)
    utils.unParseTable(string_tf_table,string_tf_path,'\t')
    print(string_tf_path)


    ##for top tfs
    #pipeline_dfci.plotCRCCorrMaps('top_tfs','/storage/cylin/home/cl6/projects/NIBR_YvsO_cyl/crc_atac/keratinocyte_combined_all/motif_beds/',tf_list_path=top_tf_path,window=50)

    ##for string tfs
    #pipeline_dfci.plotCRCCorrMaps('string_tfs','/storage/cylin/home/cl6/projects/NIBR_YvsO_cyl/crc_atac/keratinocyte_combined_all/motif_beds/',tf_list_path=string_tf_path,window=50)



    #now make a combined in out scatter
    out_path = '%s%s_EDGE_DELTA_OUT.txt' % (crc_folder,analysis_name)
    in_path = '%s%s_TF_DELTA_IN.txt' % (crc_folder,analysis_name)

    out_table = utils.parseTable(out_path,'\t')
    in_table = utils.parseTable(in_path,'\t')


    #for top tfs
    tf_list = [x[0] for x in in_table[1:] if [y[0] for y in out_table[1:]].count(x[0]) >0]
    tf_list = [tf for tf in tf_list if top_tfs.count(tf) >0]
    print(tf_list)
    print(len(tf_list))

    in_out_table = [['TF_NAME','IN_DELTA','OUT_DELTA']]
    for tf_name in tf_list:
        in_row = [line[0] for line in in_table].index(tf_name)
        out_row = [line[0] for line in out_table].index(tf_name)

        in_delta = in_table[in_row][3]
        out_delta = out_table[out_row][2]

        in_out_table.append([tf_name,in_delta,out_delta])

    in_out_path = '%s%s_IN_OUT_DELTA_TOP_TF.txt' % (crc_folder,analysis_name)
    print(in_out_path)
    utils.unParseTable(in_out_table,in_out_path,'\t')


    # for string tfs
    tf_list = [x[0] for x in in_table[1:] if [y[0] for y in out_table[1:]].count(x[0]) >0]
    tf_list = [tf for tf in tf_list if string_tfs.count(tf) >0]
    print(tf_list)
    print(len(tf_list))

    in_out_table = [['TF_NAME','IN_DELTA','OUT_DELTA','OUT_SEM']]
    for tf_name in tf_list:
        in_row = [line[0] for line in in_table].index(tf_name)
        out_row = [line[0] for line in out_table].index(tf_name)

        in_delta = in_table[in_row][3]
        out_delta = out_table[out_row][2]
        out_sem = out_table[out_row][5]
        in_out_table.append([tf_name,in_delta,out_delta,out_sem])

    in_out_path = '%s%s_IN_OUT_DELTA_STRING_TF.txt' % (crc_folder,analysis_name)
    print(in_out_path)
    utils.unParseTable(in_out_table,in_out_path,'\t')


    print('\n\n')
    print('#======================================================================')
    print('#=========================X. TF TARGET GENES===========================')
    print('#======================================================================')
    print('\n\n')



    tf_name = 'IRF2'
    filtered_signal_path = '%scrc_atac/%s/%s_EDGE_TABLE_signal_filtered.txt' % (projectFolder, analysis_name, analysis_name)
    get_tf_target_genes(tf_name,filtered_signal_path,cut_off = -0.5,output='')



#==========================================================================
#==================================THE END=================================
#==========================================================================


if __name__=="__main__":
    main()
