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

pipeline_dir = config.pipeline_dir

sys.path.append(whereAmI)
sys.path.append(pipeline_dir)

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
genePlotFolder = config.gene_plot_folder

# mask Files
maskFile = config.mask_file

# genome firectory
genomeDirectory = config.genome_dir

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

y_k27ac_list = config.get_global_name('young_h3k27ac_list')
o_k27ac_list = config.get_global_name('old_h3k27ac_list')
k27ac_list = y_k27ac_list + o_k27ac_list

y_brd4_list = config.get_global_name('young_brd4_list')
o_brd4_list = config.get_global_name('old_brd4_list')

y_atac_list = config.get_global_name('young_atac_list')
o_atac_list = config.get_global_name('old_atac_list')

chip_list_no_input = config.get_global_name('chip_list_no_input')


#for some reason needed to re-run macs for that brd4
#pipeline_dfci.callMacsSlurm(chip_dataFile,macsFolder,namesList = ['O_BRD4_1'],overwrite=False,pvalue='1e-9',useBackground =True,launch=False)

#==========================================================================
#===========================MAIN METHOD====================================
#==========================================================================


def main():


    print('main analysis for MYCN project')

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
    young_bed_path = '%sHG19_young_atac_-0_+0.bed' % (bedFolder)
    old_bed_path = '%sHG19_old_atac_-0_+0.bed' % (bedFolder)

    #combined
    combined_bed_path = '%sHG19_combined_atac_-0_+0.bed' % (bedFolder)

    #look for these to already exist
    if utils.checkOutput(young_bed_path,0.1,0.1) and utils.checkOutput(old_bed_path,0.1,0.1) and utils.checkOutput(combined_bed_path,0.1,0.1):
        print('IDENTIFIED ATAC BEDS FOR YOUNG, OLD, COMBINED:')
        print(young_bed_path,old_bed_path,combined_bed_path)

    else:
        print('MAKING ATAC PEAKS')

        #for young
        print('making young atac collection')
        young_collection = mergeMacs14(macsEnrichedFolder,y_atac_list)

        #for old
        print('making old atac collection')
        old_collection = mergeMacs14(macsEnrichedFolder,o_atac_list)

        #making beds
        young_bed = utils.locusCollectionToBed(young_collection)
        old_bed = utils.locusCollectionToBed(old_collection)

        #making combined collection
        young_loci = young_collection.getLoci()
        old_loci = old_collection.getLoci()

        combined_loci = young_loci + old_loci
        combined_collection = utils.LocusCollection(combined_loci)
        stitched_collection = combined_collection.stitchCollection()
        combined_bed= utils.locusCollectionToBed(stitched_collection)

        #writing to disk
        utils.unParseTable(young_bed,young_bed_path,'\t')
        utils.unParseTable(old_bed,old_bed_path,'\t')
        utils.unParseTable(combined_bed,combined_bed_path,'\t')



    print('\n\n')
    print('#======================================================================')
    print('#===============IV. DETECTING ACTIVE GENES IN THE SYSTEM===============')
    print('#======================================================================')
    print('\n\n')

    # assumes active genes defined in 1_map_enhancers.py
    # k27ac macs enriched peak present at the promoter (+/- 1kb tss) AND by expression
    # fpkm > 10 to define an active gene

    activity_path = '%sHG19_KERATINOCYTE_ACTIVE.txt' % (geneListFolder)

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
    analysis_name = 'keratinocyte_combined'
    enhancer_path = '%smeta_rose/combined_h3k27ac/combined_h3k27ac_AllEnhancers_ENHANCER_TO_GENE.txt' % (projectFolder)
    subpeak_path = '%sbeds/HG19_combined_atac_-0_+0.bed' % (projectFolder)
    activity_path = '%sgeneListFolder/HG19_KERATINOCYTE_ACTIVE.txt' % (projectFolder)

    call_crc(analysis_name,enhancer_path,subpeak_path,activity_path,crc_folder)
    #os.system('bash %scrc_atac/keratinocyte_combined_crc.sh' % (projectFolder))
    call_crc_script = 'bash %scrc_atac/keratinocyte_combined_crc.sh' % (projectFolder)
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
    # crc_folder = '%scrc_atac/keratinocyte_combined' % (projectFolder)
    # analysis_name = 'keratinocyte_combined'
    # tf_edge_brd4_delta(crc_folder,chip_dataFile,analysis_name,y_brd4_list,o_brd4_list)


    crc_folder = '%scrc_atac/keratinocyte_combined' % (projectFolder)
    analysis_name = 'keratinocyte_combined'
    tf_edge_brd4_delta_out(crc_folder,chip_dataFile,analysis_name,y_brd4_list,o_brd4_list)

    print('\n\n')
    print('#======================================================================')
    print('#========================VII. DELTA IN BY BRD4=========================')
    print('#======================================================================')
    print('\n\n')


    crc_folder = '%scrc_atac/keratinocyte_combined' % (projectFolder)
    analysis_name = 'keratinocyte_combined'
    tf_edge_brd4_delta_in(crc_folder,chip_dataFile,analysis_name,y_brd4_list,o_brd4_list)


    print('\n\n')
    print('#======================================================================')
    print('#========================VIII. CONNECTIVITY OF STRING TFS==============')
    print('#======================================================================')
    print('\n\n')


    #using the signal filtered edge table to identify connectivity between string TFs

    def makeDegreeTable(edge_signal_path,string_clustering_path,normalize=True,output=''):

        '''
        takes in the edge signal path, a list of TFs, if normalize then normalize to max
        possible degree (bounded 0 to 1 where 1 is all TFs)
        also annotate the cluster for each
        '''
        #get the string TFs

        string_table = utils.parseTable(string_clustering_path,'\t')
        string_tfs = utils.uniquify([line[4] for line in string_table[1:]])

        tf_cluster_dict = {}
        for line in string_table[1:]:
            tf_cluster_dict[line[4]] = [line[1],line[2]]
        print('identified %s string tfs with interactions' % (len(string_tfs)))

        #now load in the edge table
        edge_signal_table = utils.parseTable(edge_signal_path,'\t')

        edge_dict = defaultdict(list)
        all_tfs = []
        for line in edge_signal_table[1:]:
            [source,target] = line[0].split('_')
            all_tfs.append(source)

        all_tfs = utils.uniquify(all_tfs)
        for line in edge_signal_table[1:]:
            [source,target] = line[0].split('_')
            if all_tfs.count(source) > 0:
                if all_tfs.count(target) > 0:
                    edge_dict[source].append(target)



        non_string_tfs = [tf_name for tf_name in all_tfs if string_tfs.count(tf_name) ==0]

        #counting only unique edges
        for tf in all_tfs:
            edge_dict[tf] = utils.uniquify(edge_dict[tf])


        all_out_edges = []
        for tf in all_tfs:
            all_out_edges+= edge_dict[tf]


        degree_table = [['TF','CLUSTER','COLOR','IN_DEGREE','OUT_DEGREE','CLUSTER_IN_DEGREE','OUTGROUP_IN_DEGREE','CLUSTER_OUT_DEGREE','OUTGROUP_OUT_DEGREE']]
        for tf in string_tfs + non_string_tfs:

            out_degree = len(edge_dict[tf])
            in_degree = all_out_edges.count(tf)

            #now figure out which TFs are in the same cluster
            if tf_cluster_dict.keys().count(tf) > 0:
                cluster_tfs = [tf_i for tf_i in tf_cluster_dict.keys() if tf_cluster_dict[tf_i] == tf_cluster_dict[tf]]
                outgroup_tfs = [tf_j for tf_j in all_tfs if cluster_tfs.count(tf_j) == 0]
                cluster_out_degree = len([tf_k for tf_k in edge_dict[tf] if cluster_tfs.count(tf_k) >0])
                outgroup_out_degree = len([tf_l for tf_l in edge_dict[tf] if outgroup_tfs.count(tf_l) >0])

                all_cluster_out_edges = []
                for tf_m in cluster_tfs:
                    all_cluster_out_edges+= edge_dict[tf_m]
                all_outgroup_out_edges = []
                for tf_n in outgroup_tfs:
                    all_outgroup_out_edges+= edge_dict[tf_n]

                cluster_in_degree = all_cluster_out_edges.count(tf)
                outgroup_in_degree = all_outgroup_out_edges.count(tf)


            else:

                cluster_out_degree= 0
                outgroup_out_degree = 0

                cluster_in_degree= 0
                outgroup_in_degree = 0

            if normalize:
                out_degree = round(float(out_degree)/len(all_tfs),2)
                in_degree = round(float(in_degree)/len(all_tfs),2)

                if string_tfs.count(tf) > 0:
                    cluster_out_degree = round(float(cluster_out_degree)/len(cluster_tfs),2)
                    outgroup_out_degree = round(float(outgroup_out_degree)/len(outgroup_tfs),2)

                    cluster_in_degree = round(float(cluster_in_degree)/len(cluster_tfs),2)
                    outgroup_in_degree = round(float(outgroup_in_degree)/len(outgroup_tfs),2)

            if string_tfs.count(tf) > 0:
                new_line = [tf] + tf_cluster_dict[tf] + [in_degree,out_degree,cluster_in_degree,outgroup_in_degree,cluster_out_degree,outgroup_out_degree]
            else:
                new_line = [tf] + ['7','grey'] + [in_degree,out_degree,cluster_in_degree,outgroup_in_degree,cluster_out_degree,outgroup_out_degree]
            degree_table.append(new_line)
        if len(output) >0:
            utils.unParseTable(degree_table,output,'\t')
            return(output)
        else:
            return(degree_table)


    #get the string TFs
    edge_signal_path = '%scrc_atac/keratinocyte_combined/keratinocyte_combined_EDGE_TABLE_signal_filtered.txt' % (projectFolder)
    string_clustering_path = config.string_clustering_path # '/grail/string/string_MCL_clusters.tsv'
    degree_table_path = '%scrc_atac/keratinocyte_combined/keratinocyte_combined_DEGREE_TABLE_STRING_TF_signal_filtered.txt' % (projectFolder)
    degree_table_path = makeDegreeTable(edge_signal_path,string_clustering_path,normalize=True,output=degree_table_path)


    print('\n\n')
    print('#======================================================================')
    print('#========================IX. TF COR HEATMAPS===========================')
    print('#======================================================================')
    print('\n\n')

    crc_folder = '%scrc_atac/keratinocyte_combined/' % (projectFolder)
    analysis_name = 'keratinocyte_combined'
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
    string_interaction_path = config.string_interaction_path   #'/grail/string/string_interactions.tsv'
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


    #use the filtered signal table from tf_edge_brd4_delta_out function

    def get_tf_target_genes(tf_name,filtered_signal_path,cut_off = 0.5,output=''):

        '''
        looks for target genes with at least one edge satisfying the cutoff criteria
        if positive cutoff given, assume greater than, if negative assume less than
        '''

        if output == '':
            output = filtered_signal_path.replace('EDGE_TABLE_signal_filtered.txt','%s_TARGETS_filtered_%s.txt' % (tf_name,cut_off))

        target_table = [['GENE','EDGE','Y_vs_O_BRD4_LOG2']]
        filtered_signal_table = utils.parseTable(filtered_signal_path,'\t')

        for line in filtered_signal_table[1:]:
            source_tf = line[0].split('_')[0]
            target_tf = line[0].split('_')[-1]

            if source_tf == tf_name:
                delta_edge = float(line[-1])
                if cut_off > 0 and delta_edge > cut_off:
                    target_table.append([target_tf,line[1],line[-1]])
                elif cut_off < 0 and delta_edge < cut_off:
                    target_table.append([target_tf,line[1],line[-1]])

        utils.unParseTable(target_table,output,'\t')
        print('For %s identified %s target genes with cutoff of %s' % (tf_name,str(len(target_table)-1), cut_off))
        print('Writing output to %s' % (output))
        return output


    tf_name = 'IRF2'
    filtered_signal_path = '%scrc_atac/keratinocyte_combined/keratinocyte_combined_EDGE_TABLE_signal_filtered.txt' % (projectFolder)
    get_tf_target_genes(tf_name,filtered_signal_path,cut_off = -0.5,output='')



#==========================================================================
#===================SPECIFIC FUNCTIONS FOR ANALYSIS========================
#==========================================================================


#specific functions written for this analysis


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



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~DEFINING KERATINOCYTE H3K27AC ENHANCER LANDSCAPE~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def define_enhancer_landscape(projectFolder,pipeline_dir,analysis_name,chip_dataFile,names_list= []):

    '''
    defines the enhancer baseline using H3K27ac chips
    enhancers defined using auto optimized stitching of nearby regions
    w/ a 2.5kb tss exclusion
    uses the meta rose code and writes out a .sh file for reproducibility
    '''

    #For H3K27AC
    #with TSS exclusion and auto stitching

    dataDict = pipeline_dfci.loadDataTable(chip_dataFile)

    if len(names_list) == 0:
        names_list = [name for name in dataDict.keys() if name.count('H3K27AC') == 1]

    bamFileList = [dataDict[name]['bam'] for name in names_list]
    bamString = string.join(bamFileList,',')

    controlBams = [dataDict[name]['background'] for name in names_list]
    controlFileList = [dataDict[name]['bam'] for name in controlBams]
    controlBamString = string.join(controlFileList,',')

    bedFileList = [macsEnrichedFolder + dataDict[name]['enrichedMacs'] for name in names_list]
    bedString = string.join(bedFileList,',')

    roseFolder = '%smeta_rose/' % (projectFolder)
    roseFolder = utils.formatFolder(roseFolder,True)

    outputFolder = '%s%s/' % (roseFolder,analysis_name)
    bashFileName = '%s%s_meta_rose.sh' % (roseFolder,analysis_name)

    bashFile = open(bashFileName,'w')
    bashFile.write('#!/usr/bin/bash\n\n')
    bashFile.write('#SBATCH --mem=16000\n\n') #set memory for slurm
    bashFile.write('cd %s\n' % (pipeline_dir))

    metaRoseCmd = 'python %sROSE2_META.py -g hg19 -i %s -r %s -c %s -o %s -n %s -t 0 --mask %s' % (pipeline_dir,bedString,bamString,controlBamString,outputFolder,analysis_name,maskFile)

    bashFile.write(metaRoseCmd + '\n')
    bashFile.close()


    #the 4KB parameter is
    region_map_path = '%s%s/%s_SuperEnhancers.table.txt' % (roseFolder,analysis_name,analysis_name)
    return bashFileName,region_map_path


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~MAKING SUBPEAKS FROM MACS14 FILES~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def mergeMacs14(macsEnrichedFolder,names_list):

    '''
    merges a bunch of peaks across dataset names
    spits out as a locus collection
    '''
    macsEnrichedFolder = utils.formatFolder(macsEnrichedFolder,False)

    bed_loci = []
    for name in names_list:
        macs_peak_path = '%s%s_peaks.bed' % (macsEnrichedFolder,name)
        macs_peaks = utils.parseTable(macs_peak_path,'\t')
        peaks_loci = []
        for line in macs_peaks:
            peaks_loci.append(utils.Locus(line[0],line[1],line[2],'.',line[3]))

        bed_loci += peaks_loci

    bed_collection = utils.LocusCollection(bed_loci)
    stitched_collection = bed_collection.stitchCollection()

    return stitched_collection


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~DEFINING ACTIVE GENE LISTS~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


def makeActiveList(mapped_path,exp_path,annotFile,output=''):

    '''
    makes a list of active genes
    '''


    exp_table =utils.parseTable(exp_path,'\t') #BC1 = young BC=old

    exp_dict= defaultdict(float)
    for line in exp_table[1:]:
        exp_vector = [float(x) for x in line[1:]]
        exp_dict[line[0]] = max(exp_vector)

    startDict = utils.makeStartDict(annotFile)

    #now make the active gene list first by chromatin
    mapped_table = utils.parseTable(mapped_path,'\t')
    active_table = []
    for i in range(1,len(mapped_table)):
        line = mapped_table[i]
        #check old

        if max([int(line[i]) for i in [2,3]]) and max([int(line[i]) for i in [4,5,6]]) == 1:
            old_active = True
        else:
            old_active = False

        if max([int(line[i]) for i in [7,8]]) and max([int(line[i]) for i in [9,10,11]]) == 1:
            young_active = True
        else:
            young_active = False

        if old_active or young_active:
            chromatin_active = True

        gene_name = startDict[line[1]]['name']

        if exp_dict[gene_name] > 10.0:
            active_table.append([line[1],gene_name])

    if len(output) == 0:
        return active_table
    else:
        utils.unParseTable(active_table,output,'\t')
        return output


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~CALLING CRC ANALYSIS~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


def call_crc(analysis_name,enhancer_path,subpeak_path,activity_path,crc_folder=''):

    '''
    runs crc
    '''


    if len(crc_folder) == 0:

        crc_folder = utils.formatFolder('%scrc' % (projectFolder),True)
    else:
        crc_folder = utils.formatFolder(crc_folder,True)


    output_folder = utils.formatFolder('%s%s' % (crc_folder,analysis_name),True)

    crc_path = config.crc_path

    crc_cmd = '%s -e %s -g HG19 -o %s -n %s -s %s -a %s -c %s' % (crc_path ,enhancer_path,output_folder,analysis_name,subpeak_path,activity_path, '/grail/genomes/Homo_sapiens/UCSC/hg19/Sequence/Chromosomes/')
    crc_bash_path = '%s%s_crc.sh' % (crc_folder,analysis_name)

    crc_bash = open(crc_bash_path,'w')
    crc_bash.write('#!/usr/bin/bash\n\n')
    crc_bash.write('#SBATCH --mem=16000\n\n')
    crc_bash.write(crc_cmd +'\n')
    crc_bash.close()

    print('wrote crc command for %s to %s' % (analysis_name,crc_bash_path))
    return crc_bash_path


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~PLOTTING METAS ACROSS NB~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


def plot_keratinocyte_genes(plot_prefix,chip_dataFile,figure_gff_path,bed_string):

    '''
    plots all varieties and iterations of tracks for shep on data
    '''


    #first establish the plot folder
    plotFolder = utils.formatFolder('%sKERATINOCYTE_ALL/' % (genePlotFolder),True)
    #plot_prefix = 'HG19_KERATINOCYTE_GENES'

    #we also have to set the extension properly between datasets

    #go by data file
    dataDict = pipeline_dfci.loadDataTable(chip_dataFile)
    names_list = dataDict.keys()
    #initial check for consistency of read lengths
    # for name in names_list:
    #     bam = utils.Bam(dataDict[name]['bam'])
    #     read_length = bam.getReadLengths()[0]
    #     bam_extension = 200-read_length
    #     print('For dataset %s in %s using an extension of %s' % (name,chip_dataFile,bam_extension))

    bam = utils.Bam(dataDict[names_list[0]]['bam'])
    read_length = bam.getReadLengths()[0]
    bam_extension = 200-read_length
    print('For datasets in %s using an extension of %s' % (chip_dataFile,bam_extension))

    #first do individuals
    for plot_group in ['BRD4','H3K27AC']:
        plotList = [name for name in dataDict.keys() if name.upper().count(plot_group) > 0 and name.upper().count('INPUT') == 0]
        plotName = '%s_INDIVIDUAL_%s' % (plot_prefix,plot_group)
        print(plotName)
        pipeline_dfci.callBatchPlot(chip_dataFile,figure_gff_path,plotName,plotFolder,plotList,uniform=True,bed =bed_string,plotType= 'MULTIPLE',extension=bam_extension,multiPage = False,debug=False,nameString = '',rpm=True,rxGenome = '')

    #now as metas
    plotList = y_k27ac_list + o_k27ac_list
    groupList = ['YOUNG_H3K27AC']*3 + ['OLD_H3K27AC']*3
    groupString = ','.join(groupList)

    plotName = '%s_H3K27AC_META_RELATIVE' % (plot_prefix)
    pipeline_dfci.callBatchPlot(chip_dataFile,figure_gff_path,plotName,plotFolder,plotList,uniform=False,bed =bed_string,plotType= 'MERGE',extension=bam_extension,multiPage = False,debug=False,nameString = groupString,rpm=True,rxGenome = '')


    plotName = '%s_H3K27AC_META_UNIFORM' % (plot_prefix)
    pipeline_dfci.callBatchPlot(chip_dataFile,figure_gff_path,plotName,plotFolder,plotList,uniform=True,bed =bed_string,plotType= 'MERGE',extension=bam_extension,multiPage = False,debug=False,nameString = groupString,rpm=True,rxGenome = '')


    #now as metas for brd4
    plotList = y_brd4_list + o_brd4_list
    groupList = ['YOUNG_BRD4']*2 + ['OLD_BRD4']*2
    groupString = ','.join(groupList)

    plotName = '%s_BRD4_META_RELATIVE' % (plot_prefix)
    pipeline_dfci.callBatchPlot(chip_dataFile,figure_gff_path,plotName,plotFolder,plotList,uniform=False,bed =bed_string,plotType= 'MERGE',extension=bam_extension,multiPage = False,debug=False,nameString = groupString,rpm=True,rxGenome = '')


    plotName = '%s_BRD4_META_UNIFORM' % (plot_prefix)
    pipeline_dfci.callBatchPlot(chip_dataFile,figure_gff_path,plotName,plotFolder,plotList,uniform=True,bed =bed_string,plotType= 'MERGE',extension=bam_extension,multiPage = False,debug=False,nameString = groupString,rpm=True,rxGenome = '')

    #now metas for all
    plotList =  y_k27ac_list + o_k27ac_list +y_brd4_list + o_brd4_list
    groupList = ['YOUNG_H3K27AC']*3 + ['OLD_H3K27AC']*3  +['YOUNG_BRD4']*2 + ['OLD_BRD4']*2
    groupString = ','.join(groupList)

    plotName = '%s_ALL_META_RELATIVE' % (plot_prefix)
    pipeline_dfci.callBatchPlot(chip_dataFile,figure_gff_path,plotName,plotFolder,plotList,uniform=False,bed =bed_string,plotType= 'MERGE',extension=bam_extension,multiPage = False,debug=False,nameString = groupString,rpm=True,rxGenome = '')


    plotName = '%s_ALL_META_UNIFORM' % (plot_prefix)
    pipeline_dfci.callBatchPlot(chip_dataFile,figure_gff_path,plotName,plotFolder,plotList,uniform=True,bed =bed_string,plotType= 'MERGE',extension=bam_extension,multiPage = False,debug=False,nameString = groupString,rpm=True,rxGenome = '')





def plot_keratinocyte_atac(plot_prefix,atac_dataFile,figure_gff_path,bed_string):

    '''
    plots all varieties and iterations of tracks for shep on data
    '''


    #first establish the plot folder
    plotFolder = utils.formatFolder('%sKERATINOCYTE_ATAC/' % (genePlotFolder),True)
    #plot_prefix = 'HG19_KERATINOCYTE_GENES'

    #we also have to set the extension properly between datasets

    #go by data file
    dataDict = pipeline_dfci.loadDataTable(atac_dataFile)
    names_list = dataDict.keys()
    #initial check for consistency of read lengths
    # for name in names_list:
    #     bam = utils.Bam(dataDict[name]['bam'])
    #     read_length = bam.getReadLengths()[0]
    #     bam_extension = 200-read_length
    #     print('For dataset %s in %s using an extension of %s' % (name,atac_dataFile,bam_extension))

    bam = utils.Bam(dataDict[names_list[0]]['bam'])
    read_length = bam.getReadLengths()[0]
    bam_extension = 200-read_length
    print('For datasets in %s using an extension of %s' % (atac_dataFile,bam_extension))

    #first do individuals
    for plot_group in ['Y','O']:
        plotList = [name for name in dataDict.keys() if name.upper().count(plot_group) > 0 and name.upper().count('INPUT') == 0]
        plotName = '%s_INDIVIDUAL_%s' % (plot_prefix,plot_group)
        print(plotName)
        pipeline_dfci.callBatchPlot(atac_dataFile,figure_gff_path,plotName,plotFolder,plotList,uniform=True,bed =bed_string,plotType= 'MULTIPLE',extension=bam_extension,multiPage = False,debug=False,nameString = '',rpm=True,rxGenome = '')

    y_list = [name for name in dataDict.keys() if name.upper().count('Y') >0]
    o_list = [name for name in dataDict.keys() if name.upper().count('O') >0]
    #now as metas
    plotList = y_list + o_list
    groupList = ['YOUNG']*len(y_list) + ['OLD']*len(o_list)
    groupString = ','.join(groupList)

    plotName = '%s_META_RELATIVE' % (plot_prefix)
    pipeline_dfci.callBatchPlot(atac_dataFile,figure_gff_path,plotName,plotFolder,plotList,uniform=False,bed =bed_string,plotType= 'MERGE',extension=bam_extension,multiPage = False,debug=False,nameString = groupString,rpm=True,rxGenome = '')


    plotName = '%s_META_UNIFORM' % (plot_prefix)
    pipeline_dfci.callBatchPlot(atac_dataFile,figure_gff_path,plotName,plotFolder,plotList,uniform=True,bed =bed_string,plotType= 'MERGE',extension=bam_extension,multiPage = False,debug=False,nameString = groupString,rpm=True,rxGenome = '')



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~RUNNING ENHANCER PROMOTER ANALYSIS~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def wrap_enhancer_promoter(dataFile,input_path,activity_path,analysis_name,names_list,useBackground=True,top=5000):

    '''
    runs enhancer promoter on everybody with the conserved regions and union of active genes
    '''

    #hard coded paths
    tads_path ='%shESC_domains_hg19.bed' %(bedFolder)

    #setting the output folder
    ep_folder = utils.formatFolder('%senhancerPromoter/' % (projectFolder),True)



    dataDict = pipeline_dfci.loadDataTable(dataFile)

    bams_list = [dataDict[name]['bam'] for name in names_list]
    bams_string = ' '.join(bams_list)

    background_names = [dataDict[name]['background'] for name in names_list]
    background_list = [dataDict[background_name]['bam'] for background_name in background_names]
    background_string = ' '.join(background_list)


    ep_bash_path = '%s%s_enhancer_promoter.sh' % (ep_folder,analysis_name)
    ep_bash = open(ep_bash_path,'w')

    ep_bash.write('#!/usr/bin/bash\n\n\n')
    ep_bash.write('#SBATCH --mem=16000\n\n')

    ep_bash.write('#enhancer promoter analysis for %s\n\n' % (analysis_name))

    ep_bash.write('cd %s%s/\n\n' % (ep_folder,analysis_name))
    if useBackground:

        python_cmd = 'python %senhancerPromoter.py -b %s -c %s -g %s -i %s -o %s -a %s --name %s --tads %s --top %s\n\n' % (pipeline_dir,bams_string,background_string,genome.upper(),input_path,ep_folder,activity_path,analysis_name,tads_path,top)

        ep_bash.write(python_cmd)
    else:

        python_cmd = 'python %senhancerPromoter.py -b %s -g %s -i %s -o %s -a %s --name %s --tads %s --top %s\n\n' % (pipeline_dir,bams_string,genome.upper(),input_path,ep_folder,activity_path,analysis_name,tads_path,top)

        ep_bash.write(python_cmd)

    ep_bash.close()

    #generate an anticipated output to check for completeness
    output_path = '%s%s/%s_TOP_%s_ORDERED.txt' % (ep_folder,analysis_name,analysis_name,top)

    return(ep_bash_path,output_path)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~CALCULATING CHANGES IN BRD4 OUT DEGREE BY TF~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


def tf_motif_brd4_delta(crc_folder,chip_dataFile,analysis_name,y_brd4_list,o_brd4_list):

    '''
    calculates changes in brd4 out degree at each predicted motif occurrence this is by subpeaks
    '''
    crc_folder = utils.formatFolder(crc_folder,False)
    subpeak_bed_path = '%s%s_all_subpeak.bed' % (crc_folder,analysis_name)
    #direct the output to the crc folder
    signal_path = '%s%s_all_subpeak_signal.txt' % (crc_folder,analysis_name)

    #set final output
    output_path = '%s%s_brd4_subpeak_delta_out.txt' % (crc_folder,analysis_name)


    all_brd4_list = y_brd4_list + o_brd4_list
    pipeline_dfci.map_regions(chip_dataFile,[subpeak_bed_path],mappedFolder,signalFolder,all_brd4_list,True,signal_path)

    #now bring in the signal table as a dictionary using the locus line as the id
    print('making log2 y vs o signal dict at subpeaks')
    signal_table = utils.parseTable(signal_path,'\t')
    signal_dict = defaultdict(float)

    #figure out columns for young and old
    o_columns = [signal_table[0].index(name) for name in o_brd4_list]
    y_columns = [signal_table[0].index(name) for name in y_brd4_list]
    for line in signal_table[1:]:
        o_signal = numpy.mean([float(line[col]) for col in o_columns])
        y_signal = numpy.mean([float(line[col]) for col in y_columns])
        delta = numpy.log2(y_signal/o_signal)
        signal_dict[line[1]] = delta

    #loading subpeak collection
    subpeak_table = utils.parseTable(subpeak_bed_path,'\t')
    subpeak_loci = []

    for line in subpeak_table:
        line_locus = utils.Locus(line[0],int(line[1]),int(line[2]),'.')
        subpeak_loci.append(line_locus)


    #next create a dictionary of all tfs in the system
    tf_table_path = '%s%s_GENE_TF_TABLE.txt' % (crc_folder,analysis_name)

    tf_table = utils.parseTable(tf_table_path,'\t')
    tf_list = utils.uniquify([line[0] for line in tf_table[1:]])
    tf_list.sort()

    print('finding motifs for tfs')
    print(tf_list)


    out_degree_delta_dict = defaultdict(list)
    for tf_name in tf_list:
        print(tf_name)
        motif_bed_path = '%smotif_beds/%s_motifs.bed' % (crc_folder,tf_name)
        motif_bed = utils.parseTable(motif_bed_path,'\t')
        motif_loci = []
        for line in motif_bed[1:]:
            motif_locus = utils.Locus(line[0],line[1],line[2],'.')
            motif_loci.append(motif_locus)

        motif_collection = utils.LocusCollection(motif_loci)

        #now iterate through the subpeak loci
        for subpeak_locus in subpeak_loci:
            if motif_collection.getOverlap(subpeak_locus,'both'):
                out_degree_delta_dict[tf_name].append(signal_dict[subpeak_locus.__str__()])

    out_degree_table = [['TF_NAME','DELTA_MEAN','DELTA_MEDIAN','DELTA_STD']]
    for tf_name in tf_list:
        delta_line = [tf_name] + [numpy.mean(out_degree_delta_dict[tf_name]),numpy.median(out_degree_delta_dict[tf_name]),numpy.std(out_degree_delta_dict[tf_name])]
        out_degree_table.append(delta_line)


    utils.unParseTable(out_degree_table,output_path,'\t')
    print(output_path)
    return(output_path)







#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~CALCULATING CHANGES IN BRD4 OUT DEGREE BY TF CLUSTERS~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


def tf_group_brd4_delta(crc_folder,chip_dataFile,analysis_name,tf_list,tf_group_name,y_brd4_list,o_brd4_list):

    '''
    calculates changes in brd4 out degree at each predicted motif occurrence
    '''
    crc_folder = utils.formatFolder(crc_folder,False)
    subpeak_bed_path = '%s%s_all_subpeak.bed' % (crc_folder,analysis_name)
    #direct the output to the crc folder
    signal_path = '%s%s_all_subpeak_signal.txt' % (crc_folder,analysis_name)

    #set final output
    output_path = '%s%s_all_subpeak_signal_%s.txt' % (crc_folder,analysis_name,tf_group_name)

    if len(tf_list) == 0:
        #next create a dictionary of all tfs in the system
        tf_table_path = '/storage/cylin/home/cl6/projects/NIBR_YvsO_cyl/crc/keratinocyte_combined/keratinocyte_combined_GENE_TF_TABLE.txt'

        tf_table = utils.parseTable(tf_table_path,'\t')
        tf_list = utils.uniquify([line[0] for line in tf_table[1:]])
        tf_list.sort()

    #now bring in the signal table as a dictionary using the locus line as the id
    print('making log2 y vs o signal dict at subpeaks')
    signal_table = utils.parseTable(signal_path,'\t')
    signal_dict = defaultdict(float)

    #figure out columns for young and old
    o_columns = [signal_table[0].index(name) for name in o_brd4_list]
    y_columns = [signal_table[0].index(name) for name in y_brd4_list]
    for line in signal_table[1:]:
        o_signal = numpy.mean([float(line[col]) for col in o_columns])
        y_signal = numpy.mean([float(line[col]) for col in y_columns])
        delta = numpy.log2(y_signal/o_signal)
        signal_dict[line[1]] = delta


    #loading motifs
    print('making motif dictionaries')
    motif_dict = {}
    for tf_name in tf_list:
        print(tf_name)
        motif_bed_path = '%smotif_beds/%s_motifs.bed' % (crc_folder,tf_name)
        motif_bed = utils.parseTable(motif_bed_path,'\t')
        motif_loci = []
        for line in motif_bed[1:]:
            motif_locus = utils.Locus(line[0],line[1],line[2],'.')
            motif_loci.append(motif_locus)

        motif_collection = utils.LocusCollection(motif_loci)
        motif_dict[tf_name] = motif_collection


    #loading subpeak collection
    subpeak_table = utils.parseTable(subpeak_bed_path,'\t')
    #set the output table
    out_signal_table = [signal_table[0] + ['LENGTH'] + tf_list]
    for i in range(len(subpeak_table)):
        line = subpeak_table[i]
        line_locus = utils.Locus(line[0],int(line[1]),int(line[2]),'.')

        motif_line = [line_locus.len()]
        for tf_name in tf_list:
            motif_line.append(len(motif_dict[tf_name].getOverlap(line_locus)))
        out_signal_table.append([line[0]] + signal_table[(i+1)][1:] + motif_line)

    utils.unParseTable(out_signal_table,output_path,'\t')
    print(output_path)
    return(output_path)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~CALCULATING CHANGES IN BRD4 OUT DEGREE BY TF EDGES~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


def tf_edge_brd4_delta_out(crc_folder,chip_dataFile,analysis_name,y_brd4_list,o_brd4_list,output=''):

    '''
    calculates changes in brd4 out degree at each predicted motif occurrence this is by subpeaks
    '''
    crc_folder = utils.formatFolder(crc_folder,False)
    edge_path = '%s%s_EDGE_TABLE.txt' % (crc_folder,analysis_name)

    #make a gff of the edge table
    edge_table = utils.parseTable(edge_path,'\t')
    edge_gff = []
    for line in edge_table[1:]:
        gff_line = [line[2],'%s_%s' % (line[0],line[1]),'',line[3],line[4],'','.','','%s_%s' % (line[0],line[1])]
        edge_gff.append(gff_line)

    edge_gff_path = '%s%s_EDGE_TABLE.gff' % (crc_folder,analysis_name)
    utils.unParseTable(edge_gff,edge_gff_path,'\t')


    #direct the output to the crc folder
    signal_path = '%s%s_EDGE_TABLE_signal.txt' % (crc_folder,analysis_name)



    all_brd4_list = y_brd4_list + o_brd4_list
    if utils.checkOutput(signal_path,0,0) == False:
        signal_table_list = pipeline_dfci.map_regions(chip_dataFile,[edge_gff_path],mappedFolder,signalFolder,all_brd4_list,True,signal_path,extendReadsTo=100)
        print(signal_table_list)
    else:
        print('Found previous signal table at %s' % (signal_path))



    #now bring in the signal table as a dictionary using the locus line as the id
    print('making log2 y vs o signal table at edges')
    signal_table = utils.parseTable(signal_path,'\t')
    signal_dict = defaultdict(float)

    #figure out columns for young and old
    o_columns = [signal_table[0].index(name) for name in o_brd4_list]
    y_columns = [signal_table[0].index(name) for name in y_brd4_list]
    o_signal_vector = []
    y_signal_vector = []
    for line in signal_table[1:]:
        o_signal = numpy.mean([float(line[col]) for col in o_columns])
        y_signal = numpy.mean([float(line[col]) for col in y_columns])

        o_signal_vector.append(o_signal)
        y_signal_vector.append(y_signal)

    o_median = numpy.median(o_signal_vector)
    y_median = numpy.median(y_signal_vector)

    print('old median BRD4 signal')
    print(o_median)
    print('young median BRD4 signal')
    print(y_median)

    #now that we have the median, we can take edges where at least 1 edge is above the median
    #and both are above zero and generate a new table w/ the fold change

    signal_filtered_path = string.replace(signal_path,'.txt','_filtered.txt')
    if utils.checkOutput(signal_filtered_path,0,0):
        print('Found filtered signal table for edges at %s' % (signal_filtered_path))
        signal_table_filtered = utils.parseTable(signal_filtered_path,'\t')
    else:

        signal_table_filtered = [signal_table[0]+['O_BRD4_MEAN','Y_BRD4_MEAN','Y_vs_O_LOG2']]
        for line in signal_table[1:]:
            o_signal = numpy.mean([float(line[col]) for col in o_columns])
            y_signal = numpy.mean([float(line[col]) for col in y_columns])

            if (o_signal > o_median or y_signal > y_median) and min(o_signal,y_signal) >0:
                delta = numpy.log2(y_signal/o_signal)
                new_line = line + [o_signal,y_signal,delta]
                signal_table_filtered.append(new_line)

        utils.unParseTable(signal_table_filtered,signal_filtered_path,'\t')

    #now get a list of all TFs in the system
    tf_list = utils.uniquify([line[0].split('_')[0] for line in signal_table_filtered[1:]])
    tf_list.sort()
    print(tf_list)

    out_degree_table = [['TF_NAME','EDGE_COUNT','DELTA_MEAN','DELTA_MEDIAN','DELTA_STD','DELTA_SEM']]

    for tf_name in tf_list:
        print(tf_name)
        edge_vector = [float(line[-1]) for line in signal_table_filtered[1:] if line[0].split('_')[0] == tf_name]

        edge_count = len(edge_vector)
        delta_mean = round(numpy.mean(edge_vector),4)
        delta_median = round(numpy.median(edge_vector),4)
        delta_std = round(numpy.std(edge_vector),4)
        delta_sem = round(stats.sem(edge_vector),4)
        tf_out_line = [tf_name,edge_count,delta_mean,delta_median,delta_std,delta_sem]
        out_degree_table.append(tf_out_line)

    if output == '':
        #set final output
        output_path = '%s%s_EDGE_DELTA_OUT.txt' % (crc_folder,analysis_name)

    else:
        output_path = output

    utils.unParseTable(out_degree_table,output_path,'\t')
    print(output_path)
    return(output_path)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~CALCULATING CHANGES IN BRD4 IN DEGREE AT TFS~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


def tf_edge_brd4_delta_in(crc_folder,chip_dataFile,analysis_name,y_brd4_list,o_brd4_list,output=''):

    '''
    calculates changes in BRD4 at in degree edges
    '''


    crc_folder = utils.formatFolder(crc_folder,False)
    enhancer_tf_path = '%s%s_ENHANCER_TF_TABLE.txt' % (crc_folder,analysis_name)

    #make a gff of the edge table
    enhancer_tf_table = utils.parseTable(enhancer_tf_path,'\t')
    enhancer_tf_gff = []
    for line in enhancer_tf_table[1:]:
        gff_line = [line[1],line[0],'',line[2],line[3],'','.','',line[0]]
        enhancer_tf_gff.append(gff_line)

    enhancer_tf_gff_path = '%s%s_ENHANCER_TF_TABLE.gff' % (crc_folder,analysis_name)
    utils.unParseTable(enhancer_tf_gff,enhancer_tf_gff_path,'\t')


    #direct the output to the crc folder
    signal_path = '%s%s_ENHANCER_TF_TABLE_signal.txt' % (crc_folder,analysis_name)


    all_brd4_list = y_brd4_list + o_brd4_list
    if utils.checkOutput(signal_path,0,0) == False:
        signal_table_list = pipeline_dfci.map_regions(chip_dataFile,[enhancer_tf_gff_path],mappedFolder,signalFolder,all_brd4_list,True,signal_path,extendReadsTo=100)
        print(signal_table_list)
    else:
        print('Found previous signal table at %s' % (signal_path))



    #now bring in the signal table as a dictionary using the locus line as the id
    print('making an enhancer signal dict')
    signal_table = utils.parseTable(signal_path,'\t')
    y_signal_dict = defaultdict(float)
    o_signal_dict = defaultdict(float)

    #signal here is calculated as AUC

    #figure out columns for young and old
    o_columns = [signal_table[0].index(name) for name in o_brd4_list]
    y_columns = [signal_table[0].index(name) for name in y_brd4_list]

    for line in signal_table[1:]:
        region_coords = [int(x) for x in line[1].split(':')[-1].split('-')]
        region_length = region_coords[1] - region_coords[0]
        o_signal = region_length*numpy.mean([float(line[col]) for col in o_columns])
        y_signal = region_length*numpy.mean([float(line[col]) for col in y_columns])
        y_signal_dict[line[0]] = y_signal
        o_signal_dict[line[0]] = o_signal

    #now grab the gene table
    gene_tf_path = '%s%s_GENE_TF_TABLE.txt' % (crc_folder,analysis_name)
    gene_tf_table = utils.parseTable(gene_tf_path,'\t')
    y_tf_dict = defaultdict(float)
    o_tf_dict = defaultdict(float)

    for line in gene_tf_table[1:]:
        y_tf_dict[line[0]] += y_signal_dict[line[-1]]
        o_tf_dict[line[0]] += o_signal_dict[line[-1]]

    tf_list = utils.uniquify([line[0] for line in gene_tf_table[1:]])

    tf_list.sort()
    print(tf_list)

    in_degree_table = [['TF_NAME','Y_BRD4_IN','O_BRD4_IN','LOG2_Y_vs_O_BRD4']]

    for tf_name in tf_list:
        y_signal = round(y_tf_dict[tf_name],4)
        o_signal = round(o_tf_dict[tf_name],4)
        delta = round(numpy.log2(y_signal/o_signal),4)
        new_line = [tf_name,y_signal,o_signal,delta]
        in_degree_table.append(new_line)


    if output == '':
        #set final output
        output_path = '%s%s_TF_DELTA_IN.txt' % (crc_folder,analysis_name)

    else:
        output_path = output

    utils.unParseTable(in_degree_table,output_path,'\t')
    print(output_path)
    return(output_path)



#==========================================================================
#==================================THE END=================================
#==========================================================================


if __name__=="__main__":
    main()
