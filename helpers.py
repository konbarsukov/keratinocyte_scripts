

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

import pipeline_dfci
import utils
import string
import numpy
import os
from scipy import stats
import re
from collections import defaultdict
#==========================================================================
#============================PARAMETERS====================================
#==========================================================================



projectName = config.project_name
genome = config.genome_name
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


py27_path = config.py27_path #'/opt/bin/python' #'/storage/cylin/anaconda3/envs/py27_anaconda/bin/python'




#==========================================================================
#===================SPECIFIC FUNCTIONS FOR ANALYSIS========================
#==========================================================================


# script 1

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~DEFINING ACTIVE GENE LISTS~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


def makeActiveList(mapped_path,exp_path,annotFile,output=''):

    '''
    makes a list of active genes
    10fpkm expression cutoff
    bound by k27ac and brd4 in at least 1 dataset
    '''


    exp_table =utils.parseTable(exp_path,'\t') #BC1 = group1 BC=group2

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
        #check group2
        if max([int(line[i]) for i in [2,3]]) and max([int(line[i]) for i in [4,5,6]]) == 1:
            group2_active = True
        else:
            group2_active = False

        if max([int(line[i]) for i in [7,8]]) and max([int(line[i]) for i in [9,10,11]]) == 1:
            group1_active = True
        else:
            group1_active = False

        if group2_active or group1_active:
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


# script 3

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


# script 4

#specific functions written for this analysis



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


#use the filtered signal table from tf_edge_delta_out function
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
#~~~~~~~~~~~~~~~~~~~~~~~~~~CALLING CRC ANALYSIS~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


def call_crc(analysis_name,enhancer_path,subpeak_path,activity_path,crc_folder):

    '''
    runs crc
    '''

    crc_folder = utils.formatFolder(crc_folder,True)

    output_folder = utils.formatFolder('%s%s' % (crc_folder,analysis_name),True)

    crc_path = config.crc_path

    crc_cmd = '%s -e %s -g HG19 -o %s -n %s -s %s -a %s -c %s' % (crc_path ,enhancer_path,output_folder,analysis_name,subpeak_path,activity_path, config.genome_folder)
    crc_bash_path = '%s%s_crc.sh' % (crc_folder,analysis_name)

    crc_bash = open(crc_bash_path,'w')
    crc_bash.write('#!/usr/bin/bash\n\n')
    crc_bash.write('#SBATCH --mem=16000\n\n')
    crc_bash.write(crc_cmd +'\n')
    crc_bash.close()

    print('wrote crc command for %s to %s' % (analysis_name,crc_bash_path))
    return crc_bash_path



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~CALCULATING CHANGES IN BRD4 OUT DEGREE BY TF EDGES~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


def tf_edge_delta_out(crc_folder,chip_dataFile,analysis_name,group1_list,group2_list,output=''):

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

    #get a list of all chip datasets
    all_chip_list = group1_list + group2_list
    if utils.checkOutput(signal_path,0,0) == False:
        signal_table_list = pipeline_dfci.map_regions(chip_dataFile,[edge_gff_path],mappedFolder,signalFolder,all_chip_list,True,signal_path,extendReadsTo=100)
        print(signal_table_list)
    else:
        print('Found previous signal table at %s' % (signal_path))

    #now bring in the signal table as a dictionary using the locus line as the id
    print('making log2 group1 over group2 signal table at edges')
    signal_table = utils.parseTable(signal_path,'\t')
    signal_dict = defaultdict(float)

    #figure out columns for group1 and group2
    group2_columns = [signal_table[0].index(name) for name in group2_list]
    group1_columns = [signal_table[0].index(name) for name in group1_list]
    group2_signal_vector = []
    group1_signal_vector = []
    for line in signal_table[1:]:
        group2_signal = numpy.mean([float(line[col]) for col in group2_columns])
        group1_signal = numpy.mean([float(line[col]) for col in group1_columns])

        group2_signal_vector.append(group2_signal)
        group1_signal_vector.append(group1_signal)

    group2_median = numpy.median(group2_signal_vector)
    group1_median = numpy.median(group1_signal_vector)

    print('group2 median signal (rpm/bp)')
    print(group2_median)
    print('group1 median signal (rpm/bp)')
    print(group1_median)

    #now that we have the median, we can take edges where at least 1 edge is above the median
    #and both are above zero and generate a new table w/ the fold change

    signal_filtered_path = string.replace(signal_path,'.txt','_filtered.txt')
    if utils.checkOutput(signal_filtered_path,0,0):
        print('Found filtered signal table for edges at %s' % (signal_filtered_path))
        signal_table_filtered = utils.parseTable(signal_filtered_path,'\t')
    else:

        signal_table_filtered = [signal_table[0]+['GROUP2_MEAN','GROUP1_MEAN','LOG2_GROUP1_OVER_GROUP2']]
        for line in signal_table[1:]:
            group2_signal = numpy.mean([float(line[col]) for col in group2_columns])
            group1_signal = numpy.mean([float(line[col]) for col in group1_columns])

            if (group2_signal > group2_median or group1_signal > group1_median) and min(group2_signal,group1_signal) >0:
                delta = numpy.log2(group1_signal/group2_signal)
                new_line = line + [group2_signal,group1_signal,delta]
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



# def tf_edge_brd4_delta_out(crc_folder,chip_dataFile,analysis_name,y_brd4_list,o_brd4_list,output=''):
#
#     '''
#     calculates changes in brd4 out degree at each predicted motif occurrence this is by subpeaks
#     '''
#     crc_folder = utils.formatFolder(crc_folder,False)
#     edge_path = '%s%s_EDGE_TABLE.txt' % (crc_folder,analysis_name)
#
#     #make a gff of the edge table
#     edge_table = utils.parseTable(edge_path,'\t')
#     edge_gff = []
#     for line in edge_table[1:]:
#         gff_line = [line[2],'%s_%s' % (line[0],line[1]),'',line[3],line[4],'','.','','%s_%s' % (line[0],line[1])]
#         edge_gff.append(gff_line)
#
#     edge_gff_path = '%s%s_EDGE_TABLE.gff' % (crc_folder,analysis_name)
#     utils.unParseTable(edge_gff,edge_gff_path,'\t')
#
#
#     #direct the output to the crc folder
#     signal_path = '%s%s_EDGE_TABLE_signal.txt' % (crc_folder,analysis_name)
#
#
#
#     all_brd4_list = y_brd4_list + o_brd4_list
#     if utils.checkOutput(signal_path,0,0) == False:
#         signal_table_list = pipeline_dfci.map_regions(chip_dataFile,[edge_gff_path],mappedFolder,signalFolder,all_brd4_list,True,signal_path,extendReadsTo=100)
#         print(signal_table_list)
#     else:
#         print('Found previous signal table at %s' % (signal_path))
#
#
#
#     #now bring in the signal table as a dictionary using the locus line as the id
#     print('making log2 y vs o signal table at edges')
#     signal_table = utils.parseTable(signal_path,'\t')
#     signal_dict = defaultdict(float)
#
#     #figure out columns for young and old
#     o_columns = [signal_table[0].index(name) for name in o_brd4_list]
#     y_columns = [signal_table[0].index(name) for name in y_brd4_list]
#     o_signal_vector = []
#     y_signal_vector = []
#     for line in signal_table[1:]:
#         o_signal = numpy.mean([float(line[col]) for col in o_columns])
#         y_signal = numpy.mean([float(line[col]) for col in y_columns])
#
#         o_signal_vector.append(o_signal)
#         y_signal_vector.append(y_signal)
#
#     o_median = numpy.median(o_signal_vector)
#     y_median = numpy.median(y_signal_vector)
#
#     print('old median BRD4 signal')
#     print(o_median)
#     print('young median BRD4 signal')
#     print(y_median)
#
#     #now that we have the median, we can take edges where at least 1 edge is above the median
#     #and both are above zero and generate a new table w/ the fold change
#
#     signal_filtered_path = string.replace(signal_path,'.txt','_filtered.txt')
#     if utils.checkOutput(signal_filtered_path,0,0):
#         print('Found filtered signal table for edges at %s' % (signal_filtered_path))
#         signal_table_filtered = utils.parseTable(signal_filtered_path,'\t')
#     else:
#
#         signal_table_filtered = [signal_table[0]+['O_BRD4_MEAN','Y_BRD4_MEAN','Y_vs_O_LOG2']]
#         for line in signal_table[1:]:
#             o_signal = numpy.mean([float(line[col]) for col in o_columns])
#             y_signal = numpy.mean([float(line[col]) for col in y_columns])
#
#             if (o_signal > o_median or y_signal > y_median) and min(o_signal,y_signal) >0:
#                 delta = numpy.log2(y_signal/o_signal)
#                 new_line = line + [o_signal,y_signal,delta]
#                 signal_table_filtered.append(new_line)
#
#         utils.unParseTable(signal_table_filtered,signal_filtered_path,'\t')
#
#     #now get a list of all TFs in the system
#     tf_list = utils.uniquify([line[0].split('_')[0] for line in signal_table_filtered[1:]])
#     tf_list.sort()
#     print(tf_list)
#
#     out_degree_table = [['TF_NAME','EDGE_COUNT','DELTA_MEAN','DELTA_MEDIAN','DELTA_STD','DELTA_SEM']]
#
#     for tf_name in tf_list:
#         print(tf_name)
#         edge_vector = [float(line[-1]) for line in signal_table_filtered[1:] if line[0].split('_')[0] == tf_name]
#
#         edge_count = len(edge_vector)
#         delta_mean = round(numpy.mean(edge_vector),4)
#         delta_median = round(numpy.median(edge_vector),4)
#         delta_std = round(numpy.std(edge_vector),4)
#         delta_sem = round(stats.sem(edge_vector),4)
#         tf_out_line = [tf_name,edge_count,delta_mean,delta_median,delta_std,delta_sem]
#         out_degree_table.append(tf_out_line)
#
#     if output == '':
#         #set final output
#         output_path = '%s%s_EDGE_DELTA_OUT.txt' % (crc_folder,analysis_name)
#
#     else:
#         output_path = output
#
#     utils.unParseTable(out_degree_table,output_path,'\t')
#     print(output_path)
#     return(output_path)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~CALCULATING CHANGES IN BRD4 IN DEGREE AT TFS~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


def tf_edge_delta_in(crc_folder,chip_dataFile,analysis_name,group1_list,group2_list,output=''):

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

    all_chip_list = group1_list + group2_list
    if utils.checkOutput(signal_path,0,0) == False:
        signal_table_list = pipeline_dfci.map_regions(chip_dataFile,[enhancer_tf_gff_path],mappedFolder,signalFolder,all_chip_list,True,signal_path,extendReadsTo=100)
        print(signal_table_list)
    else:
        print('Found previous signal table at %s' % (signal_path))


    #now bring in the signal table as a dictionary using the locus line as the id
    print('making an enhancer signal dict')
    signal_table = utils.parseTable(signal_path,'\t')
    group1_signal_dict = defaultdict(float)
    group2_signal_dict = defaultdict(float)

    #signal here is calculated as AUC

    #figure out columns for group2 and group1
    group2_columns = [signal_table[0].index(name) for name in group2_list]
    group1_columns = [signal_table[0].index(name) for name in group1_list]

    for line in signal_table[1:]:
        region_coords = [int(x) for x in line[1].split(':')[-1].split('-')]
        region_length = region_coords[1] - region_coords[0]
        group2_signal = region_length*numpy.mean([float(line[col]) for col in group2_columns])
        group1_signal = region_length*numpy.mean([float(line[col]) for col in group1_columns])
        group1_signal_dict[line[0]] = group1_signal
        group2_signal_dict[line[0]] = group2_signal

    #now grab the gene table
    gene_tf_path = '%s%s_GENE_TF_TABLE.txt' % (crc_folder,analysis_name)
    gene_tf_table = utils.parseTable(gene_tf_path,'\t')
    group1_tf_dict = defaultdict(float)
    group2_tf_dict = defaultdict(float)

    for line in gene_tf_table[1:]:
        group1_tf_dict[line[0]] += group1_signal_dict[line[-1]]
        group2_tf_dict[line[0]] += group2_signal_dict[line[-1]]

    tf_list = utils.uniquify([line[0] for line in gene_tf_table[1:]])

    tf_list.sort()
    print(tf_list)

    in_degree_table = [['TF_NAME','GROUP1_IN','GROUP2_IN','LOG2_GROUP1_vs_GROUP2']]

    for tf_name in tf_list:
        group1_signal = round(group1_tf_dict[tf_name],4)
        group2_signal = round(group2_tf_dict[tf_name],4)
        delta = round(numpy.log2(group1_signal/group2_signal),4)
        new_line = [tf_name,group1_signal,group2_signal,delta]
        in_degree_table.append(new_line)


    if output == '':
        #set final output
        output_path = '%s%s_TF_DELTA_IN.txt' % (crc_folder,analysis_name)

    else:
        output_path = output

    utils.unParseTable(in_degree_table,output_path,'\t')
    print(output_path)
    return(output_path)


# def tf_edge_brd4_delta_in(crc_folder,chip_dataFile,analysis_name,y_brd4_list,o_brd4_list,output=''):
#
#     '''
#     calculates changes in BRD4 at in degree edges
#     '''
#
#
#     crc_folder = utils.formatFolder(crc_folder,False)
#     enhancer_tf_path = '%s%s_ENHANCER_TF_TABLE.txt' % (crc_folder,analysis_name)
#
#     #make a gff of the edge table
#     enhancer_tf_table = utils.parseTable(enhancer_tf_path,'\t')
#     enhancer_tf_gff = []
#     for line in enhancer_tf_table[1:]:
#         gff_line = [line[1],line[0],'',line[2],line[3],'','.','',line[0]]
#         enhancer_tf_gff.append(gff_line)
#
#     enhancer_tf_gff_path = '%s%s_ENHANCER_TF_TABLE.gff' % (crc_folder,analysis_name)
#     utils.unParseTable(enhancer_tf_gff,enhancer_tf_gff_path,'\t')
#
#
#     #direct the output to the crc folder
#     signal_path = '%s%s_ENHANCER_TF_TABLE_signal.txt' % (crc_folder,analysis_name)
#
#
#     all_brd4_list = y_brd4_list + o_brd4_list
#     if utils.checkOutput(signal_path,0,0) == False:
#         signal_table_list = pipeline_dfci.map_regions(chip_dataFile,[enhancer_tf_gff_path],mappedFolder,signalFolder,all_brd4_list,True,signal_path,extendReadsTo=100)
#         print(signal_table_list)
#     else:
#         print('Found previous signal table at %s' % (signal_path))
#
#
#
#     #now bring in the signal table as a dictionary using the locus line as the id
#     print('making an enhancer signal dict')
#     signal_table = utils.parseTable(signal_path,'\t')
#     y_signal_dict = defaultdict(float)
#     o_signal_dict = defaultdict(float)
#
#     #signal here is calculated as AUC
#
#     #figure out columns for young and old
#     o_columns = [signal_table[0].index(name) for name in o_brd4_list]
#     y_columns = [signal_table[0].index(name) for name in y_brd4_list]
#
#     for line in signal_table[1:]:
#         region_coords = [int(x) for x in line[1].split(':')[-1].split('-')]
#         region_length = region_coords[1] - region_coords[0]
#         o_signal = region_length*numpy.mean([float(line[col]) for col in o_columns])
#         y_signal = region_length*numpy.mean([float(line[col]) for col in y_columns])
#         y_signal_dict[line[0]] = y_signal
#         o_signal_dict[line[0]] = o_signal
#
#     #now grab the gene table
#     gene_tf_path = '%s%s_GENE_TF_TABLE.txt' % (crc_folder,analysis_name)
#     gene_tf_table = utils.parseTable(gene_tf_path,'\t')
#     y_tf_dict = defaultdict(float)
#     o_tf_dict = defaultdict(float)
#
#     for line in gene_tf_table[1:]:
#         y_tf_dict[line[0]] += y_signal_dict[line[-1]]
#         o_tf_dict[line[0]] += o_signal_dict[line[-1]]
#
#     tf_list = utils.uniquify([line[0] for line in gene_tf_table[1:]])
#
#     tf_list.sort()
#     print(tf_list)
#
#     in_degree_table = [['TF_NAME','Y_BRD4_IN','O_BRD4_IN','LOG2_Y_vs_O_BRD4']]
#
#     for tf_name in tf_list:
#         y_signal = round(y_tf_dict[tf_name],4)
#         o_signal = round(o_tf_dict[tf_name],4)
#         delta = round(numpy.log2(y_signal/o_signal),4)
#         new_line = [tf_name,y_signal,o_signal,delta]
#         in_degree_table.append(new_line)
#
#
#     if output == '':
#         #set final output
#         output_path = '%s%s_TF_DELTA_IN.txt' % (crc_folder,analysis_name)
#
#     else:
#         output_path = output
#
#     utils.unParseTable(in_degree_table,output_path,'\t')
#     print(output_path)
#     return(output_path)
