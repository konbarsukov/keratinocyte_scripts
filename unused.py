# script 1

#specific functions written for this analysis



# #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# #~~~~~~~~~~~~~~~~~~~~MAKING SUBPEAKS FROM MACS2 FILES~~~~~~~~~~~~~~~~~~~~~~
# #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# def mergeMacs2(macs2Folder,names_list):
#
#     '''
#     merges a bunch of peaks across dataset names
#     spits out as a locus collection
#     '''
#     macs2Folder = utils.formatFolder(macs2Folder,False)
#
#     bed_loci = []
#     for name in names_list:
#         macs2_peak_path = '%s%s/%s_peaks.narrowPeak' % (macs2Folder,name,name)
#         macs2_peaks = utils.parseTable(macs2_peak_path,'\t')
#         peaks_loci = []
#         for line in macs2_peaks:
#             peaks_loci.append(utils.Locus(line[0],line[1],line[2],'.',line[3]))
#
#         bed_loci += peaks_loci
#
#     bed_collection = utils.LocusCollection(bed_loci)
#     stitched_collection = bed_collection.stitchCollection()
#
#     return stitched_collection

# #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# #~~~~~~~~~~~~~~~~~~~~~~~~PLOTTING METAS ACROSS NB~~~~~~~~~~~~~~~~~~~~~~~~~~
# #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#
# def plot_keratinocyte_genes(plot_prefix,chip_dataFile,figure_gff_path,bed_string):
#
#     '''
#     plots all varieties and iterations of tracks for shep on data
#     '''
#
#
#     #first establish the plot folder
#     plotFolder = utils.formatFolder('%sKERATINOCYTE_ALL/' % (genePlotFolder),True)
#     #plot_prefix = 'HG19_KERATINOCYTE_GENES'
#
#     #we also have to set the extension properly between datasets
#
#     #go by data file
#     dataDict = pipeline_dfci.loadDataTable(chip_dataFile)
#     names_list = dataDict.keys()
#     #initial check for consistency of read lengths
#     # for name in names_list:
#     #     bam = utils.Bam(dataDict[name]['bam'])
#     #     read_length = bam.getReadLengths()[0]
#     #     bam_extension = 200-read_length
#     #     print('For dataset %s in %s using an extension of %s' % (name,chip_dataFile,bam_extension))
#
#     bam = utils.Bam(dataDict[names_list[0]]['bam'])
#     read_length = bam.getReadLengths()[0]
#     bam_extension = 200-read_length
#     print('For datasets in %s using an extension of %s' % (chip_dataFile,bam_extension))
#
#     #first do individuals
#     for plot_group in ['BRD4','H3K27AC']:
#         plotList = [name for name in dataDict.keys() if name.upper().count(plot_group) > 0 and name.upper().count('INPUT') == 0]
#         plotName = '%s_INDIVIDUAL_%s' % (plot_prefix,plot_group)
#         print(plotName)
#         pipeline_dfci.callBatchPlot(chip_dataFile,figure_gff_path,plotName,plotFolder,plotList,uniform=True,bed =bed_string,plotType= 'MULTIPLE',extension=bam_extension,multiPage = False,debug=False,nameString = '',rpm=True,rxGenome = '')
#
#     #now as metas
#     plotList = y_k27ac_list + o_k27ac_list
#     groupList = ['YOUNG_H3K27AC']*3 + ['OLD_H3K27AC']*3
#     groupString = ','.join(groupList)
#
#     plotName = '%s_H3K27AC_META_RELATIVE' % (plot_prefix)
#     pipeline_dfci.callBatchPlot(chip_dataFile,figure_gff_path,plotName,plotFolder,plotList,uniform=False,bed =bed_string,plotType= 'MERGE',extension=bam_extension,multiPage = False,debug=False,nameString = groupString,rpm=True,rxGenome = '')
#
#
#     plotName = '%s_H3K27AC_META_UNIFORM' % (plot_prefix)
#     pipeline_dfci.callBatchPlot(chip_dataFile,figure_gff_path,plotName,plotFolder,plotList,uniform=True,bed =bed_string,plotType= 'MERGE',extension=bam_extension,multiPage = False,debug=False,nameString = groupString,rpm=True,rxGenome = '')
#
#
#     #now as metas for brd4
#     plotList = y_brd4_list + o_brd4_list
#     groupList = ['YOUNG_BRD4']*2 + ['OLD_BRD4']*2
#     groupString = ','.join(groupList)
#
#     plotName = '%s_BRD4_META_RELATIVE' % (plot_prefix)
#     pipeline_dfci.callBatchPlot(chip_dataFile,figure_gff_path,plotName,plotFolder,plotList,uniform=False,bed =bed_string,plotType= 'MERGE',extension=bam_extension,multiPage = False,debug=False,nameString = groupString,rpm=True,rxGenome = '')
#
#
#     plotName = '%s_BRD4_META_UNIFORM' % (plot_prefix)
#     pipeline_dfci.callBatchPlot(chip_dataFile,figure_gff_path,plotName,plotFolder,plotList,uniform=True,bed =bed_string,plotType= 'MERGE',extension=bam_extension,multiPage = False,debug=False,nameString = groupString,rpm=True,rxGenome = '')
#
#     #now metas for all
#     plotList =  y_k27ac_list + o_k27ac_list +y_brd4_list + o_brd4_list
#     groupList = ['YOUNG_H3K27AC']*3 + ['OLD_H3K27AC']*3  +['YOUNG_BRD4']*2 + ['OLD_BRD4']*2
#     groupString = ','.join(groupList)
#
#     plotName = '%s_ALL_META_RELATIVE' % (plot_prefix)
#     pipeline_dfci.callBatchPlot(chip_dataFile,figure_gff_path,plotName,plotFolder,plotList,uniform=False,bed =bed_string,plotType= 'MERGE',extension=bam_extension,multiPage = False,debug=False,nameString = groupString,rpm=True,rxGenome = '')
#
#
#     plotName = '%s_ALL_META_UNIFORM' % (plot_prefix)
#     pipeline_dfci.callBatchPlot(chip_dataFile,figure_gff_path,plotName,plotFolder,plotList,uniform=True,bed =bed_string,plotType= 'MERGE',extension=bam_extension,multiPage = False,debug=False,nameString = groupString,rpm=True,rxGenome = '')





# def plot_keratinocyte_atac(plot_prefix,atac_dataFile,figure_gff_path,bed_string):
#
#     '''
#     plots all varieties and iterations of tracks for shep on data
#     '''
#
#
#     #first establish the plot folder
#     plotFolder = utils.formatFolder('%sKERATINOCYTE_ATAC/' % (genePlotFolder),True)
#     #plot_prefix = 'HG19_KERATINOCYTE_GENES'
#
#     #we also have to set the extension properly between datasets
#
#     #go by data file
#     dataDict = pipeline_dfci.loadDataTable(atac_dataFile)
#     names_list = dataDict.keys()
#     #initial check for consistency of read lengths
#     # for name in names_list:
#     #     bam = utils.Bam(dataDict[name]['bam'])
#     #     read_length = bam.getReadLengths()[0]
#     #     bam_extension = 200-read_length
#     #     print('For dataset %s in %s using an extension of %s' % (name,atac_dataFile,bam_extension))
#
#     bam = utils.Bam(dataDict[names_list[0]]['bam'])
#     read_length = bam.getReadLengths()[0]
#     bam_extension = 200-read_length
#     print('For datasets in %s using an extension of %s' % (atac_dataFile,bam_extension))
#
#     #first do individuals
#     for plot_group in ['Y','O']:
#         plotList = [name for name in dataDict.keys() if name.upper().count(plot_group) > 0 and name.upper().count('INPUT') == 0]
#         plotName = '%s_INDIVIDUAL_%s' % (plot_prefix,plot_group)
#         print(plotName)
#         pipeline_dfci.callBatchPlot(atac_dataFile,figure_gff_path,plotName,plotFolder,plotList,uniform=True,bed =bed_string,plotType= 'MULTIPLE',extension=bam_extension,multiPage = False,debug=False,nameString = '',rpm=True,rxGenome = '')
#
#     y_list = [name for name in dataDict.keys() if name.upper().count('Y') >0]
#     o_list = [name for name in dataDict.keys() if name.upper().count('O') >0]
#     #now as metas
#     plotList = y_list + o_list
#     groupList = ['YOUNG']*len(y_list) + ['OLD']*len(o_list)
#     groupString = ','.join(groupList)
#
#     plotName = '%s_META_RELATIVE' % (plot_prefix)
#     pipeline_dfci.callBatchPlot(atac_dataFile,figure_gff_path,plotName,plotFolder,plotList,uniform=False,bed =bed_string,plotType= 'MERGE',extension=bam_extension,multiPage = False,debug=False,nameString = groupString,rpm=True,rxGenome = '')
#
#
#     plotName = '%s_META_UNIFORM' % (plot_prefix)
#     pipeline_dfci.callBatchPlot(atac_dataFile,figure_gff_path,plotName,plotFolder,plotList,uniform=True,bed =bed_string,plotType= 'MERGE',extension=bam_extension,multiPage = False,debug=False,nameString = groupString,rpm=True,rxGenome = '')



# #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# #~~~~~~~~~~~~~~~~~~~~RUNNING ENHANCER PROMOTER ANALYSIS~~~~~~~~~~~~~~~~~~~~
# #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# def wrap_enhancer_promoter(dataFile,input_path,activity_path,analysis_name,names_list,useBackground=True,top=5000):
#
#     '''
#     runs enhancer promoter on everybody with the conserved regions and union of active genes
#     '''
#
#     #hard coded paths
#     tads_path ='%shESC_domains_hg19.bed' %(bedFolder)
#
#     #setting the output folder
#     ep_folder = utils.formatFolder('%senhancerPromoter/' % (projectFolder),True)
#
#
#
#     dataDict = pipeline_dfci.loadDataTable(dataFile)
#
#     bams_list = [dataDict[name]['bam'] for name in names_list]
#     bams_string = ' '.join(bams_list)
#
#     background_names = [dataDict[name]['background'] for name in names_list]
#     background_list = [dataDict[background_name]['bam'] for background_name in background_names]
#     background_string = ' '.join(background_list)
#
#
#     ep_bash_path = '%s%s_enhancer_promoter.sh' % (ep_folder,analysis_name)
#     ep_bash = open(ep_bash_path,'w')
#
#     ep_bash.write('#!/usr/bin/bash\n\n\n')
#     ep_bash.write('#SBATCH --mem=16000\n\n')
#
#     ep_bash.write('#enhancer promoter analysis for %s\n\n' % (analysis_name))
#
#     ep_bash.write('cd %s%s/\n\n' % (ep_folder,analysis_name))
#     if useBackground:
#
#         python_cmd = 'python %senhancerPromoter.py -b %s -c %s -g %s -i %s -o %s -a %s --name %s --tads %s --top %s\n\n' % (pipeline_dir,bams_string,background_string,genome.upper(),input_path,ep_folder,activity_path,analysis_name,tads_path,top)
#
#         ep_bash.write(python_cmd)
#     else:
#
#         python_cmd = 'python %senhancerPromoter.py -b %s -g %s -i %s -o %s -a %s --name %s --tads %s --top %s\n\n' % (pipeline_dir,bams_string,genome.upper(),input_path,ep_folder,activity_path,analysis_name,tads_path,top)
#
#         ep_bash.write(python_cmd)
#
#     ep_bash.close()
#
#     #generate an anticipated output to check for completeness
#     output_path = '%s%s/%s_TOP_%s_ORDERED.txt' % (ep_folder,analysis_name,analysis_name,top)
#
#     return(ep_bash_path,output_path)



# script 4



# #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# #~~~~~~~~~~~~~~~CALCULATING CHANGES IN BRD4 OUT DEGREE BY TF~~~~~~~~~~~~~~~
# #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#
# def tf_motif_brd4_delta(crc_folder,chip_dataFile,analysis_name,y_brd4_list,o_brd4_list):
#
#     '''
#     calculates changes in brd4 out degree at each predicted motif occurrence this is by subpeaks
#     '''
#     crc_folder = utils.formatFolder(crc_folder,False)
#     subpeak_bed_path = '%s%s_all_subpeak.bed' % (crc_folder,analysis_name)
#     #direct the output to the crc folder
#     signal_path = '%s%s_all_subpeak_signal.txt' % (crc_folder,analysis_name)
#
#     #set final output
#     output_path = '%s%s_brd4_subpeak_delta_out.txt' % (crc_folder,analysis_name)
#
#
#     all_brd4_list = y_brd4_list + o_brd4_list
#     pipeline_dfci.map_regions(chip_dataFile,[subpeak_bed_path],mappedFolder,signalFolder,all_brd4_list,True,signal_path)
#
#     #now bring in the signal table as a dictionary using the locus line as the id
#     print('making log2 y vs o signal dict at subpeaks')
#     signal_table = utils.parseTable(signal_path,'\t')
#     signal_dict = defaultdict(float)
#
#     #figure out columns for young and old
#     o_columns = [signal_table[0].index(name) for name in o_brd4_list]
#     y_columns = [signal_table[0].index(name) for name in y_brd4_list]
#     for line in signal_table[1:]:
#         o_signal = numpy.mean([float(line[col]) for col in o_columns])
#         y_signal = numpy.mean([float(line[col]) for col in y_columns])
#         delta = numpy.log2(y_signal/o_signal)
#         signal_dict[line[1]] = delta
#
#     #loading subpeak collection
#     subpeak_table = utils.parseTable(subpeak_bed_path,'\t')
#     subpeak_loci = []
#
#     for line in subpeak_table:
#         line_locus = utils.Locus(line[0],int(line[1]),int(line[2]),'.')
#         subpeak_loci.append(line_locus)
#
#
#     #next create a dictionary of all tfs in the system
#     tf_table_path = '%s%s_GENE_TF_TABLE.txt' % (crc_folder,analysis_name)
#
#     tf_table = utils.parseTable(tf_table_path,'\t')
#     tf_list = utils.uniquify([line[0] for line in tf_table[1:]])
#     tf_list.sort()
#
#     print('finding motifs for tfs')
#     print(tf_list)
#
#
#     out_degree_delta_dict = defaultdict(list)
#     for tf_name in tf_list:
#         print(tf_name)
#         motif_bed_path = '%smotif_beds/%s_motifs.bed' % (crc_folder,tf_name)
#         motif_bed = utils.parseTable(motif_bed_path,'\t')
#         motif_loci = []
#         for line in motif_bed[1:]:
#             motif_locus = utils.Locus(line[0],line[1],line[2],'.')
#             motif_loci.append(motif_locus)
#
#         motif_collection = utils.LocusCollection(motif_loci)
#
#         #now iterate through the subpeak loci
#         for subpeak_locus in subpeak_loci:
#             if motif_collection.getOverlap(subpeak_locus,'both'):
#                 out_degree_delta_dict[tf_name].append(signal_dict[subpeak_locus.__str__()])
#
#     out_degree_table = [['TF_NAME','DELTA_MEAN','DELTA_MEDIAN','DELTA_STD']]
#     for tf_name in tf_list:
#         delta_line = [tf_name] + [numpy.mean(out_degree_delta_dict[tf_name]),numpy.median(out_degree_delta_dict[tf_name]),numpy.std(out_degree_delta_dict[tf_name])]
#         out_degree_table.append(delta_line)
#
#
#     utils.unParseTable(out_degree_table,output_path,'\t')
#     print(output_path)
#     return(output_path)







# #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# #~~~~~~~~~~~~CALCULATING CHANGES IN BRD4 OUT DEGREE BY TF CLUSTERS~~~~~~~~~
# #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#
# def tf_group_brd4_delta(crc_folder,chip_dataFile,analysis_name,tf_list,tf_group_name,y_brd4_list,o_brd4_list):
#
#     '''
#     calculates changes in brd4 out degree at each predicted motif occurrence
#     '''
#     crc_folder = utils.formatFolder(crc_folder,False)
#     subpeak_bed_path = '%s%s_all_subpeak.bed' % (crc_folder,analysis_name)
#     #direct the output to the crc folder
#     signal_path = '%s%s_all_subpeak_signal.txt' % (crc_folder,analysis_name)
#
#     #set final output
#     output_path = '%s%s_all_subpeak_signal_%s.txt' % (crc_folder,analysis_name,tf_group_name)
#
#     if len(tf_list) == 0:
#         #next create a dictionary of all tfs in the system
#         tf_table_path = '/storage/cylin/home/cl6/projects/NIBR_YvsO_cyl/crc/keratinocyte_combined/keratinocyte_combined_GENE_TF_TABLE.txt'
#
#         tf_table = utils.parseTable(tf_table_path,'\t')
#         tf_list = utils.uniquify([line[0] for line in tf_table[1:]])
#         tf_list.sort()
#
#     #now bring in the signal table as a dictionary using the locus line as the id
#     print('making log2 y vs o signal dict at subpeaks')
#     signal_table = utils.parseTable(signal_path,'\t')
#     signal_dict = defaultdict(float)
#
#     #figure out columns for young and old
#     o_columns = [signal_table[0].index(name) for name in o_brd4_list]
#     y_columns = [signal_table[0].index(name) for name in y_brd4_list]
#     for line in signal_table[1:]:
#         o_signal = numpy.mean([float(line[col]) for col in o_columns])
#         y_signal = numpy.mean([float(line[col]) for col in y_columns])
#         delta = numpy.log2(y_signal/o_signal)
#         signal_dict[line[1]] = delta
#
#
#     #loading motifs
#     print('making motif dictionaries')
#     motif_dict = {}
#     for tf_name in tf_list:
#         print(tf_name)
#         motif_bed_path = '%smotif_beds/%s_motifs.bed' % (crc_folder,tf_name)
#         motif_bed = utils.parseTable(motif_bed_path,'\t')
#         motif_loci = []
#         for line in motif_bed[1:]:
#             motif_locus = utils.Locus(line[0],line[1],line[2],'.')
#             motif_loci.append(motif_locus)
#
#         motif_collection = utils.LocusCollection(motif_loci)
#         motif_dict[tf_name] = motif_collection
#
#
#     #loading subpeak collection
#     subpeak_table = utils.parseTable(subpeak_bed_path,'\t')
#     #set the output table
#     out_signal_table = [signal_table[0] + ['LENGTH'] + tf_list]
#     for i in range(len(subpeak_table)):
#         line = subpeak_table[i]
#         line_locus = utils.Locus(line[0],int(line[1]),int(line[2]),'.')
#
#         motif_line = [line_locus.len()]
#         for tf_name in tf_list:
#             motif_line.append(len(motif_dict[tf_name].getOverlap(line_locus)))
#         out_signal_table.append([line[0]] + signal_table[(i+1)][1:] + motif_line)
#
#     utils.unParseTable(out_signal_table,output_path,'\t')
#     print(output_path)
#     return(output_path)
