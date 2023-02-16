import argparse
import os
import shutil
#  from TSS_Mapping_End_command import Preprocessing
from TSS_Mapping_End_command import Coverage_mapping_singleunit
from TSS_Mapping_End_command import Coverage_mapping_multiunit

# This is the sample command line that is how user has to define the variables, and flags;
# python3 command_line_interface.py -r /Users/nehamittal/Desktop/TSS_Gene_Mapping/Reference_genomes/hg38.ncbiRefSeq.gtf -s /Users/nehamittal/Desktop/TSS_Gene_Mapping/Samples_Data/streck_cfDNA_205_S3__chr21__alignment_summary.csv 
# -o /Users/nehamittal/Desktop/Multiple_Results -c chr21 -f transcript -u 5000 -d 5000 -si 0 -ei 800 -Yax Strand)

#Set up argparse 
parser = argparse.ArgumentParser()
parser.add_argument('-r', '--ref_file',
                        help = 'Reference File',
                        required = True)
parser.add_argument('-s', '--sample_file',
                        help = 'Sample File',
                        required = True)
parser.add_argument('-o', '--output_dir',
                    help = 'Enter folder name for saving results',
                    required = True)
parser.add_argument('-c', '--chr_category',
                    help = 'Enter chromosome number to filter',
                    required = True)
parser.add_argument('-f', '--feature',
                    help = 'Enter the feature to filter(transcript/exon/CDS)',
                    required = True)
parser.add_argument('-u', '--upstream',
                    help = 'Enter the upstream to extend from tss site',
                    required = True)
parser.add_argument('-d', '--downstream',
                    help = 'Enter the downstream to extend from tss site',
                    required = True)
parser.add_argument('-i', '--index',
                    help = 'Enter the index value for the particular gene based on TSS extension file',
                    required = False)
parser.add_argument('-si', '--startindex',
                    help = 'Enter the startindex value for the particular gene based on TSS extension file that user wants in the range',
                    required = False)
parser.add_argument('-ei', '--endindex',
                    help = 'Enter the endindex value for the particular gene based on TSS extension file that user wants in the range',
                    required = False)
parser.add_argument('-Yax', '--Yaxis_category',
                    help = 'Enter the variable name on which user wants to do the classification for Y-axis at heatmap',
                    required = False)
args = parser.parse_args()

# # Raise error if reference file or sample file is missing
# if args.ref_file is None:
#     raise ValueError('Reference File is missing')
# elif args.sample_file is None:
#     raise ValueError('Sample File is missing')


#create directory folder for reports if the directory does not already exists
if os.path.isdir(args.output_dir) == False:
    os.makedirs(args.output_dir)
    # os.makedirs(args.output_dir +'/mapping_results')

else:
    option = input("Output Directory Already Existing.\n Do you wish to overwrite? (y/n) ")
    if option == 'y':
        shutil.rmtree(args.output_dir)
        os.makedirs(args.output_dir)
    elif option=='n':
        # print("ghbnm")
        new_path = input("Enter new path ")
        if os.path.isdir(new_path) == False:
            os.makedirs(new_path)
            args.output_dir = new_path
        else:
            raise ValueError('Output directory already exists. Please delete the current directory or rename')
    else:
        raise ValueError('Output directory already exists. Please delete the current directory or rename')

if (args.startindex is None) or (args.endindex is None) and (args.index is not None):
    print("Start and End index not passed")
    figure_format = input("Enter the type of figure format: png, pdf, jpeg, png and pdf, png and jpeg, pdf and jpeg, all three")
    CMS = Coverage_mapping_singleunit(ref_file = args.ref_file, sample_file = args.sample_file, output_dir = args.output_dir, 
                 chr_category = args.chr_category, feature = args.feature, upstream = args.upstream, 
                 downstream = args.downstream, index = args.index)
    
    TSS_extension = CMS.TSS_sites_extension()
    Coverage_read = CMS.Sample_coverage_calculation()
    decode_list_value = CMS.Coverage_mapping(Coverage_read, TSS_extension)
    CMS.Graph_plotting(decode_list_value, figure_format)
    
elif (args.startindex is not None) and (args.endindex is not None) and (args.index is None):
    print("Start and End index  passed")
    plot_type =  input("Enter type of plot: Heatmap/Multiple univariate/both")
    figure_format = input("Enter the type of figure format: png, pdf, jpeg, png and pdf, png and jpeg, pdf and jpeg, all three")
    CMM = Coverage_mapping_multiunit(ref_file = args.ref_file, sample_file = args.sample_file, output_dir = args.output_dir, 
                    chr_category = args.chr_category, feature = args.feature, upstream = args.upstream, 
                    downstream = args.downstream, startindex = args.startindex, endindex = args.endindex)
    
    TSS_extension = CMM.TSS_sites_extension()
    Coverage_read = CMM.Sample_coverage_calculation()

    multithread_decode_list_final = CMM.Coverage_mapping_multi(Coverage_read, TSS_extension)

    if plot_type == "Heatmap":
        CMM.get_Heatmap_multivariate_mapping(args.Yaxis_category, figure_format)
    elif plot_type == "Multiple univariate":
        CMM.get_multiple_univariate_plots(figure_format)
    elif plot_type == "both":
        CMM.get_multiple_univariate_plots(figure_format)
        CMM.get_Heatmap_multivariate_mapping(args.Yaxis_category, figure_format)
    else:
        raise ValueError("Invalid Plot type Selected")
else:
    pass

