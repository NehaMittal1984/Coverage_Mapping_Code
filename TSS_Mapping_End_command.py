import pyranges as pr
import pandas as pd
import numpy as np
import threading
import matplotlib.pyplot as plt1 #These are the python module dependencies;
import matplotlib.pyplot as plt2
import seaborn as sn

class Preprocessing:
    def __init__(self, ref_file, sample_file, output_dir, chr_category, feature, upstream, downstream):
        self.ref_file = ref_file
        self.sample_file = sample_file
        self.output_dir = output_dir
        self.chr_category = chr_category
        self.feature = feature
        self.upstream = int(upstream)
        self.downstream = int(downstream)
        # self.index = int(index)
        # self.Yax_category = Yax_category

        # TSS_extension = self.TSS_sites_extension()
        # Read = self.Sample_file()
        # self.Coverage_calculation(Read)
        # decode_list_value = self.Coverage_calculation(Read)
        # print(decode_list_value)
        # self.Plotting(decode_list_value)


   ### ***Reading Reference file, Filtration, Identification of the TSS Sites, and Extended upstream/downstream regions*** ### 
    def TSS_sites_extension(self):
        genome = pr.read_gtf(self.ref_file) #Reading the reference genome file;
        genome1 = (genome.Chromosome == self.chr_category) & (genome.Feature == self.feature) #Filtered all the transcript rows on chromosome 21 from the reference file but still having all the data columns;
        genome1 = genome[genome1].drop_duplicate_positions(strand = True, keep = 'first') #Removing the start, and end positions duplicates;
        print(genome1) #printing the genome1 variable;
        genome1 = genome1.features.tss() #to find the tss site using pyranges;
        print(genome1) #printing the tss sites;
        genome1.df.to_csv(self.output_dir+"/TSS_Sites.csv", header=True, index=True) #storing the TSS sites in to the csv format;
        TSS_sites = pd.read_csv(self.output_dir+"/TSS_Sites.csv", usecols = ['Chromosome','Start', 'End','Strand',"gene_id"]) #Reading the CSV file in pandas data frame;
        TSS_sites_chr21 = TSS_sites.drop_duplicates(subset='gene_id', keep="first", ignore_index=True) #filtered the duplicate gene IDs';
        TSS_sites_chr21.to_csv(self.output_dir+"/TSS_Sites_Remove_Gene_Duplicates.csv", header=True, index=True) #Writing in to a new CSV file;
        self.TSS1 = pd.read_csv(self.output_dir+"/TSS_Sites_Remove_Gene_Duplicates.csv", usecols = [1,2,3,4,5]) #Reading the sample file;
        # Category = TSS1["self.Yax_category"] #filter for category definition over the y-axis of heat map from the reference file on which you want classification;
        TSS2 = pr.PyRanges(self.TSS1) #Converting the Start, and End Columns in to the range using PyRanges function;
        TSS2 = TSS2.extend({"3": self.upstream, "5": self.downstream}) #Calculating the upstream and downstream 5kb region;
        # print(TSS2) #printing the TSS sites with extennsion;
        # print(len(TSS2)) #print the length of the TSS2 file with extension;
        TSS2.to_csv(self.output_dir+"/TSS_Sites_Extension.csv")
        return(TSS2)
    
    def Get_TSS1(self):
        return self.TSS1

    ### ****Sample File Reading, Filtration, and Converting in to the Pyranges start and end position*** ###
    def Sample_coverage_calculation(self):
        Read1 = pd.read_csv(self.sample_file) #Reading the sample file;
        Read1.columns = ['Chromosome', 'Start', 'End', 'Mapping Quality', 'Fragment Length'] #Adding the name in each columns;
        Read1 = pr.PyRanges(Read1) #Converting the Start, and End Columns in to the range using PyRanges function;
        Read1 = Read1.drop_duplicate_positions(strand = None, keep = 'first') #Removing the duplicated range of start and end position;
        print(Read1) #printing the data files after cleaning;
        Read1.df.to_csv(self.output_dir+"/Sample_File.csv", header=True, index=True)

        # return(Read1)

### ***Calculating the run length encodings using read data file and Pyranges RLE function (Coverage calculation)*** ###
        coverage_read1 = (Read1.to_rle()) #Calculating the coverage using read length encoding function of pyranges
        print(coverage_read1) #printing the read encoding length of the sample file;
        return(coverage_read1)

### *** This class will do mapping of the sample rle over a particular index from reference file and plot in the graph *** ###
    
class Coverage_mapping_singleunit(Preprocessing):
    def __init__(self, ref_file, sample_file, output_dir, chr_category, feature, upstream, downstream, index):
        Preprocessing.__init__(self, ref_file, sample_file, output_dir, chr_category, feature, upstream, downstream)
        self.index = int(index)

    #Mapping the sample read coverage over the specific TSS site (specific Gene)#
       
    def Coverage_mapping(self, Coverage_read, TSS_extension):

    # for i in range(len(TSS_extension)):
        site_run = Coverage_read[self.chr_category][TSS_extension.Start.values[self.index]:TSS_extension.End.values[self.index]] #mapping the coverage of read over gene(s);
        run_val = list(site_run.values) #list of value (which value) that has to repeat;
    # print(run_val)
        run_len = list(site_run.runs) #list of value that how many time above value has to repeat;
    # print(run_len)

    # **** To Do: Code Needed, For Loop for each index value or for each gene to plot gene by gene*** #

    # *******Decoding the run length encoding********** #
        decode_list = []

        for x, y in zip(run_val, run_len): #decoding the values based on length     
            z = [x]*y
            decode_list.append(z) #the values are in array here;
        decode_list = list(np.concatenate(decode_list))
        # print(decode_list)

        return(decode_list)
    
    # ****Plotting the read coverage decoded values over the Extended TSS site**** #

    def Graph_plotting(self, decode_list_value, figure_format):
        x = np.arange(-self.upstream, self.downstream+1) #plotting the coverage
        # print(x)
        plt1.plot(x, decode_list_value)
        # figure = input("Format of Figure to Save, : pdf, jpeg or png")
        if figure_format == "pdf" or figure_format == "png" or figure_format == "jpeg":
            plt1.savefig(self.output_dir+"/Univariate_Graph."+ figure_format)
        elif figure_format == "png and pdf":
            plt1.savefig(self.output_dir+"/Univariate_Graph.png")
            plt1.savefig(self.output_dir+"/Univariate_Graph.pdf")
        elif figure_format == "png and jpeg":
            plt1.savefig(self.output_dir+"/Univariate_Graph.png")
            plt1.savefig(self.output_dir+"/Univariate_Graph.jpeg") 
        elif figure_format == "pdf and jpeg":
            plt1.savefig(self.output_dir+"/Univariate_Graph.jpeg")
            plt1.savefig(self.output_dir+"/Univariate_Graph.pdf")
        elif figure_format == "all three":
            plt1.savefig(self.output_dir+"/Univariate_Graph.png")
            plt1.savefig(self.output_dir+"/Univariate_Graph.pdf")
            plt1.savefig(self.output_dir+"/Univariate_Graph.jpeg")    
        else:
             plt1.savefig(self.output_dir+"/Univariate_Graph.png")
        plt1.show()

# Coverage_mapping(ref_file = '/Users/nehamittal/Desktop/TSS_Gene_Mapping/Reference_genomes/hg38.ncbiRefSeq.gtf', 
                                # sample_file = '/Users/nehamittal/Desktop/TSS_Gene_Mapping/Samples_Data/streck_cfDNA_205_S3__chr21__alignment_summary.csv', output_dir = '/Users/nehamittal/Desktop/Test_results', chr_category = 'chr21', 
                                # feature = 'transcript', upstream = '5000', downstream = '5000', index = '2')


## *** This class will do mapping of the sample rle over a multiple indices using multithreading approach from reference file and plot the graphs *** ##

class Coverage_mapping_multiunit(Preprocessing):
    def __init__(self, ref_file, sample_file, output_dir, chr_category, feature, upstream, downstream, startindex, endindex):
        Preprocessing.__init__(self, ref_file, sample_file, output_dir, chr_category, feature, upstream, downstream)
        self.startindex = int(startindex)
        self.endindex = int(endindex)
        self.chr_category = chr_category
        self.upstream = int(upstream)
        self.downstream = int(downstream)
        self.output_dir = output_dir

    def Coverage_mapping_multi(self, Coverage_read, TSS_extension):
        print("The length of the TSS extension file", len(TSS_extension)) #printing the TSS sites with extennsion using TSS2 variable;
        if (self.startindex < len(TSS_extension)) and (self.endindex < len(TSS_extension)):
            pass
        elif self.startindex >= self.endindex:
            raise ValueError('Start Index should be less than the End Index')
        elif self.startindex >= len(TSS_extension):
            raise ValueError('Start Index should be less than the length of the TSS extension site file')
        else:
            self.endindex = len(TSS_extension)
        print("Printing the end index value after adjustment", self.endindex)
        
        self.decode_dict = {}
        multithread_decode_list = []
        # decode_list_array = np.array([[]])
        # print(type(decode_list_array))

        def Multiple_univariate_mapping(i, num):
        
            while(i < num):
                site_run = Coverage_read[self.chr_category][TSS_extension.Start.values[i]:TSS_extension.End.values[i]] #mapping the coverage of read over gene(s);
                # print(site_run)
                site_run_val = list(site_run.values) #list of value (which value) that has to repeat;
            # print(site_run_val)
                site_run_len = list(site_run.runs) #list of value that how many time above value has to repeat;
            # print(site_run_len)
            # print(np.repeat(site_run_val, site_run_len, axis = 1))
                
                decode_list = []
                for x, y in zip(site_run_val, site_run_len): #decoding the values based on length     
                    z = [x]*y
                    decode_list.append(z) #the values are in array here;
                decode_list = list(np.concatenate(decode_list))
                # print(decode_list)

                self.decode_dict[i] = decode_list
                
                i = i+1


        if __name__ == "TSS_Mapping_End_command":

            per_thread_capacity = 100 #whether it should be user defined or by us the maximum processing capacity of each thread;

            num = self.endindex - self.startindex
            print('total plots', num)

            num_of_threads = num // per_thread_capacity #gives the quotient of division
            index_remaining = num % per_thread_capacity #gives the remainder of division
            print('num_of_threads ', num_of_threads)
            print('index_remaining ', index_remaining)

            if index_remaining == 0:
                pass
            else:
                num_of_threads = num_of_threads + 1 #add 1 more thread to process the remaining index
                print('num_of_threads', num_of_threads)
            
            thread_list = []
            index_pos = self.startindex
            for thr in range(num_of_threads):       

                if index_remaining != 0 and thr == (num_of_threads-1): #executes only for the last loop if there is any remaining index
                    print('remaining_index_thread')
                    thread = threading.Thread(target = Multiple_univariate_mapping, args=(index_pos,index_pos+index_remaining,))
                
                else:
                    # print('main_ thread')
                    thread = threading.Thread(target = Multiple_univariate_mapping, args=(index_pos,index_pos+per_thread_capacity,))
                thread_list.append(thread)
                thread_list[thr].start()
                index_pos = index_pos + per_thread_capacity
                print('thread number ', thr)

            for thread in thread_list:
                thread.join()

        # multithread_decode_list = []
        j = self.startindex
        # x = np.arange(-self.upstream, self.downstream+1)

        # figure = input("Format of Figure to Save, : pdf, jpeg or png")

        while j < self.endindex:
            # print(decode_dict[j])
            # plt1.plot(x, self.decode_dict[j])
            # plt1.xticks(np.arange(-self.upstream, self.downstream+1, 1000))
            
            # if figure == "pdf":
            #     plt1.savefig(self.output_dir+"/Univariate_Graph"+str(j)+ ".pdf")
            # elif figure == "jpeg":
            #     plt1.savefig(self.output_dir+"/Univariate_Graph"+str(j)+ ".jpeg")
            # else:
            #     plt1.savefig(self.output_dir+"/Univariate_Graph"+str(j)+ ".png")

            # plt1.clf()
            # # plt1.show()
            
            multithread_decode_list.append(self.decode_dict[j])
            j = j+1
        # plt1.clf()
        # return multithread_decode_list
    

    def get_Heatmap_multivariate_mapping(self, Yaxis_category, figure_format):
        category = self.Get_TSS1()
        # print(category)

        multithread_decode_list_final = []
        j = self.startindex

        while j < self.endindex:
            
            multithread_decode_list_final.append(self.decode_dict[j])
            j = j+1
        Heatmap_data = pd.DataFrame(multithread_decode_list_final, columns = np.arange(-self.upstream, self.downstream+1))
        Yaxis_category = category[Yaxis_category] #filtering the category (particular column) needs for classification at the Y-axis;
        axis_range = self.endindex - self.startindex
        Yaxis_category[0:axis_range] = Yaxis_category[self.startindex:self.endindex]
        # print(Yaxis_category)

        # creating a dictionary for category:color
        lut = dict(zip(Yaxis_category[0:axis_range].unique(), sn.color_palette("mako", Yaxis_category[0:axis_range].unique().size)))
        row_colors = Yaxis_category[0:axis_range].map(lut)   #corresponding color per column

        # lut = dict(zip(Yaxis_category.unique(), sn.color_palette("mako", Yaxis_category.unique().size)))
        # row_colors = Yaxis_category.map(lut)   #corresponding color per column
        
        heat_map = sn.clustermap(Heatmap_data, xticklabels=500, yticklabels = False, cmap="RdBu_r", figsize = (10, 7), row_cluster = True,
                    col_cluster = False, vmin = 0, vmax = 100, row_colors = row_colors) 
        
        heat_map.fig.suptitle('TSS +/- 5kb', fontsize = 20, fontname = 'Times New Roman', weight = 'bold')
        heat_map.fig.subplots_adjust(right = 0.8)
        heat_map.ax_cbar.set_position((0.9, .4, .03, .4))

        # figure = input("Format of Figure to Save, : pdf, jpeg or png")
        if figure_format == "pdf" or figure_format == "png" or figure_format == "jpeg":
            plt1.savefig(self.output_dir+"/Univariate_Graph."+ figure_format)
        elif figure_format == "png and pdf":
            plt1.savefig(self.output_dir+"/Univariate_Graph.png")
            plt1.savefig(self.output_dir+"/Univariate_Graph.pdf")
        elif figure_format == "png and jpeg":
            plt1.savefig(self.output_dir+"/Univariate_Graph.png")
            plt1.savefig(self.output_dir+"/Univariate_Graph.jpeg") 
        elif figure_format == "pdf and jpeg":
            plt1.savefig(self.output_dir+"/Univariate_Graph.jpeg")
            plt1.savefig(self.output_dir+"/Univariate_Graph.pdf")
        elif figure_format == "all three":
            plt1.savefig(self.output_dir+"/Univariate_Graph.png")
            plt1.savefig(self.output_dir+"/Univariate_Graph.pdf")
            plt1.savefig(self.output_dir+"/Univariate_Graph.jpeg")    
        else:
             plt1.savefig(self.output_dir+"/Univariate_Graph.png")
        plt1.show()


    def get_multiple_univariate_plots(self, figure_format):
        j = self.startindex
        x = np.arange(-self.upstream, self.downstream+1)

        # figure1 = input("Format of Figure to Save, : pdf, jpeg or png")
        
        fig = plt1.figure()

        while j < self.endindex:
            
            plt1.plot(x, self.decode_dict[j])
            plt1.xticks(np.arange(-self.upstream, self.downstream+1, 1000))

            if figure_format == "pdf" or figure_format == "png" or figure_format == "jpeg":
                plt1.savefig(self.output_dir+"/Univariate_Graph"+ str(j)+ "." + figure_format)
            elif figure_format == "png and pdf":
                plt1.savefig(self.output_dir+"/Univariate_Graph" + str(j)+ ".png")
                plt1.savefig(self.output_dir+"/Univariate_Graph" + str(j)+ ".pdf")
            elif figure_format == "png and jpeg":
                plt1.savefig(self.output_dir+"/Univariate_Graph" + str(j)+ ".png")
                plt1.savefig(self.output_dir+"/Univariate_Graph" + str(j)+ ".jpeg") 
            elif figure_format == "pdf and jpeg":
                plt1.savefig(self.output_dir+"/Univariate_Graph" + str(j)+ ".jpeg")
                plt1.savefig(self.output_dir+"/Univariate_Graph" + str(j)+ ".pdf")
            elif figure_format == "all three":
                plt1.savefig(self.output_dir+"/Univariate_Graph" + str(j)+ ".png")
                plt1.savefig(self.output_dir+"/Univariate_Graph" + str(j)+ ".pdf")
                plt1.savefig(self.output_dir+"/Univariate_Graph" + str(j)+ ".jpeg")    
            else:
                plt1.savefig(self.output_dir+"/Univariate_Graph" + str(j)+ ".png")
        
            plt1.clf()
            # plt1.show()
            
            # multithread_decode_list.append(self.decode_dict[j])
            j = j+1
        plt1.close(fig)



    
            
