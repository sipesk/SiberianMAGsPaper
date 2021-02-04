'''this will take the outputs from bowtie and make a final table to be used in R (or other)
for abundance. 
Inputs will be the .txt files of the names in each of the
contigs for each of the MAGs, the total # of reads in the metagenome,
and to calculate the RPKM (reads per kilobase per million reads'''
import os, argparse, sys
import pandas

def finding(data_dir, metagenome_sample_info):
    '''function is to look through the the data dir and find 1) total reads in the metagenomes 2) all the bin files(.txt)
    and 3) all the sample files with the info'''
    #file end with .txt into a list
    string_text = '.txt'  
    #files end with .csv into a list  
    string_csv = '.csv'    
    #dir with all the files needed 
    filelist = os.listdir(data_dir)     
    #iterate through all files in data_dir
    #start empty csv list
    csvs = []                   
    #start empty text list
    texts = []                  
    sample_info = False
    #iterate through the files in the dir and
    for f in filelist:   
        #if the file it the inputed metagenome sample info       
        if f == metagenome_sample_info: 
             # then change Boliean to True                
            sample_info = True   
            # if the files end with .csv                      
        elif f.endswith(string_csv):  
             # add those files to the above empty csvs list                  
            csvs.append(f)   
             # if the files end with .txt                          
        elif f.endswith(string_text):  
            # add those files to the aboce empty txt list                
            texts.append(f) 
     # if not the metagenome sample info file the return                             
    if not sample_info:                                
        print("Your {} was not found in {}, check names".format(metagenome_sample_info, data_dir))
        sys.exit(1)
    else:
        #print("These are the csv files: {}, \nthese are the text files {}, ...\nthis is your metagenome sample info {}".format(csvs, texts, metagenome_sample_info))
    #the lists of csvs and test file names
        return csvs, texts                                 

#function takes the csv files, the numeric sample names and the metagenome sample info,
def sample_match(csv_file, sample_depths, number_reads):            
    print("Finding match for {}.".format(csv_file))
    # for each item and it's value in the sample depth number list
    for index, value in enumerate(sample_depths): 
        # start_index is equal to the index where the csv file name matches the sample depth number                  
        start_index = csv_file.find(value)           
        # if the value in the csv file name AND the index befor the start index is not numeric               
        if value in csv_file and not csv_file[start_index - 1].isnumeric():  
            # then print the index (the matching potition, if it doesnt match it will return -1)       
            print(index)                           
            # set target_number reads to the division of the total # reads/1million of each of the files iterated                                     
            target_num_reads = number_reads["total # reads/1000000"][index]     
            # print the information of the matching csv, start index and the sample depth number       
            print("Csv file: {}\n Start index: {}\n value is: {}\n".format(csv_file, start_index, value))  
             #return the divion and the sample depth number   
            return target_num_reads, value                                                                    

#Function to do math on all of the csv files in the data dir and metagenomic sample info
def csv_math(list_of_csvs, data_dir, number_reads): 
    # set empty list to sample depth, this will be populated with the names of the samples           
    sample_depths = []
    # Set master list as None in order to make it a populated dataframe
    master_list = None
    #set number_reads to the file for the metagenome_sample_info
    number_reads = pandas.read_csv(os.path.join(data_dir, number_reads), delimiter="\t")   
    # for each of the sample names in the metagenom..._info in the column "Sample Depth" 
    for sample in number_reads["Sample Depth"]:
        # remove the '.' in the name and replace it with nothing
        temp_sample = (str(sample).replace(".", ""))
        # Check if we added a zero decimal
        if temp_sample[-1] is "0":
            # Remove the zero
            temp_sample = temp_sample[:-1]
        # Add to the list
        sample_depths.append(temp_sample)
    for i in list_of_csvs:
        #for each of the files in the csv file list
        #call target_num_reads and the value of the sample names to be the output of the sample_match func
        target_num_reads, value = sample_match(i, sample_depths, number_reads)
        # temp file set to each csv file in the csv file list in order to add more, sep with a tab
        temp = pandas.read_csv(os.path.join(data_dir, i), delimiter="\t")
         # define temp col names
        temp.columns = ["contigs", "genelength", "maps", "null"]
        # if the master_list is none, or just starting
        if master_list is None:
            # the master_list will be a dataframe
            master_list = pandas.DataFrame()
             # insert each infomation from each csv file into the master_list at the first pos, named whats in the "" and the input is the final argument
            master_list.insert(0, "Reference Contig", temp["contigs"])
        result = temp["maps"]/((temp["genelength"]/1000)*target_num_reads)
    # set result to be the equation of the maps in the temp file / (genelength / thousand) * target num reads - because that value is defined above
    #Insert the result are the end of the length of the master list
        master_list.insert(len(master_list.columns), "{}".format(value), result)       
    #print the master_list to the output
    print(master_list)
        master_list.to.csv(os.path.join(data_dir, "master_list.csv")                                                 
    return master_list

# function to match the contigs that exist in a bin with the master list
def contig_match(texts, master_list, input_data_dir):
    ''' Get depth values of contigs in text files from master list and sum '''       
    #Open text files
    text_file = open(text_file, 'r')
    #Geth the lines of the text file
    lines = text_file.readlines()
    #create a list to store the indices of matching contigs
    matches = []
    for line in lines:
        match = (master_list.loc[master_list["Reference Contig"] == line])
        print("Line from text file is {}\n Match from master list is {}\n".format(line,match))
        break

    for text_file in texts: 
        # for each bin in the text file,                  
        text_file = open(os.path.join(input_data_dir, text_file), 'r')                       
         # open the txt file
        lines = text_file.readlines()                 

#main function to pass arguments to the above funtions
# a and b are set outputs for the finding function which takes the data input dir and the metagenome sample info file
# csv math function is called on the input data dir and the input metagenome sample info file
def main(data_dir, metagenome_sample_info):                     
    '''Entry point for application'''                   
    a, b = finding(data_dir, metagenome_sample_info)             
    master_list = csv_math(a, data_dir, metagenome_sample_info)
    for text in text_files:
        text = os.path.join(input_data_dir, text)   
        temp_file = contig_match(text, master_list, data_dir)


# parser = argparse.ArgumentParser()
# parser.add_argument('data_dir', 
#                     help= "Directory that holds all of the Bowtie2 outputs")
# parser.add_argument('metagenome_sample_info', 
#                     help= "All reads w/in each of your inputed metagenome samples")
# args=parser.parse_args()
input_data_dir = '/Users/ksipes/Documents/UTK/DimentionsPermafrost/bowtie/NEW-bowtie'
input_metagenome_sample_info = 'Total_Reads_ALL_Sample.txt'
main(input_data_dir, input_metagenome_sample_info)