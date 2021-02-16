""" Script to process data and generate a text file 
Author: Nick Sipes
Email: ntsipes@gmail.com
Date: 02-15-2021
"""

import os
import argparse
import pandas as pd
import logging
from datetime import datetime


def find_master_file(list_of_files):
    """ Verify total_reads_all_sample file exists and remove from text file list"""
    # Iterate through all files
    for f in list_of_files:
        # Look for the file, ignoring case
        if "total_reads_all_sample.txt" in f.lower(): 
            logging.info("total reads file found! \n")
            # Store the filepath
            total_reads_fp = f
            logging.debug("\ttotal reads fp {}".format(total_reads_fp))
            # Remove the file from the list
            list_of_files.remove(f)
            # Exit the loop
            break
    # Return the list
    return list_of_files, total_reads_fp

def match_csvs_to_sample_depths(list_of_csvs, total_reads_fp):
    """ Get sample depths from the total reads file """
    # Store csvs with matching total number reads value
    filenames_with_total_num_reads = []
    # Read the file as a csv with tab delimiter
    data = pd.read_csv(total_reads_fp, delimiter="\t")
    # Iterate through the list of csvs
    while len(list_of_csvs) != 0:
        f = list_of_csvs.pop()
        # Get just the filename
        path, filename = os.path.split(f)
        # Find the matching sample depth by filename
        depth, idx = get_depth_for_filename(filename, total_reads_fp)
        logging.debug("Match found!\n\tDepth: {} Filename: {}\n".format(depth, filename))
        # Add filename, depth, and total number of reads to list as a tuple (immutable list)
        filenames_with_total_num_reads.append((f, str(depth), data["total # reads/1000000"][idx]))
    # Return the matches
    logging.debug("Filenames with their matching total number of reads:\n\t{}\n"
                        .format(filenames_with_total_num_reads))
    return filenames_with_total_num_reads

def generate_master_csv(list_of_csvs_with_data):
    """ Perform math on csvs matched to their total number of reads """
    # Define column names for all CSVs
    col_names = ["contigs", "genelength", "maps", "null"]
    # Create the master CSV
    master_csv = None 
    # Iterte through CSV's with their matches
    for f in list_of_csvs_with_data:
        # Read the CSV
        temp = pd.read_csv(f[0], delimiter="\t")
        # Apply the column names
        temp.columns = col_names
        if master_csv is None:
            # Make the master_csv a dataFrame
            master_csv = pd.DataFrame()
            # Create the first column in the master_csv
            master_csv.insert(0, "Reference Contig", temp["contigs"])
        # Perform math using data from csv and the total number reads from the related sample depth
        math_result = temp["maps"] / ((temp["genelength"] / 1000.0) * f[2])
        # Insert the result into the master_csv
        master_csv.insert(len(master_csv.columns), "{}".format(f[1]), math_result)
    # Create a CSV from the master_csv dataFrame object
    data_directory_path = os.path.split(list_of_csvs_with_data[0][0])
    master_csv.to_csv(os.path.join(data_directory_path[0], "master_list.csv"))
    # Return the master_csv dataFrame
    return master_csv    

def find_file_type_in_dir(directory, extension_to_find, exclusions=[]):
    """ Find files that end with extension """
    # List all files in directory
    list_of_files = os.listdir(directory)
    # Create an empty list to store results
    matches = []
    # Iterate through all files
    for f in list_of_files:
        # Find where the file extension starts in the file name
        extension_start = f.find(".")
        # Compare the file extension to extension_to_find
        if f[extension_start:] in extension_to_find and f not in exclusions:
            # Append with the full path
            abs_filepath = os.path.join(directory, f)
            matches.append(abs_filepath)
            # Display matches if debugging
            logging.debug("\tMatch found: {}".format(abs_filepath))
    # Show the filenames for visbility
    display_filenames(matches)
    # Return all the matches
    return matches

def display_filenames(list_of_filenames):
    """ Display filenames from a list of filepaths """
    for f in list_of_filenames:
        result = os.path.split(f)
        logging.info("\tfilename is {}".format(result[1]))

def get_depth_for_filename(filename, total_reads_fp):
    """ Extract the depth from the filename """
    # Get all available depths
    data = pd.read_csv(total_reads_fp, delimiter="\t")
    available_depths = data["Sample Depth"]
    # Look for available depths in the filename
    for idx, depth in enumerate(available_depths):
        depth_str = (str(depth).replace(".", ""))
        # Check if a trailing 0 was added and remove it
        if depth_str[-1] == "0":
            depth_str = depth_str[:-1]
        # Find where the depth starts in the filename
        start_idx = filename.find(depth_str)
        # Check the depth_str is in the filename and the character before the depth is not numeric
        if depth_str in filename and not filename[start_idx - 1].isnumeric():
            #logging.debug("Match found! {} {}".format(depth_str, filename))
            return depth, idx

def aggregate_text_files(list_of_text_files, total_reads_fp):
    """ Combine all text files into a list to make matching more efficient """
    # Empty list to store the mega list in
    text_file_mega_list = []
    for f in list_of_text_files:
        # Get the filename
        path, name = os.path.split(f)
        # Read lines from the file
        lines = open(f, 'r').readlines()
        # Append to the mega list as a tuple
        for line in lines:
            depth, _ = get_depth_for_filename(name, total_reads_fp)
            text_file_mega_list.append((name, depth, line.strip("\n")))
    logging.debug("Mega list is \t {}".format(text_file_mega_list))
    return text_file_mega_list

def match_contigs_and_sum(agg_text_file, master_dataframe):      
    """ Sum the depth values of contigs in text file from master list"""
    ref_contig_col = master_dataframe["Reference Contig"]
    # Create a list to hold the dictionaries for each text file
    text_file_totals = {}
    for idx, item in enumerate(ref_contig_col):
        for line in agg_text_file:
            # Find the row that matches the NODE from the line
            if line[2] in item:
                logging.info("\tMatch found! {} and {}".format(line[1], item))
                # Get the row of the NODE as a series
                node_data = master_dataframe.iloc[ref_contig_col[ref_contig_col == item].index[0]]
                # Check if the text file total is being tracked
                if line[0] in text_file_totals.keys():
                    # Accumulate the total
                    text_file_totals[line[0]] = text_file_totals[line[0]].add(node_data)
                else:
                    # Start tracking the total
                    text_file_totals.update({line[0] : node_data})
    return text_file_totals

def output_files(results):
    """ Clean the process results to generate a text and csv file """
    data = []
    stop_enabled = True
    for n in results.keys():
        # Get the current line 
        line_results = results[n]
        # Rename it based on the text file
        line_results["Reference Contig"] = n
        # Make the Series into a Dataframe to transpose
        line_results = line_results.to_frame().transpose()
        # Print the result
        print(f"Key: {n} Data: {line_results}")
        data.append(line_results)
        if stop_enabled and input("Press 's' + ENTER to skip line-by-line results") == 's':
            stop_enabled = False
    # Concatinate all the results and sort columns in ascending order
    combined_result = pd.concat(data).sort_index(axis=1)
    # Move the text file column to the front of the dataframe
    text_file_column = combined_result["Reference Contig"]
    combined_result.drop(labels=["Reference Contig"], axis=1, inplace=True)
    combined_result.insert(0, "Reference Contig", text_file_column)
    # Generate a timestamp
    timestamp = datetime.now().strftime("%Y%m%d_%H%M")
    # Save the output csv to the current working directory
    output_fp = os.path.join(os.getcwd(), f"{timestamp}_process_results.csv")
    combined_result.to_csv(output_fp)
    logging.info("\t*****-----PROCESS COMPLETE-----******")
    logging.info(f"Output data can be found at: {output_fp}")
    return output_fp

# Set our logger visibility level
logging.basicConfig(level=logging.DEBUG)

# Get the input filepath
args = argparse.ArgumentParser(description="Data to process")
args.add_argument("--data-dir", type=str,
                    help="Full filepath of directory containing all txt and csv data")
args = args.parse_args()

### READ THIS ###
# If you want to ignore files in your data directory, add them to this list #
files_to_ignore =[ "master_list.csv", "siberiangenomes__sed.txt" ]


logging.info("Path to access data is {}".format(args.data_dir))
# Search for csv files
csv_files = find_file_type_in_dir(args.data_dir, ".csv", exclusions=files_to_ignore)
# Search for text files
txt_files = find_file_type_in_dir(args.data_dir, ".txt", exclusions=files_to_ignore)
# Find the total reads file and get it's path
txt_files, total_reads_fp = find_master_file(txt_files)
# Match CSV's to their sample depth information
csvs_with_sample_depths = match_csvs_to_sample_depths(csv_files, total_reads_fp)
# Create the master csv from all the csvs
master_csv_dataframe = generate_master_csv(csvs_with_sample_depths)
# Aggregate the text files to make searching faster
aggregated_text_file_lines = aggregate_text_files(txt_files, total_reads_fp)
# Match text file contents to sum the proper values
results = match_contigs_and_sum(aggregated_text_file_lines, master_csv_dataframe)
# Convert the result into output files
output_filepath = output_files(results)
