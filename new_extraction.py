# importing files
import csv
import sqlite3
import subprocess
import concurrent.futures
import threading
import os
import time
import pandas as pd
from datetime import datetime

# initializing variables
name = "ITS_RefSeq_Fungi"
extn = "tar.gz"
db_name = f"./extracted_files/{name}/{name}"
output_file = f"extracted_files/{name}/{name}.tsv"
blastn_file = f"extracted_files/{name}/{name}-blastn.tsv"
sql_file = f"./extracted_files/{name}/taxonomy4blast.sqlite3"
blastn_lock = threading.Lock()

# creating directory in the current directory
def creating_directory():
    print(f"Creating directory named {name}\n")
    create_dir = f"mkdir -p ./extracted_files/{name}"
    result = subprocess.run(create_dir, shell=True, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"Error running command: {result.stderr}")
        exit(1)
    print(f"Successfully created a directory named {name}\n")

# extracting the the compressed file in compressed_files to the newly created directory
def extraction():
    print(f"Extracting {name}.{extn} from compressed_files to current directory\n")
    extract = f"tar -xvzf ./compressed_files/{name}.{extn} -C ./extracted_files/{name}"
    result = subprocess.run(extract, shell=True, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"Error running blast: {result.stderr}")
        exit(1)
    print(f"Successfully extracted {name}.{extn}\n")

# running blastdbcmd command extract data into a temporary file sample.tsv
def blasting():
    print("Running Blast Query to enter data in TSV file\n")
    command = f'blastdbcmd -db {db_name} -entry all -outfmt \'%o,%a,%i,"%t",%s,%g,%l,%h,%T,%X,%e,%L,%C,%S,%N,%B,%K,%P\' > sample.tsv'
    result = subprocess.run(command, shell=True, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"Error running blast: {result.stderr}")
        exit(1)
    print("Done. Entered into sample.tsv\n")

# adding parent-taxid into new field and creating final tsv file
def adding_parent():
    print("Adding parent field to TSV file\n")
    conn = sqlite3.connect(sql_file)
    cursor = conn.cursor()
    print("Loading taxid info into memory for faster lookup...")
    cursor.execute("SELECT taxid, parent FROM TaxidInfo")
    taxid_dict = {str(taxid): str(parent) for taxid, parent in cursor.fetchall()}
    conn.close()
    
    with open("sample.tsv", mode="r", newline="") as infile, open(output_file, mode="w", newline="") as outfile:
        reader = csv.reader(infile)
        writer = csv.writer(outfile, delimiter='\t')
        header = ["ordinal_number", "accession", "sequence_id", "sequence_title", "sequence", "gi", "sequence_length", "sequence_hash_value", "taxid", "taxid_leaf", "membership_integer", "common_taxonomic_name", "common_taxonomic_name_leaf", "scientific_name", "scientific_name_leaf", "blast_name", "taxonomic_super_kingdom", "pig", "taxid_parent"]
        writer.writerow(header)
        next(reader)
        for row in reader:
            taxid = row[9]
            parent = taxid_dict.get(taxid, "")
            row.append(parent)
            writer.writerow(row)
    
    print(f"Added parent field. The new TSV file is {output_file}\n")

# creating header for the blastn file
def blastn_file_creation():
    print(f"Creating {blastn_file} with header\n")
    with open(blastn_file, mode="w", newline="") as write_file:
        writer = csv.writer(write_file, delimiter='\t')
        header = ["query_sequence_id", "query_gi", "query_accession", "query_accession_version", "query_sequence_length", "subject_sequence_id", "subject_all_sequence_id", "subject_gi", "subject_all_gi", "subject_accession", "subject_accession_version", "subject_all_accession", "subject_sequence_length", "query_start", "query_end", "subject_start", "subject_end", "query_sequence", "subject_sequence", "expect_value", "bit_score", "raw_score", "alignment_length", "percentage_identity", "number_of_identical_matches", "number_of_mismatches", "number_of_positive_scoring_matches", "number_of_gap_opens", "number_of_gaps", "percentage_of_positive_scoring_matches", "query/subject_frame", "query_frames", "subject_frames", "blast_traceback_operations", "subject_taxid", "subject_scientific_name", "subject_common_name", "subject_blast_name", "subject_super_kingdom", "subject_all_taxids", "subject_all_scientific_names", "subject_all_common_names", "subject_all_blast_names", "subject_all_super_kingdoms", "subject_strand", "query_coverage_per_subject", "query_coverage_per_hsp", "query_coverage_per_unique_subject", "subject_title", "subject_all_titles"]
        writer.writerow(header)
    print(f"Successfully created {blastn_file} with header\n")

# worker function for blastn()
def run_blastn_for_accession(row):
    accession = row[1]
    blastdbcmd_cmd = f"blastdbcmd -db {db_name} -entry {accession}"
    blastn_cmd = f'blastn -db {db_name} -outfmt "6 qseqid qgi qacc qaccver qlen sseqid sallseqid sgi sallgi sacc saccver sallacc slen qstart qend sstart send qseq sseq evalue bitscore score length pident nident mismatch positive gapopen gaps ppos frames qframe sframe btop staxid ssciname scomname sblastname sskingdom staxids sscinames scomnames sblastnames sskingdoms sstrand qcovs qcovhsp qcovus stitle salltitles" -max_target_seqs 10'
    try:
        print(f"Running blastn for accession: {accession}")
        blastdbcmd_output = subprocess.run(
            blastdbcmd_cmd,
            shell=True,
            check=True,
            capture_output=True,
            text=True
        )
        blastn_output = subprocess.run(
            blastn_cmd,
            input=blastdbcmd_output.stdout,
            shell=True,
            check=True,
            capture_output=True,
            text=True
        )
        with blastn_lock:
            with open(blastn_file, 'a', newline='') as f:
                f.write(blastn_output.stdout)
    except subprocess.CalledProcessError as e:
        print(f"Error processing accession {accession}: {e.stderr}")

def blastn():
    print("Running blastn for each accession and appending to blastn TSV file\n")
    blastn_file_creation()
    with open(output_file, mode="r", newline="") as read_file:
        reader = csv.reader(read_file, delimiter='\t')
        accession_rows = list(reader)[1:]
    
    max_workers = 14
    with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = [executor.submit(run_blastn_for_accession, row) for row in accession_rows]
        concurrent.futures.wait(futures)

    print("All blastn tasks completed successfully!\n")

# moving the compressed file from Downloads to compressed_files
def moving():
    print(f"Moving {name}.{extn} from Downloads to compressed_files\n")
    source_path = os.path.expanduser(f"~/Downloads/{name}.{extn}")
    dest_dir = "./compressed_files"
    dest_path = os.path.join(dest_dir, f"{name}.{extn}")
    if os.path.exists(source_path):
        try:
            import shutil
            shutil.move(source_path, dest_path)
            print(f"Successfully moved {name}.{extn} from Downloads to compressed_files\n")
        except Exception as e:
            print(f"Error moving file: {e}")
            exit(1)
    else:
        print(f"Warning: Source file {source_path} not found. Skipping move operation.")

# removing the temporary files
def removing_file():
    print("Removing sample.tsv file")
    rm_file = "rm sample.tsv"
    result = subprocess.run(rm_file, shell=True, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"Error running blast: {result.stderr}")
        exit(1)
    print("Successfully removed sample.tsv file")

# removing directory containing extracted files
def removing_directory():
    print(f"Removing {name} directory\n")
    rm_dir = f"rm -r {name}"
    result = subprocess.run(rm_dir, shell=True, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"Error running blast: {result.stderr}")
        exit(1)
    print(f"Successfully removed {name} directory\n")

# --- Main Execution ---
if __name__ == "__main__":
    start_datetime = datetime.now()
    
    creating_directory()
    extraction()
    blasting()
    adding_parent()
    blastn()
    # moving()
    # removing_file()
    # removing_directory()
    
    end_datetime = datetime.now()
    total_duration = end_datetime - start_datetime

    # Format timedelta to show days, hours, minutes, seconds, milliseconds
    total_seconds = total_duration.total_seconds()
    days, remainder = divmod(total_seconds, 86400)
    hours, remainder = divmod(remainder, 3600)
    minutes, remainder = divmod(remainder, 60)
    seconds, milliseconds = divmod(remainder, 1)

    duration_str = []
    if int(days) > 0:
        duration_str.append(f"{int(days)} days")
    if int(hours) > 0:
        duration_str.append(f"{int(hours)} hours")
    if int(minutes) > 0:
        duration_str.append(f"{int(minutes)} minutes")
    if int(seconds) > 0:
        duration_str.append(f"{int(seconds)} seconds")
    if int(milliseconds * 1000) > 0:
        duration_str.append(f"{int(milliseconds * 1000)} milliseconds")
    
    formatted_duration = ", ".join(duration_str) if duration_str else "less than a millisecond"
    
    print("\n" + "="*40)
    print("       Script Execution Summary       ")
    print("="*40)
    print(f"Started on: {start_datetime.strftime('%A, %B %d, %Y %H:%M:%S')}")
    print(f"Ended at:   {end_datetime.strftime('%A, %B %d, %Y %H:%M:%S')}")
    print(f"Total time taken: {formatted_duration}")
    print("="*40 + "\n")
    print("All tasks completed successfully!")