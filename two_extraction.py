# importing files
import csv
import sqlite3
import subprocess

# initializing variables
name = "ITS_RefSeq_Fungi"
extn = "tar.gz"
db_name = f"{name}/{name}"    # Enter the path to your DB
output_file = f"csv_files/{name}.tsv"   # Changed to .tsv
blastn_file = f"blastn/{name}-blastn.tsv"	# Changed to .tsv
sql_file = f"./{name}/taxonomy4blast.sqlite3"

# creating directory in the current directory
def creating_directory():
	print(f"Creating directory named {name}\n")
	create_dir = f"mkdir {name}"
	result = subprocess.run(create_dir, shell=True, capture_output=True, text=True)
	if result.returncode != 0:
		print(f"Error running blast: {result.stderr}")
		exit(1)
	print(f"Successfully created a  directory named {name}\n")

# extracting the the compressed file in Downloads to the newly created directory
def extraction():
	print(f"Extracting {name}.{extn} from Downloads to current directory\n")
	extract = f"tar -xvzf ~/Downloads/{name}.{extn} -C ./{name}"
	result = subprocess.run(extract, shell=True, capture_output=True, text=True)
	if result.returncode != 0:
		print(f"Error running blast: {result.stderr}")
		exit(1)
	print(f"Successfully extracted {name}.{extn} from Downloads to current directory\n")

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

	with open("sample.tsv", mode="r", newline="") as infile, open(output_file, mode="w", newline="") as outfile:
		reader = csv.reader(infile)
		writer = csv.writer(outfile, delimiter='\t')  # Changed to tab delimiter
		header = ["ordinal_number", "accession", "sequence_id", "sequence_title", "sequence", "gi", "sequence_length", "sequence_hash_value", "taxid", "taxid_leaf", "membership_integer", "common_taxonomic_name", "common_taxonomic_name_leaf", "scientific_name", "scientific_name_leaf", "blast_name", "taxonomic_super_kingdom", "pig", "taxid_parent"]
		writer.writerow(header)

		for row in reader:
			taxid = row[9]
			query = "SELECT parent FROM TaxidInfo WHERE taxid = ?"
			cursor.execute(query, (taxid,))
			result = cursor.fetchone()
			parent = str(result[0]) if result else ""
			row.append(parent)
			writer.writerow(row)
	
	conn.close()
	print(f"Added parent field. The new TSV file is {output_file}\n")

def blastn_file_creation():
	print(f"Creating {blastn_file} with header\n")

	with open(blastn_file, mode="w", newline="") as write_file:
		writer = csv.writer(write_file, delimiter='\t')  # Changed to tab delimiter
		header = ["query_sequence_id", "query_gi", "query_accession", "query_accession_version", "query_sequence_length", "subject_sequence_id", "subject_all_sequence_id", "subject_gi", "subject_all_gi", "subject_accession", "subject_accession_version", "subject_all_accession", "subject_sequence_length", "query_start", "query_end", "subject_start", "subject_end", "query_sequence", "subject_sequence", "expect_value", "bit_score", "raw_score", "alignment_length", "percentage_identity", "number_of_identical_matches", "number_of_mismatches", "number_of_positive_scoring_matches", "number_of_gap_opens", "number_of_gaps", "percentage_of_positive_scoring_matches", "query/subject_frame", "query_frames", "subject_frames", "blast_traceback_operations", "subject_taxid", "subject_scientific_name", "subject_common_name", "subject_blast_name", "subject_super_kingdom", "subject_all_taxids", "subject_all_scientific_names", "subject_all_common_names", "subject_all_blast_names", "subject_all_super_kingdoms", "subject_strand", "query_coverage_per_subject", "query_coverage_per_hsp", "query_coverage_per_unique_subject", "subject_title", "subject_all_titles"]
		writer.writerow(header)

	print(f"Successfully created {blastn_file} with header\n")

# runnning blastn command
def blastn():
	print("Running blastn for each accession and appending to blastn TSV file\n")

	with open(output_file, mode="r", newline="") as read_file, open(blastn_file, mode="a", newline="") as write_file:
		reader = csv.reader(read_file, delimiter='\t')  # Changed to tab delimiter
		write_blastn = csv.writer(write_file, delimiter='\t')  # Changed to tab delimiter
		is_header = True

		for row in reader:
			if is_header:
				is_header = False
				continue
			
			print(f"Running blastn for accession: {row[1]}\n")
			accession = row[1]

			blast = f"blastdbcmd -db {name}/{name} -entry {accession} -out query.fasta"
			result = subprocess.run(blast, shell=True, capture_output=True, text=True)

			if result.returncode != 0:
				print(f"Error running blast: {result.stderr}")
				exit(1)
			
			blastn = f'blastn -query query.fasta -db {name}/{name} -outfmt "6 qseqid qgi qacc qaccver qlen sseqid sallseqid sgi sallgi sacc saccver sallacc slen qstart qend sstart send qseq sseq evalue bitscore score length pident nident mismatch positive gapopen gaps ppos frames qframe sframe btop staxid ssciname scomname sblastname sskingdom staxids sscinames scomnames sblastnames sskingdoms sstrand qcovs qcovhsp qcovus stitle salltitles" -max_target_seqs 10 -out bn.tsv'
			blastn_result = subprocess.run(blastn, shell=True, capture_output=True, text=True)

			if blastn_result.returncode != 0:
				print(f"Error running blast: {blastn_result.stderr}")
				exit(1)
			
			print(f"Appending results of accession: {row[1]} to {blastn_file}\n")
			with open("bn.tsv", mode="r", newline="") as read_blastn_file:
				reader = csv.reader(read_blastn_file, delimiter='\t')  # Fixed filename and delimiter
				data = list(reader)
				write_blastn.writerows(data)

# moving the compressed file from Downloads to compressed_files
def moving():
	print(f"Moving {name}.{extn} from Downloads to compressed_files\n")
	moving = f"mv ~/Downloads/{name}.{extn} ./compressed_files"
	result = subprocess.run(moving, shell=True, capture_output=True, text=True)
	if result.returncode != 0:
		print(f"Error running blast: {result.stderr}")
		exit(1)
	print(f"Successfully moved {name}.{extn} from Downloads to compressed_files\n")

# removing the temporary files
def removing_file():
	print("Removing sample.tsv file")  # Updated print message
	rm_file = "rm sample.tsv"  # Fixed file name
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


creating_directory()
extraction()
blasting()
adding_parent()
blastn_file_creation()
blastn()
moving()
removing_file()
removing_directory()
print("All tasks completed successfully!")