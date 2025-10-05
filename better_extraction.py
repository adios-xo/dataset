import csv
import subprocess
import concurrent.futures
import os

# ----------------------------------------------------------------
# Function: fetch_sequences
# Purpose: For a list of accession IDs, fetch corresponding sequences
#          from the local BLAST database and write them to a FASTA file.
# ----------------------------------------------------------------
def fetch_sequences(batch_accessions, db_path, fasta_file):
    with open(fasta_file, "w") as f_out:
        for acc in batch_accessions:
            try:
                # Using blastdbcmd to fetch sequence by accession ID
                result = subprocess.run(
                    ['blastdbcmd', '-db', db_path, '-entry', acc, '-outfmt', '%f'],
                    capture_output=True, text=True, check=True
                )
                # Append to the batch FASTA file
                f_out.write(result.stdout)
            except subprocess.CalledProcessError:
                # If sequence not found or any error, skip it
                print(f"[Warning] Could not fetch sequence for accession: {acc}")

# ----------------------------------------------------------------
# Function: run_blast_and_parse
# Purpose: Run blastn on the given FASTA file and return parsed results.
#          The result is a list of rows, where each row is a list of BLAST fields.
# ----------------------------------------------------------------
def run_blast_and_parse(fasta_file, db_path):
    # These are the fields we want in the BLAST output (format 6 = tab-delimited)
    outfmt_fields = [
        "qseqid", "qgi", "qacc", "qaccver", "qlen", "sseqid", "sallseqid", "sgi", "sallgi",
        "sacc", "saccver", "sallacc", "slen", "qstart", "qend", "sstart", "send", "qseq",
        "sseq", "evalue", "bitscore", "score", "length", "pident", "nident", "mismatch",
        "positive", "gapopen", "gaps", "ppos", "frames", "qframe", "sframe", "btop", "staxid",
        "ssciname", "scomname", "sblastname", "sskingdom", "staxids", "sscinames", "scomnames",
        "sblastnames", "sskingdoms", "sstrand", "qcovs", "qcovhsp", "qcovus", "stitle", "salltitles"
    ]

    outfmt_string = "6 " + " ".join(outfmt_fields)

    # Run blastn with the given output format
    result = subprocess.run(
        ['blastn', '-query', fasta_file, '-db', db_path, '-outfmt', outfmt_string],
        capture_output=True, text=True, check=True
    )

    # Split the output into lines and then fields
    lines = result.stdout.strip().split('\n')
    parsed_rows = [line.split('\t') for line in lines if line.strip()]

    # Wrap the last two columns (stitle, salltitles) in quotes for clean CSV formatting
    for row in parsed_rows:
        if len(row) >= 2:
            row[-2] = f'"{row[-2]}"'
            row[-1] = f'"{row[-1]}"'

    return parsed_rows

# ----------------------------------------------------------------
# Function: process_batch
# Purpose: Fetch sequences, run BLAST, parse results for a single batch
# ----------------------------------------------------------------
def process_batch(batch_accessions, db_path, batch_num):
    # Temporary FASTA file name for this batch
    fasta_file = f"batch_{batch_num}.fasta"

    print(f"[Batch {batch_num}] Fetching {len(batch_accessions)} sequences...")
    fetch_sequences(batch_accessions, db_path, fasta_file)

    print(f"[Batch {batch_num}] Running BLAST...")
    blast_results = run_blast_and_parse(fasta_file, db_path)

    # Clean up temporary file
    try:
        os.remove(fasta_file)
    except OSError:
        pass

    print(f"[Batch {batch_num}] Completed. Found {len(blast_results)} hits.")
    return blast_results

# ----------------------------------------------------------------
# Function: main
# Purpose: Read accessions, batch them, process in parallel, write to CSV
# ----------------------------------------------------------------
def main():
    # --- Configuration ---
    input_csv = "accessions.csv"           # CSV with one accession ID per row
    output_csv = "blast_results.csv"       # Final output file
    blast_db_path = "your_blast_db"        # Path to local BLAST database
    batch_size = 50                        # How many accessions per batch
    max_workers = 14                       # Number of parallel threads to use

    # --- Step 1: Load accession IDs from CSV ---
    with open(input_csv, newline='') as f:
        reader = csv.reader(f)
        accession_ids = [row[0].strip() for row in reader if row]

    print(f"[Info] Loaded {len(accession_ids)} accession IDs.")

    # --- Step 2: Split into batches ---
    batches = [accession_ids[i:i+batch_size] for i in range(0, len(accession_ids), batch_size)]
    print(f"[Info] Created {len(batches)} batches of size ~{batch_size}.")

    # --- Step 3: Define BLAST output fields (same as in parsing) ---
    header_fields = [
        "qseqid", "qgi", "qacc", "qaccver", "qlen", "sseqid", "sallseqid", "sgi", "sallgi",
        "sacc", "saccver", "sallacc", "slen", "qstart", "qend", "sstart", "send", "qseq",
        "sseq", "evalue", "bitscore", "score", "length", "pident", "nident", "mismatch",
        "positive", "gapopen", "gaps", "ppos", "frames", "qframe", "sframe", "btop", "staxid",
        "ssciname", "scomname", "sblastname", "sskingdom", "staxids", "sscinames", "scomnames",
        "sblastnames", "sskingdoms", "sstrand", "qcovs", "qcovhsp", "qcovus", "stitle", "salltitles"
    ]

    # Write header to CSV output file
    with open(output_csv, 'w', newline='', encoding='utf-8') as f_out:
        writer = csv.writer(f_out)
        writer.writerow(header_fields)

    # --- Step 4: Process batches in parallel using ThreadPool ---
    # We use threads because BLAST runs in subprocesses (I/O-bound),
    # and ThreadPool works better in this case than ProcessPool on many systems.
    with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = {
            executor.submit(process_batch, batch, blast_db_path, idx): idx
            for idx, batch in enumerate(batches, start=1)
        }

        # Open output CSV in append mode to add results as they come in
        with open(output_csv, 'a', newline='', encoding='utf-8') as f_out:
            writer = csv.writer(f_out)
            for future in concurrent.futures.as_completed(futures):
                batch_num = futures[future]
                try:
                    results = future.result()
                    writer.writerows(results)
                    print(f"[Batch {batch_num}] Results written to CSV.")
                except Exception as e:
                    print(f"[Error] Batch {batch_num} failed: {e}")

    print("[Done] All batches processed and results saved.")

# Entry point
if __name__ == "__main__":
    main()
