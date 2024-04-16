"""Module to handle blast requests"""

import os
import subprocess
import shutil
import tempfile

from .config import db_config_array, db_alt_config_array

class Blast:
    """Class for Blast arguments"""
    def __init__(self, db, threads):
        self.db = db
        self.threads = threads

    def parse_hits(self):
        # Parse BLAST results and extract relevant information
        # You can use Biopython's BLAST parser or custom parsing logic here
        pass

    def determine_rfb(self):
        # Implement logic to determine the rfb type
        pass

    def determine_mlst(self):
        # Implement logic to determine MLST based on BLAST results
        pass

    def determine_flic(self):
        # Implement logic to determine fliC type
        pass

    def determine_crispr(self):
        # Implement logic to determine CRISPR type
        pass

    def combine_data(self):
        # Combine the results to infer the serotype
        pass

    def write_out(self, count_dict, outdir, input_filename, db_marker):
        with open(os.path.join(outdir, input_filename, f"{db_marker}_hits.txt"), 'w') as fout:
            for gene, count in count_dict.items():
                fout.write(f"{gene};{count}\n")

    def coverage_hits(self, outdir, input_filename, db_marker):
        hits_fpath = os.path.join(outdir, input_filename, f"{db_marker}_hits.txt")
        rfb_hits_count_fpath = os.path.join(self.db, 'RFB_hits_count.csv')
        hitscov_fpath = os.path.join(outdir, input_filename, f"{db_marker}_hitscoverage.txt")
        with open(hits_fpath, 'r') as fin1, \
            open(rfb_hits_count_fpath, 'r') as fin2, \
            open(hitscov_fpath, 'w') as fout:
            for line in fin1:
                parts1 = line.rstrip().split(';')
                for line2 in fin2:
                    parts2 = line2.rstrip().split(';')
                    if parts1[0] == parts2[0]:
                        coverage = float(parts1[1]) / float(parts2[1]) * 100
                        fout.write(f"{line.rstrip()};{parts2[1]};{coverage}\n")
                        break

    def run_blastn(self, outdir, input_filename, input_fasta_file, db_fasta, db_marker, identity, coverage):
        blast_cmd = f"blastn -db {self.db}/{db_fasta} -query {input_fasta_file} -out {os.path.join(outdir, input_filename, f'{db_marker}_blastout.txt')} -num_threads {self.threads} -num_alignments 10000 -outfmt 6 -word_size 11 -dust no"
        subprocess.run(blast_cmd, shell=True, check=True)
        self.filter_blast(outdir, input_filename, db_marker, identity, coverage)

    def filter_blast(self, outdir, input_filename, db_marker, identity, coverage):
        with open(os.path.join(outdir, input_filename, f'{db_marker}_blastout.txt'), 'r') as blastfile, \
            open(os.path.join(outdir, input_filename, f'{db_marker}_allrecords.txt'), 'w') as fout:
            for line in blastfile:
                parts = line.rstrip().split('\t')
                length = float(parts[3])
                q_start = float(parts[6])
                q_end = float(parts[7])
                q_len = q_end - q_start + 1
                if float(parts[2]) >= identity and (length / q_len) * 100 >= coverage:
                    fout.write(line)

    def get_hits(self, outdir, input_filename, db_marker):
        gene_counts = {}
        with open(os.path.join(outdir, input_filename, f"{db_marker}_allrecords.txt"), 'r') as fin:
            for line in fin:
                parts = line.strip().split('\t')
                gene = parts[1].split('_')[0]
                gene_counts[gene] = gene_counts.get(gene, 0) + 1
        sorted_gene_counts = dict(sorted(gene_counts.items(), key=lambda item: item[1], reverse=True))
        self.write_out(sorted_gene_counts, outdir, input_filename, db_marker)
    
    def make_blast_db(self):
        for root, dirs, files in os.walk(self.db):
            for filename in files:
                if filename.endswith('.fasta'):
                    subprocess.run(f"makeblastdb -dbtype nucl -in {os.path.join(root, filename)}", shell=True, check=True)

    def run(self, list_file, outdir, keep_files):
        with open(list_file, 'r') as f:
            input_files = f.read().splitlines()

            for filepath in input_files:
                input_filename = os.path.splitext(os.path.basename(filepath))[0]
                tempdir = tempfile.mkdtemp(prefix="ShigaPass_")
                input_fasta_file = os.path.join(tempdir, input_filename + "_parsed.fasta")
                shutil.copy(filepath, input_fasta_file)
                new_outdir = os.path.join(outdir, input_filename)
                if not os.path.exists(new_outdir):
                    os.makedirs(new_outdir)

                for search_parameters in db_config_array:
                    self.run_blastn(outdir, input_filename, input_fasta_file, search_parameters["db_fasta"],
                                    search_parameters["db_marker"], search_parameters["identity"],
                                    search_parameters["coverage"])
                    self.get_hits(outdir, input_filename, search_parameters["db_marker"])
                    self.coverage_hits(outdir, input_filename, search_parameters["db_marker"])

                shutil.rmtree(tempdir)
                if not keep_files:
                    shutil.rmtree(new_outdir)
