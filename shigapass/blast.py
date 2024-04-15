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

    def write_out(self, count_dict, outdir, filename, query)
        with open(os.path.join(outdir, filename, f"{query}_hits.txt"), 'w') as f:
            for gene, count in count_dict.items():
                f.write(f"{gene};{count}\n")

    def coverage_hits(self, outdir, input_filename, query):
        with open(os.path.join(outdir, input_filename, f"{query}_hits.txt"), 'r') as f1, \
            open(os.path.join(self.db, 'RFB_hits_count.csv'), 'r') as f2, \
            open(os.path.join(outdir, input_filename, f"{query}_hitscoverage.txt"), 'w') as f3:
            for line in f1:
                parts1 = line.strip().split(';')
                for line2 in f2:
                    parts2 = line2.strip().split(';')
                    if parts1[0] == parts2[0]:
                        coverage = float(parts1[1]) / float(parts1[2]) * 100
                        f3.write(f"{line.strip()};{coverage}\n")
                        break

    def run_blastn(self, outdir, input_filename, input_fasta_file, db_fasta, db_marker, identity, coverage):
        blast_cmd = f"blastn -db {self.db}/{db_fasta} -query {input_fasta_file} -out {os.path.join(outdir, input_filename, f'{db_marker}_blastout.txt')} -num_threads {self.threads} -num_alignments 10000 -outfmt 6 -word_size 11 -dust no"
        subprocess.run(blast_cmd, shell=True, check=True)
        self.filter_blast(identity, coverage)

    def filter_blast(self, outdir, input_filename, db_marker, identity, coverage):
        with open(os.path.join(outdir, input_filename, f'{db_marker}_blastout.txt'), 'r') as blastfile, \
            open(os.path.join(outdir, input_filename, f'{db_marker}_allrecords.txt'), 'w') as fout:
            for line in blastfile:
                parts = line.rstrip().split('\t')
                if float(parts[2]) >= identity and (float(parts[3]) / float(parts[1])) * 100 >= coverage:
                    fout.write(line)

    def get_hits(self, outdir, input_filename, db_marker):
        with open(os.path.join(outdir, input_filename, f"{db_marker}_allrecords.txt"), 'r') as f:
            lines = f.readlines()
            gene_counts = {}
            for line in lines:
                parts = line.strip().split('\t')
                gene = parts[1].split('_')[0]
                gene_counts[gene] = gene_counts.get(gene, 0) + 1
    
    def make_blast_db(self):
        for root, dirs, files in os.walk(self.db):
            for filename in files:
                if filename.endswith('.fasta'):
                    subprocess.run(f"makeblastdb -dbtype nucl -in {os.path.join(root, filename)}", shell=True, check=True)

    def run(self, list_file, outdir):
        with open(list_file, 'r') as f:
            input_files = f.read().splitlines()

            for filepath in input_files:
                input_filename = os.path.splitext(os.path.basename(filepath))[0]
                tempdir = tempfile.mkdtemp(prefix="ShigaPass_")
                input_fasta_file = os.path.join(tempdir, input_filename + "_parsed.fasta")
                shutil.copy(filepath, input_fasta_file)

                for search_parameters in db_config_array:
                    self.run_blastn(outdir, input_filename, input_fasta_file, search_parameters["db_fasta"],
                                    search_parameters["db_marker"], search_parameters["identity"],
                                    search_parameters["coverage"])
                    self.get_hits(outdir, input_filename, search_parameters["db_marker"])
                    self.coverage_hits(outdir, input_filename, search_parameters["db_marker"])

                shutil.rmtree(tempdir)
