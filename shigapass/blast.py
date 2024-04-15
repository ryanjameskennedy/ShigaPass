"""Module to handle blast requests"""

import os
import subprocess
import shutil
import tempfile

class Blast:
    """Class for Blast arguments"""
    def __init__(self, db, threads):
        self.db = db
        self.threads = threads

    def write_out(self, count_dict, outdir, filename, query)
        with open(os.path.join(outdir, filename, f"{query}_hits.txt"), 'w') as f:
            for gene, count in count_dict.items():
                f.write(f"{gene};{count}\n")

    def coverage_hits(self, outdir, filename, query):
        with open(os.path.join(outdir, filename, f"{query}_hits.txt"), 'r') as f1, \
                open(os.path.join(self.db, 'RFB_hits_count.csv'), 'r') as f2, \
                open(os.path.join(outdir, filename, f"{query}_hitscoverage.txt"), 'w') as f3:
            for line in f1:
                parts1 = line.strip().split(';')
                for line2 in f2:
                    parts2 = line2.strip().split(';')
                    if parts1[0] == parts2[0]:
                        coverage = float(parts1[1]) / float(parts1[2]) * 100
                        f3.write(f"{line.strip()};{coverage}\n")
                        break

    def run_blastn(self, outdir, filename, query, query, id_val, cov_val):
        blast_cmd = f"blastn -db {self.db}/{query} -query {query} -out {os.path.join(outdir, filename, f'{query}_blastout.txt')} -num_threads {self.threads} -num_alignments 10000 -outfmt 6 -word_size 11 -dust no"
        subprocess.run(blast_cmd, shell=True, check=True)
        with open(os.path.join(outdir, filename, f'{query}_blastout.txt'), 'r') as f:
            lines = f.readlines()
        with open(os.path.join(outdir, filename, f'{query}_allrecords.txt'), 'w') as f:
            for line in lines:
                parts = line.strip().split('\t')
                if float(parts[2]) >= id_val and (float(parts[3]) / float(parts[1])) * 100 >= cov_val:
                    f.write(line)

    def get_hits(self, outdir, filename, query):
        with open(os.path.join(outdir, filename, f"{query}_allrecords.txt"), 'r') as f:
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
            files = f.read().splitlines()

            for filepath in files:
                filename = os.path.splitext(os.path.basename(filepath))[0]
                tempdir = tempfile.mkdtemp(prefix="ShigaPass_")
                shutil.copy(filepath, os.path.join(tempdir, filename + "_parsed.fasta"))
                fasta_file = os.path.join(tempdir, filename + "_parsed.fasta")

                for query in ["ipaH_150-mers.fasta", "RFB_serotypes_AtoC_150-mers_v2.fasta"]:
                    self.run_blastn(outdir, filename, fasta_file, query, 98, 95)
                    self.get_hits(outdir, filename, query.split('.')[0])
                    self.coverage_hits(outdir, filename, query.split('.')[0])

                shutil.rmtree(tempdir)
