"""Module to handle blast requests"""

import re
import os
import subprocess
import shutil
import tempfile
import logging

from .config import db_config_array, db_alt_config_array, db_mlst_config_array, flex_score_array
from .utils import dict_to_str, write_out, read_file, get_val_from_file

LOG = logging.getLogger(__name__)

class Blast:
    """Class for Blast arguments"""
    def __init__(self, db, threads):
        self.db = db
        self.threads = threads

    def get_db_array(self, db_marker, rfb_filter=None):
        """Retrieve db array from the config"""
        return [db for db in db_config_array if db.get("db_marker") == db_marker and db.get("rfb") == rfb_filter]

    def determine_ipah(self, outdir, input_filename, input_fasta_file):
        """Determine the ipah type"""
        LOG.info("Determing ipaH  for: %s" % input_filename)
        search_parameters = self.get_db_array("ipaH")
        sorted_gene_counts, hits_coverage = self.search(outdir, input_filename, input_fasta_file, search_parameters)
        if sorted_gene_counts:
            ipah = "ipaH+"
            ipah_hits = sorted_gene_counts.get("ipaH")
            ipah_coverage = hits_coverage
        else:
            ipah = "ipaH-"
            ipah_hits = 0
            ipah_coverage = 0
        return ipah, ipah_hits, ipah_coverage

    def determine_alt_rfb(self, outdir, input_filename, input_fasta_file, rfb, rfb_hits, rfb_coverage):
        """Determine alternative rfb types"""
        for db_array in db_alt_config_array:
            if rfb == db_array["rfb"]:
                search_parameters = [db_array]
                sorted_gene_counts, hits_coverage = self.search(outdir, input_filename, input_fasta_file, search_parameters)
                if sorted_gene_counts:
                    LOG.info("Rfb has changed from %s to %s" % (rfb, db_array["new_rfb"]))
                    rfb = db_array["new_rfb"]
                    if rfb != "C6":
                        rfb_hits = sorted_gene_counts.get(rfb)
                        rfb_coverage = hits_coverage
                else:
                    rfb = db_array["rfb"]
                    LOG.info("Rfb has remained the same: %s" % rfb)
                break
        return rfb, rfb_hits, rfb_coverage

    def determine_rfb(self, outdir, input_filename, input_fasta_file):
        """Determine the rfb type"""
        LOG.info("Determining rfb type for: %s" % input_filename)
        search_parameters = self.get_db_array("rfb")
        sorted_gene_counts, hits_coverage = self.search(outdir, input_filename, input_fasta_file, search_parameters)
        rfb, rfb_hits, rfb_coverage, flexserotype, comments = [""] * 5
        if sorted_gene_counts:
            rfb = list(sorted_gene_counts.keys())[0]
            rfb_hits = sorted_gene_counts[rfb]
            rfb_coverage = hits_coverage
            rfb, rfb_hits, rfb_coverage = self.determine_alt_rfb(outdir, input_filename, input_fasta_file, rfb, rfb_hits, rfb_coverage)
            if rfb == "A3b":
                rfb_hits_fpath = os.path.join(outdir, input_filename, "rfb_hits.txt")
                rfb = get_val_from_file(rfb_hits_fpath, 0, 0)
                if rfb == "A3a":
                    rfb = "A3"
                    LOG.info("Rfb has changed to %s; hits detected are unique for %s" % (rfb, rfb))
                else:
                    rfb = "A3/A16"
                    LOG.info("Hits detected are common with A3 and Aprov97-10607")
            elif rfb == "A3a":
                rfb = "A3"
                LOG.info("Rfb has changed to %s; hits detected are unique for %s" % (rfb, rfb))
            elif rfb == "C1":
                search_parameters = self.get_db_array("additionalrfb", "C1")
                sorted_gene_counts, hits_coverage = self.search(outdir, input_filename, input_fasta_file, search_parameters)
                add_rfb_hits_fpath = os.path.join(outdir, input_filename, "additionalrfb_hits.txt")
                rfb = get_val_from_file(add_rfb_hits_fpath, 1, 0)
                if rfb:
                    rfb == "C1"
                    LOG.info("Rfb has remained %s" % rfb)
                else:
                    rfb == "C20"
                    LOG.info("Rfb has changed to %s" % rfb)
            elif rfb == "B1-5":
                LOG.info("Determining phage and plasmid encoded O-antigen modification genes")
                search_parameters = self.get_db_array("POAC")
                sorted_gene_counts, hits_coverage = self.search(outdir, input_filename, input_fasta_file, search_parameters)
                flexserotype = self.get_poac(outdir, input_filename)
            elif not rfb:
                search_parameters = self.get_db_array("additionalrfb", "C1")
                sorted_gene_counts, hits_coverage = self.search(outdir, input_filename, input_fasta_file, search_parameters)
                add_rfb_hits_fpath = os.path.join(outdir, input_filename, "additionalrfb_hitscoverage.txt")
                rfb_hits = get_val_from_file(add_rfb_hits_fpath, 3, 1)
                rfb_coverage = get_val_from_file(add_rfb_hits_fpath, 3, 3)
                if rfb_hits >= 30:
                    rfb = "D"
                    LOG.info("Rfb has changed to %s" % rfb)
                else:
                    rfb_hits = 0
                    rfb_coverage = 0
            else:
                LOG.info("Rfb blast completed")
            comments = self.check_multiple_rfbs(outdir, input_filename, "rfb", rfb)
        return rfb, rfb_hits, rfb_coverage, flexserotype, comments

    def determine_mlst(self, outdir, input_filename, input_fasta_file):
        """Determine MLST based on BLAST results"""
        LOG.info("Determining MLST type for: %s" % input_filename)
        for mlst_array in db_mlst_config_array:
            gene = mlst_array["db_marker"]
            self.run_blastn(outdir, input_filename, input_fasta_file,
                             mlst_array["db_fasta"], gene,
                             mlst_array["identity"], mlst_array["coverage"])
            mlst_blastfpath = os.path.join(outdir, input_filename, f"{gene}_allrecords.txt")
            mlst_blast = read_file(mlst_blastfpath)
            try:
                mlst = mlst_blast.split("\t")[1].split(":")[0]
                gene_st = f"ST_${mlst}"
            except KeyError:
                gene_st = f"ST_${gene}"
            output = f"{gene}:{gene_st}\n"
            mlst_outfpath = os.path.join(outdir, input_filename, "mlst_alleles.txt")
            write_out(output, mlst_outfpath, "a")

    def determine_flic(self):
        # Implement logic to determine fliC type
        pass

    def determine_crispr(self):
        # Implement logic to determine CRISPR type
        pass

    def combine_data(self):
        # Combine the results to infer the serotype
        pass

    def check_multiple_rfbs(self, outdir, input_filename, db_marker, rfb):
        comments = ""
        hitscov_fpath = os.path.join(outdir, input_filename, f"{db_marker}_hitscoverage.txt")
        content = read_file(hitscov_fpath)
        if len(content) >= 3:
            comments = f"More than one rfb is detected: {len(content)}"
            LOG.info("Mutliple rfbs")
        elif len(content) == 2:
            if rfb == "AprovBEDP02-5104" or rfb == "A16" or  rfb == "A3":
                comments = ""
            else:
                comments = f"More than one rfb is detected: {len(content)}"
        return comments

    def get_poac(self, outdir, input_filename):
        # Read and modify the file content
        poac_hits_fpath = os.path.join(outdir, input_filename, "POAC_hits.txt")
        content = read_file(poac_hits_fpath)

        # Perform the replace using re
        replacements = {
            "gtrX": "32", "gtrII": "4", "gtrIC": "2", "gtrIV": "8", "gtrV": "16",
            "gtrI": "1", "oac1b": "128", "oac": "64", "optII": "256"
        }
        for pattern, repl in replacements.items():
            content = re.sub(pattern, repl, content)
        write_out(content, poac_hits_fpath)
        # Get score
        lines = content.rstrip().split("\n")
        score = sum(int(line.split(";")[0]) for line in lines)
        score_outfpath = os.path.join(outdir, input_filename, "score.txt")
        write_out(f"score={score}\n", score_outfpath)
        LOG.info("This is the score: %s" % score)

        # Find the corresponding flexserotype
        flexserotype = "Unknown"
        for flex_array in flex_score_array:
            if score == flex_array["score"]:
                flexserotype = flex_array["flex"]
                break
        LOG.info("This is the flexserotype: %s" % flexserotype)

        # Sort, and concatenate phages
        phages = sorted(set(line.split(";")[0] for line in lines))
        phages_str = ";".join(phages)

        # Write the summary to a CSV file
        summary_outfpath = os.path.join(outdir, "ShigaPass_Flex_summary.csv")
        write_out(f"{input_filename};{phages_str};{flexserotype}\n", summary_outfpath, "a")
        return flexserotype

    def get_hits(self, outdir, input_filename, db_marker):
        gene_counts = {}
        filtered_records_fpath = os.path.join(outdir, input_filename, f"{db_marker}_allrecords.txt")
        with open(filtered_records_fpath, "r") as fin:
            for line in fin:
                parts = line.strip().split("\t")
                gene = parts[1].split("_")[0]
                gene_counts[gene] = gene_counts.get(gene, 0) + 1
        sorted_gene_counts = dict(sorted(gene_counts.items(), key=lambda item: item[1], reverse=True))
        outfpath = os.path.join(outdir, input_filename, f"{db_marker}_hits.txt")
        content = dict_to_str(sorted_gene_counts)
        write_out(content, outfpath)
        LOG.info("The gene hits output was written to: %s" % outfpath)
        return sorted_gene_counts

    def coverage_hits(self, outdir, input_filename, db_marker):
        hits_fpath = os.path.join(outdir, input_filename, f"{db_marker}_hits.txt")
        rfb_hits_count_fpath = os.path.join(self.db, "RFB_hits_count.csv")
        hitscov_fpath = os.path.join(outdir, input_filename, f"{db_marker}_hitscoverage.txt")
        hits_coverage = ""
        with open(hits_fpath, "r") as fin1, \
            open(rfb_hits_count_fpath, "r") as fin2, \
            open(hitscov_fpath, "w") as fout:
            for line in fin1:
                parts1 = line.rstrip().split(";")
                for line2 in fin2:
                    parts2 = line2.rstrip().split(";")
                    if parts1[0] == parts2[0]:
                        hits_coverage = float(parts1[1]) / float(parts2[1]) * 100
                        fout.write(f"{line.rstrip()};{parts2[1]};{hits_coverage}\n")
                        break
            if not hits_coverage:
                fout.write("")
            LOG.info("The hits coverage was written to: %s" % hitscov_fpath)
        return hits_coverage

    def run_blastn(self, outdir, input_filename, input_fasta_file, db_fasta, db_marker, identity, coverage):
        blast_cmd = f"blastn -db {self.db}/{db_fasta} -query {input_fasta_file} -out {os.path.join(outdir, input_filename, f"{db_marker}_blastout.txt")} -num_threads {self.threads} -num_alignments 10000 -outfmt 6 -word_size 11 -dust no"
        subprocess.run(blast_cmd, shell=True, check=True)
        blast_outfpath = os.path.join(outdir, input_filename, f"{db_marker}_blastout.txt")
        filtered_blast_outfpath = os.path.join(outdir, input_filename, f"{db_marker}_allrecords.txt")
        self.filter_blast(blast_outfpath, filtered_blast_outfpath, identity, coverage)

    def filter_blast(self, blast_outfpath, filtered_blast_outfpath, identity, coverage):
        with open(blast_outfpath, "r") as blastfile, \
            open(filtered_blast_outfpath, "w") as fout:
            for line in blastfile:
                parts = line.rstrip().split("\t")
                length = float(parts[3])
                try:
                    q_len = float(parts[1].split("_")[-1])
                except ValueError:
                    q_len = float(parts[1].split(":")[-1])
                if float(parts[2]) >= identity and (length / q_len) * 100 >= coverage:
                    fout.write(line)
        LOG.info("The blast filtered output was written to: %s" % filtered_blast_outfpath)

    def search(self, outdir, input_filename, input_fasta_file, db_array):
        search_parameters = db_array[0]
        self.run_blastn(outdir, input_filename, input_fasta_file, search_parameters["db_fasta"],
            search_parameters["db_marker"], search_parameters["identity"],
            search_parameters["coverage"])
        sorted_gene_counts = self.get_hits(outdir, input_filename, search_parameters["db_marker"])
        hits_coverage = self.coverage_hits(outdir, input_filename, search_parameters["db_marker"])
        return sorted_gene_counts, hits_coverage

    def make_blast_db(self):
        for root, dirs, files in os.walk(self.db):
            for filename in files:
                if filename.endswith(".fasta"):
                    subprocess.run(f"makeblastdb -dbtype nucl -in {os.path.join(root, filename)}", shell=True, check=True)

    def run(self, list_file, outdir, keep_files):
        with open(list_file, "r") as f:
            input_files = f.read().splitlines()

            for filepath in input_files:
                input_filename = os.path.splitext(os.path.basename(filepath))[0]
                tempdir = tempfile.mkdtemp(prefix="ShigaPass_")
                input_fasta_file = os.path.join(tempdir, input_filename + "_parsed.fasta")
                shutil.copy(filepath, input_fasta_file)
                new_outdir = os.path.join(outdir, input_filename)
                if not os.path.exists(new_outdir):
                    os.makedirs(new_outdir)

                ipah, ipah_hits, ipah_coverage = self.determine_ipah(outdir, input_filename, input_fasta_file)
                if ipah == "ipaH-":
                    serotype="Not Shigella/EIEC"
                    rfb="ND"
                    mlst="ND"
                    flic="ND"
                    crispr="ND"
                    hit="ND"
                    RFB_coverage=0
                else:
                    rfb, rfb_hits, rfb_coverage, flexserotype, comments = self.determine_rfb(outdir, input_filename, input_fasta_file)
                    self.determine_mlst(outdir, input_filename, input_fasta_file)

                shutil.rmtree(tempdir)
                if not keep_files:
                    shutil.rmtree(new_outdir)
