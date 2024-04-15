import os
import sys

from .__version__ import VERSION
from .cli import get_args
from .main import OptionsParser


def print_help():
    """Print help output"""
    print(f"""
                    ...::: ShigaPass v{VERSION} :::...
Author(s): Iman Yassine
          Ryan James Kennedy

Description: This tool is used to predict Shigella serotypes.

Usage: shigapass [options]
Options:
    -l, --list         List of input file(s) (FASTA) with their path(s) (mandatory)
    -o, --outdir       Output directory (mandatory)
    -p, --db           Path to databases directory (mandatory)
    -t, --threads      Number of threads (optional, default: 2)
    -u, --mkdb         Call the makeblastdb utility for databases initialisation (optional, but required when running the script for the first time)
    -k, --keep         Do not remove subdirectories (optional)
    -v, --version      Display the version and exit
    -h, --help         Display this help and exit

Example: shigapass -l list_of_fasta.txt -o ShigaPass_Results -p ShigaPass/ShigaPass_DataBases -t 4 -u -k

NOTE: The -u option should be used when running the script for the first time and after databases updates.
""")

def main():
    """Main function for initiating software"""
    args = None
    if len(sys.argv) == 1:
        print_help()
        sys.exit(0)
    elif sys.argv[1] in {"-V", "--version"}:
        print(f"ShigaPass version {VERSION}")
        sys.exit(0)
    elif sys.argv[1] in {"-h", "--h", "-help", "--help"}:
        print_help()
        sys.exit(0)
    else:
        args = get_args().parse_args()
        sp_parser = OptionsParser(VERSION)
        sp_parser.parse_options(args)
    try:
        print("Done")
    except SystemExit:
        print("Controlled exit resulting from early termination.")
        sys.exit(1)
    except KeyboardInterrupt:
        print("Controlled exit resulting from interrupt signal.")
        sys.exit(1)
    except Exception as e:
        error_message = "Uncontrolled exit resulting from an unexpected error.\n\n"
        error_message += "-" * 80 + "\n"
        error_message += f"EXCEPTION: {type(e).__name__}\n"
        error_message += f"MESSAGE: {e}\n"
        error_message += "-" * 80 + "\n\n"
        print(error_message)
        sys.exit(1)
    """
    LIST = "/path/to/list_of_fasta_files"
    OUTDIR = "/path/to/output_directory"

    # Iterate over files in a directory
    for fasta_file in os.listdir(LIST):
        if fasta_file.endswith(".fasta"):
            # Perform analysis for each FASTA file
            namedir = os.path.basename(fasta_file).replace(".fasta", "")
            outdir_namedir = os.path.join(OUTDIR, namedir)

            if not os.path.exists(outdir_namedir):
                os.makedirs(outdir_namedir)
            else:
                # Clear existing files in the directory
                files = os.listdir(outdir_namedir)
                for file in files:
                    os.remove(os.path.join(outdir_namedir, file))

            # Copy and modify the FASTA file
            with open(fasta_file, "r") as infile:
                with open(os.path.join(outdir_namedir, f"{namedir}_parsed.fasta"), "w") as outfile:
                    for record in SeqIO.parse(infile, "fasta"):
                        record.id = record.id.replace("_", "~")
                        SeqIO.write(record, outfile, "fasta")

            # Run BLAST searches
            run_blast(os.path.join(outdir_namedir, f"{namedir}_parsed.fasta"), "ipaH_150-mers.fasta", os.path.join(outdir_namedir, "ipaH_blast.txt"))
            run_blast(os.path.join(outdir_namedir, f"{namedir}_parsed.fasta"), "RFB_serotypes_AtoC_150-mers_v2.fasta", os.path.join(outdir_namedir, "rfb_blast.txt"))
            # Run other BLAST searches and analysis functions
            
            # Determine and write results to output files
            """

if __name__ == "__main__":
    main()
