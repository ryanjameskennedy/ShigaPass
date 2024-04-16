db_config_array = [
    {"db_fasta": "IPAH/ipaH_150-mers.fasta", "db_marker": "ipaH", "identity": 98, "coverage": 95},
    {"db_fasta": "RFB/RFB_serotypes_AtoC_150-mers_v2.fasta", "db_marker": "rfb", "identity": 98, "coverage": 95},
    {"db_fasta": "RFB/galF_SB1.fasta", "db_marker": "additionalrfb", "identity": 98, "coverage": 95},
    {"db_fasta": "RFB/RFB_serogroup_D_150-mers.fasta", "db_marker": "additionalrfb", "identity": 100, "coverage": 100},
    {"db_fasta": "RFB/POAC-genes_150-mers.fasta", "db_marker": "POAC", "identity": 98, "coverage": 95},
    {"db_fasta": "FLIC/fliC_Shigella_v1.fasta", "db_marker": "flic", "identity": 98, "coverage": 95},
    {"db_fasta": "CRISPR/CRISPR_spacers.fasta", "db_marker": "crispr", "identity": 100, "coverage": 100}
]

db_alt_config_array = {
    "additionalrfb": {
        "rfbs": ["A2", "A3a", "A3b", "C10"],
        "files": ["RFB_AprovBEDP02-5104_150-mers.fasta", "RFB_A16_150-mers_v2.fasta", "RFB_A16_150-mers_v2.fasta", "taurine_SB6.fasta"],
        "new_rfbs": ["AprovBEDP02-5104", "A16", "A16", "C6"],
    }
}