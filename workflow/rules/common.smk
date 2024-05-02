import pandas as pd
import os.path
import glob
from pathlib import Path

#### Defining global variables ####
BAMS = Path(config["data"])
REFS = Path(config["references"])
OUTDIR = Path(config["outdir"])
SAMPLES_DIR= OUTDIR / "samples"
REFS_DIR = OUTDIR / "references"
DATASET_DIR = OUTDIR / "dataset"

CHROM_NAMES = config["chrom_names"]

SAMPLEFILE=config["sample_table"]
SAMPLETABLE=(pd.read_csv(config["sample_table"], sep=","))
SAMPLES=list(set(SAMPLETABLE["sample"]))


d={'sample': SAMPLETABLE["sample"],
    'group': SAMPLETABLE["group"]}

SAMPLE_REFERENCE = pd.DataFrame(data=d).set_index("sample", drop=False)

#### Defining which final output files are being requested ####
def get_final_output():
    final_output = expand(SAMPLES_DIR / "cnv" / "{sample}" / "copy_number_variants.tsv",sample=SAMPLES)
    final_output.append(DATASET_DIR / "cnv" / "copy_number_variants_dataset.tsv")
    return final_output
