"""
GTEx: germline nonsense variants in the tissues of GTEx v8, Supplementary Data 2 of Kim et al. 2024 (Commun Biol
7:1461).

The table has one row per variant and tissue. NMD_efficiency is its column NMD_efficiency, and ASE_altered is the
variant allele fraction in the RNA. For ASE_altered 0, the table caps NMD_efficiency at 10, so build keeps the rows with
ASE_altered above 0. In 4,102 of the 17,816 rows with ASE_altered above 0, NMD_efficiency differs from
-log2(ASE_altered / 0.5). The difference moves the Spearman correlations of scripts/validate_model.py.py by only about
0.001.
"""

import logging
from pathlib import Path

import polars as pl
from pyfaidx import Fasta

from nmd_scanner.benchmark.sources import fetch, finish_table, hgvsc_alleles

logger = logging.getLogger(__name__)

NAME = "GTEx"
# Supplementary Data 1 to 12 of Kim et al. 2024, one sheet each
URL = (
    "https://static-content.springer.com/esm/art%3A10.1038%2Fs42003-024-07136-y/MediaObjects/"
    "42003_2024_7136_MOESM3_ESM.xlsx"
)
SHA256 = "637c92b17095ef6147b74ae0ab6549f674adc6a67ee04c06d483fa7331ff34c6"
SHEET = "Supplementary Data 2"
# The GTEx v8 names (SMTSD) of the abbreviations in the column Tissue_type. AbExp uses the same names.
TISSUE_NAMES = {
    "ADPSBQ": "Adipose - Subcutaneous",
    "ADPVSC": "Adipose - Visceral (Omentum)",
    "ADRNLG": "Adrenal Gland",
    "ARTAORT": "Artery - Aorta",
    "ARTCRN": "Artery - Coronary",
    "ARTTBL": "Artery - Tibial",
    "BLDDER": "Bladder",
    "BREAST": "Breast - Mammary Tissue",
    "BRNACC": "Brain - Anterior cingulate cortex (BA24)",
    "BRNAMY": "Brain - Amygdala",
    "BRNCDT": "Brain - Caudate (basal ganglia)",
    # As in GTEx. Supplementary Data 9 of Kim et al. 2024 swaps the names of BRNCHA and BRNCHB.
    "BRNCHA": "Brain - Cerebellar Hemisphere",
    "BRNCHB": "Brain - Cerebellum",
    "BRNCTXA": "Brain - Cortex",
    "BRNCTXB": "Brain - Frontal Cortex (BA9)",
    "BRNHPP": "Brain - Hippocampus",
    "BRNHPT": "Brain - Hypothalamus",
    "BRNNCC": "Brain - Nucleus accumbens (basal ganglia)",
    "BRNPTM": "Brain - Putamen (basal ganglia)",
    "BRNSNG": "Brain - Substantia nigra",
    "BRNSPC": "Brain - Spinal cord (cervical c-1)",
    "CLNSGM": "Colon - Sigmoid",
    "CLNTRN": "Colon - Transverse",
    "CVSEND": "Cervix - Endocervix",
    "CVXECT": "Cervix - Ectocervix",
    "ESPGEJ": "Esophagus - Gastroesophageal Junction",
    "ESPMCS": "Esophagus - Mucosa",
    "ESPMSL": "Esophagus - Muscularis",
    "FIBRBLS": "Cells - Cultured fibroblasts",
    "FLLPNT": "Fallopian Tube",
    "HRTAA": "Heart - Atrial Appendage",
    "HRTLV": "Heart - Left Ventricle",
    "KDNCTX": "Kidney - Cortex",
    "KDNMDL": "Kidney - Medulla",
    "LCL": "Cells - EBV-transformed lymphocytes",
    "LIVER": "Liver",
    "LUNG": "Lung",
    "MSCLSK": "Muscle - Skeletal",
    "NERVET": "Nerve - Tibial",
    "OVARY": "Ovary",
    "PNCREAS": "Pancreas",
    "PRSTTE": "Prostate",
    "PTTARY": "Pituitary",
    "SKINNS": "Skin - Not Sun Exposed (Suprapubic)",
    "SKINS": "Skin - Sun Exposed (Lower leg)",
    "SLVRYG": "Minor Salivary Gland",
    "SNTTRM": "Small Intestine - Terminal Ileum",
    "SPLEEN": "Spleen",
    "STMACH": "Stomach",
    "TESTIS": "Testis",
    "THYROID": "Thyroid",
    "UTERUS": "Uterus",
    "VAGINA": "Vagina",
    "WHLBLD": "Whole Blood",
}


def download(downloads: Path) -> pl.DataFrame:
    """Read the sheet of Supplementary Data 2 with fetch. polars reads it with fastexcel, of the extra scripts."""
    # The first line of the sheet is its title
    return pl.read_excel(fetch(downloads, URL, SHA256), sheet_name=SHEET, read_options={"header_row": 1})


def build(raw_df: pl.DataFrame, transcripts_df: pl.DataFrame, fasta: Fasta) -> pl.DataFrame:
    """
    The benchmark table of raw_df, the table of download, with the alleles of hgvsc_alleles and the column tissue. It has
    the rows with ASE_altered above 0. tissue is the name of Tissue_type in TISSUE_NAMES. build raises an
    InvalidOperationError at an abbreviation that TISSUE_NAMES lacks.
    """
    table_df = raw_df.with_row_index("source_row").filter(pl.col("ASE_altered") > 0)
    logger.info("%s: %d of %d rows lack ASE_altered or have 0", NAME, len(raw_df) - len(table_df), len(raw_df))
    table_df = table_df.with_columns(tissue=pl.col("Tissue_type").replace_strict(TISSUE_NAMES))
    return finish_table(hgvsc_alleles(table_df, transcripts_df, fasta), transcripts_df, fasta)
