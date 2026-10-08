# /// script
# requires-python = "==3.12.*"
# dependencies = [
#     "nmd-scanner==0.5.0",
#     "pandas==3.0.6",
#     "numpy==2.5.3",
#     "scipy==1.18.1",
#     "pyarrow==24.0.0",
#     "pyfaidx==0.9.0.4",
#     "onnxruntime==1.30.0",
#     "openpyxl==3.1.5",
# ]
# [tool.uv]
# exclude-newer = "2026-10-07T16:00:00Z"
# exclude-newer-package = { nmd-scanner = "2026-10-07T20:00:00Z" }
# ///
"""
Compute the numbers of the section "Model validation" of the main README: how well nmd_efficiency_rf.onnx and
nmd_escape alone rank the measured NMD efficiency in 4 cohorts.

Run it with:

    uv run scripts/validate_model.py --gff3 gencode.v42.annotation.gff3.gz --fasta GRCh38.fa --train-dir train/ \
        --out-dir out/

--train-dir is the --out-dir of the run of scripts/train_model.py that trained nmd_efficiency_rf.onnx.

The script pins nmd-scanner 0.5.0. Version 0.5.0 computes the same features as 0.4.0, on which the model was trained.
PyPI got 0.5.0 after the date of exclude-newer, so exclude-newer-package sets a later date for nmd-scanner alone.

Cohorts:
- TCGA: the out-of-fold predictions in oof_predictions.parquet of --train-dir.
- Geuvadis: germline nonsense variants in lymphoblastoid cell lines (Iha et al. 2025), Zenodo record 16666299. Their
  NMD efficiency is -log2(VAF_mean / 0.5). Only the variants of the study's accurate annotation count: no rescued MNV,
  no nontranslating gene, no mixed NMD isoforms.
- GTEx: germline variants in the tissues of GTEx v8, Supplementary Data 2 of Kim et al. 2024. Their NMD efficiency is
  -log2(VAF_RNA / 0.5). The table has one row per variant and tissue (variant-tissue rows). The validation uses the
  mean of each variant over its tissues.
- MMRF-TARGET: somatic variants in tumors, MMRF_TARGET_dataset.csv of NMDEff (https://github.com/hjkng/nmdeff).

Steps:
1. Download the tables of Geuvadis, GTEx and MMRF-TARGET and check their sha256. Zenodo has the Geuvadis tables only
   inside one zip of 3.2 GB, so fetch just the 2 members by HTTP range requests. With --geuvadis-dir, read the 2
   members from a local directory instead, with the same check.
2. GTEx and MMRF-TARGET give the alleles in HGVSc on the transcript strand. Convert them to the forward strand with
   add_alleles of make_benchmark_vcfs.py. Geuvadis gives genomic alleles. Fail if a REF does not match the FASTA.
3. Run nmd_scanner.annotate on the unique variants of the 3 cohorts, and score its rows with the ONNX model. Join each
   measured row to the NMD-Scanner row of its transcript and variant.
4. Keep the rows with nmd_model_status "ok" and a variant allele fraction above 0 in the RNA.
5. Per cohort, compute the Spearman correlation of the predictions and of (not nmd_escape) with the measured NMD
   efficiency. The 95% CIs are percentiles of N_RESAMPLES bootstrap resamples of the rows. The gain, model minus
   rules, has a paired CI from the same resamples. Each cohort draws from its own generator, seeded with SEED.
6. Per cohort, compute the mean measured NMD efficiency of the rows with nmd_escape False. For GTEx, also compute it
   over the variant-tissue rows, in total and per tissue with at least MIN_TISSUE_ROWS of them.

Outputs in --out-dir: variants.vcf (the input of nmd_scanner.annotate) and validation.json.
"""

import argparse
import io
import json
import logging
import urllib.request
import zipfile
from pathlib import Path

import numpy as np
import onnxruntime
import pandas as pd
from make_benchmark_vcfs import KEY, add_alleles, add_key, check_sha256, download, read_nmdeff, read_strands, write_vcf
from pyfaidx import Fasta
from scipy.stats import rankdata

import nmd_scanner
from nmd_scanner.schema import MODEL_INPUTS

logger = logging.getLogger("validate_model")

ONNX = Path(__file__).resolve().parent.parent / "nmd_efficiency_rf.onnx"
# Supplementary Data 1 to 12 of Kim et al. 2024 (Commun Biol 7:1461), one sheet each
GTEX_URL = (
    "https://static-content.springer.com/esm/art%3A10.1038%2Fs42003-024-07136-y/MediaObjects/"
    "42003_2024_7136_MOESM3_ESM.xlsx"
)
GTEX_SHA256 = "637c92b17095ef6147b74ae0ab6549f674adc6a67ee04c06d483fa7331ff34c6"
GTEX_SHEET = "Supplementary Data 2"
# Zenodo record 16666299 of Iha et al. 2025 has one file, NMD-rules.zip, of ZENODO_SIZE bytes
ZENODO_URL = "https://zenodo.org/api/records/16666299/files/NMD-rules.zip/content"
ZENODO_SIZE = 3_189_701_188
NONSENSE_LIST = "NMD-rules/Datasets/Nonsense_list.txt"
ISOFORMS = "NMD-rules/Datasets/IsoformLevel.NMDannotation.txt"
# zipfile checks the CRC32 of each member, and read_member checks the sha256
GEUVADIS_SHA256 = {
    NONSENSE_LIST: "5a248be88a6b2a287c2baf12e542d638257a06daadc21a4c440af812518bce26",
    ISOFORMS: "bc16f2b869149f07e0a5d7ce4173df4fbca9ea9c4fab2dd9f3c07a82c4011a7d",
}
# A Geuvadis variant with one of these flags in IsoformLevel.NMDannotation.txt is not in the accurate annotation
ACCURATE_ANNOTATION_FLAGS = ["resMNV", "nonTranslating", "Multiisoforms"]
# VAF_DNA of a heterozygous germline variant
GERMLINE_VAF_DNA = 0.5

SEED = 0
N_RESAMPLES = 2000
MIN_TISSUE_ROWS = 100


class RemoteFile:
    """
    A read-only file at a URL that answers HTTP range requests. Each read is one request, so zipfile can read single
    members of a remote zip.
    """

    def __init__(self, url: str, size: int) -> None:
        self.url = url
        self.size = size
        self.position = 0

    def seekable(self) -> bool:
        return True

    def tell(self) -> int:
        return self.position

    def seek(self, offset: int, whence: int = io.SEEK_SET) -> int:
        self.position = offset + {io.SEEK_SET: 0, io.SEEK_CUR: self.position, io.SEEK_END: self.size}[whence]
        return self.position

    def read(self, n: int = -1) -> bytes:
        end = self.size if n < 0 else min(self.position + n, self.size)
        if end <= self.position:
            return b""
        request = urllib.request.Request(self.url, headers={"Range": f"bytes={self.position}-{end - 1}"})
        with urllib.request.urlopen(request, timeout=120) as response:
            # A server that ignores the range answers 200 and sends the whole file
            if response.status != 206:
                raise ValueError(f"{self.url} answered a range request with HTTP {response.status}")
            data = response.read()
        self.position += len(data)
        return data


def read_member(member: str, geuvadis_dir: Path | None) -> io.BytesIO:
    """Read a member of the Zenodo zip from geuvadis_dir, or else by HTTP range requests, and check its sha256."""
    if geuvadis_dir is None:
        with zipfile.ZipFile(RemoteFile(ZENODO_URL, ZENODO_SIZE)) as archive:
            data = archive.read(member)
    else:
        data = (geuvadis_dir / Path(member).name).read_bytes()
    return io.BytesIO(check_sha256(data, GEUVADIS_SHA256[member], member))


def read_geuvadis(geuvadis_dir: Path | None, fasta: Fasta) -> pd.DataFrame:
    """
    Read the nonsense SNVs of Nonsense_list.txt, one row per variant and MANE Select transcript. Add the flags of
    IsoformLevel.NMDannotation.txt, which has the evaluated variants. REF and ALT are genomic. Fail if a REF does not
    match the FASTA.
    """
    table = pd.read_csv(read_member(NONSENSE_LIST, geuvadis_dir), sep="\t").rename(
        columns={"CHROM": "chrom", "POS": "pos", "REF": "ref", "ALT": "alt", "TranscriptID": "transcript"}
    )
    fasta_ref = [
        fasta[chrom][pos - 1 : pos - 1 + len(ref)].seq.upper()
        for chrom, pos, ref in zip(table["chrom"], table["pos"], table["ref"], strict=True)
    ]
    if (table["ref"] != fasta_ref).any():
        raise ValueError(f"rows of {NONSENSE_LIST} have a REF that does not match the FASTA")
    flags = pd.read_csv(
        read_member(ISOFORMS, geuvadis_dir), sep="\t", usecols=["Variant_ID", *ACCURATE_ANNOTATION_FLAGS]
    )
    return table.merge(flags, on="Variant_ID", how="left", validate="many_to_one")


def score(features: pd.DataFrame) -> pd.DataFrame:
    """
    Return KEY, nmd_model_status and nmd_escape of each NMD-Scanner row, and the prediction of the model for the rows
    with nmd_model_status "ok". Fail if the model does not take MODEL_INPUTS.
    """
    session = onnxruntime.InferenceSession(ONNX, providers=["CPUExecutionProvider"])
    if json.loads(session.get_modelmeta().custom_metadata_map["feature_names"]) != list(MODEL_INPUTS):
        raise ValueError(f"{ONNX} does not take MODEL_INPUTS")
    ok = (features["nmd_model_status"] == "ok").to_numpy()
    prediction = np.full(len(features), np.nan)
    prediction[ok] = session.run(None, {"input": features[MODEL_INPUTS][ok].astype("float64").to_numpy()})[0].ravel()
    return features[[*KEY, "nmd_model_status", "nmd_escape"]].assign(prediction=prediction)


def scored_rows(measured: pd.DataFrame, scored: pd.DataFrame, vaf_rna: str) -> pd.DataFrame:
    """
    Join each measured row to the NMD-Scanner row of its transcript and variant. Keep the rows with nmd_model_status
    "ok" and a variant allele fraction above 0 in the RNA, which is the column vaf_rna.
    """
    rows = measured.merge(scored, on=KEY, how="inner", validate="many_to_one")
    return rows[(rows["nmd_model_status"] == "ok") & (rows[vaf_rna] > 0)].astype({"nmd_escape": bool})


def spearman(x: np.ndarray, y: np.ndarray) -> float:
    return float(np.corrcoef(rankdata(x), rankdata(y))[0, 1])


def percentiles(values: np.ndarray) -> list[float]:
    return [float(v) for v in np.percentile(values, [2.5, 97.5])]


def correlations(rows: pd.DataFrame) -> dict:
    """
    Spearman correlation of the prediction and of the rules with NMD_efficiency, and their difference, with bootstrap
    CIs. The rules score is (not nmd_escape), so that both scores rank high NMD efficiency first.
    """
    y = rows["NMD_efficiency"].to_numpy(dtype="float64")
    model = rows["prediction"].to_numpy(dtype="float64")
    rules = (~rows["nmd_escape"]).to_numpy(dtype="float64")
    rng = np.random.default_rng(SEED)
    model_resampled = np.empty(N_RESAMPLES)
    rules_resampled = np.empty(N_RESAMPLES)
    for resample in range(N_RESAMPLES):
        i = rng.integers(0, len(y), len(y))
        model_resampled[resample] = spearman(model[i], y[i])
        rules_resampled[resample] = spearman(rules[i], y[i])
    model_value = spearman(model, y)
    rules_value = spearman(rules, y)
    return {
        "n": len(y),
        "model": model_value,
        "model_ci": percentiles(model_resampled),
        "rules": rules_value,
        "rules_ci": percentiles(rules_resampled),
        "gain": model_value - rules_value,
        "gain_ci": percentiles(model_resampled - rules_resampled),
    }


def mean_nmd_triggering(rows: pd.DataFrame) -> float:
    """The mean NMD_efficiency of the rows with nmd_escape False."""
    return float(rows.loc[~rows["nmd_escape"], "NMD_efficiency"].mean())


def estimate(value: float, ci: list[float]) -> str:
    return f"{value:.2f} ({ci[0]:.2f} to {ci[1]:.2f})"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--gff3", required=True, type=Path, help="GENCODE GFF3 of GRCh38, optionally gzip-compressed")
    parser.add_argument("--fasta", required=True, type=Path, help="GRCh38 FASTA with chr names and a .fai index")
    parser.add_argument(
        "--train-dir",
        required=True,
        type=Path,
        help="--out-dir of the train_model.py run that trained the model; it has the out-of-fold predictions",
    )
    parser.add_argument("--out-dir", required=True, type=Path, help="directory for the outputs")
    parser.add_argument(
        "--geuvadis-dir",
        type=Path,
        help="local directory with Nonsense_list.txt and IsoformLevel.NMDannotation.txt of NMD-rules.zip; skips the"
        " range requests to Zenodo. The sha256 of both files is still checked.",
    )
    args = parser.parse_args()
    logging.basicConfig(level=logging.WARNING, format="%(asctime)s %(name)s %(levelname)s %(message)s")
    logger.setLevel(logging.INFO)
    args.out_dir.mkdir(parents=True, exist_ok=True)
    trained = args.train_dir / "nmd_efficiency_rf.onnx"
    if trained.read_bytes() != ONNX.read_bytes():
        raise ValueError(f"{trained} differs from {ONNX}, so its out-of-fold predictions are of another model")

    strands = read_strands(args.gff3)
    fasta = Fasta(str(args.fasta))
    geuvadis = read_geuvadis(args.geuvadis_dir, fasta)
    # The first line of the sheet is its title
    gtex = pd.read_excel(io.BytesIO(download(GTEX_URL, GTEX_SHA256)), sheet_name=GTEX_SHEET, header=1)
    gtex = add_alleles(gtex, strands, fasta).rename(columns={"Tissue_type": "tissue"})
    mmrf = add_alleles(read_nmdeff("MMRF_TARGET_dataset.csv"), strands, fasta)
    vcf = args.out_dir / "variants.vcf"
    write_vcf(
        pd.concat([table[["chrom", "pos", "ref", "alt"]] for table in (geuvadis, gtex, mmrf)]).drop_duplicates(), vcf
    )
    scored = score(add_key(nmd_scanner.annotate(vcf, args.gff3, args.fasta, sequences=False)))

    geuvadis = scored_rows(geuvadis, scored, "VAF_mean")
    accurate = ~geuvadis[ACCURATE_ANNOTATION_FLAGS].astype("boolean").fillna(False).any(axis=1)
    geuvadis = geuvadis[accurate].assign(NMD_efficiency=lambda rows: -np.log2(rows["VAF_mean"] / GERMLINE_VAF_DNA))
    # ASE_altered is VAF_RNA. For VAF_RNA 0, the table caps NMD_efficiency at 10.
    gtex = scored_rows(gtex, scored, "ASE_altered")
    gtex_variants = gtex.groupby(KEY).agg(
        NMD_efficiency=("NMD_efficiency", "mean"),
        prediction=("prediction", "first"),
        nmd_escape=("nmd_escape", "first"),
    )
    cohorts = {
        # Every row of tcga_dataset.csv has VAF_RNA above 0
        "tcga": pd.read_parquet(args.train_dir / "oof_predictions.parquet").astype({"nmd_escape": bool}),
        "geuvadis": geuvadis,
        "gtex": gtex_variants,
        "mmrf_target": scored_rows(mmrf, scored, "VAF_RNA"),
    }

    validation = {}
    for cohort, rows in cohorts.items():
        result = {**correlations(rows), "mean_nmd_triggering": mean_nmd_triggering(rows)}
        logger.info(
            "%s, %d rows: model %s, rules %s, gain %s, NMD-triggering mean %.2f",
            cohort,
            result["n"],
            estimate(result["model"], result["model_ci"]),
            estimate(result["rules"], result["rules_ci"]),
            estimate(result["gain"], result["gain_ci"]),
            result["mean_nmd_triggering"],
        )
        validation[cohort] = result
    counts = gtex["tissue"].value_counts()
    tissue_means = {
        tissue: mean_nmd_triggering(gtex[gtex["tissue"] == tissue])
        for tissue in counts.index[counts >= MIN_TISSUE_ROWS]
    }
    tissue_means = dict(sorted(tissue_means.items(), key=lambda item: item[1]))
    validation["gtex_variant_tissue_rows"] = {
        "n": len(gtex),
        "tissues": gtex["tissue"].nunique(),
        "mean_nmd_triggering": mean_nmd_triggering(gtex),
        "min_rows_per_tissue": MIN_TISSUE_ROWS,
        "mean_nmd_triggering_per_tissue": tissue_means,
    }
    logger.info(
        "gtex, %d variant-tissue rows: NMD-triggering mean %.2f, from %.2f to %.2f between the %d tissues with at"
        " least %d rows",
        len(gtex),
        mean_nmd_triggering(gtex),
        min(tissue_means.values()),
        max(tissue_means.values()),
        len(tissue_means),
        MIN_TISSUE_ROWS,
    )
    (args.out_dir / "validation.json").write_text(json.dumps(validation, indent=2) + "\n")


if __name__ == "__main__":
    main()
