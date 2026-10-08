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
Validate the NMD efficiency model nmd_efficiency_rf.onnx on 4 cohorts, and compute the table of the section "Model
validation" of the main README.

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
1. Download the tables of Geuvadis, GTEx and MMRF-TARGET to <out-dir>/downloads/ and check their sha256. Zenodo has
   the Geuvadis tables only inside one zip of 3.2 GB, so fetch just the 2 members by HTTP range requests. With
   --geuvadis-dir, read the 2 members from a local directory instead, with the same check.
2. GTEx and MMRF-TARGET give the alleles in HGVSc on the transcript strand. Convert them to the forward strand with
   add_alleles of scripts/make_benchmark_vcfs.py. Geuvadis gives genomic alleles. Fail if a REF does not match the
   FASTA. Write the unique variants of the 3 cohorts to <out-dir>/variants.vcf.
3. Run nmd_scanner.annotate on that VCF. It runs add_features_and_rules. Join each measured row to the NMD-Scanner
   row of its transcript and variant, on the transcript ID without version, chrom, pos, ref and alt. NMD-Scanner's
   start is 0-based, so start + 1 is pos.
4. Keep the rows with nmd_model_status "ok" and a variant allele fraction above 0 in the RNA. Score them with the
   ONNX model.
5. Per cohort, compute the Spearman correlation of the predictions and of (not nmd_escape) with the measured NMD
   efficiency. The 95% CIs are percentiles of N_RESAMPLES bootstrap resamples of the rows. The gain, model minus
   rules, has a paired CI from the same resamples. Each cohort draws from its own generator, seeded with SEED.
6. Per cohort, compute the mean measured NMD efficiency of the rows with nmd_escape False. For GTEx, also compute it
   over the variant-tissue rows, in total and per tissue with at least MIN_TISSUE_ROWS of them.

Outputs in --out-dir: downloads/, variants.vcf, features.parquet (the annotate result without sequences), the scored
rows tcga_rows.parquet, geuvadis_rows.parquet, gtex_rows.parquet (the variant-tissue rows),
gtex_variants.parquet and mmrf_target_rows.parquet, and validation.json and validation.md.
"""

import argparse
import dataclasses
import hashlib
import importlib.metadata
import io
import json
import logging
import time
import urllib.request
import zipfile
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import onnxruntime
import pandas as pd
import pyarrow.parquet as pq
from make_benchmark_vcfs import add_alleles, read_strands
from pyfaidx import Fasta
from scipy.stats import rankdata

import nmd_scanner
from nmd_scanner.schema import MODEL_INPUTS

logger = logging.getLogger("validate_model")

DEFAULT_ONNX = Path(__file__).resolve().parent.parent / "nmd_efficiency_rf.onnx"

# The commit that make_benchmark_vcfs.py pins, the tag v1.0 of NMDEff
MMRF_URL = (
    "https://raw.githubusercontent.com/hjkng/nmdeff/08c92768fcb689236a833db6a2f2d9bcbe919f12/MMRF_TARGET_dataset.csv"
)
MMRF_SHA256 = "c063914da054292c6c6b867eb58650d5378c74d5c852f54c874eab8e51cb4a19"
# Supplementary Data 1 to 12 of Kim et al. 2024 (Commun Biol 7:1461), one sheet each
GTEX_URL = (
    "https://static-content.springer.com/esm/art%3A10.1038%2Fs42003-024-07136-y/MediaObjects/"
    "42003_2024_7136_MOESM3_ESM.xlsx"
)
GTEX_SHA256 = "637c92b17095ef6147b74ae0ab6549f674adc6a67ee04c06d483fa7331ff34c6"
GTEX_SHEET = "Supplementary Data 2"
# Zenodo record 16666299 of Iha et al. 2025 has one file, NMD-rules.zip
ZENODO_URL = "https://zenodo.org/api/records/16666299/files/NMD-rules.zip/content"
ZENODO_SIZE = 3_189_701_188
# zipfile checks the CRC32 of each member, and zenodo_members checks the sha256
GEUVADIS_MEMBERS = {
    "NMD-rules/Datasets/Nonsense_list.txt": "5a248be88a6b2a287c2baf12e542d638257a06daadc21a4c440af812518bce26",
    "NMD-rules/Datasets/IsoformLevel.NMDannotation.txt": (
        "bc16f2b869149f07e0a5d7ce4173df4fbca9ea9c4fab2dd9f3c07a82c4011a7d"
    ),
}
# A Geuvadis variant with one of these flags in IsoformLevel.NMDannotation.txt is not in the accurate annotation
ACCURATE_ANNOTATION_FLAGS = ["resMNV", "nonTranslating", "Multiisoforms"]
# VAF_DNA of a heterozygous germline variant
GERMLINE_VAF_DNA = 0.5

KEY = ["transcript", "chrom", "pos", "ref", "alt"]
SEED = 0
N_RESAMPLES = 2000
MIN_TISSUE_ROWS = 100
GTEX_ROWS_UNIT = "variant-tissue rows"

# The rows of the README table: cohort, label, variants
COHORTS = [
    ("tcga", "TCGA (Kim et al. 2024), out of fold, folds grouped by chromosome", "somatic, tumors"),
    (
        "geuvadis",
        "Geuvadis (Iha et al. 2025), lymphoblastoid cell lines, the study's accurate annotation",
        "germline",
    ),
    ("gtex", "GTEx v8 (Kim et al. 2024), 54 tissues, mean per variant", "germline"),
    ("mmrf_target", "MMRF-TARGET (Kim et al. 2024)", "somatic, tumors"),
]


@dataclass(frozen=True)
class StepCount:
    step: str
    rows: int


@dataclass(frozen=True)
class CohortRows:
    """The scored rows of a cohort, with the measured NMD efficiency y, the prediction and nmd_escape."""

    rows: pd.DataFrame
    # What one row is, for example "variants"
    unit: str
    row_counts: list[StepCount]


@dataclass(frozen=True)
class Correlations:
    n: int
    model: float
    model_ci: list[float]
    rules: float
    rules_ci: list[float]
    gain: float
    gain_ci: list[float]


@dataclass(frozen=True)
class Level:
    unit: str
    n: int
    n_nmd_triggering: int
    # The mean measured NMD efficiency of the rows with nmd_escape False
    mean_nmd_triggering: float


@dataclass(frozen=True)
class CohortResult:
    cohort: str
    label: str
    variants: str
    correlations: Correlations
    level: Level
    row_counts: list[StepCount]


class RemoteFile:
    """
    A read-only file at a URL that answers HTTP range requests. Each read is one request, so zipfile can read single
    members of a remote zip.
    """

    def __init__(self, url: str) -> None:
        self.url = url
        self.position = 0
        _, content_range = self._fetch(0, 0)
        self.size = int(content_range.rsplit("/", 1)[1])

    def _fetch(self, start: int, end: int) -> tuple[bytes, str]:
        request = urllib.request.Request(self.url, headers={"Range": f"bytes={start}-{end}"})
        with urllib.request.urlopen(request, timeout=120) as response:
            # A server that ignores the range answers 200 and sends the whole file
            if response.status != 206:
                raise ValueError(f"{self.url} answered a range request with HTTP {response.status}")
            return response.read(), response.headers["Content-Range"]

    def seekable(self) -> bool:
        return True

    def tell(self) -> int:
        return self.position

    def seek(self, offset: int, whence: int = io.SEEK_SET) -> int:
        base = {io.SEEK_SET: 0, io.SEEK_CUR: self.position, io.SEEK_END: self.size}[whence]
        self.position = base + offset
        return self.position

    def read(self, n: int = -1) -> bytes:
        end = self.size if n < 0 else min(self.position + n, self.size)
        if end <= self.position:
            return b""
        data, _ = self._fetch(self.position, end - 1)
        self.position += len(data)
        return data


def parse_args() -> argparse.Namespace:
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
    parser.add_argument("--onnx", type=Path, default=DEFAULT_ONNX, help="model to validate (default: %(default)s)")
    parser.add_argument(
        "--geuvadis-dir",
        type=Path,
        help="local directory with Nonsense_list.txt and IsoformLevel.NMDannotation.txt of NMD-rules.zip; skips the"
        " range requests to Zenodo. The sha256 of both files is still checked.",
    )
    parser.add_argument(
        "--features",
        type=Path,
        help="features.parquet of an earlier run with the same inputs; skips nmd_scanner.annotate",
    )
    return parser.parse_args()


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def check_sha256(path: Path, expected: str) -> None:
    digest = sha256(path)
    if digest != expected:
        raise ValueError(f"{path} has sha256 {digest}, expected {expected}")


def download(url: str, expected_sha256: str, path: Path) -> Path:
    if not path.exists():
        logger.info("Downloading %s", url)
        with urllib.request.urlopen(url, timeout=120) as response:
            path.write_bytes(response.read())
    check_sha256(path, expected_sha256)
    return path


def zenodo_members(download_dir: Path, geuvadis_dir: Path | None) -> list[Path]:
    """
    Return the paths of the GEUVADIS_MEMBERS, in that order, and check their sha256. Read them from geuvadis_dir if it
    is given. Otherwise fetch the missing ones from the Zenodo zip to download_dir.
    """
    directory = download_dir if geuvadis_dir is None else geuvadis_dir
    paths = {member: directory / member.rsplit("/", 1)[1] for member in GEUVADIS_MEMBERS}
    missing = [member for member, path in paths.items() if not path.exists()]
    if missing and geuvadis_dir is None:
        logger.info("Reading %s from %s", ", ".join(missing), ZENODO_URL)
        remote = RemoteFile(ZENODO_URL)
        if remote.size != ZENODO_SIZE:
            raise ValueError(f"{ZENODO_URL} has {remote.size} bytes, expected {ZENODO_SIZE}")
        with zipfile.ZipFile(remote) as archive:
            for member in missing:
                paths[member].write_bytes(archive.read(member))
    for member, path in paths.items():
        check_sha256(path, GEUVADIS_MEMBERS[member])
    return list(paths.values())


def check_train_dir(train_dir: Path, onnx_path: Path) -> None:
    """Fail if the random forest of the training run differs from the model to validate."""
    trained = train_dir / "nmd_efficiency_rf.onnx"
    if sha256(trained) != sha256(onnx_path):
        raise ValueError(f"{trained} differs from {onnx_path}, so its out-of-fold predictions are of another model")


def read_mmrf(path: Path, strands: dict[str, str], fasta: Fasta) -> pd.DataFrame:
    return add_alleles(pd.read_csv(path), strands, fasta)


def read_gtex(path: Path, strands: dict[str, str], fasta: Fasta) -> pd.DataFrame:
    # The first line of the sheet is its title
    table = pd.read_excel(path, sheet_name=GTEX_SHEET, header=1)
    return add_alleles(table, strands, fasta).rename(columns={"Tissue_type": "tissue"})


def read_geuvadis(nonsense_list: Path, isoforms: Path, fasta: Fasta) -> pd.DataFrame:
    """
    Read the nonsense SNVs of Nonsense_list.txt, one row per variant and MANE Select transcript. Add the flags of
    IsoformLevel.NMDannotation.txt, which has the evaluated variants. REF and ALT are genomic. Fail if a REF does not
    match the FASTA.
    """
    table = pd.read_csv(nonsense_list, sep="\t")
    fasta_ref = [
        str(fasta[chrom][pos - 1 : pos - 1 + len(ref)]).upper()
        for chrom, pos, ref in zip(table["CHROM"], table["POS"], table["REF"], strict=True)
    ]
    mismatches = table["REF"].to_numpy() != np.array(fasta_ref)
    if mismatches.any():
        raise ValueError(f"{mismatches.sum()} rows of {nonsense_list.name} have a REF that does not match the FASTA")
    flags = pd.read_csv(isoforms, sep="\t", usecols=["Variant_ID", *ACCURATE_ANNOTATION_FLAGS])
    table = table.merge(flags, on="Variant_ID", how="left", validate="many_to_one")
    return table.rename(
        columns={"CHROM": "chrom", "POS": "pos", "REF": "ref", "ALT": "alt", "TranscriptID": "transcript"}
    )


def write_vcf(tables: list[pd.DataFrame], path: Path) -> int:
    variants = pd.concat([table[["chrom", "pos", "ref", "alt"]] for table in tables]).drop_duplicates()
    with open(path, "w") as handle:
        handle.write("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for chrom, pos, ref, alt in variants.itertuples(index=False):
            handle.write(f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\t.\t.\n")
    return len(variants)


def annotate(vcf: Path, gff3: Path, fasta: Path, out_path: Path) -> pd.DataFrame:
    started = time.time()
    features = nmd_scanner.annotate(vcf, gff3, fasta, sequences=False)
    logger.info("nmd_scanner.annotate: %d rows in %.0f s", len(features), time.time() - started)
    pq.write_table(nmd_scanner.to_arrow(features), out_path)
    return features


def score(features: pd.DataFrame, onnx_path: Path) -> pd.DataFrame:
    """
    Return KEY, nmd_model_status and nmd_escape of each NMD-Scanner row, and the prediction of the model for the rows
    with nmd_model_status "ok". Fail if the model does not take MODEL_INPUTS, or if 2 rows have the same KEY.
    """
    session = onnxruntime.InferenceSession(onnx_path, providers=["CPUExecutionProvider"])
    feature_names = json.loads(session.get_modelmeta().custom_metadata_map["feature_names"])
    if feature_names != list(MODEL_INPUTS):
        raise ValueError(f"{onnx_path} takes {feature_names}, not MODEL_INPUTS")
    scored = pd.DataFrame(
        {
            "transcript": features["transcript_id"].astype(str).str.split(".").str[0].to_numpy(),
            "chrom": features["chrom"].astype(str).to_numpy(),
            "pos": features["start"].astype("int64").to_numpy() + 1,
            "ref": features["ref"].astype(str).to_numpy(),
            "alt": features["alt"].astype(str).to_numpy(),
            "nmd_model_status": features["nmd_model_status"].astype(str).to_numpy(),
            "nmd_escape": features["nmd_escape"].to_numpy(),
        }
    )
    if scored.duplicated(KEY).any():
        raise ValueError("NMD-Scanner returned more than one row for a transcript and variant")
    ok = (scored["nmd_model_status"] == "ok").to_numpy()
    x = features[MODEL_INPUTS][ok].astype("float64").to_numpy()
    prediction = np.full(len(scored), np.nan)
    prediction[ok] = session.run(None, {"input": x})[0].ravel()
    return scored.assign(prediction=prediction)


def status_ok(joined: pd.DataFrame) -> pd.DataFrame:
    return joined[joined["nmd_model_status"] == "ok"].astype({"nmd_escape": bool})


def tcga_rows(train_dir: Path) -> CohortRows:
    # Every row of tcga_dataset.csv has VAF_RNA above 0
    rows = pd.read_parquet(train_dir / "oof_predictions.parquet").rename(columns={"NMD_efficiency": "y"})
    row_counts = [StepCount("out-of-fold rows, folds grouped by chromosome", len(rows))]
    rows = rows.astype({"nmd_escape": bool})
    return CohortRows(rows[["tcga_row", *KEY, "y", "prediction", "nmd_escape"]], "rows", row_counts)


def geuvadis_rows(measured: pd.DataFrame, scored: pd.DataFrame) -> CohortRows:
    # One row per variant and transcript
    joined = measured.merge(scored, on=KEY, how="inner", validate="one_to_one")
    ok = status_ok(joined)
    evaluated = ok[ok["VAF_mean"].notna()]
    positive = evaluated[evaluated["VAF_mean"] > 0]
    flagged = positive[ACCURATE_ANNOTATION_FLAGS].astype("boolean").fillna(False).any(axis=1)
    accurate = positive[~flagged]
    row_counts = [
        StepCount("Nonsense_list.txt rows", len(measured)),
        StepCount("rows joined to an NMD-Scanner row of their transcript", len(joined)),
        StepCount("nmd_model_status ok", len(ok)),
        StepCount("evaluated (VAF_mean present)", len(evaluated)),
        StepCount("VAF_mean > 0", len(positive)),
        StepCount("accurate annotation (no flag in IsoformLevel.NMDannotation.txt)", len(accurate)),
    ]
    rows = accurate.assign(y=-np.log2(accurate["VAF_mean"] / GERMLINE_VAF_DNA))
    return CohortRows(rows[[*KEY, "Variant_ID", "VAF_mean", "y", "prediction", "nmd_escape"]], "variants", row_counts)


def gtex_rows(measured: pd.DataFrame, scored: pd.DataFrame) -> tuple[CohortRows, pd.DataFrame]:
    """
    Return the GTEx variants, with the mean NMD efficiency over their tissues, and the variant-tissue rows that the
    means come from.
    """
    joined = measured.merge(scored, on=KEY, how="inner", validate="many_to_one")
    ok = status_ok(joined)
    # ASE_altered is VAF_RNA. For VAF_RNA 0, the table caps NMD_efficiency at 10.
    positive = ok[ok["ASE_altered"] > 0]
    rows = positive[[*KEY, "tissue", "ASE_altered", "NMD_efficiency", "prediction", "nmd_escape"]].rename(
        columns={"NMD_efficiency": "y"}
    )
    variants = (
        rows.groupby(KEY)
        .agg(
            y=("y", "mean"),
            prediction=("prediction", "first"),
            nmd_escape=("nmd_escape", "first"),
            tissues=("y", "size"),
        )
        .reset_index()
    )
    row_counts = [
        StepCount(f"{GTEX_SHEET} variant-tissue rows", len(measured)),
        StepCount("rows joined to an NMD-Scanner row of their transcript", len(joined)),
        StepCount("nmd_model_status ok", len(ok)),
        StepCount("ASE_altered (VAF_RNA) > 0", len(rows)),
        StepCount("variants, mean over their tissues", len(variants)),
    ]
    return CohortRows(variants, "variants, mean over their tissues", row_counts), rows


def mmrf_rows(measured: pd.DataFrame, scored: pd.DataFrame) -> CohortRows:
    joined = measured.merge(scored, on=KEY, how="inner", validate="many_to_one")
    ok = status_ok(joined)
    positive = ok[ok["VAF_RNA"] > 0]
    row_counts = [
        StepCount("MMRF_TARGET_dataset.csv rows", len(measured)),
        StepCount("rows joined to an NMD-Scanner row of their transcript", len(joined)),
        StepCount("nmd_model_status ok", len(ok)),
        StepCount("VAF_RNA > 0", len(positive)),
    ]
    rows = positive.rename(columns={"NMD_efficiency": "y"})
    return CohortRows(rows[[*KEY, "Project", "VAF_RNA", "y", "prediction", "nmd_escape"]], "rows", row_counts)


def spearman(x: np.ndarray, y: np.ndarray) -> float:
    return float(np.corrcoef(rankdata(x), rankdata(y))[0, 1])


def percentiles(values: np.ndarray) -> list[float]:
    return [float(v) for v in np.percentile(values, [2.5, 97.5])]


def correlations(rows: pd.DataFrame) -> Correlations:
    """
    Spearman correlation of the prediction and of the rules with y, and their difference, with bootstrap CIs. The
    rules score is (not nmd_escape), so that both scores rank high NMD efficiency first.
    """
    y = rows["y"].to_numpy(dtype="float64")
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
    return Correlations(
        n=len(y),
        model=model_value,
        model_ci=percentiles(model_resampled),
        rules=rules_value,
        rules_ci=percentiles(rules_resampled),
        gain=model_value - rules_value,
        gain_ci=percentiles(model_resampled - rules_resampled),
    )


def level(rows: pd.DataFrame, unit: str) -> Level:
    triggering = rows.loc[~rows["nmd_escape"], "y"]
    return Level(unit=unit, n=len(rows), n_nmd_triggering=len(triggering), mean_nmd_triggering=float(triggering.mean()))


def tissue_levels(rows: pd.DataFrame) -> dict[str, Level]:
    """
    The level of each GTEx tissue with at least MIN_TISSUE_ROWS variant-tissue rows, sorted by mean_nmd_triggering.
    """
    counts = rows["tissue"].value_counts()
    levels = {
        tissue: level(rows[rows["tissue"] == tissue], GTEX_ROWS_UNIT)
        for tissue in counts.index[counts >= MIN_TISSUE_ROWS]
    }
    return dict(sorted(levels.items(), key=lambda item: item[1].mean_nmd_triggering))


def estimate(value: float, ci: list[float]) -> str:
    return f"{value:.2f} ({ci[0]:.2f} to {ci[1]:.2f})"


def markdown(results: list[CohortResult], gtex_level: Level, tissues: dict[str, Level]) -> str:
    lines = [
        "# Validation of nmd_efficiency_rf.onnx",
        "",
        'Spearman correlation with the measured NMD efficiency, with 95% bootstrap CIs. The column "Canonical'
        ' rules" ranks the same rows by `nmd_escape` alone.',
        "",
        "| Cohort | Variants | Rows | Model | Canonical rules |",
        "| --- | --- | ---: | --- | --- |",
    ]
    lines += [
        f"| {r.label} | {r.variants} | {r.correlations.n:,} | {estimate(r.correlations.model, r.correlations.model_ci)}"
        f" | {estimate(r.correlations.rules, r.correlations.rules_ci)} |"
        for r in results
    ]
    lines += [
        "",
        f"{N_RESAMPLES:,} bootstrap resamples of the rows per cohort, seed {SEED}. The gain is model minus rules, with"
        " a paired CI.",
        "",
        "| Cohort | Gain |",
        "| --- | --- |",
    ]
    lines += [f"| {r.label} | {estimate(r.correlations.gain, r.correlations.gain_ci)} |" for r in results]
    lines += [
        "",
        "Mean measured NMD efficiency of the rows with `nmd_escape` False (NMD-triggering).",
        "",
        "| Cohort | Unit | All | `nmd_escape` False | Mean |",
        "| --- | --- | ---: | ---: | ---: |",
    ]
    levels = [(r.label, r.level) for r in results] + [("GTEx v8 (Kim et al. 2024)", gtex_level)]
    lines += [
        f"| {label} | {cohort_level.unit} | {cohort_level.n:,} | {cohort_level.n_nmd_triggering:,} |"
        f" {cohort_level.mean_nmd_triggering:.2f} |"
        for label, cohort_level in levels
    ]
    lowest, highest = next(iter(tissues.items())), next(reversed(tissues.items()))
    lines += [
        "",
        f"Between the {len(tissues)} GTEx tissues with at least {MIN_TISSUE_ROWS} {GTEX_ROWS_UNIT}, the mean ranges"
        f" from {lowest[1].mean_nmd_triggering:.2f} ({lowest[0]}) to {highest[1].mean_nmd_triggering:.2f}"
        f" ({highest[0]}).",
        "",
        "| Cohort | Step | Rows |",
        "| --- | --- | ---: |",
    ]
    lines += [f"| {r.cohort} | {c.step} | {c.rows:,} |" for r in results for c in r.row_counts]
    return "\n".join(lines) + "\n"


def main() -> None:
    args = parse_args()
    logging.basicConfig(level=logging.WARNING, format="%(asctime)s %(name)s %(levelname)s %(message)s")
    logger.setLevel(logging.INFO)
    out_dir = args.out_dir
    download_dir = out_dir / "downloads"
    download_dir.mkdir(parents=True, exist_ok=True)
    check_train_dir(args.train_dir, args.onnx)

    mmrf_path = download(MMRF_URL, MMRF_SHA256, download_dir / "MMRF_TARGET_dataset.csv")
    gtex_path = download(GTEX_URL, GTEX_SHA256, download_dir / "42003_2024_7136_MOESM3_ESM.xlsx")
    nonsense_list, isoforms = zenodo_members(download_dir, args.geuvadis_dir)

    strands = read_strands(args.gff3)
    fasta = Fasta(str(args.fasta))
    geuvadis = read_geuvadis(nonsense_list, isoforms, fasta)
    gtex = read_gtex(gtex_path, strands, fasta)
    mmrf = read_mmrf(mmrf_path, strands, fasta)
    vcf_path = out_dir / "variants.vcf"
    n_variants = write_vcf([geuvadis, gtex, mmrf], vcf_path)
    logger.info("unique variants in %s: %d", vcf_path.name, n_variants)

    if args.features is None:
        features = annotate(vcf_path, args.gff3, args.fasta, out_dir / "features.parquet")
    else:
        features = pd.read_parquet(args.features)
    scored = score(features, args.onnx)

    gtex_variants, gtex_variant_tissue_rows = gtex_rows(gtex, scored)
    cohort_rows = {
        "tcga": tcga_rows(args.train_dir),
        "geuvadis": geuvadis_rows(geuvadis, scored),
        "gtex": gtex_variants,
        "mmrf_target": mmrf_rows(mmrf, scored),
    }
    outputs = {
        "tcga_rows": cohort_rows["tcga"].rows,
        "geuvadis_rows": cohort_rows["geuvadis"].rows,
        "gtex_rows": gtex_variant_tissue_rows,
        "gtex_variants": cohort_rows["gtex"].rows,
        "mmrf_target_rows": cohort_rows["mmrf_target"].rows,
    }
    for name, rows in outputs.items():
        rows.to_parquet(out_dir / f"{name}.parquet", index=False)

    results = []
    for cohort, label, variants in COHORTS:
        data = cohort_rows[cohort]
        for count in data.row_counts:
            logger.info("%s: %s: %d", cohort, count.step, count.rows)
        result = CohortResult(
            cohort=cohort,
            label=label,
            variants=variants,
            correlations=correlations(data.rows),
            level=level(data.rows, data.unit),
            row_counts=data.row_counts,
        )
        logger.info(
            "%s: model %.3f, rules %.3f, gain %.3f",
            cohort,
            result.correlations.model,
            result.correlations.rules,
            result.correlations.gain,
        )
        results.append(result)
    gtex_level = level(gtex_variant_tissue_rows, GTEX_ROWS_UNIT)
    tissues = tissue_levels(gtex_variant_tissue_rows)

    validation = {
        "versions": {
            name: importlib.metadata.version(name)
            for name in ("nmd-scanner", "onnxruntime", "pandas", "numpy", "scipy")
        },
        "inputs": {
            "onnx": args.onnx.name,
            "onnx_sha256": sha256(args.onnx),
            "mmrf_url": MMRF_URL,
            "mmrf_sha256": MMRF_SHA256,
            "gtex_url": GTEX_URL,
            "gtex_sha256": GTEX_SHA256,
            "gtex_sheet": GTEX_SHEET,
            "zenodo_url": ZENODO_URL,
            "zenodo_members_sha256": GEUVADIS_MEMBERS,
            "gff3": args.gff3.name,
            "fasta": args.fasta.name,
        },
        "unique_variants_annotated": n_variants,
        "bootstrap": {"resamples": N_RESAMPLES, "seed": SEED, "ci": "percentiles 2.5 and 97.5"},
        "cohorts": [dataclasses.asdict(r) for r in results],
        "gtex_variant_tissue_rows": {
            "level": dataclasses.asdict(gtex_level),
            "tissues": int(gtex_variant_tissue_rows["tissue"].nunique()),
            "min_rows_per_tissue": MIN_TISSUE_ROWS,
            "levels_per_tissue": {tissue: dataclasses.asdict(tissue_level) for tissue, tissue_level in tissues.items()},
        },
    }
    (out_dir / "validation.json").write_text(json.dumps(validation, indent=2) + "\n")
    (out_dir / "validation.md").write_text(markdown(results, gtex_level, tissues))
    logger.info("Wrote %s", out_dir / "validation.md")


if __name__ == "__main__":
    main()
