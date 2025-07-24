# for faster access to the fasta sequence

from functools import lru_cache

# Global cache to map id -> actual fasta object
_fasta_cache = {}

@lru_cache(maxsize=None)
def fetch_seq_cached(chrom, start, end, fasta_id):
    fasta = _fasta_cache[fasta_id]  # lookup the actual fasta object
    return fasta[chrom][start:end].seq.upper()

def add_exon_cds_sequence(df, fasta, chrom_col="Chromosome", start_col="Start", end_col="End", new_col="Exon_CDS_seq"):
    fasta_id = id(fasta)
    _fasta_cache[fasta_id] = fasta  # store for later access

    def extract_seq(row):
        return fetch_seq_cached(row[chrom_col], row[start_col], row[end_col], fasta_id)

    df[new_col] = df.apply(extract_seq, axis=1)
    return df