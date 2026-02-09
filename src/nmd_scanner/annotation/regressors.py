from importlib.resources import files
import joblib
import numpy as np

path = files('nmd_scanner.assets').joinpath('nmdeff_regressor.pkl')

def nmdeff_predict(start_loss: np.ndarray, stop_loss: np.ndarray, total_exon_count: np.ndarray, ptc_less_than_150nt_to_start: np.ndarray,
                   nmd_long_exon_rule: np.ndarray, nmd_start_proximal_rule: np.ndarray, nmd_single_exon_rule: np.ndarray, nmd_escape: np.ndarray,
                   downstream_exon_count: np.ndarray, nmd_last_exon_rule: np.ndarray, ptc_to_start_codon: np.ndarray, stop_codon_distance: np.ndarray,
                   ptc_exon_length: np.ndarray, ptc_to_intron: np.ndarray, upstream_exon_count: np.ndarray, nmd_50nt_penultimate_rule: np.ndarray,
                   utr5_length: np.ndarray, utr3_length: np.ndarray, transcript_length: np.ndarray) -> np.ndarray:
    
    with open(path, 'rb') as f:
        model = joblib.load(f)
        X = np.column_stack((start_loss, stop_loss, total_exon_count, ptc_less_than_150nt_to_start,
                             nmd_long_exon_rule, nmd_start_proximal_rule, nmd_single_exon_rule, nmd_escape,
                             downstream_exon_count, nmd_last_exon_rule, ptc_to_start_codon, stop_codon_distance,
                             ptc_exon_length, ptc_to_intron, upstream_exon_count, nmd_50nt_penultimate_rule,
                             utr5_length, utr3_length, transcript_length))
        return model.predict(X)
    