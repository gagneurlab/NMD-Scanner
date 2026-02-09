from importlib.resources import files
import joblib
import numpy as np
import pandas as pd

path = files('nmd_scanner.assets').joinpath('nmdeff_regressor.pkl')

def nmdeff_predict(start_loss: np.ndarray, stop_loss: np.ndarray, total_exon_count: np.ndarray, ptc_less_than_150nt_to_start: np.ndarray,
                   nmd_long_exon_rule: np.ndarray, nmd_start_proximal_rule: np.ndarray, nmd_single_exon_rule: np.ndarray, nmd_escape: np.ndarray,
                   downstream_exon_count: np.ndarray, nmd_last_exon_rule: np.ndarray, ptc_to_start_codon: np.ndarray, stop_codon_distance: np.ndarray,
                   ptc_exon_length: np.ndarray, ptc_to_intron: np.ndarray, upstream_exon_count: np.ndarray, nmd_50nt_penultimate_rule: np.ndarray,
                   utr5_length: np.ndarray, utr3_length: np.ndarray, transcript_length: np.ndarray) -> np.ndarray:
    
    with open(path, 'rb') as f:
        model = joblib.load(f)
        
        # Create DataFrame with proper feature names to avoid sklearn warning
        X = pd.DataFrame({
            'start_loss': start_loss.flatten(),
            'stop_loss': stop_loss.flatten(),
            'total_exon_count': total_exon_count.flatten(),
            'ptc_less_than_150nt_to_start': ptc_less_than_150nt_to_start.flatten(),
            'nmd_long_exon_rule': nmd_long_exon_rule.flatten(),
            'nmd_start_proximal_rule': nmd_start_proximal_rule.flatten(),
            'nmd_single_exon_rule': nmd_single_exon_rule.flatten(),
            'nmd_escape': nmd_escape.flatten(),
            'downstream_exon_count': downstream_exon_count.flatten(),
            'nmd_last_exon_rule': nmd_last_exon_rule.flatten(),
            'ptc_to_start_codon': ptc_to_start_codon.flatten(),
            'stop_codon_distance': stop_codon_distance.flatten(),
            'ptc_exon_length': ptc_exon_length.flatten(),
            'ptc_to_intron': ptc_to_intron.flatten(),
            'upstream_exon_count': upstream_exon_count.flatten(),
            'nmd_50nt_penultimate_rule': nmd_50nt_penultimate_rule.flatten(),
            'utr5_length': utr5_length.flatten(),
            'utr3_length': utr3_length.flatten(),
            'transcript_length': transcript_length.flatten()
        })
        
        return model.predict(X)
    