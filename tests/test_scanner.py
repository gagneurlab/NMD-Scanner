import nmd_scanner

def test_annotate_nmd_minimal(vcf_path, gtf_path, fasta_path, tmp_path):
    """Test minimal output mode with just IDs and NMD predictions"""
    output = str(tmp_path)
    df = nmd_scanner.annotate_nmd(vcf_path, gtf_path, fasta_path, output, detailed=False)
    
    # Check DataFrame structure
    assert df is not None
    assert df.shape[0] > 0
    
    # Check minimal columns are present
    assert 'variant_id' in df.columns
    assert 'transcript_id' in df.columns
    assert 'gene_id' in df.columns
    assert 'chromosome' in df.columns
    assert 'alt_is_premature' in df.columns
    assert 'nmd_escape' in df.columns
    assert 'nmd_last_exon_rule' in df.columns
    assert 'nmd_50nt_penultimate_rule' in df.columns
    assert 'nmd_long_exon_rule' in df.columns
    assert 'nmd_start_proximal_rule' in df.columns
    assert 'nmd_single_exon_rule' in df.columns
    
    # Check that transcript features ARE present in minimal mode
    assert 'utr5_length' in df.columns
    assert 'utr3_length' in df.columns
    assert 'total_exon_count' in df.columns
    
    # Check that detailed CDS columns are NOT present
    assert 'ref_cds_seq' not in df.columns
    assert 'alt_cds_seq' not in df.columns
    
    print(f"Minimal output shape: {df.shape}")
    print(df.head())

def test_annotate_nmd_detailed(vcf_path, gtf_path, fasta_path, tmp_path):
    """Test detailed output mode with all CDS, transcript, and NMD fields"""
    output = str(tmp_path)
    df = nmd_scanner.annotate_nmd(vcf_path, gtf_path, fasta_path, output, detailed=True)
    
    # Check DataFrame structure
    assert df is not None
    assert df.shape[0] > 0
    
    # Check that all detailed columns are present
    assert 'variant_id' in df.columns
    assert 'transcript_id' in df.columns
    assert 'gene_id' in df.columns
    assert 'ref_cds_seq' in df.columns
    assert 'alt_cds_seq' in df.columns
    assert 'utr5_length' in df.columns
    assert 'utr3_length' in df.columns
    assert 'total_exon_count' in df.columns
    assert 'nmd_escape' in df.columns
    
    print(f"Detailed output shape: {df.shape}")
    print(df.head())

def test_nmd_escape_logic(vcf_path, gtf_path, fasta_path):
    """Test that NMD escape logic produces boolean values"""
    df = nmd_scanner.annotate_nmd(vcf_path, gtf_path, fasta_path, detailed=False)
    
    # Check that NMD escape columns have boolean values
    assert df['nmd_escape'].dtype == bool or df['nmd_escape'].dtype == 'object'
    assert df['nmd_last_exon_rule'].dtype == bool or df['nmd_last_exon_rule'].dtype == 'object'
    
    # Check that at least some variants have NMD predictions
    assert df['alt_is_premature'].notna().any()
    
    print(f"Premature variants: {df['alt_is_premature'].sum()}")
    print(f"NMD escape: {df['nmd_escape'].sum()}")