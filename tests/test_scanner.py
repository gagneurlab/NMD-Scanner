import nmd_scanner
import pandas as pd
import os

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


def test_output_format_csv(vcf_path, gtf_path, fasta_path, tmp_path):
    """Test CSV output format"""
    output_file = tmp_path / "test_output.csv"
    df = nmd_scanner.annotate_nmd(vcf_path, gtf_path, fasta_path, str(output_file), output_format='csv')
    
    # Check that file was created
    assert output_file.exists()
    
    # Read back and verify
    df_read = pd.read_csv(output_file)
    assert df_read.shape == df.shape
    assert list(df_read.columns) == list(df.columns)


def test_output_format_parquet(vcf_path, gtf_path, fasta_path, tmp_path):
    """Test parquet output format"""
    output_file = tmp_path / "test_output.parquet"
    df = nmd_scanner.annotate_nmd(vcf_path, gtf_path, fasta_path, str(output_file), output_format='parquet')
    
    # Check that file was created
    assert output_file.exists()
    
    # Read back and verify
    df_read = pd.read_parquet(output_file)
    assert df_read.shape == df.shape
    assert list(df_read.columns) == list(df.columns)


def test_parquet_format_simple(tmp_path):
    """Test parquet write/read with simple data"""
    # Create a simple dataframe
    test_df = pd.DataFrame({
        'variant_id': ['var1', 'var2'],
        'nmd_escape': [True, False],
        'value': [1.5, 2.5]
    })
    
    # Write as parquet
    output_file = tmp_path / "simple_test.parquet"
    test_df.to_parquet(output_file, index=False)
    
    # Read back
    df_read = pd.read_parquet(output_file)
    
    assert df_read.shape == test_df.shape
    assert list(df_read.columns) == list(test_df.columns)
    assert df_read['variant_id'].tolist() == ['var1', 'var2']


def test_output_format_inferred_from_extension(vcf_path, gtf_path, fasta_path, tmp_path):
    """Test that output format is inferred from file extension"""
    # Test parquet inference
    parquet_file = tmp_path / "test_inferred.parquet"
    df1 = nmd_scanner.annotate_nmd(vcf_path, gtf_path, fasta_path, str(parquet_file))
    assert parquet_file.exists()
    df_parquet = pd.read_parquet(parquet_file)
    assert df_parquet.shape == df1.shape
    
    # Test csv inference
    csv_file = tmp_path / "test_inferred.csv"
    df2 = nmd_scanner.annotate_nmd(vcf_path, gtf_path, fasta_path, str(csv_file))
    assert csv_file.exists()
    df_csv = pd.read_csv(csv_file)
    assert df_csv.shape == df2.shape


def test_output_format_directory_default_csv(vcf_path, gtf_path, fasta_path, tmp_path):
    """Test that directory output defaults to CSV"""
    df = nmd_scanner.annotate_nmd(vcf_path, gtf_path, fasta_path, str(tmp_path))
    
    # Should create CSV file in directory
    csv_files = list(tmp_path.glob("*.csv"))
    assert len(csv_files) == 1
    assert csv_files[0].exists()


def test_output_format_directory_with_format_option(vcf_path, gtf_path, fasta_path, tmp_path):
    """Test specifying format when output is a directory"""
    df = nmd_scanner.annotate_nmd(vcf_path, gtf_path, fasta_path, str(tmp_path), output_format='parquet')
    
    # Should create parquet file in directory
    parquet_files = list(tmp_path.glob("*.parquet"))
    assert len(parquet_files) == 1
    assert parquet_files[0].exists()
