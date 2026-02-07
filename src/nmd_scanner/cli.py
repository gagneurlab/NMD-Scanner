# Import dependencies
import click
import os
from . import scanner


def is_valid_output_path(path):
    """
    Validate if the output path exists or is creatable.
    Allows a file path (where the parent directory must exist) or a directory.
    """

    if os.path.exists(path):
        return True
    parent = os.path.dirname(path)
    return os.path.isdir(parent) if parent else False


@click.command()
@click.option('--vcf', required=True, type=click.Path(exists=True),
              help='Path to VCF file')
@click.option('--gtf', required=True, type=click.Path(exists=True),
              help='Path to GTF file')
@click.option('--fasta', required=True, type=click.Path(exists=True),
              help='Path to FASTA file')
@click.option('--output', required=True,
              help='Path to output file or output directory')
@click.option('--reassign-exons', is_flag=True,
              help='Recompute exon numbers (recommended for hg19; may be slow)')
def cli(vcf, gtf, fasta, output, reassign_exons):
    """Run NMD pipeline to analyze variants for nonsense-mediated decay."""
    
    # Check that the output path is valid
    if not is_valid_output_path(output):
        raise click.BadParameter(f"Invalid output path: {output}")

    # Run the main pipeline
    scanner.annotate_nmd(vcf, gtf, fasta, output, reassign_exons=reassign_exons)


if __name__ == '__main__':
    cli()
