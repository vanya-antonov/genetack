
from io import StringIO
import os
import pandas as pd
from pprint import pprint
import subprocess
import tempfile

from Bio import SeqIO


TEST_DATA_DIR = 'data'


def test_genetack_gm_run(tmpdir):
    fsmark_dir = os.path.join(tmpdir, 'fsmark')
    prm_str = ' '.join([
        '--output_filtered_fs',
        '--save_fsmark_files', fsmark_dir,
        '--gm_mod_fn', os.path.join(TEST_DATA_DIR, 'S_griseus.gm_mod.txt'),
        '--fs_mod_fn', os.path.join(TEST_DATA_DIR, 'S_griseus.fs_mod.txt')
    ])

    fna_fn = os.path.join(TEST_DATA_DIR, 'S_griseus.50kb.fasta')
    cmd_str = 'genetack_gm.pl ' + prm_str + ' ' + fna_fn
    gtgm_out = subprocess.run(
        cmd_str, shell=True, stdout=subprocess.PIPE,
        universal_newlines=True, stderr=subprocess.DEVNULL
    ).stdout

    # Create Bio.Seq object
    seq = SeqIO.read(fna_fn, 'fasta').seq

    all_fs_df = pd.read_csv(
        StringIO(gtgm_out), sep='\s+')

    # Validate the number of predctions
    assert all_fs_df.shape[0] == 27

    # Test that the file exists
    assert os.path.exists(fsmark_dir)

    fsmark_fn = os.path.join(
        fsmark_dir, str(all_fs_df.iloc[0].Fragment_Left) + '.fsmark')
    assert os.path.isfile(fsmark_fn)

    _validate_zero_based_coordinates(all_fs_df, fsmark_dir, seq)

def _validate_zero_based_coordinates(all_fs_df, fsmark_dir, seq):
    for index, fs in all_fs_df.iterrows():
        fsmark_fn = os.path.join(
            fsmark_dir, str(fs.Fragment_Left) + '.fsmark')
        assert os.path.isfile(fsmark_fn)

        nt_seq_from_coords = seq[fs.Fragment_Left:fs.Fragment_Right]
        if fs.Strand == '-':
            nt_seq_from_coords = nt_seq_from_coords.reverse_complement()

        fsmark_df = pd.read_csv(fsmark_fn, sep='\s+')
        nt_seq_from_fsmark = ''.join(fsmark_df.Letter.to_list())
        assert nt_seq_from_coords.upper() == nt_seq_from_fsmark.upper()

        gene_seq = seq[fs.Gene_Left:fs.Gene_Right]
        if fs.Strand == '-':
            gene_seq = gene_seq.reverse_complement()

        stop_codon = gene_seq[-3:]
        assert stop_codon.upper() in ['TGA', 'TAG', 'TAA']

        start_codon = gene_seq[0:3]
        assert start_codon.upper() in ['ATG', 'GTG', 'TTG']


def test_fs_coord_after_codon():
    fsmark_dir = os.path.join(TEST_DATA_DIR, 'S_griseus.50kb.fsmark_dir')
    gtgm_fn = os.path.join(TEST_DATA_DIR, 'S_griseus.50kb.genetack_gm')
    fna_fn = os.path.join(TEST_DATA_DIR, 'S_griseus.50kb.fasta')

    # Create Bio.Seq object
    seq = SeqIO.read(fna_fn, 'fasta').seq

    # https://stackoverflow.com/a/22605281/310453
    # https://stackoverflow.com/a/31324373/310453
    all_fs_df = pd.read_csv(gtgm_fn, sep='\s+')

    _validate_single_fs_genes(all_fs_df, seq)

    for index, fs in all_fs_df.iterrows():
        fsmark_fn = os.path.join(fsmark_dir, str(fs.Fragment_Left) + '.fsmark')

        # Test that the file exists
        assert os.path.exists(fsmark_fn)

        _validate_fshift(fs, fsmark_fn)

def _validate_single_fs_genes(all_fs_df, seq):
    # Remove genes with >1 fs
    # https://stackoverflow.com/a/34272155/310453
    df = all_fs_df.drop_duplicates(
        subset=['Gene_Left', 'Gene_Right'], keep=False)
    for index, fs in df.iterrows():
        # Get the gene seq before (upstream) FS
        # Need to compute num_upstream_codons because gene border
        # doesn't correspond to end of the start codon
        if fs.Strand == '+':
            num_upstream_codons = int((fs.FS_coord_adj - fs.Gene_Left) / 3)
            coord_1 = fs.FS_coord_adj - 3*num_upstream_codons
            coord_2 = fs.FS_coord_adj
        else:
            num_upstream_codons = int((fs.Gene_Right - fs.FS_coord_adj) / 3)
            coord_1 = fs.FS_coord_adj
            coord_2 = fs.FS_coord_adj + 3*num_upstream_codons

        # Make sure we have complete codons
        assert (coord_2 - coord_1) % 3 == 0

        upstream_nt = seq[coord_1:coord_2]
        if fs.Strand == '-':
            upstream_nt = upstream_nt.reverse_complement()

        # Make sure there is no stop-codons
        upstream_aa = upstream_nt.translate()
        assert('*' not in upstream_aa)


def _validate_fshift(fs, fsmark_fn):
    fsmark_df = pd.read_csv(fsmark_fn, sep='\t')
    if fs.Strand == '+':
        coord_in_fragment = fs.FS_coord_adj - fs.Fragment_Left
    else:
        coord_in_fragment = fs.Fragment_Right - fs.FS_coord_adj
    # state:  0 1 2 0 1 2
    #         A T G T T C
    #        | | | | | | |
    #        0 1 2 3 4 5 6
    # Correct fs-coord is 3 because the letter before it has the last codon position
    assert fsmark_df.iloc[coord_in_fragment-1].Emission == '2_wo_stop'

