"""Tests for PAF and GFF3 output functionality."""

import os
import tempfile

from tsplit.utils import resolve_output_path, write_gff3, write_paf


def test_resolve_output_path_none():
    """Test that None filepath returns None."""
    assert resolve_output_path(None, None) is None
    assert resolve_output_path(None, '/some/dir') is None


def test_resolve_output_path_absolute():
    """Test that absolute paths are returned unchanged."""
    abs_path = '/absolute/path/to/file.txt'
    assert resolve_output_path(abs_path, None) == abs_path
    assert resolve_output_path(abs_path, '/some/dir') == abs_path


def test_resolve_output_path_relative_no_outdir():
    """Test that relative paths without outdir are returned unchanged."""
    rel_path = 'relative/path/file.txt'
    assert resolve_output_path(rel_path, None) == rel_path


def test_resolve_output_path_relative_with_outdir():
    """Test that relative paths are joined with outdir."""
    with tempfile.TemporaryDirectory() as tmpdir:
        rel_path = 'output.txt'
        result = resolve_output_path(rel_path, tmpdir)
        expected = os.path.join(os.path.abspath(tmpdir), rel_path)
        assert result == expected

        # Test with a more complex relative path
        rel_path2 = 'subdir/output.txt'
        result2 = resolve_output_path(rel_path2, tmpdir)
        expected2 = os.path.join(os.path.abspath(tmpdir), rel_path2)
        assert result2 == expected2


def test_write_paf():
    """Test that PAF output is correctly formatted."""
    # Sample alignment data
    alignment_data = [
        {
            'qry_name': 'seq1',
            'qry_length': 1000,
            'qry_start': 100,
            'qry_end': 200,
            'strand': '+',
            'ref_name': 'seq1',
            'ref_length': 1000,
            'ref_start': 0,
            'ref_end': 100,
            'num_matches': 95,
            'aln_block_length': 100,
            'mapping_quality': 255,
        }
    ]

    with tempfile.TemporaryDirectory() as tmpdir:
        paf_file = os.path.join(tmpdir, 'test.paf')
        write_paf(alignment_data, paf_file)

        # Check file exists
        assert os.path.exists(paf_file)

        # Check content
        with open(paf_file) as f:
            lines = f.readlines()
            assert len(lines) == 1
            fields = lines[0].strip().split('\t')
            assert len(fields) == 12
            assert fields[0] == 'seq1'
            assert fields[1] == '1000'
            assert fields[2] == '100'
            assert fields[3] == '200'
            assert fields[4] == '+'
            assert fields[5] == 'seq1'
            assert fields[6] == '1000'
            assert fields[7] == '0'
            assert fields[8] == '100'
            assert fields[9] == '95'
            assert fields[10] == '100'
            assert fields[11] == '255'


def test_write_gff3():
    """Test that GFF3 output is correctly formatted."""
    # Sample feature data
    feature_data = [
        {
            'seqid': 'seq1',
            'type': 'TIR',
            'start': 1,
            'end': 100,
            'score': '95.5',
            'strand': '+',
            'phase': '.',
            'attributes': 'ID=seq1_TIR_L;Name=Left_TIR;length=100',
        },
        {
            'seqid': 'seq1',
            'type': 'TIR',
            'start': 901,
            'end': 1000,
            'score': '95.5',
            'strand': '-',
            'phase': '.',
            'attributes': 'ID=seq1_TIR_R;Name=Right_TIR;length=100',
        },
    ]

    with tempfile.TemporaryDirectory() as tmpdir:
        gff_file = os.path.join(tmpdir, 'test.gff3')
        write_gff3(feature_data, gff_file, source='tSplit_test')

        # Check file exists
        assert os.path.exists(gff_file)

        # Check content
        with open(gff_file) as f:
            lines = f.readlines()
            # Header + 2 features
            assert len(lines) == 3
            assert lines[0].strip() == '##gff-version 3'

            # Check first feature
            fields1 = lines[1].strip().split('\t')
            assert len(fields1) == 9
            assert fields1[0] == 'seq1'
            assert fields1[1] == 'tSplit_test'
            assert fields1[2] == 'TIR'
            assert fields1[3] == '1'
            assert fields1[4] == '100'
            assert fields1[5] == '95.5'
            assert fields1[6] == '+'
            assert fields1[7] == '.'
            assert 'ID=seq1_TIR_L' in fields1[8]

            # Check second feature
            fields2 = lines[2].strip().split('\t')
            assert len(fields2) == 9
            assert fields2[0] == 'seq1'
            assert fields2[2] == 'TIR'
            assert fields2[3] == '901'
            assert fields2[4] == '1000'
            assert fields2[6] == '-'


def test_write_empty_paf():
    """Test that empty PAF files are created correctly."""
    with tempfile.TemporaryDirectory() as tmpdir:
        paf_file = os.path.join(tmpdir, 'empty.paf')
        write_paf([], paf_file)

        assert os.path.exists(paf_file)
        with open(paf_file) as f:
            content = f.read()
            assert content == ''


def test_write_empty_gff3():
    """Test that empty GFF3 files still have header."""
    with tempfile.TemporaryDirectory() as tmpdir:
        gff_file = os.path.join(tmpdir, 'empty.gff3')
        write_gff3([], gff_file)

        assert os.path.exists(gff_file)
        with open(gff_file) as f:
            lines = f.readlines()
            assert len(lines) == 1
            assert lines[0].strip() == '##gff-version 3'
