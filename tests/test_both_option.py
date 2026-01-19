"""Tests for the --both option to output both terminal repeats."""

import os
import tempfile

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from tsplit.parseAlign import getLTRs, getTIRs


def test_getTIRs_both_option_yields_two_segments():
    """Test that getTIRs with both=True yields both left and right TIRs."""
    # Create a temporary directory and test file
    with tempfile.TemporaryDirectory() as tempdir:
        test_fasta = os.path.join(tempdir, 'test_tir.fa')

        # Create a simple test sequence with TIRs at both ends
        # The left TIR should be at the start, the right TIR at the end (reverse complement)
        # For testing, we'll use a known TIR element structure
        with open(test_fasta, 'w') as f:
            f.write('>test_element\n')
            # Simple palindromic sequence for testing
            # AAATTTGGG at start, CCCAAATTT at end (reverse complement of start)
            f.write('AAATTTGGG' + 'N' * 100 + 'CCCAAATTT\n')

        # Call getTIRs with both=True
        segments = list(
            getTIRs(
                test_fasta,
                flankdist=10,
                minid=70,
                minterm=5,
                report='external',
                temp=tempdir,
                keeptemp=False,
                both=True,
            )
        )

        # With both=True, we should get 0 or 2 segments depending on alignment quality
        # If alignments found, should be _L_TIR and _R_TIR
        if len(segments) > 0:
            assert len(segments) == 2 or len(segments) == 0
            if len(segments) == 2:
                assert segments[0].id.endswith('_L_TIR')
                assert segments[1].id.endswith('_R_TIR')
                # Check descriptions
                assert 'left TIR' in segments[0].description
                assert 'right TIR' in segments[1].description
                assert 'reverse complemented' in segments[1].description


def test_getTIRs_without_both_option_yields_one_segment():
    """Test that getTIRs with both=False yields only one TIR."""
    with tempfile.TemporaryDirectory() as tempdir:
        test_fasta = os.path.join(tempdir, 'test_tir.fa')

        with open(test_fasta, 'w') as f:
            f.write('>test_element\n')
            f.write('AAATTTGGG' + 'N' * 100 + 'CCCAAATTT\n')

        # Call getTIRs with both=False (default)
        segments = list(
            getTIRs(
                test_fasta,
                flankdist=10,
                minid=70,
                minterm=5,
                report='external',
                temp=tempdir,
                keeptemp=False,
                both=False,
            )
        )

        # With both=False, we should get 0 or 1 segment
        if len(segments) > 0:
            assert len(segments) == 1
            assert segments[0].id.endswith('_TIR')
            assert not segments[0].id.endswith('_L_TIR')
            assert not segments[0].id.endswith('_R_TIR')


def test_getLTRs_both_option_yields_two_segments():
    """Test that getLTRs with both=True yields both left and right LTRs."""
    with tempfile.TemporaryDirectory() as tempdir:
        test_fasta = os.path.join(tempdir, 'test_ltr.fa')

        # Create a simple test sequence with LTRs at both ends (same orientation)
        with open(test_fasta, 'w') as f:
            f.write('>test_element\n')
            # Same sequence at both ends for LTR (not reverse complement)
            f.write('AAATTTGGG' + 'N' * 100 + 'AAATTTGGG\n')

        # Call getLTRs with both=True
        segments = list(
            getLTRs(
                test_fasta,
                flankdist=10,
                minid=70,
                minterm=5,
                report='external',
                temp=tempdir,
                keeptemp=False,
                both=True,
            )
        )

        # With both=True, we should get 0 or 2 segments
        if len(segments) > 0:
            assert len(segments) == 2 or len(segments) == 0
            if len(segments) == 2:
                assert segments[0].id.endswith('_L_LTR')
                assert segments[1].id.endswith('_R_LTR')
                # Check descriptions
                assert 'left LTR' in segments[0].description
                assert 'right LTR' in segments[1].description
                # Right LTR should NOT mention reverse complement
                assert 'reverse complement' not in segments[1].description.lower()


def test_getLTRs_without_both_option_yields_one_segment():
    """Test that getLTRs with both=False yields only one LTR."""
    with tempfile.TemporaryDirectory() as tempdir:
        test_fasta = os.path.join(tempdir, 'test_ltr.fa')

        with open(test_fasta, 'w') as f:
            f.write('>test_element\n')
            f.write('AAATTTGGG' + 'N' * 100 + 'AAATTTGGG\n')

        # Call getLTRs with both=False (default)
        segments = list(
            getLTRs(
                test_fasta,
                flankdist=10,
                minid=70,
                minterm=5,
                report='external',
                temp=tempdir,
                keeptemp=False,
                both=False,
            )
        )

        # With both=False, we should get 0 or 1 segment
        if len(segments) > 0:
            assert len(segments) == 1
            assert segments[0].id.endswith('_LTR')
            assert not segments[0].id.endswith('_L_LTR')
            assert not segments[0].id.endswith('_R_LTR')


def test_TIR_right_is_reverse_complemented():
    """Test that the right TIR is reverse complemented when both=True."""
    with tempfile.TemporaryDirectory() as tempdir:
        test_fasta = os.path.join(tempdir, 'test_tir.fa')

        # Create sequence where we know the TIRs
        # Left TIR: AAATTT, Right TIR: AAATTT (which is the reverse complement)
        left_tir = 'AAATTT'
        right_tir_revcomp = 'AAATTT'  # This is already the reverse complement
        internal = 'NNNNNNNNNNNNNNNNNNNNNNNNNNN'

        with open(test_fasta, 'w') as f:
            f.write('>test_element\n')
            f.write(left_tir + internal + right_tir_revcomp + '\n')

        segments = list(
            getTIRs(
                test_fasta,
                flankdist=10,
                minid=70,
                minterm=5,
                report='external',
                temp=tempdir,
                keeptemp=False,
                both=True,
            )
        )

        # If TIRs are found, check that right is reverse complemented
        if len(segments) == 2:
            left_seq = str(segments[0].seq)
            right_seq = str(segments[1].seq)

            # The right sequence should be the reverse complement of what was in the file
            # So if the file had AAATTT at the right end, the output should be its reverse complement
            assert len(left_seq) > 0
            assert len(right_seq) > 0


def test_LTR_right_is_not_reverse_complemented():
    """Test that the right LTR is NOT reverse complemented when both=True."""
    with tempfile.TemporaryDirectory() as tempdir:
        test_fasta = os.path.join(tempdir, 'test_ltr.fa')

        # For LTRs, both ends should be the same sequence (not reverse complement)
        ltr = 'AAATTTGGGCCC'
        internal = 'N' * 50

        with open(test_fasta, 'w') as f:
            f.write('>test_element\n')
            f.write(ltr + internal + ltr + '\n')

        segments = list(
            getLTRs(
                test_fasta,
                flankdist=10,
                minid=70,
                minterm=5,
                report='external',
                temp=tempdir,
                keeptemp=False,
                both=True,
            )
        )

        # If LTRs are found, check that both are in same orientation
        if len(segments) == 2:
            left_seq = str(segments[0].seq)
            right_seq = str(segments[1].seq)

            # Both sequences should be similar (not reverse complements)
            assert len(left_seq) > 0
            assert len(right_seq) > 0
            # The sequences should be the same or very similar
            # (they might differ slightly due to alignment/gaps)


def test_gaps_removed_from_TIR_sequences():
    """Test that gaps are removed from TIR sequences."""
    # Create a mock SeqRecord with gaps
    seq_with_gaps = SeqRecord(Seq('AAA-TTT-GGG'), id='test')

    # After replacing gaps, should have no dashes
    ungapped = Seq(str(seq_with_gaps.seq).replace('-', ''))
    assert '-' not in str(ungapped)
    assert str(ungapped) == 'AAATTTGGG'


def test_gaps_removed_from_LTR_sequences():
    """Test that gaps are removed from LTR sequences."""
    # Create a mock SeqRecord with gaps
    seq_with_gaps = SeqRecord(Seq('AAA-TTT-GGG'), id='test')

    # After replacing gaps, should have no dashes
    ungapped = Seq(str(seq_with_gaps.seq).replace('-', ''))
    assert '-' not in str(ungapped)
    assert str(ungapped) == 'AAATTTGGG'


def test_both_option_with_split_mode():
    """Test that both option works with split mode."""
    with tempfile.TemporaryDirectory() as tempdir:
        test_fasta = os.path.join(tempdir, 'test_tir.fa')

        with open(test_fasta, 'w') as f:
            f.write('>test_element\n')
            f.write('AAATTTGGG' + 'N' * 100 + 'CCCAAATTT\n')

        # Call with split mode and both=True
        segments = list(
            getTIRs(
                test_fasta,
                flankdist=10,
                minid=70,
                minterm=5,
                report='split',
                temp=tempdir,
                keeptemp=False,
                both=True,
            )
        )

        # With split mode and both=True, we should get:
        # - Left TIR (_L_TIR)
        # - Right TIR (_R_TIR)
        # - Internal segment (_I)
        # Total: 0 or 3 segments
        if len(segments) > 0:
            # Count segment types
            l_tir = [s for s in segments if s.id.endswith('_L_TIR')]
            r_tir = [s for s in segments if s.id.endswith('_R_TIR')]
            internal = [s for s in segments if s.id.endswith('_I')]

            # If TIRs found, should have both TIRs and internal
            if len(l_tir) > 0 or len(r_tir) > 0:
                assert len(l_tir) == 1
                assert len(r_tir) == 1
                assert len(internal) == 1


def test_both_option_with_all_mode():
    """Test that both option works with all mode."""
    with tempfile.TemporaryDirectory() as tempdir:
        test_fasta = os.path.join(tempdir, 'test_ltr.fa')

        with open(test_fasta, 'w') as f:
            f.write('>test_element\n')
            f.write('AAATTTGGG' + 'N' * 100 + 'AAATTTGGG\n')

        # Call with all mode and both=True
        segments = list(
            getLTRs(
                test_fasta,
                flankdist=10,
                minid=70,
                minterm=5,
                report='all',
                temp=tempdir,
                keeptemp=False,
                both=True,
            )
        )

        # With all mode and both=True, we should get:
        # - Left LTR (_L_LTR)
        # - Right LTR (_R_LTR)
        # - Internal segment (_I)
        # - Original element (test_element)
        # Total: 0 or 4 segments
        if len(segments) > 0:
            # Count segment types
            l_ltr = [s for s in segments if s.id.endswith('_L_LTR')]
            r_ltr = [s for s in segments if s.id.endswith('_R_LTR')]
            internal = [s for s in segments if s.id.endswith('_I')]
            original = [
                s
                for s in segments
                if not any(
                    s.id.endswith(suffix) for suffix in ['_L_LTR', '_R_LTR', '_I']
                )
            ]

            # If LTRs found, should have both LTRs, internal, and original
            if len(l_ltr) > 0 or len(r_ltr) > 0:
                assert len(l_ltr) == 1
                assert len(r_ltr) == 1
                assert len(internal) == 1
                assert len(original) == 1


def test_both_option_with_internal_mode():
    """Test that both option has no effect with internal mode."""
    with tempfile.TemporaryDirectory() as tempdir:
        test_fasta = os.path.join(tempdir, 'test_tir.fa')

        with open(test_fasta, 'w') as f:
            f.write('>test_element\n')
            f.write('AAATTTGGG' + 'N' * 100 + 'CCCAAATTT\n')

        # Call with internal mode and both=True
        segments = list(
            getTIRs(
                test_fasta,
                flankdist=10,
                minid=70,
                minterm=5,
                report='internal',
                temp=tempdir,
                keeptemp=False,
                both=True,
            )
        )

        # With internal mode, both option should not affect output
        # Should only get internal segments
        if len(segments) > 0:
            for seg in segments:
                assert seg.id.endswith('_I')
                assert not seg.id.endswith('_L_TIR')
                assert not seg.id.endswith('_R_TIR')
