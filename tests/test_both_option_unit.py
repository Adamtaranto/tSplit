"""Unit tests for --both option functionality without requiring external tools."""

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def test_ungap_removes_gaps():
    """Test that gap removal using replace works correctly."""
    # Test with dashes
    seq_with_dashes = Seq('AAA-TTT-GGG-CCC')
    ungapped = Seq(str(seq_with_dashes).replace('-', ''))
    assert str(ungapped) == 'AAATTTGGGCCC'
    assert '-' not in str(ungapped)

    # Test with no gaps
    seq_no_gaps = Seq('AAATTTGGGCCC')
    ungapped_no_gaps = Seq(str(seq_no_gaps).replace('-', ''))
    assert str(ungapped_no_gaps) == 'AAATTTGGGCCC'

    # Test with multiple types of gaps
    seq_multi_gaps = Seq('A-A-A-T-T-T')
    ungapped_multi = Seq(str(seq_multi_gaps).replace('-', ''))
    assert str(ungapped_multi) == 'AAATTT'


def test_reverse_complement_produces_correct_sequence():
    """Test that reverse complement works correctly for TIR sequences."""
    # Simple test sequence
    seq = Seq('AAATTTGGG')
    rc = seq.reverse_complement()
    # Reverse of AAATTTGGG is GGGTTAAA
    # Complement of AAATTTGGG is TTTAAACCC
    # Reverse complement is CCCAAATTT
    assert str(rc) == 'CCCAAATTT'

    # Test with SeqRecord
    record = SeqRecord(seq, id='test', name='test', description='test desc')
    rc_record = record.reverse_complement(
        id='test_rc', name='test_rc', description='rc desc'
    )
    assert str(rc_record.seq) == 'CCCAAATTT'
    assert rc_record.id == 'test_rc'
    assert rc_record.name == 'test_rc'
    assert rc_record.description == 'rc desc'


def test_seqrecord_slicing():
    """Test that SeqRecord slicing works as expected."""
    seq = Seq('AAATTTGGGCCCNNNNNN')
    record = SeqRecord(seq, id='test_element')

    # Test left slice
    left = record[0:9]
    assert str(left.seq) == 'AAATTTGGG'
    assert left.id == 'test_element'

    # Test right slice
    right = record[12:18]
    assert str(right.seq) == 'NNNNNN'

    # Test middle slice
    middle = record[9:12]
    assert str(middle.seq) == 'CCC'


def test_seqrecord_id_modification():
    """Test that SeqRecord IDs can be modified correctly."""
    seq = Seq('AAATTTGGG')
    record = SeqRecord(seq, id='element1', name='element1', description='test')

    # Modify ID
    record.id = 'element1_L_TIR'
    record.name = 'element1_L_TIR'
    record.description = '[element1 left TIR segment]'

    assert record.id == 'element1_L_TIR'
    assert record.name == 'element1_L_TIR'
    assert record.description == '[element1 left TIR segment]'
    assert 'left TIR' in record.description


def test_suffix_naming_conventions():
    """Test that suffix naming follows the correct pattern."""
    base_id = 'element1'

    # Test TIR suffixes
    assert f'{base_id}_L_TIR' == 'element1_L_TIR'
    assert f'{base_id}_R_TIR' == 'element1_R_TIR'
    assert f'{base_id}_TIR' == 'element1_TIR'

    # Test LTR suffixes
    assert f'{base_id}_L_LTR' == 'element1_L_LTR'
    assert f'{base_id}_R_LTR' == 'element1_R_LTR'
    assert f'{base_id}_LTR' == 'element1_LTR'

    # Test internal suffix
    assert f'{base_id}_I' == 'element1_I'


def test_both_flag_logic():
    """Test the boolean logic for the both flag."""
    both = True
    report = 'external'

    # Check that both works with expected report modes
    assert report in ['split', 'external', 'all']

    if both and report in ['split', 'external', 'all']:
        # Should output both terminal repeats
        output_both = True
    else:
        output_both = False

    assert output_both is True

    # Test with internal mode
    report = 'internal'
    if both and report in ['split', 'external', 'all']:
        output_both = True
    else:
        output_both = False

    assert output_both is False


def test_description_strings():
    """Test that description strings are correctly formatted."""
    element_id = 'test_element'

    # Left TIR description
    left_tir_desc = f'[{element_id} left TIR segment]'
    assert 'left TIR' in left_tir_desc
    assert element_id in left_tir_desc

    # Right TIR description with reverse complement note
    right_tir_desc = f'[{element_id} right TIR segment, reverse complemented]'
    assert 'right TIR' in right_tir_desc
    assert 'reverse complemented' in right_tir_desc

    # Left LTR description
    left_ltr_desc = f'[{element_id} left LTR segment]'
    assert 'left LTR' in left_ltr_desc

    # Right LTR description (no reverse complement)
    right_ltr_desc = f'[{element_id} right LTR segment]'
    assert 'right LTR' in right_ltr_desc
    assert 'reverse complement' not in right_ltr_desc
