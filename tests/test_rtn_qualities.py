"""rtn_qualities.py: base quality 2 at chosen positions for rtn, original qualities back afterwards."""
import rtn_qualities as rq


def sam(pos, cigar, seq, qual=None, flag=0, name='r', mapq='60'):
    return '\t'.join([name, str(flag), 'chrM', str(pos), mapq, cigar, '*', '0', '0', seq, qual or 'I' * len(seq)]) + '\n'


def quals(line):
    return line.split('\t')[10].rstrip('\n')


def test_lowered_where_the_read_is_aligned_to_a_position():
    assert quals(rq.lower(sam(100, '10M', 'A' * 10), {102, 109})) == 'II#IIIIII#'


def test_insertions_and_soft_clips_move_the_read_not_the_reference():
    # 2 soft-clipped, 3 aligned (100-102), 1 inserted, 4 aligned (103-106)
    assert quals(rq.lower(sam(100, '2S3M1I4M', 'A' * 10), {102, 103})) == 'IIII#I#III'


def test_a_deletion_skips_its_reference_position():
    # 100-101 aligned, 102 deleted, 103-104 aligned
    assert quals(rq.lower(sam(100, '2M1D2M', 'A' * 4), {102, 103})) == 'II#I'


def test_the_appended_copy_of_the_start_counts_as_the_start():
    assert quals(rq.lower(sam(16568, '4M', 'A' * 4), {1})) == 'II#I'


def test_unmapped_reads_and_headers_stay():
    assert rq.lower('@HD\tVN:1.6\n', {1}) == '@HD\tVN:1.6\n'
    line = sam(100, '4M', 'AAAA', flag=4)
    assert rq.lower(line, {101}) == line


def test_restore_gives_back_the_qualities_and_keeps_rtn_mapq():
    orig = {('r', 0, '100'): 'ABCD'}
    out = rq.restore(sam(100, '4M', 'AAAA', qual='A#C#', mapq='0'), orig)
    assert quals(out) == 'ABCD' and out.split('\t')[4] == '0'


def test_restore_reports_a_record_without_its_original():
    assert rq.restore(sam(101, '4M', 'AAAA'), {('r', 0, '100'): 'ABCD'}) is None
