"""mark_removed_reads.py: every mapped read of p09 and p10, the removed ones at MAPQ 0, each with its reason in ZF."""
import mark_removed_reads as mr


def rec(name, pos=100, cigar='10M', mapq='60', flag=0, tags=()):
    return [name, str(flag), 'chrM', str(pos), mapq, cigar, '*', '0', '0', 'A' * 10, 'I' * 10, *tags]


def line(f):
    return '\t'.join(f) + '\n'


def zf(out_line):
    f = out_line.rstrip('\n').split('\t')
    return f[4], f[-1]


def run_p09(inputs, rtn, clipped=True):
    rtn_records = {mr.key(f): f for f in rtn}
    out, removed = mr.p09(['@HD\tVN:1.6\n'] + [line(f) for f in inputs], rtn_records, 30, clipped)
    return out, removed


def test_p09_kept_read_keeps_its_mapq_and_rtn_tags():
    out, removed = run_p09([rec('a')], [rec('a', tags=('ZH:i:1', 'ZN:i:2'))])
    assert out[0] == '@HD\tVN:1.6\n'
    assert out[1].rstrip('\n').split('\t')[11:] == ['ZH:i:1', 'ZN:i:2', 'ZF:Z:kept']
    assert zf(out[1]) == ('60', 'ZF:Z:kept') and removed == {}


def test_p09_rtn_mapq_0_is_numt():
    out, removed = run_p09([rec('a')], [rec('a', mapq='0')])
    assert zf(out[1]) == ('0', 'ZF:Z:numt') and removed == {'a': 'numt'}


def test_p09_rtn_setting_a_low_bwa_mapq_to_0_is_numt():
    out, removed = run_p09([rec('a', mapq='10')], [rec('a', mapq='0')])
    assert zf(out[1]) == ('0', 'ZF:Z:numt') and removed == {'a': 'numt'}


def test_p09_low_bwa_mapq_that_rtn_left_alone_or_out_is_mapq():
    out, removed = run_p09([rec('a', mapq='10'), rec('b', mapq='2')], [rec('a', mapq='10')])
    assert [zf(x) for x in out[1:]] == [('0', 'ZF:Z:mapq'), ('0', 'ZF:Z:mapq')] and removed == {'a': 'mapq', 'b': 'mapq'}


def test_p09_clipped_read_only_with_the_clip_rule():
    inp, rtn = [rec('a', cigar='3S7M')], [rec('a', cigar='3S7M')]
    assert zf(run_p09(inp, rtn)[0][1]) == ('0', 'ZF:Z:clipped')
    assert zf(run_p09(inp, rtn, clipped=False)[0][1]) == ('60', 'ZF:Z:kept')


def test_p09_unmapped_read_is_left_out_but_its_reason_goes_to_p10():
    out, removed = run_p09([rec('a', pos=0, cigar='*', mapq='0', flag=4)], [])
    assert out == ['@HD\tVN:1.6\n'] and removed == {'a': 'unmapped'}


def test_p09_record_rtn_left_out_takes_its_primary_reason():
    primary = rec('a', cigar='4S6M')
    supplementary = rec('a', pos=5000, cigar='6H4M', flag=2048)
    out, removed = run_p09([primary, supplementary], [rec('a', cigar='4S6M')])
    assert [zf(x)[1] for x in out[1:]] == ['ZF:Z:clipped', 'ZF:Z:clipped'] and removed == {'a': 'clipped'}


def test_p09_record_rtn_left_out_of_a_kept_read_is_numt():
    out, _ = run_p09([rec('a'), rec('a', pos=5000, cigar='6H4M', flag=2048)], [rec('a')])
    assert [zf(x)[1] for x in out[1:]] == ['ZF:Z:kept', 'ZF:Z:numt']


def test_p09_primary_rtn_left_out_is_numt():
    out, removed = run_p09([rec('a')], [])
    assert zf(out[1]) == ('0', 'ZF:Z:numt') and removed == {'a': 'numt'}


def test_p10_reasons():
    kept, removed = {'k', 'low'}, {'n': 'numt', 'c': 'clipped', 'u': 'unmapped'}
    reason = lambda f: mr.p10_reason(f, kept, removed, 30)
    assert reason(rec('k')) == 'kept'
    assert reason(rec('low', mapq='12')) == 'mapq'
    assert reason(rec('n')) == 'numt'
    assert reason(rec('c')) == 'clipped'
    assert reason(rec('t')) == 'trimmed'
    assert reason(rec('u')) == 'unmapped'


def test_mark_sets_mapq_0_only_for_removed_reads():
    assert zf(mr.mark(rec('a'), 'kept')) == ('60', 'ZF:Z:kept')
    assert zf(mr.mark(rec('a'), 'trimmed')) == ('0', 'ZF:Z:trimmed')
