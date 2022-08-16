import os
import sys
import re
import gzip
import pandas as pd
from subprocess import getstatusoutput

from biofrost.utils import to_str




def index_annotation(gtf):
    """
    Generate binned index for element in gtf
    """
    from .utils import tree

    gtf_index = defaultdict(dict)
    intron_index = defaultdict(dict)
    splice_site_index = tree()

    last_exon = None
    with open(gtf, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            content = line.rstrip().split('\t')
            # only include gene and exon feature for now
            if content[2] not in ['gene', 'exon']:
                continue

            parser = GTFParser(content)

            # Extract splice site
            if content[2] == 'exon':
                splice_site_index[parser.contig][parser.start][parser.strand]['start'] = 1
                splice_site_index[parser.contig][parser.end][parser.strand]['end'] = 1

                # Load intron
                if last_exon is not None and last_exon.attr['transcript_id'] == parser.attr['transcript_id']:
                    intron_start = last_exon.end if last_exon.strand == '+' else last_exon.start
                    intron_end = parser.start if parser.strand == '+' else parser.end
                    intron_strand = parser.strand

                    intron_start, intron_end = min(intron_start, intron_end), max(intron_start, intron_end)
                    start_div, end_div = intron_start // 500, intron_end // 500
                    for i in range(start_div, end_div + 1):
                        intron_index[parser.contig].setdefault(i, []).append((intron_start, intron_end, intron_strand))

                last_exon = parser

            # Binned index
            start_div, end_div = parser.start // 500, parser.end // 500
            for i in range(start_div, end_div + 1):
                gtf_index[parser.contig].setdefault(i, []).append(parser)

    return gtf_index, intron_index, splice_site_index


def revcomp(seq):
    """
    Convert sequence to reverse complementary
    """
    trantab = str.maketrans("ATCG", "TAGC")
    return seq.translate(trantab)[::-1]


def get_bsj(seq, bsj):
    """Return transformed sequence of given BSJ"""
    return seq[bsj:] + seq[:bsj]


def get_n50(sequence_lengths):
    """
    Get n50 of sequence lengths
    """
    sequence_lengths = sorted(sequence_lengths, reverse=True)
    total_bases = sum(sequence_lengths)
    target_bases = total_bases * 0.5
    bases_so_far = 0
    for sequence_length in sequence_lengths:
        bases_so_far += sequence_length
        if bases_so_far >= target_bases:
            return sequence_length
    return 0

def get_mm_exons(hit):
    """
    Get blocks of aligned segments
    :param hit:
    :return:
    """
    r_start = hit.r_st
    r_end = hit.r_st
    r_block = []
    for length, operation in hit.cigar:
        if operation == 0:
            r_end += length
        elif operation == 1:
            pass
        elif operation == 2:
            r_end += length
        elif operation == 3:
            r_block.append([r_start, r_end, r_end - r_start + 1])
            r_start = r_end + length
            r_end = r_start
        elif operation == 4:
            pass
    if r_end > r_start:
        r_block.append([r_start, r_end, r_end - r_start + 1])
    return r_block


def get_exons(hit):
    """
    Get exons from pysam alignment results
    """
    r_start, r_end = hit.reference_start, hit.reference_start
    q_start, q_end = hit.query_alignment_start, hit.query_alignment_start

    r_block = []
    for operation, length in hit.cigar:
        if operation == 0:
            r_end += length
            q_end += length
        elif operation == 1:
            q_end += length
        elif operation == 2:
            r_end += length
        elif operation == 3:
            r_block.append([r_start, r_end, q_start, q_end])
            r_start = r_end + length
            r_end = r_start
            q_start = q_end
        elif operation == 4:
            pass

    if r_end > r_start:
        r_block.append([r_start, r_end, q_start, q_end])
    else:
        pass

    return r_block





def iter_mm_paf(paf_file, is_cigar=False):
    with open(paf_file, 'r') as f:
        last_id, last_len = None, None
        hits = []
        for line in f:
            content = line.rstrip().split('\t')
            qname, rname = content[0], content[5] 
            qlen, rlen, rstart, rend = int(content[1]), int(content[6]), int(content[7]), int(content[8])
            cigar = content[-1].split(":")[-1] if is_cigar else None
            if qname != last_id:
                if last_id is not None:
                    yield last_id, last_len, hits
                last_id, last_len = qname, qlen
                hits = []
            hits.append([rname, rstart, rend, cigar])
        yield last_id, last_len, hits


def load_idr_bed(infile, as_pyranges=False):
    import pyranges as pr
    idr_header = [
        'chrom', 'chromStart', 'chromEnd', 'name', 'score', 'strand', 'signalValue', 'p-value', 'q-value', 'summit', 'localIDR', 'globalIDR', 
        'rep1_chromStart', 'rep1_chromEnd', 'rep1_signalValue', 'rep1_summit',
        'rep2_chromStart', 'rep2_chromEnd', 'rep2_signalValue', 'rep2_summit',
    ]
    df = pd.read_csv(infile, sep="\t", header=None)
    if df.shape[0] <= 20:
        return None

    df.columns = idr_header
    df['IDR_score'] = 2 ** (df['score'] / -125)

    if as_pyranges:
        peaks = pr.PyRanges(df.rename({"chrom": "Chromosome", "chromStart": "Start", "chromEnd": "End", "strand": "Strand"}, axis=1))
        return peaks

    return df


def load_narrowPeak(infile, as_pyranges=False):
    import pyranges as pr
    narrowPeak_header = [
        'chrom', 'chromStart', 'chromEnd', 'name', 'score', 'strand', 'signalValue', 'pValue', 'qvalue', 'peak'
    ]
    df = pd.read_csv(infile, sep="\t", header=None)
    df.columns = narrowPeak_header
    
    if as_pyranges:
        peaks = pr.PyRanges(df.rename({"chrom": "Chromosome", "chromStart": "Start", "chromEnd": "End", "strand": "Strand"}, axis=1))
        return peaks

    return df