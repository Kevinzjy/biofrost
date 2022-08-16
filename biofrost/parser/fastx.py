"""Functions for parsing FASTA/FASTX format"""
import gzip
from collections import namedtuple

import pysam

from biofrost.utils import to_str, is_gzipped

__all__ = [
    "Faidx",
    "is_fastq",
    "read_fastx",
    "yield_fastx",
]

Seq = namedtuple("Seq", "header seq qual")


class Faidx(object):
    """
    Convert pysam.FastaFile to unify API like mappy::Aligner

    Attributes
    ----------
    faidx : pysam.FastaFile
        fasta file handle
    contig_len : dict
        length of each contig

    Methods
    -------
    info(additional=""):
        Prints the person's name and age.

    """

    def __init__(self, infile):
        """Load contig information from faidx"""
        import pysam

        if not os.path.exists(infile + '.fai'):
            sys.stderr.write('Generating faidx for input sequences ...')
        try:
            faidx = pysam.FastaFile(infile)
        except ValueError:
            raise ValueError('Cannot load Faidx, index file is missing')
        except IOError:
            raise IOError('Cannot generate Faidx, file could not be opened')

        self.faidx = faidx
        self.contig_len = {contig: faidx.get_reference_length(contig) for contig in faidx.references}

    def seq(self, contig, start, end):
        """Return sequence of given coordinate"""
        return self.faidx.fetch(contig, start, end)

    def close(self):
        self.faidx.close()


def is_fastq(f_name, is_gz=False):
    """Determine is file have quality line (FASTQ) or not (FASTA)

    Args:
        f_name (str): Path of FASTA/FASTQ file
        is_gz (bool, optional): Whether input file is gzipped. Defaults to False.

    Raises:
        KeyError: If the first character is not "@" or ">"

    Returns:
        bool: True if file is in FASTQ format, otherwise False
    """
    f = gzip.open(f_name, 'rb') if is_gz else open(f_name, 'r')
    line = to_str(f.readline()).rstrip()

    is_qual = False
    if line.startswith('>'):
        is_qual = False
    elif line.startswith('@'):
        is_qual = True
    else:
        raise KeyError(f"Cannot determine FASTA/FASTQ format for {f_name} line {line}")

    f.close()
    return is_qual


def _read_fasta(fname, is_gz=False):
    sequences = {}
    seq_header = None
    seq_id = None
    seq = ''
    f = gzip.open(fname, 'rb') if is_gz else open(fname, 'r')
    for line in f:
        context = to_str(line).rstrip()
        if context.startswith('>'):
            if seq_id is not None:
                sequences[seq_id] = Seq(seq_header, seq, None)
            seq_id = context.lstrip('>').split()[0]
            seq_header = context.split()[1:]
            seq = ''
        else:
            seq += context
    sequences[seq_id] = Seq(seq_header, seq, None)
    f.close()
    return sequences


def _read_fastq(fname, is_gz=False):
    sequences = {}
    f = gzip.open(fname, 'rb') if is_gz else open(fname, 'r')
    for line in f:
        header = to_str(line).rstrip().lstrip('@')
        seq_id = header.split()[0]
        seq = to_str(f.readline()).rstrip()
        sep = to_str(f.readline()).rstrip()
        qual = to_str(f.readline()).rstrip()
        sequences[seq_id] = Seq(header.split()[1:], seq, qual)
    f.close()
    return sequences


def read_fastx(fa_file):
    """Read FASTA/Q file into dict

    Args:
        fa_file (str): Path of input FASTA/Q

    Returns:
        dict: dict of seq_id -> Seq(seq, qual)
    """
    is_gz = is_gzipped(fa_file)
    is_fq = is_fastq(fa_file, is_gz)
    if is_fq:
        return _read_fastq(fa_file, is_gz)
    else:
        return _read_fasta(fa_file, is_gz)


def _yield_fasta(fname, is_gz=False):
    seq_header = None
    seq_id = None
    seq = ''
    f = gzip.open(fname, 'rb') if is_gz else open(fname, 'r')
    for line in f:
        context = to_str(line).rstrip()
        if context.startswith('>'):
            if seq_id is not None:
                yield seq_id, Seq(seq_header, seq, None)
            seq_id = context.lstrip('>').split()[0]
            seq_header = context.split()[1:]
            seq = ''
        else:
            seq += context
    yield seq_id, Seq(seq_header, seq, None)
    f.close()


def _yield_fastq(fname, is_gz=False):
    f = gzip.open(fname, 'rb') if is_gz else open(fname, 'r')
    for line in f:
        header = to_str(line).rstrip().lstrip('@')
        seq_id = header.split()[0]
        seq = to_str(f.readline()).rstrip()
        sep = to_str(f.readline()).rstrip()
        qual = to_str(f.readline()).rstrip()
        yield seq_id, Seq(header.split()[1:], seq, qual)
    f.close()


def yield_fastx(fa_file):
    """Iter through FASTA/FASTQ records

    Args:
        fa_file (str): Path of input FASTA/FASTQ

    Yields:
        tuple: (read_id, Seq(header, seq, qual))

    Examples
    --------
    >>> from biofrost.parser import yield_fastx
    >>> for seq_id, seq in yield_fastx(test_dir / "parser" / "test.fa"):
    >>>     print(seq_id, seq.seq)
    seq1 AAAAATTTTT
    seq2 ACGATCGAC
    """
    is_gz = is_gzipped(fa_file)
    is_fq = is_fastq(fa_file, is_gz)
    if is_fq:
        yield from _yield_fastq(fa_file, is_gz)
    else:
        yield from _yield_fasta(fa_file, is_gz)
