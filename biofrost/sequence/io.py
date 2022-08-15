import gzip
import pysam
from biofrost.utils import to_str


class Faidx(object):
    """
    Convert pysam.FastaFile to unify API like mappy2::Aligner

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


class Fasta(object):
    """
    load fasta file into memory
    """
    def __init__(self, infile):
        faidx = Faidx(infile)
        self.contig_len = faidx.contig_len
        self.genome = {ctg: faidx.seq(ctg, 0, size) for ctg, size in self.contig_len.items()}
        faidx.close()

    def seq(self, contig, start, end):
        if contig not in self.genome:
            return None
        return self.genome[contig][start:end]


def load_fasta(fname, is_gz=False):
    sequences = {}
    seq_id = None
    seq = ''
    f = gzip.open(fname, 'rb') if is_gz else open(fname, 'r')
    for line in f:
        if to_str(line).startswith('>'):
            if seq_id is not None:
                sequences[seq_id] = seq
            seq_id = to_str(line).rstrip().split()[0].lstrip('>')
            seq = ''
        else:
            seq += to_str(line).rstrip()
    sequences[seq_id] = seq
    f.close()
    return sequences


def yield_fasta(fname, is_gz=False):
    sequences = {}
    seq_id = None
    seq = ''
    f = gzip.open(fname, 'rb') if is_gz else open(fname, 'r')
    for line in f:
        if to_str(line).startswith('>'):
            if seq_id is not None:
                yield (seq_id, seq)
            seq_id = to_str(line).rstrip().lstrip('>')
            seq = ''
        else:
            seq += to_str(line).rstrip()
    yield (seq_id, seq)
    f.close()


def load_fastq(fname, is_gz=False):
    sequences = {}
    seq_id = None
    seq = ''
    f = gzip.open(fname, 'rb') if is_gz else open(fname, 'r')
    for line in f:
        seq_id = to_str(line).rstrip().lstrip('@').split(' ')[0]
        seq = to_str(f.readline()).rstrip()
        sep = to_str(f.readline()).rstrip()
        qual = to_str(f.readline()).rstrip()
        sequences[seq_id] = (seq, qual)
    f.close()
    return sequences


def yield_fastq(fname, is_gz=False):
    f = gzip.open(fname, 'rb') if is_gz else open(fname, 'r')
    for line in f:
        read_id = to_str(line).rstrip().lstrip('@')
        seq = to_str(f.readline()).rstrip()
        sep = to_str(f.readline()).rstrip()
        qual = to_str(f.readline()).rstrip()
        yield (read_id, seq, sep, qual)
    f.close()