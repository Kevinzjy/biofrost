"""Functions for parsing alignment results"""
import pandas as pd
__all__ = [
    "PAFParser",
    "read_paf",
    "yield_paf",
    "batch_yield_paf",
    "read_event",
    "batch_yield_event",
]

STRANDNESS = {"+": 1, "-": -1, ".": 0}


class PAFParser(object):
    """
    A unified API to parse minimap2 paf-format output like mappy.Alignment

    Parameters
    ----------
    content : list
        A paf line in list format

    Attributes
    ----------
    ctg : str
        Name of the reference sequence the query is mapped to
    ctg_len : int
        Total length of the reference genome
    r_st : int
        Start position on the reference
    r_en : int
        End position on the reference
    query : str
        Name of the query sequence
    q_len : int
        Length of the query sequence
    q_st : int
        Start position on the query
    q_en : int
        End position on the query
    strand : int
        +1 if on the forward strand; -1 if on the reverse strand; 0 if unknown
    mlen : int
        Length of the matching bases in the alignment, excluding ambiguous base matches
    blen : int
        Length of the alignment, including both alignment matches and gaps but excluding ambiguous bases
    tags : dict
        Dict of additional tag information in the SAM format
    """

    def __init__(self, content):
        self.query = content[0]
        self.q_len = int(content[1])
        self.q_st = int(content[2])
        self.q_en = int(content[3])
        self.strand = STRANDNESS[content[4]]
        self.ctg = content[5]
        self.ctg_len = int(content[6])
        self.r_st = int(content[7])
        self.r_en = int(content[8])
        self.mlen = int(content[9])
        self.blen = int(content[10])
        self.mapq = int(content[11])
        self._tag_str = content[12:]

    def __repr__(self):
        """Convert to string"""
        return "\t".join([str(x) for x in [self.q_st, self.q_en, self.strand, self.ctg, self.ctg_len, self.r_st, self.r_en, self.mlen, self.blen, self.mapq, ]])

    @property
    def tags(self):
        """Get the dictionary of tags"""
        tag_dict = {}
        for x in self._tag_str:
            tag_name, tag_type, tag_str = x[:2], x[3], x[5:]
            if tag_type == "i":
                tag_dict[tag_name] = int(tag_str)
            elif tag_type == "A":
                tag_dict[tag_name] = tag_str
            elif tag_type == "f":
                tag_dict[tag_name] = float(tag_str)
            elif tag_type == "Z":
                tag_dict[tag_name] = tag_str
            else:
                raise KeyError(f"Unknown tag format in PAF: {x}")
        return tag_dict

    def to_series(self, tags=[]):
        """Convert PAFParser to pd.Series

        Args:
            tags (list): list of keys for wanted tags

        Returns:
            pd.Series: series in ordinary PAF format
        """
        tmp_row = [self.q_len, self.q_st, self.q_en, self.strand, self.ctg, self.ctg_len, self.r_st, self.r_en, self.mlen, self.blen, self.mapq]
        tmp_cols = ["qlen", "qstart", "qend", "strand", "rname", "rlen", "rstart", "rend", "mlen", "blen", "qual"]
        for i in tags:
            tmp_cols.append(i)
            if i in self.tags:
                tmp_row.append(self.tags[i])
            else:
                tmp_row.append(pd.NA)
        return pd.Series(data=tmp_row, index=tmp_cols)



def read_paf(paf_file, index_col=None):
    """Load minimap2 paf into dataframe

    Args:
        paf_file (str): Path of minimap2 paf output
        index_col (int, optional): Index of columns to use as dataframe index. Defaults to None.

    Returns:
        pd.DataFrame: paf dataframe
    """
    paf_col = ["qname", "qlen", "qstart", "qend", "strand", "rname", "rlen", "rstart", "rend", "mlen", "blen", "qual", "tag"]

    paf_data = []
    with open(paf_file, 'r') as f:
        for line in f:
            content = line.rstrip().split("\t")
            paf_data.append(content[:12] + [content[13:], ])
    paf_data = pd.DataFrame(paf_data, columns = paf_col)
    paf_data['qlen'] = paf_data['qlen'].astype(int)
    paf_data['qstart'] = paf_data['qstart'].astype(int)
    paf_data['qend'] = paf_data['qend'].astype(int)
    paf_data['rlen'] = paf_data['rlen'].astype(int)
    paf_data['rstart'] = paf_data['rstart'].astype(int)
    paf_data['rend'] = paf_data['rend'].astype(int)
    paf_data['mlen'] = paf_data['mlen'].astype(int)
    paf_data['blen'] = paf_data['blen'].astype(int)
    if index_col is not None:
        paf_data = paf_data.set_index(paf_col[index_col])

    return paf_data


def yield_paf(paf_file):
    """Iter through minimap2 PAF output

    Args:
        paf_file (str): Path of minimap2 paf output

    Yields:
        PAFParser: PAFParser object for each hit

    Examples
    --------
    >>> #Iter through PAF
    >>> paf_file = "minimap2_output.paf"
    >>> for hit in yield_paf(paf_file):
    >>>     print(hit)
    >>>     break
    1	434	1	CECB01001999.1114.3745	2632	729	1162	433	433	0
    """
    with open(paf_file, 'r') as f:
        for line in f:
            content = line.rstrip().split("\t")
            parser = PAFParser(content)
            yield parser


def _convert_chunk(chunk, tags=[]):
    """Convert chunk to pd.DataFrame"""
    df = pd.DataFrame([hit.to_series() for hit in chunk])
    for t in tags:
        df[t] = [hit.tags[t] if t in hit.tags else pd.NA for hit in chunk]
    return df


def batch_yield_paf(paf_file, tags=[]):
    """Iter through minimap2 PAF by query name

    Args:
        paf_file (str): Path of input PAF file
        tags (list, optional): Optional tags to append to the dataframe. Defaults to [].

    Yields:
        (str, pd.DataFrame): The query sequence name and PAF formatted dataframe of hits of each query sequence

    Examples
    --------
    >>> from biofrost.parser.alignment import batch_yield_paf
    >>> for read_id, chunk in batch_yield_paf('/data/zhangjy/git/biofrost/tests/parser/test.paf', tags=["NM", "cs", ]):
    >>>     break
    >>> chunk.head()
    qlen	qstart	qend	strand	rname	rlen	rstart	rend	mlen	blen	qual	NM	cs
    0	828	1	434	1	CECB01001999.1114.3745	2632	729	1162	433	433	0	0	:433
    1	828	1	434	1	AGXZ01000038.90296.92802	2507	350	783	433	433	0	0	:433
    2	828	1	434	1	CDZI01065205.2.2427	2426	476	909	433	433	0	0	:433
    3	828	1	434	1	CEAH01070560.1.2109	2109	478	911	433	433	0	0	:433
    4	828	1	434	1	MKXP01000142.34.2826	2793	718	1151	433	433	0	0	:433
    """
    chunk = []
    read_id = None
    for hit in yield_paf(paf_file):
        if hit.query != read_id:
            if read_id is not None:
                yield read_id, _convert_chunk(chunk, tags)
            read_id = hit.query
            chunk = [hit, ]
        else:
            chunk.append(hit)
    yield read_id, _convert_chunk(chunk, tags)


# def _parse_cs(read):
#     r_st = read.reference_start
#     cs_str = dict(read.tags)['cs']
#     cs_tuple = re.findall(cs_pattern, cs_str)

#     x = r_st
#     base_alignment = defaultdict(list)
#     read_alignment = []
#     for cs in cs_tuple:
#         op, regex = cs[0], cs[1:]
#         if op == "=":  # Identical sequence (long form)
#             mlen = len(regex)
#             for i in np.arange(x, x + mlen):
#                 base_alignment[i].append("=")
#                 read_alignment.append("=")
#             x += mlen
#         elif op == ":":  # Identical sequence length
#             mlen = int(regex)
#             for i in np.arange(x, x + mlen):
#                 base_alignment[i].append("=")
#                 read_alignment.append("=")
#             x += mlen
#         elif op == "*":  # Substitution: ref to query
#             ref, alt = regex[0], regex[1]
#             base_alignment[x].append(alt)
#             read_alignment.append(alt)
#             x += 1
#         elif op == "+":  # Insertion to the reference
#             base_alignment[x].append(f"+{regex}")
#             base_alignment[x + 1].append(f"{regex}+")
#             read_alignment.append(f"+{regex}")
#         elif op == "-":  # Deletion from the reference
#             dlen = len(regex)
#             for i in np.arange(x, x + dlen):
#                 base_alignment[i].append("-")
#                 read_alignment.append("-")
#             x += dlen
#         elif op == "~":  # Intron length and splice signal
#             ilen = int(regex[2:-2])
#             x += ilen
#         else:
#             raise CSError(f"Unknown cs string: {cs}")
#     return base_alignment, read_alignment


def read_event(event_file):
    """Read nanopolish eventalign result

    Args:
        event_file (str): Path of eventalign output file

    Returns:
        pd.DataFrame: Aligned events
    """
    event_data = []
    with open(event_file, 'r') as f:
        header = f.readline().rstrip().split('\t')
        for line in f:
            content = line.rstrip().split('\t')
            event_data.append(content)
    event_data = pd.DataFrame(event_data, columns=header)
    return event_data


def batch_yield_event(event_file):
    """Iter through nanopolish eventalign results by read index

    Args:
        event_file (str): Path of eventalign output file

    Yields:
        ((str, str, str), pd.DataFrame): ((transcript_id, read_index, strand),aligned events in pd.DataFrame format)
    """
    chunk = []
    read_index = None
    with open(event_file, 'r') as f:
        header = f.readline().rstrip().split('\t')
        for line in f:
            content = line.rstrip().split('\t')
            if (content[0], content[3], content[4]) != read_index:
                if read_index is not None:
                    yield read_index, pd.DataFrame(chunk, columns=header)
                read_index = (content[0], content[3], content[4])
                chunk = [content, ]
            else:
                chunk.append(content)
        yield read_index, pd.DataFrame(chunk, columns=header)
