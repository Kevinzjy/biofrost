"""Functions for processing alignment results"""
from io import StringIO
import numpy as np
import pandas as pd
from collections import Counter
from subprocess import getstatusoutput
from ..utils import mk_random_dir


def expectation_maximization(reads, classifications, noise_threshold=0, max_iter=1000, max_delta=1e-10, verbose=False):
    """Expectation-Maximization algorithm for gene/taxon quantification
    Args:
        reads (np.array): array of query reads
        classifications (np.array): array of aligned genes / taxons
        noise_threshold (float, optional): Minimum abundance for kept during EM iteration [0-1]. Defaults to 0.
        max_iter (int, optional): Maximum EM iteration numbers. Defaults to 1000.
        max_delta (_type_, optional): Maximum detla of estimated gene / taxon abundance. Defaults to 1e-10.
        verbose (bool, optional): Verbose output. Defaults to False.
    Returns:
        (pd.Series, int): Estimated abundace and number of iteration steps
    """
    g = np.unique(classifications)
    a_i = pd.Series(data=np.zeros(g.shape) + 1/g.shape[0], index=g, dtype=np.float64)

    df = pd.DataFrame({"Read": reads, "Gene": classifications})
    df = df.drop_duplicates()

    # Seperate uniquely mapped reads
    hits = df.groupby("Read").apply(lambda x: x.shape[0])
    unique_reads = hits[hits == 1].index
    unique_genes = pd.Series(Counter(df.loc[df['Read'].isin(unique_reads)]['Gene']), dtype=np.float64)
    unique_genes.index.name = "Gene"
    unique_genes = unique_genes.reindex(g).fillna(0)

    # Assign multi-mapped reads
    ambiguous_reads = df.loc[~df['Read'].isin(unique_reads)]
    if ambiguous_reads.shape[0] == 0:
        return unique_genes, 1
    n_ambiguous = ambiguous_reads['Read'].unique().shape[0]
    n_unique = unique_reads.shape[0]

    by_read = ambiguous_reads.groupby(ambiguous_reads.index)
    by_gene = ambiguous_reads.groupby('Gene')

    n_step = 0
    d_i = 0
    trace = []
    while 1:
        # Expectation
        p = by_read.apply(lambda x: a_i[x['Gene']].sum())
        n_j = by_gene.apply(lambda x: (a_i[x.name] / p[x.index]).sum())

        # Maximization
        a_j = pd.Series(n_j, dtype=np.float64).reindex(g).fillna(0) 
        a_j = (a_j / a_j.sum() * n_ambiguous + unique_genes / unique_genes.sum() * n_unique) / (n_ambiguous + n_unique)
        a_j = a_j / a_j.sum()

        # De-noise
        a_j[a_j < noise_threshold] = 0
        a_j = a_j / a_j.sum()

        # Output Delta
        delta = (a_j - a_i).abs().sum()
        if delta <= max_delta or n_step >= max_iter:
            break

        d_j = len(str(int(1/delta)))
        if verbose and (n_step % int(max_iter/100) == 0 or d_j > d_i):
            print(f"Iter: {n_step}, delta: {delta:.12f}")

        d_i = d_j
        a_i = a_j
        n_step += 1
        trace.append(delta)

    # Final abundance
    abundance = a_j

    return abundance, n_step


def run_blastn(queries, subjects, task="blastn", blast_dir="/data/public/software/ncbi-blast-2.13.0+/bin"):
    """Run blastn of queries and subjects

    Args:
        queries (dict): dict of query sequences
        subjects (dict): dict of subject sequences
        task (str, optional): task to execut. Permissible values: blastn / blastn-short / dc-megablast / megablast / rmblastn. Defaults to "blastn".

    Raises:
        ValueError: _description_
        ValueError: _description_

    Returns:
        _type_: _description_
    """
    assert task in ['blastn', 'blastn-short', 'dc-megablast', 'megablast', 'rmblastn']

    # Create output directory
    tmp_dir = mk_random_dir()

    # Prepare input
    queries_fa = tmp_dir / "queries.fa"
    subjects_fa = tmp_dir / "subjects.fa"

    with open(queries_fa, "w") as out:
        for i, j in queries.items():
            out.write(f">{i}\n{j}\n")

    with open(subjects_fa, "w") as out:
        for i, j in subjects.items():
            out.write(f">{i}\n{j}\n")

    # Run commands
    if blast_dir is None:
        status, output = getstatusoutput("which makeblastdb")
        if status != 0:
            raise OSError("makeblastdb not found")
        status, output = getstatusoutput("which blastn")
        if status != 0:
            raise OSError("makeblastdb not found")
        blast_dir = Path(output).parent

    blast_dir = "/data/public/software/ncbi-blast-2.13.0+/bin"

    db_cmd = f"{blast_dir}/makeblastdb -in {subjects_fa} -dbtype nucl"
    run_cmd = f"{blast_dir}/blastn -query {queries_fa} -db {subjects_fa} -task {task} -outfmt 6"

    status, output = getstatusoutput(db_cmd)
    if status != 0:
        raise ValueError(f"Failed run to command: {db_cmd}\n{output}\n")

    status, output = getstatusoutput(run_cmd)
    if status != 0:
        raise ValueError(f"Failed run to command: {run_cmd}\n{output}\n")

    # Parse output
    tabular_out = pd.read_csv(StringIO(output), sep="\t", header=None)
    tabular_out.columns = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore']

    return tabular_out
