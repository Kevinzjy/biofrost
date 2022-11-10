"""Functions for processing alignment results"""
import numpy as np
import pandas as pd
from collections import Counter


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

    # Seperate uniquely mapped reads
    hits = df.groupby("Read").apply(lambda x: x.shape[0])
    unique_reads = hits[hits == 1].index
    unique_genes = pd.Series(Counter(df.loc[df.index.isin(unique_reads)]['Gene']), dtype=np.float64)
    unique_genes.index.name = "Gene"
    unique_genes = unique_genes.reindex(g).fillna(0)

    # Assign multi-mapped reads
    ambiguous_reads = df.loc[~df.index.isin(unique_reads)]
    if ambiguous_reads.shape[0] == 0:
        return unique_genes, 1

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
        a_j = pd.Series(n_j, dtype=np.float64).reindex(g).fillna(0) + unique_genes
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