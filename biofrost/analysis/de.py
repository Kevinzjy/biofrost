import numpy as np
from numpy import errstate

def calcNormCPM(count_df):
    """Generate normalized cpm dataframe

    Args:
        count_df (pd.DataFrame): Raw read count dataframe

    Returns:
        pd.DataFrame: Normalized CPM dataframe
    """
    factor = calcNormFactors(count_df)
    cpm_df = count_df.apply(lambda x: 1000 * 1000 * x / (x.sum() * factor[x.name]), axis=0)
    return cpm_df

def calcNormFactors(df):
    """Calculate TMM normalization factors of given read count matrix

    Args:
        df (pd.DataFrame): read count matrix, where rows are genes and columns are samples

    Returns:
        pd.Series: library size factor of each sample
    """
    f75 = df.apply(lambda x: np.quantile(x, .75), axis=0)
    ref_idx = abs(f75-np.mean(f75)).idxmin()

    f = df.apply(lambda x: calcFactorTMM(x, df[ref_idx]), axis=0)

    f = f / np.exp(np.mean(np.log(f)))
    return f


def calcFactorTMM(obs, ref, libsize_obs=None, libsize_ref=None, logratioTrim=.3, sumTrim=0.05, doWeighting=True, Acutoff=-1e10):
    """Calculate TMM factor for given pair of samples

    Args:
        obs (pd.Seires): observation sample
        ref (pd.Series): reference sample
        libsize_obs (int, optional): Library size of observation sample. Defaults to sum of read counts.
        libsize_ref (int, optional): LIbrary size of reference sample. Defaults to sum of read counts.
        logratioTrim (float, optional): _description_. Defaults to .3.
        sumTrim (float, optional): _description_. Defaults to 0.05.
        doWeighting (bool, optional): _description_. Defaults to True.
        Acutoff (_type_, optional): _description_. Defaults to -1e10.

    Returns:
        _type_: _description_
    """
    # TMM between two libraries
    obs = obs.astype(float)
    ref = ref.astype(float)

    nO = obs.sum() if libsize_obs is None else libsize_obs
    nR = ref.sum() if libsize_ref is None else libsize_ref

    with errstate(divide='ignore'):
        logR = np.log2((obs/nO)/(ref/nR))               # log ratio of expression, accounting for library size
        absE = (np.log2(obs/nO) + np.log2(ref/nR)) / 2  # absolute expression
        v = (nO-obs)/nO/obs + (nR-ref)/nR/ref           # estimated asymptotic variance

    # remove infinite values, cutoff based on A
    fin = np.isfinite(logR) & np.isfinite(absE) & (absE > Acutoff)

    logR = logR[fin]
    absE = absE[fin]
    v = v[fin]

    # If difference is too small
    # if max(abs(logR)) < 1e-6:
    #     return 1

    # taken from the original mean() function
    n = logR.shape[0]
    loL = np.floor(n * logratioTrim) + 1
    hiL = n + 1 - loL
    loS = np.floor(n * sumTrim) + 1
    hiS = n + 1 - loS

    # keep <- (rank(logR) %in% loL:hiL) & (rank(absE) %in% loS:hiS)
    # a fix from leonardo ivan almonacid cardenas, since rank() can return
    # non-integer values when there are a lot of ties
    keep = set(logR.sort_values().index[int(loL):int(hiL)]) & set(absE.sort_values().index[int(loS):int(hiS)])

    if doWeighting:
        f = sum(logR[keep]/v[keep]) / sum(1/v[keep])
    else:
        f = np.mean(logR[keep])

    if np.isnan(f):
        f = 0

    return 2 ** f