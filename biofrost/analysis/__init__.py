"""Collections of useful analysis functions & algorithms"""
from .alignment import expectation_maximization, run_blastn
from .bioinfo import get_n50
from .de import calcNormFactors, calcFactorTMM, calcNormCPM
from .ontology import clusterProfiler

__all__ = [
    "expectation_maximization",
    "run_blastn",

    "get_n50",

    "calcNormFactors",
    "calcFactorTMM",
    "calcNormCPM",

    "clusterProfiler",
]