"""Collections of useful analysis functions & algorithms"""
from .alignment import expectation_maximization
from .bioinfo import get_n50
from .de import calcNormFactors, calcFactorTMM, calcNormCPM

__all__ = [
    "expectation_maximization",
    "get_n50",
    "calcNormFactors",
    "calcFactorTMM",
    "calcNormCPM"
]