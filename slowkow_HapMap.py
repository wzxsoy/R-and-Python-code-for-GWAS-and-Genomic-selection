#!/usr/bin/env python
"""
HapMap.py
Kamil Slowikowski
February 25, 2014

# Quickly read the data and automatically compute allele frequencies.
    >>> from HapMap import HapMap
    >>> ceu = HapMap("CEU", "/media/blue/data/EPI511/HapMap3/")

# Show a quick view of the data.
    >>> ceu

    CEU pop 718848 SNPs 112 individuals 57/55 M/F

    /media/blue/data/EPI511/HapMap3/CEU.ind
    ind sex  pop
    0  NA06989   F  CEU
    1  NA11891   M  CEU
    2  NA11843   M  CEU
    3  NA12341   F  CEU
    4  NA06984   M  CEU

    /media/blue/data/EPI511/HapMap3/CEU.geno
    ind         NA06989  NA11891  NA11843  NA12341  NA06984
    snp
    rs3131972         1        1        2        1        2
    rs3131969         1        2        2        1        2
    rs3131967         1        2        2        1        2
    rs1048488         1        1        2        1        2
    rs12562034        1        2        2        2        1

    /media/blue/data/EPI511/HapMap3/CEU.snp
    snp  chrom  cM     pos ref alt
    0   rs3131972      1   0  742584   G   A
    1   rs3131969      1   0  744045   G   A
    2   rs3131967      1   0  744197   C   T
    3   rs1048488      1   0  750775   T   C
    4  rs12562034      1   0  758311   G   A

# geno, ind, and snp are DataFrame objects.
    >>> ceu.snp

    <class 'pandas.core.frame.DataFrame'>
    Int64Index: 718848 entries, 0 to 718847
    Data columns (total 6 columns):
    snp      718848  non-null values
    chrom    718848  non-null values
    cM       718848  non-null values
    pos      718848  non-null values
    ref      718848  non-null values
    alt      718848  non-null values
    dtypes: float64(1), int64(2), object(3)

# Get the first 10 rows of SNP information.
    >>> ceu.snp.head()

              snp  chrom  cM     pos ref alt
    0   rs3131972      1   0  742584   G   A
    1   rs3131969      1   0  744045   G   A
    2   rs3131967      1   0  744197   C   T
    3   rs1048488      1   0  750775   T   C
    4  rs12562034      1   0  758311   G   A

# Minor allele frequencies.
    >>> ceu.maf.head()

    snp
    rs3131972     0.165179
    rs3131969     0.133929
    rs3131967     0.117117
    rs1048488     0.166667
    rs12562034    0.102679
    dtype: float64

# Reference alleles.
    >>> ceu.snp.ref.head()

    0    G
    1    G
    2    C
    3    T
    4    G
    Name: ref, dtype: object

# Count males and females.
    >>> ceu.ind.sex.value_counts()

    M    57
    F    55
    dtype: int64

# Mean genotype for SNP 42.
    >>> ceu.geno.ix[42].mean()
    1.1517857142857142

# Allele counts of SNP 42.
    >>> ceu.geno.ix[42].value_counts()

    1    55
    2    37
    0    20
    dtype: int64

# Genotypes for SNP 42.
    >>> ceu.geno.ix[42]

    ind
    NA06989    2
    NA11891    0
    NA11843    1
    ...
    NA12815    0
    NA12043    2
    NA12264    0
    Name: rs9442373, Length: 112, dtype: float64

# Name of SNP 42.
    >>> ceu.geno.index[42]
    'rs9442373'

# Get the position of SNP 42.
    >>> ceu.snp.ix[42]['pos']
    1052501

# Select a slice of rows, and columns by name.
    >>> ceu.snp.ix[42:52, ['chrom', 'pos']]

        chrom      pos
    42      1  1052501
    43      1  1084601
    44      1  1089205
    45      1  1096336
    46      1  1096647
    47      1  1109721
    48      1  1120069
    49      1  1125105
    50      1  1145994
    51      1  1153667
    52      1  1163474
"""


# if future: then 1 / 3 == 0.33 else 1 / 3 == 0
from __future__ import division


import os
import sys
from itertools import combinations, izip, tee, chain
from functools import partial
import operator as op
from copy import deepcopy
from time import time, asctime
from contextlib import contextmanager

import pandas as pd
import numpy as np


class HapMap(object):
    """
    A class to wrap around data stored in the following formats:

    ==> HapMap3/ASW.geno <==
    0210121000012112112002101111101011110100101111012
    0210121000012112122002111111101911110100102111022
    0210121001012112122002111111101119110200102111022
    0210121000022112112002111211101111120100202111022
    2122212222222212222222222222221212122222221122222
    1210121000022112122102111211101111121120212111022
    0112112210211211222111110222122111021220221222011
    1102221112121201222220201111211112122112112111111
    0222111111222121110202221012211121101221210122111
    2211121222110211122120121211112222122211012101221

    ==> HapMap3/ASW.ind <==
                 NA19916 M        ASW
                 NA19835 F        ASW
                 NA20282 F        ASW
                 NA19703 M        ASW
                 NA19901 F        ASW
                 NA19908 M        ASW
                 NA19914 F        ASW
                 NA20287 F        ASW
                 NA19713 F        ASW
                 NA19904 M        ASW

    ==> HapMap3/ASW.snp <==
                          rs3131972     1    0.0     742584 G A
                          rs3131969     1    0.0     744045 G A
                          rs3131967     1    0.0     744197 C T
                          rs1048488     1    0.0     750775 T C
                         rs12562034     1    0.0     758311 G A
                          rs4040617     1    0.0     769185 A G
                          rs4970383     1    0.0     828418 C A
                          rs4475691     1    0.0     836671 C T
                          rs1806509     1    0.0     843817 C A
                          rs7537756     1    0.0     844113 A G
    """
    def __init__(self, pop, root="."):
        """Create a HapMap object.

        Parameters
        ----------
        pop     Prefix of the .geno .ind .snp files like "CEU".
        root    Directory with the .geno .ind .snp files.
        """
        self.pop = pop
        # File names.
        self.ind_file = os.path.join(root, pop + ".ind")
        self.snp_file = os.path.join(root, pop + ".snp")
        self.geno_file = os.path.join(root, pop + ".geno")
        # DataFrames of individuals, snps, and genotypes.
        self.ind = self._ind_dataframe(self.ind_file)
        self.snp = self._snp_dataframe(self.snp_file)
        self.geno = self._geno_dataframe(self.geno_file)
        # Compute reference and minor allele frequencies.
        self.compute_freqs()
        # Convenient property access to indices for each chromosome.
        for c in self.snp.chrom.drop_duplicates():
            idx = self.snp.index[self.snp.chrom == c]
            setattr(self, 'chr{}'.format(c), idx)

    def chrom(self, c):
        """Convenience method to get indices for a single chromosome."""
        return self.snp.index[self.snp.chrom == c]

    def idxsnp(self, snp):
        """Get the index for a SNP.

        Parameters
        ----------
        snp     SNP identifier
        """
        return self.snp.index[self.snp.snp == snp][0]

    def window(self, chrom, start=0, end=np.inf):
        """Get the index for an inclusive window on a chromosome like this:

        Parameters
        ----------
        chrom : int [1-22, X=23]
            A chromosome.
        start : int
        end : int
            Start and end of window.

        Example
        -------
            >>> ASW = HapMap("ASW")
            >>> ASW.window(1, 1e7, 2e7).shape
            (2730,)
        """
        # Fix if start and end are in the wrong order.
        start, end = sorted([start, end])
        # Indices for SNPs on this chromosome.
        idx = getattr(self, 'chr{}'.format(chrom))
        # Indices for SNPs with positions in the window.
        return idx[(start <= self.snp.pos[idx]) & (self.snp.pos[idx] <= end)]

    def windows(self, chrom, size=1e7):
        """Split a chromosome into equal-sized non-overlapping windows and
        return a grouped ``pd.DataFrame`` object for convenient iteration.

        Parameters
        ----------
        chrom : int [1-22, X=23]
            A chromosome.
        size : int
            Size of the window.
        """
        att = 'chr{}'.format(chrom)
        idx = getattr(self, att)
        win = self.snp.pos[idx] 
        return self.snp.ix[idx].groupby(win // size)

    def compute_freqs(self):
        """Update the allele counts and frequencies in the ``geno`` DataFrame:

            ref         Reference allele count for each SNP.
            var         2 - ref (or variant allele count).
            ref_freq    ref / (ref + var)
            maf         min(ref_freq, 1 - ref_freq)
        """
        # Reference and alternate allele counts and frequencies.
        self.ref = self.geno.sum(axis=1)
        self.var = (2 - self.geno).sum(axis=1)
        self.ref_freq = self.ref / (self.ref + self.var)
        # Minor allele frequencies.
        self.maf = self.ref_freq.apply(lambda x: min(x, 1 - x))
        return self

    def _geno_dataframe(self, geno_file, check=False):
        """Read the POP.geno file that contains a matrix of genotypes.
        
        Parameters
        ----------
        geno_file   Full path to the POP.geno file.
        check       Ensure the file contains only these values: 0, 1, 2, NaN
        """
        # Translate character codes '2', '1', '0' to integers 2, 1, 0.
        d = np.frombuffer(open(geno_file).read(), dtype=np.int8) - 48
        # Remove newlines.
        d = d[d != -38]
        # Number of individuals.
        n = self.ind.shape[0]
        d = d.reshape((len(d) / n, n))
        # Column and row names are taken from the other files.
        d = pd.DataFrame(d, columns=self.ind.ind, index=self.snp.snp)
        # 0, 1, 2 are valid values. 9 is null.
        d = d.replace(9, np.nan)
        # Sanity check. Ensure the geno file has no invalid values.
        if check:
            invalid_values = (d != 0) & (d != 1) & (d != 2) & ~np.isnan(d)
            assert 0 == invalid_values.sum().sum(), "Bad values in geno file!"
        return d

    def _ind_dataframe(self, x):
        cols = ["ind", "sex", "pop"]
        return pd.read_table(x, sep=" +", header=None, names=cols)

    def _snp_dataframe(self, x):
        cols = ["snp", "chrom", "cM", "pos", "ref", "alt"]
        return pd.read_table(x, sep=" +", header=None, names=cols)

    def pairwise_r2(self, chrom):
        """Compute r^2 (Pearson correlation coefficient) for all pairs of
        adjacent SNPs along a chromosome.
        """
        g = self.geno
        pairs = pairwise(self.snp.index[self.snp.chrom == chrom])
        return pd.Series([r2(g.values[i], g.values[j]) for i, j in pairs])

    def copy(self):
        """Return a deep copy of this object and all contained data."""
        return deepcopy(self)

    def __repr__(self):
        f = '{} pop {} SNPs {} individuals {}/{} M/F\n\n' \
            '{}\n{}\n\n' \
            '{}\n{}\n\n' \
            '{}\n{}\n\n'
        return f.format(self.pop,
                        self.geno.shape[0], self.geno.shape[1],
                        self.ind.sex.value_counts()['M'],
                        self.ind.sex.value_counts()['F'],
                        self.ind_file,
                        self.ind.head(5).to_string(),
                        self.geno_file,
                        self.geno.iloc[:5, :5].to_string(),
                        self.snp_file,
                        self.snp.head(5).to_string())


def admixture_chisq(n, g, t):
    """Case-only admixture association chisq (1 dof). Week 4 Slide 133.

        2 [  2 N g log(g / t)  +  2 N (1-g) log( (1-g) / (1-t) )  ]

    Parameters
    ----------
    n : int
        Number of cases.
    g : float
        Local ancestry.
    t : float
        Genome-wide ancestry.

    Example
    -------
    >>> admixture_chisq(100, 0.12, 0.08)
    3.8153047642583502

    """
    a = 2 * n * g * np.log(g / t)
    b = 2 * n * (1 - g) * np.log((1 - g) / (1 - t))
    return 2 * (a + b)


def local_ancestry(pop, pop1, pop2, chroms=[], windowsize=1e7,
                   precision=0.001, alphas=[0, 0.5, 1.0]):
    """Estimate local ancestry along a chromosome in non-overlapping windows
    for each individual in ``pop`` as a linear combination of the individuals
    from two other populations ``pop1`` and ``pop2``.

    Parameters
    ----------
    pop : HapMap
        Calculate local ancestry estimates for individuals in this population.
    pop1 : HapMap
    pop2 : HapMap
        Estimate the most likely linear combination of these populations.
    chroms : list of ints
        List of chromosomes [1-22, X=23] or all by default.
    windowsize : int
        Number of base pairs in each non-overlapping window.
    precision : float
        Estimate alpha (linear combination parameter) to this precision.
    alphas : list of floats
        If provided, only test these values of alpha

    Example
    -------
    >>> local_ancestry(ASW, CEU, YRI, 22, alphas=[])
            Ind  Chrom  Window  SNPs     CEU     YRI
    0   NA19916     22       2  1182  0.0760  0.9240
    1   NA19916     22       3  2495  0.0760  0.9240
    2   NA19916     22       4  2575  0.0420  0.9580
    3   NA19916     22       5  3195  0.0005  0.9995
    """
    # Compute all chromosomes by default.
    if type(chroms) != list:
        chroms = [chroms]
    if chroms == []:
        chroms = range(1, 24)
    result = []
    for ind in pop.ind.ind:
        for chrom in chroms:
            for j, window in pop.windows(chrom, windowsize):
                # The individual's genotypes in this window.
                g = pop.geno.ix[window.index, ind].values
                # Reference allele frequencies for this chromosome.
                p1 = pop1.ref_freq[window.index].values
                p2 = pop2.ref_freq[window.index].values
                # Find the maximum likelihood estimate of alpha.
                alpha = mle_alpha(g, p1, p2, precision, alphas)
                result.append([ind, chrom, j, len(g), alpha, 1 - alpha])
    cols = ['Ind', 'Chrom', 'Window', 'SNPs', pop1.pop, pop2.pop]
    return pd.DataFrame(result, columns=cols)


def compare_populations(pops):
    """Compute Fst between all combinations of given populations.

    Parameters
    ----------
    pops : list of HapMap objects

    Example
    -------
    >>> CEU = HapMap("CEU", "/path/to/data")
    >>> CHB = HapMap("CHB", "/path/to/data")
    >>> YRI = HapMap("YRI", "/path/to/data")
    >>> print compare_populations([CEU, CHB, YRI])
    Pops     Chrom   Fst1   Fst2
    CEU,CHB  1       0.103  0.112
    CEU,YRI  1       0.139  0.151
    CHB,YRI  1       0.156  0.171
    ...
    CEU,CHB  22      0.100  0.107
    CEU,YRI  22      0.141  0.152
    CHB,YRI  22      0.164  0.179
    CEU,CHB  Genome  0.102  0.109
    CEU,YRI  Genome  0.141  0.153
    CHB,YRI  Genome  0.159  0.173
    """
    out = []
    for c in range(1, 24): # X=23
        # Fst for all pairs of populations between chromosomes.
        for a, b in combinations(pops, 2):
            fa, fb = [p.ref_freq[p.chrom(c)] for p in (a, b)]
            na, nb = [p.ind.shape[0] for p in (a, b)]
            f1, f2 = [f(fa, fb, na, nb) for f in (Fst1, Fst2)]
            out.append(["{},{}".format(a.pop, b.pop), c, f1, f2])
    for a, b in combinations(pops, 2):
        na, nb = [p.ind.shape[0] for p in (a, b)]
        # Fst for all pairs of populations between whole genomes.
        f1, f2 = [f(a.ref_freq, b.ref_freq, na, nb) for f in (Fst1, Fst2)]
        out.append(["{},{}".format(a.pop, b.pop), "Genome", f1, f2])
    return pd.DataFrame(out, columns=["Pops", "Chrom", "Fst1", "Fst2"])


def compare_chromosomes(pops, chroms=2):
    """Compute F_{ST} between all combinations of given populations using
    a single chromosome.

    Parameters
    ----------
    pops : list of HapMap objects
    chroms : int or list [1-22, X=23]
    """
    if type(chroms) != list:
        chroms = [chroms]
    out = []
    for chrom in chroms:
        for a, b in combinations(pops, 2):
            na, nb = [p.ind.shape[0] for p in (a, b)]
            fa, fb = [p.ref_freq[p.chrom(chrom)] for p in (a, b)]
            f1, f2 = [f(fa, fb, na, nb) for f in [Fst1, Fst2]]
            out.append(["{},{}".format(a.pop, b.pop), chrom, f1, f2])
    return pd.DataFrame(out, columns=["Pops", "Chrom", "Fst1", "Fst2"])


def match_pops(pops1, pops2, chrom):
    """For each /individual/ in the second set of populations, find the best
    matching /population/ from the first set of populations. The best match
    has allele frequencies with the greatest likelihood to generate the
    individual's genotypes.

    Parameters
    ----------
    pops1 : list of HapMap objects
    pops2 : list of HapMap objects
    chrom : int within [1-22, X=23]
        Use genotypes from this chromosome.
    """
    # Select SNP ids on a single chromosome.
    for p in flatten([pops1] + [pops2]):
        p.chrom = p.snp.index[p.snp.chrom == chrom]
    result = {}
    for p1 in pops1:
        ps = p1.ref_freq[p1.chrom]
        result[p1.pop] = []
        for p2 in pops2:
            gs = p2.geno.ix[p2.chrom]
            for i in p2.ind.ind:
                result[p1.pop].append(ll(gs[i].values, ps.values))
    d = pd.DataFrame(result)
    d.index = flatten([p.ind.ind.tolist() for p in pops2])
    d['BestMatch'] = d.apply(lambda r: r.idxmax(), axis=1)
    d['Population'] = flatten([[p.pop] * p.geno.shape[1] for p in pops2])
    return d


def cluster(pops, n, chroms=[]):
    """Cluster individuals from the given populations into n clusters using
    the genotypes from one chromosome.

    Parameters
    ----------
    pops : list of HapMap objects
    n : int
        Cluster the individuals into this many groups.
    chroms : int or list of ints [1-22, X=23]
        Use genotypes from these chromosomes.

    Example
    -------
    >>> cluster([TSI, JPT, LWK), 3, 22)
    Step 0
    TSI 1002120211 2020222120 2000201210 2002120120 0212012122
        1211221100 2010111210 1201022210 21001210
    JPT 1100101100 0020020202 1101121000 2011212000 2221021000
        2220021110 1222101221 2022121002 101110
    LWK 1212221201 1012110000 2201000112 2120022011 1212211221
        2121111010 0102222222 0101120101 2022020100
    Step 1
    TSI 1002120211 2020222120 2000201210 2002220120 0212212122
        1211221102 2010112210 1202022212 21001211
    JPT 0000000000 0000000000 0100001000 0001000000 0000000000
        0000000000 0000000000 0000000000 000000
    LWK 1212221211 1112111121 2211221112 2121121111 1212211221
        2121111211 2111222222 1111122111 2122121122
    Step 2
    TSI 2222222222 2222222222 2222222222 2222222222 2222222222
        2222222222 2222222222 2222222222 22222222
    JPT 0000000000 0000000000 0000000000 0000000000 0000000000
        0000000000 0000000000 0000000000 000000
    LWK 1111111111 1111111111 1111111111 1111111111 1111111111
        1111111111 1111111111 1111111111 1111111111
    """
    # Get names of all individuals.
    inds = np.array(flatten([p.ind.ind.tolist() for p in pops]))
    # Select genotypes on multiple chromosomes.
    if type(chroms) == list and len(chroms):
        idx = np.concatenate([pops[0].chrom(c) for c in chroms])
    # Select genotypes on a single chromosome.
    else:
        idx = pops[0].chrom(chroms)
    genotypes = (p.geno.ix[idx] for p in pops)
    genotypes = pd.concat(genotypes, axis=1, join='inner')
    # Cluster the individuals with an expectation maximization algorithm.
    return em(genotypes.values, n, freq)


def mle_alpha(g, p1, p2, precision=0.001, alphas=[]):
    """Estimate an individual's fractional ancestry in two populations by
    computing the likelihood of generating the individual's genotypes from a
    linear combination of the two ancestral populations. Ignore NaNs.

    Parameters
    ----------
    g : np.array
        Genotypes for a single individual.
    p1 : np.array
    p2 : np.array
        Reference allele frequencies for each population.
    precision : float
        Estimate alpha (the linear combination parameter) to this precision.
    alphas : list of floats
        Restrict the possible values of alpha to this list.
    Example
    -------
    >>> mle_alpha(ASW.geno.values[:, 1], CEU, YRI)
      CEU      Ind   YRI
    0.079  NA19916 0.921
    """
    alpha = np.nan
    if len(alphas) > 1:
        scores = map(partial(ll_alpha, g, p1, p2), alphas)
        best = np.argmax(scores)
        return alphas[best]
    # Find the argument that produces the maximum value of the function.
    return maximize(partial(ll_alpha, g, p1, p2), precision=precision)


def ll_alpha(g, p1, p2, a):
    """Compute the log-likelihood to generate an individual's genotypes given
    the allele frequencies of two populations and the fractional ancestry of
    the individual. Ignore NaNs.
    See: Model-based clustering when allele frequencies in ancestral
    population are known. Slide 26 week 3.

        sum[   log[ a p1 + (1-a) p2 ] g  
             + log[ a (1-p1) + (1-a) (1-p2) ] (2-g)   ]

    Parameters
    ----------
    g : np.array
        Genotypes for one individual.
    p1 : np.array
        Reference allele frequencies in the first population.
    p2 : np.array
        Reference allele frequencies in the second population
    a : float within [0, 1]
        The individual's fractional ancestry from the first population.
    """
    g, p1, p2 = drop_nan(g, p1, p2)
    x = np.log(a * p1       + (1 - a) * p2      ) * g
    y = np.log(a * (1 - p1) + (1 - a) * (1 - p2)) * (2 - g)
    return np.nansum(x + y)


def random_geno(p, n=56):
    """Generate genotypes for a single SNP in Hardy-Weinberg equilibrium.
    
    Parameters
    ----------
    p : float
        Allele frequency of the reference allele.
    n : int
        Generate this many genotypes.

     Genotype         Frequency
    -------------------------------
        2                p^2
        1           2 * p * (1 - p)
        0             (1 - p)^2 
    """
    ps = [(1 - p) ** 2, 2 * p * (1 - p), p ** 2]
    return np.random.choice(a=[0, 1, 2], size=n, p=ps)


def case_control_association(pop1, pop2, cases, controls,
                             chrom=14, e1=np.nan, e2=np.nan):
    """Perform a case/control association study with cases and controls
    sampled from two populations. Optionally adjust the study for population
    stratification.

    Parameters
    ----------
    pop1 : HapMap
    pop2 : HapMap
        HapMap objects for two populations.
    cases : pair of lists of integers
    controls : pair of lists of integers
        Indices of individuals in pop1 and pop2.
    chrom : [1-22, X=23]
        Only use SNPs from this chromosome.
    e1 : float
    e2 : float
        Adjust the case-control study for population stratification using
        these eigenvalues.
    """
    # Select a single chromosome.
    chr_idx = pop1.snp.index[pop1.snp.chrom == chrom]
    # Concatenate the cases and controls into a single dataframe.
    geno = pd.concat([pop1.geno.ix[chr_idx, cases[0]],
                      pop2.geno.ix[chr_idx, cases[1]],
                      pop1.geno.ix[chr_idx, controls[0]],
                      pop2.geno.ix[chr_idx, controls[1]]],
                     axis=1).astype(float)
    # (number of pop1 cases, number of pop2 cases)
    n_cases = map(len, cases)
    # (number of pop1 controls, number of pop2 controls)
    n_controls = map(len, controls)
    # Case/control status that matches the constructed genotype dataframe.
    case = phen = np.repeat([True, False], [sum(n_cases), sum(n_controls)])
    # Size of each population.
    n_pop1 = len(cases[0]) + len(controls[0])
    n_pop2 = len(cases[1]) + len(controls[1])
    print "\n# Case/control association study for {} {} and {} {}".format(
          n_pop1, pop1.pop, n_pop2, pop2.pop)
    print "{}/{} cases/controls".format(sum(n_cases), sum(n_controls)),
    print " {}/{} {}".format(n_cases[0], n_controls[0], pop1.pop),
    # Ancestry status that matches the constructed genotype dataframe.
    ancestry = np.repeat([True, False, True, False], n_cases + n_controls)
    # What is the F_{ST} between pop1 and pop2?
    pop1_freq = freq(geno.ix[:, ancestry].values)
    pop2_freq = freq(geno.ix[:, ~ancestry].values)
    pp_fst = Fst2(pop1_freq, pop2_freq, n_pop1, n_pop2)
    # Calculate the admixture difference between cases and controls.
    admix = n_cases[0] / sum(n_cases) - n_controls[0] / sum(n_controls)
    print " {} admixture difference".format(sig(admix)),
    # Calculate effective study size as per Alkes' notes on the board.
    effective_size = 4 / (1 / sum(n_cases) + 1 / sum(n_controls))
    print " Neff = {}".format(sig(effective_size))
    # What is the F_{ST} between cases and controls? Week 3, slide 47.
    cc_fst = pp_fst * admix ** 2
    print "Fst({},{})        = {}".format(pop1.pop, pop2.pop, sig(pp_fst))
    print "Fst(cases,controls) = {}".format(sig(cc_fst))
    # Degrees of freedom.
    k = 0
    # Correct for population stratification.
    if ~np.isnan(e1) and ~np.isnan(e2):
        print "** Correcting for population stratification with 1 PC **"
        # An eigenvector with 2 unique values that distinguishes ancestry.
        eig = ancestry * e1 + ~ancestry * e2
        # Normalize so sum of squares is 1.
        eig = eig / np.sqrt(np.dot(eig, eig))
        geno = geno - np.apply_along_axis(proj, 1, geno, eig)
        phen = phen - proj(phen, eig)
        # We adjusted for one vector, so we decrease the degrees of freedom.
        k = 1
    # Compute Armitage Trend Test for each SNP.
    atts = geno.apply(lambda g: att(g.values, phen, k), axis=1)
    # Are the statistics inflated? What is \lambda_{GC}?
    expected_lambda = 1 + effective_size * cc_fst
    lambda_gc = atts.median() / 0.455
    after_gc = atts.median() / lambda_gc / 0.455
    print "lambda  Expected  Observed  After Genomic Control"
    form ="        {:>8}  {:>8}  {:>21}"
    print form.format(*map(sig, [expected_lambda, lambda_gc, after_gc]))


def Fst1(p1, p2, n1, n2):
    """Estimate F_{ST} as the expected value of a ratio. Slide 79 week 1.

                                     1        1  
                 (p1 - p2)^2  -  ( ----  +  ---- ) p (1-p)
                                   2 n1     2 n2
        Mean[    -----------------------------------------    ]
                                2 p (1-p)

    Parameters
    ----------
    p1 : np.array
    p2 : np.array
        Reference allele frequencies for each population.
    n1 : int
    n2 : int
        Number of individuals in each population.
    """
    p = (p1 + p2) / 2
    num = np.power(p1 - p2, 2) - (1 / (2 * n1) + 1 / (2 * n2)) * p * (1 - p)
    den = 2 * p * (1 - p)
    return np.nanmean(num / den)


def Fst2(p1, p2, n1, n2):
    """Estimate F_{ST} as the ratio of sums. Slide 79 week 1.

                                 1        1  
       Sum[  (p1 - p2)^2  -  ( ----  +  ---- ) p (1-p)  ]
                               2 n1     2 n2
       --------------------------------------------------
       Sum[                 2 p (1-p)                   ]

    Parameters
    ----------
    p1 : np.array
    p2 : np.array
        Reference allele frequencies for each population.
    n1 : int
    n2 : int
        Number of individuals in each population.
    """
    p = (p1 + p2) / 2
    num = np.power(p1 - p2, 2) - (1 / (2 * n1) + 1 / (2 * n2)) * p * (1 - p)
    den = 2 * p * (1 - p)
    return np.nansum(num) / np.nansum(den)


def r2(a, b):
    """Pearson correlation coefficient, ignoring NaNs. Slide 125 week 1.
    
        [ Mean(a b) -  Mean(a) Mean(b) ]^2

    Parameters
    ----------
    a : np.array
    b : np.array
    """
    a, b = drop_nan(a, b)
    num = np.power(np.mean(a * b) - np.mean(a) * np.mean(b), 2)
    den = np.var(a) * np.var(b)
    return num / den


def att(a, b, k=0):
    """Armitage Trend Test chisq statistic, ignoring NaNs. Slide 142 week 1.

        (N - K - 1) Corr(a, b)

    Parameters
    ----------
    a : np.array
    b : np.array
    k : int
        Adjust the degrees of freedom.
    """
    a, b = drop_nan(a, b)
    return (len(a) - k - 1) * r2(a, b)


def em(vs, k, fun, iterations=10):
    """Cluster the columns in a matrix into multiple groups with an
    expectation maximization (EM) algorithm.

    Parameters
    ----------
    vs : np.array
        Matrix with columns to cluster into multiple groups.
    k : int
        Cluster the columns into this many groups.
    fun : function with single parameter
        Compute probability vectors of each group with this function.
    iterations : int
        Limit the maximum number of iterations for the EM algorithm.
    """
    n = vs.shape[1]
    # Maintain a summary vector for each group.
    us = np.zeros((k, vs.shape[0]))
    # Initialize group assignments randomly.
    gs = np.random.choice(k, n)
    gs_ = gs.copy()
    results = np.zeros((iterations, n)).astype(int)
    # Expectation maximization algorithm with limited iterations.
    for rep in range(iterations):
        # Update the group summaries.
        for k_ in range(k):
            us[k_] = fun(vs[:, gs == k_])
        # Reassign each individual to the group with maximum likelihood.
        for n_ in range(n):
            v = vs[:, n_]
            lls = [ll(v, us[k_]) for k_ in range(k)]
            gs[n_] = np.argmax(lls)
        # Record the previous group assignments.
        results[rep] = gs_.copy()
        # Check if EM has converged.
        if np.all(gs == gs_):
            break
        # Check if the EM got stuck lumping everyone into one group.
        elif len(np.unique(gs)) == 1:
            gs = np.random.choice(k, n)
        gs_ = gs.copy()
    # Return the history of group assignments for each column.
    # The last row has the final assignment of each column to a group.
    return pd.DataFrame(results[:rep + 1])


def ll(g, p):
    """Log-likelihood to observe counts of multiple events with weighted
    probabilities, ignoring NaNs.

    E.g. Log-likelihood to generate genotypes with given the allele
    frequencies. See: Supervised clustering without fractional ancestry.
    Week 3 slide 25.

        Sum[  log(p) g  +  log(1-p) (2-g)  ]

    Parameters
    ----------
    g : np.array
        Genotypes for one individual.
    p : np.array
        Reference allele frequencies in one population.
    """
    g, p = drop_nan(g, p)
    return np.nansum(np.log(p) * g + np.log(1 - p) * (2 - g))


def maximize(fun, interval=[0,1], precision=0.001):
    """"Maximize a function over the given interval with specified precision.

    Parameters
    ----------
    fun : single-parameter function that returns a float
        Find the argument that maximizes the result of this function.
    interval : list of [start, end]
        Inclusive range to search for the maximum.
    precision : float
        Continue searching until we achieve desired precision.
    """
    # Find alpha within this interval which maximizes the function.
    amin, amax = interval[:2]
    alpha = np.nan
    # Limit iterations based on the desired precision.
    for rep in range(int(np.abs(np.log10(precision)))):
        # Sample 10 points in the interval and find the maximum point.
        step = (amax - amin) / 10.0
        alphas = np.arange(amin, amax + step, step)
        best = np.argmax(map(fun, alphas))
        # Narrow the interval containing the solution.
        amin = alphas[max(0, best - 1)]
        amax = alphas[min(len(alphas) - 1, best + 1)]
        # Stop estimating when we achieve desired precision.
        alpha = (amax + amin) / 2.0
        if amax - amin < precision:
            break
    return alpha


def freq(geno):
    """Compute allele frequencies for a genotype matrix, ignoring NaNs.

        ref / (ref + var)

    Parameters
    ----------
    geno : np.array
        Matrix of SNPs (rows) and individual genotypes (columns).
    """
    ref = np.nansum(geno, axis=1)
    var = np.nansum(2 - geno, axis=1)
    return ref / (ref + var)


def normalize(geno):
    """Normalize a genotype matrix, ignoring NaNs. Subtract each SNP's mean
    (across individuals) and then normalize each SNP by sqrt(p * (1 - p))
    where p is the mean genotype, ignoring NaNs. See slide 121.

    Parameters
    ----------
    geno : np.array
        Matrix of SNPs (rows) and individual genotypes (columns).
    """
    p = np.nanmean(geno, axis=1) / 2
    q = np.sqrt(p * (1 - p))
    return np.apply_along_axis(lambda x: (x - 2 * p) / q, 0, geno)


def cov(x, nan=np.nan):
    """Compute covariance of all pairs of rows in a matrix. Produces different
    results than np.cov(), ignoring NaNs.

    Parameters
    ----------
    x : np.array
    nan : float
        If not np.nan (default), then replace NaNs with this value.
    """
    x[np.isnan(x)] = nan
    s = np.zeros((x.shape[1], x.shape[1]))
    for i in range(x.shape[1]):
        # Subtract the mean from each column.
        x[:, i] = x[:, i] - np.nanmean(x[:, i])
        # Compute the diagonal entries of the covariance matrix.
        s[i, i] = np.nansum(x[:, i] * x[:, i])
    for (i, j) in combinations(range(x.shape[1]), 2):
        # Compute the off-diagonal entries of the covariance matrix.
        s[j, i] = s[i, j] = np.nansum(x[:, i] * x[:, j])
    return s / (2 * x.shape[0])


# Helper functions.


def drop_nan(*xs):
    """With a list of vectors, drop indices where any vector has a NaN value.
    Return a copy of the vectors without NaNs.

    Parameters
    ----------
    *xs : one or more np.array objects   OR   one np.ndarray object
        Return arrays with NaN positions omitted.

    Example
    -------
        >>> a = np.array([np.nan, 1, 2])
        >>> b = np.array([np.nan, np.nan, 3])
        >>> drop_nan(a, b)
        array([[ 2.],
               [ 3.]])
    """
    if len(xs) == 1 and type(xs[0]) == np.ndarray and len(xs[0].shape) > 1:
        xs = xs[0]
    f = lambda x: ~np.isnan(x)
    idx = np.bitwise_and.reduce(np.apply_along_axis(f, 1, xs), 0)
    return np.apply_along_axis(lambda x: x[idx], 1, xs)


def test_drop_nan():
    # Test with two np.arrays.
    a = np.random.random(10)
    b = np.random.random(10)
    a[5] = np.nan
    b[3] = np.nan
    n = ~np.isnan(a) & ~np.isnan(b)
    x, y = drop_nan(a, b)
    assert np.all(x == a[n]) 
    assert np.all(y == b[n])
    # Test with one multidimensional np.array.
    a = np.random.random((5, 5))
    a[0, 1] = np.nan
    a[2, 3] = np.nan
    a = drop_nan(a)
    assert np.all(a.flatten() == a[~np.isnan(a)])
    assert a.shape == (5, 3)


def pairwise(iterable):
    "s -> (s0,s1), (s1,s2), (s2, s3), ..."
    a, b = tee(iterable)
    next(b, None)
    return izip(a, b)


def flatten(*xs):
    "[[a,b,c],[d,e,f]] -> [a,b,c,d,e,f]"
    return list(chain.from_iterable(*xs))


def batches(xs, n=3):
    "[a,b,c,d,e,f,g] -> [[a,b,c],[d,e,f],[g]]"
    xs = iter(xs)
    done = False
    while not done:
        batch = []
        for _ in range(n):
            try:
                batch.append(next(xs))
            except:
                done = True
                break
        if batch:
            yield batch


def subtract_projections(a, b, axis=1):
    """Subtract from each vector a[i] its orthogonal projection onto the line
    spanned by the vector b[j] for all pairs of i and j.
    Week 3, slide 112.

    Parameters
    ----------
    a : np.array
    b : np.array
        Matrices with matching height or with matching width.
    axis : 0 or 1
        Operate on rows by default. Use axis=0 to operate on columns.
    """
    for i in range(b.shape[1 - axis]):
        a = a - np.apply_along_axis(proj, axis, a, b[i])
    return a


def proj(v, u):
    """Project v orthogonally onto the line spanned by u, ignoring NaNs.

    Example
    -------
        >>> a, b = np.random.random(5), np.random.random(5)
        >>> a = a - proj(a, b)
        >>> np.dot(a, b)
        -4.9960036108132044e-16
    """
    return u * np.nansum(u * v) / np.nansum(u * u)


def norm(v):
    """Normalize a vector so mean = 0 and sum of squares = 1, ignoring NaNs.

    Example
    -------
        >>> a = np.random.choice([np.nan] + range(10), 20)
        >>> b = norm(a)
        >>> np.nansum(b * b)
        1.0000000000000002
        >>> np.nanmean(b)
        -1.3829601044592879e-17
    """
    u = v - np.nanmean(v)
    return u / np.sqrt(np.nansum(u * u))


def gram_schmidt(x, axis=1):
    """Perform Gram-Schmidt orthonormalization of a matrix to produce a
    matrix where each row has sum of squares = 1 and mean = 0 and all rows
    are orthogonal to each other so ``np.dot(r[i], r[j]) == 0``.

    Parameters
    ----------
    x : np.array
    axis : 0 or 1
        Normalize all rows by default. Set axis=0 to normalize columns.
    """
    x = np.apply_along_axis(norm, axis, x)
    # Skip the first vector.
    for i in range(1, x.shape[1 - axis]):
        # Subtract from row i the projection of i onto each previous row j.
        for j in range(i):
            x[i] = x[i] - proj(x[i], x[j])
    return np.apply_along_axis(norm, axis, x)


@contextmanager
def log(msg, stream=sys.stderr):
    """Use `with log('message'):` to log and time an event."""
    start = time()
    stream.write("{} # {}\n".format(asctime(), msg))
    yield
    elapsed = int(time() - start + 0.5)
    stream.write("{} # done in {} s\n".format(asctime(), elapsed))
    stream.flush()


def dataframe_string(df, flt=None, **kwargs):
    if flt is None:
        flt = sig
    return df.to_string(float_format=flt, **kwargs)


def sig(x, p=3):
    """Round a number to significant figures, retaining trailing zeros.
    Ported from JS to Python by Randle Taylor. Modified by Kamil Slowikowski.
    Based on the webkit javascript implementation taken from here:
        https://chromium.googlesource.com/external/Webkit/+/Safari-3-branch/JavaScriptCore/kjs/number_object.cpp

    Parameters
    ----------
    x : float
        Format this number.
    p : int
        Display this many significant figures.
    """
    import math
    x = float(x)
    # Zero.
    if x == 0.:
        return "0." + "0" * (p - 1)
    out = []
    # Negative numbers.
    if x < 0:
        out.append("-")
        x *= -1
    # Powers of ten.
    e = int(math.log10(x))
    tens = math.pow(10, e - p + 1)
    n = math.floor(x / tens)
    if n < math.pow(10, p - 1):
        e -= 1
        tens = math.pow(10, e - p + 1)
        n = math.floor(x / tens)

    if abs((n + 1.) * tens - x) <= abs(n * tens - x):
        n += 1

    if n >= math.pow(10,p):
        n /= 10.
        e += 1

    m = "%.*g" % (p, n)

    if e < -6 or e >= p:
        out.append(m[0])
        if p > 1:
            out.append(".")
            out.extend(m[1:p])
        out.append('e')
        if e > 0:
            out.append("+")
        out.append(str(e))
    elif e == (p - 1):
        out.append(m)
    elif e >= 0:
        out.append(m[:e + 1])
        if e + 1 < len(m):
            out.append(".")
            out.extend(m[e + 1:])
    else:
        out.append("0.")
        out.extend(["0"] * -(e + 1))
        out.append(m)

    return "".join(out)

