import datetime
import time
import pickle
import concurrent.futures
import numpy as np
import os
from pyteomics import mgf, fasta
import shutil
import json
import logging
from collections import deque
import itertools as it
import re

std_aa_mass = {
    'G': 57.02146,
    'A': 71.03711,
    'S': 87.03203,
    'P': 97.05276,
    'V': 99.06841,
    'T': 101.04768,
    'C': 103.00919,
    'L': 113.08406,
    'I': 113.08406,
    'N': 114.04293,
    'D': 115.02694,
    'Q': 128.05858,
    'K': 128.09496,
    'E': 129.04259,
    'M': 131.04049,
    'H': 137.05891,
    'F': 147.06841,
    'U': 150.95364,
    'R': 156.10111,
    'Y': 163.06333,
    'W': 186.07931,
    'O': 237.14773,
}
expasy_rules = {
    'arg-c': r'R',
    'asp-n': r'\w(?=D)',
    'bnps-skatole': r'W',
    'caspase 1': r'(?<=[FWYL]\w[HAT])D(?=[^PEDQKR])',
    'caspase 2': r'(?<=DVA)D(?=[^PEDQKR])',
    'caspase 3': r'(?<=DMQ)D(?=[^PEDQKR])',
    'caspase 4': r'(?<=LEV)D(?=[^PEDQKR])',
    'caspase 5': r'(?<=[LW]EH)D',
    'caspase 6': r'(?<=VE[HI])D(?=[^PEDQKR])',
    'caspase 7': r'(?<=DEV)D(?=[^PEDQKR])',
    'caspase 8': r'(?<=[IL]ET)D(?=[^PEDQKR])',
    'caspase 9': r'(?<=LEH)D',
    'caspase 10': r'(?<=IEA)D',
    'chymotrypsin high specificity': r'([FY](?=[^P]))|(W(?=[^MP]))',
    'chymotrypsin low specificity':
        r'([FLY](?=[^P]))|(W(?=[^MP]))|(M(?=[^PY]))|(H(?=[^DMPW]))',
    'clostripain': r'R',
    'cnbr': r'M',
    'enterokinase': r'(?<=[DE]{3})K',
    'factor xa': r'(?<=[AFGILTVM][DE]G)R',
    'formic acid': r'D',
    'glutamyl endopeptidase': r'E',
    'granzyme b': r'(?<=IEP)D',
    'hydroxylamine': r'N(?=G)',
    'iodosobenzoic acid': r'W',
    'lysc': r'K',
    'ntcb': r'\w(?=C)',
    'pepsin ph1.3': r'((?<=[^HKR][^P])[^R](?=[FL][^P]))|'
                    r'((?<=[^HKR][^P])[FL](?=\w[^P]))',
    'pepsin ph2.0': r'((?<=[^HKR][^P])[^R](?=[FLWY][^P]))|'
                    r'((?<=[^HKR][^P])[FLWY](?=\w[^P]))',
    'proline endopeptidase': r'(?<=[HKR])P(?=[^P])',
    'proteinase k': r'[AEFILTVWY]',
    'staphylococcal peptidase i': r'(?<=[^E])E',
    'thermolysin': r'[^DE](?=[AFILMV])',
    'thrombin': r'((?<=G)R(?=G))|'
                r'((?<=[AFGILTVM][AFGILTVWA]P)R(?=[^DE][^DE]))',
    'trypsin': r'([KR](?=[^P]))|((?<=W)K(?=P))|((?<=M)R(?=P))',
    'trypsin_exception': r'((?<=[CD])K(?=D))|((?<=C)K(?=[HY]))|((?<=C)R(?=K))|((?<=R)R(?=[HR]))',
}
_nterm_mod = r'[^-]+-$'
_cterm_mod = r'-[^-]+$'
std_nterm = 'H-'
std_cterm = '-OH'


def mass_calculation(sequence, aa_mass):
    """calculate the mono-isotopic mass given the specific peptide"""
    gen = [a for a in range(len(sequence)) if str.isupper(sequence[a])]
    if sequence[0:gen[0] + 1] not in aa_mass:  # mass of water and the first amino acid
        _mass = 18.01055 + aa_mass[sequence[0:gen[0]]] + aa_mass[sequence[gen[0]]]
    else:
        _mass = 18.01055 + aa_mass[sequence[gen[0]]]
    for _index in range(1, len(gen)):
        if sequence[gen[_index - 1] + 1:gen[_index] + 1] in aa_mass:
            _mass += aa_mass[sequence[gen[_index]]]
        else:
            _mass += aa_mass[sequence[gen[_index - 1] + 1:gen[_index]]] + aa_mass[sequence[gen[_index]]]
    return _mass


def cleave(sequence, rule, missed_cleavages=0, min_length=1, Methionine_drop=True):
    """cleave the protein into peptides and specify their relative positions"""
    rule = expasy_rules.get(rule, rule)
    peptides = []
    ml = missed_cleavages + 2
    trange = range(ml)
    cleavage_sites = deque([0], maxlen=ml)
    cl = 1

    for ii in it.chain([x.end() for x in re.finditer(rule, sequence)], [None]):
        cleavage_sites.append(ii)
        if cl < ml:
            cl += 1
        for jj in trange[:cl - 1]:
            seq = sequence[cleavage_sites[jj]:cleavage_sites[-1]]
            if seq and len(seq) >= min_length:
                peptides.append((seq, cleavage_sites[jj]))  # peptide sequence and position (index from 0)
    if Methionine_drop:
        '''add methionine dropped peptides'''
        meth_group = [ele[1:] for ele, index in peptides if (ele[0] == 'M' and index == 0 and len(ele) > min_length)]
        for _ele in meth_group:
            peptides.append((_ele, 1))
    return peptides


def database_generate(tuple_bag):
    """primary amine, peptide n/c and protein n/c needed"""
    fas_list, parse_rule, num_max_mod, max_length, min_length, miss_cleavage,\
        link_site, var_mod, fix_mod = tuple_bag
    db_peptides = set()
    for description, sequence in fas_list:

        sequence = sequence.replace('L', 'I')  # replace L to I

        if all(word not in ['B', 'J', 'X', 'Z', 'O', 'U'] for word in sequence):
            # form_description = description.split(' ')[0]
            new_peptides = cleave(sequence, parse_rule, missed_cleavages=miss_cleavage, min_length=min_length)

            '''special attention for xl-ms/regular peptide digestion'''
            if link_site:  # xl-ms peptides digestion
                new_peptides = (new_peptide for new_peptide in new_peptides if len(new_peptide[0]) <= max_length and
                                not set(new_peptide[0][0:-1]).isdisjoint(link_site))
            else:  # else link_site is empty, activate linear peptide search
                new_peptides = (new_peptide for new_peptide in new_peptides if len(new_peptide[0]) <= max_length)

            for element, index in new_peptides:
                base_peptide = str(element)
                for fix in fix_mod.keys():
                    base_peptide = base_peptide.replace(fix, fix_mod[fix])
                forms = [base_peptide]
                for i in range(len(base_peptide)):
                    if base_peptide[i] in var_mod.keys():
                        for j in range(len(forms)):
                            forms.append(forms[j][:i] + var_mod[base_peptide[i]] + forms[j][i+1:])
                forms = [i for i in forms if sum([i.count(j) for j in var_mod.values()]) <= num_max_mod]

                db_peptides.update(forms)

    return db_peptides


if __name__ == '__main__':
    '''import configuration'''

    Parse_rule = 'trypsin'
    Fasta_path = 'BSA.fasta'
    Max_length = 50
    Min_length = 5
    Miss_cleavage = 2
    Num_max_mod = 2  # exclude fixed modifications
    Link_site = ['K']  # support multiple sites such as ['K', 'R']
    Fix_mod = {'C': [103.00919 + 57.021464, 'c']}
    Var_mod = {'M': [131.04049 + 15.99, 'm']}

    starttime = datetime.datetime.now()
    '''add modification masses into the dict'''

    Fasta_list = list(fasta.read(Fasta_path))  # without decoy database


    '''create the new fix_mod and var_mod mapping'''
    Fix_mod = {key: Fix_mod[key][1] for key in Fix_mod.keys()}
    Var_mod = {key: Var_mod[key][1] for key in Var_mod.keys()}

    div_data = 1  # divide fasta by every # number
    section = [i for i in range(len(Fasta_list)) if i % div_data == 0]
    section_length = len(section)
    args = []
    for i in range(section_length):
        args.append((Fasta_list[section[i]:section[i] + div_data], Parse_rule, Num_max_mod, Max_length, Min_length,
                     Miss_cleavage, Link_site, Var_mod, Fix_mod))
    print('start multiprocessing')
    with concurrent.futures.ProcessPoolExecutor() as executor:
        result = executor.map(database_generate, args)
    for i in result:
        print(i)