import time
import pickle
from process_mgf import gen_mzML
import sys
import json
from database import database_generate
import concurrent.futures
from pyteomics import fasta
from process_mgf import aa_combination
import datetime
import math
import pyteomics.mgf as mgf
from match_ptm import match
import csv
from scipy.stats import zscore
import numpy as np


def write_raw_result(_raw_res, file_name):
    with open("{}_raw_result.csv".format(file_name.split('\\')[-1].split('.')[0]), 'w', newline='') as csvfile:
        fieldnames = ['Raw_result']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        for _spectrum in _raw_res:
            writer.writerow({'Raw_result': _spectrum})


def write_final_result(_final_res, file_name):
    if _final_res:
        min_z_score = _final_res[-1][0]
        excellent_z = 2 * abs(min_z_score)
        good_z = 1.5 * abs(min_z_score)
        mild_z = 1.2 * abs(min_z_score)
    with open("{}_final_result.csv".format(file_name.split('\\')[-1].split('.')[0]), 'w', newline='') as csvfile:
        fieldnames = ['Mod_mass', 'Site', 'Score', 'z-score', 'Confidence', 'Constitution']

        writer = csv.writer(csvfile)
        writer.writerow(fieldnames)
        for _spectrum in _final_res:
            total = sum(_spectrum[1][2].values())
            each = [(k1, round(v1/total, 3)) for k1, v1 in _spectrum[1][2].items()]
            confi_grade = 'unsure' if _spectrum[0] < mild_z else 'moderate' if _spectrum[0] < good_z else 'good' if \
                _spectrum[0] < excellent_z else 'significant'
            writer.writerow([_spectrum[1][0], '; '.join(_spectrum[1][2].keys()),
                            _spectrum[1][1], _spectrum[0], confi_grade, str(each)[1:-1]])

        # writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        # writer.writeheader()
        # for _spectrum in _final_res:
        #     writer.writerow({'Mod_mass': _spectrum[1][0], 'Site': '; '.join(_spectrum[1][2]), 'Score': _spectrum[1][1],
        #                      'z-score': _spectrum[0]})


def run(tuple_bag):
    sub_spectra, resolution, xl_m, aa_dict, concatenated_seq, unimod, link_sites, min_mod, max_mod, unknown = tuple_bag
    res = []
    for spectrum in sub_spectra:
        res.append(match(spectrum, resolution, aa_dict, concatenated_seq, link_sites, xl_m, unimod, min_mod, max_mod,
                         max_unknown=unknown))
    return res


if __name__ == "__main__":
    '''...............................................paramenters...................................................'''


    if len(sys.argv) != 1:
        params = sys.argv[1]
    else:
        params = r'params.txt'
    with open(params, 'r') as f:  # need to replace 'N-term', 'C-term' with '[', ']'
        params_json = json.load(f)
    Parse_rule = params_json['Parse_rule']
    Fasta_path = params_json['Fasta_path']
    print('Database is', Fasta_path)
    Data_path = params_json['Data_path']
    print('File is', Data_path)
    Max_length = params_json['Max_length']
    Min_length = params_json['Min_length']
    Miss_cleavage = params_json['Miss_cleavage']
    Num_max_mod = params_json['Num_max_mod'] # exclude fixed modifications
    Num_max_unknown = params_json['Num_max_unknown']  # exclude fixed modifications
    Link_site = params_json['Link_site']  # support multiple sites such as ['K', 'R']

    if not Link_site:
        print('Attention! Activate PTM screening in linear peptides!')

    xl_mass = params_json['xl_mass']  # CBDPS
    workers = params_json['workers'] # multithreading workers
    Fix_mod = params_json['Known_fix_mod']  # {'C': [103.00919 + 57.021464, 'c']}
    Var_mod = params_json['Known_var_mod']  # {'M': [131.04049 + 15.99, 'm']}
    std_aa_mass = params_json['std_aa_mass']
    tolerance_da = params_json['resolution']
    upper_ptm_mass = params_json['upper_ptm_mass']
    lower_ptm_mass = params_json['lower_ptm_mass']
    '''update std_aa_mass by Fix_mod and Var_mod'''
    for key, value in Fix_mod.items():
        std_aa_mass.pop(key)
        std_aa_mass[value[1]] = round(value[0], 5)
    for key, value in Var_mod.items():
        std_aa_mass[value[1]] = round(value[0], 5)

    '''....................................generate mzml for non-cleavable search...................................'''
    # gen_mzML(Data_path, Data_path.replace('mgf', 'mzML'))

    '''......................................generate database and save to pickle file..............................'''

    starttime = datetime.datetime.now()
    Fasta_list = list(fasta.read(Fasta_path))  # without decoy database

    '''create the new fix_mod and var_mod mapping'''
    Fix_mod = {key: Fix_mod[key][1] for key in Fix_mod.keys()}
    Var_mod = {key: Var_mod[key][1] for key in Var_mod.keys()}

    interval = math.ceil(len(Fasta_list) / workers)
    args = []
    for i in range(workers):
        args.append((Fasta_list[i * interval: (i + 1) * interval], Parse_rule, Num_max_mod, Max_length, Min_length,
                     Miss_cleavage, Link_site, Var_mod, Fix_mod))
    print('Database in silico digestion...')
    with concurrent.futures.ProcessPoolExecutor(max_workers=workers) as executor:
        result = executor.map(database_generate, args)
    db_peptide = set()
    for i in result:
        db_peptide.update(i)

    '''Concatenate peptide connceted with $ sign'''
    db_peptide = '$'.join(db_peptide)
    print('Peptides number:', db_peptide.count('$'), '; Time spent:', datetime.datetime.now() - starttime)

    # with open('db_peptide', 'wb') as pickout:
    #     pickle.dump(db_peptide, pickout)
    # print(db_peptide)

    # with open('db_peptide', 'rb') as pickout:
    #     db_peptide = pickle.load(pickout)
    # print(db_peptide)
    '''...........................................allocate spectra and multiprocess...............................'''
    t1 = datetime.datetime.now()
    std_aa_mass = aa_combination(std_aa_mass)  # combination

    with open('db_unimod', 'r') as f:  # need to replace 'N-term', 'C-term' with '[', ']'
        db_unimod = json.load(f)

    mass_unimod = dict()

    for key, value in db_unimod.items():
        masses = set([i[0] for i in value if upper_ptm_mass > i[0] > lower_ptm_mass])
        if key == 'N-term':
            name = '['
        elif key == 'C-term':
            name = ']'
        else:
            name = key
        if masses:
            mass_unimod[name] = sorted(masses)
    spectra = list(mgf.read(Data_path))
    # interval = math.ceil(len(spectra) / workers)
    spectra_div = 100
    section = [i for i in range(len(spectra)) if i % spectra_div == 0]
    section_length = len(section)
    args = []
    for i in range(section_length):
        if i != len(section) - 1:
            args.append((spectra[section[i]: section[i + 1]], tolerance_da, xl_mass, std_aa_mass, db_peptide,
                         mass_unimod, Link_site, lower_ptm_mass, upper_ptm_mass, Num_max_unknown))
        else:
            args.append((spectra[section[i]:], tolerance_da, xl_mass, std_aa_mass, db_peptide,
                         mass_unimod, Link_site, lower_ptm_mass, upper_ptm_mass, Num_max_unknown))

    print('ptm screening...')
    with concurrent.futures.ProcessPoolExecutor(max_workers=workers) as executor:
        result = executor.map(run, args)
    print('start writing result...')
    raw_res = []
    final_res = []
    '''write_raw_result(result, Data_path)'''
    for i in result:
        for j in i:
            ele = j.split('_')[1]
            if ele:
                raw_res.append(j)
                score = 0
                ptm = []
                for sub in ele.split(';'):
                    info = sub.split('$')
                    score += float(info[1])
                    for mod in info[3:]:
                        mod_site = int(mod.split(':')[0])
                        mod_mass = float(mod.split(':')[1])
                        mod_aa = info[0][mod_site]
                        if mod_site == 0:
                            types = 'n-term'
                        elif mod_site == len(info[0]) - 1:
                            types = 'c-term'
                        else:
                            types = 'regular'
                        ptm.append((mod_mass, mod_aa, types))
                for p in ptm:
                    final_res.append((p[0], p[1], p[2], score))  # (mod mass, mod aa, type, score)

    final_set = dict()
    final_score = dict()
    for i in final_res:
        # print(i)
        final_set.setdefault(i[0], list())
        final_score.setdefault(i[0], 0)
        final_score[i[0]] += i[3]
        final_set[i[0]].append((i[1], i[2], i[3]))
    temp = []
    majority = 0.7  # if specific aa occupies major nterm or cterm, it is then aa specific mod
    for k, v in final_set.items():
        nterm = dict()
        cterm = dict()
        regular = dict()
        for vv in v:
            if vv[1] == 'c-term':
                cterm.setdefault(vv[0], 0)
                cterm[vv[0]] += vv[2]
            elif vv[1] == 'n-term':
                nterm.setdefault(vv[0], 0)
                nterm[vv[0]] += vv[2]
            else:
                regular.setdefault(vv[0], 0)
                regular[vv[0]] += vv[2]

        cterm_total = sum([i for i in cterm.values()])
        nterm_total = sum([i for i in nterm.values()])
        for nterm_k, nterm_v in nterm.items():
            if nterm_v > majority * nterm_total:
                nterm = {nterm_k: nterm_v}
                break

        for cterm_k, cterm_v in cterm.items():
            if cterm_v > majority * cterm_total:
                cterm = {cterm_k: cterm_v}
                break

        if len(nterm) > 1:
            regular['n-term'] = sum(nterm.values())
        elif len(nterm) == 1:
            regular.setdefault(list(nterm)[0], 0)
            regular[list(nterm)[0]] += nterm[list(nterm)[0]]

        if len(cterm) > 1:
            regular['c-term'] = sum(cterm.values())
        elif len(cterm) == 1:
            regular.setdefault(list(cterm)[0], 0)
            regular[list(cterm)[0]] += cterm[list(cterm)[0]]
        temp.append((k, final_score[k], regular))

    confidence = zscore(np.log([i[1] for i in temp]))
    output = sorted(list(zip(confidence, temp)), reverse=True)
    write_raw_result(raw_res, Data_path)
    write_final_result(output, Data_path)
    # print(output)
    print('Screening finished, time spent:', datetime.datetime.now() - t1)

