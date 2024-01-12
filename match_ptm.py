import ctypes
import time
from process_mgf import spectrum_process
from graph import find_all_paths, graph_path_filter
import matplotlib.pyplot as plt
import pyteomics.mgf as mgf
import pickle
import json

'''#####################################################'''


def match(spectrum, resolution, aa_dict, db_peptide, link_sites, xl_mass, db_unimod, min_ptm, max_ptm, max_unknown=2):
    # print('working on {}'.format(spectrum['params']['scans']))
    clibrary = ctypes.CDLL(r".\cpplibrary.so")
    # aa_dict should contain some two AAs, db_unimod should be K:[mass, mass]
    PROTON = 1.007276
    pre_mass = (spectrum['params']['pepmass'][0] - PROTON) * spectrum['params']['charge'][0]

    '''process spectrum, reorder, merge nearby peaks'''
    # print('i am here process spectrum')
    mz_py, intensity_py = spectrum_process(list(spectrum['m/z array']), list(spectrum['intensity array']),
                                           pre_mass, atol=resolution)
    mz = (ctypes.c_double * len(mz_py))()
    intensity = (ctypes.c_double * len(intensity_py))()
    resize = int(mz_py[-1] / resolution) + 1
    digit_spectrum = (ctypes.c_double * resize)()  # digitized spectrum
    for i in range(len(mz_py)):
        mz[i] = mz_py[i]
        intensity[i] = intensity_py[i]

    '''digitize spectrum in SEQUEST way'''
    # print('i am here sequest')
    digitizeSpectrum = clibrary.digitizeSpectrum
    digitizeSpectrum.argtypes = [ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double), ctypes.c_int,
                                 ctypes.POINTER(ctypes.c_double), ctypes.c_double]
    digitizeSpectrum(mz, intensity, len(mz), digit_spectrum, resolution)

    '''tag construction'''
    # print('i am here tag constrction')
    length = 10000  # suppose at most #length tag1
    start_idx = (ctypes.c_int * length)()
    end_idx = (ctypes.c_int * length)()
    char_name = (ctypes.c_char_p * length)()

# std_aa_mass = aa_combination(std_aa_mass)  # combination

    aa_length = len(aa_dict)
    aa_mass = (ctypes.c_double * aa_length)()
    aa_name = (ctypes.c_char_p * aa_length)()
    n = 0
    for key, value in aa_dict.items():
        aa_name[n] = key.encode()
        aa_mass[n] = value
        n += 1

    constructTag = clibrary.constructTag
    constructTag.argtypes = [ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double), ctypes.c_int,
                             ctypes.c_double, ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int),
                             ctypes.POINTER(ctypes.c_char_p), ctypes.c_int, ctypes.POINTER(ctypes.c_char_p),
                             ctypes.POINTER(ctypes.c_double), ctypes.c_int]
    constructTag(mz, intensity, len(mz), resolution, start_idx, end_idx, char_name, length, aa_name, aa_mass, aa_length)

    '''construct graph given the above discovered tag1'''
    # n = 1
    '''build directed graph in dict'''
    # print('i am here graph')
    graph = dict()
    start_points = []
    end_points = []
    n = 0
    for i in range(len(start_idx)):
        if end_idx[i] != 0:
            start_points.append(start_idx[i])
            end_points.append(end_idx[i])
            graph.setdefault(start_idx[i], [])
            graph[start_idx[i]].append(tuple([end_idx[i], char_name[i].decode()]))
            # print('new is '.format(n), start_idx[i], end_idx[i], char_name[i])

            n += 1
    # print(start_points)
    # print(end_points)
    start_points = set(start_points) - set(end_points)
    # print(start_points)
    all_paths = []
    # t1 = time.perf_counter()
    for i in start_points:
        paths = find_all_paths(graph, i)
        all_paths = paths + all_paths
    # t2 = time.perf_counter()
    # print('time find all path is ', t2-t1)
    # for i in all_paths:
    #     print('original path is ', i)
    # print('original lenthis', len(all_paths))
    # timer1 = time.perf_counter()
    all_paths = graph_path_filter(all_paths, intensity)
    # timer2 = time.perf_counter()
    # print('total tag markers to search is', len(all_paths), 'time is', timer2 - timer1)
    # print([i for i in mz])
    # print([i for i in intensity])
    # print(len(all_paths))

    '''start fuzzy string match using bitap algorithm'''
    # print('i am here fuzzy')
    bitap = clibrary.bitap

    candidates = dict()  # matched peptides and the corresponding matched position in set(),
    # used to solve the unchange position
    coarse_input = []  # [[mass_1st, b/y, peak_pos, peptide_seq, alias, xl_pos, pattern_len], [], []...]

    text = db_peptide.encode()

    for i in all_paths:
        mass_1st = mz[i[0][0]]
        check_by_ion = 0 if i[1].replace('(', '').replace(')', '') == i[2] else 1
        # print([mz[ii] for ii in i[0]], i[0], i[1], i[2], i[3])
        # print(i)
        pat = i[2]
        size_res = 5000  # at most 5000 peptide candidates per tagN
        pattern = pat.encode()
        res = (ctypes.c_int64 * size_res)(*[-1] * size_res)  # 64 bit in case of very large database, initialize to -1

        bitap.argtypes = [ctypes.c_char_p, ctypes.c_char_p, ctypes.c_int, ctypes.POINTER(ctypes.c_int64),
                          ctypes.c_int, ctypes.c_int]
        bitap(text, pattern, len(pattern), res, size_res, i[3])

        for j in res:
            break_flag = 0  # avoid pattern is SEIRHA, result is SEIR$A
            pos = j
            start = j
            end = j
            if j == -1:
                break
            while start > 0 and db_peptide[start - 1] != '$':
                start -= 1
            while end < len(db_peptide) and db_peptide[end] != '$':
                end += 1
                break_flag += 1
            if break_flag < len(pat):
                continue
            candidates.setdefault(db_peptide[start: end], set())

            alias = db_peptide[start: pos] + pat + db_peptide[pos + len(pat): end]
            if check_by_ion == 0:
                align_pos = pos - start
            else:
                align_pos = pos - start + len(pat)

            coarse_input.append((mass_1st, check_by_ion, align_pos, db_peptide[start: end], alias,
                                 [ii for ii in range(len(alias)-1) if alias[ii] in link_sites], len(pat)))
            for k in range(len(pat)):
                if pat[k] == db_peptide[k + pos]:
                    candidates[db_peptide[start: end]].add(k + pos - start)
            # print('added ', db_peptide[start: end], pos - start)

    '''coarse scoring function to pre-filter the peptides from only tag information, return top N results'''
    # tc1 = time.perf_counter()
    coarse_input_len = len(coarse_input)
    coarse_score_ls = (ctypes.c_double * coarse_input_len)()

    mass_1st_ls = (ctypes.c_double * coarse_input_len)(*[i[0] for i in coarse_input])
    check_by_ion_ls = (ctypes.c_int * coarse_input_len)(*[i[1] for i in coarse_input])
    align_pos_ls = (ctypes.c_int * coarse_input_len)(*[i[2] for i in coarse_input])
    alias_ls = (ctypes.c_char_p * coarse_input_len)(*[i[4].encode() for i in coarse_input])
    sites_ls = (ctypes.c_char_p * coarse_input_len)(*['$'.join([str(ii) for ii in i[5]]).encode() for i in coarse_input])
    pattern_len_ls = (ctypes.c_int * coarse_input_len)(*[i[6] for i in coarse_input])

    coarseScore = clibrary.coarseScore
    coarseScore.argtypes = [ctypes.POINTER(ctypes.c_double), ctypes.c_int, ctypes.POINTER(ctypes.c_double),
                            ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_char_p),
                            ctypes.POINTER(ctypes.c_char_p), ctypes.POINTER(ctypes.c_int), ctypes.c_int,
                            ctypes.c_double, ctypes.POINTER(ctypes.c_char_p), ctypes.POINTER(ctypes.c_double),
                            ctypes.c_int, ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double,
                            ctypes.POINTER(ctypes.c_double)]

    coarseScore(digit_spectrum, resize, mass_1st_ls, check_by_ion_ls, align_pos_ls, alias_ls, sites_ls, pattern_len_ls,
                coarse_input_len, pre_mass, aa_name, aa_mass, aa_length, xl_mass, resolution, min_ptm, max_ptm,
                coarse_score_ls)
    topN = 10
    index = sorted(range(len(coarse_score_ls)), key=lambda ii: coarse_score_ls[ii], reverse=True)[:topN]
    new_candidates = dict()
    for i in index:
        # print(i, coarse_input[i])
        if coarse_score_ls[i] > 0:
            new_candidates[coarse_input[i][3]] = candidates[coarse_input[i][3]]
    candidates = new_candidates
    # tc2 = time.perf_counter()
    # print('time spend on coarse', tc2 - tc1)
    # print('total candidate', coarse_input_len)

    '''discover ptm from the given tag, peptide information'''
    # print('i am here ptm')
    pep_len = len(candidates)  # total number of peptide candidates
    pep_seq = (ctypes.c_char_p * pep_len)()  # peptide sequence char**
    pep_sites = (ctypes.c_char_p * pep_len)()  # cross-linked sites for each peptides
    pep_unchange_pos = (ctypes.c_char_p * pep_len)()  # postions cannot be modified
    n = 0
    # print('length of candidates ', len(candidates))
    # mark1 = False
    # mark2 = False
    for i, j in sorted(candidates.items()):
        pep_seq[n] = i.encode()
        pep_unchange_pos[n] = '$'.join([str(ii) for ii in sorted(j)]).encode()
        pep_sites[n] = '$'.join([str(ii) for ii in range(len(i)-1) if i[ii] in link_sites]).encode()
        # exclude the last aa
        # print(pep_seq[n], pep_unchange_pos[n], pep_sites[n])

        n += 1
    # print('mark',mark1, mark2)
    # with open('db_unimod', 'r') as f:  # need to replace 'N-term', 'C-term' with '[', ']'
    #     db_unimod = json.load(f)

    ptm_len = len(db_unimod)
    ptm_site = (ctypes.c_char_p * ptm_len)()
    ptm_mass = (ctypes.c_char_p * ptm_len)()
    n = 0
    for i, j in db_unimod.items():
        # try -50 to 200 as ptm range
        # masses = set([ele[0] for ele in j if 100 > ele[0] > -50])
        ptm_mass[n] = '$'.join([str(ele) for ele in j]).encode()
        # ptm_mass[n] = '$'.join([str(ele) for ele in sorted(list(masses))]).encode()
        # if i == 'N-term':
        #     name = '['
        # elif i == 'C-term':
        #     name = ']'
        # else:
        #     name = i
        ptm_site[n] = i.encode()
        # print(ptm_site[n], ptm_mass[n])
        n += 1

    # str_res = (ctypes.c_char_p * 2)()  # at most two returned result
    # t1 = time.perf_counter()
    # print('this is time consuming')
    # print('sequence is ', len(pep_seq))
    # tp1 = time.perf_counter()
    matchPTM = clibrary.matchPTM
    matchPTM.argtypes = [ctypes.POINTER(ctypes.c_double), ctypes.c_int, ctypes.POINTER(ctypes.c_char_p),
                         ctypes.POINTER(ctypes.c_char_p), ctypes.POINTER(ctypes.c_char_p), ctypes.c_int,
                         ctypes.c_double, ctypes.POINTER(ctypes.c_char_p), ctypes.POINTER(ctypes.c_double),
                         ctypes.c_int, ctypes.POINTER(ctypes.c_char_p), ctypes.POINTER(ctypes.c_char_p),
                         ctypes.c_int, ctypes.c_double, ctypes.c_double, ctypes.c_int]

    matchPTM.restype = ctypes.c_char_p

    res = matchPTM(digit_spectrum, resize, pep_seq, pep_sites, pep_unchange_pos, pep_len, pre_mass, aa_name, aa_mass, aa_length,
                   ptm_site, ptm_mass, ptm_len, xl_mass, resolution, max_unknown)

    # tp2 = time.perf_counter()
    # print('time spend ptm', tp2 - tp1)
    # print(spectrum['params']['scans'], res.decode())
    # print('finished on {}'.format(spectrum['params']['scans']))
    return spectrum['params']['scans'] + "_" + res.decode()


if __name__ == '__main__':
    spectra = mgf.read('demo.mgf')
    for spectrum in spectra:
        print(spectrum)
