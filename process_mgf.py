from pyteomics import mgf
from psims.mzml.writer import MzMLWriter
import numpy as np
import time


def gen_mzML(path, out_file):
    """
    Convert mgf to mzML file for the non-cleavable XL-MS search engine.
    :param path: input mgf file path, e.g., input.mgf
    :param out_file: output mzML, e.g., output.mzML
    """
    scans = mgf.read(path, convert_arrays=0)
    with MzMLWriter(open(out_file, 'wb'), close=True) as out:
        # Add default controlled vocabularies
        out.controlled_vocabularies()

        out.file_description([  # the list of file contents terms
            "MS1 spectrum",
            "MSn spectrum"
        ])

        out.software_list([
            {"id": "psims-writer",
             "version": "0.1.2",
             "params": [
                 "python-psims",
             ]}
        ])
        source = out.Source(1, ["electrospray ionization", "electrospray inlet"])
        analyzer = out.Analyzer(2, [
            "quadrupole"
        ])
        detector = out.Detector(3, ["inductive detector"])
        config = out.InstrumentConfiguration(id="IC1", component_list=[source, analyzer, detector],
                                             params=["LTQ-FT"])
        out.instrument_configuration_list([config])

        methods = []
        methods.append(
            out.ProcessingMethod(
                order=1, software_reference="psims-writer", params=[
                    "MS:1000035",  # peak picking
                    "Conversion to mzML"
                ]))
        processing = out.DataProcessing(methods, id='DP1')
        out.data_processing_list([processing])

        # Open the run and spectrum list sections
        with out.run(id="my_analysis"):
            spectrum_count = len(scans)
            with out.spectrum_list(count=spectrum_count):
                for scan in scans:
                    zipped = sorted(zip(scan['m/z array'], scan['intensity array']))
                    mz_array, intensity_array = zip(*zipped)
                    bpi = max(intensity_array)
                    bpmz = mz_array[intensity_array.index(bpi)]
                    out.write_spectrum(
                        mz_array, intensity_array,
                        id='controllerType=0 controllerNumber=1 scan={}'.format(scan['params']['scans']),
                        # id=scan['params']['scans'],
                        params=[
                            "MSn Spectrum",
                            {"ms level": 2},
                            {"base peak m/z": bpmz},
                            {"base peak intensity": bpi, 'unitName': "number of detector counts"},
                            {"total ion current": sum(intensity_array)}
                        ],
                        # Include precursor information
                        precursor_information={
                            "mz": scan['params']['pepmass'][0],
                            "intensity": scan['params']['pepmass'][1],
                            "charge": int(scan['params']['charge'][0]),
                            # "scan_id": 99999,
                            "activation": ["N/A", {"collision energy": 0}],
                            "isolation_window": [scan['params']['pepmass'][0] - 1, scan['params']['pepmass'][0],
                                                 scan['params']['pepmass'][0] + 1]
                        })


def spectrum_process(mz, intensity, pre_mass, atol=0.02):
    """generate b y starting peaks, sort peaks and merge peaks"""
    PROTON = 1.007276
    WATER = 18.010564
    # mz_new1 = mz.copy()
    # intensity_new1 = intensity.copy()
    # for idx in range(len(mz)):
    #     intensity_new1[idx] = 1.0 * intensity[idx]
    #     mz_new1.append(round(pre_mass + 2 * PROTON - mz[idx], 5))
    #     intensity_new1.append(round(1.0 * intensity[idx], 5))

    zipped = list(zip(mz, intensity))
    zipped.append((PROTON, 1.0))  # add b0 ion
    zipped.append((PROTON + WATER, 1.0))  # add y0 ion
    zipped.append((pre_mass + PROTON - WATER, 1.0))  # add bN ion
    zipped.append((pre_mass + PROTON, 1.0))  # add yN ion
    zipped = sorted(zipped)
    # print(zipped)
    zipped_new = [zipped[0]]
    for idx in range(1, len(zipped)):
        if zipped[idx][0] < pre_mass + 10.0:  # make sure remove irrelevant peaks
            if abs(zipped_new[-1][0] - zipped[idx][0]) < atol:
                weighted_mass = (zipped_new[-1][0] * zipped_new[-1][1] + zipped[idx][0] *
                                 zipped[idx][1]) / (zipped_new[-1][1] + zipped[idx][1])
                sum_intensity = round(zipped_new[-1][1] + zipped[idx][1], 5)
                zipped_new[-1] = (round(weighted_mass, 5), sum_intensity)
            else:
                zipped_new.append((zipped[idx][0], zipped[idx][1]))

    return zip(*zipped_new)  # sorted


def aa_combination(aa_dict):
    """generate AB combination for any A,B in aa_dict"""
    aa_dict_new = aa_dict.copy()
    max_val = max([i for i in aa_dict.values()])
    for k1, v1 in aa_dict.items():
        for k2, v2 in aa_dict.items():
            if v1 + v2 < max_val + 10.0:
                aa_dict_new['(' + k1 + k2 + ')'] = v1 + v2
    return aa_dict_new


if __name__ == '__main__':
    """Serves to reorder the peaks in the mgf and generate mzML file for searching engine."""
    gen_mzML('ams_00252_mascot.mgf', 'ams_00252_mascot.mzML')


    # a=[1,3,5,7,9]*100
    # b=[1,1,1,1,1]*100
    # t1 = time.perf_counter()
    # c,d = spectrum_process(a, b, 10-2.014552, tol=1e-3)
    # t2 = time.perf_counter()
    # print(c,d)
    # print(t2-t1)


    # m = mgf.read('demo.mgf')
    # m = next(m)
    # # print(m)
    # t1 = time.perf_counter()
    # pre_mass = (m['params']['pepmass'][0] - 1.007276) * m['params']['charge'][0]
    # a, b = spectrum_process(list(m['m/z array']), list(m['intensity array']), pre_mass, 0.02)
    # t2 = time.perf_counter()
    # print(a,b)
    # print(t2 - t1)