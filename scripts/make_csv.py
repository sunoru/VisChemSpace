import csv
import numpy as np
import sys

from chemspace import load_nmr, load_ir


def main_nmr(input, output):
    raw = load_nmr(input)
    results = []
    for each in raw:
        p = [x for x in enumerate(each.nmr_vector) if x[1] > 0]
        tmp = [each.smiles, each.bin_width, len(p), len(each.nmr_vector)]
        tmp += p
        results.append(tmp)
    with open(output, 'w', newline='') as fo:
        writer = csv.writer(fo)
        writer.writerows(results)


def main_ir(input, output):
    raw = load_ir(input)
    results = []
    for each in raw.data:
        tmp = [each, raw.wn_range, raw.bin_size, len(raw.data[each][0]), raw.num_bins]
        tmp += list(np.array(raw.data[each]).transpose())
        results.append(tmp)
    with open(output, 'w', newline='') as fo:
        writer = csv.writer(fo)
        writer.writerows(results)


if __name__ == '__main__':
    {
        'nmr': main_nmr,
        'ir': main_ir
    }[sys.argv[1]](sys.argv[2], sys.argv[3])
