import argparse
import os

import numpy as np
import pandas as pd

from chemspace import Fingerprints, load_data, NMRVector, IRVector, ChemicalSpaceGraph


def parse_args():
    parser = argparse.ArgumentParser(description='The complete pipeline.')
    parser.add_argument('nmr_input', metavar='NMR', type=str, default='',
                        help='NMR input file (pickle format).')
    parser.add_argument('ir_input', metavar='IR', type=str, default = '',
                        help='IR (loader) input file (pickle format).')
    parser.add_argument('--nmr_width', type=float, default=0.1,
                        help='Width of each bin in the fingerprints for NMR spectra.')
    parser.add_argument('--ir_width', type=float, default=0.2,
                        help='Width of each bin in the fingerprints for IR spectra.')
    parser.add_argument('--threshold', type=float, default=3000,
                        help='Threshold for similarity.')
    parser.add_argument('--fingerprints-cache', type=str, default='cache.h5',
                        help='Use the cache file for fingerprints.')
    parser.add_argument('--graph-cache', type=str, default='cache.p',
                        help='Use the cache file for graph.')
    parser.add_argument('--gephi-workspace', type=str, default='chemspace',
                        help='The workspace name in Gephi.')
    return parser.parse_args()


def main(args):
    if os.path.exists(args.fingerprints_cache) or os.path.exists(args.graph_cache):
        graph = ChemicalSpaceGraph.from_file(args.fingerprints_cache, args.graph_cache)
    else:
        nmr, ir = load_data(args.nmr_input, args.ir_input)
        nmr_vectors = [NMRVector.from_old(x) for x in nmr]
        ir_vectors = [IRVector.from_old(ir, x.smiles) for x in nmr_vectors if x.smiles in ir.data]
        nmr_vectors = [x for x in nmr_vectors if x.smiles in ir.data]
        # because maybe they are not consistent.
        assert all(nmr.smiles == ir.smiles for nmr, ir in zip(nmr_vectors, ir_vectors))
        length = len(nmr_vectors)
        print('Number of molecules: %d' % length)
        all_fingerprints = []
        for i, (nmr, ir) in enumerate(zip(nmr_vectors, ir_vectors)):
            print('Generating fingerprints for %s... (%d/%d)\r' % (nmr.smiles, i + 1, length), end='')
            all_fingerprints.append(Fingerprints(nmr, ir, nmr_width=args.nmr_width, ir_width=args.ir_width))
        print()
        graph = ChemicalSpaceGraph()
        graph.set_fingerprints(all_fingerprints)
        print('Saving fingerprints to %s...' % args.fingerprints_cache)
        graph.save_fingerprints(args.fingerprints_cache)
    if not graph.edges:
        graph.generate_graph(threshold=args.threshold, log=True)
        print('Saving graph to %s...' % args.graph_cache)
        graph.save_graph(args.graph_cache)
    distances = pd.Series(list(graph.edges.values()))
    print(distances.describe())
    print('Starting Gephi streaming...')
    graph.gephi_start(args.threshold, args.gephi_workspace)


if __name__ == '__main__':
    main(parse_args())
