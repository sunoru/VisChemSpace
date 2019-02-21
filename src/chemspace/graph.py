import os

import numpy as np
import pandas as pd
import pickle
from gephistreamer import graph, streamer
from itertools import combinations
# from sklearn.metrics import cohen_kappa_score
from chemspace import Fingerprints


class ChemicalSpaceGraph:
    fingerprints_df: pd.DataFrame = None
    nodes: [str] = None
    edges: {(str, str): float} = None
    stream = None
    base: str

    def __init__(self):
        pass

    @staticmethod
    def similarity(x: pd.Series, y: pd.Series):
        # TODO: Find a good method.
        # It's now just Euclidean distance.
        return np.linalg.norm(x - y)

    def set_fingerprints(self, fingerprints: [Fingerprints], base=None):
        self.fingerprints_df = pd.DataFrame({
            x.smiles: x.data for x in fingerprints
        })
        if not base:
            base = fingerprints[0].smiles
        self.base = base
        assert self.base in self.fingerprints_df

    def generate_graph(self, threshold=28, log=False):
        self.nodes = list(self.fingerprints_df.columns)
        self.edges = {}
        possible_edges = list(combinations(self.nodes, 2))
        length = len(possible_edges)
        for i, (v1, v2) in enumerate(possible_edges):
            x = self.fingerprints_df[v1]
            y = self.fingerprints_df[v2]
            if log and i % 5000 == 0:
                print('Calculating the distance between %s and %s... (%d/%d)' % (v1, v2, i + 1, length), end='')
            d = self.similarity(x, y)
            if log and i % 5000 == 0:
                print(' %f\r' % d, end='')
            if d <= threshold:
                self.edges[(v1, v2)] = d
        if log:
            print('\nGraph generated: %d nodes and %d edges.' % (len(self.nodes), len(self.edges)))

    def save_fingerprints(self, filename):
        self.fingerprints_df.to_hdf(filename, 'fingerprints')

    def save_graph(self, filename):
        with open(filename, 'wb') as fo:
            pickle.dump({
                'nodes': self.nodes,
                'edges': self.edges
            }, fo)

    @classmethod
    def from_file(cls, fingerprints, edges, base=None):
        x = cls()
        if fingerprints and os.path.exists(fingerprints):
            x.fingerprints_df = pd.read_hdf(fingerprints, 'fingerprints')
        if edges and os.path.exists(edges):
            with open(edges, 'rb') as fi:
                data = pickle.load(fi)
                x.nodes = data['nodes']
                x.edges = data['edges']
        if not base:
            base = x.nodes[0]
        x.base = base
        return x

    def gephi_start(self, threshold=28, workspace='chemspace'):
        self.stream = streamer.Streamer(streamer.GephiWS(workspace=workspace))

        base = self.fingerprints_df[self.base]
        distances = [[ChemicalSpaceGraph.similarity(self.fingerprints_df[x], base), x] for x in self.nodes]
        distances.sort()
        rankings = {
            x[1]: [i, x[0]]
            for i, x in enumerate(distances)
        }
        nodes = [
            graph.Node(
                x,
                distance=rankings[x][1],
                ranking=rankings[x][0]
            ) for x in self.nodes
        ]
        self.stream.add_node(*nodes)
        similarities = pd.Series(list(self.edges.values()))
        m = similarities.max()
        edges = [
            graph.Edge(x, y, directed=False, weight=1 - self.edges[(x, y)] / m, label=self.edges[(x, y)])
            for x, y in self.edges
            if self.edges[(x, y)] <= threshold
        ]
        self.stream.add_edge(*edges)
