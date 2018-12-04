import os

import numpy as np
import pandas as pd
import pickle
from gephistreamer import graph, streamer
from itertools import combinations
from sklearn.metrics import cohen_kappa_score
from chemspace import Fingerprints


class ChemicalSpaceGraph:
    fingerprints_df: pd.DataFrame = None
    nodes: [str] = None
    edges: {(str, str): float} = None
    stream = None

    def __init__(self):
        pass

    @staticmethod
    def similarity(x: pd.Series, y: pd.Series):
        # <del>A method based on Tanimoto Similarity
        # DOI: 10.1021/jm401411z</del>
        # TODO: Find a good method.
        # Cohen's Kappa
        xx = np.ceil(x * 10)
        yy = np.ceil(y * 10)
        return cohen_kappa_score(xx, yy)
        # It's now just Euclidean distance.
        # return 1 / np.linalg.norm(x - y)

    def set_fingerprints(self, fingerprints: [Fingerprints]):
        self.fingerprints_df = pd.DataFrame({
            x.smiles: x.data for x in fingerprints
        })

    def generate_graph(self, threshold=0.2, log=False):
        self.nodes = list(self.fingerprints_df.columns)  # [:200]  # TODO
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
            if d >= threshold:
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
    def from_file(cls, fingerprints, edges):
        x = cls()
        if fingerprints and os.path.exists(fingerprints):
            x.fingerprints_df = pd.read_hdf(fingerprints, 'fingerprints')
        if edges and os.path.exists(edges):
            with open(edges, 'rb') as fi:
                data = pickle.load(fi)
                x.nodes = data['nodes']
                x.edges = data['edges']
        return x

    def gephi_start(self, threshold=0.2, workspace='chemspace'):
        self.stream = streamer.Streamer(streamer.GephiWS(workspace=workspace))
        nodes = [
            graph.Node(x)
            for x in self.nodes
        ]
        self.stream.add_node(*nodes)
        # similarities = pd.Series(list(self.edges.values()))
        # m = similarities.min()
        # k = 1 / (similarities.max() - m)
        edges = [
            graph.Edge(x, y, directed=False, weight=self.edges[(x, y)])
            # graph.Edge(x, y, directed=False, weight=k * (self.edges[(x, y)] - m))
            for x, y in self.edges
            if self.edges[(x, y)] >= threshold
        ]
        self.stream.add_edge(*edges)
