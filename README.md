# Visualization for Chemical Multi-Spectral Space

Scripts for exploring the chemical space with spectral (NMR, IR) data. See the [complete report](./report/report.pdf).

## Requirements

- Python 3.x
- Dependent Libraries (See [requirements.txt](./requirements.txt)):
    - GephiStreamer
    - NumPy
    - Pandas
- [Gephi](https://gephi.org/)

You can install the dependencies by:

```bash
$ pip install -r requirements.txt
```

## Installation

In Gephi, install a plugin called [Graph Streaming](https://gephi.org/plugins/#/plugin/graphstreaming).

Download or clone this repository.

## Usage

The scripts work with Gephi. The main pipeline is as following:

1. In Gephi, keep a workspace open and set a name (e.g. `sample`) for it.

2. Start a stream server for that workspace. (Just right click `Master Server` under the `Streaming` tab and start.)

3. Run the scripts (visualize the sample):

```bash
$ export PYTHONPATH=./src  # to import `chemspace` module.
$ python3 scripts/visualize_main.py ./data/sample/m150_nmr.p ./data/sample/m150_ir.p --gephi-workspace sample
```

It is recommended to provide arguments to specify the cache files:

```bash
$ python3 scripts/visualize_main.py ./data/sample/m150_nmr.p ./data/sample/m150_ir.p --gephi-workspace sample \
  --fingerprints-cache cache-sample.h5 --graph-cache cache-sample.p 
```

More usage information can be found by:

```bash
$ python3 scripts/visualize_main.py --help
```

4. Under the `Layout` tab in Gephi, choose `ForceAtlas 2` and `Run` with the default parameters.
Stop when the nodes in the graph are distributed in an equilibrium.

5. Under the `Appearance-Nodes-Ranking` tab, select `d` and set expected colors. Go to `Preview`, and `Refresh`
with whatever parameters you want. It would also be better if the background is set to black or dark gray.

The Gephi file for the sample can be found in the [release page](https://github.com/sunoru/VisChemSpace/releases).

## License

[MIT](./LICENSE.md)
