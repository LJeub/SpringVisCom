# SpringVisCom
SpringVisCom produces network layouts emphasising pre-specified community structure

## Installation

Download the latest release from [here](https://github.com/LJeub/SpringVisCom/releases/latest) and add the 
code to your MATLAB/OCTAVE search path (e.g. using `addpath(genpath(SpringVisCom))` from the parent directory).

## Usage

### Single-layer networks

To obtain node coordinates for a single-layer network use
```Octave
xy = SpringVisCom(A, S)
```
The functions assume that the network is given as an adjacency matrix `A` such that `A(i, j)=w` if nodes `i` and `j`
are connected by an edge of weight `w` (`w=1` for an unweighted network) and `A(i, j)=0` if nodes `i` and `j` are 
not connected. It additionally takes community assingments `S` for the nodes (e.g. obtained using 
[GenLouvain](https://github.com/GenLouvain/GenLouvain)) as input. If no community assingments are available, set `S=[]`.

To plot the network use
```Octave
GraphPlot(xy, A)
```

.

For more information and available options use `help('SpringVisCom')` and `help('GraphPlot')` from the MATLAB/OCTAVE prompt
and have a look at the [usage examples](https://ljeub.github.io/code/netvis.html).

### Multilayer networks

This package also includes support for laying out and plotting multilayer networks provided as a cell array of adjacency matrices. 
See `help('MultilayerSpringVisCom')` and `help('MultilayerGraphPlot')` for more information.
