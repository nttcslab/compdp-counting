# CompDP: A Framework for Simultaneous Subgraph Counting under Connectivity Constraints

This repository includes the codes to reproduce the results of experiments in the paper "CompDP: A Framework for Simultaneous Subgraph Counting under Connectivity Constraints."

## Requirements

We use [TdZdd](https://github.com/kunisura/TdZdd) and [boost C++ libraries](https://www.boost.org/) for implementing both proposed and baseline methods. Before building our code, you must place the header files of [TdZdd](https://github.com/kunisura/TdZdd) into this directory. More specifically, all header files of [TdZdd/include/tdzdd](https://github.com/kunisura/TdZdd/tree/master/include/tdzdd) must be placed on `tdzdd/` directory; e.g. `tdzdd/DdEval.hpp`, `tdzdd/DdSpec.hpp`, etc. Additionally, you should rewrite `path/to/boost/` in Line 9 of `CMakeLists.txt` so that this is the boost root directory of your environment.

After that, if your environment has CMake version >=3.8, you can build all codes with the following commands:

```shell
(moving to src/ directory)
mkdir release
cd release
cmake -DCMAKE_BUILD_TYPE=Release ..
make
```

After running this, all binaries are generated in `release/` directory.

### Verified environments

We verified that the building process of our codes and the commands presented below worked fine in the following macOS and Linux environments:

- macOS Big Sur 11.2.1 + Apple clang 12.0.0
- CentOS 7.9 + gcc 4.8.5

## How to reproduce experimental results

After building our code, the following binaries are generated: `main` and `baseline`. `main` implements our proposed method, while `baseline` implements the baseline method.

All data used in our experiments are in `data.tar.gz`. It includes the data of the following graphs:
- Grid graphs: `grid8x8.txt`, `grid8x16.txt`, `grid8x24.txt`, `grid8x32.txt`, `grid9x9.txt`, `grid10x10.txt`, `grid11x11.txt`, `grid12x12.txt`, `grid13x13.txt`
- Rocketfuel: `1221.dat.txt` (Rocketfuel-1221), `1755.dat.txt` (Rocketfuel-1755), `6461.dat.txt` (Rocketfuel-6461)
- [Romegraph](http://graphdrawing.org/download/rome-graphml.tgz): `grafo<id>.100.graphml.txt`

These graphs are described as an edgelist, and the edges are ordered according to the edge ordering described in the paper. For each graph, the following data are included:

- `<graphname>.txt`: the edgelist of the graph. Used as `[graph_file]`.
- `<graphname>.txt.1`: the center vertex used for `[vertex_file]` of cycle problem.
- `<graphname>.txt.2`: the most distant pair of vertices used for `[vertex_file]` of path problem.
- `<graphname>.txt.4`: the most distant four vertices used for `[vertex_file]` of Steiner tree and rooted spanning forest (RSF) problems.

### Experiments

The proposed method can be executed by the following command:

```shell
./main [mode] [graph_file] [vertex_file]
```

After execution, the program solves the subgraph counting for every vertex and output each count value. The four problem settings in the experiment can be treated by setting `[mode]`;

- `./main 0 <graphname>.txt <graphname>.txt.2` solves the path counting.
- `./main 1 <graphname>.txt <graphname>.txt.1` solves the cycle counting.
- `./main 2 <graphname>.txt <graphname>.txt.4` solves the Steiner tree counting.
- `./main 3 <graphname>.txt <graphname>.txt.4` solves the RST counting.

The baseline method can also be executed in the same way as follows:

```shell
./baseline [mode] [graph_file] [vertex_file]
```

## License

This software is released under the NTT license, see `LICENSE.pdf`.
