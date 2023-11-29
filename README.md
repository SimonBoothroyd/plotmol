# Interactive Plots with Molecular Annotations

[![Test Status](https://github.com/simonboothroyd/plotmol/actions/workflows/ci.yaml/badge.svg?branch=main)](https://github.com/simonboothroyd/plotmol/actions/workflows/ci.yaml)
[![codecov](https://codecov.io/gh/SimonBoothroyd/plotmol/branch/main/graph/badge.svg?token=Aa8STE8WBZ)](https://codecov.io/gh/SimonBoothroyd/plotmol)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

This framework provides a very simple set of functionalities for creating interactive plots that are annotated with 
molecular structures using the [`bokeh` library](https://docs.bokeh.org/en/latest/index.html).

## Installation

The package and its dependencies can be installed using the `conda` package manager:

```shell
conda install -c conda-forge plotmol
```

## Getting Started

See the [examples](examples) directory for an example of this framework in use.

## Development

A development conda environment can be created and activated by running:

```shell
make env
conda activate plotmol
```

The environment will include all example and development dependencies, including linters and testing apparatus.

Unit / example tests can be run using:

```shell
make test
make test-examples
```

The codebase can be formatted by running:

```shell
make format
```

or checked for lint with:

```shell
make lint
```

## License

The main package is release under the [MIT license](LICENSE). 

## Copyright

Copyright (c) 2023, Simon Boothroyd
