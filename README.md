[![Build Status][ci-badge]][ci-link]
[![Coverage Status][cov-badge]][cov-link]
[![Docs status][docs-badge]][docs-link]
[![PyPI version][pypi-badge]][pypi-link]

# aiida-hydrogen-restorer

AiiDA plugin package with workflows to detect missing hydrogens in crystal structures from external databases and use DFT simulations to restore them in the correct positions.

## Installation

```shell
pip install aiida-hydrogen-restorer
```

## Development

```shell
git clone https://github.com/epfl_theos/aiida-hydrogen-restorer .
cd aiida-hydrogen-restorer
pip install --upgrade pip
pip install -e .[pre-commit,testing]  # install extra dependencies
pre-commit install  # install pre-commit hooks
pytest -v  # discover and run all tests
```

See the [developer guide](http://aiida-hydrogen-restorer.readthedocs.io/en/latest/developer_guide/index.html) for more information.

## License

MIT

## Contact

mbercx@gmail.com


[ci-badge]: https://github.com/epfl_theos/aiida-hydrogen-restorer/workflows/ci/badge.svg?branch=master
[ci-link]: https://github.com/epfl_theos/aiida-hydrogen-restorer/actions
[cov-badge]: https://coveralls.io/repos/github/epfl_theos/aiida-hydrogen-restorer/badge.svg?branch=master
[cov-link]: https://coveralls.io/github/epfl_theos/aiida-hydrogen-restorer?branch=master
[docs-badge]: https://readthedocs.org/projects/aiida-hydrogen-restorer/badge
[docs-link]: http://aiida-hydrogen-restorer.readthedocs.io/
[pypi-badge]: https://badge.fury.io/py/aiida-hydrogen-restorer.svg
[pypi-link]: https://badge.fury.io/py/aiida-hydrogen-restorer
