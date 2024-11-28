# Minkowski

Minkowski is a Python library for computing dense lattice packings of equal hyperspheres. 

## Installation
`Python >= 3.8` is required. Note that **only Linux is supported at this time**. Please contact the authors if you would like to see support for additional systems (or submit a pull request)!

```
pip install minkowski
```

## Usage

```python
import minkowski as mx

o = mx.HeuristicSearch("hkz")
o.cube(3, 999, 1)
```

Please refer to the Wiki for full documentation.

## Contributing

Pull requests are welcome. Alternatively, to request a feature or report a bug, please feel free to open an issue.

## References
<a id="1">[1]</a> 
Edmund Hlawka, Rudolf Taschner, and Johannes Schoißengeier. Geometric and Analytic Number Theory. Universitext. Springer Berlin Heidelberg, Berlin, Heidelberg, 1991.

