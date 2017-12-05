# GEO

Library for working with Gene Expression Omnibus (GEO) :gem:

## Dependencies

```sh
pip install pandas
pip install geoparse
```

## Get started

```python
from geo.get_and_parse_geo_data import get_and_parse_geo_data

geo_dict = get_and_parse_geo_data('GSE26375', directory_path='.')

print(geo_dict.keys())

print(geo_dict['id_x_sample'])
print(geo_dict['gene_x_sample'])
print(geo_dict['information_x_sample'])
```

Look at `notebook/geo.ipynb` for demonstration.

## Development

If you find a bug or have any trouble with this library, please [submit an issue](https://github.com/KwatME/geo/issues) and I'll take care of it.
