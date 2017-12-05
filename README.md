# GEO

Library for working with Gene Expression Omnibus (GEO) :gem:

Example

```python
geo_dict = get_and_parse_geo_data('GSE26375', directory_path='.')

print(geo_dict.keys())

print(geo_dict['id_x_sample'])
print(geo_dict['gene_x_sample'])
print(geo_dict['information_x_sample'])
```
