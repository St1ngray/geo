from os.path import abspath, expanduser

import GEOparse
from pandas import concat


def get_and_parse_geo_data(geo_id, directory_path='.'):
    """
    Get and parse GEO data.
    Arguments:
        geo_id (str):
        directory_path (str):
    Returns:
        dict: GEO dict
            {
                information_x_sample: DataFrame (n_information, n_samples),
                id_x_sample: DataFrame (n_ids, n_samples),
                id_to_gene_symbol: dict (n_ids),
                gene_x_sample: DataFrame (n_genes, n_samples),
            }
    """

    directory_path = abspath(expanduser(directory_path))
    print('Downloading GEO data in {} ...'.format(directory_path))
    gse = GEOparse.get_GEO(geo=geo_id, destdir=directory_path)

    print('Title: {}'.format(gse.get_metadata_attribute('title')))
    print('N samples: {}'.format(len(gse.get_metadata_attribute('sample_id'))))

    geo_dict = {}

    # Make ID-x-sample
    values = []
    for sample_id, gsm in gse.gsms.items():
        print('{} ...'.format(sample_id))

        sample_table = gsm.table
        sample_table.columns = sample_table.columns.str.lower().str.replace(
            ' ', '_')
        sample_values = sample_table.set_index('id_ref').squeeze()
        sample_values.name = sample_id
        values.append(sample_values)

    geo_dict['id_x_sample'] = concat(
        values, axis=1).sort_index().sort_index(axis=1)
    print('id_x_sample.shape: {}'.format(geo_dict['id_x_sample'].shape))

    # Make ID-to-gene-symbol dict
    id_to_gene_symbol = None
    for platform_id, gpl in gse.gpls.items():
        print('{} ...'.format(platform_id))

        platform_table = gpl.table
        platform_table.columns = platform_table.columns.str.lower(
        ).str.replace(' ', '_')
        platform_table.set_index('id', inplace=True)

        if 'gene_symbol' not in platform_table.columns:
            # Make gene_symbol column

            if 'gene_assignment' in platform_table.columns:
                gene_symbols = []
                for a in platform_table['gene_assignment']:
                    try:
                        gene_symbols.append(a.split('//')[1].strip())
                    except IndexError as e:
                        print('{}: {}'.format(e, a))
                        gene_symbols.append('NO GENE NAME')
                platform_table['gene_symbol'] = gene_symbols

            elif 'oligoset_genesymbol' in platform_table.columns:
                platform_table['gene_symbol'] = platform_table[
                    'oligoset_genesymbol']

        if 'gene_symbol' in platform_table:
            id_to_gene_symbol = platform_table['gene_symbol'].dropna()
            print('id_to_gene_symbol:\n{}'.format(id_to_gene_symbol))
            geo_dict['id_to_gene_symbol'] = id_to_gene_symbol.to_dict()

            # Make gene-x-sample
            gene_x_sample = geo_dict['id_x_sample'].copy()
            id_to_gene_symbol = geo_dict['id_to_gene_symbol']

            gene_x_sample.index = geo_dict['id_x_sample'].index.map(
                lambda i: id_to_gene_symbol.get(str(i), 'NO GENE NAME'))

            gene_x_sample.drop('NO GENE NAME', inplace=True)

            gene_x_sample.index.name = 'gene_symbol'

            geo_dict['gene_x_sample'] = gene_x_sample.sort_index().sort_index(
                axis=1)
            print('gene_x_sample.shape: {}'.format(geo_dict['gene_x_sample']
                                                   .shape))

        else:
            geo_dict['id_to_gene_symbol'] = None
            geo_dict['gene_x_sample'] = None
            print(
                '\tgene_symbol is not a GPL column ({}); IDs may be already gene symbols.'.
                format(', '.join(platform_table.columns)))

        # Make information_x_sample
        geo_dict['information_x_sample'] = gse.phenotype_data.T
        print('information:\n\t{}'.format('\n\t'.join(geo_dict[
            'information_x_sample'].index)))

    return geo_dict
