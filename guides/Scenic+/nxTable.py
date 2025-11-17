import pandas as pd
import anndata
import scanpy as sc

def _format_df_nx(df, key, var):
    """
    A helper function to format differential test results
    """
    df.index = df['names']
    df = pd.DataFrame(df['logfoldchanges'])
    df.columns = [var+'_Log2FC_'+key]
    df.index.name = None
    return df

def _get_log2fc_nx(scplus_obj: 'SCENICPLUS',
                  variable,
                  features,
                  contrast = "gene" #: Optional[str] = 'gene'
                  ):
    """
    A helper function to derive log2fc changes
    """
    if contrast == 'gene':
        adata = anndata.AnnData(X=scplus_obj.X_EXP, obs=pd.DataFrame(
            index=scplus_obj.cell_names), var=pd.DataFrame(index=scplus_obj.gene_names))
    if contrast == 'region':
        adata = anndata.AnnData(X=scplus_obj.X_ACC.T, obs=pd.DataFrame(
            index=scplus_obj.cell_names), var=pd.DataFrame(index=scplus_obj.region_names))
    adata.obs = pd.DataFrame(scplus_obj.metadata_cell[variable])
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    adata = adata[:, features]
    sc.tl.rank_genes_groups(
        adata, variable, method='wilcoxon', corr_method='bonferroni')
    groups = adata.uns['rank_genes_groups']['names'].dtype.names
    diff_list = [_format_df_nx(sc.get.rank_genes_groups_df(
        adata, group=group), group, variable) for group in groups]
    return pd.concat(diff_list, axis=1)

def create_nx_tables(scplus_obj: 'SCENICPLUS',
                     eRegulon_metadata_key = "eRegulon_metadata", #: str ='eRegulon_metadata',
                     subset_eRegulons = None, #,: List = None,
                     subset_regions = None, #: List = None,
                     subset_genes = None, #: List = None,
                     add_differential_gene_expression = False, #: bool = False,
                     add_differential_region_accessibility = False, #: bool = False,
                     differential_variable = [] ): #: List =[]):
    """
    A function to format eRegulon data into tables for plotting eGRNs.
    
    Parameters
    ---------
    scplus_obj: SCENICPLUS
        A SCENICPLUS object with eRegulons
    eRegulon_metadata_key: str, optional
        Key where the eRegulon metadata dataframe is stored
    subset_eRegulons: list, optional
        List of eRegulons to subset
    subset_regions: list, optional
        List of regions to subset
    subset_genes: list, optional
        List of genes to subset
    add_differential_gene_expression: bool, optional
        Whether to calculate differential gene expression logFC for a given variable
    add_differential_region_accessibility: bool, optional
        Whether to calculate differential region accessibility logFC for a given variable
    differential_variable: list, optional
        Variable to calculate differential gene expression or region accessibility.
    Return
    ---------
    A dictionary with edge feature tables ('TF2G', 'TF2R', 'R2G') and node feature tables ('TF', 'Gene', 'Region')
    """
    import pandas as pd
    er_metadata = scplus_obj.uns[eRegulon_metadata_key].copy()
    if subset_eRegulons is not None:
        # subset_eRegulons = [x + '_[^a-zA-Z0-9]' for x in subset_eRegulons]
        er_metadata = er_metadata[er_metadata['Region_signature_name'].str.contains(
            '|'.join(subset_eRegulons))]
    if subset_regions is not None:
        er_metadata = er_metadata[er_metadata['Region'].isin(subset_regions)]
    if subset_genes is not None:
        er_metadata = er_metadata[er_metadata['Gene'].isin(subset_genes)]
    nx_tables = {}
    nx_tables['Edge'] = {}
    nx_tables['Node'] = {}
    # Generate edge tables
    r2g_columns = [x for x in er_metadata.columns if 'R2G' in x]
    tf2g_columns = [x for x in er_metadata.columns if 'TF2G' in x]
    nx_tables['Edge']['TF2R'] = er_metadata[er_metadata.columns.difference(
        r2g_columns + tf2g_columns)].drop('Gene', axis=1).drop_duplicates()
    nx_tables['Edge']['TF2R'] = nx_tables['Edge']['TF2R'][['TF', 'Region'] +
                                                          nx_tables['Edge']['TF2R'].columns.difference(['TF', 'Region']).tolist()]
    nx_tables['Edge']['R2G'] = er_metadata[er_metadata.columns.difference(
        tf2g_columns)].drop('TF', axis=1).drop_duplicates()
    nx_tables['Edge']['R2G'] = nx_tables['Edge']['R2G'][['Region', 'Gene'] +
                                                        nx_tables['Edge']['R2G'].columns.difference(['Region', 'Gene']).tolist()]
    nx_tables['Edge']['TF2G'] = er_metadata[er_metadata.columns.difference(
        r2g_columns)].drop('Region', axis=1).drop_duplicates()
    nx_tables['Edge']['TF2G'] = nx_tables['Edge']['TF2G'][['TF', 'Gene'] +
                                                          nx_tables['Edge']['TF2G'].columns.difference(['TF', 'Gene']).tolist()]
    # Generate node tables
    tfs = list(set(er_metadata['TF']))
    nx_tables['Node']['TF'] = pd.DataFrame(
        'TF', index=tfs, columns=['Node_type'])
    nx_tables['Node']['TF']['TF'] = tfs
    genes = list(set(er_metadata['Gene']))
    genes = [x for x in genes if x not in tfs]
    nx_tables['Node']['Gene'] = pd.DataFrame(
        'Gene', index=genes, columns=['Node_type'])
    nx_tables['Node']['Gene']['Gene'] = genes
    regions = list(set(er_metadata['Region']))
    nx_tables['Node']['Region'] = pd.DataFrame(
        'Region', index=regions, columns=['Node_type'])
    nx_tables['Node']['Region']['Region'] = regions
    # Add gene logFC
    if add_differential_gene_expression is True:
        for var in differential_variable:
            nx_tables['Node']['TF'] = pd.concat([nx_tables['Node']['TF'], _get_log2fc_nx(
                scplus_obj, var, nx_tables['Node']['TF'].index.tolist(), contrast='gene')], axis=1)
            nx_tables['Node']['Gene'] = pd.concat([nx_tables['Node']['Gene'], _get_log2fc_nx(
                scplus_obj, var, nx_tables['Node']['Gene'].index.tolist(), contrast='gene')], axis=1)
    if add_differential_region_accessibility is True:
        for var in differential_variable:
            nx_tables['Node']['Region'] = pd.concat([nx_tables['Node']['Region'], _get_log2fc_nx(
                scplus_obj, var, nx_tables['Node']['Region'].index.tolist(), contrast='region')], axis=1)
    return nx_tables
