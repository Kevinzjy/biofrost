from biofrost.analysis import TaxonDB


def test_TaxonDB():
    taxon_db = TaxonDB(
        nodes_dmp="/data/public/database/taxonomy/nodes.dmp",
        names_dmp="/data/public/database/taxonomy/names.dmp",
    )
    print(taxon_db.taxon_nodes.loc[9606])
    print(taxon_db.taxon_names[9606])
    paths = taxon_db.fetch_levels(9606, "scientific name", True)
    print(paths)
