"""Functions for taxonomy related analysis"""
from collections import defaultdict, namedtuple, Counter
from biofrost.analysis import expectation_maximization

import pandas as pd
__all__ = [
    "TaxonDB",
    "SilvaDB",
]

Name = namedtuple("Name", "name_txt unique_name")
Levels = ['superkingdom', 'clade', 'phylum', 'class', 'order', 'family', 'genus', 'species']


class TaxonDB(object):
    """
    A class to parse NCBI taxonomy (https://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz) ids to levels

    Parameters
    ----------
    nodes_dmp : str
        Path of nodes.dmp downloaded from NCBI taxonomy database.
    names_dmp : str
        Path of names.dmp downloaded from NCBI taxonomy database.

    Attributes
    ----------
    taxon_nodes : pd.DataFrame
        information for each taxonomy node
    taxon_names : pd.DataFrame
        family name of the person

    Examples
    --------
    >>> from biofrost.analysis import TaxonDB
    >>> taxon_db = TaxonDB(
    >>>     nodes_dmp="/data/public/database/taxonomy/nodes.dmp",
    >>>     names_dmp="/data/public/database/taxonomy/names.dmp",
    >>> )
    >>> taxon_db.taxon_nodes.loc[9606]
    parent_tax_id        9605
    rank              species
    embl_code              HS
    div_id                  5
    i_div_flag              1
    gc_id                   1
    i_gc_flag               1
    mgc_id                  2
    i_mgc_flag              1
    h_genbank_flag          1
    h_subtree_flag          0
    comments             None
    Name: 9606, dtype: object
    >>> taxon_db.taxon_names[9606]
    {'authority': [Name(name_txt='Homo sapiens Linnaeus, 1758', unique_name='')], 'scientif
    ic name': [Name(name_txt='Homo sapiens', unique_name='')], 'genbank common name': [Name
    (name_txt='human', unique_name='')]}
    >>> taxon_db.fetch_levels(9606, "scientific name", True)
    [('no rank', 'root'), ('no rank', 'cellular organisms'), ('superkingdom', 'Eukaryota'),
     ('clade', 'Opisthokonta'), ('kingdom', 'Metazoa'), ('clade', 'Eumetazoa'), ('clade', '
    Bilateria'), ('clade', 'Deuterostomia'), ('phylum', 'Chordata'), ('subphylum', 'Craniat
    a <chordates>'), ('clade', 'Vertebrata <vertebrates>'), ('clade', 'Gnathostomata <verte
    brates>'), ('clade', 'Teleostomi'), ('clade', 'Euteleostomi'), ('superclass', 'Sarcopte
    rygii'), ('clade', 'Dipnotetrapodomorpha'), ('clade', 'Tetrapoda'), ('clade', 'Amniota'
    ), ('class', 'Mammalia'), ('clade', 'Theria <mammals>'), ('clade', 'Eutheria'), ('clade
    ', 'Boreoeutheria'), ('superorder', 'Euarchontoglires'), ('order', 'Primates'), ('subor
    der', 'Haplorrhini'), ('infraorder', 'Simiiformes'), ('parvorder', 'Catarrhini'), ('sup
    erfamily', 'Hominoidea'), ('family', 'Hominidae'), ('subfamily', 'Homininae'), ('genus'
    , 'Homo'), ('species', 'Homo sapiens')]
    """

    def __init__(self, nodes_dmp, names_dmp):
        """Init TaxonDB"""
        self.taxon_nodes = self._init_nodes(nodes_dmp)
        self.taxon_names = self._init_names(names_dmp)

    def _init_nodes(self, nodes_dmp):
        """Load NCBI taxonomy nodes.dmp

        Args:
            nodes_dmp (str): path of nodes.dmp downloaded from NCBI taxonomy database

        Returns:
            taxon_nodes: pd.DataFrame

        Description
        -----
        General information.
        Field terminator is "\t|\t"
        Row terminator is "\t|\n"

        nodes.dmp file consists of taxonomy nodes. The description for each node includes the following
        fields:
                tax_id                                  -- node id in GenBank taxonomy database
                parent tax_id                           -- parent node id in GenBank taxonomy database
                rank                                    -- rank of this node (superkingdom, kingdom, ...)
                embl code                               -- locus-name prefix; not unique
                division id                             -- see division.dmp file
                inherited div flag  (1 or 0)            -- 1 if node inherits division from parent
                genetic code id                         -- see gencode.dmp file
                inherited GC  flag  (1 or 0)            -- 1 if node inherits genetic code from parent
                mitochondrial genetic code id           -- see gencode.dmp file
                inherited MGC flag  (1 or 0)            -- 1 if node inherits mitochondrial gencode from parent
                GenBank hidden flag (1 or 0)            -- 1 if name is suppressed in GenBank entry lineage
                hidden subtree root flag (1 or 0)       -- 1 if this subtree has no sequence data yet
                comments                                -- free-text comments and citations
        """
        taxon_nodes = []
        with open(nodes_dmp, 'r') as f:
            for line in f:
                content = line.rstrip("\t|\n").split("\t|\t")
                taxon_nodes.append(content)
        taxon_nodes = pd.DataFrame(taxon_nodes)
        taxon_nodes.columns = ['tax_id', 'parent_tax_id', 'rank', 'embl_code', 'div_id', 'i_div_flag', 'gc_id', 'i_gc_flag', 'mgc_id', 'i_mgc_flag', 'h_genbank_flag', 'h_subtree_flag', 'comments']
        taxon_nodes = taxon_nodes.set_index('tax_id')
        taxon_nodes.index = taxon_nodes.index.astype(int)
        taxon_nodes['parent_tax_id'] = taxon_nodes['parent_tax_id'].astype(int)
        return taxon_nodes

    def _init_names(self, names_dmp):
        """Load NCBI taxonomy nodes.dmp

        Args:
            names_dmp (str): path of names.dmp

        Returns:
            dict: tax_id -> name_class -> [(name_txt, uniq_name)]

        Description
        -----
        Taxonomy names file (names.dmp):
                tax_id                                  -- the id of node associated with this name
                name_txt                                -- name itself
                unique name                             -- the unique variant of this name if name not unique
                name class                              -- (synonym, common name, ...)
        """
        taxon_names = defaultdict(dict)
        with open(names_dmp, 'r') as f:
            for line in f:
                content = line.rstrip("\t|\n").split("\t|\t")
                tax_id, name_txt, uniq_name, name_class = content
                taxon_names[int(tax_id)].setdefault(name_class, []).append(Name(name_txt, uniq_name))
        return taxon_names

    def _backtrace_node(self, leaf_node):
        tax_id = leaf_node.name
        parent_id = leaf_node['parent_tax_id']
        rank = leaf_node['rank']
        if tax_id not in self.taxon_names:
            raise KeyError(f"{tax_id} not found!")

        tax_name = self.taxon_names[tax_id]

        if tax_id == 1:
            return [(rank, tax_id, tax_name), ]

        parent_node = self.taxon_nodes.loc[parent_id]
        return self._backtrace_node(parent_node) + [(rank, tax_id, tax_name), ]

    def fetch_taxon_name(self, tax_id, name_class="scientific name", use_uniq=False):
        """Get name of taxonomy id

        Args:
            tax_id (int): query taxonomy id
            name_class (str, optional): type of name (authority, synonym, blast name, scientific name, common name, taxonomy id). Defaults to "scientific name".
            use_uniq (bool, optional): use the unique variant of taxonomy name. Defaults to False.

        Returns:
            str / list: taxonomy name(s)
        """

        tax_name = self.taxon_names[tax_id]
        tmp_names = []
        for name in tax_name[name_class]:
            if use_uniq and name.uniq_name != "":
                tmp_names.append(name.unique_name)
            else:
                tmp_names.append(name.name_txt)

        if len(tmp_names) == 1:
            return tmp_names[0]
        return tmp_names

    def fetch_taxon_path(self, tax_id, name_class="scientific name", use_uniq=False):
        """Get all levels of query taxonomy id

        Args:
            tax_id (int): query taxonomy id
            name_class (str, optional): Type of name (authority, synonym, blast name, scientific name, common name, taxonomy id). Defaults to "scientific name".
            use_uniq (bool, optional): Use the unique variant of taxonomy name. Defaults to False.

        Returns:
            generator object: [(level, name)] if name is unique, otherwise [(level, list(names))]
        """
        leaf_node = self.taxon_nodes.loc[tax_id]
        paths = self._backtrace_node(leaf_node)

        if name_class in ["authority", "synonym", "blast name", "scientific name", "common name"]:
            for p_lvl, p_id, p_name in paths:
                # No name class found
                if name_class not in p_name:
                    yield p_lvl, ""

                # Collect name txt
                tmp_names = []
                for name in p_name[name_class]:
                    if use_uniq and name.unique_name != "":
                        tmp_names.append(name.unique_name)
                    else:
                        tmp_names.append(name.name_txt)

                # Deal with unique name_class
                if len(tmp_names) == 1:
                    yield p_lvl, tmp_names[0]
                else:
                    yield p_lvl, tmp_names
        elif name_class == "taxonomy id":
            for p_lvl, p_id, _ in paths:
                yield p_lvl, p_id
        else:
            raise KeyError(f"Name class: {name_class} not supported!")

    def fetch_taxon_level(self, tax_id, rank="genus", name_class="scientific name", use_uniq=False):
        """Fetch specific taxonomy level of given id

        Args:
            tax_id (int): Query taxonomy id
            rank (str, optional): Rank to query {superkingdom, kingdom, phylum, order, family, genus, species}. Defaults to "genus".
            name_class (str, optional): Type of name (synonym, common name, ...). Defaults to "scientific name".
            use_uniq (bool, optional): Use the unique variant of taxonomy name. Defaults to False. Defaults to False.

        Returns:
            str: Name of specific rank, None if not exist.
        """
        assert rank in ['superkingdom', 'kingdom', 'phylum', 'order', 'family', 'genus', 'species']
        paths = self.fetch_taxon_path(tax_id, name_class, use_uniq)
        for name_rank, name_txt in paths:
            if name_rank == rank:
                return name_txt
        return None

    def batch_fetch_path(self, taxon_ids, levels=Levels):
        """Fetch taxonomy path for given ids

        Args:
            taxon_ids (iter): Iterator of query taxonomy ids

        Returns:
            pd.DataFrame: dataframe of taxnomy levels
        """
        taxon_names = pd.DataFrame({x: dict(self.fetch_taxon_path(x, name_class="scientific name")) for x in taxon_ids}).T
        taxon_col = [i for i in levels if i in taxon_names.columns]
        taxon_names = taxon_names.loc[:, taxon_col]
        return taxon_names


class SilvaDB(object):
    """
    A class to parse SILVA 16S/18S database

    Parameters
    ----------
    silva_txt : str
        Path of tax_slv_[ls]su_VERSION.txt.gz
    silva_map : str
        Path of tax_slv_[ls]su_VERSION.map.gz
    silva_acc : str
        Path of tax_slv_[ls]su_VERSION.acc_taxid.gz

    Attributes
    ----------
    slv_accid : pd.DataFrame
        Information of SILVA accession id

    Examples
    --------
    >>> from biofrost.analysis import SilvaDB
    >>> silva_db = SilvaDB(
    >>>     silva_txt="/data/public/database/SILVA/tax_slv_lsu_138.1.txt.gz",
    >>>     silva_map="/data/public/database/SILVA/tax_slv_lsu_138.1.map.gz",
    >>>     silva_acc="/data/public/database/SILVA/tax_slv_lsu_138.1.acc_taxid.gz",
    >>> )
    >>> silva_db.seq2tax('AAAA02037088.378.3887')
    taxid                                                29156
    name                                           Hypocreales
    path     Eukaryota;Amorphea;Obazoa;Opisthokonta;Nucletm...
    rank                                                 order
    Name: AAAA02037088.378.3887, dtype: object
    """
    def __init__(self, silva_txt, silva_map, silva_acc):
        self.slv_txt   = self._init_rank(silva_txt)
        self.slv_map   = self._init_map(silva_map)
        self.slv_accid = self._init_acc_taxid(silva_acc)

        self.slv_accid['name'] = self.slv_map.loc[self.slv_accid['taxid']]['name'].values
        self.slv_accid['path'] = self.slv_txt.loc[self.slv_accid['taxid']]['path'].values
        self.slv_accid['rank'] = self.slv_txt.loc[self.slv_accid['taxid']]['rank'].values

    def _init_rank(self, silva_txt):
        """Load silva txt

        tax_slv_[ls]su_VERSION.txt
        -------------------------
        These files contain taxonomic rank designations for all taxonomic paths
        used in the SILVA taxonomies. Additionally, a unique numeric identifier is
        assigned to each taxon (path). These identifiers will be mostly stable in
        upcoming SILVA releases.

        IDs used in the SSU and LSU files do not match.

        Field description:
        path:
            The full taxonomic path including the name of the group itself.
            Segments are separated with ";"
        taxid:
            numerical identifier
        rank:
            The rank designation.
        remark:
            Can be empty ('') or a or w.
            a: Marks taxa of environmental origin. That is, taxa containing no
            sequence coming from a cultivated organism.
            w: Marks taxa scheduled for revision in the next release.
        release:
            The SILVA release version
        """
        slv_txt = pd.read_csv(silva_txt, compression="gzip", sep="\t", header=None)
        slv_txt.columns = ['path', 'taxid', 'rank', 'remark', 'release']
        slv_txt = slv_txt.set_index('taxid')
        return slv_txt

    def _init_map(self, silva_map):
        """Load silva taxonomy map"""
        slv_map = pd.read_csv(silva_map, compression="gzip", sep="\t", header=None)
        slv_map.columns = ['taxid', 'name', 'mark', 'parent']
        slv_map = slv_map.set_index('taxid')
        return slv_map

    def _init_acc_taxid(self, silva_acc):
        """Load silva acc taxid

        tax_TAXNAME_[ls]su_VERSION.acc_taxid
        ----------------------------
        Mapping of 'SILVA' sequence IDs (<INSDC primary accession>.<start>.<stop>)
        used in FASTA files to the numeric SILVA taxid (MEGAN compatible).
        """
        slv_accid = pd.read_csv(silva_acc, compression="gzip", sep="\t", header=None)
        slv_accid.columns = ['sequence', 'taxid']
        slv_accid = slv_accid.set_index('sequence')
        return slv_accid

    def seq2tax(self, sequence_name):
        """Get taxonomy information for SILVA sequence name

        Args:
            sequence_name (str): query sequence id

        Returns:
            pd.Series: (taxid, name, path, rank)
        """
        return self.slv_accid.loc[sequence_name]


def load_centrifuge(clf_file, taxon_db, lvl_name="genus", aln_perc=.1, use_EM=False):
    """Load centrifuge classification output

    Args:
        clf_file (_type_): path to centrifuge classification output
        taxon_db (TaxonDB): see biofrost.parser.TaxonDB
        lvl_name (str, optional): taxonomy level (superkingdom, clade, phylum, class, order, family, genus, species). Defaults to "genus".
        aln_perc (float, optional): threshold for filter high confidence alignment. Defaults to .1.
        use_EM (bool, optional): whether use EM to calculate taxonomy abundance. Defaults to False.

    Returns:
        pd.Series: abundance of each taxonomy ids
    """

    # Load centrifuge output
    df = pd.read_csv(clf_file, sep="\t")
    df = df.loc[df['hitLength']/df['queryLength'] >= aln_perc]
    df = df.loc[df['taxID'].isin(taxon_db.taxon_nodes.index)]
    df['rank'] = taxon_db.taxon_nodes.loc[df['taxID']]['rank'].values

    # Convert taxonomy IDs to full level paths
    taxon_levels = pd.DataFrame({taxID: dict(taxon_db.fetch_taxon_path(taxID, name_class="taxonomy id")) for taxID in df['taxID'].unique()}).T
    taxon_levels = taxon_levels.loc[:, taxon_levels.columns.isin(Levels)]

    # Merge dataframe
    df = pd.concat([df.reset_index(drop=True), taxon_levels.loc[df['taxID']].reset_index(drop=True)], axis=1)

    # Calculate abundance
    lvl_df = df.loc[~df[lvl_name].isna()]
    lvl_df = lvl_df[['readID', lvl_name]].drop_duplicates()
    if use_EM:
        abundance, _ = expectation_maximization(lvl_df['readID'].values, lvl_df[lvl_name].values, noise_threshold=1/lvl_df.shape[0], verbose=True)
    else:
        n_hits = Counter(lvl_df['readID'])
        abundance = defaultdict(int)
        for r, c in zip(lvl_df['readID'], lvl_df[lvl_name]):
            abundance[c] += 1 / n_hits[r]
        abundance = pd.Series(abundance).sort_values(ascending=False)

    return abundance


def load_diamond2(clf_file, taxon_db, lvl_name="genus", max_e=1e-6, use_EM=False):
    """Load diamond2 alignment format 102 output

    Args:
        clf_file (_type_): path of diamond2 output
        taxon_db (TaxonDB): see biofrost.parser.TaxonDB
        lvl_name (str, optional): taxonomy level (superkingdom, clade, phylum, class, order, family, genus, species). Defaults to "genus".
        max_e (float, optional): evalue threshold. Defaults to 1e-6.
        use_EM (bool, optional): whether use EM to calculate taxonomy abundance. Defaults to False.

    Returns:
        pd.Series: abundance of each taxonomy ids
    """

    df = pd.read_csv(clf_file, sep="\t", header=None)
    df.columns = ['readID', 'taxID', 'evalue']
    df = df.loc[(df['evalue'] <= max_e) & (df['taxID'] != 0)]
    df['rank'] = taxon_db.taxon_nodes.loc[df['taxID']]['rank'].values

    # Convert taxonomy IDs to full level paths
    taxon_levels = pd.DataFrame({taxID: dict(taxon_db.fetch_taxon_path(taxID, name_class="taxonomy id")) for taxID in df['taxID'].unique()}).T
    taxon_levels = taxon_levels.loc[:, taxon_levels.columns.isin(Levels)]

    # Merge dataframe
    df = pd.concat([df.reset_index(drop=True), taxon_levels.loc[df['taxID']].reset_index(drop=True)], axis=1)

    # Calculate abundance
    lvl_df = df.loc[~df[lvl_name].isna()]
    lvl_df = lvl_df[['readID', lvl_name]].drop_duplicates()
    if use_EM:
        abundance, _ = expectation_maximization(lvl_df['readID'].values, lvl_df[lvl_name].values, noise_threshold=1/lvl_df.shape[0], verbose=True)
    else:
        n_hits = Counter(lvl_df['readID'])
        abundance = defaultdict(int)
        for r, c in zip(lvl_df['readID'], lvl_df[lvl_name]):
            abundance[c] += 1 / n_hits[r]
        abundance = pd.Series(abundance).sort_values(ascending=False)
    return abundance


def load_kraken2(clf_file, taxon_db, lvl_name="genus"):
    """Load kraken2 kreport output

    Args:
        clf_file (Path): path of kraken2 kreport output
        taxon_db (TaxonDB): see biofrost.parser.TaxonDB
        lvl_name (str, optional): taxonomy level (superkingdom, clade, phylum, class, order, family, genus, species). Defaults to "genus".

    Returns:
        pd.Series: abundance of each taxonomy ids
    """

    df = pd.read_csv(clf_file, sep="\t", header=None)
    df.columns = ["Abundance", "Reads", "Assigned", "Class", "taxID", "Name"]
    df = df.loc[df['taxID'] != 0]

    df['rank'] = taxon_db.taxon_nodes.loc[df['taxID']]['rank'].values

    taxon_levels = pd.DataFrame({taxID: dict(taxon_db.fetch_taxon_path(taxID, name_class="taxonomy id")) for taxID in df['taxID'].unique()}).T
    taxon_levels = taxon_levels.loc[:, taxon_levels.columns.isin(Levels)]

    df = pd.concat([df.reset_index(drop=True), taxon_levels.loc[df['taxID']].reset_index(drop=True)], axis=1)

    lvl_df = df.loc[~df[lvl_name].isna()]
    lvl_df = lvl_df.loc[lvl_df['rank'] == lvl_name]

    abundance = lvl_df.set_index('taxID')["Abundance"]

    return abundance