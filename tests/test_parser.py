from pathlib import Path
test_dir = Path(__file__).parent


def test_PAF():
    from biofrost.parser import yield_paf
    paf_file = test_dir / "parser" / "test.paf"
    for hit in yield_paf(paf_file):
        print(hit)
        break


def test_batch_PAF():
    from biofrost.parser import batch_yield_paf
    paf_file = test_dir / "parser" / "test.paf"
    for read_id, chunk in batch_yield_paf(paf_file):
        print(read_id)
        print(chunk)


def test_gtf():
    from biofrost.parser import yield_gff
    gtf_file = test_dir / "parser" / "chrM.gtf"
    print(len([record for record in yield_gff(gtf_file)]))


def test_gff():
    from biofrost.parser import yield_gff
    gff_file = test_dir / "parser" / "chrM.gff3"
    print(len([record for record in yield_gff(gff_file, is_gtf=False)]))


# def test_TaxonDB():
#     from biofrost.parser import TaxonDB
#     taxon_db = TaxonDB(
#         nodes_dmp="/data/public/database/taxonomy/nodes.dmp",
#         names_dmp="/data/public/database/taxonomy/names.dmp",
#     )
#     print(taxon_db.taxon_nodes.loc[9606])
#     print(taxon_db.taxon_names[9606])
#     paths = taxon_db.fetch_levels(9606, "scientific name", True)
#     print(paths)


# def test_SilvaDB():
#     from biofrost.parser import SilvaDB
#     silva_db = SilvaDB(
#         silva_txt="/data/public/database/SILVA/tax_slv_lsu_138.1.txt.gz",
#         silva_map="/data/public/database/SILVA/tax_slv_lsu_138.1.map.gz",
#         silva_acc="/data/public/database/SILVA/tax_slv_lsu_138.1.acc_taxid.gz",
#     )
#     print(silva_db.seq2tax('AAAA02037088.378.3887'))
#     assert False


def test_read_fastx():
    from biofrost.parser import read_fastx
    assert read_fastx(test_dir / "parser" / "test.fa")['seq1'].seq == "AAAAATTTTT"
    assert read_fastx(test_dir / "parser" / "test.fa.gz")['seq1'].seq == "AAAAATTTTT"
    assert read_fastx(test_dir / "parser" / "test.fq")['seq2'].seq == "ACGATCGAC"
    assert read_fastx(test_dir / "parser" / "test.fq.gz")['seq2'].seq == "ACGATCGAC"
