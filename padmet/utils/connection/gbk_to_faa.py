# -*- coding: utf-8 -*-
"""
Description:
    convert GBK to FAA with Bio package

::

    usage:
        padmet gbk_to_faa    --gbk=FILE --output=FILE [--qualifier=STR] [-v]

    option:
        -h --help    Show help.
        --gbk=FILE    path to the gbk file.
        --output=FILE    path to the output, a FAA file.
        --qualifier=STR    the qualifier of the gene id [default: locus_tag].
        -v   print info
"""
import docopt

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet.IUPAC import protein


def command_help():
    """
    Show help for analysis command.
    """
    print(docopt.docopt(__doc__))


def gbk_to_faa_cli(command_args):
    args = docopt.docopt(__doc__, argv=command_args)
    gbk_file = args["--gbk"]
    output = args["--output"]
    qualifier = args["--qualifier"]
    verbose = args["-v"]
    gbk_to_faa(gbk_file, output, qualifier, verbose)


def gbk_to_faa(gbk_file, output, qualifier='locus_tag', verbose=True):
    """
    convert GBK to FAA with Bio package

    Parameters
    ----------
    gbk_file: str
        path to the gbk file
    output: str
        path to the output, a FAA file
    qualifier: str
        he qualifier of the gene id
    verbose: bool
        if True print information
    """
    fasta_records = []
    fasta_ids = {}
    with open(gbk_file, "r") as gbk:
        for seq_record in SeqIO.parse(gbk, "genbank"):
            seq_feature_cds = (seq_feature for seq_feature in seq_record.features if seq_feature.type == "CDS")
            for seq_feature in seq_feature_cds:
                try:
                    fasta_id = seq_feature.qualifiers[qualifier][0]
                    if fasta_id not in fasta_ids:
                        fasta_record = SeqRecord(Seq(seq_feature.qualifiers['translation'][0], protein), id=fasta_id, description=fasta_id)
                        fasta_records.append(fasta_record)
                        fasta_ids[fasta_id] = 1
                    else:
                        fasta_ids[fasta_id] += 1
                        isoform_id = fasta_id + '_isoform' + str(fasta_ids[fasta_id])
                        fasta_record = SeqRecord(Seq(seq_feature.qualifiers['translation'][0], protein), id=isoform_id, description=fasta_id)
                        fasta_records.append(fasta_record)
                except KeyError:
                    if verbose:
                        try:
                            print("locus without Translation: "+seq_feature.qualifiers['locus_tag'][0])
                        except KeyError:
                            print("locus without Translation: "+seq_feature.qualifiers.get('gene',["Unknown"])[0])

    SeqIO.write(fasta_records, output, "fasta")


