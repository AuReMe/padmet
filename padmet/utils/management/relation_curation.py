#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
:

    usage:
        padmet relation_curation    --padmet=FILE --id_in=STR [--type=STR] [-v]
        padmet relation_curation    --padmet=FILE --id_out=STR [--type=STR] [-v]
        padmet relation_curation    --padmet=FILE --id_in=STR [--type=STR] --to-remove==STR --output=FILE [-v]
        padmet relation_curation    --padmet=FILE --id_in=STR --id_out=STR [--type=STR] --to-remove==STR --output=FILE [-v]

    option:
        -h --help    Show help.
        --model_metabolic=FILE    pathname to the metabolic network of the model (sbml).
        --study_metabolic=FILE    ****.
        --inp=FILE    ****.
        --omcl=FILE    ****.
        --output=FILE    ****.
        -v   print info.
"""
import docopt

from padmet.classes import PadmetSpec


def command_help():
    """
    Show help for analysis command.
    """
    print(docopt.docopt(__doc__))


def relation_curation_cli(command_args):
    args = docopt.docopt(__doc__, argv=command_args)
    padmet_file = args["--padmet"]
    id_in = args["--id_in"]
    id_out = args["--id_out"]
    _type = args["--type"]
    output = args["--output"]
    verbose = args["-v"]
    to_remove = args["--to-remove"]
    padmet = PadmetSpec(padmet_file)
    get_relations(padmet=padmet, id_in=id_in, id_out=id_out, _type=_type, verbose=verbose)

    if to_remove:
        if to_remove != "all":
            to_remove = to_remove.split(";")
        get_relations(padmet=padmet, id_in=id_in, id_out=id_out, _type=_type, to_remove=to_remove, output=output, verbose=False)


def get_relations(padmet, id_in=None, id_out=None, _type=None, to_remove=None, output=None, verbose=False):
    """
    """
    if id_in:
        if id_in not in padmet.dicOfNode.keys() :
            raise KeyError("id %s not in padmet" %id_in)
        if id_in not in padmet.dicOfRelationIn.keys() :
            raise KeyError("id %s is in padmet but is not associated to any relation as id_in" %id_in)
        
        if _type:
            rlts = [rlt for rlt in padmet.dicOfRelationIn[id_in] if rlt.type == _type]
            print("%s relations found with id_in: %s and type: %s" %(len(rlts),id_in, _type))
        else:
            rlts = padmet.dicOfRelationIn[id_in]
            print("%s relations found with id_in: %s" %(len(rlts),id_in))
        [print("INDEX: %s; RELATION: %s" %(k,v.toString())) for k,v in (enumerate(rlts))]

    elif id_out:
        if id_out not in padmet.dicOfNode.keys() :
            raise KeyError("id %s not in padmet" %id_out)
        if id_out not in padmet.dicOfRelationOut.keys() :
            raise KeyError("id %s is in padmet but is not associated to any relation as id_out" %id_in)

        if _type:
            rlts = [rlt for rlt in padmet.dicOfRelationOut[id_out] if rlt.type == _type]
            print("%s relations found with id_in: %s and type: %s" %(len(rlts),id_out, _type))
        else:
            rlts = padmet.dicOfRelationIn[id_out]
            print("%s relations found with id_in: %s" %(len(rlts),id_out))
        [print("INDEX: %s; RELATION: %s" %(k,v.toString())) for k,v in (enumerate(rlts))]


    if not to_remove:
        return rlts
    elif to_remove == "all":
        print("Removing all relations")
        for rlt in rlts:
            padmet._delRelation(rlt)
        padmet.generateFile(output)
    else:
        for index, rlt in enumerate(rlts):
            if str(index) in to_remove:
                print("Removing relation %s" %rlt.toString())
                padmet._delRelation(rlt)
        padmet.generateFile(output)
    




if __name__ == "__main__":
    main()
