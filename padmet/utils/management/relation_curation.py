#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
"""

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
