# -*- coding: utf-8 -*-
"""
Description:
    From a list of genes, get from the linked reactions the list of products.

    R1 is linked to G1, R1 produces M1 and M2.  output: M1,M2. Takes into account reversibility

"""
    
def gene_to_targets(padmet, genes_file, output, verbose=False):
    """
    From a list of genes, get from the linked reactions the list of products.
    R1 is linked to G1, R1 produces M1 and M2.  output: M1,M2. Takes into account reversibility

    Parameters
    ----------
    padmet: padmet.classes.PadmetSpec
        padmet to explore
    genes_file: str
        path of genes file, 1 gene id by line
    output: str
        pathname of the output file
    verbose: bool
        if True print information
    """
    
    with open(genes_file,'r') as f:
        all_genes = f.read().splitlines()
        nb_genes = len(all_genes)
    all_targets = set()
    count = 0
    for gene_id in all_genes:
        count += 1
        try:
            #check if gene id in padmet
            padmet.dicOfNode[gene_id]
            if verbose: print("%s/%s %s" % (count, nb_genes, gene_id))
        except KeyError:
            if verbose: print("Gene %s not found in padmetSpec" %gene_id)
            continue

            #get all reactions linked to the gene
        try:
            reactions_linked = [padmet.dicOfNode[rlt.id_in] 
            for rlt in padmet.dicOfRelationOut[gene_id]
            if rlt.type == "is_linked_to"]
        except KeyError:
            if verbose: print("the gene %s is not linked to any reaction" %gene_id)
            continue
        if verbose: print("\t%s reactions linked" %len(reactions_linked))

        for reaction_node in reactions_linked:
            reaction_node_id = reaction_node.id
            if verbose: print("\t\t"+reaction_node_id)
            #all products
            products = set([rlt.id_out 
            for rlt in padmet.dicOfRelationIn.get(reaction_node_id, None)
            if rlt.type == "produces"])
            #add the products as targets
            all_targets = all_targets.union(products)
            #if reaction is reversible (reversible or UNKNOWN direction) add reactants
            if reaction_node.misc["DIRECTION"][0] != "LEFT-TO-RIGHT":
                #all reactants
                reactants = set([rlt.id_out 
                for rlt in padmet.dicOfRelationIn.get(reaction_node_id, None)
                if rlt.type == "consumes"])
                all_targets = all_targets.union(reactants)
    
    if len(all_targets) != 0:
        with open(output, 'w') as f:
            f.write("\n".join(all_targets))
        
                
