#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
This file is part of padmet.

padmet is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

padmet is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with padmet. If not, see <http://www.gnu.org/licenses/>.

@author: Meziane AITE, meziane.aite@inria.fr
Description:
Contains all necessary functions to generate wikiPages from a padmet file and update 
a wiki online. Require WikiManager module (with wikiMate,Vendor)
"""
from padmetRef import PadmetRef
from padmetSpec import PadmetSpec
from multiprocessing import Pool, cpu_count
import os
import shutil

global gene_template, metabolite_template, reaction_template, pathway_template
"""
Theses templates are used to create the pages for a wiki.
"""

gene_template = ['[[Category:Gene]]\n', 
'== Gene gene_id ==\n', 
'* Common name:\n', 
'* Synonyms:\n', 
'== Reactions Associated == \n', 
'== Encoded proteins == \n',
'== Other data == \n',
'== External links ==\n']

reaction_template = ['[[Category:Reaction]]\n', 
'== Reaction [DB_LINK_ACCESS=reaction_id reaction_id] == \n', 
'* Common name:\n', 
'* EC number:\n', 
'* Synonyms:\n', 
'== Reaction Formula ==\n', 
'\n', 
'reaction_formula_id\n', '\n', 
'reaction_formula_cname\n',
'== Pathways ==\n', 
'== Genes associated with this reaction  ==\n', 
'Genes have been associated with this reaction based on different elements listed below.\n',
'== Original source(s) ==\n',
'== External links ==\n']

metabolite_template = ['[[Category:Metabolite]]\n', 
'== Metabolite [DB_LINK_ACCESS=metabolite_id metabolite_id] ==\n', 
'* Common name:\n', 
'* Synonyms:\n', 
'== Reactions known to consume the compound: ==\n', 
'== Reactions known to produce the compound: ==\n', 
'== Reactions of unknown directionality: ==\n', 
'== External links ==\n', 
'\n']

pathway_template = ['[[Category:Pathway]]\n', 
'== Pathway [DB_LINK_ACCESS=pathway_id pathway_id] ==\n', 
'* Common name:\n', 
'* Synonyms:\n', 
'== Reactions found == {{#set:number of reactions found=0}}\n', 
'== Reaction NOT found == {{#set:number of reactions NOT found=0}}\n', 
'== External links ==\n']

source_template = ['{{#ask: [[Category:Reaction]] [[source::SRC]]\n',
'| ?name\n',
'| ?ec-number\n',
'| ?is in pathway\n',
'| ?comment\n',
'}}\n']

category_gene = ['{{#ask: [[Category:Gene]]\n',
'| ?associated to reaction\n',
'}}']

category_metabolite = ['{{#ask: [[Category:Metabolite]]\n',
'| ?name\n',
'| ?consumed by\n',
'| ?produced by\n',
'| ?produced or consumed by\n',
'}}']

category_pathway = ['{{#ask: [[Category:Pathway]]\n',
'| ?name\n',
'| ?number of reactions found\n',
'| ?number of reactions NOT found\n',
'}}']

category_reaction = ['{{#ask: [[Category:Reaction]]\n',
'| ?name\n',
'| ?ec-number\n',
'| ?source\n',
'| ?is in pathway\n',
'| ?compartment\n',
'}}']

main_page = ['=== genome-scale metabolic network ===\n']

sidebar = ['* navigation\n',
'** mainpage|mainpage-description\n',
'** Special:Ask|Semantic search\n',
'* Categories\n',
'** Category:Reaction|Reactions\n',
'** Category:Gene|Genes\n',
'** Category:Metabolite|Metabolites\n',
'** Category:Pathway|Pathways\n',
'* Origin of reactions\n']

def create_all_wikiPages(padmetSpec_file, output_dir, padmetRef_file = None, verbose = False):
    """
    Main function to generete wikiPages for a padmet. The padmetRef is used to add extra
    information from the full database (ex: the numbers of total reactions in a pathway)
    Steps:
    0/ In output_dir, create 4 folder: metabolites, reactions, genes, pathways, navigation.
    1/ Get all reactions, metabolites in reactions, pathways and genes
    2/ for each element, create a file from the specific template in wikiCode

    @param padmetRef_file: the pathname to the file padmet of reference
    @param padmetSpec_file: the pathname to the file padmet to convert to wiki
    @param output_dir: the pathname to folder where to store
    @param verbose: print informations
    @type padmetRef_file, padmetSpec_file, output_dir, verbose: str
    @return: _
    @rtype: None
    """
    global padmetSpec, padmetRef, mp_output_dir, DB_LINK_ACCESS, all_reactions_id, all_genes_id, all_pathways_id, all_metabolites_id 
    
    if not output_dir.endswith("/"): output_dir += "/"
    mp_output_dir = output_dir
    padmetSpec = PadmetSpec(padmetSpec_file)
    #check if bigg or metacyc
    try:
        db = padmetSpec.info["DB_info"]["DB"].upper()
        if db == "METACYC":
            DB_LINK_ACCESS = "http://metacyc.org/META/NEW-IMAGE?object"
        elif db == "BIGG":
            DB_LINK_ACCESS = "http://bigg.ucsd.edu/search?query"
        else:
            DB_LINK_ACCESS = "UNKNOWN"
    except KeyError:
        DB_LINK_ACCESS = "UNKNOWN"
        
    if padmetRef_file is not None:
        padmetRef = PadmetRef(padmetRef_file)
    else:
        padmetRef = None

    # create all the directory required
    if verbose: print("Creating directory...\n")
    createDirectory(mp_output_dir, verbose)

    # def the reactions, metabolites, pathways and genes
    if verbose: print("Loading all entry ...\n")

    all_genes_id = [node.id for node in padmetSpec.dicOfNode.itervalues() 
    if node.type == "gene"]
    nbGenes = len(all_genes_id)
    if verbose: print("%s genes's wikiPages to create" %nbGenes)

    #get all reactions    
    all_reactions_id = [node.id for node in padmetSpec.dicOfNode.itervalues() 
    if node.type == "reaction"]
    nbReactions = len(all_reactions_id)
    if verbose: print("%s reactions's wikiPages to create" %nbReactions)

        #for each reactions, recovering consumed and produced compounds
    all_metabolites_id = set([rlt.id_out for rlt in padmetSpec.getAllRelation()
    if rlt.type in ["consumes","produces"]])
    #bug with E- need to pop
    all_metabolites_id = list(all_metabolites_id)
    try:
        all_metabolites_id.remove("E-")
    except ValueError:
        pass
    nbMetabolites = len(all_metabolites_id)
    if verbose: print("%s Metabolites wikiPages to create" %nbMetabolites)
    
    #get all pathways    
    all_pathways_id = set([rlt.id_out for rlt in padmetSpec.getAllRelation() 
    if rlt.type == "is_in_pathway" and padmetSpec.dicOfNode[rlt.id_in].type == "reaction"])
    nbPathways = len(all_pathways_id)
    if verbose: print("%s Pathways wikiPages to create" %nbPathways)

    #create all wiki pages
    if verbose: print("Starting multiprocessing wikiPage Generation...")
    p = Pool(processes=cpu_count())
    for entryGenerator, entryFunction in [(all_genes_id, mp_createWikiPageGene), (all_reactions_id, mp_createWikiPageReaction),
   (all_pathways_id, mp_createWikiPagePathway), (all_metabolites_id, mp_createWikiPageMetabolite)]:
       resultats = p.map(entryFunction, entryGenerator)
       [None for _ in resultats]

    createDefaultPage()
    if verbose: print("wikiPages successfully created\n")

def createDirectory(dirPath, verbose = False):
    """
    create the folders genes, reactions, metabolites, pathways in the folder dirPath/
    if already exist, it will replace old folders (and delete old files)
    """
    #simple check that dirPath is a dir:
    dirNames = ["genes","reactions","metabolites","pathways","navigation"]
    #creatings the directory which will contains the wiki pages
    for d in dirNames:
        if not os.path.exists(dirPath+d):
            if verbose: print("Creating directory: "+dirPath+d)
            os.makedirs(dirPath+d)
        else:
            if verbose: print("The directory "+dirPath+d+" already exist. Old pages will be deleted")
            shutil.rmtree(dirPath+d)
            os.makedirs(dirPath+d)

def mp_createWikiPageGene(gene_id):
    """
    multiProcessing version of the function.
    @param gene_id: the id of the gene to create a wiki page
    @type gene_id: str
    @return: _
    @rtype: None
    """
    pageInArray = createWikiPageGene(padmetSpec, gene_template, gene_id)
    gene_id = gene_id.replace("/",".")
    fileName = mp_output_dir+"genes/"+gene_id+".txt"
    createWikiFile(pageInArray, fileName)

def createWikiPageGene(padmetSpec, gene_template, gene_id):
    """
    create a file with all the wikicode to create the page of the given Gene
    @param padmetSpec: the Padmet instance of the network
    @param gene_template: the template gene page
    @param gene_id: gene id
    @type padmetSpec: PadmetSpec
    @type gene_template: list
    @type gene_id: str
    @return: pageInArray corresponding to the wikiPage
    @rtype: list
    """
    #try te retrieve the gene node for gene_id     
    try:
        gene_node = padmetSpec.dicOfNode[gene_id]
        nodeClass = gene_node.type
        if nodeClass != "gene":
            raise TypeError("The given arguments doesn't refer to a gene: "+gene_id)
            exit()
    except KeyError:
        print("Unable to find the node with given arguments: "+gene_id)
        exit()

    pageInArray = list(gene_template)
    #2nd line is where the gene name is defined
    pageInArray[1] = pageInArray[1].replace("gene_id",gene_id)
    pageInArray[1] = pageInArray[1].replace("DB_LINK_ACCESS",DB_LINK_ACCESS)


    if "COMMON_NAME" in gene_node.misc.keys():
        cname_index = pageInArray.index("* Common name:\n")+1
        cnames = ["** "+cname+"{{#set:name="+cname+"}}\n" 
        for cname in gene_node.misc["COMMON_NAME"]]
        pageInArray[cname_index:cname_index] = cnames

    #check synonymous ('has_name' relations) / generator of synonymous
    try:
        synonymes = [padmetSpec.dicOfNode[rlt.id_out].misc["LABEL"][0] 
        for rlt in padmetSpec.dicOfRelationIn[gene_id]
        if rlt.type == "has_name"]
        if len(synonymes) != 0:
            syns_index = pageInArray.index("* Synonyms:\n")+1
            syns = [ "** "+syn+"{{#set:name="+syn+"}}\n" for syn in synonymes]
            pageInArray[syns_index:syns_index] = syns
    except KeyError:
        pass

    #Recovering all associations
    rxn_index = pageInArray.index("== Reactions Associated == \n") + 1
    try:
        #k=gene_id, v = list of ASSIGNMENT
        assoc_rlts = [(rlt.id_in, rlt.misc.get("ASSIGNMENT",["UNKNOWN"])) for rlt in padmetSpec.dicOfRelationOut[gene_id]
        if rlt.type == "is_linked_to"]
        dict_assign_gene = {}
        for rxn_id, all_assignments in assoc_rlts:
            for assign in all_assignments:
                try:
                    dict_assign_gene[assign].append(rxn_id)
                except KeyError:
                    dict_assign_gene[assign] = [rxn_id]
            
        #Create a list of string relative to reaction association
        assoc_gene_rxn = []
        for assignment, list_of_rxn in dict_assign_gene.iteritems():
            type_of_assoc = "* "+assignment+"\n"
            type_of_assoc = "* '''"+assignment+"''' {{#set:evidence="+assignment+"}}\n"
            assoc_gene_rxn.append(type_of_assoc)
            assoc_gene_rxn += ["** [["+rxn_id+"]]{{#set:associated to reaction="+rxn_id+"}}\n" for rxn_id in list_of_rxn]
    except KeyError:
        assoc_gene_rxn = ["* NONE\n"]
    #Insert the latter in the pageInArray
    pageInArray[rxn_index:rxn_index] = assoc_gene_rxn

    #Recover all encoded proteins
    prot_index = pageInArray.index("== Encoded proteins == \n") + 1
    try:
        prot_encoded = ["* "+rlt.id_out for rlt in padmetSpec.dicOfRelationIn[gene_id]
         if rlt.type == "codes_for"]
        prot_encoded = "\n".join(prot_encoded)
    except KeyError:
        prot_encoded = ["* NONE\n"]
    #Insert the latter in the pageInArray
    pageInArray[prot_index:prot_index] = prot_encoded

    if "TF" in gene_node.misc.keys():
        tf_index = pageInArray.index("== Transcription factors == \n")+1
        tfs = ["* "+tf+"{{#set:TF="+tf+"}}\n" 
        for tf in gene_node.misc["TF"]]
        pageInArray[tf_index:tf_index] = tfs

    #External links
    try:
        xrefs = []
        xref_index = pageInArray.index("== External links ==\n")+1
        xrefs_nodes = [padmetSpec.dicOfNode[rlt.id_out] for rlt in padmetSpec.dicOfRelationIn[gene_id] if rlt.type == "has_xref"]
        for node in xrefs_nodes:
            xrefs.append("* "+node.misc["DB"][0]+": "+node.misc["ID"][0]+"{{#set:xref="+node.misc["ID"][0]+"}}\n")
            for k,v in node.misc.iteritems():
                if k not in ["DB","ID"]:
                    for i in v:
                        xrefs.append("** "+str(k)+": "+str(i)+"\n")
                   
        pageInArray[xref_index:xref_index] = xrefs
    except KeyError:
        pass

    return pageInArray

    
def mp_createWikiPageReaction(reaction_id):
    """
    multiProcessing version of the function.
    @param reactop,_id: the id of the reaction to create a wiki page
    @type reaction_id: str
    @return: _
    @rtype: None
    """
    pageInArray = createWikiPageReaction(padmetSpec, reaction_template, reaction_id, padmetRef)
    reaction_id = reaction_id.replace("/",".")
    fileName = mp_output_dir+"reactions/"+reaction_id+".txt"
    createWikiFile(pageInArray,fileName)

def createWikiPageReaction(padmetSpec, reaction_template, reaction_id, padmetRef):
    """
    create a file with all the wikicode to create the page of the given Reaction
    @param padmetSpec: the Padmet instance of the network
    @param reaction_template: the template reaction page
    @param reaction_id: reaction id
    @param padmetRef: the Padmet instance of the reference database
    @type padmetSpec: PadmetSpec
    @type reaction_template: list
    @type reaction_id: str
    @type padmetRef: PadmetRef
    @return: pageInArray corresponding to the wikiPage
    @rtype: list
    """
    #try te retrieve the gene node for gene_id     
    try:
        reaction_node = padmetSpec.dicOfNode[reaction_id]
        nodeClass = reaction_node.type
        if nodeClass != "reaction":
            raise TypeError("The given arguments doesn't refer to a reaction: "+reaction_id)
            exit()
    except KeyError:
        print("Unable to find the node with given arguments: "+reaction_id)
        exit()
    
    pageInArray = list(reaction_template)
    #2nd line is where the reaction name is defined
    pageInArray[1] = pageInArray[1].replace("reaction_id", reaction_id)
    pageInArray[1] = pageInArray[1].replace("DB_LINK_ACCESS",DB_LINK_ACCESS)


    if "COMMON_NAME" in reaction_node.misc.keys():
        cname_index = pageInArray.index("* Common name:\n")+1
        cnames = ["** "+cname+"{{#set:name="+cname+"}}\n" 
        for cname in reaction_node.misc["COMMON_NAME"]]
        pageInArray[cname_index:cname_index] = cnames
    #check if ec number associated
    if "EC_NUMBER" in reaction_node.misc.keys():
        ec = reaction_node.misc["EC_NUMBER"][0]
        ec = ec.replace("EC-","")
        i = pageInArray.index("* EC number:\n")
        pageInArray[i] = "EC number: [http://enzyme.expasy.org/EC/"+ec+" "+ec+"]{{#set:ec-number="+ec+"}}\n"

    #check synonymous ('has_name' relations) / generator of synonymous
    try:
        synonymes = [padmetSpec.dicOfNode[rlt.id_out].misc["LABEL"][0] 
        for rlt in padmetSpec.dicOfRelationIn[reaction_id]
        if rlt.type == "has_name"]
        if len(synonymes) != 0:
            syns_index = pageInArray.index("* Synonyms:\n")+1
            syns = [ "** "+syn+"{{#set:name="+syn+"}}\n" for syn in synonymes]
            pageInArray[syns_index:syns_index] = syns
    except KeyError:
        pass
    """
    #check if hmm association
    i = pageInArray.index("* HMM model:\n")
    try:
        #recuperer tt les relations catalyses et si ya assignement qlq chose
        HMM = [rlt for rlt in padmetSpec.getRelation(nodeId,"out") if rlt.type == "catalyses"]
        pageInArray[i] = "* HMM model: Yes\n"
    except TypeError:
        pageInArray[i] = "* HMM model: No\n"
        HMM = None
    """

    # Recovering the formula
    direction = reaction_node.misc["DIRECTION"][0]
    if direction == "UNKNOWN":
        direction = " '''=>/<=>''' "
    elif direction == "REVERSIBLE":
        direction = " '''<=>''' "
    elif direction == "LEFT-TO-RIGHT":
        direction = " '''=>''' "
    
    #generator ['stoic cpd_id', ...]
    reactants = [rlt.misc["STOICHIOMETRY"][0]+" [["+rlt.id_out+"]]" 
    for rlt in padmetSpec.dicOfRelationIn.get(reaction_id, None) if rlt.type == "consumes"]
    products = [rlt.misc["STOICHIOMETRY"][0]+" [["+rlt.id_out+"]]"
    for rlt in padmetSpec.dicOfRelationIn.get(reaction_id, None) if rlt.type == "produces"]
    formula_id = " '''+''' ".join(reactants)+direction+" '''+''' ".join(products)+"\n"
        
    try:
        reactants = [rlt.misc["STOICHIOMETRY"][0]+" "+padmetSpec.dicOfNode[rlt.id_out].misc["COMMON_NAME"][0] 
        for rlt in padmetSpec.dicOfRelationIn.get(reaction_id, None) if rlt.type == "consumes"]
        products = [rlt.misc["STOICHIOMETRY"][0]+" "+padmetSpec.dicOfNode[rlt.id_out].misc["COMMON_NAME"][0]
        for rlt in padmetSpec.dicOfRelationIn.get(reaction_id, None) if rlt.type == "produces"]
        formula_cname = " '''+''' ".join(reactants)+direction+" '''+''' ".join(products)+"\n"
    except KeyError:
        formula_cname = "\n"
    
    i = pageInArray.index("reaction_formula_id\n")
    pageInArray[i] = pageInArray[i].replace("reaction_formula_id\n",formula_id)
    i = pageInArray.index("reaction_formula_cname\n")
    pageInArray[i] = pageInArray[i].replace("reaction_formula_cname\n",formula_cname)

    #recovering pathways associated
    try:
        ptw_index = pageInArray.index("== Pathways ==\n") + 1
        pathways = [padmetSpec.dicOfNode[rlt.id_out] for rlt in padmetSpec.dicOfRelationIn[reaction_id] 
        if rlt.type == "is_in_pathway"]
        for ptw in pathways:
            ptw_id = ptw.id
            #recovering the nb of reactions associated to the pathway
            if padmetRef is not None:
                try:
                    nbReactionsTotal = len([rlt for rlt in padmetRef.dicOfRelationOut[ptw_id] if rlt.type == "is_in_pathway"])
                # If keyError: pathway not in padmetRef, pathway added manualy
                except KeyError: 
                    nbReactionsTotal = "NA"
            else:
                nbReactionsTotal = "NA"
            nbReactionsFound = len([rlt for rlt in padmetSpec.dicOfRelationOut[ptw_id] if rlt.type == "is_in_pathway"])
            pageInArray.insert(ptw_index,"{{#set:is in pathway="+ptw_id+"}}\n")
            ptw_index += 1
            try:
                pageInArray.insert(ptw_index,"* [["+ptw_id+"]]: "+ptw.misc["COMMON_NAME"][0]+", [http://metacyc.org/META/NEW-IMAGE?object="+ptw_id+" view in MetaCyc]\n")
            except KeyError:
                pageInArray.insert(ptw_index,"* [["+ptw_id+"]]: "+" ["+DB_LINK_ACCESS+"="+ptw_id+" view in MetaCyc]\n")
            ptw_index += 1
            pageInArray.insert(ptw_index,"** '''"+str(nbReactionsFound)+"''' reactions found over '''"+str(nbReactionsTotal)+"''' reactions in the full pathway as defined in MetaCyc\n")
            ptw_index += 1
    except KeyError:
        pass
    #recovering genes associated 
    assoc_index = pageInArray.index("== Genes associated with this reaction  ==\n") + 2
    try:
        #k=gene_id, v = list of ASSIGNMENT
        assoc_rlts = [(rlt.id_out, rlt.misc.get("ASSIGNMENT",["UNKNOWN"])) for rlt in padmetSpec.dicOfRelationIn[reaction_id]
        if rlt.type == "is_linked_to"]
        dict_assign_gene = {}
        for gene_id, all_assignments in assoc_rlts:
            for assign in all_assignments:
                try:
                    dict_assign_gene[assign].append(gene_id)
                except KeyError:
                    dict_assign_gene[assign] = [gene_id]
            
        #Create a list of string relative to reaction association
        assoc_gene_rxn = []
        for assignment, list_of_genes in dict_assign_gene.iteritems():
            type_of_assoc = "* "+assignment+"\n"
            type_of_assoc = "* '''"+assignment+"''' {{#set:evidence="+assignment+"}}\n"
            assoc_gene_rxn.append(type_of_assoc)
            assoc_gene_rxn += ["** [["+gene_id+"]] {{#set:gene association="+gene_id+" ("+assignment+")}}\n" 
        for gene_id in list_of_genes]
    except KeyError:
        assoc_gene_rxn = ["* NONE\n"]

    #recovering orginal sources used to add this reaction
    if "SOURCE" in reaction_node.misc.keys():
        source_index = pageInArray.index('== Original source(s) ==\n') + 1
        try:
            comment = reaction_node.misc["COMMENT"][0]
        except KeyError:
            comment = None
        for src in reaction_node.misc["SOURCE"]:
            if comment:
                sources = ["* "+src+"{{#set:source="+src+"}}\n", 
                "** "+comment+"{{#set:comment="+comment+"}}\n"]
            else:
                sources = ["* "+src+"{{#set:source="+src+"}}\n"]
            pageInArray[source_index:source_index] = sources
    else:
        source_index = pageInArray.index('== Original source(s) ==\n') + 1
        sources = []
        for rlt in padmetSpec.dicOfRelationIn[reaction_id]:
            if rlt.type == "has_suppData":
                src = padmetSpec.dicOfNode[rlt.id_out].misc["ORIGIN_FILE"][0]
                sources.append("* "+src+"{{#set:source="+src+"}}\n") 
        pageInArray[source_index:source_index] = sources

        
    #Insert the latter in the pageInArray
    pageInArray[assoc_index:assoc_index] = assoc_gene_rxn
    """            
    # enzyme associated
    i = pageInArray.index("== MetaCyc Enzymes associated with this reaction ==\n")+1
    enzymes = (padmetRef.dicOfNode[rlt.id_in] for rlt in padmetRef.getRelation(nodeId,"out") if rlt.type == "catalyses")
    try:
        for e in enzymes:
            pageInArray.insert(i,"* "+e.misc["common_name"][0]+", [http://metacyc.org/META/NEW-IMAGE?object="+e.misc["metacyc_id"][0]+" "+e.misc["metacyc_id"][0]+"\n")
            i += 1
    except StopIteration:
        pass
    
    #HMM hits
    if HMM is not None:
        i = pageInArray.index("== Non validated associations (HMM hits) ==\n")+1
        for h in HMM:
            enzymeMetacycID = padmetSpec.dicOfNode[h.id_in].misc["metacyc_id"][0]
            pageInArray.insert(i,"* "+enzymeMetacycID+", "+h.misc["HMM"]+" (["+enzymeMetacycID+"]))")
            i += 1
    """    

    #External links
    try:
        xref_index = pageInArray.index("== External links ==\n")+1
        xrefs_nodes = [padmetSpec.dicOfNode[rlt.id_out] for rlt in padmetSpec.dicOfRelationIn[reaction_id] if rlt.type == "has_xref"]
        xrefs = [xrefLink(xrefNode) for xrefNode in xrefs_nodes]        
        pageInArray[xref_index:xref_index] = xrefs
    except KeyError:
        pass
    return(pageInArray)

def mp_createWikiPagePathway(pathway_id):
    """
    multiProcessing version of the function.
    @param pathway_id: the id of the pathway to create a wiki page
    @type pathway_id: str
    @return: _
    @rtype: None
    """
    pageInArray = createWikiPagePathway(padmetSpec, pathway_template, pathway_id, padmetRef)
    pathway_id = pathway_id.replace("/",".")
    fileName = mp_output_dir+"pathways/"+pathway_id+".txt"
    createWikiFile(pageInArray,fileName)
    
def createWikiPagePathway(padmetSpec, pathway_template, pathway_id, padmetRef = None ):
    """
    create a file with all the wikicode to create the page of the given Reaction
    @param padmetSpec: the Padmet instance of the network
    @param pathway_template: the template pathway page
    @param pathway_id: pathway id
    @type padmetSpec: PadmetSpec
    @type pathway_template: list
    @type pathway_id: str
    @return: pageInArray corresponding to the wikiPage
    @rtype: list
    """
    #try te retrieve the pathway node for pathway_id     
    try:
        pathway_node = padmetSpec.dicOfNode[pathway_id]
        nodeClass = pathway_node.type
        if nodeClass != "pathway":
            raise TypeError("The given arguments doesn't refer to a pathway: "+pathway_id)
            exit()
    except KeyError:
        print("Unable to find the node with given arguments: "+pathway_id)
        exit()
        
    pageInArray = list(pathway_template)
    #2nd line is where the pathway name is defined
    pageInArray[1] = pageInArray[1].replace("pathway_id", pathway_id)
    pageInArray[1] = pageInArray[1].replace("DB_LINK_ACCESS",DB_LINK_ACCESS)

    
    if "COMMON_NAME" in pathway_node.misc.keys():
        cname_index = pageInArray.index("* Common name:\n")+1
        cnames = ["** "+cname+"{{#set:name="+cname+"}}\n" 
        for cname in pathway_node.misc["COMMON_NAME"]]
        pageInArray[cname_index:cname_index] = cnames
    
    #check synonymous ('has_name' relations) / generator of synonymous
    try:
        synonymes = [padmetSpec.dicOfNode[rlt.id_out].misc["LABEL"][0] 
        for rlt in padmetSpec.dicOfRelationIn[pathway_id]
        if rlt.type == "has_name"]
        if len(synonymes) != 0:
            syns_index = pageInArray.index("* Synonyms:\n")+1
            syns = [ "** "+syn+"{{#set:name="+syn+"}}\n" for syn in synonymes]
            pageInArray[syns_index:syns_index] = syns
    except KeyError:
        pass
    #check nb reactions found in padmetRef and in current network
    rxn_found_index = pageInArray.index("== Reactions found == {{#set:number of reactions found=0}}\n")
    reactionsFound = [padmetSpec.dicOfNode[rlt.id_in] for rlt in padmetSpec.dicOfRelationOut[pathway_id] if rlt.type == "is_in_pathway" 
    and padmetSpec.dicOfNode[rlt.id_in].type == "reaction"]
    pageInArray[rxn_found_index] = pageInArray[rxn_found_index].replace("0",str(len(reactionsFound)))
    #r_id: r_name view in Metacyc
    rxn_found_index += 1
    for r_node in reactionsFound:
        r_id = r_node.id
        try:
            pageInArray.insert(rxn_found_index,"* [["+r_id+"]]: "+r_node.misc["COMMON_NAME"][0]+", [http://metacyc.org/META/NEW-IMAGE?object="+r_id+" view in MetaCyc]\n")
        except KeyError:
            pageInArray.insert(rxn_found_index,"* [["+r_id+"]]: "+" ["+DB_LINK_ACCESS+"="+r_id+" view in MetaCyc]\n")
        rxn_found_index += 1

    rxn_not_found_index = pageInArray.index("== Reaction NOT found == {{#set:number of reactions NOT found=0}}\n")
    if padmetRef is not None:
        reactionsTotal = [padmetRef.dicOfNode[rlt.id_in] for rlt in padmetRef.dicOfRelationOut[pathway_id] if rlt.type == "is_in_pathway"
        and padmetRef.dicOfNode[rlt.id_in].type == "reaction"]
        reactionsIdsNotFound = set([node.id for node in reactionsTotal]).difference(set([node.id for node in reactionsFound]))
        if len(reactionsIdsNotFound) != 0:
            pageInArray[rxn_not_found_index] = pageInArray[rxn_not_found_index].replace("0",str(len(reactionsIdsNotFound)))
            reactionsNotFound = [padmetRef.dicOfNode[rxn_id] for rxn_id in reactionsIdsNotFound] 
            rxn_not_found_index += 1
            for r_node in reactionsNotFound:
                r_id = r_node.id
                try:
                    pageInArray.insert(rxn_not_found_index,"* [["+r_id+"]]: "+r_node.misc["COMMON_NAME"][0]+", [http://metacyc.org/META/NEW-IMAGE?object="+r_id+" view in MetaCyc]\n")
                except KeyError:
                    pageInArray.insert(rxn_not_found_index,"* [["+r_id+"]]: "+" ["+DB_LINK_ACCESS+"="+r_id+" view in MetaCyc]\n")
                rxn_not_found_index += 1
    else:
        pageInArray[rxn_not_found_index] = pageInArray[rxn_not_found_index].replace("0","NA")

    #External links
    try:
        xref_index = pageInArray.index("== External links ==\n")+1
        xrefs_nodes = [padmetSpec.dicOfNode[rlt.id_out] for rlt in padmetSpec.dicOfRelationIn[pathway_id] if rlt.type == "has_xref"]
        xrefs = [xrefLink(xrefNode) for xrefNode in xrefs_nodes]        
        pageInArray[xref_index:xref_index] = xrefs
    except KeyError:
        pass
    
    return(pageInArray)

def mp_createWikiPageMetabolite(metabolite_id):
    """
    multiProcessing version of the function.
    @param metabolite_id: the id of the metabolite to create a wiki page
    @type metabolite_id: str
    @return: _
    @rtype: None
    """
    pageInArray = createWikiPageMetabolite(padmetSpec, metabolite_template, metabolite_id)
    metabolite_id = metabolite_id.replace("/",".")    
    fileName = mp_output_dir+"metabolites/"+metabolite_id+".txt"
    createWikiFile(pageInArray, fileName)

def createWikiPageMetabolite(padmetSpec, metabolite_template, metabolite_id):
    """
    create a file with all the wikicode to create the page of the given metabolite
    @param padmetSpec: the Padmet instance of the network
    @param metabolite_template: the template metabolite page
    @param metabolite_id: metabolite id
    @type padmetSpec: PadmetSpec
    @type metabolite_template: list
    @type metabolite_id: str
    @return: pageInArray corresponding to the wikiPage
    @rtype: list
    """
    #try te retrieve the pathway node for metabolite_id     
    try:
        metabolite_node = padmetSpec.dicOfNode[metabolite_id]
        nodeClass = metabolite_node.type
    except KeyError:
        print("Unable to find the node with given arguments: "+metabolite_id)
        exit()
        
    pageInArray = list(metabolite_template)
    #2nd line is where the metabolite name is defined
    pageInArray[1] = pageInArray[1].replace("metabolite_id", metabolite_id)
    pageInArray[1] = pageInArray[1].replace("DB_LINK_ACCESS",DB_LINK_ACCESS)

    
    if "COMMON_NAME" in metabolite_node.misc.keys():
        cname_index = pageInArray.index("* Common name:\n")+1
        cnames = ["** "+cname+"{{#set:name="+cname+"}}\n" 
        for cname in metabolite_node.misc["COMMON_NAME"]]
        pageInArray[cname_index:cname_index] = cnames
    
    #check synonymous ('has_name' relations) / generator of synonymous
    try:
        synonyms = [padmetSpec.dicOfNode[rlt.id_out].misc["LABEL"][0] 
        for rlt in padmetSpec.dicOfRelationIn[metabolite_id]
        if rlt.type == "has_name"]
        if len(synonyms) != 0:
            syns_index = pageInArray.index("* Synonyms:\n")+1
            syns = [ "** "+syn+"{{#set:name="+syn+"}}\n" for syn in synonyms]
            pageInArray[syns_index:syns_index] = syns    
    except KeyError:
        pass
    #reactions that consume or produce the compound
    produce = (padmetSpec.dicOfNode[rlt.id_in] for rlt in padmetSpec.dicOfRelationOut[metabolite_id] if rlt.type == "produces")
    consume = (padmetSpec.dicOfNode[rlt.id_in] for rlt in padmetSpec.dicOfRelationOut[metabolite_id] if rlt.type == "consumes")
    unknown = set()
    i = pageInArray.index("== Reactions known to consume the compound: ==\n")+1
    try:
        for r_node in consume:
            r_id = r_node.id
            try:
                if r_node.misc["DIRECTION"][0] != "REVERSIBLE":
                    pageInArray.insert(i,"* [["+r_id+"]]{{#set:consumed by="+r_id+"}}\n")
                    i += 1
                else:
                    unknown.add(r_id)
            except KeyError:
                pass
    #0 reactions that consume
    except StopIteration:
        pass
    i = pageInArray.index("== Reactions known to produce the compound: ==\n")+1
    try:
        for r_node in produce:
            r_id = r_node.id
            try:
                if r_node.misc["DIRECTION"][0] != "REVERSIBLE":
                    pageInArray.insert(i,"* [["+r_id+"]]{{#set:produced by="+r_id+"}}\n")
                    i += 1
                else:
                    unknown.add(r_id)
            except KeyError:
                pass
    #0 reactions that produce
    except StopIteration:
        pass
    i = pageInArray.index("== Reactions of unknown directionality: ==\n")+1
    for r_id in unknown:
        pageInArray.insert(i,"* [["+r_id+"]]{{#set:produced or consumed by="+r_id+"}}\n")
        i += 1

    #External links
    try:
        xref_index = pageInArray.index("== External links ==\n")+1
        xrefs_nodes = [padmetSpec.dicOfNode[rlt.id_out] for rlt in padmetSpec.dicOfRelationIn[metabolite_id] if rlt.type == "has_xref"]
        #xrefs = ["* "+node.misc["DB"][0]+": "+node.misc["ID"][0]+"\n" for node in xrefs_nodes]
        xrefs = [xrefLink(xrefNode) for xrefNode in xrefs_nodes]        
        pageInArray[xref_index:xref_index] = xrefs
    except KeyError:
        pass
    
    return(pageInArray)

def assoc_correspondance(type_of_assoc):
    """
    The association gene-reaction can be convert to something more meaningful than a
    simple tag. the dictionnary dictOfAssoc make the link for some of the associations.
    @param type_of_assoc: the evidence allowing to link a gene to a reaction
    @type type_of_assoc: str
    @return: the meaningful correspondence
    @rtype: str
    """
    dictOfAssoc = {"EC-NUMBER":"EC Number (pathway Tools)", 
    "AUTOMATED-NAME-MATCH":"Functional description (Pathway Tools)",
    "GAP-FILLING":"Gap-Filling (Meneco)",
    "MANUAL":"MANUAL"}
    result = dictOfAssoc.get(type_of_assoc, type_of_assoc)    
    return result


def xrefLink(xrefNode):
    if xrefNode.misc["DB"][0] == "METACYC":
        toInsert = "* pubchem : [http://pubchem.ncbi.nlm.nih.gov/summary/summary.cgi?cid="+xrefNode.misc["ID"][0]+" "+xrefNode.misc["ID"][0]+"]\n"

    elif xrefNode.misc["DB"][0] == "UNIPROT":
        toInsert = "* uniprot : [http://www.uniprot.org/uniprot/"+xrefNode.misc["ID"][0]+" "+xrefNode.misc["ID"][0]+"] {{#set:uniprot id="+xrefNode.misc["ID"][0]+"}}\n"

    elif xrefNode.misc["DB"][0] in ["KEGG","LIGAND-RXN","LIGAND-CPD"]:
        toInsert = "* kegg : [http://www.genome.jp/dbget-bin/www_bget?"+xrefNode.misc["ID"][0]+" "+xrefNode.misc["ID"][0]+"] {{#set:kegg id="+xrefNode.misc["ID"][0]+"}}\n"

    elif xrefNode.misc["DB"][0] == "RHEA":
        toInsert = "* RHEA : [http://www.ebi.ac.uk/rhea/reaction.xhtml?id="+xrefNode.misc["ID"][0]+" "+xrefNode.misc["ID"][0]+"]\n"

    elif xrefNode.misc["DB"][0] == "WIKIPEDIA":
        toInsert = "* wikipedia : [http://en.wikipedia.org/wiki/"+xrefNode.misc["ID"][0]+" "+xrefNode.misc["ID"][0]+"]\n"

    elif xrefNode.misc["DB"][0] == "CHEBI":
        toInsert = "* chebi : [http://www.ebi.ac.uk/chebi/searchId.do?chebiId="+xrefNode.misc["ID"][0]+" "+xrefNode.misc["ID"][0]+"] {{#set:chebi id="+xrefNode.misc["ID"][0]+"}}\n"

    elif xrefNode.misc["DB"][0] == "PUBCHEM":
        toInsert = "* pubchem : [http://pubchem.ncbi.nlm.nih.gov/summary/summary.cgi?cid="+xrefNode.misc["ID"][0]+" "+xrefNode.misc["ID"][0]+"]\n"

    elif xrefNode.misc["DB"][0] == "ECOCYC":
        toInsert = "* ecocyc : [http://metacyc.org/ECOLI/NEW-IMAGE?object="+xrefNode.misc["ID"][0]+" "+xrefNode.misc["ID"][0]+"]\n"

    elif xrefNode.misc["DB"][0] == "CHEMSPIDER":
        toInsert = "* CHEMSPIDER : [http://www.chemspider.com/Chemical-Structure."+xrefNode.misc["ID"][0]+".html "+xrefNode.misc["ID"][0]+"]\n"

    elif xrefNode.misc["DB"][0] == "umbbd-compounds":
        toInsert = "* umbbd-compounds : [http://umbbd.ethz.ch/servlets/pageservlet?ptype=c&compID="+xrefNode.misc["ID"][0]+" "+xrefNode.misc["ID"][0]+"]\n"

    elif xrefNode.misc["DB"][0] == "ARACYC":
        toInsert = "* ARACYC : [http://metacyc.org/ARA/NEW-IMAGE?object="+xrefNode.misc["ID"][0]+" "+xrefNode.misc["ID"][0]+"]\n"

    elif xrefNode.misc["DB"][0] == "PIR":
        toInsert = "* PIR : [http://pir.georgetown.edu/cgi-bin/nbrfget?uid="+xrefNode.misc["ID"][0]+" "+xrefNode.misc["ID"][0]+"]\n"

    elif xrefNode.misc["DB"][0] == "NCI":
        toInsert = "* NCI : [http://cactus.nci.nih.gov/ncidb2.2/?nsc="+xrefNode.misc["ID"][0]+" "+xrefNode.misc["ID"][0]+"]\n"

    elif xrefNode.misc["DB"][0] == "knapsacK":
        toInsert = "* knapsacK : [http://kanaya.naist.jp/knapsack_jsp/information.jsp?word="+xrefNode.misc["ID"][0]+" "+xrefNode.misc["ID"][0]+"]\n"

    else:
        toInsert = "* "+xrefNode.misc["DB"][0]+" : "+xrefNode.misc["ID"][0]+"\n"
                
    return(toInsert)                

def createDefaultPage():
    """
    Default pages are the pages used in a wiki to navigate: the category, the sidebar, the mainpage.
    Sidebar: in the sidebar allows to navigate between categories of nodes and origins of reactions.
    The origin of a reaction is recovered from node.misc['source']. if no source from the suppData
    associated to this reaction suppdata_node.misc['origin_file]
    """
    #all source represent: all rxn.misc["source"] and the value of origin_file for suppData of all rxn
    all_sources = set()
    for rxn_id in all_reactions_id:
        rxn_node = padmetSpec.dicOfNode[rxn_id]
        if "SOURCE" in rxn_node.misc.keys():
            for src in rxn_node.misc["SOURCE"]:
                all_sources.add(src)
        else:
            for rlt in padmetSpec.dicOfRelationIn[rxn_id]:
                if rlt.type == "has_suppData":
                    src = os.path.splitext(padmetSpec.dicOfNode[rlt.id_out].misc["ORIGIN_FILE"][0])[0]
                    all_sources.add(src)
    for source in all_sources:
        pageInArray = list(source_template)
        pageInArray[0] = pageInArray[0].replace("SRC",source)
        fileName = mp_output_dir+"navigation/"+source+".txt"
        createWikiFile(pageInArray, fileName)
    #Category pages
    createWikiFile(category_gene, mp_output_dir+"navigation/"+"Category:Gene.txt")
    createWikiFile(category_metabolite, mp_output_dir+"navigation/"+"Category:Metabolite.txt")
    createWikiFile(category_pathway, mp_output_dir+"navigation/"+"Category:Pathway.txt")
    createWikiFile(category_reaction, mp_output_dir+"navigation/"+"Category:Reaction.txt")
    #Navigation pages
    for src in all_sources:
        sidebar.append("** "+src+"|"+src+"\n")
    createWikiFile(sidebar, mp_output_dir+"navigation/"+"MediaWiki:Sidebar.txt")
    
def createWikiFile(pageInArray,fileName):
    with open(fileName,'w') as f:
        for line in pageInArray:
            f.write(line)


