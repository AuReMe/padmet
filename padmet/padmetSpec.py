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
PadmetSpec is an object representing the metabolic network of a species(organism)
based on a reference database PadmetRef.
"""
from policy import Policy
from node import Node
from relation import Relation
from padmetRef import PadmetRef
import sbmlPlugin
try:
    from libsbml import *
except:
    print("package libsbml needed, use this cmd:\n pip install "
          + "python-libsbml")
    exit()


class PadmetSpec:
    """
    PadmetSpec is an object representing the metabolic network of a species(organism)
    based on a reference database PadmetRef.
    contains <Policy>, <Node> and <Relation>
    The policy defines the way Node and Relation are associated
    A node is an Object that contains information about an element of the network 
    (can be a pathway, reaction...).
    A realtion defines how two nodes are connected. In a relation there is 
    a node "in" and a node "out". (reactionX'in' consumes metaboliteX'out')
    PadmetSpec contains 3 attributs: 
        dicOfNode: a dictionary of node: key=Node's unique id / value = <Node>
        dicOfRelationIn: a dictionnary of relation with: key= nodeIN id / value = list of <relation>
        dicOfRelationOut: a dictionnary of relation with: key= nodeOut id / value = list of <relation>
        policy: a <policy>
        info: a dictionnary of informations about the network, the database used...
        This dictionnary is always represented in the header of a padmet file
    """
    def __init__(self, padmetSpec_file = None):
        """
        if None, initializes an empty <PadmetSpec>
        @param padmetSpec_file: pathname of the padmet file
        @type padmetSpec_file: str
        @return: _
        @rtype: None
        """
        if padmetSpec_file is not None:        
            self.loadGraph(padmetSpec_file)
        else:
            self.dicOfRelationIn = {}
            self.dicOfRelationOut = {}
            self.dicOfNode = {}
            self.policy = Policy()
            self.info = {}         
#==============================================================================
# Constructor / getter            
#==============================================================================

    def setInfo(self,source):
        """
        All the information printed in the header of a padmet stocked in a dict.
        {"metacyc":{version:XX,...},"ecocyc":{...}...}
        set Info from a dictionnary or copying from an other padmet
        @param source: may be a dict or an other padmet from where will be copied the info
        @tyme source: dict or Padmet
        """
        if type(source) is dict:
            self.info = source
        else: 
            self.info = source.info

    def setPolicy(self,source):
        """
        Set policy from a list or copying from an other padmet
        @param source: may be a list or an other padmet from where will be copied the policy
        @type source: list or Padmet
        """
        if type(source) is list:
            self.policy = Policy(source)
        else:
            self.policy = source.policy
        

    def setDicOfNode(self,source):
        """
        Set dicOfNode from a dict or copying from an other padmet
        @param source: may be a dict or an other padmet from where will be copied the dicOfNode
        @type source: dict or Padmet
        """
        if type(source) is dict:
            self.dicOfNode = source
        else:
            self.dicOfNode = source.dicOfNode


    def setdicOfRelationIn(self,source):
        """
        Set dicOfRelationIn from a dict or copying from an other padmet
        @param source: may be a dict or an other padmet from where will be copied the dicOfRelationIn
        @type source: dict or Padmet
        @rtype: None
        """
        if type(source) is dict:
            self.dicOfRelationIn = source
        else:
            self.dicOfRelationIn = source.dicOfRelationIn


    def setdicOfRelationOut(self,source):
        """
        Set dicOfRelationOut from a dict or copying from an other padmet
        @param source: may be a dict or an other padmet from where will be copied the dicOfRelationIn
        @type source: dict or Padmet
        @rtype: None
        """
        if type(source) is dict:
            self.dicOfRelationOut = source
        else:
            self.dicOfRelationOut = source.dicOfRelationOut

    
    def getAllRelation(self):
        """
        return a set of all relations
        @rtype: set
        """
        all_relation = set()
        for list_rlt in self.dicOfRelationIn.itervalues():
            for rlt in list_rlt:
                all_relation.add(rlt)
        for list_rlt in self.dicOfRelationIn.itervalues():
            for rlt in list_rlt:
                all_relation.add(rlt)
        return all_relation


    def loadGraph(self, padmet_file):
        """
        Allow to recover all the informations of the padmet file.
        The section Data Base informations corresponds to the self.info
        The section Nodes corresponds to the data of each nodes in self.dicOfNode, sep ="\t"
        The section Relations corresponds to the data of each relations in self.dicOfRelationIn/Out, sep ="\t"
        @param padmet_file: the pathname of the padmet file to load.
        @return: _
        @rtype: None
        """
        with open(padmet_file, 'r') as f:        
            padmet_in_array = [line for line in f.read().splitlines() if len(line) != 0]

        self.policy = Policy()        
        self.info = {}
        self.dicOfNode = {}
        self.dicOfRelationIn = {}
        self.dicOfRelationOut = {}

        try:
            info_index = padmet_in_array.index("Data Base informations")
        except ValueError:
            info_index = None
        #Index where start policy
        policy_index = padmet_in_array.index("Policy")
        #Index where start Nodes
        node_index = padmet_in_array.index("Nodes")
        #Index where start Relations       
        relation_index = padmet_in_array.index("Relations")
        
        #
        if info_index is not None:
            info_section = padmet_in_array[info_index+1:policy_index]
            #
            for line in info_section:
                if "\t" not in line:
                    line = line.replace(":","")
                    self.info[line] = {}
                    current_data = line
                else:
                    line = line.replace("\t","").split(":")
                    self.info[current_data][line[0]] = line[1]

        #
        policy_in_array = [line.split("\t") for line in padmet_in_array[policy_index+1:node_index]]
        self.setPolicy(policy_in_array)


        # generator of data to input in Node()
        node_section = (line.split("\t") for line in padmet_in_array[node_index+1:relation_index])
        #
        for data in node_section:
            #Create a new Node object (cf node.py)
            node_type = data[0]
            node_id = data[1]
            node_misc = {}
            if len(data) > 2:
                i = 2
                #add the misc data in a dict
                while (i+1) < len(data): 
                    #in case diff information for the same key
                    try:
                        node_misc[data[i]].append(data[i+1])
                    except KeyError:
                        node_misc[data[i]] = [data[i+1]]
                    i += 2
            node = Node(node_type, node_id, node_misc)
            #add the node into the dictionnay
            self.dicOfNode[node.id] = node 

        # generator of data to input in Relation()
        relation_section = (line.split("\t") for line in padmet_in_array[relation_index+1:])
        #instantiate a new Relation object (cf relation.py)
        for data in relation_section:
            rlt_id_in = data[0]
            rlt_type = data[1]
            rlt_id_out = data[2]
            rlt_misc = {}
            if len(data) > 3:
                i = 3
                while (i+1) < len(data): 
                    #in case diff information for the same key               
                    try:
                        rlt_misc[data[i]].append(data[i+1])
                    except KeyError:
                        rlt_misc[data[i]] = [data[i+1]]
                    i+=2

            relation = Relation(rlt_id_in, rlt_type, rlt_id_out, rlt_misc)
            try:
                self.dicOfRelationIn[rlt_id_in].append(relation)
            except KeyError:
                self.dicOfRelationIn[rlt_id_in] = [relation]
            try:
                self.dicOfRelationOut[rlt_id_out].append(relation)
            except KeyError:
                self.dicOfRelationOut[rlt_id_out] = [relation]

    def updateFromSbml(self, sbml_file, padmetRef = None, assocIdOriginRef = None, verbose = False, force = False):
        """
        Copy data from sbml to padmet. If the id of the sbml are 
        different with the padmetRef, need to add a dictionnary of associations
        between the id in the sbml and the id of the database of reference.
        if no padmetRef, will create reactions and compounds based on the given sbml
        @param sbml_file: pathname of the sbml file
        @param padmetRef: <PadmetRef> padmet of reference.
        @param assocIdOriginRef (opt): pathname of the file containing the 
        association original_id ref_id, sep ="\t"
        @param verbose: print info
        @param force: if true, allow to copy reactions that are not in padmetRef by creating new reaction
        @type sbml_file, assocIdOriginRef: str
        @type padmetRef: PadmetRef
        @type verbose, force: bool
        @return: _
        @rtype: None
        """
        file_name = os.path.basename(sbml_file)
        #using libSbml to read sbml_file
        if verbose: print ("loading sbml file: %s" %sbml_file)
        reader = SBMLReader()
        document = reader.readSBML(sbml_file)
        for i in range(document.getNumErrors()):
            print (document.getError(i).getMessage())

        model = document.getModel()
        listOfSpecies = model.getListOfSpecies()
        listOfReactions = model.getListOfReactions()
        nbReactions = len(listOfReactions)
        nbSpecies = len(listOfSpecies)
        if verbose:
            print("nb species: %s" %nbSpecies)
            print("nb reactions: %s" %nbReactions) 
        dicOfAssocReactions = {}
        dicOfAssocSpecies = {}
        #reading assocIdOriginRef line by line, one line = "origin_id\tref_id\n"
        if assocIdOriginRef is not None:
            if verbose: print("Parsing %s" %assocIdOriginRef)
            with open(assocIdOriginRef,'r') as f: #open the file               
                all_assoc = (line.split("\t") for line in f.read().splitlines()
                if not line.startswith("#"))
                for id_origin, id_ref in all_assoc:
                    if id_ref in padmetRef.dicOfNode.keys():
                        nodeType = padmetRef.dicOfNode[id_ref].type
                        #in case of reaction
                        if nodeType == "reaction":
                            dicOfAssocReactions.update({id_origin:id_ref})
                        #in case of species
                        else:
                            dicOfAssocSpecies.update({id_origin:id_ref})
                    else:
                        print("id ref %s from dict not found in padmetRef" %id_ref)
        else:
            for species in listOfSpecies:
                dicOfAssocSpecies[species.getId()] = sbmlPlugin.convert_from_coded_id(species.getId())[0]
            for reaction in listOfReactions:
                dicOfAssocReactions[reaction.getId()] = sbmlPlugin.convert_from_coded_id(reaction.getId())[0]

        if verbose: print("start with reactions")
        #dictionnary used to stock association gene:reaction
        dicOfGeneReaction = {}
        count = 0
        for (idOrigin, idRef) in dicOfAssocReactions.iteritems():
            count += 1
            #recovere the reaction element from sbml
            reactionSBML = listOfReactions.getElementBySId(idOrigin)
            nextStep = False
            if verbose:
                print("%s/%s" % (count, len(dicOfAssocReactions.keys())))
                print("idOrigin: %s \t idRef: %s" % (idOrigin, idRef))
            try:
                #check if reaction already in the padmet
                self.dicOfNode[idRef]
                nextStep = True
            except KeyError:
                if padmetRef is not None and idRef in padmetRef.dicOfNode.keys():
                    #copy the node reaction and dependency relations/nodes cf copyNode()
                    if verbose: print("Try to copy node from padmetRef")
                    self.copyNode(padmetRef, idRef)
                    self.dicOfNode[idRef]
                    nextStep = True
                elif force or padmetRef is None:
                    #force will allow to create not found reactions
                    cname = reactionSBML.getName()
                    if verbose: print("creating new reaction: %s" %idRef)
                    if reactionSBML.getReversible():
                        reaction_dir = "REVERSIBLE"
                    else:
                        reaction_dir = "LEFT-TO-RIGHT"
                    if len(cname) != 0:
                        reaction_node = Node("reaction", idRef, {"DIRECTION":[reaction_dir], "COMMON_NAME":[cname]})
                    else:
                        reaction_node = Node("reaction", idRef, {"DIRECTION":[reaction_dir]})
                    self.dicOfNode[idRef] = reaction_node
                    reactants = reactionSBML.getListOfReactants()
                    for reactant in reactants: #convert ids
                        try:
                            reactant_id = dicOfAssocSpecies[reactant.getSpecies()]
                        except KeyError:
                            reactant_id = sbmlPlugin.convert_from_coded_id(reactant.getSpecies().getId())[0]
                        reactant_compart = model.getElementBySId(reactant.getSpecies()).getCompartment()
                        reactant_stoich = reactant.getStoichiometry()
                        if reactant_id not in self.dicOfNode.keys():
                            if padmetRef is not None and reactant_id in padmetRef.dicOfNode.keys():
                                self._copyNodeExtend(padmetRef, reactant_id)
                            else:
                                if verbose: print("creating new compound: %s" %reactant_id)
                                reactant_cname = listOfSpecies.getElementBySId(reactant.getSpecies()).getName()
                                if reactant_cname:
                                    reactant_node = Node("compound", reactant_id, {"COMMON_NAME":[reactant_cname]})
                                else:
                                    reactant_node = Node("compound", reactant_id)
                                self.dicOfNode[reactant_id] = reactant_node
                        if reactant_compart is None:
                            if verbose: print("%s has no compart, set to 'c'" %reactant)
                            reactant_compart = "c"
                        consumes_rlt = Relation(idRef,"consumes",reactant_id, {"STOICHIOMETRY":[reactant_stoich],"COMPARTMENT":[reactant_compart]})
                        self._addRelation(consumes_rlt)
                    products = reactionSBML.getListOfProducts()
                    for product in products: #convert ids
                        try:
                            product_id = dicOfAssocSpecies[product.getSpecies()]
                        except KeyError:
                            product_id = sbmlPlugin.convert_from_coded_id(product.getSpecies().getId())[0]
                        product_compart = model.getElementBySId(product.getSpecies()).getCompartment()
                        product_stoich = product.getStoichiometry()
                        if product_id not in self.dicOfNode.keys():
                            if padmetRef is not None and product_id in padmetRef.dicOfNode.keys():
                                self._copyNodeExtend(padmetRef, product_id)
                            else:
                                print("creating new compound: %s" %product_id)
                                product_cname = listOfSpecies.getElementBySId(product.getSpecies()).getName()
                                if product_cname:
                                    product_node = Node("compound", product_id, {"COMMON_NAME":[product_cname]})
                                else:
                                    product_node = Node("compound", product_id)
                                self.dicOfNode[product_id] = product_node
                        if product_compart is None:
                            if verbose: print("%s has no compart, set to 'c'" %product)
                            product_compart = "c"
                        produces_rlt = Relation(idRef, "produces", product_id, {"STOICHIOMETRY":[product_stoich],"COMPARTMENT":[product_compart]})
                        self._addRelation(produces_rlt)
                    nextStep = True
                else:
                        print("Error: idRef %s not found" %idRef)
            if nextStep:
                #Reaction was found in current network
                #update SOURCE in misc
                try:
                    if file_name not in self.dicOfNode[idRef].misc["SOURCE"]:
                        self.dicOfNode[idRef].misc["SOURCE"].append(file_name)
                except KeyError:
                    self.dicOfNode[idRef].misc["SOURCE"] = [file_name]
                #Extracting all data to create the supplementary data node
                #Using sbmlPlugin to recovere the formula from the sbml
                formula = sbmlPlugin.extractFormula(reactionSBML)
                #Using sbmlPlugin to recover the note section from the sbml
                notes = sbmlPlugin.parseNotes(reactionSBML)
                #data will be stored in a suppData node
                data = {"ORIGIN_FILE":[file_name],"ORIGIN_ID":[str(idOrigin)],
                "NAME":[reactionSBML.getName()],"REVERSIBLE":[str(reactionSBML.getReversible())],
                "FORMULA":[formula]}
                #add notes to data                
                data.update(notes)
                #create the node suppData and the relation has_suppData
                suppData = self.createNode("suppData",data,[[idRef,"has_suppData","_self"]])
                if verbose: print("Creating suppData %s" %suppData[0])
                                
                #parses gene_association and create gene node or update already existing genes
                if "GENE_ASSOCIATION" in notes.keys():
                    #Using sbmlPlugin to recover all genes associated to the reaction
                    listOfGenes = sbmlPlugin.parseGeneAssoc(notes["GENE_ASSOCIATION"][0])
                    if len(listOfGenes) != 0:
                        for gene in listOfGenes:
                            try:
                                dicOfGeneReaction[gene].add(idRef)
                            except KeyError:
                                dicOfGeneReaction[gene] = set([idRef])
        if verbose: print("Creates gene's nodes:")
        count = 0
        nbGenes = len(dicOfGeneReaction.keys())
        for gene_id, rxns_associated in dicOfGeneReaction.iteritems():
            count += 1
            if verbose: print("%s/%s %s" % (count, nbGenes, gene_id))
            try:
                #check if gene already in the padmet
                gene_node = self.dicOfNode[gene_id]
            except KeyError:
                gene_node = Node("gene",gene_id)
                self.dicOfNode[gene_id] = gene_node
            for rxn_id in rxns_associated:
                try:
                    linked_rlt = [rlt for rlt in self.dicOfRelationIn[rxn_id] if rlt.type == "is_linked_to"
                    and rlt.id_out == gene_id][0]
                    try:
                        linked_rlt.misc["ASSIGNMENT"].append(file_name)
                    except KeyError:
                        linked_rlt.misc["ASSIGNMENT"] = [file_name]
                except IndexError:
                    linked_rlt = Relation(rxn_id, "is_linked_to", gene_id,{"ASSIGNMENT":[file_name]})
                    self._addRelation(linked_rlt)
                        
        if verbose: print("Case species:")
        count = 0
        #in case of species
        for (idOrigin, idRef) in dicOfAssocSpecies.iteritems():
            count += 1
            speciesSBML = listOfSpecies.getElementBySId(idOrigin)
            if verbose:
                print("%s/%s" % (count, len(dicOfAssocSpecies.keys())))
                print("idOrigin: %s \t %s" % (idOrigin, idRef))
            try:
                #check if species is already in the padmet
                self.dicOfNode[idRef]
                #Extracting all data to create the supplementary data node
                if verbose: print("Creating suppData")
                notes = sbmlPlugin.parseNotes(speciesSBML)
                try:
                    #skipping INCHI because of length
                    notes.pop("INCHI")
                except KeyError:
                    pass
                data = {"ORIGIN_FILE":[file_name],"ORIGIN_ID":[str(idOrigin)],
                "NAME":[speciesSBML.getName()],"CHARGE":[str(speciesSBML.getCharge())],
                "COMPARTMENT":[speciesSBML.getCompartment()]}
                data.update(notes)
                self.createNode("suppData",data,[[idRef,"has_suppData","_self"]])
            except KeyError:
                pass

            
    def generateFile(self, output):
        """
        Allow to create a padmet file to stock all the data.
        @param output: pathname of the padmet file to create
        @return: _
        @rtype: None
        """
        # Order the dictionary of node by unique id and the tuple of relation
        # by the node in id.
        dicKeys = self.dicOfNode.keys()
        dicKeys.sort()
        orderedTuplOfRelation = tuple(sorted(self.getAllRelation(),
                                           key=lambda x: x.id_in, reverse=False))
        with open(output,'w') as f:
            if len(self.info) != 0:
                f.write("Data Base informations\n")
                f.write("\n")
                for k,data in self.info.iteritems():
                    f.write(k+":\n")
                    for k,v in data.iteritems():
                       f.write("\t"+k+":"+v+"\n")
                f.write("\n")
            # It writes the policy of the TGDBP file.
            f.write("Policy\n")
            f.write("\n")
            policyInArray = self.policy.getPolicyInArray()
            for line in policyInArray:
                line = "\t".join(line)
                line = line + "\n"
                f.write(line)
            f.write("\n")
            # It writes the nodes of the TGDBP file.
            f.write("Nodes\n")
            f.write("\n")
            for nodeId in dicKeys:
                node = self.dicOfNode[nodeId]
                line = node.toString()+"\n"
                f.write(line)
            f.write("\n")
            # It writes the relations of the TGDBP file.
            f.write("Relations\n")
            f.write("\n")
            for rlt in orderedTuplOfRelation:
                line = rlt.toString()+"\n"
                f.write(line)

    def extract_pathway(self, node_id, padmetRef_file, output, sbml = None): 
        """
        Allow to extract a pathway in a csv file. Need a padmetRef to check the total number
        of reactions in the pathway.
        Header = "Reactions (metacyc_id)", "Reactions (common_name)", "EC-Number",
                      "Formula (metacyc_id)", "Formula (common_name)", "Found in the network"
        @param node_id: id of the pathway
        @param padmetRef_ile: pathname of the padmet ref file
        @param output: pathname of the output to create
        @param sbml: if true, create a sbml file of this pathway
        @return: _
        @rtype: None
        """
        #Retriving the ID and the node of the pathway.
        try:
            node = self.dicOfNode[node_id]
        except KeyError:
            print("%s not found" %node_id)
            return
        node_type = node.type
        if node_type != "pathway":
            print("%s is not a pathway" %node_id)
            return
        
        padmetRef = PadmetRef(padmetRef_file)
        #Retriving all the reactions associated to the pathway from the padmetRef
        total_reactions = [padmetRef.dicOfNode[rlt.id_in] for rlt in padmetRef.dicOfRelationOut.get(node_id, None)
        if rlt.type == "is_in_pathway"
        and padmetRef.dicOfNode[rlt.id_in].type == "reaction"]
        #Retriving the reactions (ID) present in the network
        found_reactions = [self.dicOfNode[rlt.id_in].id for rlt in self.dicOfRelationOut.get(node_id, None)
        if rlt.type == "is_in_pathway"
        and self.dicOfNode[rlt.id_in].type == "reaction"]

        if len(found_reactions) == 0:
            print("0 reactions associated to this pathway in the network")
            return
        
        with open(output,'w') as f:
            #Define header
            header = ["Reactions (metacyc_id)", "Reactions (common_name)", "EC-Number",
                      "Formula (metacyc_id)", "Formula (common_name)", "Found in the network"]
            header = "\t".join(header)+"\n"
            f.write(header)
            for rNode in total_reactions:
                rId = rNode.id
                metacyc_id = rId
                if rId in found_reactions:
                    in_network = "yes"
                else:
                    in_network = "no"
                try:
                    ec = rNode.misc["EC_NUMBER"][0]
                except KeyError:
                    ec = "Unknown"
                try:
                    common_name = rNode.misc["COMMON_NAME"][0]
                except KeyError:
                    common_name = "Unknown"
                try:
                    direction = rNode.misc["DIRECTION"][0]
                except KeyError:
                    direction = " =>/<=> "
                if direction == "REVERSIBLE":
                    direction = " <=> "
                elif direction == "LEFT-TO-RIGHT":
                    direction = " => "
        
                reactants = [rlt.misc["STOICHIOMETRY"][0]+" "+rlt.id_out 
                for rlt in padmetRef.dicOfRelationIn.get(rId, None) if rlt.type == "consumes"]
                products = [rlt.misc["STOICHIOMETRY"][0]+" "+rlt.id_out
                for rlt in padmetRef.dicOfRelationIn.get(rId, None) if rlt.type == "produces"]
                metIdFormula = " + ".join(reactants)+direction+" + ".join(products)
                
                try:
                    reactants = [rlt.misc["STOICHIOMETRY"][0]+" "+padmetRef.dicOfNode[rlt.id_out].misc["COMMON_NAME"][0] 
                    for rlt in padmetRef.dicOfRelationIn.get(rId, None) if rlt.type == "consumes"]
                    products = [rlt.misc["STOICHIOMETRY"][0]+" "+padmetRef.dicOfNode[rlt.id_out].misc["COMMON_NAME"][0]
                    for rlt in padmetRef.dicOfRelationIn.get(rId, None) if rlt.type == "produces"]
                    cnameFormula = " + ".join(reactants)+direction+" + ".join(products)
                except KeyError:
                    cnameFormula = "Unknown"
                    
                line = "\t".join([metacyc_id,common_name,ec,metIdFormula,cnameFormula,in_network])
                line = line+"\n"
                f.write(line)
        if sbml is not None:
        #generate an sbml for the extracted data
            print("Generating SBML file...")
            #TODO

        
    def network_report(self, output_dir, padmetRef_file = None, verbose = False):
        """
        Summurizes the network in a folder (output_dir) of 4 files.
        all_pathways.csv: report on the pathways of the network. PadmetRef is used to recover the total reactions of a pathways. (sep = "\t")
        line = dbRef_id, Common name, Number of reactions found, 
            Total number of reaction, Ratio (Reaction found / Total)
        all_reactions.csv: report on the reactions of the network.  (sep = "\t")
        line = dbRef_id, Common name, formula (with id), 
            formula (with common name), in pathways, associated genes, sources
        all_metabolites.csv: report on the metabolites of the network. (sep = "\t")
        line = dbRef_id, Common name, Produced (p), Consumed (c), Both (cp)
        all_genes.csv: report on the genes of the network. (sep= "\t")
        line = "id", "Common name", "linked reactions"
        @param padmetRef_file: pathname of the padmet of reference
        @param output_dir: pathname of the folder where to create the reports
        @param verbose: if true print info.
        @type padmetRef_file, output_dir: str
        @type verbose: bool
        @return: _
        @rtype: None
        """
        os.system("mkdir -p "+output_dir)
        if not output_dir.endswith("/"):
            output_dir += "/"
        all_pathways = output_dir+"all_pathways.csv"
        all_reactions = output_dir+"all_reactions.csv"
        all_metabolites = output_dir+"all_metabolites.csv"
        all_genes = output_dir+"all_genes.csv"
        
        if padmetRef_file is not None:
            padmetRef = PadmetRef(padmetRef_file)
            pathways = set([self.dicOfNode[rlt.id_out] for rlt in self.getAllRelation() 
            if rlt.type == "is_in_pathway" and self.dicOfNode[rlt.id_in].type == "reaction"])
            nb_pathways = len(pathways)

        reactions = [node for node in self.dicOfNode.values()
        if node.type == "reaction"]
        nb_reactions = len(reactions)

        metabolites = set([self.dicOfNode[rlt.id_out] 
        for rlt in self.getAllRelation()
        if rlt.type == "consumes" or rlt.type == "produces"])
        nb_metabolites = len(metabolites)
        
        genes = [node for node in self.dicOfNode.values()
        if node.type == "gene"]
        nb_genes = len(genes)

        if padmetRef_file is not None:
            with open(all_pathways,'w') as f:
                header = ["dbRef_id", "Common name", "Number of reaction found", 
                "Total number of reaction", "Ratio (Reaction found / Total)"]
                header = "\t".join(header)+"\n"
                f.write(header)
    
                count = 1
                for pnode in pathways:
                    pnode_id = pnode.id
                    if verbose: print("Pathway : %s %s/%s" % (pnode_id, count, nb_pathways))
                    # Recover the first common_name
                    try:
                        common_name = pnode.misc["COMMON_NAME"][0]
                    except KeyError:
                        common_name = "Unknown"
                    
                    number_of_reaction_found = len([rlt for rlt in self.dicOfRelationOut[pnode_id]
                    if rlt.type == "is_in_pathway"
                    and self.dicOfNode[rlt.id_in].type == "reaction"])
    
                    try:
                        total_number_of_reaction = len([rlt for rlt in padmetRef.dicOfRelationOut[pnode_id]
                        if rlt.type == "is_in_pathway"
                        and padmetRef.dicOfNode[rlt.id_in].type == "reaction"])
                        # If keyError: pathway not in padmetRef, pathway added manualy
                    except KeyError: 
                        total_number_of_reaction = "NA"
                        
                    try:
                        if type(total_number_of_reaction) is int:
                            ratio = "%.2f" % (float(number_of_reaction_found)/total_number_of_reaction)
                        else:
                            ratio = "NA"
                        line = [pnode_id,common_name, str(number_of_reaction_found), 
                                str(total_number_of_reaction), str(ratio)]
                        line = "\t".join(line)+"\n"
                        f.write(line)
                    except ZeroDivisionError:
                        pass
                    count += 1
                
        with open(all_reactions,'w') as f:
            header = ["dbRef_id", "Common name", "formula (with id)", 
            "formula (with common name)", "in pathways", "associated genes", "sources"]
            header = "\t".join(header)+"\n"
            f.write(header)

            count = 0
            for rnode in reactions:
                count += 1
                rnode_id = rnode.id
                if verbose: print("Reaction : %s %s/%s" % (rnode_id, count, nb_reactions))
                # Recover the first common_name
                try:
                    common_name = rnode.misc["COMMON_NAME"][0]
                except KeyError:
                    common_name = "Unknown"
                # Recovering pathways associated
                try:
                    in_pathways = ";".join([rlt.id_out for rlt in self.dicOfRelationIn.get(rnode_id, None) 
                    if rlt.type == "is_in_pathway"])
                except TypeError:
                    in_pathways = ""
                # Recovering the formula
                direction = rnode.misc["DIRECTION"][0]
                if direction == "UNKNOWN":
                    direction = " =>/<=> "
                elif direction == "REVERSIBLE":
                    direction = " <=> "
                elif direction == "LEFT-TO-RIGHT":
                    direction = " => "
                
                reactants = [rlt.misc["STOICHIOMETRY"][0]+" "+rlt.id_out 
                for rlt in self.dicOfRelationIn[rnode_id] if rlt.type == "consumes"]
                products = [rlt.misc["STOICHIOMETRY"][0]+" "+rlt.id_out
                for rlt in self.dicOfRelationIn[rnode_id] if rlt.type == "produces"]
                formula_id = " + ".join(reactants)+direction+" + ".join(products)
                
                try:
                    reactants = [rlt.misc["STOICHIOMETRY"][0]+" "+self.dicOfNode[rlt.id_out].misc["COMMON_NAME"][0] 
                    for rlt in self.dicOfRelationIn[rnode_id] if rlt.type == "consumes"]
                    products = [rlt.misc["STOICHIOMETRY"][0]+" "+self.dicOfNode[rlt.id_out].misc["COMMON_NAME"][0]
                    for rlt in self.dicOfRelationIn[rnode_id] if rlt.type == "produces"]
                    formula_cname = " + ".join(reactants)+direction+" + ".join(products)
                except KeyError:
                    formula_cname = ""
                
                # Recovering genes associated
                try:
                    gene_assoc = ";".join([rlt.id_out 
                    for rlt in self.dicOfRelationIn[rnode_id]
                    if rlt.type == "is_linked_to"])
                except TypeError:
                    gene_assoc = ""
                
                sources = ""
                if "SOURCE" in rnode.misc.keys():
                    sources = rnode.misc["SOURCE"][0]
                else:
                    sources = []
                    for rlt in self.dicOfRelationIn[rnode_id]:
                        if rlt.type == "has_suppData":
                            src = os.path.splitext(self.dicOfNode[rlt.id_out].misc["ORIGIN_FILE"][0])[0]
                            sources.append(src)
                    sources = ";".join(sources)

                line = "\t".join([rnode_id, common_name, formula_id, formula_cname, in_pathways, gene_assoc, sources])
                line = line+"\n"
                f.write(line)

        with open(all_metabolites,'w') as f:
            header = ["dbRef_id", "Common name", "Produced (p), Consumed (c), Both (cp)"]
            header = "\t".join(header)+"\n"
            f.write(header)

            count = 0
            for mnode in metabolites:
                count += 1
                mnode_id = mnode.id
                is_consumed = False
                is_produced = False
                if verbose: print("Metabolite : %s %s/%s" % (mnode_id, count, nb_metabolites))
                # Recover the first common_name
                try:
                    common_name = mnode.misc["COMMON_NAME"][0]
                except KeyError:
                    common_name = "Unknown"
                rlts_c_and_p = [(rlt.type,rlt.id_in) 
                for rlt in self.dicOfRelationOut[mnode_id]
                if rlt.type == "consumes" or rlt.type == "produces"]

                if "REVERSIBLE" in [self.dicOfNode[rxn_id].misc["DIRECTION"][0] 
                for (rlt_type, rxn_id) in rlts_c_and_p]:
                    statut = "cp"
                else:
                    for rlt_type, rxn_id in rlts_c_and_p:
                        if rlt_type == "consumes":
                            is_consumed = True
                        elif rlt_type == "produces":
                            is_produced = True
                    if is_consumed and is_produced:
                        statut = "cp"
                    elif is_consumed:
                        statut = "c"
                    elif is_produced:
                        statut = "p"
                    else:
                        statut = "NA"

                line = "\t".join([mnode_id, common_name, statut])
                line = line+"\n"
                f.write(line)

        with open(all_genes,'w') as f:
            header = ["id", "Common name", "linked reactions"]
            header = "\t".join(header)+"\n"
            f.write(header)

            count = 0
            for gnode in genes:
                count += 1
                gnode_id = gnode.id
                if verbose: print("Gene : %s %s/%s" % (gnode_id, count, nb_genes))
                # Recover the first common_name
                try:
                    common_name = gnode.misc["COMMON_NAME"][0]
                except KeyError:
                    common_name = "Unknown"
                
                # Recovering linked reactions
                try:
                    rxn_assoc = ";".join([rlt.id_in 
                    for rlt in self.dicOfRelationOut.get(gnode_id, None) 
                    if rlt.type == "is_linked_to"])
                except TypeError:
                    rxn_assoc = ""

                line = "\t".join([gnode_id, common_name, rxn_assoc])
                line = line+"\n"
                f.write(line)
            
                
#==============================================================================
# For Nodes 
#==============================================================================
            
    def _addNode(self,node):
        """
        Allows to add a node, only if the id is not already used.
        @param node: the node to add
        @type node: Node
        @return: True if added, False if no.        
        @rtype: Bool
        """
        if node.id not in self.dicOfNode.keys():
            self.dicOfNode[node.id] = node
            return True
        else:
            return False

        
    def _copyNodeExtend(self, padmet, node_id):
        """
        Allows to copy a node from an other padmet with the first childrens only.
        Recursive function, call itself for the relations where the node is "in"
        NB: particular case: we dont want to recovere the relations "prot catalyses reaction"
        do nothing for the relations where the node is "out"
        @param padmet: Padmet from where to copy the node.
        @param node_id: the id of the node to copy.
        @type padmetRef: PadmetSpec/Ref
        @type node_id: str
        @return: _
        @rtype: None
        """
        node = padmet.dicOfNode[node_id]
        # If true, the node does not exist and is susccesffuly added.
        if self._addNode(node):
            try:
                relations_in = (rlt for rlt in padmet.dicOfRelationIn.get(node_id, None)
                if rlt.type != "catalyses") 
                # For each relation, we add it with addRelation(), then
                # call _copyNodeExtend() for the id 'out'.
                for rlt in relations_in:
                    self._addRelation(rlt)
                    node_out_id = rlt.id_out
                    node_out_type = padmet.dicOfNode[node_out_id].type
                    if (node_out_type == "xref" or node_out_type == "name" 
                    or node_out_type == "suppData"):
                        self._addNode(padmet.dicOfNode[node_out_id])
                    else:
                        self._copyNodeExtend(padmet, node_out_id)                        
            # TypeError is raise if tgdb.getRealtion() is none, not relation
            # is found.
            except TypeError:
                pass

            
    def copyNode(self, padmet, node_id):
        """
        copyNode() allows to copy a node from an other padmetSpec or padmetRef. It copies all 
        the relations 'in' and 'out' and it calls the function 
        _copyNodeExtend() to recover the associated node.
        @param Padmet: PadmetSpec/Ref from where to copy the node
        @param node_id: the id of the node to copy
        @type padmetRef: PadmetSpec/Ref
        @type node_id: str
        @return: _
        @rtype: None
        """
        try:
            nodeToCopy = padmet.dicOfNode[node_id]
        except KeyError:
            raise TypeError("Unable to get the node(s) with id: %s" %node_id)

        
        # True => node succeffuly added in dictionary of node.
        if self._addNode(nodeToCopy):
            # Recover all relations of the node
            try:
                relations_in = (rlt for rlt in padmet.dicOfRelationIn.get(node_id, None))
                # Add each relation with _addRelation, then call
                # _copyNodeExtend() for the id 'out'.
                for rlt in relations_in:
                    self._addRelation(rlt)
                    node_out_id = rlt.id_out
                    node_out_type = padmet.dicOfNode[node_out_id].type
                    if (node_out_type == "xref" or node_out_type == "name" 
                    or node_out_type == "suppData"):
                        self._addNode(padmet.dicOfNode[node_out_id])
                    else:
                        self._copyNodeExtend(padmet, node_out_id)        
            # TypeError is raise if padmet.getRealtion() is none, no relation
            # found. (NoneType not iterable)
            except TypeError:
                pass                      

            try:
                # Recover all relations where the node is out.
                relations_out = (rlt for rlt in padmet.dicOfRelationOut.get(node_id, None)
                if rlt.type != "catalyses")
                # For each relation, we add it with addRelation, then call
                # _copyNodeExtend() for the id 'in'
                for rlt in relations_out:
                    self._addRelation(rlt)
                    node_in_id = rlt.id_in
                    self._copyNodeExtend(padmet, node_in_id)
            # TypeError is raised if tgdb.getRealtion() is none, not relation
            # found.
            except TypeError:
                pass

                
    def delNode(self, node_id):
        """
        Allows to delete a node, the relations associated to the node, and for 
        some relations, delete the associated node. 
        For relations where the node to del is 'in': 
            if rlt type in ['has_xref','has_name','has_suppData']: delNode out
        For relations where the node to del is 'out':
            if rlt type in ['consumes','produces']
        @param node_id: id of node to delete
        @type node_id: str
        @return: True if node successfully deleted, False if node not in dicOfNode
        @rtype: Bool
        """
        # Delete the node from dicOfNode.
        try:
            self.dicOfNode.pop(node_id)
        except KeyError:
            print("The id %s doesnt exist. Unable to delete" %node_id)
            return False
            
        # If exist delete the relations 'in' and 'out'
        try:
            # Recover the relations where the node is "in".
            relationsIn = [rlt for rlt in self.dicOfRelationIn.get(node_id, None)]
            for rltIn in relationsIn:
                #print(rltIn.toString())
                self._delRelation(rltIn)
                if (rltIn.type in ["has_xref","has_name","has_suppData"]):
                    self.delNode(rltIn.id_out)
        except TypeError:
            pass
        try:
            # Recover the relations where the node is "out"
            relationsOut = [rlt for rlt in self.dicOfRelationOut.get(node_id, None)]
            for rltOut in relationsOut:
                #print(rltOut.toString())
                self._delRelation(rltOut)
        except TypeError:
            pass
        return True

    
    def updateNode(self, node_id, data, action, verbose = False):
        """
        Allows to update miscellaneous data of a Node.
        @param node_id: the id of node to update
        @param data: tuple with data[0] refere to the miscellaneous data key 
        (ex: common_name, direction ...), data[1] is a list of value to add / update.
        data[1] can be None if the action is to pop the key
        @param action: if == "add": the list data[1] wil be added (ex: adding a common_name)
        if == "remove": if data[1] is not None, the list data[1] will be removed (ex: remove just one specifique common_name)
        if == "update":data[1] is the new value of the key data[0]
        ex: updateNode('RXN-5',('direction',['LEFT-TO-RIGHT']),update). The 
        reaction' direction will be change to left2right
        @param node_id: the id of the node to update
        @param data: tuple of data to update, data[0] is the key, data[1] is a value, list or None
        @param action: action in ['add','pop','remove','update']. Check description for more information
        @param verbose: print more info
        @type node_id, action: str
        @type data: list or None
        @type verbose: Bool
        @return: True if successfully updated, False if no
        @rtype: Bool
        """
        try:
            node = self.dicOfNode[node_id]
        except KeyError:
            print("The id %s doesnt exist. Unable to update" %node_id)
            return False

        if action == "add":
            try:
                oldData = node.misc[data[0]]
                oldData = oldData + data[1]
                oldData = set(oldData)
                newData = list(oldData)
                node.misc[data[0]] = newData
            except KeyError:
                node.misc[data[0]] = data[1]
    
        elif action == "remove":
            try:
                newData = [x for x in oldData if x not in data[0]]
                if len(newData) != 0:
                    node.misc[data[0]] = newData
                elif len(newData) == oldData:
                    print("None of the data to deleted is present")
                else:
                    node.misc.pop(data[0])
            except KeyError:
                print("No data %s associated to this node" %data[0])
                return False
        
        elif action == "update":
            node.misc[data[0]] = data[1]
            
        else:
            print("Unknown action: %s, please check docstring" %action)
            return False
        self.dicOfNode[nodeId] = node
        return True
            
#==============================================================================
# For Relations:     
#==============================================================================

    def _delRelation(self, relation):
        """
        Delete a relation from dicOfRelationIn and out
        @param relation: the relation to delete
        @type relation: Relation
        @return: True if succesfully deleted
        @rtype: Bool
        """
        delete = False
        idIn = relation.id_in
        idOut = relation.id_out
        try:
            for rlt in self.dicOfRelationIn[idIn]:
                if relation.compare(rlt):
                    self.dicOfRelationIn[idIn].remove(rlt)
                    delete = True
                    break
        except KeyError: pass
        try:
            for rlt in self.dicOfRelationOut[idOut]:
                if relation.compare(rlt):
                    self.dicOfRelationOut[idOut].remove(rlt)
                    delete = True
                    break
        except KeyError: pass
        if not delete:
            print("Unable to delete this relation, doesn't exist.")
        return delete

            
    def _addRelation(self,relation):
        """
        AddRelation() allows to add a relation if not already in allRelations.
        @param relation: the relation to add
        @type relation: Relation
        @return: true if relation was successfully added
        @rtype: bool
        """
        idIn = relation.id_in
        idOut = relation.id_out
        rlt_to_compare = self.dicOfRelationIn.get(idIn,[])
        rlt_to_compare += self.dicOfRelationOut.get(idOut,[])
        for rlt in rlt_to_compare:
            if relation.compare(rlt):
                return False
        try:
            self.dicOfRelationIn[idIn].append(relation)    
        except KeyError:
            self.dicOfRelationIn[idIn] = [relation]    
        try:
            self.dicOfRelationOut[idOut].append(relation)
        except KeyError:
            self.dicOfRelationOut[idOut] = [relation]
        return True

#==============================================================================
# manipulating de novo node:     
#==============================================================================
    
    def _basicNode(self, _type):
        """
        For padmetSpec, when creating a new node, a new id is creating. This id
        start with 'SPE_' for padmetSpec (META_ for padmetRef).
        This function generate a new id and an empty node with this id.
        @param _type: the type of the node to create
        @return: (newId,newNode)
        @rtype: tuple(str, Node)
        """
        #check if the type is define in the policy
        if _type in self.policy.getClassOfNode():
            #define a new unique_id:
            #generator of the int part of ids that contains the species tag
            max_local_id = self.getMaxLocalID()
            new_id = "SPE_"+str(max_local_id + 1)
            new_node = Node(_type, new_id)
            return(new_id, new_node)
        else:
            raise TypeError("the given type is not define in the policy")


    def getMaxLocalID(self):
        """
        For padmetSpec, when creating a new node, a new id is creating. This id
        start with 'SPE_' for padmetSpec (META_ for padmetRef) + an incremented int.
        This function extracts the max int (or max local id)
        @return: the max local id
        @rtype: int
        """
        try:
            max_local_id = max([int(node_id.replace(("SPE_"),"")) 
            for node_id in self.dicOfNode.iterkeys() 
            if node_id.startswith("SPE_")])
        except ValueError:
            max_local_id  = 0
        return max_local_id

            
    def createNode(self, _type, dicOfMisc, listOfRelation = None):
        """
        Creation of new node to add in the network.
        use ._basicNode first then completes the node with more informations
        @param _type: type of node (gene, reaction...)
        @param dicOfMisc: dictionnary of miscellaneous data
        @param listOfRelation: list of list of data needed to create a relation. (id_in,type,id_out,misc)
        @type _type: str
        @type dicOfMisc: dict
        @type listOfRelation: def = None, List.
        @return: (new_id,new_node)
        @rtype: tuple(str, Node)
        """
        (new_id, new_node) = self._basicNode(_type)
        #add the new node in the tgdbp
        new_node.misc = dicOfMisc
        self._addNode(new_node)
        if listOfRelation is not None:
            for data in listOfRelation:
                try:
                    #'_self' in a relation is where to put the id of the current new node (if needed)
                    self_index = data.index("_self")
                    data[self_index] = new_id
                except ValueError:
                    pass
                try:
                    new_relation = Relation(data[0],data[1],data[2],data[3])
                except IndexError:
                    new_relation = Relation(data[0],data[1],data[2])
                self._addRelation(new_relation)

        return (new_id, new_node)