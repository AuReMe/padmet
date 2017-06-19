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
PadmetRef is an object representing a DATABASE of metabolic network.
"""
from policy import Policy
from node import Node
from relation import Relation
from sbmlPlugin import convert_from_coded_id

try:
    from libsbml import *
except:
    print("package libsbml needed, use this cmd:\n pip install "
          + "python-libsbml")
    exit()

class PadmetRef:
    """
    PadmetRef is an object representing a DATABASE of metabolic network.
    Contains <Policy>, <Node> and <Relation>
    The policy defines the way Node and Relation are associated
    A node is an Object that contains information about an element of the network 
    (can be a pathway, reaction...).
    A realtion defines how two nodes are connected. In a relation there is 
    a node "in" and a node "out". (reactionX'in' consumes metaboliteX'out')
    PadmetRef contains 3 attributs: 
        dicOfNode: a dictionary of node: key=Node's unique id / value = <Node>
        dicOfRelationIn: a dictionnary of relation with: key= nodeIN id / value = list of <relation>
        dicOfRelationOut: a dictionnary of relation with: key= nodeOut id / value = list of <relation>
        policy: a <policy>
        info: a dictionnary of informations about the network, the database used...
        This dictionnary is always represented in the header of a padmet file
    """
    def __init__(self, padmetRef_file = None):
        """
        if None, initializes an empty <PadmetRef>
        @param padmetRef_file: pathname of the padmet file
        @type padmetRef_file: str
        @return: _
        @rtype: None
        """
        if padmetRef_file is not None:        
            self.loadGraph(padmetRef_file)
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
        
    def initFromSbml(self, sbml_file, verbose = False):
        """
        Initialize a padmetRef from sbml. Copy all species, convert id with sbmlPlugin
        stock name in COMMON NAME, stock original info in suppData. Copy all reactions,
        convert id with sbmlPlugin, stock name in common name, stock compart and stoichio data relative
        to reactants and products in the misc of consumes/produces relations
        @param sbml_file: pathname of the sbml file
        @param verbose: <bool> if True print supp info
        @type sbml_file: str
        """
        #using libSbml to read sbml_file
        if verbose: print ("loading sbml file: %s" %sbml_file)
        reader = SBMLReader()
        document = reader.readSBML(sbml_file)
        for i in range(document.getNumErrors()):
            print(document.getError(i).getMessage())
        
        model = document.getModel()
        listOfSpecies = model.getListOfSpecies()
        listOfReactions = model.getListOfReactions()
        nbReactions = len(listOfReactions)
        nbSpecies = len(listOfSpecies)
        if verbose:
            print("nb species: %s" %nbSpecies)
            print("nb reactions: %s" %nbReactions) 
            
        if verbose: print("creating species")
        for specie in listOfSpecies:
            specie_id_encoded = specie.getId()
            specie_id = convert_from_coded_id(specie_id_encoded)[0]
            if verbose: print("specie: %s, uncoded: %s" % (specie_id_encoded, specie_id))
            try:
                self.dicOfNode[specie_id]
                if verbose: print("already in padmetRef")
            except KeyError:
                specie_name = specie.getName()
                if specie_name:
                    specie_node = Node("compound", specie_id, {"COMMON_NAME": [specie_name]})
                else:
                    specie_node = Node("compound", specie_id)
                self.dicOfNode[specie_id] = specie_node
    
        if verbose: print("creating reactions")
        for reaction in listOfReactions:
            reaction_id_encoded = reaction.getId()
            reaction_id = convert_from_coded_id(reaction_id_encoded)[0]
            if verbose: print("reaction: %s, uncoded: %s" % (reaction_id_encoded, reaction_id))
            try:
                self.dicOfNode[reaction_id]
                if verbose: print("already in padmetRef")
            except KeyError:
                reaction_name = reaction.getName()
                if reaction.getReversible():
                    reaction_dir = "REVERSIBLE"
                else:
                    reaction_dir = "LEFT-TO-RIGHT"
                if reaction_name:
                    reaction_node = Node("reaction", reaction_id, {"COMMON_NAME": [reaction_name], "DIRECTION": [reaction_dir]})
                else:
                    reaction_node = Node("reaction", reaction_id, {"DIRECTION": [reaction_dir]})
                self.dicOfNode[reaction_id] = reaction_node
                reactants = reaction.getListOfReactants()
                for reactant in reactants:
                    reactant_id, x, reactant_compart = convert_from_coded_id(reactant.getSpecies())
                    if reactant_compart is None:
                        if verbose: print("%s has no compart, set to 'c'" %reactant)
                        reactant_compart = "c"
                    reactant_stoich = reactant.getStoichiometry()
                    consumes_rlt = Relation(reaction_id, "consumes", reactant_id, {"STOICHIOMETRY": [reactant_stoich], "COMPARTMENT": [reactant_compart]})
                    self._addRelation(consumes_rlt)
                products = reaction.getListOfProducts()
                for product in products:
                    product_id, x, product_compart = convert_from_coded_id(product.getSpecies())
                    if product_compart is None:
                        if verbose: print("%s has no compart, set to 'c'" %product)
                        product_compart = "c"
                    product_stoich = product.getStoichiometry()
                    produces_rlt = Relation(reaction_id, "produces", product_id, {"STOICHIOMETRY": [product_stoich], "COMPARTMENT": [product_compart]})
                    self._addRelation(produces_rlt)


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
        orderedRelation = tuple(sorted(self.getAllRelation(), key=lambda x: x.id_in, reverse=False))
        with open(output, 'w') as f:
            if len(self.info) != 0:
                f.write("Data Base informations\n")
                for k,data in self.info.iteritems():
                    f.write(k+":\n")
                    for k,v in data.iteritems():
                       f.write("\t"+k+":"+v+"\n")
                f.write("\n")
            # It writes the policy of the TGDBP file.
            f.write("Policy\n")
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
            for rlt in orderedRelation:
                line = rlt.toString()+"\n"
                f.write(line)

    def extract_data(self,output_directory, verbose = False):
        """
        extracting data on rections and compounds in flate files.(used for samifier)
        """
        
        all_reactions = output_directory+"all_reactions.ts"
        all_metabolites = output_directory+"all_metabolites.tsv"
        
        #Recovering metabolites
        if verbose: print("Recovering metabolites")
        metabolites_nodes = set([self.dicOfNode[rlt.id_out] for rlt in self.getAllRelation()
        if rlt.type in ["consumes","produces"]])
        
        #recovere all the compounds that are in a class which is involved in a reaction
        metabolites_sub_class = []
        for node in metabolites_nodes:
            if node.type == "class":
                node_id = node.id
                metabolites_sub_class += [self.dicOfNode[rlt.id_in] for rlt in self.dicOfRelationOut.get(node_id, None)
                if rlt.type == "is_a_class"]
        metabolites_nodes = metabolites_nodes.union(set(metabolites_sub_class))
        
        #Recovering reactions
        if verbose: print("Recovering reactions")
        reactions_nodes = [node for node in self.dicOfNode.itervalues()
        if node.type == "reaction"]
        
        if verbose: print("Metabolites...")
        count = 0
        with open(all_metabolites,'w') as f:
            header = "\t".join(["METACYC", "SYNONYMS","XREF"])+"\n"
            f.write(header)
            for node in metabolites_nodes:
                count += 1
                if verbose: print("Metabolite "+str(count)+"/"+str(len(metabolites_nodes)))
                
                metacyc_id = node.id

                try:
                    common_name = node.misc["COMMON_NAME"]
                except KeyError:
                    common_name = []

                try:
                    names = [self.dicOfNode[rlt.id_out].misc["LABEL"][0] for rlt in self.dicOfRelationIn.get(metacyc_id, None) 
                    if rlt.type == "has_name"]
                except TypeError:
                    names = []
                synonyms = ";".join(common_name+names)

                try:
                    INCHIKEY = node.misc["INCHI_KEY"][0]
                    INCHIKEY = [INCHIKEY.replace("InChIKey=","INCHIKEY:")]
                except KeyError:
                    INCHIKEY = []

                try:
                    SMILES = node.misc["SMILES"][0]
                    SMILES = ["SMILES:"+SMILES]
                except KeyError:
                    SMILES = []

                try:
                    xrefs_node = [self.dicOfNode[rlt.id_out] for rlt in self.dicOfRelationIn.get(metacyc_id, None) 
                    if rlt.type == "has_xref"]
                    xrefs = [":".join((node.misc["DB"][0].strip(),node.misc["ID"][0])) for node in xrefs_node]
                except TypeError:
                    xrefs = []
                all_xrefs = ";".join(INCHIKEY+SMILES+xrefs)
                line = "\t".join([metacyc_id, synonyms, all_xrefs])+"\n"
                f.write(line)

        if verbose: print("Reactions...")
        count = 0
        with open(all_reactions,'w') as f:
            header = "\t".join(["METACYC", "SYNONYMS","XREF", "REAGS", "PRODS", "REV"])+"\n"
            f.write(header)
            for node in reactions_nodes:
                count += 1
                if verbose: print("Reaction "+str(count)+"/"+str(len(reactions_nodes)))
                
                metacyc_id = node.id

                try:
                    common_name = node.misc["COMMON_NAME"]
                except KeyError:
                    common_name = []

                try:
                    names = [self.dicOfNode[rlt.id_out].misc["LABEL"][0] for rlt in self.dicOfRelationIn.get(metacyc_id, None) 
                    if rlt.type == "has_name"]
                except TypeError:
                    names = []
                synonyms = ";".join(common_name+names)

                try:
                    ec = node.misc["EC_NUMBER"]
                except KeyError:
                    SMILES = []
                
                rev = node.misc["DIRECTION"][0]
                if rev == "LEFT-TO-RIGHT":
                    rev_bool = "0"
                else:
                    rev_bool = "1"
                
                try:
                    xrefs_node = [self.dicOfNode[rlt.id_out] for rlt in self.dicOfRelationIn.get(metacyc_id, None) 
                    if rlt.type == "has_xref"]
                    xrefs = [":".join((node.misc["DB"][0].strip(),node.misc["ID"][0])) for node in xrefs_node]
                except TypeError:
                    xrefs = []
                all_xrefs = ";".join(ec+xrefs)
                
                try:
                    reactants = [rlt.id_out for rlt in self.dicOfRelationIn.get(metacyc_id, None)
                    if rlt.type == "consumes"]
                    reactants = ";".join(reactants)
                except TypeError:
                    reactants = ""

                try:
                    products = [rlt.id_out for rlt in self.dicOfRelationIn.get(metacyc_id, None)
                    if rlt.type == "produces"]
                    products = ";".join(products)
                except TypeError:
                    products = ""
                line = "\t".join([metacyc_id, synonyms, all_xrefs, reactants, products, rev_bool ])+"\n"
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
            
           
#==============================================================================
# For Relations:     
#==============================================================================

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
        For padmetRef, when creating a new node, a new id is creating. This id
        start with 'META_' for padmetRef (SPE_ for padmetSpe).
        This function generate a new id and an empty node with this id.
        @param _type: the type of the node to create
        @return: (newId,newNode)
        @rtype: tuple(str, Node)
        """
        #check if the class is define in the policy
        if _type in self.policy.getClassOfNode():
            #define a new unique_id:
            #generator of the int part of ids that contains the species tag
            max_local_id = self.getMaxLocalID()
            #localPrefix = self.getLocalPrefix()
            new_id = "META_"+str(max_local_id + 1)
            new_node = Node(_type, new_id)
            return(new_id, new_node)
        else:
            raise TypeError("the type given is not define in the policy")


    def getMaxLocalID(self):
        """
        For padmetRef, when creating a new node, a new id is creating. This id
        start with 'META_' for padmetRef (SPE_ for padmetSpec) + an incremented int.
        This function extracts the max int (or max local id)
        @return: the max local id
        @rtype: int
        """
        try:
            #localPrefix = self.getLocalPrefix()
            max_local_id = max([int(node_id.replace(("META_"),"")) 
            for node_id in self.dicOfNode.iterkeys() 
            if node_id.startswith("META_")])
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