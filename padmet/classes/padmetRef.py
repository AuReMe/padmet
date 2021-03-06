#!/usr/bin/env python
# -*- coding: utf-8 -*-
from padmet.classes.policy import Policy
from padmet.classes.node import Node
from padmet.classes.relation import Relation
from padmet.utils import sbmlPlugin
import libsbml
import os

__all__ = ["PadmetRef"]


class PadmetRef:
    """
    PadmetRef is an object representing a DATABASE of metabolic network.
    
    Contains <Policy>, <Node> and <Relation>:
    
    The policy defines the way Node and Relation are associated.

    A node is an Object that contains information about an element of the network 
    (can be a pathway, reaction...).

    A realtion defines how two nodes are connected. In a relation there is 
    a node "in" and a node "out". (reactionX'in' consumes metaboliteX'out').

    PadmetRef contains 3 attributs: 
        dicOfNode: a dictionary of node: key=Node's unique id / value = <Node>.

        dicOfRelationIn: a dictionnary of relation with: key= nodeIN id / value = list of <relation>.

        dicOfRelationOut: a dictionnary of relation with: key= nodeOut id / value = list of <relation>.

        policy: a <policy>.

        info: a dictionnary of informations about the network, the database used...
        This dictionnary is always represented in the header of a padmet file.
    """

    def __init__(self, padmetRef_file=None):
        """
        if None, initializes an empty <PadmetRef>
        
        Parameters
        ----------
        padmetRef_file: str
            pathname of the padmet file
        
        """
        if padmetRef_file is not None:
            if not os.path.exists(padmetRef_file):
                raise FileNotFoundError("No Padmet Ref file accessible at " + padmetRef_file)
            self.loadGraph(padmetRef_file)
        else:
            self.dicOfRelationIn = {}
            self.dicOfRelationOut = {}
            self.dicOfNode = {}
            self.policy = Policy()
            self.info = {}

    # ==============================================================================
    # Constructor / getter
    # ==============================================================================

    def setInfo(self, source):
        """
        All the information printed in the header of a padmet stocked in a dict.
        {"metacyc":{version:XX,...},"ecocyc":{...}...}
        set Info from a dictionnary or copying from an other padmet

        Parameters
        ----------
        source: dict or padmet.classes.PadmetRef
            may be a dict or an other padmet from where will be copied the info
        """
        if type(source) is dict:
            self.info = source
        else:
            self.info = source.info

    def setPolicy(self, source):
        """
        Set policy from a list or copying from an other padmet

        Parameters
        ----------
        source: list or padmet.classes.PadmetRef
            may be a list or an other padmet from where will be copied the policy
        """
        if type(source) is list:
            self.policy = Policy(source)
        else:
            self.policy = source.policy

    def setDicOfNode(self, source):
        """
        Set dicOfNode from a dict or copying from an other padmet

        Parameters
        ----------
        source: dict or padmet.classes.PadmetRef
            may be a dict or an other padmet from where will be copied the dicOfNode
        """
        if type(source) is dict:
            self.dicOfNode = source
        else:
            self.dicOfNode = source.dicOfNode

    def setdicOfRelationIn(self, source):
        """
        Set dicOfRelationIn from a dict or copying from an other padmet

        Parameters
        ----------
        source: dict or padmet.classes.PadmetRef
            may be a dict or an other padmet from where will be copied the dicOfRelationIn
        """
        if type(source) is dict:
            self.dicOfRelationIn = source
        else:
            self.dicOfRelationIn = source.dicOfRelationIn

    def setdicOfRelationOut(self, source):
        """
        Set dicOfRelationOut from a dict or copying from an other padmet
        
        Parameters
        ----------
        source: dict or padmet.classes.PadmetRef
            may be a dict or an other padmet from where will be copied the dicOfRelationOut
        """
        if type(source) is dict:
            self.dicOfRelationOut = source
        else:
            self.dicOfRelationOut = source.dicOfRelationOut

    def getAllRelation(self):
        """

        Returns
        -------
        set
            return a set of all relations
        """
        all_relation = set()
        for list_rlt in self.dicOfRelationIn.values():
            for rlt in list_rlt:
                all_relation.add(rlt)
        for list_rlt in self.dicOfRelationIn.values():
            for rlt in list_rlt:
                all_relation.add(rlt)
        return all_relation

    def getCompounds(self, get_id=True):
        """
        Get all the compounds from a padmet instance.

        Parameters
        ----------
        self: padmet
            padmet instance
        get_id: bool
            True for the ID and False for the node corresponding to the compounds

        Returns
        -------
        list
            list of compound IDs or nodes
        """
        total_cpd_id = set()

        all_rxns = [node for node in self.dicOfNode.values() if node.type == "reaction"]

        for rxn_node in all_rxns:
            total_cpd_id.update([rlt.id_out for rlt in self.dicOfRelationIn[rxn_node.id] if rlt.type in ["consumes","produces"]])

        if get_id is True:
            all_cpds = [node_id for (node_id, node) in self.dicOfNode.items() if node_id in total_cpd_id]
        else:
            all_cpds = [node for (node_id, node) in self.dicOfNode.items() if node_id in total_cpd_id]

        return all_cpds

    def getReactions(self, get_id=True):
        """
        Get all the reactions from a padmet instance.

        Parameters
        ----------
        self: padmet
            padmet instance
        get_id: bool
            True for the ID and False for the node corresponding to the reactions

        Returns
        -------
        list
            list of reactions IDs or nodes
        """
        if get_id is True:
            return [node.id for node in self.dicOfNode.values() if node.type == "reaction"]
        else:
            return [node for node in self.dicOfNode.values() if node.type == "reaction"]

    def getGenes(self, get_id=True):
        """
        Get all the genes from a padmet instance.

        Parameters
        ----------
        self: padmet
            padmet instance
        get_id: bool
            True for the ID and False for the node corresponding to the genes

        Returns
        -------
        list
            list of genes IDs or nodes
        """
        if get_id is True:
            return [node.id for node in self.dicOfNode.values() if node.type == "gene"]
        else:
            return [node for node in self.dicOfNode.values() if node.type == "gene"]

    def getPathways(self, get_id=True):
        """
        Get all the pathways from a padmet instance.

        Parameters
        ----------
        self: padmet
            padmet instance
        get_id: bool
            True for the ID and False for the node corresponding to the pathways

        Returns
        -------
        list
            list of pathways IDs or nodes
        """
        total_pwy_id = set()

        all_rxns = [node for node in self.dicOfNode.values() if node.type == "reaction"]

        for rxn_node in all_rxns:
            pathways_ids = set([rlt.id_out for rlt in self.dicOfRelationIn[rxn_node.id] if rlt.type == "is_in_pathway"])
            total_pwy_id.update(pathways_ids)

        if get_id is True:
            all_pwys = [node_id for (node_id, node) in self.dicOfNode.items() if node_id in total_pwy_id]
        else:
            all_pwys = [node for (node_id, node) in self.dicOfNode.items() if node_id in total_pwy_id]

        return all_pwys

    def getPathwaysReactions(self):
        """
        Get all the pathways from a padmet instance with their reactions.

        Parameters
        ----------
        self: padmet
            padmet instance

        Returns
        -------
        dict
            dictionary of pathways IDs as key with list of reactions IDs as value
        """
        # Iterate through node of type "pathway", then for each node found the relations linking the pathway to the reaction it contains.
        return {node.id: [
                            rlt.id_in for rlt in self.dicOfRelationOut[node.id]
                            if rlt.type == "is_in_pathway" and self.dicOfNode[rlt.id_in].type == 'reaction'
                            ]
                        for node in self.dicOfNode.values() if node.type == "pathway"}

    def loadGraph(self, padmet_file):
        """
        Allow to recover all the informations of the padmet file.
        The section Data Base informations corresponds to the self.info
        The section Nodes corresponds to the data of each nodes in self.dicOfNode, sep ="\t"
        The section Relations corresponds to the data of each relations in self.dicOfRelationIn/Out, sep ="\t"

        Parameters
        ----------
        padmet_file: str
            the pathname of the padmet file to load.
        """
        with open(padmet_file, "r", encoding="utf8") as f:
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
        # Index where start policy
        policy_index = padmet_in_array.index("Policy")
        # Index where start Nodes
        node_index = padmet_in_array.index("Nodes")
        # Index where start Relations
        relation_index = padmet_in_array.index("Relations")

        #
        if info_index is not None:
            info_section = padmet_in_array[info_index + 1 : policy_index]
            #
            for line in info_section:
                if "\t" not in line:
                    line = line.replace(":", "")
                    self.info[line] = {}
                    current_data = line
                else:
                    line = line.replace("\t", "").split(":")
                    self.info[current_data][line[0]] = line[1]

        #
        policy_in_array = [
            line.split("\t") for line in padmet_in_array[policy_index + 1 : node_index]
        ]
        self.setPolicy(policy_in_array)

        # generator of data to input in Node()
        node_section = (
            line.split("\t")
            for line in padmet_in_array[node_index + 1 : relation_index]
        )
        #
        for data in node_section:
            # Create a new Node object (cf node.py)
            node_type = data[0]
            node_id = data[1]
            node_misc = {}
            if len(data) > 2:
                i = 2
                # add the misc data in a dict
                while (i + 1) < len(data):
                    # in case diff information for the same key
                    try:
                        node_misc[data[i]].append(data[i + 1])
                    except KeyError:
                        node_misc[data[i]] = [data[i + 1]]
                    i += 2
            node = Node(node_type, node_id, node_misc)
            # add the node into the dictionnay
            self.dicOfNode[node.id] = node

        # generator of data to input in Relation()
        relation_section = (
            line.split("\t") for line in padmet_in_array[relation_index + 1 :]
        )
        # instantiate a new Relation object (cf relation.py)
        for data in relation_section:
            rlt_id_in = data[0]
            rlt_type = data[1]
            rlt_id_out = data[2]
            rlt_misc = {}
            if len(data) > 3:
                i = 3
                while (i + 1) < len(data):
                    # in case diff information for the same key
                    try:
                        rlt_misc[data[i]].append(data[i + 1])
                    except KeyError:
                        rlt_misc[data[i]] = [data[i + 1]]
                    i += 2

            relation = Relation(rlt_id_in, rlt_type, rlt_id_out, rlt_misc)
            try:
                self.dicOfRelationIn[rlt_id_in].append(relation)
            except KeyError:
                self.dicOfRelationIn[rlt_id_in] = [relation]
            try:
                self.dicOfRelationOut[rlt_id_out].append(relation)
            except KeyError:
                self.dicOfRelationOut[rlt_id_out] = [relation]

    def updateFromSbml(self, sbml_file, verbose=False):
        """
        Initialize a padmetRef from sbml. Copy all species, convert id with sbmlPlugin
        stock name in COMMON NAME. Copy all reactions,
        convert id with sbmlPlugin, stock name in common name, stock compart and stoichio data relative
        to reactants and products in the misc of consumes/produces relations

        Parameters
        ----------
        sbml_file: str
            pathname of the sbml file
        verbose: bool
            if True print supp info
        """
        if not os.path.exists(sbml_file):
            raise FileNotFoundError("No SBML file accessible at " + sbml_file)

        # using libSbml to read sbml_file
        if verbose:
            print("loading sbml file: %s" % sbml_file)
        reader = libsbml.SBMLReader()
        document = reader.readSBML(sbml_file)
        for i in range(document.getNumErrors()):
            print(document.getError(i).getMessage())

        model = document.getModel()
        listOfSpecies = model.getListOfSpecies()
        listOfReactions = model.getListOfReactions()
        nbReactions = len(listOfReactions)
        nbSpecies = len(listOfSpecies)
        if verbose:
            print("nb species: %s" % nbSpecies)
            print("nb reactions: %s" % nbReactions)

        if verbose:
            print("creating species")
        for specie in listOfSpecies:
            specie_id_encoded = specie.getId()
            specie_id = sbmlPlugin.convert_from_coded_id(specie_id_encoded)[0]
            if verbose:
                print("specie: %s, uncoded: %s" % (specie_id_encoded, specie_id))
            try:
                self.dicOfNode[specie_id]
                if verbose:
                    print("already in padmetRef")
            except KeyError:
                specie_name = specie.getName()
                if specie_name:
                    specie_node = Node(
                        "compound", specie_id, {"COMMON-NAME": [specie_name]}
                    )
                else:
                    specie_node = Node("compound", specie_id)
                self.dicOfNode[specie_id] = specie_node

        if verbose:
            print("creating reactions")
        for reaction in listOfReactions:
            reaction_id_encoded = reaction.getId()
            reaction_id = sbmlPlugin.convert_from_coded_id(reaction_id_encoded)[0]
            if verbose:
                print("reaction: %s, uncoded: %s" % (reaction_id_encoded, reaction_id))
            try:
                self.dicOfNode[reaction_id]
                if verbose:
                    print("already in padmetRef")
            except KeyError:
                reaction_name = reaction.getName()
                if reaction.getReversible():
                    reaction_dir = "REVERSIBLE"
                else:
                    reaction_dir = "LEFT-TO-RIGHT"
                if reaction_name:
                    reaction_node = Node(
                        "reaction",
                        reaction_id,
                        {"COMMON-NAME": [reaction_name], "DIRECTION": [reaction_dir]},
                    )
                else:
                    reaction_node = Node(
                        "reaction", reaction_id, {"DIRECTION": [reaction_dir]}
                    )
                self.dicOfNode[reaction_id] = reaction_node
                reactants = reaction.getListOfReactants()
                for reactant in reactants:
                    reactant_id, x, reactant_compart = sbmlPlugin.convert_from_coded_id(
                        reactant.getSpecies()
                    )
                    if reactant_compart is None:
                        if verbose:
                            print("%s has no compart, set to 'c'" % reactant)
                        reactant_compart = "c"
                    reactant_stoich = reactant.getStoichiometry()
                    consumes_rlt = Relation(
                        reaction_id,
                        "consumes",
                        reactant_id,
                        {
                            "STOICHIOMETRY": [reactant_stoich],
                            "COMPARTMENT": [reactant_compart],
                        },
                    )
                    self._addRelation(consumes_rlt)
                products = reaction.getListOfProducts()
                for product in products:
                    product_id, x, product_compart = sbmlPlugin.convert_from_coded_id(
                        product.getSpecies()
                    )
                    if product_compart is None:
                        if verbose:
                            print("%s has no compart, set to 'c'" % product)
                        product_compart = "c"
                    product_stoich = product.getStoichiometry()
                    produces_rlt = Relation(
                        reaction_id,
                        "produces",
                        product_id,
                        {
                            "STOICHIOMETRY": [product_stoich],
                            "COMPARTMENT": [product_compart],
                        },
                    )
                    self._addRelation(produces_rlt)

    def generateFile(self, output):
        """
        Allow to create a padmet file to stock all the data.
        @param output: pathname of the padmet file to create

        Parameters
        ----------
        output: str
            path to output file
        """
        # Order the dictionary of node by unique id and the tuple of relation
        # by the node in id.
        dicKeys = list(self.dicOfNode.keys())
        dicKeys.sort()
        orderedRelation = tuple(
            sorted(self.getAllRelation(), key=lambda x: x.id_in, reverse=False)
        )
        with open(output, "w", encoding="utf8") as f:
            if len(self.info) != 0:
                f.write("Data Base informations\n")
                for k, data in self.info.items():
                    f.write(k + ":\n")
                    for k, v in data.items():
                        f.write("\t" + k + ":" + v + "\n")
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
                line = node.toString() + "\n"
                f.write(line)
            f.write("\n")
            # It writes the relations of the TGDBP file.
            f.write("Relations\n")
            f.write("\n")
            for rlt in orderedRelation:
                line = rlt.toString() + "\n"
                f.write(line)

    # ==============================================================================
    # For Nodes
    # ==============================================================================

    def _addNode(self, node):
        """
        Allows to add a node, only if the id is not already used.

        Parameters
        ----------
        node:  padmet.classes.Node
            the node to add

        Returns
        -------
        bool
            True if added, False if no.        
        """
        if node.id not in list(self.dicOfNode.keys()):
            self.dicOfNode[node.id] = node
            return True
        else:
            return False

    def delNode(self, node_id):
        """
        Allows to delete a node, the relations associated to the node, and for 
        some relations, delete the associated node. 
        For relations where the node to del is 'in': 
            if rlt type in ['has_xref','has_name','has_suppData']: delNode out
        For relations where the node to del is 'out':
            if rlt type in ['consumes','produces']

        Parameters
        ----------
        node_id: str
            id of node to delete

        Returns
        -------
        bool
            True if node successfully deleted, False if node not in dicOfNode
        """
        # Delete the node from dicOfNode.
        try:
            self.dicOfNode.pop(node_id)
        except KeyError:
            print("The id %s doesnt exist. Unable to delete" % node_id)
            return False

        # If exist delete the relations 'in' and 'out'
        try:
            # Recover the relations where the node is "in".
            relationsIn = [rlt for rlt in self.dicOfRelationIn.get(node_id, None)]
            for rltIn in relationsIn:
                self._delRelation(rltIn)
                if rltIn.type in [
                    "has_xref",
                    "has_name",
                    "has_suppData",
                    "has_reconstructionData",
                ]:
                    self.delNode(rltIn.id_out)
        except TypeError:
            pass
        try:
            # Recover the relations where the node is "out"
            relationsOut = [rlt for rlt in self.dicOfRelationOut.get(node_id, None)]
            for rltOut in relationsOut:
                self._delRelation(rltOut)
        except TypeError:
            pass
        return True

    # ==============================================================================
    # For Relations:
    # ==============================================================================

    def _delRelation(self, relation):
        """
        Delete a relation from dicOfRelationIn and out

        Parameters
        ----------
        relation: padmet.classes.Relation
            the relation to delete

        Returns
        -------
        bool
            True if succesfully deleted
        """
        delete = False
        idIn = relation.id_in
        idOut = relation.id_out
        try:
            for rlt in self.dicOfRelationIn[idIn]:
                if relation.compare(rlt):
                    self.dicOfRelationIn[idIn].remove(rlt)
                    if len(self.dicOfRelationIn[idIn]) == 0:
                        self.dicOfRelationIn.pop(idIn)
                    delete = True
                    break
        except KeyError:
            pass
        try:
            for rlt in self.dicOfRelationOut[idOut]:
                if relation.compare(rlt):
                    self.dicOfRelationOut[idOut].remove(rlt)
                    if len(self.dicOfRelationOut[idOut]) == 0:
                        self.dicOfRelationOut.pop(idOut)
                    delete = True
                    break
        except KeyError:
            pass
        if not delete:
            print("Unable to delete this relation, doesn't exist.")
        return delete

    def _addRelation(self, relation):
        """
        AddRelation() allows to add a relation if not already in allRelations.

        Parameters
        ----------
        relation: padmet.classes.Relation
            the relation to add

        Returns
        -------
        bool
            true if relation was successfully added
        """
        idIn = relation.id_in
        idOut = relation.id_out
        rlt_to_compare = self.dicOfRelationIn.get(idIn, [])
        rlt_to_compare += self.dicOfRelationOut.get(idOut, [])
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

    # ==============================================================================
    # manipulating de novo node:
    # ==============================================================================

    def createNode(self, _type, _id, dicOfMisc={}, listOfRelation=None):
        """
        Creation of new node to add in the network.

        Parameters
        ----------
        _type: str
            type of node (gene, reaction...)
        _id: str
            id of the node
        dicOfMisc: dict
            dictionnary of miscellaneous data
        listOfRelation: list or None
            list of relation

        Returns
        -------
        padmet.classes.Node
            new_node
        """
        # add the new node in the tgdbp
        if _id in list(self.dicOfNode.keys()):
            raise ValueError("%s is already used, unable to create Node" % (_id))
        new_node = Node(_type, _id, dicOfMisc)
        self.dicOfNode[_id] = new_node
        if listOfRelation is not None:
            for rlt in listOfRelation:
                self._addRelation(rlt)

        return new_node
