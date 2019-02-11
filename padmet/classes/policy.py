#!/usr/bin/env python
# -*- coding: utf-8 -*-
class Policy:
    """
    A Policy define the types of relations and nodes of a network.

    A policy contains 3 attributes:

        policy_in_array: Is a list of list of relations

        e.g: [['reaction','consumes','compounds'],['reaction','produces','compounds']]

        class_of_node: Is a set of all the type of nodes represented in the network

        e.g: set(['reaction', 'compound'])

        type_of_arc: Is a dictionary of all the types of arcs represented in the network

        (e.g: {'reaction':['consumes','compounds']})
    """

    def __init__(self, policy_in_array=None):
        """
        Parameters
        ----------
        policy_in_array: list
            Is a list of list of relations

            e.g: [['reaction','consumes','compounds'],['reaction','produces','compounds']]

            (the default value is None)
        """
        if policy_in_array is not None:
            self.setPolicyInArray(policy_in_array)
        else:
            #set an empty policy
            self.policy_in_array = []
            self.class_of_node = set()
            self.type_of_arc = {}

    def setPolicyInArray(self, policy_in_array):
        """
        From policy_in_array, set class_of_node and type_of_arc
        
        Parameters
        ----------
        policy_in_array: list
            Is a list of list of arcs. e.g: [['reaction','consumes','compounds'],['reaction','produces','compounds']]
        """
        for relation in policy_in_array:
            if len(relation) < 3:
                raise ValueError("Array given to set the PolicyInArray is uncorrect: \
                " + str(relation))
        self.policy_in_array = policy_in_array
        self._setClassOfNode()
        self._setTypeOfArc()

    def getPolicyInArray(self):
        """
        Returns
        -------
        list
            return policy_in_array
        """
        return self.policy_in_array

    def _setClassOfNode(self):
        """
        From self.policy_in_array set class_of_node
        """
        self.class_of_node = set()
        if len(self.policy_in_array) != 0:
            for line in self.policy_in_array:
                self.class_of_node.add(line[0])
                self.class_of_node.add(line[2])
        else:
            raise ValueError("PolicyInArray is not defined")

    def getClassOfNode(self):
        """
        Returns
        -------
        set
            return class_of_node
        """
        return self.class_of_node

    def _setTypeOfArc(self):
        """
        From self.policy_in_array and self.class_of_node set type_of_arc
        """
        self.type_of_arc = {}
        if len(self.class_of_node) != 0:
            for _class in self.class_of_node:
                self.type_of_arc[_class] = []

            for line in self.policy_in_array:
                arc = line[1:]
                self.type_of_arc[line[0]].append(arc)
        else:
            raise ValueError("PolicyInArray and/or classOfNode are not defined")

    def getTypeOfArc(self):
        """
        Returns
        -------
        dict
            return type_of_arc
        """
        return self.type_of_arc

