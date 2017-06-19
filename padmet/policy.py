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
Define a policy in padmet object.
"""

class Policy:
    """
    A Policy define the types of relations, nodes in a network.
    A policy contains 3 attributs:
        policy_in_array: Is a list of list of arcs (eg: [['reaction','consumes','compounds'],['reaction','produces','compounds']])
        class_of_node: Is a set of all the type of nodes represented in the network (eg: set(reaction, compound))
        type_of_arc: Is a dictionnary of all the types of arcs represented in the network (eg: {reaction:[consumes,compounds]})
    """

    def __init__(self,policy_in_array = None):
        """
        @param policy_in_array: Is a list of list of arcs (eg: [['reaction','consumes','compounds'],['reaction','produces','compounds']])
        @type policy_in_array: list
        @return: _
        @rtype: None
        """
        if policy_in_array is not None:
            self.setPolicyInArray(policy_in_array)
        else:
            #set an empty policy
            self.policy_in_array = []
            self.class_of_node = set()
            self.type_of_arc = {}
            
    def setPolicyInArray(self,policy_in_array):
        """
        From policy_in_array, set class_of_node and type_of_arc
        @param policy_in_array: Is a list of list of arcs (eg: [['reaction','consumes','compounds'],['reaction','produces','compounds']])
        @type policy_in_array: list
        @return: _
        @rtype: None
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
        return policy_in_array
        @return: self.policy_in_array
        @rtype: list
        """
        return self.policy_in_array
            
    def _setClassOfNode(self):
        """
        From self.policy_in_array set class_of_node
        @return: _
        @rtype: None
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
        return class_of_node
        @return: self.class_of_node
        @rtype: set
        """
        return self.class_of_node

    def _setTypeOfArc(self):
        """
        From self.policy_in_array and self.class_of_node set type_of_arc
        @return: _
        @rtype: None
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
        return class_of_node
        @return: self.class_of_node
        @rtype: set
        """
        return(self.type_of_arc)

