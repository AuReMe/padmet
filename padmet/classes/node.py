#!/usr/bin/env python
# -*- coding: utf-8 -*-
#pylint: disable=too-few-public-methods
class Node:
    """
    A Node represent an element in a metabolic network
    
    e.g: compound, reaction.
    """
    def __init__(self, _type, _id, misc=None):
        """
        Parameters
        ----------
        _type: str
            The type of the node ('reaction','pathway')
        _id: str
            the identifier of the node ('rxn-45)
        misc: dict
            A dictionary of miscellaneous data ({'DIRECTION':[REVERSIBLE]})
            (the default value is None)
        """
        self.type = _type
        self.id = _id
        #if misc is None so misc = {}
        self.misc = misc or {}

    def toString(self):
        """
        This function is used to stock the information relative to the node
        in a padmet file.

        Returns
        -------
        str
            string with all data sep by tab' ex: reaction\tRXN0..
        """
        sep = "\t"
        line = sep.join([self.type, self.id])
        if len(self.misc) != 0:
            for k, n in self.misc.items():
                if len(n) == 1:
                    line += sep + sep.join([str(k), str(n[0])])
                else:
                    for i in range(len(n)):
                        line += sep + sep.join([str(k), str(n[i])])
        return line
