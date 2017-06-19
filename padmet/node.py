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
Defines the class Node used in padmet.
"""

class Node:
    """
    A Node represent an element in a metabolic network (eg: compound, reaction)
    A Node contains 3 attributs:
        _type: The type of the node (eg: 'reaction' or 'pathway')
        _id: the identifier of the node (eg: 'rxn-45)
        _misc: A dictionnary of miscellaneous data, k = tag of the data, v = list of values
        (eg: {'DIRECTION':[REVERSIBLE]})
    """
    def __init__(self, _type, _id, misc = None):
        """
        @param _type: The type of the node ('reaction','pathway')
        @param _id: the identifier of the node ('rxn-45)
        @param _misc: A dictionnary of miscellaneous data ({'DIRECTION':[REVERSIBLE]})
        @type _type, _id: str
        @type misc: dict
        @return: _
        @rtype: None
        """
        self.type = _type
        self.id = _id
        #if misc is None so misc = {}
        self.misc = misc or {}
    
    def toString(self):
        """
        This function is used to stock the information relative to the node
        in a padmet file.
        @return: string with all data sep by tab' ex: reaction\tRXN0..
        @rtype: str
        """
        sep = "\t"
        line = sep.join([self.type, self.id])
        if len(self.misc) != 0:
            for k,n in self.misc.iteritems():
                if len(n) == 1:
                    line += sep + sep.join([str(k),str(n[0])])
                else:
                    for i in range(len(n)):
                        line += sep + sep.join([str(k),str(n[i])])
        return line