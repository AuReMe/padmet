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
Define the class Relation used in padmet.
"""

class Relation:
    """
    A Relation represent a link between two elements (node) in a metabolic network
    e.g: RXN-1 consumes CPD-1
    A Relation contains 4 attributes:
        _type: The type of the relation (e.g: 'consumes' or 'produces')
        id_in: the identifier of the node corresponding to the subject of the relation (e.g: 'RXN-1)
        id_out: the identifier of the node corresponding to the object of the relation (e.g: 'CPD-1)
        misc: A dictionary of miscellaneous data, k = tag of the data, v = list of values
        (e.g: {'STOICHIOMETRY':[1.0]})
    """
    def __init__(self, id_in, _type, id_out, misc=None):
        """
        Parameters
        ----------
        id_in: str
            the identifier of the node corresponding to the subject of the relation ('RXN-1)
        _type: str
            The type of the relation (e.g: 'consumes' or 'produces')
        id_out: str
            the identifier of the node corresponding to the object of the relation ('CPD-1)
        _misc: dict
            A dictionary of miscellaneous data (e.g: {'STOICHIOMETRY':[1.0]})
        """
        self.id_in = id_in
        self.type = _type
        self.id_out = id_out
        self.misc = misc or {}

    def toString(self):
        """
        This function is used to stock the information relative to the node
        in a padmet file.

        Returns
        -------
        str
            string with all data sep by tab' ex: RXN0\tconsumes\tCPD-A..
        """
        sep = "\t"
        line = sep.join([self.id_in, self.type, self.id_out])
        if len(self.misc) != 0:
            for k, n in self.misc.iteritems():
                if len(n) == 1:
                    line += sep + sep.join([str(k), str(n[0])])
                else:
                    for i in range(len(n)):
                        line += sep + sep.join([str(k), str(n[i])])
        return line

    def compare(self, relation):
        """
        compare 2 relations. First check if ids and type are the same, then check
        the misc dictionary.
        Parameters
        ----------
        relation: padmet.Relation
            the relation to compare

        Returns
        -------
        bool
            Return True if relation are the same, False if not        
        """
        if set([relation.id_in, relation.type, relation.id_out]) == set([self.id_in, self.type, self.id_out]):
            if relation.misc == self.misc:
                return True
        return False
                        