========
Tutorial
========

PADMet format
-------------

Padmet is an object representing the metabolic network of a species (organism) based on a reference database.

This format contains Node, Policy and Relation:

- a Node is an Object that contains information about an element of the network (can be a pathway, reaction...).

- a Policy defines the way Node and Relation are associated.

- a Relation defines how two nodes are connected. In a relation there is a node "in" and a node "out". (reactionX (**node in**) consumes metaboliteX (**node out**))

A Padmet object contains 3 attributes:

- dicOfNode: a dictionary of node: key=Node's unique id / value = <Node>

- dicOfRelationIn: a dictionary of relation with: key= nodeIN id / value = list of <relation>

- dicOfRelationOut: a dictionary of relation with: key= nodeOut id / value = list of <relation>

Node
----

A Node represent an element in a metabolic network.

A Node has a type ('reaction','pathway'), an ID and a dictionary of miscellaneous data. Node can be find in the dicOfNode, this dicitonary contains Node IDs as key and Node object as values.

As an example we use the padmet_1.padmet from `test data <https://github.com/AuReMe/padmet/tree/master/tests/test_data/padmet>`__:

.. code:: python

    >>> from padmet.classes import PadmetSpec

    # Creation of a padmet object from the padmet file.
    >>> padmet_object = PadmetSpec('test_data/padmet/padmet_1.padmet')

    # Selection of the Node corresponding to the reaction 'ENOYL-COA-HYDRAT-RXN'
    >>> enoyl_coa_node = padmet_object.dicOfNode['ENOYL-COA-HYDRAT-RXN']

    >>> enoyl_coa_node.id 
    'ENOYL-COA-HYDRAT-RXN'

    >>> enoyl_coa_node.type
    'reaction'

    >>> enoyl_coa_node.misc
    {'EC-NUMBER': ['EC-4.2.1.17'], 'DIRECTION': ['LEFT-TO-RIGHT']}

Relation
--------

A Relation represent a link between two elements (node) in a metabolic network.

A Relation contains 4 attributes:

- type: The type of the relation (e.g: 'consumes' or 'produces')
- id_in: the identifier of the node corresponding to the subject of the relation (e.g: 'RXN-1')
- id_out: the identifier of the node corresponding to the object of the relation (e.g: 'CPD-1')
- misc: A dictionary of miscellaneous data, k = tag of the data, v = list of values (e.g: {'STOICHIOMETRY':[1.0]})

We will use the same example as the Node section:

.. code:: python

    >>> from padmet.classes import PadmetSpec

    # Creation of a padmet object from the padmet file.
    >>> padmet_object = PadmetSpec('test_data/padmet/padmet_1.padmet')

    # Selection of the Node corresponding to the reaction 'ENOYL-COA-HYDRAT-RXN'
    >>> enoyl_coa_relation = padmet_object.dicOfRelationIn['ENOYL-COA-HYDRAT-RXN']

    >>> len(enoyl_coa_relation)
    6
    # Node 'ENOYL-COA-HYDRAT-RXN' is an ID in for 6 Relations

    # What is the type of the first relation when this Node is In
    >>> padmet_object.dicOfRelationIn['ENOYL-COA-HYDRAT-RXN'][0].type
    'is_in_pathway'

    >>> padmet_object.dicOfRelationIn['ENOYL-COA-HYDRAT-RXN'][0].id_in
    'ENOYL-COA-HYDRAT-RXN'

    >>> padmet_object.dicOfRelationIn['ENOYL-COA-HYDRAT-RXN'][0].id_out
    'FAO-PWY'

    >>> padmet_object.dicOfRelationIn['ENOYL-COA-HYDRAT-RXN'][0].toString()
    'ENOYL-COA-HYDRAT-RXN\tis_in_pathway\tFAO-PWY'

Create padmet
-------------

Padmet files can be created using the padmet ``padmet pgdb_to_padmet`` command. This command takes as input the result of Pathway Tools software. You can get these files by using the `mpwt <https://github.com/AuReMe/mpwt>`__ package.

    .. code::

        usage:
            padmet pgdb_to_padmet --pgdb=DIR --output=FILE [--version=V] [--db=ID] [--padmetRef=FILE] [--source=STR] [-v] [--enhance]
            padmet pgdb_to_padmet --pgdb=DIR --output=FILE --extract-gene [--no-orphan] [--keep-self-rxn] [--version=V] [--db=ID] [--padmetRef=FILE] [--source=STR] [-v] [--enhance]

        options:
            -h --help     Show help.
            --version=V    Xcyc version [default: N.A].
            --db=ID    Biocyc database corresponding to the pgdb (metacyc, ecocyc, ...) [default: N.A].
            --output=FILE    padmet file corresponding to the DB.
            --pgdb=DIR    directory containg all the .dat files of metacyc (data).
            --padmetRef=FILE    padmet of reference.
            --source=STR    Tag associated to the source of the reactions, used to ensure traceability [default: GENOME].
            --enhance    use the metabolic-reactions.xml file to enhance the database.
            --extract-gene    extract genes from genes_file (use if its a specie's pgdb, if metacyc, do not use).
            --no-orhpan    remove reactions without gene associaiton (use if its a specie's pgdb, if metacyc, do not use).
            --keep-self-rxn    remove reactions with no reactants (use if its a specie's pgdb, if metacyc, do not use).
            -v   print info.


The input folder (--pgdb) should be like:

    ::

        input_folder
        ├── classes.dat
        ├── compound-links.dat
        ├── compounds.dat
        ├── dnabindsites.dat
        ├── enzrxns.dat
        ├── gene-links.dat
        ├── genes.dat
        ├── pathway-links.dat
        ├── pathways.dat
        ├── promoters.dat
        ├── protein-features.dat
        ├── protein-links.dat
        ├── proteins.dat
        ├── protligandcplxes.dat
        ├── pubs.dat
        ├── reaction-links.dat
        ├── reactions.dat
        ├── regulation.dat
        ├── regulons.dat
        ├── rnas.dat
        ├── species.dat
        ├── terminators.dat
        └── transunits.dat

This command will create a padmet file.

Explore padmet
--------------

List of reactions, compounds, pathways and genes of a metabolic network can be extracted using getter on the padmet:

.. code:: python

    >>> from padmet.classes import PadmetSpec

    # Creation of a padmet object from the padmet file.
    >>> padmet_object = PadmetSpec('test_data/padmet/padmet_1.padmet')

    # Get the list of reaction IDs in the padmet.
    >>> list_reaction_ids = padmet_object.getReactions()

    # Get the list of compounds in the padmet.
    >>> list_compound_ids = padmet_object.getCompounds()

    # Get the list of genes in the padmet.
    >>> list_gene_ids = padmet_object.getGenes()

    # Get the list of pathways in the padmet.
    >>> list_pathway_ids = padmet_object.getPathways()

It is also possible to get a list of nodes:

.. code:: python

    # Get the list of reaction IDs in the padmet.
    >>> list_reaction_nodes = padmet_object.getReactions(get_id=False)