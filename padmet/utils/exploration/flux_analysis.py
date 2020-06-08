    # -*- coding: utf-8 -*-
"""
Description:
    1./ Run flux balance analyse with cobra package on an already defined reaction.
    Need to set in the sbml the value 'objective_coefficient' to 1.
    If the reaction is reachable by flux: return the flux value and the flux value
    for each reactant of the reaction.
    If not: only return the flux value for each reactant of the reaction.
    If a reactant has a flux of '0' this means that it is not reachable by flux
    (and maybe topologically). To unblock the reaction it is required to fix the
    metabolic network by adding/removing reactions until all reactant are reachable.
    
    2./If seeds and targets given as sbml files with only compounds.
    Will also try to use the Menetools library to make a topologicall analysis.
    Topological reachabylity of the targets compounds from the seeds compounds.
    
    3./ If --all_species: will test flux reachability of all the compounds in the
    metabolic network (may take several minutes)

::

    usage:
        padmet flux_analysis --sbml=FILE
        padmet flux_analysis --sbml=FILE --seeds=FILE --targets=FILE [--all_species]
        padmet flux_analysis --sbml=FILE --all_species

    option:
        -h --help    Show help.
        --sbml=FILE    pathname to the sbml file to test for fba and fva.
        --seeds=FILE    pathname to the sbml file containing the seeds (medium).
        --targets=FILE    pathname to the sbml file containing the targets.
        --all_species    allow to make FBA on all the metabolites of the given model.
"""
import docopt
import os

from padmet.utils.sbmlPlugin import convert_from_coded_id
from cobra import Reaction
from cobra import flux_analysis as cobra_flux_analysis
from cobra.io.sbml import read_sbml_model


def command_help():
    """
    Show help for analysis command.
    """
    print(docopt.docopt(__doc__))


def flux_analysis_cli(command_args):
    args = docopt.docopt(__doc__, argv=command_args)
    sbml_file = args["--sbml"]
    seeds_file = args["--seeds"]
    targets_file = args["--targets"]
    all_species = args["--all_species"]
    flux_analysis(sbml_file, seeds_file, targets_file, all_species)


def flux_analysis(sbml_file, seeds_file = None, targets_file = None, all_species = False):
    """
    1./ Run flux balance analyse with cobra package on an already defined reaction.
    Need to set in the sbml the value 'objective_coefficient' to 1.
    If the reaction is reachable by flux: return the flux value and the flux value
    for each reactant of the reaction.
    If not: only return the flux value for each reactant of the reaction.
    If a reactant has a flux of '0' this means that it is not reachable by flux
    (and maybe topologically). To unblock the reaction it is required to fix the
    metabolic network by adding/removing reactions until all reactant are reachable.
    
    2./If seeds and targets given as sbml files with only compounds.
    Will also try to use the Menetools library to make a topologicall analysis.
    Topological reachabylity of the targets compounds from the seeds compounds.
    
    3./ If --all_species: will test flux reachability of all the compounds in the
    metabolic network (may take several minutes)

    Parameters
    ----------
    sbml_file: str
        path to sbml file to analyse
    seeds_file: str
        path to sbml file with only compounds representing the seeds/growth medium
    targets_file: str
        path to sbml file with only compounds representing the targets to reach
    all_species: bool
        if True will try to create obj function for each compound and return which are reachable by flux.
        
    """
    if targets_file:
        if not os.path.exists(targets_file):
            raise FileNotFoundError("No target SBML file accessible at " + targets_file)
        targets = read_sbml_model(targets_file).metabolites

    if seeds_file:
        if not os.path.exists(seeds_file):
            raise FileNotFoundError("No seeds SBML file accessible at " + seeds_file)

    if not os.path.exists(sbml_file):
        raise FileNotFoundError("No target SBML file accessible at " + sbml_file)

    model=read_sbml_model(sbml_file)
    
    #nb metabolites
    real_metabolites = set([i.id.replace("_"+i.compartment,"") for i in model.metabolites])
    rxn_with_ga = [i for i in model.reactions if i.gene_reaction_rule]
    print("#############")
    print("Model summary")
    print("Number of compounds: %s" %len(real_metabolites))
    print("Number of reactions: %s" %len(model.reactions))
    print("Number of genes: %s" %len(model.genes))
    print("Ratio rxn with genes/rxns: %s%%" %(100*len(rxn_with_ga)/len(model.reactions)))
    
    # Launch a topoligical analysis if menetools is installed.
    if seeds_file and targets_file:
        print("#############")
        print("Analyzing targets")
        print("#Topological analysis")
        try:
            from menetools import run_menecheck

            menetools_result = run_menecheck(draft_sbml=sbml_file, seeds_sbml=seeds_file, targets_sbml=targets_file)
            print("Number of targets: %s" %(len(targets)))
            print("Unproductible targets: " + ",".join(menetools_result[0]))
            print("Productible targets: " + ",".join(menetools_result[1]))
        except ImportError:
            print("Menetools is not installed. Can't run topological analysis.")
        print("#Flux Balance Analysis")
        fba_on_targets(targets, model)
    if all_species:
        targets = model.metabolites
        print("#Flux Balance Analysis on all model metabolites (long process...)")
        fba_on_targets(targets, model)
        return
    try:
        biomassrxn = [rxn for rxn in model.reactions if rxn.objective_coefficient == 1.0][0]
        biomassname = biomassrxn.id
    except IndexError:
        print("Need to set OBJECTIVE COEFFICIENT to '1.0' for the reaction to test")
        exit()

    print("#############")
    print("Computing optimization")
    solution= model.optimize()
    print("Testing reaction %s" %biomassname)
    print("Growth rate: %s" %solution.objective_value)
    print("Status: %s" %solution.status)
    model.summary()
    if (solution.objective_value > 1e-5):
        blocked = cobra_flux_analysis.find_blocked_reactions(model, model.reactions)
        essRxns = cobra_flux_analysis.find_essential_reactions(model)
        essGenes = cobra_flux_analysis.find_essential_genes(model)

        print('FVA analysis:')
        print('\tBlocked reactions: %s' %len(blocked))
        print('\tEssential reactions: %s' %len(essRxns))
        [print(rxn.id) for rxn in essRxns]
        print('\tEssential genes: %s' %len(essGenes))
    
    #get biomass rxn reactants
    bms_reactants = dict([(k,v) for k,v in list(biomassrxn.metabolites.items()) if v < 0])
    bms_products = dict([(k,v) for k,v in list(biomassrxn.metabolites.items()) if v > 0])
    dict_output = {"positive":{},"negative":{}}
    #for each metabolite in reactant, create a biomass rxn with only this metabolite in reactants
    biomassrxn.objective_coefficient = 0.0
    for reactant, stoich in list(bms_reactants.items()):
        test_rxn = Reaction("test_rxn")
        test_rxn.lower_bound = 0
        test_rxn.upper_bound = 1000
        metabolitedict = dict(bms_products)
        metabolitedict.update({reactant:stoich})

        model.add_reactions([test_rxn])
        test_rxn.add_metabolites(metabolitedict)
        test_rxn.objective_coefficient = 1.0

        solution = model.optimize()
        if (solution.objective_value > 1e-5):
            dict_output["positive"][reactant] = solution.objective_value
        else:
            dict_output["negative"][reactant] = solution.objective_value
        model.remove_reactions([test_rxn])
    print("%s/%s compounds with positive flux" %(len(list(dict_output["positive"].keys())), len(bms_reactants)))
    print("%s/%s compounds without flux" %(len(list(dict_output["negative"].keys())), len(bms_reactants)))

    for k,v in list(dict_output["positive"].items()):
        print("%s // %s %s positive" %(k, convert_from_coded_id(k.id)[0]+"_"+convert_from_coded_id(k.id)[2], v))
    for k,v in list(dict_output["negative"].items()):
        print("%s // %s %s NULL" %(k, convert_from_coded_id(k.id)[0]+"_"+convert_from_coded_id(k.id)[2], v))


def fba_on_targets(allspecies, model):
    """
    for each specie in allspecies, create an objective function with the current species as only product
    and try to optimze the model and get flux.

    Parameters
    ----------
    allSpecies: list
        list of species ids to test
    model: cobra.model
        Cobra model from a sbml file
    """
    #dict_output = {"positive":{},"negative":{}}
    for species in allspecies:
        #lets create a copy of the initial model
        model2 = model.copy()
        #remove all obj coef
        for rxn in model2.reactions:
            if rxn.objective_coefficient == 1.0:
                rxn.objective_coefficient = 0.0
        #Create a new reaction that consume the given species
        FBA_rxn = Reaction("FBA_TEST")
        FBA_rxn.lower_bound = 0
        FBA_rxn.upper_bound = 1000
        model2.add_reactions([FBA_rxn])
        FBA_rxn.objective_coefficient = 1.0
        metabolitedict = {}
        metabolitedict[species]=-1.0
        FBA_rxn.add_metabolites(metabolitedict)
                
        solution = model2.optimize()
        if (solution.objective_value > 1e-5):
            print("%s // %s %s positive" %(species, convert_from_coded_id(species.id)[0]+"_"+convert_from_coded_id(species.id)[2], solution.objective_value))
        else:
            print("%s // %s %s NULL" %(species, convert_from_coded_id(species.id)[0]+"_"+convert_from_coded_id(species.id)[2], solution.objective_value))
