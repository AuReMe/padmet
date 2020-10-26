#!/usr/bin/env python3.7
# -*- coding: utf-8 -*-
"""
Description:
    Prot2Genome contains functions used for blast analysis and padmet enrichment
    
::

usage:
    padmet prot2genome --query_faa=FILE --query_ids=FILE/STR --subject_gbk=FILE --subject_fna=FILE --subject_faa=FILE --output_folder=FILE [--cpu=INT] [blastp] [tblastn] [debug]
    padmet prot2genome --query_faa=FILE --query_ids=FILE/STR --subject_gbk=FILE --subject_fna=FILE --subject_faa=FILE --output_folder=FILE --exonerate=PATH  [--cpu=INT] [blastp] [tblastn] [debug]
    padmet prot2genome --padmet=FOLDER --output=FOLDER
    padmet prot2genome --studied_organisms=FOLDER --output=FOLDER
    padmet prot2genome --run=FOLDER --padmetRef=FILE [--cpu=INT] [debug]

    From aucome run fromAucome():
        -1. Extract specifique reactions in spec_reactions folder with extractReactions()
        -2. Extract genes from spec_reactions files with extractGenes()
        -3. Run tblastn + exonerate with runAllAnalysis()
    
options:
    --query_faa=FILE #TODO.
    --query_ids=FILE/STR #TODO.
    --subject_gbk=FILE #TODO.
    --subject_fna=FILE #TODO.
    --subject_faa=FILE #TODO.
    --output_folder=FILE #TODO.
    --cpu=INT     Number of cpu to use for the multiprocessing (if none use 1 cpu). [default: 1]
    blastp #TODO.
    tblastn #TODO.
    debug #TODO.
"""

import csv
import docopt
import itertools
import os
import shutil
import subprocess
import sys

from Bio.Blast.Applications import NcbiblastpCommandline, NcbitblastnCommandline
from Bio.SeqFeature import FeatureLocation
from Bio import SearchIO
from Bio import SeqIO
from padmet.classes import PadmetSpec
from padmet.utils.management import manual_curation
from multiprocessing import Pool


def command_help():
    """
    Show help for analysis command.
    """
    print(docopt.docopt(__doc__))


def prot2genome_cli(command_args):
    """
    """
    global exonerate_path
    args = docopt.docopt(__doc__, argv=command_args)

    print("Test running exonerate...")
    subprocess.run(["exonerate"], shell = True)

    run_folder = args["--run"]
    cpu =  int(args["--cpu"])
    padmetRef = args["--padmetRef"]
    debug = args["debug"]

    if run_folder:
        fromAucome(run_folder, cpu, padmetRef, blastp=True, tblastn=True, exonerate=True, debug=debug)


#Utilise gbk to faa, ajouter une option pour faire un fichier fasta une sequence.
#pident a ajouter
#TODO:
#Seq to gbk, match to faa

def fromAucome(run_folder, cpu, padmetRef, blastp=True, tblastn=True, exonerate=True, keep_tmp=False, debug=False):
    """
    This function fit an AuCoMe run. Select a aucome run folder and then the function will:
    1./ For each couple of studied organisms, extract specific reactions
        ex: For org A and org B, extract reactions in org A but not in org B and vice versa
    2./ Then for each specific reactions, extract genes associated and run blastp, tblastn and exonerate
    3./ For each reaction, for all genes associated, if no blastp match but tblastn and exonerate hit select the reaction as a hit
    4./ Create a new padmet file with the new reactions to add within
    
    Parameters
    ----------
    run_folder: str
        path to aucome run folder
    cpu: int
        number of cpu to use for multiprocessing steps
    padmetRef: str
        path to padmetRef from where to extract and add the new reactions to create new padmet files
    blastp: bool
        If true run blastp during analysis
    tblastn: bool
        If true run tblastn during analysis
    exonerate: bool
        If true run exonerate during analysis, tblastn must also be True
    keep_tmp: bool
        If true keep temporary files of analysis (with predicted gene sequence)
    debug: bool
        if true, print all raw informations of analysis
    """
    prot2genome_folder = os.path.join(run_folder,"structural_check")
    padmet_folder = os.path.join(run_folder,"orthology_based","3_padmet_filtered")
    studied_organisms_folder = os.path.join(run_folder,"studied_organisms")
    spec_reactions_folder = os.path.join(prot2genome_folder, "0_specifics_reactions")
    blast_result_folder = os.path.join(prot2genome_folder, "1_blast_results")
    blast_analysis_folder = os.path.join(blast_result_folder, "analysis")
    tmp_folder = os.path.join(blast_result_folder, "tmp")
    reactions_to_add_folder = os.path.join(prot2genome_folder, "2_reactions_to_add")
    prot2genome_padmet_folder = os.path.join(prot2genome_folder, "3_PADMETs")

    if len(os.listdir(padmet_folder)) == 0:
        padmet_folder = os.path.join(run_folder,"orthology_based","2_padmet_orthology")
        if len(os.listdir(padmet_folder)) == 0:
            sys.exit('No input padmets inside either ' + os.path.join(run_folder,"orthology_based","3_padmet_filtered") + ' or ' + os.path.join(run_folder,"orthology_based","2_padmet_orthology"))

    pool = Pool(cpu)

    print("Extracting specific reactions...")
    mp_extractReactions(padmet_folder, spec_reactions_folder, pool)
    print("Running blast analysis...")
    mp_runAnalysis(spec_reactions_folder, studied_organisms_folder, blast_analysis_folder, tmp_folder, pool, blastp, tblastn, exonerate, keep_tmp, debug)
    print("Extracting reactions to add...")
    extractAnalysis(blast_analysis_folder, spec_reactions_folder, reactions_to_add_folder)
    print("Creating padmet files...")
    mp_createPadmet(reactions_to_add_folder, padmet_folder, prot2genome_padmet_folder, padmetRef, pool, verbose=True)

    pool.close()
    pool.join()

##### create padmets #####
def mp_createPadmet(reactions_to_add_folder, padmet_folder, output_folder, padmetRef, pool, verbose=False):
    """
    Update all padmet in padmet_folder with reactions to add from file in reactiosn_to_add_folder, the informations of the reactions are extracted from padmetRef as unique source
    ex: for padmet_folder/org_a.padmet, select reactions_to_add_folder/org_a.tsv, add each reactions listed in this file based on padmetRef to create  output_folder/org_a.padmet
    Create the padmet files in multiprocess, the more cpu the more new padmet files will be created faster

    Parameters
    ----------
    reactions_to_add_folder: str
        path folder with all files containing reactions to add for each studied organism
    padmet_folder: str
        path to folder with all padmet files of studied organism
    output_folder: str
        path to output folder where to create new padmet files
    padmetRef: str
        path to padmetRef from where to extract and add the new reactions to create new padmet files
    pool: Pool object
        pool object of multiprocessing
    verbose: bool
        verbose
    """
    all_padmets = next(os.walk(padmet_folder))[2]
    all_dict_args = []

    for padmet_file in all_padmets:
        if '.padmet' in padmet_file:
            dict_args = {}
            org_id = os.path.splitext(padmet_file)[0]
            padmet_to_update = os.path.join(padmet_folder, padmet_file)
            reactions_to_add_path = os.path.join(reactions_to_add_folder, org_id+".tsv")
            output = os.path.join(output_folder, padmet_file)
            dict_args = {"reactions_to_add_path": reactions_to_add_path, "padmet_to_update": padmet_to_update, "output":output, "padmetRef":padmetRef, "verbose": verbose}
            all_dict_args.append(dict_args)
        
    pool.map(createPadmet, all_dict_args)

    
    
def createPadmet(dict_args):
    """
    function used in mp_createPadmet by each worker the Pool
    padmet are updated using funciton add_delete_rxn from padmet.utils.connection.manual_curation
    """
    reactions_to_add_path = dict_args["reactions_to_add_path"]
    padmet_to_update = PadmetSpec(dict_args["padmet_to_update"])
    output = dict_args["output"]
    padmetRef = PadmetSpec(dict_args["padmetRef"])
    verbose = dict_args["verbose"]
    manual_curation.add_delete_rxn(reactions_to_add_path, padmet_to_update, output, padmetRef=padmetRef, category="MANUAL", verbose=verbose)
    



##### extract reactions #####
    
def mp_extractReactions(padmet_folder, output_folder, pool):
    """
    From a folder of padmet files, create all dual combination and extract specific reactions to create a file in output_folder
    ex: in padmet_folder: org_a.padmet, org_b.padmet, create: output_folder: org_a_vs_org_b.tsv and org_b_vs_org_a.tsv

    Parameters
    ----------
    padmet_folder: str
        path to folder with all padmet files of studied organism
    output_folder: str
        path to output folder where to extract specific reactions
    pool: Pool object
        pool object of multiprocessing
    """
    all_padmets = next(os.walk(padmet_folder))[2]
    all_padmets = [padmet for padmet in all_padmets if '.padmet' in padmet]
    all_combi = list(itertools.combinations(all_padmets, 2))
    #test_combi = [(org_a, org_b) for org_a,org_b in all_combi if org_a in ["Ectocarpus_siliculosus.padmet", "Ectocarpus_subulatus.padmet"] and org_b in ["Ectocarpus_siliculosus.padmet", "Ectocarpus_subulatus.padmet"]]
    all_dict_args = []
    for (org_a, org_b) in all_combi:
        path_a = os.path.join(padmet_folder,org_a)
        path_b = os.path.join(padmet_folder,org_b)
        output_a = os.path.join(output_folder, "%s_VS_%s.tsv"%(os.path.splitext(org_a)[0], os.path.splitext(org_b)[0]))
        output_b = os.path.join(output_folder, "%s_VS_%s.tsv"%(os.path.splitext(org_b)[0], os.path.splitext(org_a)[0]))
        if not os.path.isfile(output_a) and not os.path.isfile(output_b):
            dict_args = {"org_a": org_a, "path_a": path_a, "output_a": output_a, "org_b": org_b, "path_b": path_b, "output_b": output_b}
            all_dict_args.append(dict_args)

    pool.map(extractReactions, all_dict_args)


def extractReactions(dict_args):
    """
    function used in mp_cextractReactions by each worker the Pool
    for org_a.padmet and org_b.padmet:
        1./ extract reactions and specific reactiosn (not in a, not in b)
        2./ extract genes associated to specific reactions
        3./ Select only reactions if they are from annotation
        rxn-1 in org_a but not in org_b, if rxn-1 doesn't come from org_a annotation, skip the reaction
        4./ create output file: header = ["reaction_id", "genes_ids", "sources"]
    """
    org_a = dict_args["org_a"]
    path_a = dict_args["path_a"]
    output_a = dict_args["output_a"]

    org_b = dict_args["org_b"]
    path_b = dict_args["path_b"]
    output_b = dict_args["output_b"]

    padmet_a = PadmetSpec(path_a)
    padmet_b = PadmetSpec(path_b)

    # Get reactions from organism without spontaneous reactions.
    a_rxns = set([node.id for node in padmet_a.dicOfNode.values() if node.type == "reaction" and 'SPONTANEOUS' not in node.misc])
    b_rxns = set([node.id for node in padmet_b.dicOfNode.values() if node.type == "reaction" and 'SPONTANEOUS' not in node.misc])

    rxn_spec_to_a = a_rxns.difference(b_rxns)
    rxn_spec_to_b = b_rxns.difference(a_rxns)

    genes_assoc_spec_to_a = set()
    for rxn_id in rxn_spec_to_a:
        genes_assoc = set([rlt.id_out for rlt in padmet_a.dicOfRelationIn[rxn_id] if rlt.type == "is_linked_to"])
        genes_assoc_spec_to_a.update(genes_assoc)    

    genes_assoc_spec_to_b = set()
    for rxn_id in rxn_spec_to_a:
        genes_assoc = set([rlt.id_out for rlt in padmet_a.dicOfRelationIn[rxn_id] if rlt.type == "is_linked_to"])
        genes_assoc_spec_to_b.update(genes_assoc)    

    genes_ids_for_blast = set()
    nb_rxn_with_no_annot_source = 0
    with open(output_a,"w") as csvfile:
        header = ["reaction_id", "genes_ids", "sources"]
        dict_writer = csv.DictWriter(csvfile, fieldnames=header, delimiter="\t")
        dict_writer.writeheader()
        for rxn_id in rxn_spec_to_a:
            genes_assoc = set([rlt.id_out for rlt in padmet_a.dicOfRelationIn[rxn_id] if rlt.type == "is_linked_to"])
            recData_dict = {"ANNOTATION":None, "ORTHOLOGY": set()}
            for recData in [padmet_a.dicOfNode[rlt.id_out].misc for rlt in padmet_a.dicOfRelationIn[rxn_id] if rlt.type == "has_reconstructionData"]:
                category = recData["CATEGORY"][0]
                if category == "ANNOTATION":
                    recData_dict["ANNOTATION"] = recData["SOURCE"][0]
                elif category == "ORTHOLOGY":
                    org_info = recData["SOURCE"][0].replace("OUTPUT_ORTHOFINDER_FROM_","")
                    recData_dict["ORTHOLOGY"].add(org_info)
            if recData_dict["ANNOTATION"] is None:
                nb_rxn_with_no_annot_source += 1
                #print("reaction %s from %s has no annotation source" %(rxn_id, org_a))
            else:
                genes_ids_for_blast.update(genes_assoc)
                line = {"reaction_id": rxn_id, "genes_ids": ";".join(genes_assoc), "sources": "ANNOTATION+"+";".join(recData_dict["ORTHOLOGY"])}
                dict_writer.writerow(line)
    print("%s reaction in %s but not in %s have no annotation source" %(nb_rxn_with_no_annot_source, org_a, org_b))

    genes_ids_for_blast = set()
    nb_rxn_with_no_annot_source = 0
    with open(output_b,"w") as csvfile:
        header = ["reaction_id", "genes_ids", "sources"]
        dict_writer = csv.DictWriter(csvfile, fieldnames=header, delimiter="\t")
        dict_writer.writeheader()
        for rxn_id in rxn_spec_to_b:
            genes_assoc = set([rlt.id_out for rlt in padmet_b.dicOfRelationIn[rxn_id] if rlt.type == "is_linked_to"])
            recData_dict = {"ANNOTATION":None, "ORTHOLOGY": set()}
            for recData in [padmet_b.dicOfNode[rlt.id_out].misc for rlt in padmet_b.dicOfRelationIn[rxn_id] if rlt.type == "has_reconstructionData"]:
                category = recData["CATEGORY"][0]
                if category == "ANNOTATION":
                    recData_dict["ANNOTATION"] = recData["SOURCE"][0]
                elif category == "ORTHOLOGY":
                    org_info = recData["SOURCE"][0].replace("OUTPUT_ORTHOFINDER_FROM_","")
                    recData_dict["ORTHOLOGY"].add(org_info)
            if recData_dict["ANNOTATION"] is None:
                nb_rxn_with_no_annot_source += 1
                #print("reaction %s from %s has no annotation source" %(rxn_id, org_b))
            else:
                genes_ids_for_blast.update(genes_assoc)
                line = {"reaction_id": rxn_id, "genes_ids": ";".join(genes_assoc), "sources": "ANNOTATION+"+";".join(recData_dict["ORTHOLOGY"])}
                dict_writer.writerow(line)
    print("%s reaction in %s but not in %s have no annotation source" %(nb_rxn_with_no_annot_source, org_b, org_a))


##### run all analysis in multiprocess #####
def mp_runAnalysis(spec_reactions_folder, studied_organisms_folder, output_folder, tmp_folder, pool, blastp, tblastn, exonerate, keep_tmp, debug):
    """
    Run different blast analysis based on files representing specific reactions of 2 padmet files.
    For each specific reaction file in spec_reactions_folder (ex: org_a_vs_org_b.tsv):
        1./ search for:
            faa file of org_a (studied_organisms_folder/org_a/org_a.faa)
            gbk file of org_b (studied_organisms_folder/org_b/org_b.gbk)
            faa file of org_b (studied_organisms_folder/org_b/org_b.faa)
            fna file of org_b (studied_organisms_folder/org_b/org_b.fna)
                if fna doesn't exist create it
        2./ if output file (blast_analysis_folder/org_a_VS_org_b.tsv) doesn't already exist run analysis
        3./ extracts all genes ids from specific reaction file with fct extractGenes()
        4./ Run blastp, tblastn, exonerate on gene_id.faa vs target.faa / fna with runAllAnalysis()
        5./ Create analysis output
        The analysis create a lot of temp files, all are in tmp_folder wich is cleaned after all loop

    Parameters
    ----------
    spec_reactions_older: str
        path folder with all files containing specific reactions
    studied_organisms_folder: str
        path to folder with all data of studied organisms. Folder contains 1 folder by org with name as org name, in each: org.gbk,org.faa,org.fna
    output_folder: str
        path to output folder where to extract blast analysis
    tmp_folder: str
        path to tmp folder where to create faa of each gene to analyse
    pool: Pool object
        pool object of multiprocessing
    blastp: bool
        If true run blastp during analysis
    tblastn: bool
        If true run tblastn during analysis
    exonerate: bool
        If true run exonerate during analysis, tblastn must also be True
    keep_tmp: bool
        If true keep temporary files of analysis (with predicted gene sequence)
    debug: bool
        if true, print all raw informations of analysis
    """
    for rxn_file in [os.path.join(spec_reactions_folder, i) for i in next(os.walk(spec_reactions_folder))[2]]:
        org_a, org_b = os.path.basename(rxn_file).replace(".tsv","").split("_VS_")
        analysis_output = os.path.join(output_folder, "%s_VS_%s.tsv"%(org_a, org_b))
        query_faa = os.path.join(studied_organisms_folder, org_a, "%s.faa"%org_a)
        subject_gbk = os.path.join(studied_organisms_folder, org_b, "%s.gbk"%org_b)
        subject_faa = os.path.join(studied_organisms_folder, org_b, "%s.faa"%org_b)
        subject_fna = os.path.join(studied_organisms_folder, org_b, "%s.fna"%org_b)
        #If GBK given and no fna, create fna from genbank
        if not os.path.exists(subject_fna):
            SeqIO.convert(subject_gbk, "genbank", subject_fna, "fasta")

        if not os.path.isfile(analysis_output):
            all_query_seq_ids = extractGenes(rxn_file)
            #Create a list of dict, each dict is the arg given to runAllAnalysis fct
            print("Extracting all data for %s vs %s" %(org_a, org_b))
            print("%s query ids to search" %len(all_query_seq_ids))
    
            all_dict_args = []
            for query_seq_id in all_query_seq_ids:
                dict_args = {"query_seq_id": query_seq_id, "query_faa": query_faa, "subject_faa": subject_faa, "subject_fna": subject_fna, "output_folder": tmp_folder, "blastp": blastp, "tblastn": tblastn, "exonerate": exonerate, "debug": debug}
                all_dict_args.append(dict_args)
            #list of dict to give to dictWritter
            all_analysis_result = []
            #Run runAllAnalysis in multiproccess.
            mp_results = pool.map(runAllAnalysis, all_dict_args)
            for _list in mp_results:
                all_analysis_result += _list
            print("Creating output analysis: %s" %analysis_output)
            #Create output file
            analysisOutput(all_analysis_result, analysis_output)
            if not keep_tmp:
                cleanTmp(tmp_folder)


def cleanTmp(tmp_folder):
    """
    Remove all files from tmp folder

    Parameters
    ----------
    tmp_folder: str
        path to tmp folder where to create faa of each gene to analyse
    """
    for tmp_folder_file in os.listdir(tmp_folder):
        folder_file_path = os.path.join(tmp_folder, tmp_folder_file)
        if os.path.isfile(folder_file_path):
            os.unlink(folder_file_path)
        if os.path.isdir(folder_file_path):
            shutil.rmtree(folder_file_path)


def extractGenes(reactions_file):
    """
    Extract genes ids and return a list from reactions_file obtained with extractReactions()

    Parameters
    ----------
    reactions_file: str
        path to reaction file
    """
    all_query_seq_ids = set()
    with open(reactions_file, 'r') as csvfile:
        reader = csv.DictReader(csvfile, delimiter="\t")
        [all_query_seq_ids.update(line["genes_ids"].split(";")) for line in reader]
    return all_query_seq_ids


def runAllAnalysis(dict_args):
    """
    For a given gene query id:
        1/ extract from query_faa the sequence and create a faa file output_folder/query_id.faa
            If isoforms found, also search for each specific isoform
        2/ if blastp, run blastp; if tblastn, run tblastn; if exonerate and tblastn has hit, run exonerate
        Run all of them and extract output as dict of data

    Returns
    -------
    list
        list of dict with all analysis output
    """
    query_seq_id = dict_args["query_seq_id"]
    query_faa = dict_args["query_faa"]
    subject_faa = dict_args["subject_faa"]
    subject_fna = dict_args["subject_fna"]
    output_folder = dict_args["output_folder"]
    blastp = dict_args["blastp"]
    tblastn = dict_args["tblastn"]
    exonerate = dict_args["exonerate"]
    debug = dict_args["debug"]

    with open(query_faa, "r") as faa:
        query_seqs = [seq_record for seq_record in SeqIO.parse(faa, "fasta") if seq_record.id.startswith(query_seq_id+"_isoform") or seq_record.id == query_seq_id]
    if len(query_seqs) > 1:
        print("/!\ Isoforms found for %s: %s"%(query_seq_id, [i.name for i in query_seqs]))
    analysis_result = list()
    for query_seq in query_seqs:
        query_seq_id = query_seq.id
        query_seq_faa = os.path.join(output_folder,query_seq_id+".faa")
        if not os.path.exists(query_seq_faa):
            SeqIO.write(query_seq, query_seq_faa, "fasta")
        current_result = {"query_seq_id": query_seq_id,}
    
        # Run BLASTP and parse the output
        if blastp:
            blastp_result = runBlastp(query_seq_faa, subject_faa, debug=debug)
            if blastp_result:
                current_result.update(blastp_result)
        # Run TBLASTN and parse the output
        if tblastn:
            tblastn_result = runTblastn(query_seq_faa, subject_fna, debug=debug)
            if tblastn_result:
                current_result.update(tblastn_result)
        
        # Run exonerate, parse output
        if exonerate:
            if not tblastn:
                IOError("must run tblastn to be able to run exonerate")
            _tblastn_hit = True
            try:
                exonerate_target_id = current_result["tblastn_sseqid"]
                start_match = int(tblastn_result['tblastn_sstart'])
                end_match = int(tblastn_result['tblastn_send'])
            except KeyError:
                _tblastn_hit = False
                if debug:
                    print("No hit from tblastn, can't run exonerate")
            if _tblastn_hit:
                output_folder = output_folder + '/' + query_seq_id + '_' + exonerate_target_id
                if not os.path.exists(output_folder):
                    os.mkdir(output_folder)
                sseq_seq_faa = os.path.join(output_folder,exonerate_target_id+".fna")
                createSeqFromTblastn(subject_fna, sseq_seq_faa, exonerate_target_id, start_match, end_match)
                exonerate_output = os.path.join(output_folder, "exonerate_output_%s_vs_%s.txt"%(query_seq_id, exonerate_target_id))
                exonerate_result = runExonerate(query_seq_faa, sseq_seq_faa, exonerate_output, debug=debug)
                exonerate_sequence = os.path.join(output_folder, "%s_vs_%s.fasta"%(query_seq_id, exonerate_target_id))
                extract_sequence(exonerate_output, exonerate_sequence)
                current_result.update(exonerate_result)
                if 'exonerate_hit_range' in current_result:
                    current_result['exonerate_hit_range'] = '-'.join([str(hit_range) for hit_range in current_result['exonerate_hit_range']])
        analysis_result.append(current_result)

    return analysis_result


def extract_sequence(exonerate_output, exonerate_sequence):
    """
    Extract protein sequence from exonerate ouput.
    """
    with open(exonerate_sequence, 'w') as output_file:
        with open(exonerate_output, 'r') as input_file:
            for chunk in input_file.read().split('Fasta_seq'):
                if '-- completed exonerate analysis' in chunk:
                    fasta_seq = chunk.split('-- completed exonerate analysis')[0]
                    output_file.write(fasta_seq)


def runSearchOnProteome(proteome_orgA, genome_orgB, output_folder, proteome_orgB=None):
    """
    From a proteome of OrgA search for missing structural annotation in genome of OrgB.
    First launch Blastp between proteome of OrgA and proteome of OrgB.
    Then launch tBlastn between proteome of OrgA and genome of OrgB to find matches.
    Use the best match to extract a region from the genome of OrgB.
    Then launch Exonerate on this region using the sequence of OrgA.

    Parameters
    ----------
    proteome_orgA: str
        path to fasta file of proteome of OrgA
    genome_orgB: str
        path to fasta file of genome of OrgB
    output_folder: str
        path to output folder
    proteome_orgB: str
        path to fasta file of proteome of OrgB
    """
    if not os.path.exists(output_folder):
        os.mkdir(output_folder)

    analysis_result = []

    if proteome_orgB:
        current_result = {"query_seq_id": proteome_orgB,}
        blastp_result = runBlastp(proteome_orgA, proteome_orgB, header=["sseqid", "evalue", "bitscore"], debug=False)
        if blastp_result:
            current_result.update(blastp_result)
    else:
        current_result = {}

    tblast_hit = None
    tblastn_result = runTblastn(proteome_orgA, genome_orgB)
    if tblastn_result:
        current_result.update(tblastn_result)
    try:
        exonerate_target_id = tblastn_result['tblastn_sseqid']
        start_match = int(tblastn_result['tblastn_sstart'])
        end_match = int(tblastn_result['tblastn_send'])
        tblast_hit = True
    except KeyError:
        print("No hit from tblastn, can't run exonerate")

    if tblast_hit:
        sseq_seq_faa = output_folder + "/" + exonerate_target_id + ".fna"
        createSeqFromTblastn(genome_orgB, sseq_seq_faa, exonerate_target_id, start_match, end_match)

        exonerate_output = output_folder + "/exonerate_output.txt"
        exonerate_result = runExonerate(proteome_orgA, sseq_seq_faa, exonerate_output, debug=False)
        current_result.update(exonerate_result)
        analysis_result.append(current_result)

        prot2genomes_output = output_folder + "/prot2genomes_result.tsv"
        analysisOutput(analysis_result, prot2genomes_output)


def runBlastp(query_seq_faa, subject_faa, header=["sseqid", "evalue", "bitscore"], debug=False):
    """
    Run blastp on querry_seq vs subectj faa and return output based on header
    Use NcbiblastpCommandline fct and extract output
    Extract 1st best hit based on bitscore

    Parameters
    ----------
    query_seq_faa: str
        path to query fasta sequence
    subject_faa: str
        path to subject fasta sequence
    header: list
        output format of blastp
    debug: bool
        if true print all raw blastp output

    Returns
    -------
    dict
        dict of the best blastp hit, add 'blastp_' tag, or empty dict if no hit
    """
    if debug:
        print("\tRunning Blastp %s vs %s" %(os.path.basename(query_seq_faa), os.path.basename(subject_faa)))
    outfmt_arg = '"%s %s"'%(6, " ".join(header))
    output = NcbiblastpCommandline(query=query_seq_faa, subject=subject_faa, evalue=1e-20, outfmt=outfmt_arg)()[0]
    output = [line.split("\t") for line in output.splitlines()]
    blastp_result = {}
    count = 0
    for line in output:
        count += 1
        blastp_result[count] = {}
        for (index, h) in enumerate(header):
            blastp_result[count]["blastp_"+h] = line[index]
    try:
        max_bitscore = max([float(hsp["blastp_bitscore"]) for hsp in blastp_result.values()])
    except ValueError:
        max_bitscore = "N.A"
    if debug:
        print("\t\tAll hsp from blastp:")
        if blastp_result.values():
            for hsp in blastp_result.values():
                print(hsp)
            print("\t\tMax biscore: %s"%max_bitscore)
        else:
            print("\t\tNo HPS")
    try:
        result = [hsp for hsp in blastp_result.values() if float(hsp["blastp_bitscore"]) == max_bitscore][0]
    except IndexError:
        result = {}
    return result


def runTblastn(query_seq_faa, subject_fna, header=["sseqid", "evalue", "bitscore", "sstart", "send"], debug=False):
    """
    Run tblastn on querry_seq vs subectj fna and return output based on header
    Use NcbitblastnCommandline fct and extract output
    Extract 1st best hit based on bitscore

    Parameters
    ----------
    query_seq_faa: str
        path to query fasta sequence
    subject_fna: str
        path to subject fna sequence
    header: list
        output format of tblastn
    debug: bool
        if true print all raw tblastn output

    Returns
    -------
    dict
        dict of the best tblastn hit, add 'tblastn_' tag, or empty dict if no hit
    """    
    print("\tRunning tBlastn %s vs %s" %(os.path.basename(query_seq_faa), os.path.basename(subject_fna)))
    outfmt_arg = '"%s %s"'%(6, " ".join(header))
    output = NcbitblastnCommandline(query=query_seq_faa, subject=subject_fna, evalue=1e-20, outfmt=outfmt_arg)()[0]
    output = [line.split("\t") for line in output.splitlines()]
    tblastn_result = {}
    count = 0
    for line in output:
        count += 1
        tblastn_result[count] = {}
        for (index, h) in enumerate(header):
            tblastn_result[count]["tblastn_"+h] = line[index]
    try:
        max_bitscore = max([float(hsp["tblastn_bitscore"]) for hsp in tblastn_result.values()])
    except ValueError:
        max_bitscore = "N.A"
    if debug:
        print("\t\tAll hsp from tblastn:")
        for hsp in tblastn_result.values():
            print(hsp)
        print("\t\tMax biscore: %s"%max_bitscore)
    try:           
        result = [hsp for hsp in tblastn_result.values() if float(hsp["tblastn_bitscore"]) == max_bitscore][0]
    except IndexError:
        result = {}
    return result


def createSeqFromTblastn(subject_fna, sseq_seq_faa, exonerate_target_id, start_match, end_match):
    """
    Use the result from the tBlastn to extract a region from the subject genome.
    The region extracted corresponds to the match region and 10kb before and 10kb after.

    Parameters
    ----------
    subject_fna: str
        path to subject fasta sequence (genome)
    sseq_seq_faa: str
        path to output fasta sequence
    exonerate_target_id: str
        ID of the contig/scaffold/chromosome where a match has been found
    start_match: int
        start of the match
    end_match: int
        end of the match
    """
    if not os.path.exists(sseq_seq_faa):
        with open(subject_fna, "r") as fna:
            sseq_seq = [seq_record for seq_record in SeqIO.parse(fna, "fasta") if seq_record.id == exonerate_target_id][0]
            sseq_seq.description = "tblastn identified sequence"
            sseq_seq.id = exonerate_target_id + "_" + str(start_match) + "_" + str(end_match)
            if start_match > end_match:
                start_match = start_match + 10000
                end_match = end_match - 10000
                if start_match > len(sseq_seq.seq):
                    start_match = len(sseq_seq.seq)
                if end_match < 0:
                    end_match = 0
                seq_location = FeatureLocation(end_match, start_match)
                sseq_seq.seq = seq_location.extract(sseq_seq.seq)
            elif start_match < end_match:
                start_match = start_match - 10000
                end_match = end_match + 10000
                if start_match < 0:
                    start_match = 0
                if end_match > len(sseq_seq.seq):
                    end_match = len(sseq_seq.seq)
                seq_location = FeatureLocation(start_match, end_match)
                sseq_seq.seq = seq_location.extract(sseq_seq.seq)
            SeqIO.write(sseq_seq, sseq_seq_faa, "fasta")


def runExonerate(query_seq_faa, sseq_seq_faa, output, debug=False):
    """
    Run exonerate on querry_seq vs subject faa
    Exonerate must be installed, and the global var PATH must be update with the exonerate/bin/
    command 'exonerate' should work from shell
    sseq_seq_faa is obtained after tblastn run based on tblastn_sseqid value
    
    Parameters
    ----------
    query_seq_faa: str
        path to query fasta sequence
    sseq_seq_faa: str
        path to subject faa sequence
    output: str
        path to exonerate output
    debug: bool
        if true print all raw exonerate output

    Returns
    -------
    dict
        dict of the best exonerate hit, add 'exonerate_' tag, or empty dict if no hit
    """
    exonerate_path = "exonerate"
    print("\tRunning Exonerate %s vs %s" %(os.path.basename(query_seq_faa), os.path.basename(sseq_seq_faa)))
    exonerate_result = {}
    cmd_args = '{0} --model protein2genome {1} {2} --score 500 --bestn 1 --percent 50 --showquerygff True  \
    --ryo ">%qi length=%ql alnlen=%qal\n>%ti length=%tl alnlen=%tal\nFasta_seq\n>%ti\n%qs" >> {3}'.format(exonerate_path, query_seq_faa, sseq_seq_faa, output)
    subprocess.run([cmd_args], shell = True)
    try:
        exonerate_raw_output = list(SearchIO.parse(output, 'exonerate-text'))
       
        best_hsp = exonerate_raw_output[0][0][0]
        for qresult in exonerate_raw_output:
            for hit in qresult:
                for hsp in hit:
                    if hsp.score > best_hsp.score:
                        best_hsp = hsp
        if debug:
            print("\t\tAll hsp from exonerate:")
            for qresult in exonerate_raw_output:
                for hit in qresult:
                    for hsp in hit:
                        print(hsp)
            print("\t\tMax score: %s"%best_hsp.score)
        exonerate_result = {"exonerate_score": best_hsp.score, "exonerate_hit_range": best_hsp.hit_range}
    except (RuntimeError, IndexError) as e:
        if debug:
            print("\t\tNo HSP")
    return exonerate_result


##### extract analysis results #####
def extractAnalysis(blast_analysis_folder, spec_reactions_folder, output_folder):
    """
    For each analysis output in blast analysis folder, obtained with runAllAnalysis()
        1./ Extract orthologues hit
        2./ For each specific reactions from spec_reactions_folder, if all genes of a reactions got ortho hit
            add reaction to reactions_to_add

    Parameters
    ----------
    blast_analysis_folder: str
        path folder with all blast analysis output files
    spec_reactions_folder: str
        path folder with all files containing specific reactions
    output_folder: str
        path folder where to extract all reactions to add

    """
    # {org_id: {org_id:set(genes_ids,...), }, }
    orthologue_dict = {}
    for analysis_file in [os.path.join(blast_analysis_folder, i) for i in next(os.walk(blast_analysis_folder))[2]]:
        org_b, org_a = os.path.basename(analysis_file).replace(".tsv","").split("_VS_")
        try:
            orthologue_dict[org_a][org_b] = set()
        except KeyError:
            orthologue_dict[org_a] = {org_b: set()}

        with open(analysis_file, 'r') as csvfile:
            reader = csv.DictReader(csvfile, delimiter="\t")
            for line in reader:
                if not line['blastp_bitscore'] and (line['exonerate_score'] and line['tblastn_bitscore']):
                    orthologue_dict[org_a][org_b].add(line['query_seq_id'])

    # {org_id: {reaction_id:{org_id:{total_genes: set(genes_ids), orthologues: set(genes_ids)}}}
    reactions_dict = {}
    for rxn_file in [os.path.join(spec_reactions_folder, i) for i in next(os.walk(spec_reactions_folder))[2]]:
        org_b, org_a = os.path.basename(rxn_file).replace(".tsv","").split("_VS_")
        if org_a not in reactions_dict.keys():
            reactions_dict[org_a] = dict()
        all_orthologues = orthologue_dict[org_a][org_b]
        with open(rxn_file, 'r') as csvfile:
            reader = csv.DictReader(csvfile, delimiter="\t")
            for line in reader:
                reaction_id = line['reaction_id']
                genes_ids = set(line['genes_ids'].split(";"))
                orthologues = genes_ids.intersection(all_orthologues)
                try:
                    reactions_dict[org_a][reaction_id][org_b] = {'total_genes': genes_ids, 'orthologues': orthologues}
                except KeyError:
                    reactions_dict[org_a][reaction_id] = {org_b: {'total_genes': genes_ids, 'orthologues': orthologues}}
                    
    for org_a, org_dict in reactions_dict.items():
        output_file = os.path.join(output_folder, "%s.tsv"%org_a)
        with open(output_file, 'w') as csvfile:
            fieldnames = ['idRef', 'Comment', 'Genes', 'Action']
            writer = csv.DictWriter(csvfile,fieldnames, delimiter="\t")
            writer.writeheader()
            for reaction_id, reaction_dict in org_dict.items():
                total_org = set(reaction_dict.keys())
                org_with_orthologues = set([org_b for org_b, org_b_dict in reaction_dict.items() if org_b_dict['total_genes'] == org_b_dict['orthologues']])
                rate = round(float(len(org_with_orthologues)/len(total_org)),2)
                if rate > 0.0:
                    line = {'idRef': reaction_id, 'Comment': 'Has Orthologues with %s'%";".join(org_with_orthologues),'Action':'add'}
                    writer.writerow(line)


def analysisOutput(analysis_result, analysis_output):
    """
    """
    analysis_header = ["query_seq_id", "blastp_sseqid", "blastp_evalue", "blastp_bitscore",  "tblastn_sseqid", "tblastn_evalue", "tblastn_bitscore", "tblastn_sstart", "tblastn_send", "exonerate_score", "exonerate_hit_range"]
    with open(analysis_output,"w") as csvfile:
        dict_writer = csv.DictWriter(csvfile, fieldnames=analysis_header, delimiter="\t")
        dict_writer.writeheader()            
        dict_writer.writerows(analysis_result)    
