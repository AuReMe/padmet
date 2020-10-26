#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pkg_resources

from padmet.classes import PadmetRef, PadmetSpec
from datetime import datetime


def instantiate_padmet(padmet_type, padmetRef_file=None, padmet_id=None, db='NA', version='NA', verbose=None):
    now = datetime.now()
    today_date = now.strftime("%Y-%m-%d")
    padmet_version = pkg_resources.require("padmet")[0].version

    POLICY_IN_ARRAY = [['class','is_a_class','class'], ['class','has_name','name'], ['class','has_xref','xref'], ['class','has_suppData','suppData'],
                    ['compound','is_a_class','class'], ['compound','has_name','name'], ['compound','has_xref','xref'], ['compound','has_suppData','suppData'],
                    ['gene','is_a_class','class'], ['gene','has_name','name'], ['gene','has_xref','xref'], ['gene','has_suppData','suppData'], ['gene','codes_for','protein'],
                    ['pathway','is_a_class','class'], ['pathway','has_name','name'], ['pathway','has_xref','xref'], ['pathway','is_in_pathway','pathway'],
                    ['protein','is_a_class','class'], ['protein','has_name','name'], ['protein','has_xref','xref'], ['protein','has_suppData','suppData'], ['protein','catalyses','reaction'],
                    ['protein','is_in_species','class'],
                    ['reaction','is_a_class','class'], ['reaction','has_name','name'], ['reaction','has_xref','xref'], ['reaction','has_suppData','suppData'], ['reaction','has_reconstructionData','reconstructionData'], ['reaction','is_in_pathway','pathway'],
                    ['reaction','consumes','class','STOICHIOMETRY','X','COMPARTMENT','Y'], ['reaction','produces','class','STOICHIOMETRY','X','COMPARTMENT','Y'],
                    ['reaction','consumes','compound','STOICHIOMETRY','X','COMPARTMENT','Y'], ['reaction','produces','compound','STOICHIOMETRY','X','COMPARTMENT','Y'],
                    ['reaction','consumes','protein','STOICHIOMETRY','X','COMPARTMENT','Y'], ['reaction','produces','protein','STOICHIOMETRY','X','COMPARTMENT','Y'],
                    ['reaction','is_linked_to','gene','SOURCE:ASSIGNMENT','X:Y']]
    dbNotes = {"PADMET":{"creation":today_date,"version":padmet_version},"DB_info":{"DB":db,"version":version}}

    if padmet_type == "PadmetRef":
        padmet = PadmetRef()
    elif padmet_type == "PadmetSpec":
        padmet = PadmetSpec()
    else:
        print("Wrong padmet_type given, it must be a string with either: PadmetRef or PadmetSpec.")

    # If there is a padmetref, use it as a template.
    if padmetRef_file:
        padmetRef = PadmetRef(padmetRef_file)
        version = padmetRef.info["DB_info"]["version"]
        db = padmetRef.info["DB_info"]["DB"]
        dbNotes = {"PADMET":{"creation":today_date,"version":padmet_version},"DB_info":{"DB":db,"version":version}}
        padmet.setInfo(dbNotes)

        padmet.setPolicy(padmetRef)

    # Otherwise use the information given by the user.
    else:
        if verbose: print("setting policy")
        padmet.setPolicy(POLICY_IN_ARRAY)
        if verbose: print("setting dbInfo")
        padmet.setInfo(dbNotes)

    if padmet_id:
        padmet.info["DB_info"]["ID"] = padmet_id

    return padmet