# -*- coding: utf-8 -*-

from __future__ import absolute_import

import pandas as pd
import numpy as np
import cobra

def met_info_from_model(model, verbose=True):
    '''
    This function looks for all the reactions in the model and returns their respective information.

    Parameters
    ----------
    model : cobra.core.Model.Model
        A cobra model.

    verbose : boolean, True by default.
        A variable to enable or disable the printings of this function.

    Returns
    -------
    met_rxn_gene_association : pandas.DataFrame
        A pandas dataframe containing the information retrieved. The columns are :
        'MetID', GeneID', 'RxnID', 'RxnName', 'SubSystem', 'RxnFormula'
    '''
    if verbose:
        print('Getting information for all metabolites in the model.')

    met_association = []
    for met in model.metabolites:
        for rxn in met.reactions:
            if len(rxn.genes) != 0:
                for gene in rxn.genes:
                    met_association.append((str(met.id), str(gene.id), rxn.id, rxn.name, rxn.subsystem, rxn.reaction))
            else:
                met_association.append((str(met.id), '', rxn.id, rxn.name, rxn.subsystem, rxn.reaction))
    labels = ['MetID', 'GeneID', 'RxnID', 'RxnName', 'SubSystem', 'RxnFormula']
    met_association = pd.DataFrame.from_records(met_association, columns=labels)
    if verbose:
        print('Information correctly obtained.')
    return met_association

def get_metabolites_from_type(metabolite_list,type):
    '''
    This function returns a list of metabolites of a certain type.

    Parameters
    ----------
    metabolite_list : list of cobra.core.metabolite.Metabolite
        A list of cobra metabolites.

    type : string
        A string specifying the type of metabolites to be retrieved. 
        Possible types are: 'aminoacids','lipids','carbohydrates','DNA','RNA'

    Returns
    -------
    filtered_metabolite_list : list
        A list containing the metabolites of a certain type contained in the model.
    '''
    import re
    filtered_metabolite_list = []

    if type == 'aminoacids':
        aminoacid_identifier = '__L_'
        exclude = ['hom','sbt','srb','lac']
        for met in metabolite_list:
            if aminoacid_identifier in met.id \
                and len(met.id.split(aminoacid_identifier)[0]) == 3 \
                and met.id.split(aminoacid_identifier)[0] not in exclude \
                or met.id == 'gly_c':
                filtered_metabolite_list.append(met)

    # if type == 'lipids':
    #     lipid_identifiers = ['triglyc','mgdg','dgdg','tag','dag','mag','ino',\
    #                             'zymst','ergst']
    #     for met in metabolite_list:
    #             for lipid_id in lipid_identifiers:
    #                 if lipid_id in met.id:
    #                     filtered_metabolite_list.append(met)
    #                     break

    if type == 'DNA':
        DNA_identifiers = ['da.p','dg.p','dt.p','dc.p']
        for met in metabolite_list:
                for DNA_id in DNA_identifiers:
                    DNA_regex = re.compile(DNA_id)
                    if re.match(DNA_regex,met.id):
                        filtered_metabolite_list.append(met)
                        break

    if type == 'RNA':
        RNA_identifiers = ['a.p','g.p','t.p','c.p','u.p']
        for met in metabolite_list:
                for RNA_id in RNA_identifiers:
                    RNA_regex = re.compile(RNA_id)
                    if re.match(RNA_regex,met.id) \
                        and len(met.id) == 5:

                        filtered_metabolite_list.append(met)
                        break

    if type == 'storage_carbohydrates':
        carb_identifiers = ['glycogen','.+ose']
        for met in metabolite_list:
                for carb_id in carb_identifiers:
                    carb_regex = re.compile(carb_id)
                    if re.match(carb_regex,met.id) or re.match(carb_regex,met.name):
                        filtered_metabolite_list.append(met)
                        break
                        
    if type == 'cell_wall':
        CW_identifiers = ['mannan','.+glcn','pa_','pc_','pe_','ps_']
        for met in metabolite_list:
            for CW_id in CW_identifiers:
                CW_regex = re.compile(CW_id)
                if re.match(CW_regex,met.id) or re.match(CW_regex,met.name):
                    filtered_metabolite_list.append(met)
                    break

    return filtered_metabolite_list

def classify_metabolites_by_type(metabolite_list):
    '''
    This function returns a dictionary with the classification of metabolites
    by type.

    Parameters
    ----------
    metabolite_list : list of cobra.core.metabolite.Metabolite
        A list of cobra metabolites.

    Returns
    -------
    classification_dict : dict
        A dictionary containing the metabolites classified by type.
    '''

    types = ['aminoacids','storage_carbohydrates','DNA','RNA','cell_wall']

    classification_dict = dict()
    all_mets =[]
    [all_mets.append(met.id) for met in metabolite_list]
    included_mets = []
    for type in types:
        classification_dict[type] =[]
        met_list = get_metabolites_from_type(metabolite_list,type=type)
        for met in met_list:
            included_mets.append(met.id)
            classification_dict[type].append(met.id)

    classification_dict['other'] = list(np.setdiff1d(all_mets,included_mets))

    return classification_dict

def get_molecular_weight(model,met):
    met = model.metabolites.get_by_id(met.id)
    formula = cobra.core.formula.Formula(met.formula)
    met_weight = formula.weight

    return met_weight