# -*- coding: utf-8 -*-

from __future__ import absolute_import

import pandas as pd
import cobra
from cobra_utils.query.met_info import classify_metabolites_by_type


def rxn_info_from_metabolites(model, metabolites, verbose=True):
    '''
    This function looks for all the reactions where the metabolites in the list participate. Also, it retrieves the genes
    associated to those reactions.

    Parameters
    ----------
    model : cobra.core.Model.Model
        A cobra model.

    metabolites : array-like
        An iterable object containing a list of metabolite ids present in the model.

    verbose : boolean, True by default.
        A variable to enable or disable the printings of this function.

    Returns
    -------
    rxn_gene_association : pandas.DataFrame
        A pandas dataframe containing the information retrieved. The columns are :
        'MetName', 'MetID', 'RxnID', 'RxnName', 'GeneID', 'Subsystem', 'RxnFormula'
    '''
    if verbose:
        print('Using list of metabolites to get reactions where they participate. Also, getting genes of those reactions.')

    rxn_gene_association = []
    for metabolite in metabolites:
        met = model.metabolites.get_by_id(metabolite)
        for rxn in met.reactions:
            if len(rxn.genes) != 0:
                for gene in rxn.genes:
                    rxn_gene_association.append(
                        (rxn.id, rxn.name, str(gene.id), rxn.subsystem, rxn.reaction, met.id, met.name))
            else:
                rxn_gene_association.append(
                    (rxn.id, rxn.name, '', rxn.subsystem, rxn.reaction, met.id, met.name))

    labels = ['RxnID', 'RxnName', 'GeneID', 'Subsystem', 'RxnFormula', 'MetID', 'MetName']
    rxn_gene_association = pd.DataFrame.from_records(rxn_gene_association, columns=labels)
    if verbose:
        print('Information correctly obtained.')
    return rxn_gene_association


def rxn_info_from_reactions(model, reactions, verbose=True):
    '''
    This function looks for all the reactions and genes that are associated from a list of reactions ids.

    Parameters
    ----------
    model : cobra.core.Model.Model
        A cobra model.

    reactions : array-like
        An iterable object containing a list of reaction ids present in the model.

    verbose : boolean, True by default.
        A variable to enable or disable the printings of this function.

    Returns
    -------
    rxn_gene_association : pandas.DataFrame
        A pandas dataframe containing the information retrieved. The columns are :
        'GeneID', 'RxnID', 'RxnName', 'SubSystem', 'RxnFormula'
    '''
    if verbose:
        print('Using list of reactions to get their information and genes associated.')

    rxn_gene_association = []
    for reaction in reactions:
        rxn = model.reactions.get_by_id(reaction)
        if len(rxn.genes) != 0:
            for gene in rxn.genes:
                rxn_gene_association.append((str(gene.id), rxn.id, rxn.name, rxn.subsystem, rxn.reaction))
        else:
            rxn_gene_association.append(('', rxn.id, rxn.name, rxn.subsystem, rxn.reaction))
    labels = ['GeneID', 'RxnID', 'RxnName', 'SubSystem', 'RxnFormula']
    rxn_gene_association = pd.DataFrame.from_records(rxn_gene_association, columns=labels)
    if verbose:
        print('Information correctly obtained.')
    return rxn_gene_association


def rxn_info_from_genes(model, genes, verbose=True):
    '''
    This function looks for all the reactions and genes that are associated from a list of gene ids.

    Parameters
    ----------
    model : cobra.core.Model.Model
        A cobra model.

    genes : array-like
        An iterable object containing a list of gene ids present in the model.

    verbose : boolean, True by default.
        A variable to enable or disable the printings of this function.

    Returns
    -------
    rxn_gene_association : pandas.DataFrame
        A pandas dataframe containing the information retrieved. The columns are :
        'GeneID', 'RxnID', 'RxnName', 'SubSystem', 'RxnFormula'
    '''
    if verbose:
        print('Using list of genes to get the reactions associated and their information.')

    rxn_gene_association = []
    for gene in genes:
        g = model.genes.get_by_id(gene)
        for rxn in g.reactions:
                rxn_gene_association.append((str(g.id), rxn.id, rxn.name, rxn.subsystem, rxn.reaction))
    labels = ['GeneID', 'RxnID', 'RxnName', 'SubSystem', 'RxnFormula']
    rxn_gene_association = pd.DataFrame.from_records(rxn_gene_association, columns=labels)
    if verbose:
        print('Information correctly obtained.')
    return rxn_gene_association


def rxn_info_from_model(model, verbose=True):
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
    rxn_gene_association : pandas.DataFrame
        A pandas dataframe containing the information retrieved. The columns are :
        'GeneID', 'RxnID', 'RxnName', 'SubSystem', 'RxnFormula'
    '''
    if verbose:
        print('Getting information for all reactions in the model.')

    rxn_gene_association = []
    for rxn in model.reactions:
        if len(rxn.genes) != 0:
            for gene in rxn.genes:
                rxn_gene_association.append((str(gene.id), rxn.id, rxn.name, rxn.subsystem, rxn.reaction))
        else:
            rxn_gene_association.append(('', rxn.id, rxn.name, rxn.subsystem, rxn.reaction))
    labels = ['GeneID', 'RxnID', 'RxnName', 'SubSystem', 'RxnFormula']
    rxn_gene_association = pd.DataFrame.from_records(rxn_gene_association, columns=labels)
    if verbose:
        print('Information correctly obtained.')
    return rxn_gene_association


def get_objective_function(model):
    '''
    This function returns the reaction set as the objective function for FBA.

    Parameters
    ----------
    model : cobra.core.Model.Model
        A cobra model.

    Returns
    -------
    reaction : cobra.core.reaction.Reaction
        A cobra reaction.
    '''
    
    for rxn in model.reactions:
        if rxn.objective_coefficient:
            obj_rxn = rxn
            print(rxn.id)
    return obj_rxn


def get_reaction_stoichiometry(reaction):
    '''
    This function returns the stoichiometry of a reaction

    Parameters
    ----------
    reaction : cobra.core.reaction.Reaction
        A cobra reaction.

    Returns
    -------
    reaction_stoichiometry : dict
        A dictionary containing the metabolite stoichiometry of the reaction.
    '''

    reaction_stoichiometry = dict()
    for met in reaction.metabolites:
        stoich = reaction.get_coefficient(met.id)
        reaction_stoichiometry[met.id] = stoich
    return reaction_stoichiometry

def biomass_breakdown(model,input_info, input_mode = 'reaction'):

    if input_mode == 'reaction':
        biomass_rxn = input_info
        stoich = get_reaction_stoichiometry(biomass_rxn)
        class_dict = classify_metabolites_by_type(biomass_rxn.metabolites)
    elif input_mode == 'dict':
        stoich = input_info
        metabolites = []
        for met_id in stoich.keys():
            met = model.metabolites.get_by_id(met_id)
            metabolites.append(met)
        class_dict = classify_metabolites_by_type(metabolites)

    exclude = ['adp','h2o','h','pi','ppi']
    contributions = dict()

    for met_type in class_dict.keys():
        contributions[met_type] = 0
        for met_id in class_dict[met_type]:
            
            if met_id.split('_')[0] not in exclude:
                met = model.metabolites.get_by_id(met_id)
                formula = cobra.core.formula.Formula(met.formula)
                met_weight = formula.weight
                if met_type == 'aminoacids':
                    met_weight -= 18 ## Taking out water after polymerization
                if met_id == 'atp_c':
                    met_stoich = stoich['atp_c'] + stoich['adp_c'] # Non-energy-related atp
                else:
                    met_stoich = stoich[met_id]
                met_mass = met_weight * met_stoich
                contributions[met_type] += abs(met_mass)/1000
    contributions

    return contributions