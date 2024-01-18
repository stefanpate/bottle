from multiprocessing import Process, Queue, Pool
from equilibrator_assets.local_compound_cache import LocalCompoundCache
from equilibrator_cache.compound_cache import CompoundCache
from equilibrator_api import ComponentContribution, Q_
import cvxpy
import numpy as np
from pathlib import Path
import sqlalchemy

from sqlalchemy.orm import joinedload, contains_eager, aliased, subqueryload
from equilibrator_cache.models import Compound, CompoundIdentifier

from minedatabase.utils import get_compound_hash

pwd = Path(__file__).parent
EQ_URI = open(pwd / "../scripts/Database/eq_uris.txt").read().strip("\n")

lc = LocalCompoundCache()
lc.ccache = CompoundCache(sqlalchemy.create_engine(EQ_URI))

cc = ComponentContribution(ccache=lc.ccache)

def rxn_constructor_in_accession_IDs(reactants_list: list,
                                     products_list: list):
    """
    Construct a reaction string of the form substrate + cofactor = product + cofactor
    except instead of SMILES, eQuilibrator accession IDs are used instead for all species
    this output can directly be fed into eQuilibrator to calculate reaction dG

    :param reactants_list (list): list of reactants on the LHS (includes cofactors)
    :param products_list (list): list of products on the RHS (includes cofactors)
    :return rxn_str_in_db_accessions (str): reaction string consisting of eQuilibrator accession IDs
    """
    # Replace with [eq_dict(reacant_id) for reactant_id in reactants_list]
    reactant_cpd_objects = lc.get_compounds(reactants_list)  # create eQuilibrator compound objects for all LHS species
    product_cpd_objects = lc.get_compounds(products_list)  # create eQuilibrator compound objects for all RHS species
    rxn_str_in_db_accessions = ''  # initialize empty string

    # try:
    for reactant in reactant_cpd_objects:

        # get eQuilibrator accession IDs for all LHS species (if decomposable)
        reactant_accession = reactant.get_accession()
        rxn_str_in_db_accessions += reactant_accession
        rxn_str_in_db_accessions += ' + '

    # remove extra ' + ' sign on the right
    rxn_str_in_db_accessions = rxn_str_in_db_accessions.rstrip(' + ')

    # add ' = ' sign before switching over to products side
    rxn_str_in_db_accessions += ' = '

    for product in product_cpd_objects:

        # get eQuilibrator accession IDs for all RHS species (if decomposable)
        product_accession = product.get_accession()
        rxn_str_in_db_accessions += product_accession
        rxn_str_in_db_accessions += ' + '

    # remove extra ' + ' sign on the right
    rxn_str_in_db_accessions = rxn_str_in_db_accessions.rstrip(' + ')

    # return reaction string in terms of eQuilibrator accession IDs
    return rxn_str_in_db_accessions

    # except Exception as e:
    #     print(e)
    #
    #     # if even one species on either LHS or RHS cannot be decomposed, return None
    #     return None

def calc_dG_frm_rxn_str(new_rxn_str: str):
    """
    Calculate dG value and estimation errors at physiological and standard conditions

    :param new_rxn_str: reaction string in eQuilibrator accession IDs

    :return rxn_object: eQuilibrator reaction object
    :return phys_dG_value (float): physiological reaction dG value (1 mM concentrations)
    :return phys_dG_error (float): error in estimating dG value at physiological conditions
    :return std_dG_value (float): standard reaction dG value (1 M concentrations)
    :return std_dG_error (float): error in estimating dG value at standard conditions
    """

    rxn_object = cc.parse_reaction_formula(new_rxn_str, )

    PhasedReaction.parse_formula(self.ccache.get_compound, formula)

    cc.p_h = Q_(7.4)
    cc.p_mg = Q_(3.0)
    cc.ionic_strength = Q_("0.25M")
    cc.temperature = Q_("298.15K")

    phys_dG_value = float(str(cc.physiological_dg_prime(rxn_object).value).rstrip(' kilojoule/mole'))
    phys_dG_error = float(str(cc.physiological_dg_prime(rxn_object).error).rstrip(' kilojoule/mole'))
    std_dG_value = float(str(cc.standard_dg_prime(rxn_object).value).rstrip(' kilojoule/mole'))
    std_dG_error = float(str(cc.standard_dg_prime(rxn_object).error).rstrip(' kilojoule/mole'))

    return rxn_object, phys_dG_value, phys_dG_error, std_dG_value, std_dG_error


def get_mdf_frm_rxn_objects(pathway_rxn_objects_list: list, lb: float, ub: float):

    """
    Calculating pathway mdf from eQuilibrator reaction objects

    :param pathway_rxn_objects_list: list of eQuilibrator reaction objects
    :param lb: lower bound on metabolite concentration
    :param ub: upper bound on metabolite concentration

    :return max_df ( float or None ):
    :return best_case_rxn_energies ( list or None ):
    :return best_case_concentrations ( list or None ):
    :return eq_compound_ids ( list or None ):
    """

    standard_dgr_prime_mean, standard_dgr_Q = cc.standard_dg_prime_multi(pathway_rxn_objects_list,
                                                                         uncertainty_representation="fullrank")

    S = cc.create_stoichiometric_matrix_from_reaction_objects(pathway_rxn_objects_list)
    Nc, Nr = S.shape

    ln_conc = cvxpy.Variable(shape=Nc, name="metabolite log concentration")  # vector
    B = cvxpy.Variable()  # scalar
    dg_prime = -(standard_dgr_prime_mean.m_as("kJ/mol") + cc.RT.m_as("kJ/mol") * S.values.T @ ln_conc)

    constraints = [
        np.log(np.ones(Nc) * lb) <= ln_conc,  # lower bound on concentrations
        ln_conc <= np.log(np.ones(Nc) * ub),  # upper bound on concentrations
        np.ones(Nr) * B <= dg_prime]

    prob_max = cvxpy.Problem(cvxpy.Maximize(B), constraints)
    prob_max.solve()
    max_df = prob_max.value

    best_case_rxn_energies = list(-dg_prime.value)
    best_case_concentrations = list(np.exp(ln_conc.value))
    eq_compound_ids = list(S.index)
    return max_df, best_case_rxn_energies, best_case_concentrations, eq_compound_ids

def calc_pathway_mdf(qin: Queue,
                     qout: Queue,
                     lb: float,
                     ub: float,
                     pred_rxns):
    """

    :param q: Shared multiprocessing queue
    :param lb: lower bound for metabolite concentration
    :param ub: upper bound for metabolite concentration

    :return:
    """

    while not qin.empty():
        path = qin.get()

        pathway_rxn_objects_list = []
        pathway_phys_dG_errors = []
        pathway_std_dG_errors = []
        for elt in path.rhashes:
            rxn = pred_rxns[elt]
            pathway_rxn_objects_list.append(rxn.eQ_rxn)
            pathway_phys_dG_errors.append(rxn.dG_phys[1])
            pathway_std_dG_errors.append(rxn.dG_std[1])

        # if dG estimation errors are within reason
        if all( abs(i) < 500 for i in pathway_phys_dG_errors ) and all( abs(j) < 500 for j in pathway_std_dG_errors):

            # then calculate mdf for pathway
            max_df, best_case_rxn_energies, \
            best_case_concentrations, eq_compound_ids = get_mdf_frm_rxn_objects(pathway_rxn_objects_list,
                                                                                               lb, ub)

        else:
            max_df, best_case_rxn_energies, best_case_concentrations, eq_compound_ids = None, None, None, None

        path.mdf = max_df
        path.dG_opt = best_case_rxn_energies
        path.conc_opt = best_case_concentrations

        qout.put(path)
        print(path, path.mdf)


def rxn_thermo_calcs(qin, qout):
    while not qin.empty():
        rxn_id, rxn = qin.get()                
        reactants, products = [side.split('.') for side in rxn.smarts.split('>>')]
        rxn_accession = rxn_constructor_in_accession_IDs(reactants, products)
        eQ_rxn, phys_dG_value, phys_dG_error, std_dG_value, std_dG_error = calc_dG_frm_rxn_str(rxn_accession)
        rxn.dG_std = (std_dG_value, std_dG_error)
        rxn.dG_phys = (phys_dG_value, phys_dG_error)
        rxn.eQ_rxn = eQ_rxn
        qout.put((rxn_id, rxn))
        print(rxn_id, qin.qsize())


def pool_f(elt):
    rxn_id, rxn = elt                
    print(rxn_id)
    reactants, products = [side.split('.') for side in rxn.smarts.split('>>')]
    rxn_accession = rxn_constructor_in_accession_IDs(reactants, products)
    eQ_rxn, phys_dG_value, phys_dG_error, std_dG_value, std_dG_error = calc_dG_frm_rxn_str(rxn_accession)
    rxn.dG_std = (std_dG_value, std_dG_error)
    rxn.dG_phys = (phys_dG_value, phys_dG_error)
    rxn.eQ_rxn = eQ_rxn
    return (rxn_id, rxn)

if __name__ == '__main__':
    import pickle
    import time
    '''
    Params
    '''
    starters = 'succinate'
    targets = 'mvacid'
    st_pair = ('succinate', 'mvacid')
    generations = 4
    n_workers = 2 # Cores to use
    ################################

    expansion_dir = '../data/processed_expansions/'
    fn = f"{starters}_to_{targets}_gen_{generations}_tan_sample_1_n_samples_1000.pk" # Expansion file name
    rxns_path = expansion_dir + 'predicted_reactions_' + fn
    paths_path = expansion_dir + 'paths_' + fn

    # Load reactions and paths
    with open(rxns_path, 'rb') as f:
        pred_rxns = pickle.load(f)

    with open(paths_path, 'rb') as f:
        paths = pickle.load(f)

    rxn_smarts = [rxn.smarts for rxn in pred_rxns.values()]
    cpd_smiles = set()
    for smarts in rxn_smarts:
        cpd_smiles.add(smarts)

    cpd_ids = [get_compound_hash(smiles) for smiles in cpd_smiles]

    # Query and make eq_dict
    eq_dict = dict()
    # Assuming cpd_chunk is a list of CompoundIdentifier accessions
    compounds_with_matching_identifiers = (
        lc.ccache.session.query(Compound)
        .options(
            joinedload(Compound.magnesium_dissociation_constants),
            joinedload(Compound.microspecies),
            subqueryload(Compound.identifiers).subqueryload(CompoundIdentifier.registry),
        )
        .join(Compound.identifiers)
        .join(CompoundIdentifier.registry)
        .filter(
            CompoundIdentifier.registry_id == 2,
            CompoundIdentifier.accession.in_(cpd_smiles)
        )
    ).all()

    for cpd_entry in compounds_with_matching_identifiers:
        cpd_entry.identifiers = [entry for entry in cpd_entry.identifiers if entry.registry_id == 2]

    eq_dict.update({cpd_entry.identifiers[0].accession.strip("pk_"): cpd_entry for cpd_entry in compounds_with_matching_identifiers})

    this_paths = paths[st_pair]

    # # initialize multiprocessing Queue
    # qin = Queue()
    # qout = Queue()

    # for k, v in list(pred_rxns.items())[:4]:
    #     qin.put((k, v))

    # workers = [Process(target=rxn_thermo_calcs, args=(qin, qout)) for i in range(n_workers)]

    # for elt in workers:
    #     elt.start()

    # for elt in workers:
    #     elt.join()

    # ctr = 0
    # while not qout.empty():
    #     rxn_id, rxn = qout.get()
    #     pred_rxns[rxn_id] = rxn
    #     ctr += 1
        
    # Pool
    # pool = Pool(processes=n_workers)
    # pred_rxns = pool.map(pool_f, pred_rxns.items())
        
    # In series
    start = time.perf_counter()
    pred_rxns = dict(map(pool_f, list(pred_rxns.items())[:100]))
    print(f"it took {time.perf_counter() - start}")

    # assert ctr == len(pred_rxns)

    # with open(rxns_path, 'wb') as f:
    #     pickle.dump(pred_rxns, f)
        
    # with open('../data/test_dgrs.pk', 'wb') as f:
    #     pickle.dump(pred_rxns, f)

    # print('saved')

    fake_path = list(pred_rxns.values())[:1]

            
    pathway_rxn_objects_list = []
    pathway_phys_dG_errors = []
    pathway_std_dG_errors = []
    lb = 1e-6 # 1 micro-mol
    ub = 500e-6 # 500 micro-mol

    for rxn in fake_path:
        # rxn = pred_rxns[elt]
        pathway_rxn_objects_list.append(rxn.eQ_rxn)
        pathway_phys_dG_errors.append(rxn.dG_phys[1])
        pathway_std_dG_errors.append(rxn.dG_std[1])

        # if dG estimation errors are within reason
        if all( abs(i) < 500 for i in pathway_phys_dG_errors ) and all( abs(j) < 500 for j in pathway_std_dG_errors):

            # then calculate mdf for pathway
            max_df, best_case_rxn_energies, \
            best_case_concentrations, eq_compound_ids = get_mdf_frm_rxn_objects(pathway_rxn_objects_list,
                                                                                                lb, ub)

        else:
            max_df, best_case_rxn_energies, best_case_concentrations, eq_compound_ids = None, None, None, None

    print('done')
    '''
    # rxns time
    10 14
    '''