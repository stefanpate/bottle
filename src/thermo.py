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
from equilibrator_api.phased_reaction import PhasedReaction

from minedatabase.utils import get_compound_hash

pwd = Path(__file__).parent
EQ_URI = open(pwd / "../scripts/Database/eq_uris.txt").read().strip("\n")

lc = LocalCompoundCache()
lc.ccache = CompoundCache(sqlalchemy.create_engine(EQ_URI))

cc = ComponentContribution(ccache=lc.ccache)

def rxn_constructor_in_accession_IDs(reactants_list: list,
                                     products_list: list,
                                     eq_dict: dict):
    """
    Construct a reaction string of the form substrate + cofactor = product + cofactor
    except instead of SMILES, eQuilibrator accession IDs are used instead for all species
    this output can directly be fed into eQuilibrator to calculate reaction dG

    :param reactants_list (list): list of reactants on the LHS (includes cofactors)
    :param products_list (list): list of products on the RHS (includes cofactors)
    :return rxn_str_in_db_accessions (str): reaction string consisting of eQuilibrator accession IDs
    """
    # Replace with [eq_dict(reacant_id) for reactant_id in reactants_list]
    # Commented from Yash
    # reactant_cpd_objects = lc.get_compounds(reactants_list)  # create eQuilibrator compound objects for all LHS species
    # product_cpd_objects = lc.get_compounds(products_list)  # create eQuilibrator compound objects for all RHS species
    try:
        reactant_cpd_objects = [eq_dict[reactant_id] for reactant_id in reactants_list]
        product_cpd_objects = [eq_dict[product_id] for product_id in products_list]
        rxn_str_in_db_accessions = ''  # initialize empty string
    except KeyError:
        return "FAILED"

    # try:
    for reactant in reactant_cpd_objects:

        # get eQuilibrator accession IDs for all LHS species (if decomposable)
        reactant_accession = [identifier.accession.split("_")[1] for identifier in reactant.identifiers if identifier.registry_id == 2][0]
        rxn_str_in_db_accessions += reactant_accession
        rxn_str_in_db_accessions += ' + '

    # remove extra ' + ' sign on the right
    rxn_str_in_db_accessions = rxn_str_in_db_accessions.rstrip(' + ')

    # add ' = ' sign before switching over to products side
    rxn_str_in_db_accessions += ' = '

    for product in product_cpd_objects:

        # get eQuilibrator accession IDs for all RHS species (if decomposable)
        product_accession = [identifier.accession.split("_")[1] for identifier in product.identifiers if identifier.registry_id == 2][0]
        rxn_str_in_db_accessions += product_accession
        rxn_str_in_db_accessions += ' + '

    # remove extra ' + ' sign on the right
    rxn_str_in_db_accessions = rxn_str_in_db_accessions.rstrip(' + ')

    # return reaction string in terms of eQuilibrator accession IDs
    # The string looks something like C0001 + X0002 = C0002 + X0001
    return rxn_str_in_db_accessions

    # except Exception as e:
    #     print(e)
    #
    #     # if even one species on either LHS or RHS cannot be decomposed, return None
    #     return None

def calc_dG_frm_rxn_str(new_rxn_str: str, eq_dict: dict):
    """
    Calculate dG value and estimation errors at physiological and standard conditions

    :param new_rxn_str: reaction string in eQuilibrator accession IDs

    :return rxn_object: eQuilibrator reaction object
    :return phys_dG_value (float): physiological reaction dG value (1 mM concentrations)
    :return phys_dG_error (float): error in estimating dG value at physiological conditions
    :return std_dG_value (float): standard reaction dG value (1 M concentrations)
    :return std_dG_error (float): error in estimating dG value at standard conditions
    """

    # Old from yash
    # rxn_object = cc.parse_reaction_formula(new_rxn_str, )

    # Assume that everything coming in here is valid, i.e. all reaction cpds are in eq_dict
    rxn_object = PhasedReaction.parse_formula(lambda k: eq_dict[k], new_rxn_str)

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
    #TODO fix this name!
    reactants = [elt[1].pkid2smi[reactant_smi] for reactant_smi in reactants]
    products = [elt[1].pkid2smi[product_smi] for product_smi in products]
    rxn_accession = rxn_constructor_in_accession_IDs(reactants, products, eq_dict)
    if rxn_accession == "FAILED":
        return (rxn_id, "FAILED")
    
    eQ_rxn, phys_dG_value, phys_dG_error, std_dG_value, std_dG_error = calc_dG_frm_rxn_str(rxn_accession, eq_dict)
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

    import re
    pattern = "([^.>>]+)"
    cpd_ids = set()
    for rxn in pred_rxns.values():
        cpd_ids.update(rxn.pkid2smi.values())


    with open(paths_path, 'rb') as f:
        paths = pickle.load(f)    

    # Query to find Compounds with matching CompoundIdentifiers (registry_id is 2)
    eq_dict = dict()

    # Finds all cpd_id matching equilibrator objects.
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
            CompoundIdentifier.registry_id == 3,
            CompoundIdentifier.accession.in_(cpd_ids)
        )
    ).all()

    # This is key to not change database
    for compound in compounds_with_matching_identifiers:
        lc.ccache.session.expunge(compound)
        
    for cpd_entry in compounds_with_matching_identifiers:
        cpd_entry.identifiers = [entry for entry in cpd_entry.identifiers if entry.registry_id == 2]

    eq_dict = {cpd_entry.identifiers[0].accession.strip("pk_"): cpd_entry for cpd_entry in compounds_with_matching_identifiers}
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
        if rxn == "FAILED":
            print(f"{rxn[0]} failed to be constructed as an equilibrator PhasedReaction object. Exiting MDF.")

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

    # Load in your pickaxe object
    #  I think here it's actually reaction prediction objects

    # 2. Identify ALL compound ids (pk_###) and add to a set
    # 3. Use a sqlalchemy query
    #     Selects Compound objects
    #     Loads in eagerly (so we have loaded in memory, not as a DB connection)
    #     Loads microspecies, magnesium, etc. that we need to calculate dG things
    #     Filters results so you only get things in the cpd_ids set
    #     Detach these from the postgres instance so we don't delete
    #     Remove all identifiers except for the coco one with pk_### in it

    # 4. Generate eq_dict from these results
    #       {"C####": Compound, "X###": Compound}

    # 5. Using your rxn_preds, pass this and eq_dict to pool_f
    # 6. Use your code as was written, with two modifications
    #       a) Remove lc.get_compounds and use eq_dict instead (2-3x order of magnitude speed increase at best)
    #       b) Use PhasedReaction directly instead of parse reaction string. This takes a function to get values (labmda k: eq_dict[k]) as well as a reaction string, which was generated with pickaxe ids.

    # 7. Your code is exactly the same from here on out
