"""
Step 2 of the model building pipeline

Goal here is to get reaction thermodynamics using the eQuilibrator software
This was developed by Avi Flemholtz and Elad Noor and further modified by Kevin Shebek
eQuilibrator is used to calculate dG and MDF of each enzymatic reaction parsed in step 01

Calculating the minimum driving force (MDF) for each reaction is crucial
MDF of reported enzymatic reactions will be used in step 03 to decide if reaction is feasible
In calculating the MDF, attention is paid to cofactors involved and their concentration constraints
Other thermodynamic metrics calculated are standard dG and physiological dG (all 1mM)

Once thermodynamic metrics are obtained, the Mongo database from step 01 is updated
or if results were saved to .json, then .json file is updated

for now, pH = 7.4, p_mg = 3.0, ionic str = 0.25 M, T = 298.15 K
"""
from multiprocessing import Process, Queue
import numpy as np
import pandas as pd
import sqlalchemy
import pymongo
import cvxpy

x = 0
for i in range(0,10):
    x += 1

from equilibrator_assets.local_compound_cache import LocalCompoundCache
from equilibrator_api import ComponentContribution, Q_
from equilibrator_cache.compound_cache import CompoundCache

# ----------------------- Set up thermo computation parameters below -----------------------
mongo_client_ID = (
    "mongodb://ymc1840:yash_is_cool@minedatabase.ci.northwestern.edu:27017"
)

URI_EQ = "postgresql://yashchainani96:sjiclass1052009@localhost:5432/thermo_cpds"

num_processes = 1  # os.cpu_count()

num_rxns_to_start_at = 11280
num_rxns_to_stop_at = 12000

db = "YC_feasib"
col = "M24_all_other_rules"

# --------------------- Helper functions for calculating thermodynamics ---------------------
def parse_rxn(rxn_str: str):
    """
    Separates a reaction string into a list of reactants and of products
    Reaction str of the form substrate_smiles + cofactor_smiles = product_smiles + cofactor_smiles
    This will return [substrate, cofactor] and [product, cofactor]

    :param rxn_str: reaction string of the form above
    :return reactants_list (list): list of reactants involved, eg:[substrate, cofactor]
    :return products_list (list): list of products involved, eg: [product, cofactor]
    """
    reactants, products = rxn_str.split(" = ")
    reactants_list = reactants.split(" + ")
    products_list = products.split(" + ")
    return reactants_list, products_list


def rxn_constructor_in_accession_IDs(reactants_list: list, products_list: list, lc: any):
    """
    Construct a reaction string of the form substrate + cofactor = product + cofactor
    Except instead of SMILES, eQuilibrator accession IDs are used instead for all species
    This output can directly be fed into eQuilibrator to calculate reaction dG

    :param reactants_list: list of reactants on the LHS (includes cofactors)
    :param products_list: list of products on the RHS (includes cofactors)
    :param lc: eQuilibrator compounds' cache object
    :return rxn_str_in_db_accessions: reaction string in eQuilibrator accession IDs
    """

    # note that compound objects will be created even if compound cannot be decomposed

    # create eQuilibrator compound objects for all LHS species
    reactant_cpd_objects = lc.get_compounds(reactants_list)

    # create eQuilibrator compound objects for all RHS species
    product_cpd_objects = lc.get_compounds(products_list)

    # initialize empty string
    rxn_str_in_db_accessions = ""

    try:
        for reactant in reactant_cpd_objects:
            # get eQuilibrator accession IDs for all LHS species (if decomposable)
            # if compound is not decomposable into groups, code will error out
            reactant_accession = reactant.get_accession()
            rxn_str_in_db_accessions += reactant_accession
            rxn_str_in_db_accessions += " + "

        # remove extra ' + ' sign on the right
        rxn_str_in_db_accessions = rxn_str_in_db_accessions.rstrip(" + ")

        # add ' = ' sign before switching over to products side
        rxn_str_in_db_accessions += " = "

        for product in product_cpd_objects:
            # get eQuilibrator accession IDs for all RHS species (if decomposable)
            product_accession = product.get_accession()
            rxn_str_in_db_accessions += product_accession
            rxn_str_in_db_accessions += " + "

        # remove extra ' + ' sign on the right
        rxn_str_in_db_accessions = rxn_str_in_db_accessions.rstrip(" + ")

        # return reaction string in terms of eQuilibrator accession IDs
        return rxn_str_in_db_accessions

    except:
        # if even one species on either LHS or RHS cannot be decomposed, return None
        return None


def create_cpds_df_frm_rxn_chunk(rxn_chunk: list, lc: any):
    """
    Create a dataframe of all unique compounds participating in a chunk of N reactions
    Compound SMILES strings and their eQuilibrator accession IDs should be stored in this dataframe
    This dataframe will then be referenced when calculating reaction dG and mdf values

    :param rxn_chunk: chunk of N reactions handled by the current worker core/ process
    :param lc: eQuilibrator compounds' cache object
    :return compounds_df: dataframe of compound SMILES and their identifiers for all species in these reactions
    """

    # initialize empty lists to store compound information
    all_cpd_smi = []
    all_cpd_accession_id = []
    all_cpd_inchi = []
    all_cpd_id = []

    ### Start processing all compounds in this chunk of N reactions by extracting all chemical species
    for rxn_dict in rxn_chunk:

        rxn_str = rxn_dict["Reaction eq"]
        reactants_list, products_list = parse_rxn(rxn_str)

        ## for all reactants on reaction LHS
        for reactant_smi in reactants_list:

            # confirm first that reactant smiles has not already been queried
            if reactant_smi not in all_cpd_smi:

                # create eQuilibrator compound object for reactant using its smiles
                try:
                    log_df = pd.DataFrame()
                    reactant_obj = lc.get_compounds(reactant_smi, log_df=log_df) # io operation to read compounds db

                    # get compound identifiers for this eQuilibrator compound object
                    reactant_accession_id = reactant_obj.get_accession()
                    reactant_inchi = reactant_obj.inchi_key
                    reactant_id = reactant_obj.id

                    # store information about this reactant
                    all_cpd_smi.append(reactant_smi)
                    all_cpd_accession_id.append(reactant_accession_id)
                    all_cpd_inchi.append(reactant_inchi)
                    all_cpd_id.append(reactant_id)

                # OSError arises when [Mg] is present
                # e.g. '*Cc1c2c3n4c1=CC1=N5C(=Cc6c(C=C)c(C)c7n6[Mg]4~5~N4=C(C=3CC2=O)[C@@H](CCC(=O)O)[C@H](C)C4=C7)C(C)=C1CC(*)*'
                # for such compounds, save their smiles but save all other references as None
                except OSError:
                    all_cpd_smi.append(reactant_smi)
                    all_cpd_accession_id.append(None)
                    all_cpd_inchi.append(None)
                    all_cpd_id.append(None)

                # AttributeError arises if a reactant cannot be decomposed by eQuilibrator
                except AttributeError:
                    print(f"\n{reactant_smi} from reaction {rxn_str} could not be decomposed. Error log:")
                    print(log_df)

                    all_cpd_smi.append(reactant_smi)
                    all_cpd_accession_id.append(None)
                    all_cpd_inchi.append(None)
                    all_cpd_id.append(None)

        ## for all products on reaction RHS
        for product_smi in products_list:

            # confirm first that product smiles have not already been queried
            if product_smi not in all_cpd_smi:

                # create equilibrator compound object for product using its smiles
                try:
                    log_df = pd.DataFrame()
                    product_obj = lc.get_compounds(product_smi)

                    # get compound identifiers for this eQuilibrator compound object
                    product_accession_id = product_obj.get_accession()
                    product_inchi = product_obj.inchi_key
                    product_id = product_obj.id

                    # store information about this product
                    all_cpd_smi.append(product_smi)
                    all_cpd_accession_id.append(product_accession_id)
                    all_cpd_inchi.append(product_inchi)
                    all_cpd_id.append(product_id)

                # OSError arises when [Mg] is present
                # e.g. '*Cc1c2c3n4c1=CC1=N5C(=Cc6c(C=C)c(C)c7n6[Mg]4~5~N4=C(C=3CC2=O)[C@@H](CCC(=O)O)[C@H](C)C4=C7)C(C)=C1CC(*)*'
                # for such compounds, save their smiles but save all other references as None
                except OSError:
                    all_cpd_smi.append(product_smi)
                    all_cpd_accession_id.append(None)
                    all_cpd_inchi.append(None)
                    all_cpd_id.append(None)

                # AttributeError arises if a product cannot be decomposed by
                except AttributeError:
                    print(f"\n{product_smi} from reaction {rxn_str} could not be decomposed. Error log:")
                    print(log_df)

                    all_cpd_smi.append(product_smi)
                    all_cpd_accession_id.append(None)
                    all_cpd_inchi.append(None)
                    all_cpd_id.append(None)

    ### finish processing compounds - i.e. no more querying PostgreSQL and Chem-axon

    # create a lookup table to referencing compound identities via their SMILES
    compounds_df = pd.DataFrame({'SMILES': all_cpd_smi,
                                 'accession_id': all_cpd_accession_id,
                                 'inchi': all_cpd_inchi,
                                 'eQuilibrator_id': all_cpd_id})

    return compounds_df


def calc_dG_frm_rxn_str(new_rxn_str: str, cc:any):
    rxn_object = cc.parse_reaction_formula(new_rxn_str)

    cc.p_h = Q_(7.4)
    cc.p_mg = Q_(3.0)
    cc.ionic_strength = Q_("0.25M")
    cc.temperature = Q_("298.15K")

    phys_dG_value = float(
        str(cc.physiological_dg_prime(rxn_object).value).rstrip(" kilojoule/mole")
    )
    phys_dG_error = float(
        str(cc.physiological_dg_prime(rxn_object).error).rstrip(" kilojoule/mole")
    )
    std_dG_value = float(
        str(cc.standard_dg_prime(rxn_object).value).rstrip(" kilojoule/mole")
    )
    std_dG_error = float(
        str(cc.standard_dg_prime(rxn_object).error).rstrip(" kilojoule/mole")
    )

    return rxn_object, phys_dG_value, phys_dG_error, std_dG_value, std_dG_error


def pick_constraints_for_MDF(
    S: any, Nc: any, Nr: any, ln_conc: any, dg_prime: any, B: any, lb: float, ub: float
):
    ### Cofactors to track for specific constraints
    ATP_inchi_key = "ZKHQWZAMYRWXGA-UHFFFAOYSA-N"
    ADP_inchi_key = "XTWYTFMLZFPYCI-UHFFFAOYSA-N"
    AMP_inchi_key = "UDMBCSSLTHHNCD-UHFFFAOYSA-N"
    NADP_plus_inchi_key = "XJLXINKUBYWONI-UHFFFAOYSA-O"
    NADPH_inchi_key = "ACFIXJIJDZMPPO-UHFFFAOYSA-N"
    NAD_plus_inchi_key = "BAWFJGJZGIEFAR-UHFFFAOYSA-O"
    NADH_inchi_key = "BOPGDPNILDQYTO-UHFFFAOYSA-N"
    PI_inchi_key = "NBIIXXVUZAFLBC-UHFFFAOYSA-L"
    PPI_inchi_key = "XPPKVPWEQAFLFU-UHFFFAOYSA-K"
    CoA_inchi_key = "RGJOEKWQDUBAIZ-UHFFFAOYSA-N"
    NH3_inchi_key = "QGZKDVFQNNGYKY-UHFFFAOYSA-O"

    ### Store the inchi keys of all the compounds in the reaction/ pathway
    compound_inchi_keys_list = []
    for compound in list(S.index):
        compound_inchi_keys_list.append(compound.inchi_key)

    ### Note the various conditions

    # Rules 14, 15, 223, 224, 408, 409
    ATP_ADP_cofactors = [
        ATP_inchi_key in compound_inchi_keys_list,
        ADP_inchi_key in compound_inchi_keys_list,
    ]

    # Rules 49, 50, 356, 402,
    ATP_ADP_PI_cofactors = [
        ATP_inchi_key in compound_inchi_keys_list,
        ADP_inchi_key in compound_inchi_keys_list,
        PI_inchi_key in compound_inchi_keys_list,
    ]

    # Rules 170, 171
    ATP_ADP_PI_CoA_cofactors = [
        ATP_inchi_key in compound_inchi_keys_list,
        ADP_inchi_key in compound_inchi_keys_list,
        PI_inchi_key in compound_inchi_keys_list,
        CoA_inchi_key in compound_inchi_keys_list,
    ]

    ATP_AMP_cofactors = [
        ATP_inchi_key in compound_inchi_keys_list,
        AMP_inchi_key in compound_inchi_keys_list,
    ]

    # Rules 66, 67, 345, 346
    ATP_AMP_PPI_cofactors = [
        ATP_inchi_key in compound_inchi_keys_list,
        AMP_inchi_key in compound_inchi_keys_list,
        PPI_inchi_key in compound_inchi_keys_list,
    ]

    ATP_AMP_PPI_CoA_cofactors = [
        ATP_inchi_key in compound_inchi_keys_list,  # Rules 38, 39
        AMP_inchi_key in compound_inchi_keys_list,
        PPI_inchi_key in compound_inchi_keys_list,
        CoA_inchi_key in compound_inchi_keys_list,
    ]

    ATP_AMP_PPI_NH3_cofactors = [
        ATP_inchi_key in compound_inchi_keys_list,  # Rules 385, 386
        AMP_inchi_key in compound_inchi_keys_list,
        PPI_inchi_key in compound_inchi_keys_list,
        NH3_inchi_key in compound_inchi_keys_list,
    ]

    NADP_NADPH_cofactors = [
        NADP_plus_inchi_key in compound_inchi_keys_list,
        NADPH_inchi_key in compound_inchi_keys_list,
    ]

    NAD_NADH_cofactors = [
        NAD_plus_inchi_key in compound_inchi_keys_list,
        NADH_inchi_key in compound_inchi_keys_list,
    ]

    ### ATP, ADP: phosphate donor/ acceptor rxns (ATP/ ADP = 10)
    # Rules 14, 15
    if all(ATP_ADP_cofactors):
        for i, compound_inchi in enumerate(compound_inchi_keys_list):
            if compound_inchi == ATP_inchi_key:
                ATP_index = i

            if compound_inchi == ADP_inchi_key:
                ADP_index = i

        ATP_ADP_constraints = [
            np.log(np.ones(Nc) * lb) <= ln_conc,  # lower bound on concentrations
            ln_conc <= np.log(np.ones(Nc) * ub),  # upper bound on concentrations
            np.ones(Nr) * B <= dg_prime,
            (ln_conc[ATP_index] - ln_conc[ADP_index]) <= 2.303,
            2.300 <= (ln_conc[ATP_index] - ln_conc[ADP_index]),
        ]

        return ATP_ADP_constraints

    ### ATP, AMP, PPI, CoA: pyrophosphate donor/ acceptor rxns (ATP/ AMP = 10, [pyrophosphate] = 1 mM)
    # Rules 66, 67, 345, 346
    if all(ATP_AMP_PPI_cofactors):
        for i, compound_inchi in enumerate(compound_inchi_keys_list):
            if compound_inchi == ATP_inchi_key:
                ATP_index = i

            if compound_inchi == AMP_inchi_key:
                AMP_index = i

            if compound_inchi == PPI_inchi_key:
                PPI_index = i

        ATP_AMP_constraints = [
            np.log(np.ones(Nc) * lb) <= ln_conc,  # lower bound on concentrations
            ln_conc <= np.log(np.ones(Nc) * ub),  # upper bound on concentrations
            np.ones(Nr) * B <= dg_prime,
            (ln_conc[ATP_index] - ln_conc[AMP_index]) <= 2.303,
            2.303 <= (ln_conc[ATP_index] - ln_conc[AMP_index]),
            ln_conc[PPI_index] <= -4.6,
            -4.7 <= ln_conc[PPI_index],
        ]

        return ATP_AMP_constraints

    ### NADP+ - NADPH reactions (NADPH/ NADP+ = 10)
    # Rules 2,3..
    if all(NADP_NADPH_cofactors):
        for i, compound_inchi in enumerate(compound_inchi_keys_list):
            if compound_inchi == NADP_plus_inchi_key:
                NADP_plus_index = i

            if compound_inchi == NADPH_inchi_key:
                NADPH_index = i

        NADP_NADPH_constraints = [
            np.log(np.ones(Nc) * lb) <= ln_conc,  # lower bound on concentrations
            ln_conc <= np.log(np.ones(Nc) * ub),  # upper bound on concentrations
            np.ones(Nr) * B <= dg_prime,
            (ln_conc[NADPH_index] - ln_conc[NADP_plus_index]) <= 2.303,
            2.300 <= (ln_conc[NADPH_index] - ln_conc[NADP_plus_index]),
        ]

        return NADP_NADPH_constraints

    ### NAD+ - NADH reactions (NADH/ NAD+ = 0.1, i.e NAD+/ NADH = 10):
    # Rules 2,3,...
    if all(NAD_NADH_cofactors):
        for i, compound_inchi in enumerate(compound_inchi_keys_list):
            if NAD_plus_inchi_key == compound_inchi:
                NAD_plus_index = i

            if NADH_inchi_key == compound_inchi:
                NADH_index = i

        NAD_NADH_constraints = [
            np.log(np.ones(Nc) * lb) <= ln_conc,  # lower bound on concentrations
            ln_conc <= np.log(np.ones(Nc) * ub),  # upper bound on concentrations
            np.ones(Nr) * B <= dg_prime,
            (ln_conc[NAD_plus_index] - ln_conc[NADH_index]) <= 2.303,
            2.300 <= (ln_conc[NAD_plus_index] - ln_conc[NADH_index]),
        ]
        return NAD_NADH_constraints

    ### All others
    else:
        reg_constraints = [
            np.log(np.ones(Nc) * lb) <= ln_conc,  # lower bound on concentrations
            ln_conc <= np.log(np.ones(Nc) * ub),  # upper bound on concentrations
            np.ones(Nr) * B <= dg_prime,
        ]

    return reg_constraints


def calc_MDF(rxn_object, lb: float, ub: float, cc:any):
    # Setting up optimization problem with variables
    standard_dgr_prime_mean, standard_dgr_Q = cc.standard_dg_prime_multi(
        rxn_object, uncertainty_representation="fullrank"
    )

    S = cc.create_stoichiometric_matrix_from_reaction_objects(rxn_object)
    Nc, Nr = S.shape

    ln_conc = cvxpy.Variable(shape=Nc, name="metabolite log concentration")  # vector
    B = cvxpy.Variable()  # scalar
    dg_prime = -(
        standard_dgr_prime_mean.m_as("kJ/mol")
        + cc.RT.m_as("kJ/mol") * S.values.T @ ln_conc
    )

    constraints = pick_constraints_for_MDF(S, Nc, Nr, ln_conc, dg_prime, B, lb, ub)

    prob_max = cvxpy.Problem(cvxpy.Maximize(B), constraints)
    prob_max.solve()
    max_df = prob_max.value
    return max_df


def get_rxn_thermo(q: Queue, process_num: int, mongo_client_ID: str, db: str, col: str, URI_EQ:str):

    while not q.empty():

        # get the chunk of reactions
        rxn_chunk = q.get()

        # announce process commencement
        print(f"\nProcess {process_num+1} starting to process current chunk of {len(rxn_chunk)} reactions")

        # create new eQuilibrator connection to PostgreSQL for each process
        lc = LocalCompoundCache()
        lc.ccache = CompoundCache(sqlalchemy.create_engine(URI_EQ))
        cc = ComponentContribution(ccache=lc.ccache)

        # create new MongoDB connection for each process
        mongo_client = pymongo.MongoClient(mongo_client_ID)
        database_to_write = mongo_client[db]
        collection_to_write = database_to_write[col]

        compounds_df = create_cpds_df_frm_rxn_chunk(rxn_chunk, lc)

        # switch reaction strings from smiles to accession ids
        for rxn_dict in rxn_chunk:

            # will use MongoDB id to update document with thermo info
            mongo_rxn_id = rxn_dict["_id"]

            # get existing reaction equation string (of the form 'SMILES A + SMILES B = SMILES C + SMILES D)
            rxn_str = rxn_dict["Reaction eq"]

            # split existing reaction equation string into reactants list and products list
            reactants_list, products_list = parse_rxn(rxn_str)

            # select identifier type for querying compounds lookup table
            # pick between 'accession_id', 'inchi', 'eQuilibrator_id'
            identifier_type = 'accession_id'

            # initialize a skip_reaction flag and set it to False first
            skip_reaction = False

            # initailize a new reaction string
            new_rxn_str = ''

            for reactant_smi in reactants_list:
                reactant_id = list(compounds_df[compounds_df['SMILES'] == reactant_smi][identifier_type])[0]

                # if reactant smiles does not have a reference, then set skip_reaction flag to True
                if reactant_id is None:
                    skip_reaction = True
                    break

                new_rxn_str += reactant_id
                new_rxn_str += ' + '

            # remove extra ' + ' sign at the end of the last added reactant
            new_rxn_str = new_rxn_str.rstrip(' + ')

            # introduce a ' = ' before switching over from reactants to products side
            new_rxn_str += ' = '

            for product_smi in products_list:
                product_id = list(compounds_df[compounds_df['SMILES'] == product_smi][identifier_type])[0]

                if product_id is None:
                    skip_reaction = True
                    break

                new_rxn_str += product_id
                new_rxn_str += ' + '

            # remove extra ' + ' sign at the end of the last added product
            new_rxn_str = new_rxn_str.rstrip(' + ')

            if not skip_reaction:
                (
                    rxn_object,
                    phys_dG_value,
                    phys_dG_error,
                    std_dG_value,
                    std_dG_error,
                ) = calc_dG_frm_rxn_str(new_rxn_str, cc)

                # Calculate MDF for this reaction
                # eQuilibrator rxn object obtained above should be in a list for MDF
                # 0.1 millimolar to 100 millimolar
                lb = 1e-4
                ub = 1e-1
                min_DF = calc_MDF([rxn_object], lb, ub, cc)

                rxn_dict_update = {
                    "Physiological dG": phys_dG_value,
                    "Physiological dG error": phys_dG_error,
                    "Standard dG": std_dG_value,
                    "Standard dG error": std_dG_error,
                    "Minimum DF": min_DF,
                }

            else:
                rxn_dict_update = {
                    "Physiological dG": None,
                    "Physiological dG error": None,
                    "Standard dG": None,
                    "Standard dG error": None,
                    "Minimum DF": None,
                }

            collection_to_write.update_one(
                {"_id": mongo_rxn_id}, {"$set": rxn_dict_update}, upsert=False
            )

        # announce process termination
        print(f"\nProcess {process_num + 1} finished processing current chunk of {len(rxn_chunk)} reactions")
        # .cache.close
        # and then delete

# ---------------------------- Start thermodynamic computations with multiprocessing -----------------
if __name__ == '__main__':

    mongo_client = pymongo.MongoClient(mongo_client_ID)
    mongo_database = mongo_client[db]
    mongo_collection = mongo_database[col]
    data = mongo_collection.find({})

    rxns_to_process = []

    # extract all relevant reactions to process
    for i, rxn_dict in enumerate(data):
        if i >= num_rxns_to_start_at and i <= num_rxns_to_stop_at:
            rxns_to_process.append(rxn_dict)
        else:
            pass

    chunk_size = 5
    print(f"\nChunk size selected to process reactions is {chunk_size}")

    # initialize shared multiprocessing queue to store all reactions
    q = Queue()

    # create reaction chunks once all relevant reactions have been extracted
    all_rxn_chunks = [ rxns_to_process[i: i + chunk_size ] for i in range(0, len(rxns_to_process), chunk_size) ]

    # add each chunk of N reactions into the shared multiprocessing queue
    for rxn_chunk in all_rxn_chunks:
        q.put(rxn_chunk)

    # initialize empty list to store multiprocessing objects
    processes = []

    for i in range(num_processes):

        p = Process(
            target=get_rxn_thermo,
            args=(
                q,
                i,
                mongo_client_ID,
                db,
                col,
                URI_EQ
            ),
        )
        p.start()
        processes.append(p)

    for p in processes:
        p.join()
