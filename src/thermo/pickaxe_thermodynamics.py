# Code to generate eQuilibrator objects from Pickaxe objects
# Author: Kevin Shebek
# Date: 2024-01-21

from dataclasses import dataclass, field
from typing import Dict, List, NamedTuple, Optional, Union

import cvxpy
import numpy as np
import sqlalchemy
from equilibrator_api import Q_, ComponentContribution
from equilibrator_api.phased_reaction import PhasedReaction
from equilibrator_assets.local_compound_cache import LocalCompoundCache
from equilibrator_cache import Compound, CompoundIdentifier
from equilibrator_cache.models import Compound, CompoundIdentifier
from minedatabase.pickaxe import Pickaxe
from sqlalchemy.orm import joinedload, subqueryload


# make an error for a reaction not generated yet
def ReactionsNotGeneratedError(Exception, rxn_ids):
    print(f"The following reactions were not generated: {rxn_ids}")
    pass

@dataclass
class PathwayMDF:
    """
    Represents a pathway in thermodynamics analysis.

    Attributes:
        mdf_value (Optional[Q_]): The minimum driving force value.
        reaction_energies (Dict[str, Q_]): A dictionary of reaction energies.
        concentrations (Dict[str, Q_]): A dictionary of concentrations.
        status (str): The status of the pathway.
        missing_microspecies (List[str]): A list of missing microspecies.
        uncertainties (Dict[str, Q_]): A dictionary of uncertainties.
    """

    mdf_value: Optional[Q_] = None
    reaction_energies: Dict[str, Q_] = field(default_factory=dict)
    concentrations: Dict[str, Q_] = field(default_factory=dict)
    status: str = "INITIALIZED"
    missing_microspecies: List[str] = field(default_factory=list)
    uncertainties: Dict[str, Q_] = field(default_factory=dict)

class PickaxeThermodynamics:
    """
    A class representing the Pickaxe Thermodynamics object.

    Attributes:
        lc (LocalCompoundCache): The local compound cache object.
        ccache (CompoundCache): The compound cache object.
        session (sqlalchemy.session): The database session.
        ccont (ComponentContribution): The component contribution object.
        eQ_compound_dict (dict): The equilibrium compound dictionary.
        eQ_reaction_dict (dict): The equilibrium reaction dictionary.
        default_p_h (str): The default pH units.
        default_p_mg (str): The default pMg units.
        default_ionic_strength_units (str): The default ionic strength units.

    Methods:
        __init__(self, lc: LocalCompoundCache) -> None: Initialize the PickaxeThermodynamics object.
        generate_eQ_reaction_from_pickaxe(self, reaction_id: str, pk: Pickaxe) -> PhasedReaction: Generates an eQ_reaction object from a Pickaxe reaction.
        generate_eQ_compound_dict_from_pickaxe(self, pk: Pickaxe) -> Dict[str, Compound]: Generates an eQdict (equilibrium dictionary) from a Pickaxe object.
        generate_eQ_reaction_dict_from_pickaxe(self, pk: Pickaxe) -> Dict[str, PhasedReaction]: Generates an eQ_reaction dictionary from a Pickaxe object.
        query_database_for_pk_ids(self, pk_ids: List[str]) -> List[Compound]: Queries the database for compounds with matching identifiers.
        calculate_pathway_mdf(self, reaction_id_list: List[str], lb: float=1e-6, ub:float=1e-2) -> Dict[Union[Q_, None], Dict[str, Q_], Union[Dict[str, Q_], None], str, List[str], Union[Dict[str, Q_], None]]: Calculates the MDF (Max-Min Driving Force) for a given pickaxe ID pathway.
    """

    def __init__(self, lc: LocalCompoundCache) -> None:
        """
        Initialize the PickaxeThermodynamics object.

        Parameters:
        lc (LocalCompoundCache): The local compound cache object.

        Returns:
        None
        """
        self.lc = lc
        self.ccache = lc.ccache
        self.session = lc.ccache.session
        self.ccont = ComponentContribution(ccache=self.ccache)
        self.eQ_compound_dict = {}
        self.eQ_reaction_dict = {}

        self.default_p_h_units = ""
        self.default_p_mg_units = ""
        self.default_ionic_strength_units = "molar"

    @property
    def p_h(self):
        return self.ccont.p_h

    @p_h.setter
    def p_h(self, value: Union[Q_, float]):
        if not isinstance(value, Q_):
            value = Q_(value)
        self.ccont.p_h = value

    @property
    def p_mg(self):
        return self.ccont.p_mg

    @p_mg.setter
    def p_mg(self, value: Union[Q_, float]):
        if not isinstance(value, Q_):
            value = Q_(value)
        self.ccont.p_mg = value

    @property
    def ionic_strength(self):
        return self.ccont.ionic_strength

    @ionic_strength.setter
    def ionic_strength(self, value: Union[Q_, float]):
        if not isinstance(value, Q_):
            value = Q_(value, "molar")
        self.ccont.ionic_strength = value

    def _eQdict_getter(self, compound_id: str) -> Compound:
        """
        Getter function for the eQdict. Used to generate PhasedReaction objects.

        Args:
            compound_id (str): The compound identifier to get from the eQdict.

        Returns:
            Compound: The compound object from the eQdict.
        """
        return self.eQ_compound_dict[compound_id]

    def generate_eQ_reaction_from_pickaxe(
        self, reaction_id: str, pk: Pickaxe
    ) -> PhasedReaction:
        """
        Generates an eQ_reaction object from a Pickaxe reaction.

        Args:
            reaction_id (str): The ID of the reaction in the Pickaxe object.
            pk (Pickaxe): The Pickaxe object containing the reaction.
        """
        reaction = pk.reactions[reaction_id]

        # Generate the reaction string
        reactant_string = ""
        for reactant in reaction["Reactants"]:
            reactant_string += (reactant[1] + " + ") * reactant[0]
        # Remove the last "+"
        reactant_string = reactant_string[:-3]

        product_string = ""
        for product in reaction["Products"]:
            product_string += (product[1] + " + ") * product[0]
        # Remove the last "+"
        product_string = product_string[:-3]

        # Combine reactant and product strings
        reaction_string = reactant_string + " = " + product_string

        eQ_reaction = PhasedReaction.parse_formula(self._eQdict_getter, reaction_string)

        return eQ_reaction

    def generate_eQ_compound_dict_from_pickaxe(
        self, pk: Pickaxe
    ) -> Dict[str, Compound]:
        """
        Generates an eQdict (equilibrium dictionary) from a Pickaxe object and


        Args:
            pk (Pickaxe): The Pickaxe object containing the compounds.
        """
        compound_ids = [compound_id for compound_id in pk.compounds]
        compounds = self.query_database_for_pk_ids(compound_ids)

        # Disconnect the compounds from the database session and remove non coco identifiers
        for compound in compounds:
            self.session.expunge(compound)
            compound.identifiers = [
                entry for entry in compound.identifiers if entry.registry_id == 2
            ]

        eQ_compound_dict = {
            compound.identifiers[0].accession: compound for compound in compounds
        }

        self.eQ_compound_dict = eQ_compound_dict

    def generate_eQ_reaction_dict_from_pickaxe(
        self, pk: Pickaxe
    ) -> Dict[str, PhasedReaction]:
        for reaction_id in pk.reactions:
            try:
                eQreaction = self.generate_eQ_reaction_from_pickaxe(reaction_id, pk)
                self.eQ_reaction_dict[reaction_id] = eQreaction
            except BaseException as e:
                self.eQ_reaction_dict[reaction_id] = "FAILED"

    def query_database_for_pk_ids(self, pk_ids: List[str]) -> List[Compound]:
        """
        Queries the database for compounds with matching identifiers.

        Args:
            pk_ids (List[str]): List of pickaxe compound identifiers to query the database for.

        Returns:
            List[Compound]: List of compounds with matching identifiers.
        """
        compounds_with_matching_identifiers = (
            self.session.query(Compound)
            .options(
                joinedload(Compound.magnesium_dissociation_constants),
                joinedload(Compound.microspecies),
                subqueryload(Compound.identifiers).subqueryload(
                    CompoundIdentifier.registry
                ),
            )
            .join(Compound.identifiers)
            .join(CompoundIdentifier.registry)
            .filter(
                CompoundIdentifier.registry_id == 2,
                CompoundIdentifier.accession.in_(pk_ids),
            )
        ).all()

        return compounds_with_matching_identifiers

    def calculate_pathway_mdf(
        self, reaction_id_list: List[str], lb: float = 1e-6, ub: float = 1e-2
    ) -> PathwayMDF:
        """
        Calculates the MDF (Max-Min Driving Force) for a given pickaxe ID pathway.

        Args:
            reaction_id_list (List[str]): List of reaction IDs.
            lb (float, optional): Lower bound on concentrations. Defaults to 1e-6.
            ub (float, optional): Upper bound on concentrations. Defaults to 1e-2.

        Returns:
            Dict[Union[Q_, None], Dict[str, Q_], Union[Dict[str, Q_], None], str, List[str], Union[Dict[str, Q_], None]]:
            A dictionary containing the MDF value, reaction energies, concentrations, status, missing microspecies, and uncertainties.
        """
        # Define a pathway mdf results object
        pathway_mdf = PathwayMDF()

        # Get the eQ_reactions from the reaction_id_list
        pathway_eQ_reactions = [
            self.eQ_reaction_dict.get(r_id, "NONEXIST") for r_id in reaction_id_list
        ]

        # Ensure all reactions are generated
        non_exist_reactions = [
            r_id
            for r_id in reaction_id_list
            if isinstance(self.eQ_reaction_dict.get(r_id, "NONEXIST"), str)
            and self.eQ_reaction_dict.get(r_id, "NONEXIST") == "NONEXIST"
        ]
        if non_exist_reactions:
            raise ReactionsNotGeneratedError(non_exist_reactions)

        # find any self.eQ_reaction_dict[r_id] is a string instance and return if failed
        failed_reactions = {
            r_id: "FAILED"
            for r_id in reaction_id_list
            if isinstance(self.eQ_reaction_dict[r_id], str)
        }

        if failed_reactions:
            pathway_mdf.status = "FAILED"
            return pathway_mdf

        # Calculate the standard_dg_prime and standard_dg_prime_uncertainty for MDF setup
        (
            standard_dgr_prime,
            standard_dgr_uncertainty,
        ) = self.ccont.standard_dg_prime_multi(
            pathway_eQ_reactions, uncertainty_representation="fullrank"
        )

        S = self.ccont.create_stoichiometric_matrix_from_reaction_objects(
            pathway_eQ_reactions
        )

        # See if any missing microspecies, means ChemAxon bypassed
        missing_microspecies = []
        for compound in S.index:
            if not compound.microspecies:
                missing_microspecies.append(
                    [
                        identifier.accession
                        for identifier in compound.identifiers
                        if identifier.registry_id == 2
                        and identifier.accession[0] in ["X", "C", "T"]
                    ][0]
                )

        # Setup the MDF problem
        Nc, Nr = S.shape
        RT = self.ccont.RT

        ln_conc = cvxpy.Variable(
            shape=Nc, name="metabolite log concentration"
        )  # vector
        B = cvxpy.Variable()  # scalar
        dg_prime = -(
            standard_dgr_prime.m_as("kJ/mol") + RT.m_as("kJ/mol") * S.values.T @ ln_conc
        )

        # constraints = [
        #     np.log(np.ones(Nc) * lb) <= ln_conc,  # lower bound on concentrations
        #     ln_conc <= np.log(np.ones(Nc) * ub),  # upper bound on concentrations
        #     np.ones(Nr) * B <= dg_prime,
        # ]
        constraints = self.pick_constraints_for_MDF(S, Nc, Nr, ln_conc, dg_prime, B, lb, ub)

        # Solve the MDF problem
        prob_max = cvxpy.Problem(cvxpy.Maximize(B), constraints)
        prob_max.solve()
        pathway_mdf.mdf_value = prob_max.value
        pathway_mdf.status = prob_max.status

        # Convert the results to the correct units
        best_case_rxn_energies = [
            Q_(val, "kilojoule/mole") for val in list(-dg_prime.value)
        ]
        best_case_concentrations = [
            Q_(val, "mole") for val in list(np.exp(ln_conc.value))
        ]
        pathway_mdf.concentrations = {
            [
                identifier.accession
                for identifier in compound.identifiers
                if identifier.registry_id == 2
                and identifier.accession[0] in ["X", "C", "T"]
            ][0]: conc
            for compound, conc in zip(list(S.index), best_case_concentrations)
        }
        pathway_mdf.reaction_energies = {
            rxn_id: dg for rxn_id, dg in zip(reaction_id_list, best_case_rxn_energies)
        }
        # define pathway_mdf unc
        pathway_mdf.uncertainties = {
            rxn_id: uncertainty[0]
            for rxn_id, uncertainty in zip(reaction_id_list, standard_dgr_uncertainty)
        }

        return pathway_mdf

    def pick_constraints_for_MDF(
            self, S: any, Nc: any, Nr: any,
            ln_conc: any, dg_prime: any, B: any, lb: float, ub: float
    ):
        '''
        AUTHOR: YASH CHAINANI
        '''
        
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