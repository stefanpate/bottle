import hashlib
import re
from collections import defaultdict
from typing import Dict, List

import networkx
import networkx as nx
from rdkit.Chem import CanonSmiles

def path_to_generation(path, DG):
    rxns = [DG.nodes()[rxn_id] for rxn_id in path]
    rxn_types = [rxn["Type"] for rxn in rxns]

    if rxn_types == ["Biology"]:
        return "B1"
    elif rxn_types == ["Biology", "Biology"]:
        return "B2"
    elif rxn_types == ["Chemistry"]:
        return "CB"
    elif rxn_types == ["Biology", "Chemistry"]:
        return "C1"
    elif rxn_types == ["Biology", "Biology", "Chemistry"]:
        return "C2"
    else:
        return "None"

class Pathways:
    def __init__(self, target_cid, pathways, DG, feas_dict):
        self.target_id = target_cid
        self.pathways = [Pathway(target_cid, pathway, DG, feas_dict) for pathway in pathways]

class Pathway:
    def __init__(self, target_cid: str, pathway: List[str], DG: networkx.DiGraph, feas_dict: Dict[str, bool]) -> None:
        # Basic Information
        self._id = hashlib.md5(",".join(pathway).encode()).hexdigest()
        self.pathway = pathway
        self.length = len(pathway)
        self.terminal_generation = path_to_generation(pathway, DG)
        self.DG = DG
        
        # Reaction Information        
        self.reaction_types = {rxn_id: DG.nodes()[rxn_id]["Type"] for rxn_id in pathway}
        self.reaction_stoichs = {rxn_id: DG.nodes()[rxn_id]["Stoich"] for rxn_id in pathway}

        # Starting/Ending/Targets
        self.target_cid = target_cid
        self.sinks = [smiles for smiles, s in DG.nodes()[self.pathway[-1]]["Stoich"].items() if s > 0]
        self.sources = [smiles for smiles, s in DG.nodes()[self.pathway[0]]["Stoich"].items() if s < 0]

        cpd_sets = set()
        bio_sets = set()
        for rxn_id, rxn_stoich in self.reaction_stoichs.items():
            for cpd in rxn_stoich:
                cpd_sets.add(cpd)
                if self.reaction_types[rxn_id] == "Biology":
                    bio_sets.add(cpd)

        self.compounds = list(cpd_sets)
        self.bio_compounds = list(bio_sets)

        self.compound_ids = set()
        self.bio_compound_ids = set()
        for rxn in self.pathway:
            reactants = [rct[1] for rct in DG.nodes()[rxn]["Reactants"]]
            products = [prd[1] for prd in DG.nodes()[rxn]["Products"]]

            self.compound_ids.update(reactants)
            self.compound_ids.update(products)

            if self.reaction_types[rxn] == "Biology":
                self.bio_compound_ids.update(reactants)
                self.bio_compound_ids.update(products)
        
        # Feasibility Metrics
        self.reaction_feas = []
        for rxn_id in pathway:
            if DG.nodes()[rxn_id]["Type"] == "Biology":
                self.reaction_feas.append(feas_dict[rxn_id])
            else:
                self.reaction_feas.append("Chemistry")

        self.MDF = -1000
    

# class Pathways:
#     def __init__(self, pathway_data):
#         """
#         pathway data is a dictionary containing at least:
#             paths_rxn_only : dict of lists
#             DG : nx.DiGraph 
#         """
#         # Define variables first
#         self.DG = pathway_data["DG"]
#         self.pathways = [Pathway(path, self.DG) for paths in pathway_data["paths_rxn_only"].values() for path in paths]
#         self.pathways_dict = {pathway._id: pathway for pathway in self.pathways}

#         self.target_pathways = defaultdict(list)
#         self.source_pathways = defaultdict(list)
#         self.compounds = set()
#         self.reactions = set()
#         self._populate_pathway_info()
        
#     def _populate_pathway_info(self):
#         # populate the target and source dictionaries
#         for pathway in self.pathways:
#             for source in pathway.sources:
#                 # Check to see if pathway is in dict already
#                 if pathway not in self.source_pathways:
#                     self.source_pathways[source].append(pathway)
            
#             for target in pathway.sinks:
#                 if target not in self.target_pathways:
#                     self.target_pathways[target].append(pathway)
        
#         # populate unique compounds and reactions
#         # Redo just using DG
#         for pathway in self.pathways:
#             for reaction in pathway.reactions:
#                 self.reactions.add(reaction)
#                 for compound in self.DG.nodes()[reaction]["stoich"]:
#                     self.compounds.add(compound)
                    
#         self.compounds = tuple(self.compounds)
#         self.reactions = tuple(self.reactions)

def get_canon_smiles(smi):
    try:
        return CanonSmiles(smi)
    except BaseException as e:
        return smi
def create_graph(df_cpds, df_rxns):
    # Generate a directed bipartite graph
    # 1. Add Compound Nodes
    # 2. Add Reaction Nodes
    # 3. Add directed edges
    cpd_node_list = []
    rxn_node_list = []
    edge_list = []

    starting_nodes = []

    # Get compound information
    for i, cpd_row in df_cpds.iterrows():
        cpd_node_list.append(
            (
                get_canon_smiles(cpd_row.SMILES),
                {
                    "InChIKey": cpd_row.InChIKey,
                    "Type": cpd_row.Type,
                    "Generation": cpd_row.Generation,
                }
            )
        )

        if cpd_row.Type == "Starting Compound":
            starting_nodes.append(get_canon_smiles(cpd_row.SMILES))

    # Get reaction information
    for i, rxn_row in df_rxns.iterrows():
        stoich = get_stoich(rxn_row["SMILES equation"])
        rxn_node_list.append(
            (
                rxn_row["Rxn hash"],
                {
                    "Rule": rxn_row["Reaction rules"],
                    "Type": "reaction",
                    "Stoich": stoich,
                    "feasible": None
                }
            )
        )

        for smi, stoi in stoich.items():
            if stoi < 0:
                edge_list.append((get_canon_smiles(smi), rxn_row["Rxn hash"]))
            else:
                edge_list.append((rxn_row["Rxn hash"], get_canon_smiles(smi)))
        
    # Create Graph
    DG = nx.DiGraph()  
    DG.add_nodes_from(cpd_node_list, bipartite=0)
    DG.add_nodes_from(rxn_node_list, bipartite=1)
    DG.add_edges_from(edge_list)
    
    return DG


def get_stoich(smiles_reaction):
    reactants, products = smiles_reaction.split("=>")
    reactants = [pair.strip(" ").split(" ") for pair in reactants.split(" + ")]
    reactants = {get_canon_smiles(pair[1]): -1*int(pair[0].strip("(").strip(")")) for pair in reactants}
    products = [pair.strip(" ").split(" ") for pair in products.split(" + ")]
    products = {get_canon_smiles(pair[1]): int(pair[0].strip("(").strip(")")) for pair in products}

    stoich = {**reactants, **products}
    return stoich   


def clean_up_graph(G, i_max=10):
    c_2, r_2 = get_cpd_rxns(G)
    c_1, r_1 = (0, 0)

    completed = True
    i = 0
    while (c_1, r_1) != (c_2, r_2):
        if i >= i_max:
            completed = False
        else:
            i += 1

        G = remove_orphans(G)
        G = remove_imbalanced_reactions(G)
        G = remove_disconnected_compound(G)

        c_1, r_1 = c_2, r_2
        c_2, r_2 = get_cpd_rxns(G)
        
    return G

def remove_orphans(G):
    """Check for orphan nodes and remove them"""
    orphans = []
    for node in G.nodes():
        if G.degree(node) == 0:
            orphans.append(node)
    G.remove_nodes_from(orphans)
    return G

def remove_imbalanced_reactions(G):
    """Check for reactions that aren't balanced, i.e. stoich doesn't match up"""
    rxns_to_remove = []
    for node in G.nodes():
        if node.startswith("R"):
            stoich = G.nodes()[node]["stoich"]
            reactants = [i for i in stoich.values() if i < 0]
            products = [i for i in stoich.values() if i > 0]
            
            in_degree = G.in_degree(node)
            out_degree = G.out_degree(node)
            
            if (abs(sum(reactants)) != in_degree) or (sum(products) != out_degree):
                rxns_to_remove.append(node)
                
    G.remove_nodes_from(rxns_to_remove)
    return G

def remove_disconnected_compound(G):
    """Check for compounds with no input rxns and gen > 1"""
    disconnected = []
    for node in G.nodes():
        if G.nodes()[node]["Type"] != "reaction":
            if (G.in_degree(node) == 0) & (G.nodes()[node]["Type"] != "Starting Compound"):
                disconnected.append(node)
    G.remove_nodes_from(disconnected)
    return G

def get_paths_to_targets(G, starting_node, target_nodes, max_depth, remove_starting_targets=True):
    def depth_first(G, current_node, target_nodes, visited, depth):
        # Are we too deep?
        if max_depth <= depth:
            return

        # Record current positon and decide if we should go further
        visited.append(current_node)
        for out_edge in G.out_edges(current_node):
            # Have we reached our target?
            if out_edge[1] in target_nodes:
                visited.append(out_edge[1])
                found_paths.append(visited)
                return

            # Is the connection a cofactor? If so, terminate here
            elif G.nodes()[out_edge[1]]["Type"] == "Coreactant":
                return
            # If not, have we seen this node before?
            elif out_edge[1] in visited:
                return
            # Can also add in more constraints here, e.g. thermo, type, etc.
            # Finally continue traversal
            else:
                depth_first(G, out_edge[1], target_nodes, visited.copy(), depth+1)
    
    # if remove_starting_targets:
    #     target_nodes = [target for target in target_nodes if starting_node != target_nodes]
    #     if not target_nodes:
    #         return []

    found_paths = []
    depth_first(G, starting_node, target_nodes, [], 0)

    return found_paths

def get_reverse_paths_to_starting(G, begin_node, end_nodes, max_depth):
    """Get the reverse pathways from target comopunds to starting compounds.


    Uses a breadth first search to identify pathways from a specified start of compounds
    to targets. This method does a reverse search.

    Parameters
    ----------
    G : nx.DiGraph
        The DiGraph containing network information.
    begin_node : str
        The pickaxe cid of the beginning node
    end_nodes : List[str]
        A list of pickaxe cids for the desired end nodes
    max_depth : int
        The maximum depth to search.
    """
    def depth_first_reversed(G, current_node, target_nodes, visited, depth):
        # Are we too deep?
        if max_depth <= depth:
            return

        # Record current positon and decide if we should go further
        visited.append(current_node)
        for in_edge in G.in_edges(current_node):
            # Have we reached our target?
            if in_edge[0] in target_nodes:
                visited.append(in_edge[0])
                found_paths.append(visited)
                return

            # Is the connection a cofactor? If so, terminate here
            elif G.nodes()[in_edge[0]]["Type"] == "Coreactant":
                continue
            # If not, have we seen this node before?
            elif in_edge[0] in visited:
                return
            # Can also add in more constraints here, e.g. thermo, Type, etc.
            # Finally continue traversal
            else:
                depth_first_reversed(G, in_edge[0], target_nodes, visited.copy(), depth+1)

    found_paths = []
    depth_first_reversed(G, begin_node, end_nodes, [], 0)

    return found_paths

def get_cpd_rxns(G):
    compound_nodes = set()

    for n, d in G.nodes(data=True):
        try:
            if d["bipartite"] == 0:
                compound_nodes.add(n)
        except:
            continue

    reaction_nodes = G.nodes() - compound_nodes
    return compound_nodes, reaction_nodes

#########################
## To Remove
def remove_bad_compounds(G):
    # Remove the following:
    #  1. If cpd and gen != 0 and in_deg == 0
    to_be_removed = set()
    for node in G.nodes():
        if "Generation" in G.nodes()[node]:
            if (
                (
                    G.nodes()[node]["Generation"] != 0
                    & G.in_degree(node) == 0
                )
                or (G.degree(node) == 0)
            ):
                to_be_removed.add(node)

def check_degree_matches(G, rxn):
    rxn_node = G.nodes()[rxn]
    in_degree = len([i for i in rxn_node["stoich"].values() if i <= 0])
    out_degree = len([i for i in rxn_node["stoich"].values() if i >= 0])
    
    if (G.in_degree(rxn) == in_degree) & (G.out_degree(rxn) == out_degree):
        return True
    else:
        return False
    
def create_graph_from_pickaxe(pk, rxn_type):
    # Generate a directed bipartite graph
    # 1. Add Compound Nodes
    # 2. Add Reaction Nodes
    # 3. Add directed edges
    cpd_node_list = []
    rxn_node_list = []
    edge_list = []

    starting_nodes = []
    smiles_to_cid = {}

    # Get compound information
    for i, cpd in pk.compounds.items(): 
        cpd_node_list.append(
            (
                i,
                {
                    "SMILES": cpd["SMILES"],
                    "InChIKey": cpd["InChI_key"],
                    "Type": cpd["Type"],
                    "Generation": cpd["Generation"],
                    "_id": cpd["_id"]
                }
            )
        )

        if cpd["Type"] == "Starting Compound":
            starting_nodes.append(i)
        
        smiles_to_cid[cpd["SMILES"]] = i

    # Get reaction information
    for i, rxn in pk.reactions.items():
        stoich = get_stoich_pk(i, pk)
        rxn_node_list.append(
            (
                i,
                {
                    "Rule": rxn["Operators"],
                    "Type": rxn_type,
                    "Stoich": stoich,
                    "feasible": None,
                    "Reactants": rxn["Reactants"],
                    "Products": rxn["Products"]
                }
            )
        )
        
        for _, c_id in rxn["Reactants"]:
            edge_list.append((c_id, i))
        for _, c_id in rxn["Products"]:
            edge_list.append((i, c_id))
    
    # Create Graph
    DG = nx.DiGraph(smiles_to_cid=smiles_to_cid)  
    DG.add_nodes_from(cpd_node_list, bipartite=0)
    DG.add_nodes_from(rxn_node_list, bipartite=1)
    DG.add_edges_from(edge_list)
    
    return DG, rxn_node_list, edge_list

def add_nodes_from_pickaxe(pk, DG, rxn_type):
    # Add to a directed bipartite graph
    # 1. Add Compound Nodes
    # 2. Add Reaction Nodes
    # 3. Add directed edges
    cpd_node_list = []
    rxn_node_list = []
    edge_list = []

    starting_nodes = []
    smiles_to_cid = {}

    print("Getting Node Info")
    # Get compound information
    for i, cpd in pk.compounds.items():
        if cpd["_id"] not in DG:
            if cpd["Type"] == "Starting":
                cpd["Type"] = "Starting Compound"
            cpd_node_list.append(
                (
                    i,
                    {
                        "SMILES": cpd["SMILES"],
                        # "InChIKey": cpd["InChI_key"],
                        "Type": cpd["Type"] if cpd["Type"] != "Starting Compound" else "Predicted",
                        "_id": cpd["_id"]
                    }
                )
            )
        
            smiles_to_cid[cpd["SMILES"]] = i

    # Get reaction information
    for i, rxn in pk.reactions.items():
        if rxn not in DG:
            stoich = get_stoich_pk(i, pk)
            rxn_node_list.append(
                (
                    i,
                    {
                        "Rule": rxn["Operators"],
                        "Type": rxn_type,
                        "Stoich": stoich,
                        "feasible": None,
                        "Reactants": rxn["Reactants"],
                        "Products": rxn["Products"]
                    }
                )
            )

            for _, c_id in rxn["Reactants"]:
                edge_list.append((c_id, i))
            for _, c_id in rxn["Products"]:
                edge_list.append((i, c_id))
    
    print("Adding Nodes")
    # Create Graph  
    DG.add_nodes_from(cpd_node_list, bipartite=0)
    DG.add_nodes_from(rxn_node_list, bipartite=1)
    DG.add_edges_from(edge_list)
    
    DG.graph["smiles_to_cid"] = {**DG.graph["smiles_to_cid"], **smiles_to_cid}
    
    return DG, rxn_node_list, edge_list

# def create_graph_from_pickaxe(pk):
#     # Generate a directed bipartite graph
#     # 1. Add Compound Nodes
#     # 2. Add Reaction Nodes
#     # 3. Add directed edges
#     cpd_node_list = []
#     rxn_node_list = []
#     edge_list = []

#     starting_nodes = []

#     # Get compound information
#     for i, cpd in pk.compounds.items(): 
#         cpd_node_list.append(
#             (
#                 get_canon_smiles(cpd["SMILES"]),
#                 {
#                     "InChIKey": cpd["InChI_key"],
#                     "type": cpd["Type"],
#                     "Generation": cpd["Generation"],
#                     "_id": cpd["_id"]
#                 }
#             )
#         )
        
#         if 'CC1=C(C(O)CO)C=C1' == get_canon_smiles(cpd["SMILES"]):
#             print(cpd)

#         if cpd["Type"] == "Starting Compound":
#             starting_nodes.append(get_canon_smiles(cpd["SMILES"]))

#     # Get reaction information
#     for i, rxn in pk.reactions.items():
#         stoich = get_stoich_pk(i, pk)
#         rxn_node_list.append(
#             (
#                 i,
#                 {
#                     "Rule": rxn["Operators"],
#                     "type": "reaction",
#                     "Stoich": stoich,
#                     "feasible": None
#                 }
#             )
#         )

#         for smi, stoi in stoich.items():
#             if stoi < 0:
#                 edge_list.append((get_canon_smiles(smi), i))
#             else:
#                 edge_list.append((i, get_canon_smiles(smi)))
    
#     # Create Graph
#     DG = nx.DiGraph()  
#     DG.add_nodes_from(cpd_node_list, bipartite=0)
#     DG.add_nodes_from(rxn_node_list, bipartite=1)
#     DG.add_edges_from(edge_list)
    
#     return DG, rxn_node_list, edge_list

def get_stoich_pk(rxn_id, pk):
    reaction = pk.reactions[rxn_id]
    lhs, rhs = [r for r in reaction["SMILES_rxn"].split("=>")]

    lhs, rhs = [[rct.strip(" ") for rct in side.split(" + ")] for side in [lhs, rhs]]
    [rct.split(" ") for rct in lhs]

    pat = re.compile("\((\d+)\) (.+)")

    reactants = {get_canon_smiles(smiles): -1*int(stoich) for stoich, smiles in [pat.findall(rct)[0] for rct in lhs]}
    products = {get_canon_smiles(smiles): int(stoich) for stoich, smiles in [pat.findall(rct)[0] for rct in rhs]}
    reactants.update(products)

    return reactants

def count_cpds_and_rxns(G):
    compound_nodes = set()
    for n, d in G.nodes(data=True):
        try:
            if d["bipartite"] == 0:
                compound_nodes.add(n)
        except:
            continue
    reaction_nodes = G.nodes() - compound_nodes

    return compound_nodes, reaction_nodes