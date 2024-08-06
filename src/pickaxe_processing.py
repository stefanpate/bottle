import re
import networkx as nx
from rdkit.Chem import CanonSmiles
from minedatabase.pickaxe import Pickaxe
from typing import Iterable

def prune_pickaxe(pk:Pickaxe, paths:Iterable):
    '''
    Prune Pickaxe object of all but reactions and compounds in
    paths
    '''
    pruned_rxns = set()
    pruned_cpds = set()

    for p in paths:
        pruned_cpds.add(p[0]) # Starting molecule

        for i in range(1, len(p)):
            pruned_rxns.add(p[i])
            
            pk_rxn = pk.reactions[p[i]]
            compounds = pk_rxn["Reactants"] + pk_rxn["Products"]

            for _, cpd in compounds:
                pruned_cpds.add(cpd)
    
    pk.reactions = {k:pk.reactions[k] for k in pruned_rxns}
    pk.compounds = {k:pk.compounds[k] for k in pruned_cpds}

    return pk

def get_canon_smiles(smi):
    try:
        return CanonSmiles(smi)
    except BaseException as e:
        return smi

def find_paths(G, starters, target, n_gen):
    """
    Get the pathways from starter to target compounds.

    Parameters
    ----------
    G : nx.DiGraph
        The DiGraph containing network information.
    starters : List[str]
        A list of pickaxe cids for the desired end nodes
    target : str
        The pickaxe cid of the beginning node
    n_gen : int
        Number of generations.
    """
    max_depth = n_gen * 2 # Bipartite mol -> rxn -> mol graph

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
    depth_first_reversed(G, target, starters, [], 0)

    # Reverse paths, extract starting compound id and reaction ids
    if found_paths:
        found_paths = list(
            set(
                [tuple([path[0]] + path[1::2]) for path in [[*reversed(ind_path)] for ind_path in found_paths]]
            )
        )

    return found_paths

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

def pk_rhash_to_smarts(rhash, pk):
    '''
    Make reaction smarts string for
    reaction indexed by rhash in a pk
    object
    '''
    rxn_stoich = get_stoich_pk(rhash, pk)
    products = ".".join([".".join([smi]*stoich) for smi, stoich in rxn_stoich.items() if stoich >= 0])
    reactants = ".".join([".".join([smi]*abs(stoich)) for smi, stoich in rxn_stoich.items() if stoich <= 0])
    rxn_sma = ">>".join([reactants, products])
    return rxn_sma