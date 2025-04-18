import itertools

from .Smile import check_smile, remove_canonical_duplicates


def make_base(backbone_info):
    """
    Creates a list of atomic species based on the provided backbone information.

    Parameters:
    backbone_info (dict): A dictionary where keys are atomic species (str) and values are their counts (int).

    Returns:
    list: A list containing the atomic species repeated according to their counts.
    """
    base = []
    for specie, count in backbone_info.items():
        for i in range(count):
            base.append(specie)
    return base


def get_all_bases(backbone_info):
    """
    Generates all unique permutations of atomic species based on the provided backbone information.

    Parameters:
    backbone_info (dict): A dictionary where keys are atomic species (str) and values are their counts (int).

    Returns:
    list: A list of tuples, each containing a unique permutation of the atomic species.
    """
    base = make_base(backbone_info)
    return list(set(list(itertools.permutations(base))))


def get_all_backbones(backbone_info):
    """
    Generates all unique permutations of atomic species based on the provided backbone information.

    Parameters:
    backbone_info (dict): A dictionary where keys are atomic species (str) and values are their counts (int).

    Returns:
    list: A list of tuples, each containing a unique permutation of the atomic species.
    """
    base = make_base(backbone_info)
    return list(set(list(itertools.permutations(base))))


#     base = make_base(backbone_info)
#     return get_all_bases(base)


def get_all_huged_backbones(backbone, a, b):
    """
    Generates all possible backbones by inserting two specified atomic species into all possible positions.

    Parameters:
    backbone (list): A list of atomic species representing the backbone.
    a (str): The first atomic species to be inserted.
    b (str): The second atomic species to be inserted.

    Returns:
    list: A list of lists, each containing a unique backbone with the two atomic species inserted.
    """
    rng = [n for n in range(0, len(backbone) + 1)]
    out_trj = []
    for i in rng:
        for j in rng[i + 1 :]:
            bckb = list(backbone).copy()
            bckb.insert(i, a)
            bckb.insert(j + 1, b)
            out_trj.append(bckb)
    return out_trj


def get_cl_marked(backbones):
    """
    Adds a 'Cl' atom at the beginning of each backbone in the provided list of backbones.

    Parameters:
    backbones (list): A list of backbones, where each backbone is a list of atomic species.

    Returns:
    list: A list of backbones with 'Cl' added at the beginning of each.
    """
    out_trj = []
    for backbone in backbones:
        x = ["Cl"] + list(backbone)
        out_trj.append(x)
    return out_trj


def get_s1s_marked(backbones):
    """
    Adds 'S1S' and '1' atoms at various positions in each backbone in the provided list of backbones.

    Parameters:
    backbones (list): A list of backbones, where each backbone is a list of atomic species.

    Returns:
    list: A list of backbones with 'S1S' and '1' added at various positions.
    """
    out_trj = []
    for backbone in backbones:
        backbone = list(backbone)
        rng = [n for n in range(0, len(backbone) + 1)]

        for i in [0]:
            for j in rng[i + 1 :]:
                bckb = backbone.copy()
                bckb.insert(i, "S1S")
                bckb.insert(j + 1, "1")
                out_trj.append(bckb)
    return out_trj


def insert_unsaturated_bond(backbones, bond):
    """
    Inserts an unsaturated bond at various positions in each backbone in the provided list of backbones.

    Parameters:
    backbones (list): A list of backbones, where each backbone is a list of atomic species.
    bond (str): The unsaturated bond to be inserted.

    Returns:
    list: A list of backbones with the unsaturated bond inserted at various positions.
    """
    out_trj = []
    for backbone in backbones:
        backbone = list(backbone)
        rng = [n for n in range(0, len(backbone) + 1)]

        for i in [0]:
            for j in rng[i + 1 :]:
                bckb = backbone.copy()
                # bckb.insert(i, 'S1S')
                bckb.insert(j, bond)
                out_trj.append(bckb)
    return out_trj


def get_all_side_chains(backbones):
    """
    Generates all possible side chains by inserting '(' and ')' at various positions in each backbone in the provided list of backbones.

    Parameters:
    backbones (list): A list of backbones, where each backbone is a list of atomic species.

    Returns:
    list: A list of backbones with '(' and ')' inserted at various positions to represent side chains.
    """
    out_trj = []
    for backbone in backbones:
        out_trj += get_all_huged_backbones(backbone, "(", ")")
    # rng = [n for n in range(0, len(backbone)+1)]
    # out_trj = []
    # for i in rng:
    #     for j in rng[i+1:]:
    #         bckb = backbone.copy()
    #         bckb.insert(i, '(')
    #         bckb.insert(j+1, ')')
    #         out_trj.append(bckb)
    return out_trj


def get_rings(backbone, ring_marker):
    """
    Generates all possible ring structures by inserting a ring marker at various positions in the provided backbone.

    Parameters:
    backbone (list): A list of atomic species representing the backbone.
    ring_marker (str): The marker to be used for indicating ring positions.

    Returns:
    list: A list of backbones with the ring marker inserted at various positions to represent rings.
    """
    out_trj = get_all_huged_backbones(backbone, str(ring_marker), str(ring_marker))
    return out_trj


def get_all_ringed(backbones, ring_marker):
    """
    Generates all possible ring structures for each backbone in the provided list of backbones by inserting a ring marker at various positions.

    Parameters:
    backbones (list): A list of backbones, where each backbone is a list of atomic species.
    ring_marker (str): The marker to be used for indicating ring positions.

    Returns:
    list: A list of backbones with the ring marker inserted at various positions to represent rings.
    """
    out_trj = []
    for backbone in backbones:
        backbone = list(backbone)
        out_trj += get_rings(backbone, ring_marker)
    return out_trj


def make_unsaturated(backbone, brackets):
    """
    Generates all possible unsaturated structures by inserting brackets around each atomic species in the provided backbone.

    Parameters:
    backbone (list): A list of atomic species representing the backbone.
    brackets (list): A list of bracket pairs to be used for indicating unsaturation.

    Returns:
    list: A list of backbones with brackets inserted around each atomic species to represent unsaturation.
    """
    unsaturated = []
    open_bracket = brackets[0]
    for i, _ in enumerate(backbone):
        for closed_bracket in brackets[1:]:
            bckb = backbone.copy()
            if open_bracket not in bckb[i]:
                bckb[i] = f"{open_bracket}{bckb[i]}{closed_bracket}"
                unsaturated.append(bckb)

    return unsaturated


def make_all_unsaturated_backbones(all_backbones, brackets):
    """
    Generates all possible unsaturated structures for each backbone in the provided list of backbones by inserting brackets around each atomic species.

    Parameters:
    all_backbones (list): A list of backbones, where each backbone is a list of atomic species.
    brackets (list): A list of bracket pairs to be used for indicating unsaturation.

    Returns:
    list: A list of backbones with brackets inserted around each atomic species to represent unsaturation.
    """
    all_unsaturated_backbones = []
    for backbone in all_backbones:
        for _ in range(len(backbone)):
            all_unsaturated_backbones += make_unsaturated(backbone, brackets)


def get_all_unsaturated_from_backbone(backbone, brackets, return_dict=False):
    """
    Generates all possible unsaturated structures from a backbone by inserting brackets around each atomic species.

    Parameters:
    backbone (list): A list of atomic species representing the backbone.
    brackets (list): A list of bracket pairs to be used for indicating unsaturation.
    return_dict (bool): Whether to return the results as a dictionary. Default is False.

    Returns:
    list or dict: A list of backbones with brackets inserted around each atomic species to represent unsaturation.
                  If return_dict is True, returns a dictionary where keys are the number of unsaturated positions and values are the corresponding backbones.
    """
    multiple_unsaturated = {}

    for i in range(len(backbone)):
        if i == 0:
            multiple_unsaturated[i + 1] = make_unsaturated(backbone, brackets)
        else:
            unsaturated = []
            for bckb in multiple_unsaturated[i]:
                unsaturated += make_unsaturated(bckb, brackets)
                multiple_unsaturated[i + 1] = unsaturated
    if return_dict:
        return multiple_unsaturated
    else:
        lst = []
        for k, v in multiple_unsaturated.items():
            lst += v
        return lst


def get_checked_smiles(smiles_list):
    """
    Filters a list of SMILES strings, returning only those that pass a validity check.

    Parameters:
    smiles_list (list): A list of SMILES strings to be checked.

    Returns:
    list: A list of valid SMILES strings that passed the check.
    """
    checked_smiles = []
    for s in smiles_list:
        if check_smile(s):
            checked_smiles.append(s)
    return checked_smiles


def get_smiles_from_backbones(backbones):
    """
    Converts a list of backbones into their corresponding SMILES strings.

    Parameters:
    backbones (list): A list of backbones, where each backbone is a list of atomic species.

    Returns:
    list: A list of SMILES strings generated from the backbones.
    """
    smiles = []
    for backbone in backbones:
        smiles.append("".join(backbone))
    return smiles


def xx_get_special_symbols(config):
    """
    Generates a dictionary of special symbols based on the provided configuration.

    Parameters:
    config (dict): A dictionary containing configuration information, including 'backbone_info' and 'specials'.

    Returns:
    dict: A dictionary where keys are atomic symbols and values are lists of special symbols.
    """
    special_symbols = {}
    for symbol, _ in config["backbone_info"].items():
        special_symbols[symbol] = [symbol] + [
            f"[{symbol}{marker}" for marker in config["specials"]
        ]
    return special_symbols


def xx_unpack_symbols(config):
    """
    Unpacks special symbols based on the provided configuration.

    Parameters:
    config (dict): A dictionary containing configuration information, including 'backbone_info' and 'specials'.

    Returns:
    list: A list of lists, where each sublist contains unpacked symbols for each atomic species.
    """
    special_symbols = xx_get_special_symbols(config)
    print("special_symbols: ", special_symbols)

    unpacked_symbols = []
    for symbol, value in config["backbone_info"].items():
        print(symbol, value)
        if value == 0:
            continue
        if value == 1:
            # unpacked_symbols.append([tuple([s]) for s in special_symbols[symbol]])
            unpacked_symbols.append(tuple(special_symbols[symbol]))
        if value > 1:
            lst = list(itertools.product(special_symbols[symbol], repeat=value))
            print("lst: ", lst)
            # joined_lst_elements = [''.join(t) for t in lst]
            joined_lst_elements = [list(t) for t in lst]
            print("joined_lst_elements: ", joined_lst_elements)
            unpacked_symbols.append(joined_lst_elements)
            # for group in list(itertools.product(special_symbols[symbol], repeat=value)):
            #     unpacked_symbols.append(list(group))
            #     print(unpacked_symbols)
    return unpacked_symbols


def xx_get_all_backbones(config):
    """
    Generates all possible backbones based on the provided configuration.

    Parameters:
    config (dict): A dictionary containing configuration information, including 'backbone_info' and 'specials'.

    Returns:
    list: A list of tuples, each containing a unique combination of unpacked symbols for the backbones.
    """
    unpacked_symbols = xx_unpack_symbols(config)
    print(unpacked_symbols)
    print("xxxx")
    all_backbones = list(itertools.product(*unpacked_symbols))
    return all_backbones


def concat_tuples(list_of_tuples):
    """
    Concatenates a list of tuples into a single tuple.

    Parameters:
    list_of_tuples (list): A list of tuples to be concatenated.

    Returns:
    tuple: A single tuple containing all elements from the input list of tuples.
    """
    t = ()
    for x in list_of_tuples:
        t += x
    return t


def construct_smiles(config):
    """
    Constructs a list of unique SMILES strings based on the provided configuration.

    Parameters:
    config (dict): A dictionary containing configuration information, including:
        - 'backbone_info' (dict): Information about the backbone atomic species and their counts.
        - 'brackets' (list): A list of bracket pairs to be used for indicating unsaturation.
        - 'allow_intramolec_rings' (bool): Whether to allow intramolecular rings.
        - 'ring_marker' (str): The marker to be used for indicating ring positions.
        - 'make_labeled' (bool): Whether to label the backbones with special markers.
        - 'specials' (list): A list of special markers to be used for labeling.

    Returns:
    list: A list of unique SMILES strings generated based on the configuration.
    """
    basic_backbones = get_all_backbones(config["backbone_info"])

    multiple_unsaturated = []
    for backbone in basic_backbones:
        backbone = list(backbone)
        multiple_unsaturated += get_all_unsaturated_from_backbone(
            backbone, config["brackets"]
        )

    all_backbones = basic_backbones + multiple_unsaturated

    all_backbones += get_all_side_chains(all_backbones)

    if config["allow_intramolec_rings"]:
        ring_backbones = get_all_ringed(all_backbones, config["ring_marker"])
        all_backbones += ring_backbones

    if config["make_labeled"]:
        cl_backbones = get_cl_marked(all_backbones)
        s1s_backbones = get_s1s_marked(all_backbones)
        all_backbones = cl_backbones + s1s_backbones

    double = insert_unsaturated_bond(all_backbones, "=")
    triple = insert_unsaturated_bond(all_backbones, "#")
    all_backbones += double
    all_backbones += triple

    print(all_backbones[-10:])
    smiles_list = get_smiles_from_backbones(all_backbones)

    checked_smiles = get_checked_smiles(smiles_list)
    smiles = remove_canonical_duplicates(checked_smiles)

    return smiles


_example_config = {
    "backbone_info": {"C": 1, "N": 0, "O": 2},
    "allow_intramolec_rings": True,
    "ring_marker": 2,
    "side_chain": ["(", ")"],
    "brackets": ["[", "]", "H2]", "H3]", "H-]", "H+]"],  # , '-]', '--]', '---]']
    "make_labeled": True,
}


def _show_ussage():
    print("eg ussage: smiles = construct_smiles(config)")
