"""
Analyze protein-DNA contacts in a mmCIF file.

For each residue-base pair within the distance cutoff, classifies contacts as:

  H-bond      -- N/O...N/O heavy-atom pair within the H-bond cutoff (default
                 3.5 A). Covers backbone amide/carbonyl, polar sidechains,
                 and DNA base/phosphate/sugar acceptors and donors.

  Hydrophobic -- C...C or nonpolar-S...C pair within the general cutoff.
                 Covers base-stacking, aliphatic packing, and aromatic
                 interactions between protein and DNA.

  Mixed       -- residue-base pair that has both H-bond and hydrophobic
                 contacts, or only cross-type pairs (e.g. C...O).

Classification is heuristic — no hydrogen positions are required — and is
based solely on heavy-atom element identity and distance, which is appropriate
for cryo-EM structures where H-atoms are not modelled.

DNA positions are reported in center-out notation relative to the recombination
site center:
  L1, L2, L3 ... = bases going 5' from the crossover point
  R1, R2, R3 ... = bases going 3' from the crossover point

For attB, the crossover is the junction between the attB-L and attB-R chains.
For attP (continuous duplex), the center defaults to the midpoint of the top
strand and can be overridden with --attP-center.

Usage:
    python analyze_contacts.py 9iu2.cif
    python analyze_contacts.py 9iu2.cif -c 4.0 --hb-cutoff 3.5 -o contacts.txt
    python analyze_contacts.py 9iu2.cif --attP-center 26 --attP-top-chain E
    python analyze_contacts.py 9iu2.cif --detailed --detailed-output atom_pairs.txt
"""

import sys
import math
import argparse
from collections import defaultdict


# ---------------------------------------------------------------------------
# mmCIF tokenizer — correctly handles semicolon multi-line strings
# ---------------------------------------------------------------------------

def tokenize_cif(path):
    with open(path, encoding="utf-8", errors="replace") as fh:
        lines = fh.readlines()

    i = 0
    while i < len(lines):
        line = lines[i]
        stripped = line.strip()

        if line.startswith(";"):
            chunks = []
            i += 1
            while i < len(lines) and not lines[i].startswith(";"):
                chunks.append(lines[i].rstrip("\n"))
                i += 1
            i += 1
            yield ("value", "\n".join(chunks))
            continue

        if not stripped or stripped.startswith("#"):
            i += 1
            continue

        if stripped.lower().startswith("data_"):
            yield ("data", stripped)
            i += 1
            continue

        if stripped.lower() == "loop_":
            yield ("loop", "loop_")
            i += 1
            continue

        for tok in _line_tokens(stripped):
            lc = tok.lower()
            if lc.startswith("_"):
                yield ("key", lc)
            else:
                yield ("value", tok)

        i += 1


def _line_tokens(line):
    tokens = []
    i = 0
    while i < len(line):
        ch = line[i]
        if ch in (" ", "\t"):
            i += 1
        elif ch == "#":
            break
        elif ch in ("'", '"'):
            j = i + 1
            while j < len(line):
                if line[j] == ch and (j + 1 >= len(line) or line[j + 1] in (" ", "\t", "")):
                    break
                j += 1
            tokens.append(line[i + 1:j])
            i = j + 1
        else:
            j = i
            while j < len(line) and line[j] not in (" ", "\t"):
                j += 1
            tokens.append(line[i:j])
            i = j
    return tokens


# ---------------------------------------------------------------------------
# mmCIF parser
# ---------------------------------------------------------------------------

def parse_cif(path):
    data = {}
    tokens = list(tokenize_cif(path))
    i = 0

    while i < len(tokens):
        ttype, tval = tokens[i]

        if ttype == "loop":
            i += 1
            cols = []
            while i < len(tokens) and tokens[i][0] == "key":
                cols.append(tokens[i][1])
                i += 1
            for c in cols:
                data[c] = []
            while i < len(tokens) and tokens[i][0] == "value":
                for c in cols:
                    if i < len(tokens) and tokens[i][0] == "value":
                        data[c].append(tokens[i][1])
                        i += 1
                    else:
                        data[c].append("?")

        elif ttype == "key":
            key = tval
            i += 1
            if i < len(tokens) and tokens[i][0] == "value":
                data[key] = tokens[i][1]
                i += 1
            else:
                i += 1

        else:
            i += 1

    return data


# ---------------------------------------------------------------------------
# Chain ID mapping (label_asym_id <-> auth_asym_id)
# ---------------------------------------------------------------------------

def build_chain_id_maps(data):
    """
    Return (label_to_auth, auth_to_label) dicts built from _atom_site columns.

    mmCIF uses two parallel chain-ID systems:
      label_asym_id -- internal sequential IDs used in atom coordinates (A, B, ... G, H, ...)
      auth_asym_id  -- author-assigned IDs used in _entity_poly.pdbx_strand_id (A, G1, H2, ...)

    pdbx_strand_id contains author IDs, so we need this map to find the right
    atoms for each entity.
    """
    label_col = data.get("_atom_site.label_asym_id", [])
    auth_col  = data.get("_atom_site.auth_asym_id",  [])
    label_to_auth = {}
    auth_to_label = {}
    for lbl, auth in zip(label_col, auth_col):
        label_to_auth.setdefault(lbl, auth)
        auth_to_label.setdefault(auth, lbl)
    return label_to_auth, auth_to_label


# ---------------------------------------------------------------------------
# Chain classification
# ---------------------------------------------------------------------------

def classify_chains(data, auth_to_label):
    """
    Return (protein_chains, dna_chains) as sets of *label_asym_id* strings.
    pdbx_strand_id contains author chain IDs, which are translated via auth_to_label.
    """
    protein_chains = set()
    dna_chains = set()

    types   = data.get("_entity_poly.type", [])
    strands = data.get("_entity_poly.pdbx_strand_id", [])

    for etype, strand_field in zip(types, strands):
        auth_ids  = [c.strip() for c in strand_field.split(",") if c.strip()]
        label_ids = [auth_to_label.get(a, a) for a in auth_ids]
        etype_lc  = etype.lower()
        if "polypeptide" in etype_lc:
            protein_chains.update(label_ids)
        elif "deoxyribonucleotide" in etype_lc or "ribonucleotide" in etype_lc:
            dna_chains.update(label_ids)

    return protein_chains, dna_chains


def get_chain_roles(data, auth_to_label):
    """
    Return dict: chain_id -> entity description string
    e.g. {'E': 'attP', 'F': 'attP', 'G1': 'attB-L', 'G2': 'attB-R', ...}
    """
    entity_desc_map = {}
    for eid, desc in zip(data.get("_entity.id", []),
                         data.get("_entity.pdbx_description", [])):
        entity_desc_map[eid] = desc.strip("'\"")

    roles = {}
    for eid, strand_field in zip(data.get("_entity_poly.entity_id", []),
                                  data.get("_entity_poly.pdbx_strand_id", [])):
        desc = entity_desc_map.get(eid, "unknown")
        for auth_chain in [c.strip() for c in strand_field.split(",") if c.strip()]:
            label_chain = auth_to_label.get(auth_chain, auth_chain)
            roles[label_chain] = desc
    return roles


def compute_seq_ranges(dna_atoms):
    """Return dict: chain -> (min_seq_id_int, max_seq_id_int)."""
    ranges = {}
    for atom in dna_atoms:
        chain = atom["chain"]
        try:
            sid = int(atom["seq_id"])
        except ValueError:
            continue
        lo, hi = ranges.get(chain, (sid, sid))
        ranges[chain] = (min(lo, sid), max(hi, sid))
    return ranges


def build_position_labels(chain_roles, dna_chains, ranges,
                          attP_center=None, attP_top_chain=None):
    """
    Return dict: (chain, seq_id_str) -> center-out label  e.g. 'L3', 'R12'

    attB-L chains: the 3' end (highest seq_id) is L1; count outward 5'.
    attB-R chains: the 5' end (lowest seq_id)  is R1; count outward 3'.

    attP top strand:  positions 1..center -> L{center}..L1,
                      positions center+1..end -> R1..R{end-center}
    attP bottom strand (antiparallel complement of top):
                      seq_id 1 of bottom pairs with the 3' end of top (R side),
                      so labels mirror the top strand in reverse.
    """
    labels = {}

    attP_chains = sorted(c for c in dna_chains if "attP" in chain_roles.get(c, ""))
    top_chain = (attP_top_chain if attP_top_chain in attP_chains
                 else (attP_chains[0] if attP_chains else None))

    for chain in dna_chains:
        role = chain_roles.get(chain, "")
        if chain not in ranges:
            continue
        lo, hi = ranges[chain]

        if "attB-L" in role:
            # 3' end (hi) = L1, counting outward toward 5'
            for s in range(lo, hi + 1):
                labels[(chain, str(s))] = f"L{hi - s + 1}"

        elif "attB-R" in role:
            # 5' end (lo) = R1, counting outward toward 3'
            for s in range(lo, hi + 1):
                labels[(chain, str(s))] = f"R{s - lo + 1}"

        elif "attP" in role:
            top_lo, top_hi = ranges.get(top_chain, (lo, hi))
            top_len = top_hi - top_lo + 1
            center = attP_center if attP_center is not None else top_len // 2

            if chain == top_chain:
                # 5' side (pos 1..center): L{center}..L1
                # 3' side (pos center+1..end): R1..R{end-center}
                for s in range(lo, hi + 1):
                    pos = s - lo + 1  # 1-based position along chain
                    if pos <= center:
                        labels[(chain, str(s))] = f"L{center - pos + 1}"
                    else:
                        labels[(chain, str(s))] = f"R{pos - center}"
            else:
                # Bottom strand — antiparallel. seq_id lo pairs with top_hi
                # (the 3'/R side), so the first residue of the bottom strand
                # gets the same label as the last residue of the top strand.
                for s in range(lo, hi + 1):
                    paired = top_len - (s - lo)  # paired position on top (1-based)
                    if paired <= center:
                        labels[(chain, str(s))] = f"L{center - paired + 1}"
                    else:
                        labels[(chain, str(s))] = f"R{paired - center}"

    return labels


# ---------------------------------------------------------------------------
# Atom loading
# ---------------------------------------------------------------------------

def load_atoms(data, protein_chains, dna_chains):
    chain_col   = data.get("_atom_site.label_asym_id", [])
    seq_col     = data.get("_atom_site.label_seq_id", [])
    comp_col    = data.get("_atom_site.label_comp_id", [])
    atom_col    = data.get("_atom_site.label_atom_id", [])
    elem_col    = data.get("_atom_site.type_symbol", [])
    x_col       = data.get("_atom_site.cartn_x", [])
    y_col       = data.get("_atom_site.cartn_y", [])
    z_col       = data.get("_atom_site.cartn_z", [])
    group_col   = data.get("_atom_site.group_pdb", [])
    model_col   = data.get("_atom_site.pdbx_pdb_model_num", [])

    protein_atoms = []
    dna_atoms = []

    n = len(chain_col)
    for i in range(n):
        model = model_col[i] if model_col else "1"
        if model != "1":
            continue

        group = group_col[i].upper() if group_col else "ATOM"
        chain = chain_col[i]

        try:
            x = float(x_col[i])
            y = float(y_col[i])
            z = float(z_col[i])
        except (ValueError, IndexError):
            continue

        # Derive element: prefer type_symbol, fall back to first letter of atom name
        if elem_col and i < len(elem_col) and elem_col[i] not in ("?", "."):
            element = elem_col[i].upper()
        else:
            atom_name = (atom_col[i] if atom_col else "X").lstrip("0123456789")
            element = atom_name[0].upper() if atom_name else "X"

        atom = {
            "chain":   chain,
            "seq_id":  seq_col[i]  if seq_col  else "?",
            "comp_id": comp_col[i] if comp_col else "?",
            "atom_id": atom_col[i] if atom_col else "?",
            "element": element,
            "x": x, "y": y, "z": z,
        }

        if chain in protein_chains and group == "ATOM":
            protein_atoms.append(atom)
        elif chain in dna_chains and group in ("ATOM", "HETATM"):
            dna_atoms.append(atom)

    return protein_atoms, dna_atoms


# ---------------------------------------------------------------------------
# Contact classification
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# Display helpers
# ---------------------------------------------------------------------------

_AA1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    "MSE": "M",  # selenomethionine
}


def aa1(three_letter):
    return _AA1.get(three_letter.upper(), three_letter)


def base1(comp_id):
    """DA -> A, DT -> T, DG -> G, DC -> C, A -> A, etc."""
    name = comp_id.upper().lstrip("D")
    return name if name in ("A", "C", "G", "T", "U") else comp_id


# Elements that can participate in H-bonds as donor or acceptor
_HBOND_ELEMENTS = {"N", "O", "S"}

# Nonpolar elements that form hydrophobic / vdW contacts
_HYDROPHOBIC_ELEMENTS = {"C", "S"}

# Met SD and Cys SG are the only S atoms we want to count as hydrophobic;
# other S (none in standard DNA) are rare. We include S in both sets and let
# distance + context disambiguate at the residue level.


def classify_pair(p_elem, d_elem, dist, hb_cutoff):
    """
    Return the contact type for a single atom pair.

    "hbond"      : both atoms are H-bond capable (N/O/S) and within hb_cutoff
    "hydrophobic": both atoms are nonpolar (C or S) within the general cutoff
    "polar-vdw"  : one polar, one nonpolar — electrostatic fringe / weak
    """
    p_polar = p_elem in _HBOND_ELEMENTS
    d_polar = d_elem in _HBOND_ELEMENTS
    p_nonpolar = p_elem in _HYDROPHOBIC_ELEMENTS
    d_nonpolar = d_elem in _HYDROPHOBIC_ELEMENTS

    if p_polar and d_polar and dist <= hb_cutoff:
        return "hbond"
    if p_nonpolar and d_nonpolar:
        return "hydrophobic"
    return "polar-vdw"


def residue_contact_type(pair_types):
    """
    Summarise the set of per-atom-pair types into one residue-level label.

    priority: if any H-bond present → "H-bond (+ vdW)" when also hydrophobic,
              else pure "H-bond"; if only hydrophobic → "Hydrophobic";
              mixed polar-vdw-only → "Mixed (polar-vdW)".
    """
    types = set(pair_types)
    has_hb   = "hbond"       in types
    has_hpho = "hydrophobic" in types
    has_pvdw = "polar-vdw"   in types

    if has_hb and has_hpho:
        return "H-bond + Hydrophobic"
    if has_hb:
        return "H-bond"
    if has_hpho and not has_pvdw:
        return "Hydrophobic"
    if has_hpho:
        return "Hydrophobic (+ polar-vdW)"
    return "Mixed (polar-vdW)"


# ---------------------------------------------------------------------------
# Contact search
# ---------------------------------------------------------------------------

def find_contacts(protein_atoms, dna_atoms, cutoff, hb_cutoff, pos_labels):
    """
    Return dict:
        (prot_chain, prot_seq, prot_comp, dna_chain, dna_seq, dna_comp, dna_pos_label)
        -> [(prot_atom_id, dna_atom_id, distance, pair_type), ...]
    """
    cutoff2 = cutoff * cutoff
    contacts = defaultdict(list)

    for pa in protein_atoms:
        for da in dna_atoms:
            dx = pa["x"] - da["x"]
            dy = pa["y"] - da["y"]
            dz = pa["z"] - da["z"]
            d2 = dx*dx + dy*dy + dz*dz
            if d2 <= cutoff2:
                dist = math.sqrt(d2)
                ptype = classify_pair(pa["element"], da["element"], dist, hb_cutoff)
                dna_label = pos_labels.get((da["chain"], da["seq_id"]), da["seq_id"])
                key = (pa["chain"], pa["seq_id"], pa["comp_id"],
                       da["chain"], da["seq_id"], da["comp_id"], dna_label)
                contacts[key].append((pa["atom_id"], da["atom_id"], dist, ptype))

    return contacts


# ---------------------------------------------------------------------------
# Reporting
# ---------------------------------------------------------------------------

_TYPE_ABBREV = {
    "H-bond":                    "HB",
    "Hydrophobic":               "HP",
    "H-bond + Hydrophobic":      "HB+HP",
    "Hydrophobic (+ polar-vdW)": "HP+pv",
    "Mixed (polar-vdW)":         "Mixed",
}

_TYPE_ORDER = [
    "H-bond",
    "H-bond + Hydrophobic",
    "Hydrophobic",
    "Hydrophobic (+ polar-vdW)",
    "Mixed (polar-vdW)",
]


def _sort_key(k):
    # k = (prot_chain, prot_seq, prot_comp, dna_chain, dna_seq, dna_comp, dna_label)
    # Sort protein by chain then seq_id; DNA by chain then by L/R label numerically.
    dna_label = k[6]  # e.g. 'L3', 'R12', or raw seq_id if no label
    if dna_label and dna_label[0] in ("L", "R"):
        # L labels sort before R labels; within each, sort numerically
        label_order = 0 if dna_label[0] == "L" else 1
        try:
            label_num = int(dna_label[1:])
            # For L labels, higher number = farther from center; sort ascending so L1 comes last
            # Convention: display in order from center outward (L1, L2 ... then R1, R2 ...)
            dna_sort = (k[3], label_order, label_num)
        except ValueError:
            dna_sort = (k[3], 2, 0)
    else:
        try:
            dna_sort = (k[3], 2, int(dna_label))
        except (ValueError, TypeError):
            dna_sort = (k[3], 2, 0)
    try:
        return (k[0], int(k[1]), k[2]) + dna_sort
    except ValueError:
        return (k[0], 0, k[2]) + dna_sort


def report_summary(contacts, cutoff, hb_cutoff, output=None):
    fh = open(output, "w", encoding="utf-8") if output else sys.stdout

    # Pre-compute per-pair metadata
    rows = []
    for key in sorted(contacts, key=_sort_key):
        pc, ps, pcomp, dc, ds, dcomp, dlabel = key
        pairs = contacts[key]
        min_dist = min(d for _, _, d, _ in pairs)
        pair_types = [pt for _, _, _, pt in pairs]
        n_hb   = sum(1 for t in pair_types if t == "hbond")
        n_hp   = sum(1 for t in pair_types if t == "hydrophobic")
        rtype  = residue_contact_type(pair_types)
        rows.append((pc, ps, aa1(pcomp), dc, dlabel, base1(dcomp), len(pairs), min_dist,
                     n_hb, n_hp, rtype))

    # Counts by type
    type_counts = defaultdict(int)
    for r in rows:
        type_counts[r[10]] += 1

    print(f"Protein-DNA contacts within {cutoff:.1f} A  "
          f"(H-bond cutoff: {hb_cutoff:.1f} A)", file=fh)
    print("=" * 96, file=fh)
    print(f"{'Prot':<6}{'Res#':<7}{'AA':<6}"
          f"{'DNA':<6}{'Pos':<7}{'Base':<5}"
          f"{'Pairs':<7}{'HB':<5}{'HP':<5}{'MinDist(A)':<12}{'Contact type'}",
          file=fh)
    print("-" * 96, file=fh)

    for (pc, ps, pcomp, dc, dlabel, dcomp,
         npairs, min_dist, n_hb, n_hp, rtype) in rows:
        print(f"{pc:<6}{ps:<7}{pcomp:<6}{dc:<6}{dlabel:<7}{dcomp:<5}"
              f"{npairs:<7}{n_hb:<5}{n_hp:<5}{min_dist:<12.2f}{rtype}",
              file=fh)

    print(f"\nTotal unique residue-base contact pairs: {len(rows)}", file=fh)

    print("\nBreakdown by contact type:", file=fh)
    for t in _TYPE_ORDER:
        if type_counts[t]:
            print(f"  {t:<30} {type_counts[t]:>4}", file=fh)

    if output:
        fh.close()
        print(f"Summary written to {output}")


def report_detailed(contacts, cutoff, hb_cutoff, output=None):
    fh = open(output, "w", encoding="utf-8") if output else sys.stdout

    print(f"Detailed protein-DNA contacts within {cutoff:.1f} A  "
          f"(H-bond cutoff: {hb_cutoff:.1f} A)", file=fh)
    print("=" * 68, file=fh)

    type_label = {"hbond": "H-bond", "hydrophobic": "Hydrophobic", "polar-vdw": "polar-vdW"}

    for key in sorted(contacts, key=_sort_key):
        pc, ps, pcomp, dc, ds, dcomp, dlabel = key
        pairs = sorted(contacts[key], key=lambda t: t[2])
        pair_types = [pt for _, _, _, pt in pairs]
        rtype = residue_contact_type(pair_types)

        print(f"\n{aa1(pcomp)}{ps} (chain {pc})  <-->  {base1(dcomp)} {dlabel} (chain {dc})"
              f"  [{rtype}]", file=fh)
        for pa, da, dist, pt in pairs:
            label = type_label[pt]
            print(f"    {pa:<6} -- {da:<6}  {dist:.2f} A   {label}", file=fh)

    if output:
        fh.close()
        print(f"Detailed output written to {output}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Analyze and classify protein-DNA contacts in a mmCIF file.")
    parser.add_argument("cif", help="Path to mmCIF (.cif) file")
    parser.add_argument("-c", "--cutoff", type=float, default=4.0,
                        help="General distance cutoff in Angstroms (default: 4.0)")
    parser.add_argument("--hb-cutoff", type=float, default=3.5,
                        help="Heavy-atom distance cutoff for H-bonds (default: 3.5)")
    parser.add_argument("-o", "--output", default=None,
                        help="Write summary table to this file")
    parser.add_argument("-d", "--detailed", action="store_true",
                        help="Print every atom pair with its type (to stdout)")
    parser.add_argument("--detailed-output", default=None,
                        help="Write atom-level detail to this file")
    parser.add_argument("--attP-center", type=int, default=None,
                        help="Last L-arm position (1-based) on the attP top strand. "
                             "Defaults to half the top-strand length. "
                             "For a 52 nt top strand this defaults to 26, placing L1 "
                             "at position 26 and R1 at position 27.")
    parser.add_argument("--attP-top-chain", default=None,
                        help="Chain ID of the attP top strand (default: first "
                             "alphabetically among attP chains, e.g. 'E').")
    args = parser.parse_args()

    if args.hb_cutoff > args.cutoff:
        print(f"WARNING: H-bond cutoff ({args.hb_cutoff} A) exceeds general "
              f"cutoff ({args.cutoff} A) — no H-bonds will be found.")

    print(f"Parsing {args.cif} ...", flush=True)
    data = parse_cif(args.cif)

    label_to_auth, auth_to_label = build_chain_id_maps(data)
    protein_chains, dna_chains = classify_chains(data, auth_to_label)
    if not protein_chains and not dna_chains:
        print("WARNING: No protein or DNA chains found.")
        sys.exit(1)

    chain_roles = get_chain_roles(data, auth_to_label)

    def display(chains):
        return [f"{c}({label_to_auth.get(c,c)})" for c in sorted(chains)]

    print(f"Protein chains : {display(protein_chains)}")
    print(f"DNA chains     : {display(dna_chains)}")
    for chain in sorted(dna_chains):
        auth = label_to_auth.get(chain, chain)
        print(f"  {chain} ({auth}): {chain_roles.get(chain, '?')}")

    protein_atoms, dna_atoms = load_atoms(data, protein_chains, dna_chains)
    print(f"Protein atoms  : {len(protein_atoms):,}")
    print(f"DNA atoms      : {len(dna_atoms):,}")

    if not protein_atoms or not dna_atoms:
        print("ERROR: No atoms loaded for one or both molecule types.")
        sys.exit(1)

    ranges = compute_seq_ranges(dna_atoms)

    # Allow user to specify attP top chain by either author or label ID
    attP_top_arg = args.attP_top_chain
    if attP_top_arg:
        attP_top_arg = auth_to_label.get(attP_top_arg, attP_top_arg)

    pos_labels = build_position_labels(
        chain_roles, dna_chains, ranges,
        attP_center=args.attP_center,
        attP_top_chain=attP_top_arg,
    )

    # Determine which chain was chosen as attP top for display
    attP_top = (attP_top_arg
                or next((c for c in sorted(dna_chains)
                         if "attP" in chain_roles.get(c, "")), None))
    attP_ctr = args.attP_center
    if attP_ctr is None and attP_top and attP_top in ranges:
        lo, hi = ranges[attP_top]
        attP_ctr = (hi - lo + 1) // 2

    print(f"\nDNA position labelling (center-out notation):")
    for chain in sorted(dna_chains):
        auth = label_to_auth.get(chain, chain)
        lo, hi = ranges.get(chain, (1, 1))
        first = pos_labels.get((chain, str(lo)), "?")
        last  = pos_labels.get((chain, str(hi)), "?")
        role  = chain_roles.get(chain, "?")
        if "attP" in role and chain == attP_top:
            role += f" [top strand, center after position {attP_ctr}]"
        elif "attP" in role:
            role += " [bottom strand]"
        print(f"  chain {chain}/{auth} ({role}): seq {lo}={first} .. seq {hi}={last}")

    print(f"\nSearching for contacts within {args.cutoff} A "
          f"(H-bond cutoff: {args.hb_cutoff} A) ...", flush=True)
    contacts = find_contacts(protein_atoms, dna_atoms,
                             args.cutoff, args.hb_cutoff, pos_labels)

    report_summary(contacts, args.cutoff, args.hb_cutoff, args.output)

    if args.detailed or args.detailed_output:
        report_detailed(contacts, args.cutoff, args.hb_cutoff, args.detailed_output)


if __name__ == "__main__":
    main()
