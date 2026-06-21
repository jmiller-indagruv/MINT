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

Usage:
    python analyze_contacts.py 9iu2.cif
    python analyze_contacts.py 9iu2.cif -c 4.0 --hb-cutoff 3.5 -o contacts.txt
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
# Chain classification
# ---------------------------------------------------------------------------

def classify_chains(data):
    protein_chains = set()
    dna_chains = set()

    types   = data.get("_entity_poly.type", [])
    strands = data.get("_entity_poly.pdbx_strand_id", [])

    for etype, strand_field in zip(types, strands):
        chain_ids = [c.strip() for c in strand_field.split(",") if c.strip()]
        etype_lc = etype.lower()
        if "polypeptide" in etype_lc:
            protein_chains.update(chain_ids)
        elif "deoxyribonucleotide" in etype_lc or "ribonucleotide" in etype_lc:
            dna_chains.update(chain_ids)

    return protein_chains, dna_chains


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

def find_contacts(protein_atoms, dna_atoms, cutoff, hb_cutoff):
    """
    Return dict:
        (prot_chain, prot_seq, prot_comp, dna_chain, dna_seq, dna_comp)
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
                key = (pa["chain"], pa["seq_id"], pa["comp_id"],
                       da["chain"], da["seq_id"], da["comp_id"])
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
    try:
        return (k[0], int(k[1]), k[2], k[3], int(k[4]))
    except ValueError:
        return (k[0], 0, k[2], k[3], 0)


def report_summary(contacts, cutoff, hb_cutoff, output=None):
    fh = open(output, "w", encoding="utf-8") if output else sys.stdout

    # Pre-compute per-pair metadata
    rows = []
    for key in sorted(contacts, key=_sort_key):
        pc, ps, pcomp, dc, ds, dcomp = key
        pairs = contacts[key]
        min_dist = min(d for _, _, d, _ in pairs)
        pair_types = [pt for _, _, _, pt in pairs]
        n_hb   = sum(1 for t in pair_types if t == "hbond")
        n_hp   = sum(1 for t in pair_types if t == "hydrophobic")
        rtype  = residue_contact_type(pair_types)
        rows.append((pc, ps, pcomp, dc, ds, dcomp, len(pairs), min_dist,
                     n_hb, n_hp, rtype))

    # Counts by type
    type_counts = defaultdict(int)
    for r in rows:
        type_counts[r[10]] += 1

    print(f"Protein-DNA contacts within {cutoff:.1f} A  "
          f"(H-bond cutoff: {hb_cutoff:.1f} A)", file=fh)
    print("=" * 96, file=fh)
    print(f"{'Prot':<6}{'Res#':<7}{'AA':<6}"
          f"{'DNA':<6}{'Base#':<7}{'Base':<5}"
          f"{'Pairs':<7}{'HB':<5}{'HP':<5}{'MinDist(A)':<12}{'Contact type'}",
          file=fh)
    print("-" * 96, file=fh)

    for (pc, ps, pcomp, dc, ds, dcomp,
         npairs, min_dist, n_hb, n_hp, rtype) in rows:
        print(f"{pc:<6}{ps:<7}{pcomp:<6}{dc:<6}{ds:<7}{dcomp:<5}"
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
        pc, ps, pcomp, dc, ds, dcomp = key
        pairs = sorted(contacts[key], key=lambda t: t[2])
        pair_types = [pt for _, _, _, pt in pairs]
        rtype = residue_contact_type(pair_types)

        print(f"\n{pcomp}{ps} (chain {pc})  <-->  {dcomp}{ds} (chain {dc})"
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
    args = parser.parse_args()

    if args.hb_cutoff > args.cutoff:
        print(f"WARNING: H-bond cutoff ({args.hb_cutoff} A) exceeds general "
              f"cutoff ({args.cutoff} A) — no H-bonds will be found.")

    print(f"Parsing {args.cif} ...", flush=True)
    data = parse_cif(args.cif)

    protein_chains, dna_chains = classify_chains(data)
    if not protein_chains and not dna_chains:
        print("WARNING: No protein or DNA chains found.")
        sys.exit(1)

    print(f"Protein chains : {sorted(protein_chains)}")
    print(f"DNA chains     : {sorted(dna_chains)}")

    protein_atoms, dna_atoms = load_atoms(data, protein_chains, dna_chains)
    print(f"Protein atoms  : {len(protein_atoms):,}")
    print(f"DNA atoms      : {len(dna_atoms):,}")

    if not protein_atoms or not dna_atoms:
        print("ERROR: No atoms loaded for one or both molecule types.")
        sys.exit(1)

    print(f"Searching for contacts within {args.cutoff} A "
          f"(H-bond cutoff: {args.hb_cutoff} A) ...", flush=True)
    contacts = find_contacts(protein_atoms, dna_atoms, args.cutoff, args.hb_cutoff)

    report_summary(contacts, args.cutoff, args.hb_cutoff, args.output)

    if args.detailed or args.detailed_output:
        report_detailed(contacts, args.cutoff, args.hb_cutoff, args.detailed_output)


if __name__ == "__main__":
    main()
