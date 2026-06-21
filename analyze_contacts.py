"""
Analyze protein-DNA contacts in a mmCIF file.

Reports pairs of (amino acid residue, DNA base) whose atoms come within
a user-specified distance cutoff (default 4.0 A). Chains are classified
by entity type parsed from the _entity_poly block.

Usage:
    python analyze_contacts.py 9iu2.cif
    python analyze_contacts.py 9iu2.cif -c 3.5 -o contacts.txt
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
    """
    Yield (token_type, value) tuples where token_type is one of:
      'data'      - data_XXXX block header
      'loop'      - loop_ keyword
      'key'       - _category.item name
      'value'     - a scalar value (may be quoted or unquoted)
      'comment'   - # comment line (usually skipped by callers)
    """
    with open(path, encoding="utf-8", errors="replace") as fh:
        lines = fh.readlines()

    i = 0
    while i < len(lines):
        line = lines[i]
        stripped = line.strip()

        # Semi-colon multi-line string: starts at column 0
        if line.startswith(";"):
            chunks = []
            i += 1
            while i < len(lines) and not lines[i].startswith(";"):
                chunks.append(lines[i].rstrip("\n"))
                i += 1
            i += 1  # closing ";"
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

        # Tokenise remaining words on the line
        for tok in _line_tokens(stripped):
            lc = tok.lower()
            if lc.startswith("_"):
                yield ("key", lc)
            else:
                yield ("value", tok)

        i += 1


def _line_tokens(line):
    """Split a single CIF line into string tokens (handles '' and "" quoting)."""
    tokens = []
    i = 0
    while i < len(line):
        ch = line[i]
        if ch in (" ", "\t"):
            i += 1
        elif ch == "#":
            break  # inline comment
        elif ch in ("'", '"'):
            # quoted string
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
# mmCIF parser — builds category -> column -> [values] mapping
# ---------------------------------------------------------------------------

def parse_cif(path):
    """
    Return dict: lowercase_category.item -> list-of-str  (loop columns)
                 lowercase_category.item -> str           (scalar key-value)
    """
    data = {}
    tokens = list(tokenize_cif(path))
    i = 0

    while i < len(tokens):
        ttype, tval = tokens[i]

        if ttype == "loop":
            i += 1
            # Collect column names
            cols = []
            while i < len(tokens) and tokens[i][0] == "key":
                cols.append(tokens[i][1])
                i += 1
            # Initialise lists
            for c in cols:
                data[c] = []
            # Collect value rows
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
                i += 1  # malformed — skip

        else:
            i += 1

    return data


# ---------------------------------------------------------------------------
# Chain classification
# ---------------------------------------------------------------------------

def classify_chains(data):
    """Return (protein_chains, dna_chains) as sets of label_asym_id strings."""
    protein_chains = set()
    dna_chains = set()

    types = data.get("_entity_poly.type", [])
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
    """Return (protein_atoms, dna_atoms) — lists of atom dicts."""
    chain_col  = data.get("_atom_site.label_asym_id", [])
    seq_col    = data.get("_atom_site.label_seq_id", [])
    comp_col   = data.get("_atom_site.label_comp_id", [])
    atom_col   = data.get("_atom_site.label_atom_id", [])
    x_col      = data.get("_atom_site.cartn_x", [])
    y_col      = data.get("_atom_site.cartn_y", [])
    z_col      = data.get("_atom_site.cartn_z", [])
    group_col  = data.get("_atom_site.group_pdb", [])
    model_col  = data.get("_atom_site.pdbx_pdb_model_num", [])

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

        atom = {
            "chain":   chain,
            "seq_id":  seq_col[i] if seq_col else "?",
            "comp_id": comp_col[i] if comp_col else "?",
            "atom_id": atom_col[i] if atom_col else "?",
            "x": x, "y": y, "z": z,
        }

        if chain in protein_chains and group == "ATOM":
            protein_atoms.append(atom)
        elif chain in dna_chains and group in ("ATOM", "HETATM"):
            dna_atoms.append(atom)

    return protein_atoms, dna_atoms


# ---------------------------------------------------------------------------
# Contact search
# ---------------------------------------------------------------------------

def find_contacts(protein_atoms, dna_atoms, cutoff):
    """
    Return dict: (prot_chain, prot_seq, prot_comp, dna_chain, dna_seq, dna_comp)
                 -> [(prot_atom_id, dna_atom_id, distance), ...]
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
                key = (pa["chain"], pa["seq_id"], pa["comp_id"],
                       da["chain"], da["seq_id"], da["comp_id"])
                contacts[key].append((pa["atom_id"], da["atom_id"], math.sqrt(d2)))

    return contacts


# ---------------------------------------------------------------------------
# Reporting
# ---------------------------------------------------------------------------

def _sort_key(k):
    try:
        return (k[0], int(k[1]), k[2], k[3], int(k[4]))
    except ValueError:
        return (k[0], 0, k[2], k[3], 0)


def report_summary(contacts, cutoff, output=None):
    fh = open(output, "w", encoding="utf-8") if output else sys.stdout

    print(f"Protein-DNA contacts within {cutoff:.1f} A", file=fh)
    print("=" * 72, file=fh)
    print(f"{'Prot chain':<12}{'Res#':<8}{'AA':<6}"
          f"{'DNA chain':<12}{'Base#':<8}{'Base':<6}"
          f"{'Atom pairs':<12}{'Min dist (A)'}", file=fh)
    print("-" * 72, file=fh)

    for key in sorted(contacts, key=_sort_key):
        pc, ps, pcomp, dc, ds, dcomp = key
        pairs = contacts[key]
        min_dist = min(d for _, _, d in pairs)
        print(f"{pc:<12}{ps:<8}{pcomp:<6}{dc:<12}{ds:<8}{dcomp:<6}"
              f"{len(pairs):<12}{min_dist:.2f}", file=fh)

    print(f"\nTotal unique residue-base contact pairs: {len(contacts)}", file=fh)

    if output:
        fh.close()
        print(f"Summary written to {output}")


def report_detailed(contacts, cutoff, output=None):
    fh = open(output, "w", encoding="utf-8") if output else sys.stdout

    print(f"Detailed protein-DNA contacts within {cutoff:.1f} A", file=fh)
    print("=" * 60, file=fh)

    for key in sorted(contacts, key=_sort_key):
        pc, ps, pcomp, dc, ds, dcomp = key
        pairs = sorted(contacts[key], key=lambda t: t[2])
        print(f"\n{pcomp}{ps} (chain {pc})  <-->  {dcomp}{ds} (chain {dc})", file=fh)
        for pa, da, dist in pairs:
            print(f"    {pa:<6} -- {da:<6}  {dist:.2f} A", file=fh)

    if output:
        fh.close()
        print(f"Detailed output written to {output}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Analyze protein-DNA residue contacts in a mmCIF file.")
    parser.add_argument("cif", help="Path to mmCIF (.cif) file")
    parser.add_argument("-c", "--cutoff", type=float, default=4.0,
                        help="Distance cutoff in Angstroms (default: 4.0)")
    parser.add_argument("-o", "--output", default=None,
                        help="Write summary table to this file")
    parser.add_argument("-d", "--detailed", action="store_true",
                        help="Also print every atom pair (to stdout)")
    parser.add_argument("--detailed-output", default=None,
                        help="Write detailed atom-pair list to this file")
    args = parser.parse_args()

    print(f"Parsing {args.cif} ...", flush=True)
    data = parse_cif(args.cif)

    protein_chains, dna_chains = classify_chains(data)
    if not protein_chains and not dna_chains:
        print("WARNING: No protein or DNA chains found. "
              "Check that _entity_poly entries are present in the file.")
        sys.exit(1)

    print(f"Protein chains : {sorted(protein_chains)}")
    print(f"DNA chains     : {sorted(dna_chains)}")

    protein_atoms, dna_atoms = load_atoms(data, protein_chains, dna_chains)
    print(f"Protein atoms  : {len(protein_atoms):,}")
    print(f"DNA atoms      : {len(dna_atoms):,}")

    if not protein_atoms or not dna_atoms:
        print("ERROR: No atoms loaded for one or both molecule types. "
              "Cannot compute contacts.")
        sys.exit(1)

    print(f"Searching for contacts within {args.cutoff} A ...", flush=True)
    contacts = find_contacts(protein_atoms, dna_atoms, args.cutoff)

    report_summary(contacts, args.cutoff, args.output)

    if args.detailed or args.detailed_output:
        report_detailed(contacts, args.cutoff, args.detailed_output)


if __name__ == "__main__":
    main()
