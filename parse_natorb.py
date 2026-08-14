"""
parse_natorb_v2.py
---------------
Parses the NATURAL ORBITALS block from a MOLPRO output file and reports
all molecular orbitals that have contributions from user-specified atom indices.

Two modes:
  1) Atom search: find all MOs with contributions from specified atom indices.
  2) Antibonding partner search: given a reference occupied MO label, find
     the best-matching virtual orbital candidates ranked by strict opposite-sign
     matching of (atom, mu, type) triples.

Usage:
    python parse_natorb_v2.py
    (prompts interactively)
"""

import re
import sys


# ---------------------------------------------------------------------------
# Parsing
# ---------------------------------------------------------------------------

def parse_ao_contributions(text):
    """
    Parse all AO contributions from a text string.
    Returns a list of (cen, mu, typ, coeff) tuples.
    """
    pattern = re.compile(
        r'(?:^|\s)(\d+)\s+(\d+)\s+(\S+)\s+([-+]?\d+\.\d+)'
    )
    return pattern.findall(text)


def parse_natural_orbitals(filepath):
    """
    Reads a MOLPRO .out file, locates the NATURAL ORBITALS block,
    and parses every MO into a structured list.

    Returns a list of dicts:
        {
            'label':  str,   e.g. '88.1'
            'occ':    float,
            'energy': float,
            'aos':    list of (cen, mu, typ, coeff) tuples
        }
    """
    with open(filepath, 'r') as f:
        lines = f.readlines()

    # --- Find the NATURAL ORBITALS block ---
    start_idx = None
    for i, line in enumerate(lines):
        if 'NATURAL ORBITALS' in line:
            start_idx = i
            break

    if start_idx is None:
        print("ERROR: 'NATURAL ORBITALS' section not found in file.")
        print("       Make sure your MOLPRO job included a NATORB directive.")
        sys.exit(1)

    # Find the first orbital line after the header
    orb_header_pattern = re.compile(r'^\s{1,6}\d+\.\d+\s+\d+\.\d+')
    orbital_start = None
    for i in range(start_idx, min(start_idx + 10, len(lines))):
        if orb_header_pattern.match(lines[i]):
            orbital_start = i
            break

    if orbital_start is None:
        print("ERROR: Could not find first orbital line after NATURAL ORBITALS header.")
        sys.exit(1)

    orb_header = re.compile(
        r'^\s{1,6}(\d+\.\d+)\s+([\d.]+)\s+([-\d.]+)\s+(.*)'
    )
    continuation = re.compile(r'^(\s{20,})(.*)')

    orbitals = []
    current_orbital = None

    for line in lines[orbital_start:]:
        m_header = orb_header.match(line)
        if m_header:
            if current_orbital is not None:
                orbitals.append(current_orbital)
            label  = m_header.group(1)
            occ    = float(m_header.group(2))
            energy = float(m_header.group(3))
            rest   = m_header.group(4)
            current_orbital = {
                'label':  label,
                'occ':    occ,
                'energy': energy,
                'aos':    []
            }
            current_orbital['aos'].extend(parse_ao_contributions(rest))
        else:
            m_cont = continuation.match(line)
            if m_cont and current_orbital is not None:
                content = m_cont.group(2)
                current_orbital['aos'].extend(parse_ao_contributions(content))
            else:
                if current_orbital is not None:
                    orbitals.append(current_orbital)
                    current_orbital = None
                if orbitals:
                    break

    if current_orbital is not None:
        orbitals.append(current_orbital)

    return orbitals


# ---------------------------------------------------------------------------
# Shared utilities
# ---------------------------------------------------------------------------

def ao_type_signature(orb, target_atoms):
    """
    Build a dict mapping (cen, typ) -> net coefficient for all AO
    contributions from target atoms in this orbital.
    """
    sig = {}
    for (cen, mu, typ, coeff) in orb['aos']:
        if cen in target_atoms:
            key = (cen, typ)
            sig[key] = sig.get(key, 0.0) + float(coeff)
    return sig


def ao_full_signature(orb, target_atoms):
    """
    Build a dict mapping (cen, mu, typ) -> net coefficient for all AO
    contributions from target atoms in this orbital.
    """
    sig = {}
    for (cen, mu, typ, coeff) in orb['aos']:
        if cen in target_atoms:
            key = (cen, mu, typ)
            sig[key] = sig.get(key, 0.0) + float(coeff)
    return sig


# ---------------------------------------------------------------------------
# Mode 1: Atom search
# ---------------------------------------------------------------------------

def run_atom_search(orbitals, target_atoms, outfile):
    lines_out = []
    lines_out.append(
        f"Searching NATURAL ORBITALS for contributions from atom(s): "
        f"{', '.join(sorted(target_atoms, key=int))}\n"
    )
    lines_out.append(f"Parsed {len(orbitals)} natural orbitals.\n")
    lines_out.append("")

    found_any = False
    for orb in orbitals:
        matching_aos = [
            (cen, mu, typ, coeff)
            for (cen, mu, typ, coeff) in orb['aos']
            if cen in target_atoms
        ]
        if matching_aos:
            found_any = True
            lines_out.append(
                f"  Orbital {orb['label']:>6}  "
                f"occ={orb['occ']:.5f}  "
                f"energy={orb['energy']:.5f}"
            )
            for (cen, mu, typ, coeff) in matching_aos:
                lines_out.append(
                    f"    -> Atom {cen:>4}  mu={mu}  type={typ:>4}  coeff={coeff}"
                )
            lines_out.append("")

    if not found_any:
        lines_out.append("No orbitals found with contributions from the specified atom(s).")

    with open(outfile, 'w') as f:
        f.write('\n'.join(lines_out))
    print(f"Done. Results written to '{outfile}'.")


# ---------------------------------------------------------------------------
# Mode 2: Antibonding partner search
# ---------------------------------------------------------------------------

def build_antibonding_report(orbitals, ref_label, target_atoms, n_results):
    """
    Builds the antibonding-partner report for a single reference orbital.
    Returns a list of output lines (does not write to file).
    Raises ValueError if the orbital is not found or has no AO data on
    the specified target atoms.
    """
    ref_orb = next((o for o in orbitals if o['label'] == ref_label), None)
    if ref_orb is None:
        raise ValueError(f"Orbital '{ref_label}' not found.")

    # Build (cen, mu, typ) signature for reference, restricted to target_atoms
    ref_sig = {}
    for (cen, mu, typ, coeff) in ref_orb['aos']:
        if cen in target_atoms:
            key = (cen, mu, typ)
            ref_sig[key] = ref_sig.get(key, 0.0) + float(coeff)

    if not ref_sig:
        raise ValueError(f"Orbital {ref_label} has no contributions from target atoms.")

    lines_out = []
    lines_out.append(f"Reference natural orbital: {ref_label}")
    lines_out.append(f"  occ={ref_orb['occ']:.5f}  energy={ref_orb['energy']:.5f}")
    lines_out.append(f"  Target atom(s): {', '.join(sorted(target_atoms, key=int))}")
    lines_out.append(f"  Total (atom, mu, type) entries in reference: {len(ref_sig)}")
    lines_out.append("")

    # --- Strict filter: all (cen, mu, typ) opposite sign ---
    exact_matches = []
    scored = []

    for orb in orbitals:
        if orb['occ'] >= 0.05:
            continue
        cand_sig = {}
        for (cen, mu, typ, coeff) in orb['aos']:
            key = (cen, mu, typ)
            cand_sig[key] = cand_sig.get(key, 0.0) + float(coeff)

        shared = set(ref_sig.keys()) & set(cand_sig.keys())
        if not shared:
            continue

        n_opposite = sum(1 for k in shared if ref_sig[k] * cand_sig[k] < 0)
        n_same = len(shared) - n_opposite
        scored.append((n_opposite, -n_same, len(shared), orb, cand_sig, shared))

        # Check strict all-opposite condition
        if all(ref_sig[k] * cand_sig[k] < 0 for k in ref_sig if k in cand_sig) \
                and all(k in cand_sig for k in ref_sig):
            exact_matches.append((orb, cand_sig))

    scored.sort(key=lambda x: (-x[0], x[1], -x[2]))

    # Report exact matches first
    if exact_matches:
        lines_out.append(
            f"Exact matches (all {len(ref_sig)} (atom,mu,type) triples strictly opposite sign):\n"
        )
        for orb, cand_sig in exact_matches:
            lines_out.append(
                f"  Orbital {orb['label']:>6}  occ={orb['occ']:.5f}  energy={orb['energy']:.5f}"
            )
            for key in sorted(ref_sig.keys(), key=lambda k: -abs(ref_sig[k])):
                cen, mu, typ = key
                lines_out.append(
                    f"    (atom {cen:>3}, mu={mu}, {typ:>4})  "
                    f"ref={ref_sig[key]:+.5f}  cand={cand_sig[key]:+.5f}"
                )
            lines_out.append("")
    else:
        lines_out.append("No exact (atom,mu,type) all-opposite-sign matches found.\n")

    # Report top N by score
    lines_out.append(
        f"Top {n_results} virtual orbital candidates ranked by "
        f"number of opposite-sign (atom,mu,type) triples:\n"
    )
    for rank, (n_opp, neg_same, n_shared, orb, cand_sig, shared) in enumerate(scored[:n_results], 1):
        lines_out.append(
            f"  Rank {rank:>3}  |  Orbital {orb['label']:>6}  "
            f"occ={orb['occ']:.5f}  energy={orb['energy']:.5f}"
        )
        lines_out.append(
            f"           |  Shared triples: {n_shared}  "
            f"Opposite: {n_opp}  Same: {-neg_same}"
        )
        for key in sorted(shared, key=lambda k: -abs(ref_sig[k])):
            cen, mu, typ = key
            sign_str = "OPPOSITE" if ref_sig[key] * cand_sig[key] < 0 else "SAME    "
            lines_out.append(
                f"               (atom {cen:>3}, mu={mu}, {typ:>4})  "
                f"ref={ref_sig[key]:+.5f}  cand={cand_sig[key]:+.5f}  [{sign_str}]"
            )
        lines_out.append("")

    return lines_out


def run_antibonding_search(orbitals, ref_label, target_atoms, n_results, outfile):
    """
    Single-orbital convenience wrapper: builds the report for one reference
    orbital and writes it directly to outfile.
    """
    try:
        lines_out = build_antibonding_report(orbitals, ref_label, target_atoms, n_results)
    except ValueError as e:
        print(f"ERROR: {e}")
        sys.exit(1)

    with open(outfile, 'w') as f:
        f.write('\n'.join(lines_out))
    print(f"Done. Results written to '{outfile}'.")


def run_antibonding_search_multi(orbitals, ref_labels, atom_spec, n_results, outfile):
    """
    Multi-orbital antibonding search. ref_labels is a list of orbital labels.
    atom_spec is a dict mapping ref_label -> target_atoms (set of strings),
    OR a single set of target_atoms shared across all ref_labels, OR None
    (meaning: for each orbital, use all atoms present in that orbital).
    Writes one combined output file with a banner section per orbital.
    """
    all_lines = []
    banner_width = 70

    for i, ref_label in enumerate(ref_labels):
        # Resolve target atoms for this specific orbital
        if atom_spec is None:
            ref_orb = next((o for o in orbitals if o['label'] == ref_label), None)
            if ref_orb is None:
                all_lines.append(f"ERROR: Orbital '{ref_label}' not found. Skipping.\n")
                continue
            target_atoms = {cen for (cen, mu, typ, coeff) in ref_orb['aos']}
        elif isinstance(atom_spec, dict):
            target_atoms = atom_spec.get(ref_label)
            if not target_atoms:
                all_lines.append(
                    f"ERROR: No atom restriction provided for orbital '{ref_label}'. Skipping.\n"
                )
                continue
        else:
            # Shared set of atoms across all orbitals
            target_atoms = atom_spec

        if i > 0:
            all_lines.append("")
            all_lines.append("")

        all_lines.append("=" * banner_width)
        all_lines.append(f"  ANTIBONDING SEARCH FOR REFERENCE ORBITAL {ref_label}")
        all_lines.append("=" * banner_width)
        all_lines.append("")

        try:
            section_lines = build_antibonding_report(orbitals, ref_label, target_atoms, n_results)
            all_lines.extend(section_lines)
        except ValueError as e:
            all_lines.append(f"ERROR: {e}")
            all_lines.append("")

    with open(outfile, 'w') as f:
        f.write('\n'.join(all_lines))
    print(f"Done. Results for {len(ref_labels)} orbital(s) written to '{outfile}'.")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    filepath = input("Enter the path to your MOLPRO .out file: ").strip()

    print("\nSelect mode:")
    print("  1) Search for all natural MOs with contributions from specified atoms")
    print("  2) Find antibonding partner candidates for a reference occupied MO")
    mode = input("Enter 1 or 2: ").strip()

    outfile = input("Enter a name for the output file (e.g. results.txt): ").strip()
    if not outfile:
        outfile = "natorb_results2.txt"

    orbitals = parse_natural_orbitals(filepath)
    print(f"Parsed {len(orbitals)} natural orbitals.")

    if mode == '1':
        raw_atoms = input(
            "Enter target atom index/indices (comma-separated, e.g. 40,67,68): "
        ).strip()
        target_atoms = set(raw_atoms.replace(' ', '').split(','))
        run_atom_search(orbitals, target_atoms, outfile)

    elif mode == '2':
        raw_labels = input(
            "Enter the label(s) of the reference occupied (bonding) orbital(s) "
            "you want antibonding partner(s) for. Separate multiple orbitals "
            "with commas (e.g. 88.1,93.1,112.1): "
        ).strip()
        ref_labels = [lbl.strip() for lbl in raw_labels.split(',') if lbl.strip()]
        if not ref_labels:
            print("ERROR: No orbital label(s) provided.")
            sys.exit(1)

        raw_atoms = input(
            "Enter the atom index/indices to restrict the comparison to "
            "(comma-separated, e.g. 40,67,68), or press Enter to use ALL atoms "
            "in each reference orbital. If you entered multiple orbitals above, "
            "this same atom restriction will be applied to all of them: "
        ).strip()
        if raw_atoms:
            target_atoms = set(raw_atoms.replace(' ', '').split(','))
            atom_spec = target_atoms
        else:
            # None -> resolved per-orbital inside run_antibonding_search_multi
            atom_spec = None

        n_raw = input(
            "How many top candidates to report per orbital? (default 20): "
        ).strip()
        n_results = int(n_raw) if n_raw.isdigit() else 20

        if len(ref_labels) == 1:
            # Preserve old single-orbital behaviour/output format
            ref_label = ref_labels[0]
            if atom_spec is None:
                ref_orb = next((o for o in orbitals if o['label'] == ref_label), None)
                if ref_orb is None:
                    print(f"ERROR: Orbital '{ref_label}' not found.")
                    sys.exit(1)
                target_atoms = {cen for (cen, mu, typ, coeff) in ref_orb['aos']}
            else:
                target_atoms = atom_spec
            run_antibonding_search(orbitals, ref_label, target_atoms, n_results, outfile)
        else:
            run_antibonding_search_multi(orbitals, ref_labels, atom_spec, n_results, outfile)

    else:
        print("Unrecognised mode. Please enter 1 or 2.")
        sys.exit(1)


if __name__ == '__main__':
    main()
