"""
orbital_character_analyzer.py
------------------------------
Automatically determines the s/p/d orbital character of a specified
atom's basis functions within each molecular orbital (MO) of a
MOLPRO-generated .molden file.

Originally written to quantify metal d-orbital character on an
isolated metal center in CASSCF/CASPT2 calculations on Fe-containing
metal-organic framework systems, but works for any atom in any basis
set.

Earlier versions of this script required the user to manually count
basis functions atom-by-atom to work out where a given atom's
coefficients start within each MO's flat coefficient list. This version
instead parses the [GTO] block of the .molden file directly to
reconstruct each atom's basis function layout (how many s/p/d/f/g
shells it has, and how many functions each shell contributes), so no
manual counting is required -- you only need to specify which atom to
analyze.

Usage:
    python orbital_character_analyzer.py
    (prompts interactively for the .molden file path, which atom to
    analyze -- using the same 1-based numbering as the file's [Atoms]
    block -- and how many MOs to analyze)
"""

import sys


def detect_shell_dimensions(all_lines):
    """
    Molden files can use either spherical harmonic or Cartesian
    Gaussians for d/f/g shells, which changes how many basis functions
    each shell contributes (5 vs 6 for d, 7 vs 10 for f, 9 vs 15 for g).
    This is signaled by a [5D], [7F], [5D7F], [5D10F], or [9G] tag
    somewhere in the file. Defaults to Cartesian (6D/10F/15G) if no
    such tag is found.

    Only lines that are themselves bracketed section tags (e.g. a line
    that is exactly "[5D]") are checked -- NOT a blind substring search
    over the whole file. MOLPRO writes exponents in Fortran D-notation
    (e.g. "0.2740600135D-01"), and a naive substring search for "5d"
    will spuriously match inside such numbers (".....35D-01" contains
    "5d"), incorrectly flipping the file into spherical mode even when
    no such tag exists.
    """
    tags = {line.strip().lower() for line in all_lines if line.strip().startswith('[')}
    d_dim = 5 if any('5d' in tag for tag in tags) else 6
    f_dim = 7 if any('7f' in tag for tag in tags) else 10
    g_dim = 9 if any('9g' in tag for tag in tags) else 15
    return {'s': 1, 'p': 3, 'd': d_dim, 'f': f_dim, 'g': g_dim}


def parse_gto_section(all_lines):
    """
    Parses the [GTO] block of a .molden file. Returns a dict keyed by
    atom index (matching the [Atoms] block numbering), where each value
    is a dict of {shell_type: number_of_shells_of_that_type} for that
    atom, e.g. {'s': 3, 'p': 2, 'd': 1}.
    """
    start = None
    for i, line in enumerate(all_lines):
        if line.strip().lower() == '[gto]':
            start = i + 1
            break

    if start is None:
        print("ERROR: '[GTO]' section not found in molden file.")
        sys.exit(1)

    atoms = {}
    current_index = None
    i = start
    n = len(all_lines)

    while i < n:
        stripped = all_lines[i].strip()

        if stripped == '':
            i += 1
            continue

        if stripped.startswith('['):
            break  # left the [GTO] section

        tokens = stripped.split()

        # Atom header line, e.g. "1 0" (atom index, default charge)
        if len(tokens) == 2 and tokens[0].isdigit() and tokens[1].lstrip('-').isdigit():
            current_index = int(tokens[0])
            atoms[current_index] = {}
            i += 1
            continue

        # Shell header line, e.g. "s   3 1.00" (shell type, # primitives, scale)
        if len(tokens) >= 2 and tokens[0].isalpha() and current_index is not None:
            shell_type = tokens[0].lower()
            n_prim = int(float(tokens[1]))
            atoms[current_index][shell_type] = atoms[current_index].get(shell_type, 0) + 1
            i += 1 + n_prim  # skip the n_prim exponent/coefficient lines that follow
            continue

        i += 1

    return atoms


def build_atom_basis_table(all_lines):
    """
    Combines shell counts (from parse_gto_section) with shell
    dimensions (from detect_shell_dimensions) to produce, per atom,
    the total number of basis functions and the breakdown by s/p/d/f/g.
    Returns: {atom_index: {'total':.., 's':.., 'p':.., 'd':.., 'f':.., 'g':..}}
    """
    dims = detect_shell_dimensions(all_lines)
    shells_by_atom = parse_gto_section(all_lines)

    table = {}
    for atom_index, shells in shells_by_atom.items():
        s_count = shells.get('s', 0) * dims['s']
        p_count = shells.get('p', 0) * dims['p']
        d_count = shells.get('d', 0) * dims['d']
        f_count = shells.get('f', 0) * dims['f']
        g_count = shells.get('g', 0) * dims['g']
        total = s_count + p_count + d_count + f_count + g_count
        table[atom_index] = {
            'total': total, 's': s_count, 'p': p_count, 'd': d_count,
            'f': f_count, 'g': g_count,
        }
    return table


def analyze_mos(all_lines, n_preceding, total_funcs, ns, np_, nd, max_mos):
    """
    For each of the first max_mos molecular orbitals, computes the
    s/p/d character (as a percentage of the target atom's own
    coefficients) and the target atom's total contribution to that MO
    (as a percentage of the whole MO).

    Note: if the target atom has f or g shells, s/p/d percentages will
    not sum to 100% -- the remainder is f/g character, which is not
    broken out separately here.
    """
    in_mo_section = False
    mo_count = 0
    line_index = 0
    n_lines = len(all_lines)

    results = []

    print("\nMO   Energy         Occup    s(%)    p(%)    d(%)   AtomTot(%)")

    while line_index < n_lines:
        line = all_lines[line_index].strip().lower()

        if line.startswith("[mo]"):
            in_mo_section = True
            line_index += 1
            continue

        if in_mo_section and line.startswith("[") and not line.startswith("[mo]"):
            break

        if not in_mo_section:
            line_index += 1
            continue

        if line.startswith("ene="):
            mo_count += 1
            if mo_count > max_mos:
                break

            energy_line = all_lines[line_index].strip()
            mo_energy = float(energy_line.split("=", 1)[1].strip())

            mo_occup = None
            all_coeffs = []

            line_index += 1

            while line_index < n_lines:
                sub_line = all_lines[line_index].strip()
                sub_lower = sub_line.lower()

                if sub_lower.startswith("[") and not sub_lower.startswith("[mo]"):
                    break
                if sub_lower.startswith("ene="):
                    break

                if sub_lower.startswith("occup="):
                    mo_occup = float(sub_line.split("=", 1)[1].strip())
                elif sub_lower.startswith("spin="):
                    pass  # not used in the s/p/d breakdown
                else:
                    tokens = sub_line.split()
                    try:
                        if len(tokens) == 2:
                            all_coeffs.append(float(tokens[1]))
                        elif len(tokens) == 1:
                            all_coeffs.append(float(tokens[0]))
                    except ValueError:
                        pass

                line_index += 1

            if mo_occup is None:
                mo_occup = 0.0

            # Extract this atom's slice of coefficients using the offset
            atom_coeffs = all_coeffs[n_preceding: n_preceding + total_funcs]

            if len(atom_coeffs) < total_funcs:
                print(f"WARNING: MO {mo_count} has only {len(atom_coeffs)} coefficients "
                      f"in the expected atom range (expected {total_funcs}). Skipping.")
                continue

            s_sq = sum(c * c for c in atom_coeffs[0:ns])
            p_sq = sum(c * c for c in atom_coeffs[ns:ns + np_])
            d_sq = sum(c * c for c in atom_coeffs[ns + np_: ns + np_ + nd])

            # Normalize against the atom's FULL coefficient slice (not just
            # s+p+d), so any f/g character is correctly reflected as the
            # gap between 100% and (s_frac + p_frac + d_frac).
            atom_total_sq = sum(c * c for c in atom_coeffs)
            total_sq = sum(c * c for c in all_coeffs)

            if abs(atom_total_sq) < 1e-12:
                s_frac = p_frac = d_frac = 0.0
            else:
                s_frac = 100.0 * s_sq / atom_total_sq
                p_frac = 100.0 * p_sq / atom_total_sq
                d_frac = 100.0 * d_sq / atom_total_sq

            if abs(total_sq) < 1e-12:
                atom_total_frac = 0.0
            else:
                atom_total_frac = 100.0 * (atom_total_sq / total_sq)

            print(f"{mo_count:3d}  {mo_energy:12.6f}  {mo_occup:6.3f}"
                  f"  {s_frac:6.2f}  {p_frac:6.2f}  {d_frac:6.2f}"
                  f"   {atom_total_frac:6.2f}")

            results.append({
                'mo': mo_count,
                'energy': mo_energy,
                'occup': mo_occup,
                's_frac': s_frac,
                'p_frac': p_frac,
                'd_frac': d_frac,
                'atom_total_frac': atom_total_frac,
            })
        else:
            line_index += 1
            continue

    print(f"\nDone. Analyzed {len(results)} molecular orbital(s).")
    return results


def main():
    molden_file = input("Enter the path to your MOLPRO .molden file: ").strip()

    with open(molden_file, 'r') as f:
        all_lines = f.readlines()

    basis_table = build_atom_basis_table(all_lines)
    print(f"\nDetected {len(basis_table)} atom(s) in the [GTO] section.")

    atom_index = int(input(
        "Which atom would you like to analyze? (use the same numbering as "
        "the file's [Atoms] block, e.g. 1 for the first atom): "
    ))

    if atom_index not in basis_table:
        print(f"ERROR: Atom {atom_index} not found in the parsed basis table.")
        sys.exit(1)

    n_preceding = sum(
        info['total'] for idx, info in basis_table.items() if idx < atom_index
    )
    info = basis_table[atom_index]
    total_funcs, ns, np_, nd = info['total'], info['s'], info['p'], info['d']

    extra = ""
    if info['f']:
        extra += f", f={info['f']}"
    if info['g']:
        extra += f", g={info['g']}"
    print(f"\nAtom {atom_index}: {total_funcs} total basis functions "
          f"(s={ns}, p={np_}, d={nd}{extra})")

    max_mos = int(input("\nHow many molecular orbitals should be analyzed? "))

    print(f"\nAnalyzing atom {atom_index} across up to {max_mos} MO(s)...")
    analyze_mos(all_lines, n_preceding, total_funcs, ns, np_, nd, max_mos)


if __name__ == '__main__':
    main()