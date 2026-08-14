# MR-Processing-Tools

Some scripts developed during my PhD to help process `.molden` files and other outputs from multireference wavefunction-based calculations.

## parse_natorb.py

Parses the `NATURAL ORBITALS` block from a MOLPRO output file (`.out`) and reports molecular orbitals of interest, either by atom contribution or by antibonding-partner matching.

### Requirements

- Python 3
- No external dependencies (uses only the standard library: `re`, `sys`)

### Usage

Run the script and follow the interactive prompts:

```bash
python parse_natorb.py
```

You'll be asked for:
1. The path to your MOLPRO `.out` file
2. Which mode to run (see below)
3. A name for the output file (defaults to `natorb_results2.txt` if left blank)

### Modes

**1) Atom search**
Finds all natural orbitals with any AO contribution from a set of user-specified atom indices.

- Prompts for a comma-separated list of atom indices (e.g. `40,67,68`)
- Writes every matching orbital, its occupation number, orbital energy, and the individual AO contributions from those atoms to the output file

**2) Antibonding partner search**
Given a reference occupied (bonding) orbital, finds virtual orbital candidates that are likely antibonding partners, ranked by how many `(atom, mu, orbital type)` triples have strictly opposite-sign coefficients relative to the reference.

- Prompts for one or more reference orbital labels (e.g. `88.1` or `88.1,93.1,112.1`)
- Prompts for which atoms to restrict the comparison to (or press Enter to use all atoms present in the reference orbital)
- Prompts for how many top-ranked candidates to report per orbital (default 20)
- Only orbitals with occupation number below 0.05 are considered as antibonding candidates
- Reports any "exact matches" (all shared triples opposite in sign) first, followed by a ranked list of top candidates
- If multiple reference orbitals are given, the output file contains one banner-separated section per orbital

### Input format

The script expects a MOLPRO `.out` file containing a `NATURAL ORBITALS` section (i.e. the job must include a `NATORB` directive). If this section isn't found, the script exits with an error message.

### Output

A plain-text file (name specified at runtime) listing the orbitals matching the selected search mode, along with occupation numbers, orbital energies, and AO contribution breakdowns.

### Notes

- Atom indices are matched as strings, so make sure the indices you enter match the atom numbering exactly as it appears in the MOLPRO output.
- For antibonding partner search, "strict opposite sign" matching requires every `(atom, mu, type)` triple present in the reference orbital to also appear in the candidate orbital with an opposite-sign coefficient — this is a stricter criterion than the ranked "top candidates" list, which allows partial overlap.
