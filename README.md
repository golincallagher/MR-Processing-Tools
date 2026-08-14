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

## scrub_molden_sym.py

Removes `Sym=` lines from a MOLPRO-generated `.molden` file. MOLPRO writes a `Sym=` label on each molecular orbital line in `.molden` output, which some visualization software (e.g. VMD) cannot parse — this causes orbital data to appear unreadable. This script strips all `Sym=` lines and writes the result to a new file, leaving the original file untouched.

### Requirements

- Python 3
- No external dependencies (uses only the standard library: `sys`)

### Usage

```bash
python scrub_molden_sym.py <input.molden> <output.molden>
```

Example:
```bash
python scrub_molden_sym.py edgn2n2coms5cas1.molden edgn2n2coms5casorbitals1.molden
```

### Input format

Any `.molden` file. No specific MOLPRO directive is required — the script simply removes any line containing the substring `Sym=`.

### Output

A copy of the input file with all `Sym=` lines removed, written to the specified output path. The original file is not modified.

## orbital_character_analyzer.py

Quantifies the s/p/d orbital character of a specific atom's basis functions within each molecular orbital (MO) of a MOLPRO-generated `.molden` file. This was written for **ligand-field systems with a single well-defined metal center** — e.g. isolated Fe centers in metal-organic frameworks or single-atom catalysts (SACs) — where the main question of interest is how much d-orbital character a given MO has on the metal, and how that varies across the MO manifold (useful for identifying metal-based vs. ligand-based orbitals, e.g. when setting up a CASSCF active space).

It is **not** a general-purpose orbital population analysis tool — the "d-orbital character" framing assumes there's a specific atom of chemical interest (typically the metal) whose character you want to track, rather than a systematic breakdown across every atom in the system.

### Requirements

- Python 3
- No external dependencies (uses only the standard library: `sys`)
- The `.molden` file must contain a `[GTO]` section (standard for MOLPRO's `.molden` output) — the script parses this directly to determine each atom's basis function layout, so no manual counting of basis functions is required.

### Usage

```bash
python orbital_character_analyzer.py
```

You'll be prompted for:
1. The path to your `.molden` file
2. Which atom to analyze — use the same numbering as the file's `[Atoms]` block (e.g. `1` for the first atom)
3. How many molecular orbitals to analyze (starting from MO 1)

The script prints the detected basis function breakdown for the chosen atom (total functions, and the s/p/d count) before running the analysis — **check this against what you expect for that atom's basis set** before trusting the results, especially the first time you run it on a new system or basis set.

### Output

For each MO analyzed, prints:
- MO index, orbital energy, and occupation number
- The percentage of that atom's own character that is s, p, or d type
- The atom's total contribution to that MO, as a percentage of the whole MO

### How it works

`.molden` files store each MO's AO coefficients as one flat list, ordered atom-by-atom to match the `[Atoms]` block. To isolate a given atom's coefficients from that list, the script needs to know how many basis functions belong to every atom *before* it. Rather than requiring this to be entered manually, the script parses the `[GTO]` block to reconstruct, for every atom, how many s/p/d/f/g shells it has and how many basis functions each shell contributes — then computes the correct offset automatically.

d/f/g shells can be represented as either spherical harmonics or Cartesian Gaussians, which changes how many functions they contribute (5 vs. 6 for d, for example). The script checks for a `[5D]`, `[7F]`, or `[9G]` tag elsewhere in the file to determine which convention is in use, defaulting to Cartesian if no such tag is found.

### Notes

- **Only tested on systems where the metal center is the first atom listed in `[Atoms]`.** The atom-offset logic is written to generalize to any atom position, but has not yet been verified against a system where the atom of interest appears later in the atom list. If you use this on such a system, it's worth manually spot-checking the reported basis function count and one or two MOs against an independent source before trusting the output.
- If the atom of interest has f or g shells (uncommon for typical MOF/SAC transition-metal basis sets, but possible with polarized or larger basis sets), the reported s/p/d percentages will not sum to 100% — the remainder reflects f/g character, which isn't broken out as a separate column.
- The script assumes atom numbering in the `[GTO]` block matches the `[Atoms]` block exactly (standard for MOLPRO output). If you're working with a `.molden` file generated by different software, double check this assumption holds.
