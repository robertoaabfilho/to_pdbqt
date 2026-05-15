# -*- coding: utf-8 -*-
"""
Converter to PDBQT - PyRx / AutoDock Tools compatible
======================================================
Supported input formats: .pdb, .sdf, .mol2
Output: .pdbqt files compatible with PyRx, AutoDock Vina, GNINA, Smina

Ligands  -> Meeko + RDKit  (produces ADT-compatible format)
Receptors -> Open Babel    (rigid, polar H only)

Dependencies:
    pip install meeko rdkit

    Open Babel (for receptors):
      Windows: https://openbabel.org/wiki/Install  (check "Add to PATH")
      Linux:   sudo apt install openbabel
      Mac:     brew install open-babel

Usage:
    # Ligands (auto-detected or explicit)
    python3 convert_to_pdbqt.py ./ligands
    python3 convert_to_pdbqt.py ./ligands --type ligand

    # Receptor/protein
    python3 convert_to_pdbqt.py ./proteins --type receptor

    # Filter by format
    python3 convert_to_pdbqt.py ./ligands --formats sdf mol2

    # With separate output folder
    python3 convert_to_pdbqt.py ./ligands ./output
"""

import sys
import subprocess
import argparse
import logging
from pathlib import Path

logging.basicConfig(level=logging.INFO, format="%(message)s")
log = logging.getLogger(__name__)

SUPPORTED_FORMATS = {".pdb", ".sdf", ".mol2"}

OBABEL_FORMAT = {
    ".pdb":  "pdb",
    ".sdf":  "sdf",
    ".mol2": "mol2",
}

# Amino acids for auto-detection of receptors
AMINO_ACIDS = {
    "ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE",
    "LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL",
    "HID","HIE","HIP","HSE","HSP","HSD","CYX","ACE","NME",
}


# ============================================================
#  AutoDock atom type assignment (replicates ADT logic)
# ============================================================

def get_autodock_type(atom, mol):
    """
    Assign AutoDock atom type following ADT rules.
    C/A, N/NA, OA, S/SA, H/HD, F, Cl, Br, I, P
    """
    symbol = atom.GetSymbol()

    if symbol == "C":
        return "A" if atom.GetIsAromatic() else "C"

    if symbol == "N":
        if atom.GetTotalNumHs() == 0:
            return "NA"
        return "N"

    if symbol == "O":
        return "OA"

    if symbol == "S":
        for nb in atom.GetNeighbors():
            if nb.GetSymbol() == "O":
                return "S"
        return "SA"

    if symbol == "H":
        for nb in atom.GetNeighbors():
            if nb.GetSymbol() in ("O", "N", "S"):
                return "HD"
        return "H"

    if symbol == "P":
        return "P"

    if symbol in ("F", "Cl", "Br", "I"):
        return symbol

    return symbol


# ============================================================
#  Auto-detect receptor vs ligand
# ============================================================

def detect_molecule_type(filepath):
    """Returns 'receptor' or 'ligand' based on residue content."""
    if filepath.suffix.lower() in (".sdf", ".mol2"):
        return "ligand"
    aa_count = 0
    total = 0
    try:
        with open(filepath, "r", errors="ignore") as fh:
            for line in fh:
                if line.startswith(("ATOM", "HETATM")):
                    resname = line[17:20].strip().upper()
                    if resname in AMINO_ACIDS:
                        aa_count += 1
                    total += 1
                if total >= 200:
                    break
    except Exception:
        return "ligand"
    if total == 0:
        return "ligand"
    return "receptor" if (aa_count / total) > 0.5 else "ligand"


# ============================================================
#  Read molecule with RDKit
# ============================================================

def read_molecule(filepath):
    """Read PDB / SDF / MOL2 with RDKit. Returns mol with explicit Hs."""
    from rdkit import Chem
    from rdkit.Chem import AllChem

    ext = filepath.suffix.lower()

    if ext == ".sdf":
        suppl = Chem.SDMolSupplier(str(filepath), removeHs=False, sanitize=True)
        mols = [m for m in suppl if m is not None]
        mol = mols[0] if mols else None
    elif ext == ".mol2":
        mol = Chem.MolFromMol2File(str(filepath), removeHs=False, sanitize=True)
    elif ext == ".pdb":
        mol = Chem.MolFromPDBFile(str(filepath), removeHs=False, sanitize=True)
    else:
        return None

    if mol is None:
        return None

    mol = Chem.AddHs(mol, addCoords=True)

    if mol.GetNumConformers() == 0:
        params = AllChem.ETKDGv3()
        params.randomSeed = 42
        if AllChem.EmbedMolecule(mol, params) == -1:
            AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
        AllChem.MMFFOptimizeMolecule(mol)

    return mol


# ============================================================
#  LIGAND: Meeko -> post-process to ADT format
# ============================================================

def convert_ligand_meeko(filepath, output_folder, pH=7.4):
    """
    Convert ligand using Meeko, then post-process the PDBQT
    to match the ADT format expected by PyRx:
      - REMARK 4 header
      - Gasteiger charges recalculated with RDKit
      - Correct AutoDock atom types (OA, NA, A, etc.)
      - CONECT records
    """
    try:
        from rdkit.Chem import rdPartialCharges
        from meeko import MoleculePreparation, PDBQTWriterLegacy

        mol = read_molecule(filepath)
        if mol is None:
            log.error("  [FAIL] Could not read: %s", filepath.name)
            return False

        # Step 1: use Meeko to get rotatable bonds / torsion info
        preparator = MoleculePreparation(
            merge_these_atom_types=["H"],
            hydrate=False,
            flexible_amides=False,
            rigidify_bonds_smarts=[],
        )
        prepared_mols = preparator.prepare(mol)
        if not prepared_mols:
            log.error("  [FAIL] Meeko could not prepare: %s", filepath.name)
            return False

        # Step 2: compute Gasteiger charges with RDKit on the same mol
        rdPartialCharges.ComputeGasteigerCharges(mol)

        # Step 3: write PDBQT in ADT format
        output_file = output_folder / filepath.with_suffix(".pdbqt").name
        conf = mol.GetConformer()
        lines = []
        lines.append("REMARK   4 XXXX COMPLIES WITH FORMAT V. 2.0")

        atom_serial = {}
        serial = 1

        for atom in mol.GetAtoms():
            idx    = atom.GetIdx()
            symbol = atom.GetSymbol()
            pos    = conf.GetAtomPosition(idx)

            charge = atom.GetDoubleProp("_GasteigerCharge") \
                     if atom.HasProp("_GasteigerCharge") else 0.0
            if charge != charge or abs(charge) > 9.999:
                charge = 0.0

            ad_type = get_autodock_type(atom, mol)

            pdb_info = atom.GetPDBResidueInfo()
            if pdb_info and pdb_info.GetName().strip():
                atom_name = pdb_info.GetName().strip()
            else:
                atom_name = symbol

            line = "HETATM%5d %-4s UNK     0    %8.3f%8.3f%8.3f  1.00  0.00    %6.3f %s" % (
                serial,
                atom_name,
                pos.x, pos.y, pos.z,
                charge,
                ad_type,
            )
            lines.append(line)
            atom_serial[idx] = serial
            serial += 1

        # CONECT records
        adjacency = {}
        for bond in mol.GetBonds():
            si = atom_serial[bond.GetBeginAtomIdx()]
            sj = atom_serial[bond.GetEndAtomIdx()]
            adjacency.setdefault(si, []).append(sj)
            adjacency.setdefault(sj, []).append(si)

        for atom_s in sorted(adjacency.keys()):
            neighbors = sorted(adjacency[atom_s])
            for i in range(0, len(neighbors), 4):
                chunk = neighbors[i:i+4]
                conect = "CONECT%5d" % atom_s
                for n in chunk:
                    conect += "%5d" % n
                lines.append(conect)

        output_file.write_text("\r\n".join(lines) + "\r\n", encoding="utf-8")

        n_heavy   = sum(1 for a in mol.GetAtoms() if a.GetSymbol() != "H")
        torsions  = sum(1 for a in mol.GetAtoms()
                        if a.GetSymbol() != "H" and a.GetDegree() > 1)
        log.info("  [OK] %s -> %s  [%d heavy atoms]",
                 filepath.name, output_file.name, n_heavy)
        return True

    except ImportError as e:
        log.error("  [FAIL] Missing dependency: %s", e)
        log.error("         Run: pip install meeko rdkit")
        return False
    except Exception as e:
        log.error("  [FAIL] Unexpected error in %s: %s", filepath.name, e)
        return False


# ============================================================
#  RECEPTOR: Open Babel (rigid, polar H only)
# ============================================================

def convert_receptor_obabel(filepath, output_folder, pH=7.4):
    """
    Convert receptor/protein to PDBQT using Open Babel.
    Keeps structure rigid, adds polar hydrogens only.
    """
    ext = filepath.suffix.lower()
    fmt_in = OBABEL_FORMAT.get(ext, "pdb")
    output_file = output_folder / filepath.with_suffix(".pdbqt").name

    command = [
        "obabel",
        "-i" + fmt_in, str(filepath),
        "-opdbqt",
        "-O", str(output_file),
        "--partialcharge", "gasteiger",
        "-xr",
        "-xh",
        "--ph", str(pH),
        "-d",
    ]

    try:
        result = subprocess.run(command, capture_output=True, text=True, timeout=300)
        if result.returncode == 0 and output_file.exists():
            size_kb = output_file.stat().st_size // 1024
            log.info("  [OK] %s -> %s  [%d KB, rigid receptor]",
                     filepath.name, output_file.name, size_kb)
            return True
        else:
            error = result.stderr.strip() or result.stdout.strip()
            log.error("  [FAIL] Error in %s: %s", filepath.name, error)
            return False

    except FileNotFoundError:
        log.error("  [FAIL] 'obabel' not found.")
        log.error("     Windows: https://openbabel.org/wiki/Install")
        log.error("     Linux:   sudo apt install openbabel")
        log.error("     Mac:     brew install open-babel")
        return False
    except subprocess.TimeoutExpired:
        log.error("  [FAIL] Timeout processing %s", filepath.name)
        return False
    except Exception as e:
        log.error("  [FAIL] Unexpected error in %s: %s", filepath.name, e)
        return False


# ============================================================
#  File collection and folder processing
# ============================================================

def collect_files(folder, formats, recursive):
    files = []
    pattern = "**/*" if recursive else "*"
    for f in sorted(folder.glob(pattern)):
        if f.is_file() and f.suffix.lower() in formats:
            files.append(f)
    return files


def convert_folder(
    input_folder,
    output_folder=None,
    mol_type="auto",
    formats=None,
    recursive=False,
    pH=7.4,
):
    input_path = Path(input_folder).resolve()
    if not input_path.exists():
        log.error("Error: folder '%s' not found.", input_folder)
        sys.exit(1)

    output_path = Path(output_folder).resolve() if output_folder else input_path
    output_path.mkdir(parents=True, exist_ok=True)

    if formats:
        exts = {"." + f.lstrip(".") for f in formats}
        invalid = exts - SUPPORTED_FORMATS
        if invalid:
            log.error("Unsupported formats: %s  Valid: %s", invalid, SUPPORTED_FORMATS)
            sys.exit(1)
    else:
        exts = SUPPORTED_FORMATS

    files = collect_files(input_path, exts, recursive)

    if not files:
        log.warning("No files found in: %s", input_path)
        return

    count_by_fmt = {}
    for f in files:
        count_by_fmt[f.suffix.lower()] = count_by_fmt.get(f.suffix.lower(), 0) + 1

    sep = "=" * 56
    log.info("\n%s", sep)
    log.info("  Converter -> PDBQT  (PyRx / ADT compatible)")
    log.info("  Mol. type : %s", mol_type)
    log.info("  pH        : %s", pH)
    log.info("  Input     : %s", input_path)
    log.info("  Output    : %s", output_path)
    log.info("  Files     : %d  %s", len(files),
             "  ".join("%s: %d" % (e, n) for e, n in sorted(count_by_fmt.items())))
    log.info(sep)

    success, failed, skipped = 0, 0, 0

    for f in files:
        expected = output_path / f.with_suffix(".pdbqt").name
        if expected.exists() and expected.stat().st_mtime >= f.stat().st_mtime:
            log.info("  [SKIP] %s (already converted)", f.name)
            skipped += 1
            continue

        if mol_type == "auto":
            detected = detect_molecule_type(f)
        else:
            detected = mol_type

        if detected == "receptor":
            ok = convert_receptor_obabel(f, output_path, pH)
        else:
            ok = convert_ligand_meeko(f, output_path, pH)

        if ok:
            success += 1
        else:
            failed += 1

    log.info("\n%s", sep)
    summary = "  Result: [OK] %d converted" % success
    if skipped:
        summary += "  [SKIP] %d skipped" % skipped
    if failed:
        summary += "  [FAIL] %d errors" % failed
    log.info(summary)
    log.info(sep + "\n")

    if failed:
        sys.exit(1)


# ============================================================
#  Command line
# ============================================================

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Convert PDB/SDF/MOL2 to PDBQT (PyRx / ADT compatible)",
        formatter_class=argparse.RawTextHelpFormatter,
        epilog="""
Examples:
  # Ligands (Meeko, ADT format)
  python convert_to_pdbqt.py C:/Bioinfo/ligands
  python convert_to_pdbqt.py C:/Bioinfo/ligands --type ligand --formats sdf

  # Receptor/protein (Open Babel)
  python convert_to_pdbqt.py C:/Bioinfo/proteins --type receptor

  # Auto-detect type, separate output folder
  python convert_to_pdbqt.py C:/Bioinfo/input C:/Bioinfo/output

  # Recursive search
  python convert_to_pdbqt.py C:/Bioinfo/all --recursive
        """
    )

    parser.add_argument("input",
                        help="Folder with ligand/receptor files")
    parser.add_argument("output", nargs="?", default=None,
                        help="Output folder (default: same as input)")
    parser.add_argument("--type", dest="mol_type",
                        choices=["ligand", "receptor", "auto"],
                        default="auto",
                        help="Molecule type:\n"
                             "  ligand   -> small molecule, Meeko [default: auto]\n"
                             "  receptor -> protein, Open Babel\n"
                             "  auto     -> detect automatically")
    parser.add_argument("--formats", nargs="+",
                        choices=["pdb", "sdf", "mol2"],
                        default=None,
                        help="Formats to convert (default: all)\n"
                             "Example: --formats sdf mol2")
    parser.add_argument("--pH", type=float, default=7.4,
                        help="pH for protonation (default: 7.4)")
    parser.add_argument("--recursive", action="store_true",
                        help="Search subfolders as well")
    parser.add_argument("--verbose", action="store_true",
                        help="Show detailed debug messages")

    args = parser.parse_args()

    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    convert_folder(
        input_folder=args.input,
        output_folder=args.output,
        mol_type=args.mol_type,
        formats=args.formats,
        recursive=args.recursive,
        pH=args.pH,
    )
