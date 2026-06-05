# -*- coding: utf-8 -*-
"""
Converter to PDBQT - PyRx / AutoDock Tools compatible
======================================================
Supported input formats: .pdb, .sdf, .mol2, .cif (mmCIF/PDBx)
Output: .pdbqt files compatible with PyRx, AutoDock Vina, GNINA, Smina

Ligands  -> Meeko + RDKit  (produces ADT-compatible format)
Receptors -> Open Babel    (rigid, polar H only)
CIF receptors -> gemmi (CIF -> PDB) + Open Babel (PDB -> PDBQT)

Dependencies:
    pip install meeko rdkit gemmi

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

    # CIF receptor (mmCIF/PDBx format)
    python3 convert_to_pdbqt.py ./proteins --type receptor --formats cif

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

SUPPORTED_FORMATS = {".pdb", ".sdf", ".mol2", ".cif"}

OBABEL_FORMAT = {
    ".pdb":  "pdb",
    ".sdf":  "sdf",
    ".mol2": "mol2",
    ".cif":  "cif",
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
    # CIF files are almost always macromolecular receptors (PDB/AlphaFold)
    if filepath.suffix.lower() == ".cif":
        return "receptor"
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
#  RECEPTOR CIF: mmCIF/PDBx -> PDBQT
# ============================================================
#
#  Strategy (in order of preference):
#
#  1. obabel direct  : obabel -icif -> -opdbqt
#     Open Babel supports mmCIF natively (format "cif") since v2.4.
#     This is the fastest and most reliable path for PDB/AlphaFold files.
#
#  2. gemmi + obabel : CIF -> intermediate PDB -> PDBQT
#     Used when obabel fails (e.g. old obabel builds, non-standard mmCIF).
#     gemmi reads the mmCIF, writes a clean PDB, then obabel converts.
#
#  3. Error          : both paths failed.

def _obabel_cif_direct(cif_path, output_file, pH=7.4):
    """
    Try to convert CIF -> PDBQT using Open Babel's native mmCIF reader.
    Returns True on success, False on failure (caller should try fallback).
    """
    command = [
        "obabel",
        "-icif", str(cif_path),
        "-opdbqt",
        "-O", str(output_file),
        "--partialcharge", "gasteiger",
        "-xr",          # rigid receptor
        "-xh",          # add polar H only
        "--ph", str(pH),
        "-d",           # delete non-polar H first
    ]
    try:
        result = subprocess.run(command, capture_output=True, text=True, timeout=300)
        if result.returncode == 0 and output_file.exists() and output_file.stat().st_size > 0:
            return True
        log.debug("    obabel direct failed (rc=%d): %s",
                  result.returncode, (result.stderr or result.stdout).strip()[:200])
        return False
    except FileNotFoundError:
        raise   # let caller handle missing obabel
    except Exception as e:
        log.debug("    obabel direct exception: %s", e)
        return False


def _cif_to_pdb_gemmi(cif_path, pdb_path):
    """
    Convert mmCIF/PDBx -> PDB using gemmi (fallback pre-processor).

    Handles macromolecular mmCIF (PDB/AlphaFold) and small-molecule CIF
    (COD/CSD). Strips existing H so Open Babel re-adds them correctly.
    Returns True on success, False on error.
    """
    try:
        import gemmi
    except ImportError:
        log.error("  [FAIL] 'gemmi' not installed (needed as fallback).")
        log.error("         Run: pip install gemmi")
        return False

    # --- macromolecular path (most common for receptors) ---
    try:
        structure = gemmi.read_structure(str(cif_path))

        n_atoms = sum(
            len(res)
            for model in structure
            for chain in model
            for res in chain
        )
        if n_atoms == 0:
            raise ValueError("No atoms found via read_structure")

        structure.remove_hydrogens()

        options = gemmi.PdbWriteOptions()
        options.numbered_ter = True
        options.ter_ignores_missing_atoms = True
        structure.write_pdb(str(pdb_path), options)

        log.debug("    gemmi (macro) wrote PDB: %s  [%d atoms]",
                  pdb_path.name, n_atoms)
        return True

    except Exception as macro_err:
        log.debug("    gemmi macro path failed: %s", macro_err)

    # --- small-molecule CIF fallback (COD/CSD style) ---
    try:
        import gemmi

        doc   = gemmi.cif.read(str(cif_path))
        block = doc.sole_block()
        small = gemmi.make_small_structure_from_block(block)

        if len(small.sites) == 0:
            log.error("  [FAIL] gemmi found no atoms in CIF: %s", cif_path.name)
            return False

        pdb_lines = ["REMARK   CIF converted by gemmi (small-molecule path)"]
        serial = 1
        for site in small.sites:
            elem = site.element.name if site.element != gemmi.Element("X") else "C"
            frac = gemmi.Fractional(site.fract.x, site.fract.y, site.fract.z)
            pos  = small.cell.orthogonalize(frac)
            name = (site.label or elem)[:4].ljust(4)
            pdb_lines.append(
                "HETATM%5d %-4s LIG     1    %8.3f%8.3f%8.3f"
                "  1.00  0.00          %2s  "
                % (serial, name, pos.x, pos.y, pos.z, elem)
            )
            serial += 1

        pdb_lines.append("END")
        pdb_path.write_text("\n".join(pdb_lines) + "\n", encoding="utf-8")
        log.debug("    gemmi (small-mol) wrote PDB: %s", pdb_path.name)
        return True

    except Exception as small_err:
        log.error("  [FAIL] gemmi could not parse CIF '%s': %s",
                  cif_path.name, small_err)
        return False


def convert_receptor_cif(cif_path, output_folder, pH=7.4):
    """
    Convert mmCIF/PDBx receptor to PDBQT.

    Tries two strategies:
      1. obabel direct  (CIF -> PDBQT, no intermediate file)
      2. gemmi + obabel (CIF -> PDB -> PDBQT, intermediate PDB deleted after)

    Raises SystemExit only for a missing obabel installation.
    """
    output_file = output_folder / (cif_path.stem + ".pdbqt")
    log.info("  [CIF] %s", cif_path.name)

    # ── Strategy 1: obabel direct ──────────────────────────────────────────
    try:
        if _obabel_cif_direct(cif_path, output_file, pH):
            size_kb = output_file.stat().st_size // 1024
            log.info("  [OK] %s -> %s  [%d KB, rigid receptor, direct]",
                     cif_path.name, output_file.name, size_kb)
            return True
        log.debug("    obabel direct did not produce output; trying gemmi fallback")
    except FileNotFoundError:
        log.error("  [FAIL] 'obabel' not found.")
        log.error("     Windows: https://openbabel.org/wiki/Install")
        log.error("     Linux:   sudo apt install openbabel")
        log.error("     Mac:     brew install open-babel")
        return False

    # ── Strategy 2: gemmi -> PDB -> obabel ────────────────────────────────
    log.debug("    falling back to gemmi + obabel pipeline")
    tmp_pdb = output_folder / (cif_path.stem + "__tmp_cif.pdb")

    try:
        if not _cif_to_pdb_gemmi(cif_path, tmp_pdb):
            return False

        if not tmp_pdb.exists() or tmp_pdb.stat().st_size == 0:
            log.error("  [FAIL] gemmi produced an empty PDB for %s", cif_path.name)
            return False

        ok = convert_receptor_obabel(tmp_pdb, output_folder, pH)

        # obabel names output after tmp_pdb stem; rename to final name
        obabel_out = output_folder / tmp_pdb.with_suffix(".pdbqt").name
        if ok and obabel_out.exists() and obabel_out != output_file:
            obabel_out.replace(output_file)
            size_kb = output_file.stat().st_size // 1024
            log.info("  [OK] %s -> %s  [%d KB, rigid receptor, gemmi+obabel]",
                     cif_path.name, output_file.name, size_kb)
        return ok

    finally:
        if tmp_pdb.exists():
            tmp_pdb.unlink()


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
            if f.suffix.lower() == ".cif":
                ok = convert_receptor_cif(f, output_path, pH)
            else:
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

  # Receptor/protein PDB/MOL2 (Open Babel)
  python convert_to_pdbqt.py C:/Bioinfo/proteins --type receptor

  # Receptor mmCIF/PDBx (gemmi + Open Babel)
  python convert_to_pdbqt.py C:/Bioinfo/proteins --type receptor --formats cif

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
                        choices=["pdb", "sdf", "mol2", "cif"],
                        default=None,
                        help="Formats to convert (default: all)\n"
                             "Example: --formats sdf mol2 cif")
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
