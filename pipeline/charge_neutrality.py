"""
Charge Neutrality Preflight Check.

This module implements charge neutrality verification before grompp:
- Parse ITP files to extract [ atoms ] section charges
- Compute total system charge: Q = Σ(N_i * q_molecule_i) using numerically stable fsum
- Auto-correct small imbalances using safe-subset atom selection (NOT uniform smearing)
- Fail fast on larger imbalances or unsafe correction targets

CRITICAL PHYSICS SAFEGUARDS:
- Correction applies only to atoms passing safe-subset criteria (low |q|, allowed names)
- Never auto-selects targets unless from explicit allowlist
- Uses math.fsum for all charge summation to avoid catastrophic cancellation
- Refuses correction if preprocessor directives detected in [ atoms ]
- Patches exactly one moleculetype; fails if target is missing (no fallback)

Per AGENT_SKILLS/06_ITP_SANITIZER.md and manifest schema requirements.
"""

import math
import re
import os
import tempfile
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Set, Tuple


class ChargeNeutralityError(Exception):
    """Error during charge neutrality check."""
    pass


@dataclass
class MoleculeCharge:
    """Charge information for a molecule."""
    name: str
    itp_moleculetype_name: str
    itp_path: Path
    total_charge: float
    atom_count: int
    charges: List[float]
    # Enhanced fields for robustness
    atom_details: Optional[List[Dict]] = None  # [{nr, atomname, atomtype, charge}, ...]
    has_preprocessor: bool = False  # True if #ifdef/#include/#define in atoms section
    unparsed_atoms_lines: Optional[List[Dict[str, Any]]] = None  # [{line_no, excerpt, reason}, ...]


@dataclass
class DipoleResult:
    """Result of dipole moment calculation before/after charge correction."""
    dipole_before_debye: Optional[float] = None  # μ before correction
    dipole_after_debye: Optional[float] = None   # μ after correction
    delta_dipole_debye: Optional[float] = None   # avg |Δμ_i|
    max_delta_dipole_debye: Optional[float] = None  # max |Δμ_i|
    coordinates_file: Optional[str] = None       # GRO file used
    computed: bool = False                       # Whether computation succeeded
    # v4: Gating reasons for heuristic dipole check
    skip_reasons: List[str] = field(default_factory=list)  # e.g., ["giant_molecule", "percolated"]


@dataclass
class CorrectionResult:
    """Result of charge correction attempt."""
    original_q: float
    corrected: bool
    correction_method: Optional[str] = None
    effective_method: Optional[str] = None
    method_desc: Optional[str] = None
    correction_delta: Optional[float] = None
    patched_itp_path: Optional[Path] = None
    patched_atoms: Optional[List[Tuple[int, float, float]]] = None  # (idx, old_q, new_q)
    # Enhanced fields for physics validation
    max_delta_applied: Optional[float] = None    # Max |Δq| per atom actually applied
    total_delta: Optional[float] = None          # Total charge correction applied
    atoms_adjusted: int = 0                      # Number of atoms adjusted
    dipole_info: Optional[DipoleResult] = None   # Dipole shift info if computed
    target_molecule: Optional[str] = None        # Molecule type that was adjusted
    target_itp_moleculetype_name: Optional[str] = None  # Actual [ moleculetype ] name in patched ITP
    # New fields for enhanced manifest logging (Issue A/B)
    patched_atom_details: Optional[List[Dict]] = None  # [{nr, atomname, atomtype, line_no, charge_col_idx, old_q, new_q}, ...]
    patched_atom_details_truncated: int = 0
    correction_refused: bool = False             # True if correction was refused
    refuse_reason: Optional[str] = None          # Why correction was refused
    safe_subset_used: bool = False               # True if safe-subset filtering was applied
    unparsed_atoms_lines: Optional[List[Dict[str, Any]]] = None  # Skipped [ atoms ] lines if explicit override enabled
    unparsed_atoms_lines_truncated: int = 0
    # Audit fields (Task H) - run-scoped reproducibility
    run_id: Optional[str] = None                 # Run ID for output naming
    thresholds_used: Optional[Dict] = None       # warn/correct/neutrality thresholds
    neutrality_status: str = "unknown"           # "neutral", "corrected", "failed", "unknown"
    polymer_policy: Optional[Dict] = None        # Protected-polymer decision audit data
    # Rounding-only charge drift diagnostics
    rounding_only_passed: Optional[bool] = None
    rounding_moleculetype_tol: Optional[float] = None
    per_moleculetype_report: Optional[List[Dict[str, Any]]] = None
    dominant_charge_contributors: Optional[List[Dict[str, Any]]] = None
    # Atom chemistry audit for applied patches
    hetero_atoms_modified: bool = False
    hetero_atoms_modified_count: int = 0
    hetero_atoms_modified_details: Optional[List[Dict[str, Any]]] = None


# === Default allowlists for safe correction targets ===
# Common electrolyte solvents that are safe to modify for tiny charge corrections
DEFAULT_CORRECTION_TARGET_ALLOWLIST = frozenset({
    "PC", "EC", "DMC", "DEC", "EMC",  # Carbonate solvents
    "DOL", "DME",                      # Ether solvents
    "ACN", "AN",                       # Acetonitrile
    "DMSO",                            # DMSO
    "THP", "THF",                      # Cyclic ethers
})

# Conservative atom-name matching for safe correction targets.
# The intent is to avoid false positives (e.g., CL/CA/CU) even if that
# produces false negatives. Users can override via CLI if needed.
ALLOWED_H_PREFIXES = frozenset({
    # Exact H or common hydrogen classes (OPLS/CHARMM/GAFF style)
    "H", "HA", "HC", "HP", "HS", "HW", "HO", "HN",
    # "HGA" style hydrogen classes in some FFs
    "HG", "HGA", "HGB", "HGC", "HGD",
})

HETERO_ELEMENTS = frozenset({"O", "N", "S", "B", "F", "P"})


def _normalize_molname(name: str) -> str:
    """Normalize molecule name for case-insensitive matching."""
    return name.strip().casefold()

ALLOWED_C_PREFIXES = frozenset({
    # Common carbon classes (no bare "C" to avoid prefix ambiguity)
    "CT", "CB", "CC", "CD", "CE", "CG", "CH", "CI", "CJ", "CK",
    "CM", "CN", "CP", "CQ", "CR", "CS", "CX", "CY", "CZ",
})

# Explicit exclusions to avoid element/ion/halogen misclassification.
EXCLUDED_ATOM_NAMES = frozenset({
    "CL", "BR", "I", "F",
    "NA", "LI", "K", "MG", "CA", "ZN", "AL", "CU", "FE",
    "CO", "NI", "MN", "HG", "PB", "CD", "AG", "AU", "SN", "SB",
    "BA", "SR", "CS",
})

_H_NUMERIC_RE = re.compile(r"^H\d+$")
_C_NUMERIC_RE = re.compile(r"^C\d+$")
_EXCLUDED_RE = re.compile(
    r"^(?:"
    + "|".join(sorted(EXCLUDED_ATOM_NAMES))
    + r")(?:\d+)?$"
)


def _normalize_atomname(name: str) -> str:
    """Normalize atom name for matching (strip, uppercase)."""
    return name.strip().upper()


def _is_excluded_atomname(name: str) -> bool:
    """Return True if atom name is an excluded element/ion/halogen."""
    if not name:
        return False
    return bool(_EXCLUDED_RE.match(name))


def _matches_default_safe_atomname(name: str) -> bool:
    """
    Conservative matcher for safe atom names.
    - H: H, H1/H2..., or allowed H-class tokens (HC/HA/HP/HG/HGA...)
    - C: C1/C2..., or allowed carbon class tokens (CT/CB/CG...)
    - Never match excluded element-like tokens (CL/CA/CU/...)
    """
    if not name or _is_excluded_atomname(name):
        return False
    if name == "H" or _H_NUMERIC_RE.match(name):
        return True
    if name in ALLOWED_H_PREFIXES:
        return True
    # Allow H-class prefixes followed by digits (e.g., HGA1)
    for prefix in ALLOWED_H_PREFIXES:
        if prefix != "H" and name.startswith(prefix):
            suffix = name[len(prefix):]
            if suffix and suffix.isdigit():
                return True
    if _C_NUMERIC_RE.match(name):
        return True
    if name in ALLOWED_C_PREFIXES:
        return True
    return False


def _matches_custom_allowlist(name: str, allowlist: Set[str]) -> bool:
    """
    Backward-compatible allowlist matching:
    - Exact match wins
    - Prefix match allowed (legacy behavior), after exclusions are applied
    """
    for token in allowlist:
        if not token:
            continue
        if name == token or name.startswith(token):
            return True
    return False

# Atomtypes that should NEVER be adjusted (strongly charged heteroatoms)
# Keep display tokens in original form and use casefold for comparisons.
DEFAULT_DISALLOWED_ATOMTYPES = (
    "O", "O2", "OS", "OW", "OH", "OC",  # Oxygens (carbonyl, ether, water)
    "N", "N2", "N3", "NA", "NB", "NC",  # Nitrogens
    "S", "SH", "SS",                     # Sulfurs
    "F", "Cl", "Br", "I",               # Halogens
    "Li", "Na", "K", "Mg", "Ca", "Zn",  # Metal ions
)
DEFAULT_DISALLOWED_ATOMTYPES_CASEFOLD = frozenset(
    token.casefold() for token in DEFAULT_DISALLOWED_ATOMTYPES
)


class ChargeParser:
    """
    Parser for extracting charges from ITP [ atoms ] sections.
    
    GROMACS atoms section format:
    ; nr  type  resnr  residue  atom  cgnr  charge  mass
       1  ca       1     MOL     C1     1   0.1234  12.01
    
    Supports ITPs with multiple [ moleculetype ] blocks.
    Enhanced to detect preprocessor directives and validate field types.
    """
    
    # Task B: Updated regex to allow trailing comments after section header
    # Matches: [ atoms ], [ atoms ]  , [ atoms ] ; comment, [ atoms ]; comment
    SECTION_RE = re.compile(r'^\s*\[\s*(\w+)\s*\](?:\s*;.*)?\s*$')
    PREPROCESSOR_RE = re.compile(r'^\s*#\s*(ifdef|ifndef|else|endif|include|define)\b', re.IGNORECASE)
    NUMERIC_LITERAL_RE = re.compile(
        r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?"
    )

    @classmethod
    def _tokenize_atoms_line(cls, data_part: str) -> List[str]:
        """
        Tokenize an [ atoms ] data line robustly.

        Splits by whitespace, then explodes tokens containing glued numeric literals,
        e.g. ``1.000-0.001`` -> ``["1.000", "-0.001"]``.
        """
        tokens: List[str] = []
        for token in data_part.split():
            numeric_chunks = cls.NUMERIC_LITERAL_RE.findall(token)
            if len(numeric_chunks) >= 2 and "".join(numeric_chunks) == token:
                tokens.extend(numeric_chunks)
            else:
                tokens.append(token)
        return tokens

    @staticmethod
    def _split_data_comment(raw_line: str) -> Tuple[str, str]:
        """
        Split a line into data and comment fragments.

        Returns:
            (data_part, comment_part) where comment_part includes leading ';' if present.
        """
        stripped = raw_line.strip()
        if ";" in stripped:
            data_part, comment_part = stripped.split(";", 1)
            return data_part.strip(), ";" + comment_part
        return stripped, ""

    @classmethod
    def _is_atoms_data_line(cls, raw_line: str) -> bool:
        """
        True for non-empty [ atoms ] data lines (excluding comments/preprocessor).
        """
        stripped = raw_line.strip()
        if not stripped:
            return False
        if stripped.startswith(";"):
            return False
        if cls.PREPROCESSOR_RE.match(raw_line):
            return False
        data_part, _ = cls._split_data_comment(raw_line)
        return bool(data_part)

    @staticmethod
    def _line_excerpt(raw_line: str, max_len: int = 120) -> str:
        """Compact single-line preview for diagnostics/audit."""
        one_line = " ".join(raw_line.strip().split())
        if len(one_line) <= max_len:
            return one_line
        return one_line[: max_len - 3] + "..."

    @classmethod
    def _resolve_requested_moleculetype_name(
        cls,
        per_mol: Dict[str, Tuple[List[float], int, Optional[List[Dict]], bool]],
        requested_name: str,
        itp_path: Path,
    ) -> Tuple[str, str]:
        """
        Resolve request to real [ moleculetype ] name.

        Returns:
            (alias_name, real_itp_name)
        """
        available = list(per_mol.keys())
        if not available:
            return (requested_name or itp_path.stem), (requested_name or itp_path.stem)

        alias_name = requested_name or itp_path.stem
        if requested_name and requested_name in per_mol:
            return alias_name, requested_name

        req_key = _normalize_molname(requested_name) if requested_name else ""
        if requested_name:
            cf_matches = [name for name in available if _normalize_molname(name) == req_key]
            if len(cf_matches) == 1:
                return alias_name, cf_matches[0]
            if len(cf_matches) > 1:
                raise ChargeNeutralityError(
                    f"Ambiguous case-insensitive moleculetype match for '{requested_name}' in {itp_path}: "
                    f"{cf_matches}"
                )

        if len(available) == 1:
            return alias_name, available[0]

        raise ChargeNeutralityError(
            f"Molecule '{requested_name}' not found in {itp_path}.\n"
            f"Available moleculetypes: {available}\n"
            "Specify the exact moleculetype name."
        )

    @classmethod
    def _collect_unparsed_atoms_lines(
        cls,
        itp_path: Path,
        target_moleculetype_name: str,
    ) -> List[Dict[str, Any]]:
        """
        Return unparseable [ atoms ] data lines for a specific moleculetype.
        """
        if not itp_path.exists():
            return []

        try:
            content = itp_path.read_text(encoding="utf-8")
        except Exception as exc:
            raise ChargeNeutralityError(f"Cannot read ITP file {itp_path}: {exc}")

        target_key = _normalize_molname(target_moleculetype_name)
        current_section = ""
        current_moleculetype: Optional[str] = None
        in_target_atoms = False
        found_target = False
        unparsed: List[Dict[str, Any]] = []

        for line_no, line in enumerate(content.splitlines(), 1):
            section_match = cls.SECTION_RE.match(line)
            if section_match:
                section = section_match.group(1).lower()
                current_section = section
                if section == "moleculetype":
                    current_moleculetype = None
                    in_target_atoms = False
                elif section == "atoms":
                    in_target_atoms = (
                        current_moleculetype is not None
                        and _normalize_molname(current_moleculetype) == target_key
                    )
                else:
                    in_target_atoms = False
                continue

            stripped = line.strip()
            if current_section == "moleculetype":
                if not stripped or stripped.startswith(";") or cls.PREPROCESSOR_RE.match(line):
                    continue
                data_part, _ = cls._split_data_comment(line)
                if not data_part:
                    continue
                parts = data_part.split()
                if parts:
                    current_moleculetype = parts[0]
                    if _normalize_molname(current_moleculetype) == target_key:
                        found_target = True
                continue

            if not in_target_atoms:
                continue
            if not cls._is_atoms_data_line(line):
                continue

            try:
                parsed = cls._parse_atoms_line(line, line_no, strict=True)
            except ChargeNeutralityError as exc:
                parsed = None
                reason = str(exc).splitlines()[0]
            else:
                reason = "Unparseable [ atoms ] data line"

            if parsed is not None:
                continue
            unparsed.append(
                {
                    "line_no": line_no,
                    "excerpt": cls._line_excerpt(line),
                    "reason": reason,
                }
            )

        if not found_target:
            raise ChargeNeutralityError(
                f"Moleculetype '{target_moleculetype_name}' not found in {itp_path} while validating [ atoms ] lines."
            )

        return unparsed
    
    @classmethod
    def _parse_atoms_line(
        cls,
        line: str,
        line_no: int,
        strict: bool = False,
    ) -> Optional[Dict]:
        """
        Parse a single [ atoms ] line with robust column detection (Task C).
        
        GROMACS [ atoms ] format variants:
        - Standard (8 tokens): nr type resnr residue atom cgnr charge mass
        - No cgnr (7 tokens):  nr type resnr residue atom charge mass
        
        Detection strategy:
        1. Strip inline comments
        2. Tokenize
        3. Identify format by token count and validate:
           - nr (token 0) must be int
           - resnr (token 2) must be int
           - Last token should parse as float (mass)
           - Second-to-last should parse as float (charge)
        
        Args:
            line: Raw line from ITP file
            line_no: Line number (1-indexed) for error messages
            strict: If True, raise on ambiguous format; else return None
            
        Returns:
            Dict with {nr, atomtype, resnr, residue, atomname, cgnr, charge, mass, 
                      charge_col_idx, line_no, original_line} or None if invalid/ambiguous
            
        Raises:
            ChargeNeutralityError: In strict mode when format is ambiguous
        """
        # Strip inline comments
        stripped = line.strip()
        if not stripped or stripped.startswith(';'):
            return None
        
        if ';' in stripped:
            data_part, comment_part = stripped.split(';', 1)
            data_part = data_part.strip()
            comment_part = ";" + comment_part
        else:
            data_part = stripped
            comment_part = ""
        
        if not data_part:
            return None
        
        tokens = cls._tokenize_atoms_line(data_part)
        n_tokens = len(tokens)
        
        # Need at least 7 tokens for minimal valid format
        if n_tokens < 7:
            if strict:
                raise ChargeNeutralityError(
                    f"Line {line_no}: Too few tokens ({n_tokens}) in [ atoms ] line.\n"
                    f"Expected at least 7 (nr type resnr residue atom charge mass).\n"
                    f"Line: {line.rstrip()}"
                )
            return None
        
        def _is_int(idx: int) -> bool:
            if idx >= n_tokens:
                return False
            try:
                int(tokens[idx])
                return True
            except ValueError:
                return False

        def _is_float(idx: int) -> bool:
            if idx >= n_tokens:
                return False
            try:
                float(tokens[idx])
                return True
            except ValueError:
                return False

        has_core_ints = _is_int(0) and _is_int(2)
        looks_like_8 = (
            has_core_ints and
            n_tokens >= 8 and
            _is_int(5) and
            _is_float(6) and
            _is_float(7)
        )
        looks_like_7 = (
            has_core_ints and
            n_tokens >= 7 and
            _is_float(5) and
            _is_float(6)
        )

        if looks_like_8:
            charge_col_idx = 6
            mass_col_idx = 7
            cgnr = int(tokens[5])
        elif looks_like_7:
            charge_col_idx = 5
            mass_col_idx = 6
            cgnr = None
        else:
            msg = (
                f"Line {line_no}: Ambiguous [ atoms ] format; could not classify as "
                f"8-token (nr type resnr residue atom cgnr charge mass) or 7-token "
                f"(nr type resnr residue atom charge mass).\n"
                f"Line: {line.rstrip()}\n"
                f"Hint: check for glued numeric fields like '1.000-0.001' or malformed spacing."
            )
            if strict:
                raise ChargeNeutralityError(msg)
            print(f"    [WARN] {msg}")
            return None

        # Must exist by format checks above
        nr = int(tokens[0])
        resnr = int(tokens[2])
        atomtype = tokens[1]
        residue = tokens[3]
        atomname = tokens[4]
        charge = float(tokens[charge_col_idx])

        # Mass is optional in some odd formats; keep None if unparseable.
        try:
            mass = float(tokens[mass_col_idx]) if mass_col_idx < len(tokens) else None
        except ValueError:
            mass = None
        
        return {
            'nr': nr,
            'atomtype': atomtype,
            'resnr': resnr,
            'residue': residue,
            'atomname': atomname,
            'cgnr': cgnr,
            'charge': charge,
            'mass': mass,
            'charge_col_idx': charge_col_idx,
            'line_no': line_no,
            'original_line': line,
            'tokens': tokens,
            'comment': comment_part,
        }
    
    @classmethod
    def parse_atoms_charges_per_molecule(
        cls,
        itp_path: Path,
        include_details: bool = False,
        strict: bool = False,
    ) -> Dict[str, Tuple[List[float], int, Optional[List[Dict]], bool]]:
        """
        Parse charges per moleculetype from ITP file.
        
        Supports files with multiple [ moleculetype ] blocks.
        Enhanced with preprocessor detection and atom details.
        
        Args:
            itp_path: Path to ITP file
            include_details: If True, include atom details (nr, atomname, atomtype)
            strict: If True, raise on ambiguous atom-line parsing
            
        Returns:
            Dict of {mol_name: (list_of_charges, atom_count, atom_details, has_preprocessor)}
            atom_details is None if include_details=False
        """
        if not itp_path.exists():
            return {}
        
        try:
            content = itp_path.read_text(encoding="utf-8")
        except Exception as e:
            raise ChargeNeutralityError(f"Cannot read ITP file {itp_path}: {e}")
        
        lines = content.splitlines()
        
        result: Dict[str, Tuple[List[float], int, Optional[List[Dict]], bool]] = {}
        current_mol_name: Optional[str] = None
        current_charges: List[float] = []
        current_details: List[Dict] = []
        current_has_preprocessor: bool = False
        current_section = ""
        
        for line_no, line in enumerate(lines, 1):
            # Check for preprocessor directives anywhere
            if cls.PREPROCESSOR_RE.match(line):
                if current_section == 'atoms':
                    current_has_preprocessor = True
                continue
            
            # Check for section header
            match = cls.SECTION_RE.match(line)
            if match:
                section = match.group(1).lower()
                # If we're entering a new moleculetype, save the previous one
                if section == 'moleculetype' and current_mol_name:
                    details = current_details if include_details else None
                    result[current_mol_name] = (
                        current_charges.copy(),
                        len(current_charges),
                        details,
                        current_has_preprocessor,
                    )
                    current_charges = []
                    current_details = []
                    current_has_preprocessor = False
                    current_mol_name = None
                current_section = section
                continue
            
            # Parse moleculetype name
            if current_section == 'moleculetype':
                stripped = line.strip()
                if stripped and not stripped.startswith(';'):
                    data = stripped.split(';')[0].strip()
                    if data:
                        parts = data.split()
                        if parts:
                            current_mol_name = parts[0]
            
            # Parse atoms section using robust parser (Task C)
            if current_section == 'atoms' and current_mol_name:
                # Only non-empty/comment-stripped, non-preprocessor lines are atom data lines.
                if not cls._is_atoms_data_line(line):
                    continue
                parsed = cls._parse_atoms_line(line, line_no, strict=strict)
                if parsed is None:
                    continue
                current_charges.append(parsed['charge'])
                
                if include_details:
                    current_details.append({
                        'nr': parsed['nr'],
                        'atomname': parsed['atomname'],
                        'atomtype': parsed['atomtype'],
                        'charge': parsed['charge'],
                        'charge_col_idx': parsed['charge_col_idx'],
                        'cgnr': parsed['cgnr'],
                        'line_no': line_no,
                        'line_excerpt': cls._line_excerpt(parsed.get('original_line', line)),
                    })
        
        # Don't forget the last molecule
        if current_mol_name:
            details = current_details if include_details else None
            result[current_mol_name] = (
                current_charges.copy(),
                len(current_charges),
                details,
                current_has_preprocessor,
            )
        
        return result

    @classmethod
    def parse_molecule_bond_graph(
        cls,
        itp_path: Path,
        mol_name: str,
    ) -> Dict[int, Set[int]]:
        """
        Parse [ bonds ] for a specific moleculetype and return adjacency by atom nr.

        Returns:
            Dict[int, Set[int]] mapping 1-based atom numbers to bonded neighbors.
        """
        if not itp_path.exists():
            return {}
        try:
            content = itp_path.read_text(encoding="utf-8")
        except Exception:
            return {}

        current_section = ""
        current_mol: Optional[str] = None
        in_target = False
        graph: Dict[int, Set[int]] = {}

        for raw in content.splitlines():
            section_match = cls.SECTION_RE.match(raw)
            if section_match:
                current_section = section_match.group(1).lower()
                if current_section == "moleculetype":
                    current_mol = None
                    in_target = False
                continue

            stripped = raw.strip()
            if not stripped or stripped.startswith(";"):
                continue
            if cls.PREPROCESSOR_RE.match(raw):
                continue

            data = stripped.split(";", 1)[0].strip()
            if not data:
                continue

            if current_section == "moleculetype":
                parts = data.split()
                if parts:
                    current_mol = parts[0]
                    in_target = _normalize_molname(current_mol) == _normalize_molname(mol_name)
                continue

            if not in_target or current_section != "bonds":
                continue

            parts = data.split()
            if len(parts) < 2:
                continue
            if not (parts[0].isdigit() and parts[1].isdigit()):
                continue
            ai = int(parts[0])
            aj = int(parts[1])
            graph.setdefault(ai, set()).add(aj)
            graph.setdefault(aj, set()).add(ai)

        return graph
    
    @classmethod
    def parse_atoms_charges(cls, itp_path: Path) -> Tuple[List[float], int]:
        """
        Parse charges from [ atoms ] section (legacy single-molecule API).
        
        For multi-moleculetype ITPs, returns charges from first moleculetype only.
        Use parse_atoms_charges_per_molecule() for full support.
        
        Args:
            itp_path: Path to ITP file
            
        Returns:
            Tuple of (list of charges, atom count)
        """
        per_mol = cls.parse_atoms_charges_per_molecule(itp_path)
        if not per_mol:
            return [], 0
        
        # Return first molecule's charges for backward compatibility
        first_mol = list(per_mol.values())[0]
        return first_mol[0], first_mol[1]
    
    @classmethod
    def get_molecule_charge(
        cls,
        itp_path: Path,
        name: str = "",
        include_details: bool = False,
        strict: bool = False,
        expected_itp_moleculetype_name: Optional[str] = None,
        allow_unparsed_atoms_lines: bool = False,
    ) -> MoleculeCharge:
        """
        Get total charge for a molecule from its ITP file.
        
        FAIL-FAST: Refuses ambiguous/unknown moleculetype resolution and, by default,
        refuses unparseable [ atoms ] data lines in the resolved target moleculetype.
        
        Args:
            itp_path: Path to ITP file
            name: Pipeline alias / requested molecule name
            include_details: If True, include atom details
            strict: If True, enforce strict [ atoms ] parsing
            expected_itp_moleculetype_name: Optional explicit [ moleculetype ] name to resolve
            allow_unparsed_atoms_lines: Explicit unsafe override to tolerate unparsed data lines
            
        Returns:
            MoleculeCharge object
            
        Raises:
            ChargeNeutralityError: On unresolved moleculetype or unparseable [ atoms ] data lines
        """
        per_mol = cls.parse_atoms_charges_per_molecule(
            itp_path,
            include_details=include_details,
            strict=strict,
        )

        if not per_mol:
            alias_name = name or itp_path.stem
            itp_name = expected_itp_moleculetype_name or alias_name
            charges, atom_count, details, has_preproc = [], 0, None, False
            unparsed_atoms_lines: List[Dict[str, Any]] = []
        else:
            requested = expected_itp_moleculetype_name if expected_itp_moleculetype_name else name
            alias_name, itp_name = cls._resolve_requested_moleculetype_name(
                per_mol=per_mol,
                requested_name=requested,
                itp_path=itp_path,
            )
            if name:
                alias_name = name
            charges, atom_count, details, has_preproc = per_mol[itp_name]

            unparsed_atoms_lines = cls._collect_unparsed_atoms_lines(
                itp_path=itp_path,
                target_moleculetype_name=itp_name,
            )
            if unparsed_atoms_lines and not allow_unparsed_atoms_lines:
                first = unparsed_atoms_lines[0]
                raise ChargeNeutralityError(
                    "Unparseable [ atoms ] data line detected in target moleculetype "
                    f"'{itp_name}' ({itp_path}).\n"
                    f"Line {first.get('line_no')}: {first.get('excerpt')}\n"
                    "Refusing to continue to avoid silent charge under-count."
                )
        
        # Use math.fsum for numerically stable summation (Issue C)
        total = math.fsum(charges)
        
        return MoleculeCharge(
            name=alias_name,
            itp_moleculetype_name=itp_name,
            itp_path=itp_path,
            total_charge=total,
            atom_count=atom_count,
            charges=charges,
            atom_details=details,
            has_preprocessor=has_preproc,
            unparsed_atoms_lines=unparsed_atoms_lines,
        )


class ChargeNeutralityChecker:
    """
    Checks and corrects charge neutrality for a molecular system.
    
    Enhanced with physics validation (v3):
    - Safe-subset atom selection (not uniform smearing)
    - Target allowlist enforcement
    - Numerically stable summation (math.fsum)
    - Preprocessor detection and rejection
    - No silent fallbacks
    """
    
    # Threshold for warnings (Issue E: more forgiving for FF rounding)
    WARN_THRESHOLD = 1e-4
    
    # Threshold for auto-correction (Issue E: configurable)
    CORRECTION_THRESHOLD = 1e-4
    
    # Default tolerance for neutrality pass (Task I)
    DEFAULT_NEUTRALITY_TOL = 1e-6
    
    # Physics guardrails (can be overridden via constructor)
    DEFAULT_MAX_DELTA_PER_ATOM = 1e-5  # Slightly relaxed for safe-subset
    DEFAULT_MAX_DIPOLE_SHIFT_DEBYE = 0.01
    DEFAULT_MAX_ABS_Q_FOR_ADJUST = 0.15  # Only adjust atoms with |q| < this
    DEFAULT_POLYMER_NET_CHARGE_TOL = 1e-3
    DEFAULT_POLYMER_EXCLUSION_BONDS = 2
    DEFAULT_SPREAD_SAFE_MIN_ATOMS = 8
    DEFAULT_MOLECULETYPE_ROUNDING_TOL = 1e-4
    MAX_PATCHED_ATOM_DETAILS = 200
    MAX_UNPARSED_ATOMS_AUDIT = 50
    
    def __init__(
        self,
        threshold: float = CORRECTION_THRESHOLD,
        max_delta_per_atom: float = DEFAULT_MAX_DELTA_PER_ATOM,
        max_total_delta: Optional[float] = None,  # None = use threshold
        max_dipole_shift_debye: float = DEFAULT_MAX_DIPOLE_SHIFT_DEBYE,
        strict_charge_physics: bool = False,
        # New config options (Issue A)
        correction_method: str = "safe_subset",  # "safe_subset" or "uniform_all"
        allowed_atomnames: Optional[Set[str]] = None,
        disallow_atomtypes: Optional[Set[str]] = None,
        max_abs_q_for_adjust: float = DEFAULT_MAX_ABS_Q_FOR_ADJUST,
        # New config options (Issue B)
        target_allowlist: Optional[Set[str]] = None,
        require_explicit_target: bool = False,
        protected_polymers: Optional[Set[str]] = None,
        polymer_net_charge_tol: float = DEFAULT_POLYMER_NET_CHARGE_TOL,
        polymer_correction_method: str = "skip_if_small",  # skip_if_small | spread_safe
        polymer_exclusion_bonds: int = DEFAULT_POLYMER_EXCLUSION_BONDS,
        # Task E,F,G: New strict mode and ion handling options
        strict_mode: bool = False,  # Fail-fast on preprocessor, ambiguous parsing, count mismatches
        allow_ions: bool = False,  # Allow correction targeting detected ions (advisory override)
        warn_threshold: Optional[float] = None,  # Override WARN_THRESHOLD
        neutrality_tol: float = DEFAULT_NEUTRALITY_TOL,  # Tolerance for |Q| < tol → neutral
        rounding_moleculetype_tol: float = DEFAULT_MOLECULETYPE_ROUNDING_TOL,
        allow_non_rounding_correction: bool = False,
        allow_unparsed_atoms_lines: bool = False,
    ):
        """
        Initialize checker.
        
        Args:
            threshold: Maximum |Q| for auto-correction
            max_delta_per_atom: Maximum |Δq| allowed per atom
            max_total_delta: Maximum |ΔQ| total correction (defaults to threshold)
            max_dipole_shift_debye: Maximum allowed dipole shift in Debye
            strict_charge_physics: If True, dipole violations are errors; else warnings
            correction_method: "safe_subset" (default) or "uniform_all"
            allowed_atomnames: Atom names allowed for adjustment (None = use default)
            disallow_atomtypes: Atom types never adjusted (None = use default)
            max_abs_q_for_adjust: Max |q| for atom to be adjustable
            target_allowlist: Molecule names allowed as correction targets
            require_explicit_target: If True, must specify target_molecule explicitly
            protected_polymers: Molecule names classified as protected_polymer
            polymer_net_charge_tol: Acceptable per-molecule polymer charge drift (|q| <= tol)
            polymer_correction_method: "skip_if_small" or "spread_safe"
            polymer_exclusion_bonds: Exclude atoms within N bonds of hetero atoms for spread_safe
            strict_mode: If True, fail on preprocessor directives, ambiguous parsing
            allow_ions: If True, allow correction on molecules detected as ions
            warn_threshold: Threshold for charge warning (defaults to WARN_THRESHOLD)
            neutrality_tol: Tolerance for considering system neutral (|Q| < tol)
            rounding_moleculetype_tol: Max |Q_mol - round(Q_mol)| treated as rounding drift
            allow_non_rounding_correction: Explicit override to allow non-rounding charge fixes
            allow_unparsed_atoms_lines: Explicit unsafe override for unparseable [ atoms ] data lines
        """
        if isinstance(threshold, bool):
            raise ChargeNeutralityError(
                f"Charge neutrality threshold must be a float, not bool (got {threshold!r})."
            )
        self.threshold = threshold
        self.max_delta_per_atom = max_delta_per_atom
        self.max_total_delta = max_total_delta if max_total_delta is not None else threshold
        self.max_dipole_shift_debye = max_dipole_shift_debye
        self.strict_charge_physics = strict_charge_physics
        
        # Task E: Strict mode flag
        self.strict_mode = strict_mode
        
        # Task G: Ion handling
        self.allow_ions = allow_ions
        
        # Task I: Thresholds
        self.warn_threshold = warn_threshold if warn_threshold is not None else self.WARN_THRESHOLD
        self.neutrality_tol = neutrality_tol
        self.rounding_moleculetype_tol = max(0.0, rounding_moleculetype_tol)
        self.allow_non_rounding_correction = allow_non_rounding_correction
        self.allow_unparsed_atoms_lines = bool(allow_unparsed_atoms_lines)
        
        # Issue A: Safe-subset config
        self.correction_method = correction_method
        # If user provides allowlist, normalize and use legacy prefix matching.
        # Otherwise use conservative default matcher based on H/C class tokens.
        if allowed_atomnames is not None:
            self.allowed_atomnames = {
                _normalize_atomname(s) for s in allowed_atomnames if s and s.strip()
            }
        else:
            self.allowed_atomnames = None
        self.allowed_h_prefixes = ALLOWED_H_PREFIXES
        self.allowed_c_prefixes = ALLOWED_C_PREFIXES
        self.excluded_atom_names = EXCLUDED_ATOM_NAMES
        raw_disallow_tokens = (
            disallow_atomtypes if disallow_atomtypes is not None
            else DEFAULT_DISALLOWED_ATOMTYPES
        )
        self.disallow_atomtypes_tokens = tuple(
            sorted(
                {
                    token.strip()
                    for token in raw_disallow_tokens
                    if token and token.strip()
                },
                key=lambda s: s.casefold(),
            )
        )
        if disallow_atomtypes is None:
            self.disallow_atomtypes = set(DEFAULT_DISALLOWED_ATOMTYPES_CASEFOLD)
        else:
            self.disallow_atomtypes = {
                token.casefold() for token in self.disallow_atomtypes_tokens
            }
        self.max_abs_q_for_adjust = max_abs_q_for_adjust
        
        # Issue B: Target allowlist
        raw_allowlist = (
            target_allowlist if target_allowlist is not None else DEFAULT_CORRECTION_TARGET_ALLOWLIST
        )
        self.target_allowlist = {
            _normalize_molname(token) for token in raw_allowlist if token and token.strip()
        }
        self.require_explicit_target = require_explicit_target
        self.protected_polymers = {
            _normalize_molname(token) for token in (protected_polymers or set()) if token and token.strip()
        }
        self.polymer_net_charge_tol = polymer_net_charge_tol
        if polymer_correction_method not in {"skip_if_small", "spread_safe"}:
            raise ChargeNeutralityError(
                f"Invalid polymer_correction_method '{polymer_correction_method}'. "
                "Use 'skip_if_small' or 'spread_safe'."
            )
        self.polymer_correction_method = polymer_correction_method
        self.polymer_exclusion_bonds = max(0, polymer_exclusion_bonds)
        self.polymer_spread_min_atoms = self.DEFAULT_SPREAD_SAFE_MIN_ATOMS
        
        self.molecule_charges: Dict[str, MoleculeCharge] = {}
        self.unparsed_atoms_audit: List[Dict[str, Any]] = []
    
    def load_molecule_charges(
        self,
        molecule_itp_paths: Dict[str, Path],
    ) -> Dict[str, MoleculeCharge]:
        """
        Load charge information for all molecules (with atom details).
        
        Args:
            molecule_itp_paths: Dict of molecule name to ITP path
            
        Returns:
            Dict of molecule name to MoleculeCharge
        """
        self.molecule_charges = {}
        self.unparsed_atoms_audit = []
        for name, itp_path in molecule_itp_paths.items():
            self.molecule_charges[name] = ChargeParser.get_molecule_charge(
                itp_path,
                name,
                include_details=True,
                strict=self.strict_mode,
                allow_unparsed_atoms_lines=self.allow_unparsed_atoms_lines,
            )
            mol_unparsed = self.molecule_charges[name].unparsed_atoms_lines or []
            if mol_unparsed:
                print(
                    f"    [WARN] allow_unparsed_atoms_lines=True: "
                    f"{len(mol_unparsed)} unparseable [ atoms ] data line(s) in "
                    f"{itp_path} ({self.molecule_charges[name].itp_moleculetype_name})."
                )
                for item in mol_unparsed:
                    self.unparsed_atoms_audit.append(
                        {
                            "molecule_alias": name,
                            "itp_moleculetype_name": self.molecule_charges[name].itp_moleculetype_name,
                            "itp_path": str(itp_path),
                            **item,
                        }
                    )
        
        return self.molecule_charges

    def _thresholds_used(self) -> Dict[str, Any]:
        """Return the effective thresholds/limits used for this run."""
        return {
            "warn_threshold": self.warn_threshold,
            "correction_threshold": self.threshold,
            "max_delta_per_atom": self.max_delta_per_atom,
            "max_total_delta": self.max_total_delta,
            "neutrality_tol": self.neutrality_tol,
            "max_dipole_shift_debye": self.max_dipole_shift_debye,
            "max_abs_q_for_adjust": self.max_abs_q_for_adjust,
            "polymer_net_charge_tol": self.polymer_net_charge_tol,
            "polymer_correction_method": self.polymer_correction_method,
            "polymer_exclusion_bonds": self.polymer_exclusion_bonds,
            "rounding_moleculetype_tol": self.rounding_moleculetype_tol,
            "allow_non_rounding_correction": self.allow_non_rounding_correction,
            "allow_unparsed_atoms_lines": self.allow_unparsed_atoms_lines,
        }
    
    def compute_system_charge(
        self,
        molecule_counts: Dict[str, int],
        strict: bool = True,
    ) -> float:
        """
        Compute total system charge using numerically stable summation.
        
        Q = Σ(N_i * q_molecule_i)
        
        This uses float arithmetic with math.fsum. At the correction scales used
        here (1e-4 and below), this is deterministic and numerically sufficient,
        and avoids Decimal(str(float)) artifacts from freezing rounded float text.
        
        Args:
            molecule_counts: Dict of molecule name to count
            strict: If True, fail on unknown molecules with count > 0
            
        Returns:
            Total system charge
            
        Raises:
            ChargeNeutralityError: If strict=True and molecule missing charge definition
        """
        contributions: List[float] = []
        loaded_by_casefold: Dict[str, str] = {}
        for alias_name, mol in self.molecule_charges.items():
            for key in {
                _normalize_molname(alias_name),
                _normalize_molname(mol.itp_moleculetype_name),
            }:
                loaded_by_casefold.setdefault(key, alias_name)
        for name, count in molecule_counts.items():
            resolved = loaded_by_casefold.get(_normalize_molname(name))
            if resolved is not None:
                mol_q = self.molecule_charges[resolved].total_charge
                contributions.append(count * mol_q)
            elif count > 0 and strict:
                raise ChargeNeutralityError(
                    f"Missing charge definition for molecule '{name}' (count={count}).\n"
                    f"Loaded molecules: {list(self.molecule_charges.keys())}"
                )
        
        return math.fsum(contributions)
    
    @staticmethod
    def is_ion(mol_charge: MoleculeCharge, tol: float = 0.1) -> bool:
        """
        Check if a molecule appears to be an ion.
        
        Ions have integer or half-integer charges (±0.5, ±1, ±2, etc.)
        and should never be modified for charge smearing.
        
        Enhanced: also check atom count (ions typically 1-3 atoms) (Issue E).
        """
        q = mol_charge.total_charge
        # Check if close to integer or half-integer
        rounded = round(q * 2) / 2  # Round to nearest 0.5
        is_integer_charge = abs(q - rounded) < tol and abs(q) >= 0.5
        
        # Additional heuristic: small atom count (Issue E)
        is_small = mol_charge.atom_count <= 3
        
        return is_integer_charge or (is_small and abs(q) > 0.4)
    
    @staticmethod
    def is_neutral(mol_charge: MoleculeCharge, tol: float = 0.01) -> bool:
        """
        Check if a molecule is intended to be neutral (total charge ~0).
        
        Neutral molecules are safe targets for charge correction.
        """
        return abs(mol_charge.total_charge) < tol
    
    def _is_suspicious_target(self, mol_name: str, mol_charge: MoleculeCharge) -> List[str]:
        """
        Check if target looks like fragment/oligomer/unknown (Issue B).
        
        Returns list of warning messages (empty if OK).
        """
        warnings = []
        name_lower = mol_name.lower()
        
        # Fragment patterns
        fragment_patterns = ['frag', 'oligo', 'residue', 'chain', 'polymer']
        for pattern in fragment_patterns:
            if pattern in name_lower:
                warnings.append(f"Target '{mol_name}' contains '{pattern}' pattern")
        
        # Too many atoms for typical solvent
        if mol_charge.atom_count > 100:
            warnings.append(
                f"Target '{mol_name}' has {mol_charge.atom_count} atoms "
                "(unusually large for solvent)"
            )
        
        # Very few atoms but not an ion
        if mol_charge.atom_count < 5 and not self.is_ion(mol_charge):
            warnings.append(
                f"Target '{mol_name}' has only {mol_charge.atom_count} atoms "
                "(unusually small for neutral molecule)"
            )
        
        return warnings

    def _resolve_target_name(self, name: str) -> Optional[str]:
        """
        Resolve target to loaded alias key.

        Accepts either pipeline alias or real ITP [ moleculetype ] name (case-insensitive).
        """
        if name in self.molecule_charges:
            return name
        key = _normalize_molname(name)
        matches: List[str] = []
        for alias_name, mol in self.molecule_charges.items():
            if _normalize_molname(alias_name) == key:
                matches.append(alias_name)
                continue
            if _normalize_molname(mol.itp_moleculetype_name) == key:
                matches.append(alias_name)
        unique = sorted(set(matches), key=lambda s: s.casefold())
        if len(unique) == 1:
            return unique[0]
        return None

    def _is_allowlisted_target(self, mol_name: str) -> bool:
        keys = {_normalize_molname(mol_name)}
        mol = self.molecule_charges.get(mol_name)
        if mol is not None:
            keys.add(_normalize_molname(mol.itp_moleculetype_name))
        return any(key in self.target_allowlist for key in keys)

    def _is_protected_polymer(self, mol_name: str) -> bool:
        keys = {_normalize_molname(mol_name)}
        mol = self.molecule_charges.get(mol_name)
        if mol is not None:
            keys.add(_normalize_molname(mol.itp_moleculetype_name))
        return any(key in self.protected_polymers for key in keys)

    @staticmethod
    def _atom_primary_element(atomname: str, atomtype: str) -> str:
        """Best-effort element token from atomname/atomtype."""
        for token in (atomname, atomtype):
            if not token:
                continue
            match = re.match(r"^[A-Za-z]+", token.strip())
            if not match:
                continue
            letters = match.group(0).upper()
            if len(letters) >= 2 and letters[:2] in {"CL", "BR"}:
                return letters[:2]
            return letters[0]
        return ""

    def _polymer_policy_decision(
        self,
        q: float,
        molecule_counts: Dict[str, int],
        target_molecule: Optional[str],
    ) -> Dict[str, Optional[object]]:
        """
        Decide protected-polymer correction policy before attempting correction.

        Returns a dict with keys:
          - decision: one of {"no_polymer_target","acceptable_drift_skip","skip_non_strict","allow_spread_safe","strict_fail"}
          - target: target molecule name or None
          - net_q: per-molecule charge for target or None
          - tol: configured polymer tolerance
          - allowlisted: bool or None
        """
        decision = {
            "decision": "no_polymer_target",
            "target": None,
            "net_q": None,
            "tol": self.polymer_net_charge_tol,
            "allowlisted": None,
            "system_q": q,
        }

        if not self.protected_polymers:
            return decision

        count_by_casefold = {
            _normalize_molname(name): count for name, count in molecule_counts.items()
        }

        chosen_target: Optional[str] = None
        if target_molecule is not None:
            resolved = self._resolve_target_name(target_molecule)
            if resolved and self._is_protected_polymer(resolved):
                chosen_target = resolved
        else:
            contributors: List[str] = []
            for name, count in molecule_counts.items():
                loaded = self._resolve_target_name(name)
                if loaded is None or count <= 0 or not self._is_protected_polymer(loaded):
                    continue
                contrib = self.molecule_charges[loaded].total_charge * count
                if abs(contrib) > self.neutrality_tol:
                    contributors.append(loaded)
            unique = sorted(set(contributors), key=lambda s: s.casefold())
            if len(unique) == 1:
                candidate = unique[0]
                candidate_charge = self.molecule_charges[candidate]
                candidate_count = count_by_casefold.get(
                    _normalize_molname(candidate),
                    count_by_casefold.get(
                        _normalize_molname(candidate_charge.itp_moleculetype_name),
                        0,
                    ),
                )
                contrib = candidate_charge.total_charge * candidate_count
                if abs(q - contrib) <= max(self.neutrality_tol * 10.0, 1e-12):
                    chosen_target = candidate

        if chosen_target is None:
            return decision

        net_q = self.molecule_charges[chosen_target].total_charge
        allowlisted = self._is_allowlisted_target(chosen_target)
        decision.update(
            {
                "target": chosen_target,
                "net_q": net_q,
                "allowlisted": allowlisted,
            }
        )

        if abs(net_q) <= self.polymer_net_charge_tol:
            if self.polymer_correction_method == "spread_safe" and allowlisted:
                decision["decision"] = "allow_spread_safe"
            else:
                decision["decision"] = "acceptable_drift_skip"
            return decision

        if not allowlisted:
            decision["decision"] = "strict_fail" if self.strict_mode else "skip_non_strict"
            return decision

        decision["decision"] = "allow_spread_safe"
        return decision

    def _get_spread_safe_atom_indices(
        self,
        mol_charge: MoleculeCharge,
    ) -> Tuple[List[int], str]:
        """
        Select spread-safe polymer correction indices.

        Rules:
        - never include hetero atoms (O/N/S/B/F/P)
        - never include atoms within N bonds of hetero atoms
        - keep disallow atomtype / |q| guardrails
        """
        if not mol_charge.atom_details:
            return [], "No atom details available"

        bond_graph = ChargeParser.parse_molecule_bond_graph(
            mol_charge.itp_path,
            mol_charge.itp_moleculetype_name,
        )
        if not bond_graph and self.polymer_exclusion_bonds > 0:
            return [], (
                "No [ bonds ] graph found for polymer spread_safe selection. "
                "Cannot enforce hetero-neighborhood exclusion."
            )

        hetero_nrs: Set[int] = set()
        for idx, atom in enumerate(mol_charge.atom_details):
            nr = int(atom.get("nr", idx + 1))
            atomname = str(atom.get("atomname", ""))
            atomtype = str(atom.get("atomtype", ""))
            elem = self._atom_primary_element(atomname, atomtype)
            if elem in HETERO_ELEMENTS:
                hetero_nrs.add(nr)

        excluded_nrs = set(hetero_nrs)
        frontier = set(hetero_nrs)
        for _ in range(self.polymer_exclusion_bonds):
            if not frontier:
                break
            next_frontier: Set[int] = set()
            for nr in frontier:
                for nb in bond_graph.get(nr, set()):
                    if nb not in excluded_nrs:
                        excluded_nrs.add(nb)
                        next_frontier.add(nb)
            frontier = next_frontier

        safe_indices: List[int] = []
        for idx, atom in enumerate(mol_charge.atom_details):
            nr = int(atom.get("nr", idx + 1))
            if nr in excluded_nrs:
                continue

            atomname = str(atom.get("atomname", ""))
            atomtype = str(atom.get("atomtype", ""))
            charge = float(atom.get("charge", 0.0))
            elem = self._atom_primary_element(atomname, atomtype)

            if elem not in {"C", "H"}:
                continue
            if atomtype.casefold() in self.disallow_atomtypes:
                continue
            if abs(charge) > self.max_abs_q_for_adjust:
                continue

            safe_indices.append(idx)

        if len(safe_indices) < self.polymer_spread_min_atoms:
            return [], (
                f"spread_safe found only {len(safe_indices)} eligible atoms after excluding "
                f"hetero atoms and {self.polymer_exclusion_bonds}-bond neighborhoods. "
                "Add explicit polymer target allowlist and/or reduce "
                "--charge-fix-polymer-exclusion-bonds if physically justified."
            )

        return safe_indices, (
            f"spread_safe(n={len(safe_indices)}, "
            f"excluded_hetero_neighborhood_bonds={self.polymer_exclusion_bonds})"
        )
    
    def _get_safe_atom_indices(
        self,
        mol_charge: MoleculeCharge,
    ) -> Tuple[List[int], str]:
        """
        Get indices of atoms that are safe to adjust (Issue A).
        
        Returns:
            Tuple of (list of 0-based atom indices, method description)
        """
        if not mol_charge.atom_details:
            # No details available - cannot do safe-subset
            return [], "No atom details available"
        
        safe_indices = []
        for i, atom in enumerate(mol_charge.atom_details):
            atomname = atom.get('atomname', '')
            atomtype = atom.get('atomtype', '')
            charge = atom.get('charge', 0.0)
            
            # Check disallowed types
            if str(atomtype).casefold() in self.disallow_atomtypes:
                continue
            
            # Check |q| threshold
            if abs(charge) > self.max_abs_q_for_adjust:
                continue
            
            norm_name = _normalize_atomname(atomname)
            if _is_excluded_atomname(norm_name):
                continue

            # Check allowed names: user override or conservative default matcher
            if self.allowed_atomnames is not None:
                name_allowed = _matches_custom_allowlist(norm_name, self.allowed_atomnames)
            else:
                name_allowed = _matches_default_safe_atomname(norm_name)
            if not name_allowed:
                continue
            
            safe_indices.append(i)
        
        return safe_indices, f"safe_subset(n={len(safe_indices)})"

    def _build_moleculetype_rounding_report(
        self,
        molecule_counts: Dict[str, int],
    ) -> Tuple[List[Dict[str, Any]], bool, List[Dict[str, Any]]]:
        """
        Build per-moleculetype charge diagnostics and rounding-only decision.

        Returns:
            (rows, rounding_only_passed, dominant_contributors)
        """
        count_by_casefold = {
            _normalize_molname(name): int(count) for name, count in molecule_counts.items()
        }
        rows: List[Dict[str, Any]] = []
        for name in sorted(self.molecule_charges, key=lambda s: s.casefold()):
            mol = self.molecule_charges[name]
            mol_q = mol.total_charge
            nearest_int = int(round(mol_q))
            d = abs(mol_q - nearest_int)
            count = count_by_casefold.get(
                _normalize_molname(name),
                count_by_casefold.get(_normalize_molname(mol.itp_moleculetype_name), 0),
            )
            contribution = count * mol_q
            rows.append(
                {
                    "moleculetype": name,
                    "molecule_alias": name,
                    "itp_moleculetype_name": mol.itp_moleculetype_name,
                    "q_mol": mol_q,
                    "q_round": nearest_int,
                    "distance_to_integer": d,
                    "count": count,
                    "q_system_contribution": contribution,
                    "rounding_pass": d <= self.rounding_moleculetype_tol,
                }
            )

        active_rows = [r for r in rows if int(r["count"]) > 0]
        if not active_rows:
            active_rows = rows
        rounding_only_passed = all(
            float(r["distance_to_integer"]) <= self.rounding_moleculetype_tol
            for r in active_rows
        )
        dominant_contributors = sorted(
            active_rows,
            key=lambda r: abs(float(r["q_system_contribution"])),
            reverse=True,
        )[:5]
        return rows, rounding_only_passed, dominant_contributors

    def _print_moleculetype_rounding_report(
        self,
        rows: List[Dict[str, Any]],
        dominant_contributors: List[Dict[str, Any]],
    ) -> None:
        """Print a compact table-like report for charge diagnostics."""
        print("  Moleculetype charge audit:")
        print("    alias->itp           Q_mol        round(Q)     d=|Q-round|      count    Q_system_contrib")
        print("    ---------------------------------------------------------------------------------------------")
        for row in rows:
            alias_name = str(row.get("molecule_alias") or row.get("moleculetype") or "")
            itp_name = str(row.get("itp_moleculetype_name") or alias_name)
            name_label = alias_name if alias_name.casefold() == itp_name.casefold() else f"{alias_name}->{itp_name}"
            print(
                f"    {name_label[:20]:<20} "
                f"{float(row['q_mol']):+12.6f} "
                f"{int(row['q_round']):>12d} "
                f"{float(row['distance_to_integer']):>14.6e} "
                f"{int(row['count']):>10d} "
                f"{float(row['q_system_contribution']):+17.6e}"
            )
        if dominant_contributors:
            labels = ", ".join(
                f"{r['moleculetype']}({float(r['q_system_contribution']):+.3e})"
                for r in dominant_contributors
            )
            print(f"  Dominant |Q_system| contributors: {labels}")

    def _selected_hetero_atom_details(
        self,
        mol_charge: MoleculeCharge,
        selected_indices: List[int],
    ) -> List[Dict[str, Any]]:
        """Return hetero-atom details for selected indices."""
        details: List[Dict[str, Any]] = []
        if not mol_charge.atom_details:
            return details
        for i in selected_indices:
            if i < 0 or i >= len(mol_charge.atom_details):
                continue
            atom = mol_charge.atom_details[i]
            atomname = str(atom.get("atomname", ""))
            atomtype = str(atom.get("atomtype", ""))
            elem = self._atom_primary_element(atomname, atomtype)
            if elem not in HETERO_ELEMENTS:
                continue
            details.append(
                {
                    "idx": i,
                    "nr": atom.get("nr"),
                    "atomname": atomname,
                    "atomtype": atomtype,
                    "element": elem,
                }
            )
        return details

    def _hetero_adjustment_allowed(
        self,
        hetero_details: List[Dict[str, Any]],
    ) -> Tuple[bool, List[Dict[str, Any]]]:
        """
        Validate that any hetero atoms selected are explicitly allowlisted.

        Returns:
            (allowed, not_allowlisted_hetero)
        """
        if not hetero_details:
            return True, []
        if self.allowed_atomnames is None:
            return False, hetero_details
        not_allowlisted: List[Dict[str, Any]] = []
        for item in hetero_details:
            atomname = _normalize_atomname(str(item.get("atomname", "")))
            if not _matches_custom_allowlist(atomname, self.allowed_atomnames):
                not_allowlisted.append(item)
        return len(not_allowlisted) == 0, not_allowlisted
    
    def check_neutrality(
        self,
        molecule_counts: Dict[str, int],
    ) -> Tuple[bool, float]:
        """
        Check if system is charge-neutral.
        
        Args:
            molecule_counts: Dict of molecule name to count
            
        Returns:
            Tuple of (is_neutral, total_charge)
        """
        q = self.compute_system_charge(molecule_counts)
        # Task I: Use configurable neutrality tolerance
        is_neutral = abs(q) < self.neutrality_tol
        return is_neutral, q
    
    def can_auto_correct(self, q: float) -> bool:
        """
        Check if charge imbalance is small enough for auto-correction.
        
        Args:
            q: Total system charge
            
        Returns:
            True if |q| <= threshold
        """
        return abs(q) <= self.threshold
    
    def apply_correction(
        self,
        q: float,
        molecule_counts: Dict[str, int],
        output_dir: Path,
        run_id: str,
        target_molecule: Optional[str] = None,
    ) -> CorrectionResult:
        """
        Apply charge correction to neutralize the system.
        
        Strategy: Use safe-subset atom selection (default) or uniform distribution.
        
        Args:
            q: Total system charge to correct
            molecule_counts: Dict of molecule name to count
            output_dir: Directory for patched ITP
            run_id: Run identifier
            target_molecule: Molecule to adjust (REQUIRED if require_explicit_target)
            
        Returns:
            CorrectionResult
        """
        thresholds_used = self._thresholds_used()

        def _make_result(**kwargs) -> CorrectionResult:
            kwargs.setdefault("original_q", q)
            kwargs.setdefault("run_id", run_id)
            kwargs.setdefault("thresholds_used", thresholds_used)
            if "unparsed_atoms_lines" not in kwargs and self.unparsed_atoms_audit:
                capped = self.unparsed_atoms_audit[: self.MAX_UNPARSED_ATOMS_AUDIT]
                kwargs["unparsed_atoms_lines"] = capped
                kwargs.setdefault(
                    "unparsed_atoms_lines_truncated",
                    max(0, len(self.unparsed_atoms_audit) - len(capped)),
                )
            if "neutrality_status" not in kwargs:
                if kwargs.get("corrected"):
                    kwargs["neutrality_status"] = "corrected"
                elif abs(q) < self.neutrality_tol:
                    kwargs["neutrality_status"] = "neutral"
                elif kwargs.get("correction_refused"):
                    kwargs["neutrality_status"] = "failed"
                else:
                    kwargs["neutrality_status"] = "unknown"
            return CorrectionResult(**kwargs)

        if abs(q) < 1e-10:
            return _make_result(
                corrected=False,
                correction_method="none (already neutral)",
                correction_delta=0.0,
                total_delta=0.0,
                atoms_adjusted=0,
                neutrality_status="neutral",
            )
        
        # Check if imbalance is too large
        if abs(q) > self.threshold:
            # Check if imbalance resembles missing ion
            q_rounded = round(q * 2) / 2  # Round to nearest 0.5
            if abs(q - q_rounded) < 0.1 and abs(q) >= 0.4:
                return _make_result(
                    corrected=False,
                    correction_refused=True,
                    refuse_reason=(
                        f"Charge imbalance |Q|={abs(q):.4f} resembles missing ion "
                        f"(~{q_rounded:+.1f}e). Check molecule counts or add missing species."
                    ),
                    correction_method="REFUSED: resembles missing ion",
                )
            return _make_result(
                corrected=False,
                correction_refused=True,
                refuse_reason=f"Charge imbalance |Q|={abs(q):.6f} too large (threshold={self.threshold})",
                correction_method=f"REFUSED: |Q| > threshold",
            )
        
        if abs(q) > self.max_total_delta:
            return _make_result(
                corrected=False,
                correction_refused=True,
                refuse_reason=(
                    f"Total charge correction |Q|={abs(q):.6f} exceeds "
                    f"max_total_delta={self.max_total_delta:.6f}"
                ),
                correction_method="REFUSED: exceeds max_total_delta",
            )
        
        # === Target Selection (Issue B) ===
        target = target_molecule
        
        if target is None and self.require_explicit_target:
            return _make_result(
                corrected=False,
                correction_refused=True,
                refuse_reason=(
                    "No target molecule specified and require_explicit_target=True. "
                    "Specify target_molecule or set a target_allowlist."
                ),
                correction_method="REFUSED: no explicit target",
            )
        
        if target is not None:
            resolved_target = self._resolve_target_name(target)
            if resolved_target is None:
                return _make_result(
                    corrected=False,
                    correction_refused=True,
                    refuse_reason=f"Specified target '{target}' not found in loaded molecules",
                    correction_method="REFUSED: target not found",
                )
            target = resolved_target
        
        if target is None:
            # Auto-select from allowlist only (Issue B)
            candidates = []
            for name, count in molecule_counts.items():
                resolved_name = self._resolve_target_name(name)
                if resolved_name is None:
                    continue
                mol = self.molecule_charges[resolved_name]
                # Skip ions and non-neutral molecules
                if self.is_ion(mol) or not self.is_neutral(mol):
                    continue
                # Must be in allowlist (Issue B)
                if not self._is_allowlisted_target(resolved_name):
                    continue
                candidates.append((resolved_name, count, mol))
            
            if not candidates:
                return _make_result(
                    corrected=False,
                    correction_refused=True,
                    refuse_reason=(
                        "No suitable correction target found. "
                        f"Checked allowlist: {sorted(self.target_allowlist)}. "
                        "Either specify target_molecule explicitly or add molecules to allowlist."
                    ),
                    correction_method="REFUSED: no allowed target",
                )
            
            # Sort by count (highest first)
            candidates.sort(key=lambda x: -x[1])
            target = candidates[0][0]
            print(f"    Auto-selected correction target: {target} (count={candidates[0][1]})")
        
        # Verify selected target
        mol_charge = self.molecule_charges[target]
        target_itp_name = mol_charge.itp_moleculetype_name
        
        # Task G: Ion detection is advisory - warn and skip unless allow_ions=True
        if self.is_ion(mol_charge):
            if not self.allow_ions:
                # Advisory: warn and skip, not a hard error
                print(f"    [WARN] Target '{target}' appears to be an ion (q={mol_charge.total_charge:+.2f})")
                print(f"    [WARN] Skipping ion correction. Use --charge-fix-allow-ions to override.")
                return _make_result(
                    corrected=False,
                    correction_refused=True,
                    refuse_reason=f"Target '{target}' appears to be an ion (q={mol_charge.total_charge:+.2f}). Use --charge-fix-allow-ions to override.",
                    correction_method="SKIPPED: target is ion (advisory)",
                )
            else:
                # User explicitly allowed ion correction
                print(f"    [WARN] Target '{target}' appears to be an ion (q={mol_charge.total_charge:+.2f}), but --allow-ions is set")
        
        # Task E: Preprocessor directive fail-fast (strict mode) or skip (non-strict)
        if mol_charge.has_preprocessor:
            if self.strict_mode:
                # Strict mode: fail-fast with actionable error
                raise ChargeNeutralityError(
                    f"Target '{target}' ITP contains preprocessor directives (#ifdef/#include/#define) "
                    f"in or around [ atoms ] section.\n"
                    f"Cannot safely compute or correct charge in strict mode.\n"
                    f"Solutions:\n"
                    f"  1. Preprocess the ITP: grompp -pp system_preprocessed.top\n"
                    f"  2. Use non-strict mode (no --strict-charge-neutrality flag)"
                )
            else:
                # Non-strict: warn and mark as unknown/skip correction
                print(f"    [WARN] Target '{target}' ITP contains preprocessor directives - cannot safely correct")
                return _make_result(
                    corrected=False,
                    correction_refused=True,
                    refuse_reason=(
                        f"Target '{target}' ITP contains preprocessor directives (#ifdef/#include/#define) "
                        "in or around [ atoms ] section. Cannot safely auto-correct."
                    ),
                    correction_method="SKIPPED: preprocessor detected",
                    neutrality_status="unknown",
                )
        
        # Check for suspicious patterns (Issue B)
        suspicious_warnings = self._is_suspicious_target(target, mol_charge)
        for warn in suspicious_warnings:
            print(f"    [WARN] {warn}")
        
        count_by_casefold = {
            _normalize_molname(name): count for name, count in molecule_counts.items()
        }
        target_count = count_by_casefold.get(
            _normalize_molname(target),
            count_by_casefold.get(_normalize_molname(target_itp_name), 0),
        )
        
        if target_count == 0 or mol_charge.atom_count == 0:
            return _make_result(
                corrected=False,
                correction_refused=True,
                refuse_reason=f"Target {target} has count=0 or no atoms",
                correction_method="REFUSED: empty target",
            )
        
        # === Atom Selection (Issue A) ===
        hetero_selected_details: List[Dict[str, Any]] = []
        is_polymer_target = self._is_protected_polymer(target)
        effective_method = self.correction_method
        if is_polymer_target:
            # Polymer corrections (when explicitly allowed) use spread-safe distribution.
            effective_method = "spread_safe"

        if effective_method == "safe_subset":
            safe_indices, method_desc = self._get_safe_atom_indices(mol_charge)
            
            if not safe_indices:
                sample_names = []
                if mol_charge.atom_details:
                    names = {
                        _normalize_atomname(a.get('atomname', ''))
                        for a in mol_charge.atom_details
                        if a.get('atomname')
                    }
                    sample_names = sorted(n for n in names if n)[:8]
                return _make_result(
                    corrected=False,
                    correction_refused=True,
                    refuse_reason=(
                        f"No safe atoms found in target '{target}' for correction. "
                        f"Sample atom names: {sample_names if sample_names else 'N/A'}. "
                        f"Criteria: conservative H/C name matching, |q| < {self.max_abs_q_for_adjust}, "
                        f"not in disallowed atomtypes. Consider expanding the allowlist "
                        f"(--charge-fix-target-atomnames) or selecting a different target molecule."
                    ),
                    correction_method="REFUSED: no safe atoms",
                )
            
            # Distribute across safe atoms only
            total_atoms_to_adjust = target_count * len(safe_indices)
            delta_per_atom = -q / total_atoms_to_adjust
        elif effective_method == "spread_safe":
            safe_indices, method_desc = self._get_spread_safe_atom_indices(mol_charge)

            if not safe_indices:
                return _make_result(
                    corrected=False,
                    correction_refused=True,
                    refuse_reason=(
                        f"No spread-safe atoms found in protected polymer target '{target}': {method_desc}. "
                        "spread_safe excludes hetero atoms (O/N/S/B/F/P) and atoms within "
                        f"{self.polymer_exclusion_bonds} bonds of those atoms."
                    ),
                    correction_method="REFUSED: no spread_safe atoms",
                )

            total_atoms_to_adjust = target_count * len(safe_indices)
            delta_per_atom = -q / total_atoms_to_adjust

        else:  # uniform_all
            safe_indices = list(range(mol_charge.atom_count))
            method_desc = "uniform_all"
            total_atoms_to_adjust = target_count * mol_charge.atom_count
            delta_per_atom = -q / total_atoms_to_adjust

        hetero_selected_details = self._selected_hetero_atom_details(mol_charge, safe_indices)
        hetero_allowed, not_allowlisted_hetero = self._hetero_adjustment_allowed(hetero_selected_details)
        if not hetero_allowed:
            sample = [
                f"{item.get('atomname', '?')}/{item.get('atomtype', '?')}"
                for item in not_allowlisted_hetero[:8]
            ]
            return _make_result(
                corrected=False,
                correction_refused=True,
                refuse_reason=(
                    "Correction would modify hetero atoms without explicit allowlist. "
                    "Provide --charge-fix-target-atomnames including the exact hetero atom names "
                    f"if this is intentional. Examples: {sample}"
                ),
                correction_method="REFUSED: hetero atoms not allowlisted",
            )
        
        if abs(delta_per_atom) > self.max_delta_per_atom:
            return _make_result(
                corrected=False,
                correction_refused=True,
                refuse_reason=(
                    f"Per-atom correction |Δq|={abs(delta_per_atom):.6e} exceeds "
                    f"max_delta_per_atom={self.max_delta_per_atom:.6e}"
                ),
                correction_method="REFUSED: delta too large per atom",
            )
        
        # Build patched_atoms list with details
        patched_charges = []
        patched_atom_details: List[Dict[str, Any]] = []
        patched_atom_details_total = 0
        for i in safe_indices:
            orig_q = mol_charge.charges[i]
            new_q = orig_q + delta_per_atom
            patched_charges.append((i, orig_q, new_q))
            if mol_charge.atom_details:
                patched_atom_details_total += 1
                if len(patched_atom_details) >= self.MAX_PATCHED_ATOM_DETAILS:
                    continue
                atom = mol_charge.atom_details[i]
                patched_atom_details.append({
                    'nr': atom.get('nr'),
                    'atomname': atom.get('atomname'),
                    'atomtype': atom.get('atomtype'),
                    'line_no': atom.get('line_no'),
                    'charge_col_idx': atom.get('charge_col_idx'),
                    'line_excerpt': atom.get('line_excerpt'),
                    'old_q': orig_q,
                    'new_q': new_q,
                })
        patched_atom_details_truncated = max(
            0,
            patched_atom_details_total - len(patched_atom_details),
        )
        
        # Write patched ITP (Issue D: only patch target moleculetype)
        patched_itp_path = self._write_patched_itp(
            mol_charge.itp_path,
            delta_per_atom,
            safe_indices,  # Pass safe indices for selective patching
            output_dir,
            run_id,
            target_itp_name,
            target_count=target_count,
            original_system_q=q,
            effective_method=effective_method,
            method_desc=method_desc,
        )
        
        # === Verification Phase ===
        # Re-parse patched ITP and verify molecule charge matches expected value
        new_mol_charge = ChargeParser.get_molecule_charge(
            patched_itp_path,
            target,
            strict=self.strict_mode,
            expected_itp_moleculetype_name=target_itp_name,
            allow_unparsed_atoms_lines=self.allow_unparsed_atoms_lines,
        )
        expected_mol_q = mol_charge.total_charge + (delta_per_atom * len(safe_indices))
        actual_mol_q = new_mol_charge.total_charge
        
        mol_verify_tol = 1e-10
        if abs(actual_mol_q - expected_mol_q) > mol_verify_tol:
            return _make_result(
                corrected=False,
                correction_refused=True,
                refuse_reason=(
                    f"Write verification failed for molecule '{target}'. "
                    f"Expected q={expected_mol_q:.12e}, got q={actual_mol_q:.12e}. "
                    f"Check write precision or ITP format."
                ),
                correction_method="REFUSED: verification failed",
            )
        
        # Update cached molecule charge to reflect patched version
        self.molecule_charges[target] = new_mol_charge
        
        # Verify system is now neutral
        try:
            recomputed_q = self.compute_system_charge(molecule_counts, strict=True)
        except ChargeNeutralityError as e:
            return _make_result(
                corrected=False,
                correction_refused=True,
                refuse_reason=f"System charge verification failed: {e}",
                correction_method="REFUSED: system verification failed",
            )
        
        system_verify_tol = 1e-9
        if abs(recomputed_q) > system_verify_tol:
            return _make_result(
                corrected=False,
                correction_refused=True,
                refuse_reason=(
                    f"System still not neutral after correction. "
                    f"Residual Q={recomputed_q:.12e} (tolerance={system_verify_tol:.0e})"
                ),
                correction_method="REFUSED: residual charge",
            )
        
        print(f"    [OK] Verification passed: molecule q={actual_mol_q:.12e}, system Q={recomputed_q:.12e}")
        unparsed_audit = self.unparsed_atoms_audit[: self.MAX_UNPARSED_ATOMS_AUDIT]
        unparsed_audit_truncated = max(
            0,
            len(self.unparsed_atoms_audit) - len(unparsed_audit),
        )
        correction_method_text = (
            f"{effective_method}: distributed {-q:.12e} across "
            f"{total_atoms_to_adjust} atoms of {target} "
            f"(itp={target_itp_name}; {method_desc})"
        )
        
        return _make_result(
            corrected=True,
            correction_method=correction_method_text,
            effective_method=effective_method,
            method_desc=method_desc,
            correction_delta=delta_per_atom,
            patched_itp_path=patched_itp_path,
            patched_atoms=patched_charges,
            max_delta_applied=abs(delta_per_atom),
            total_delta=abs(q),
            atoms_adjusted=total_atoms_to_adjust,
            target_molecule=target,
            target_itp_moleculetype_name=target_itp_name,
            patched_atom_details=patched_atom_details if patched_atom_details else None,
            patched_atom_details_truncated=patched_atom_details_truncated,
            safe_subset_used=(effective_method == "safe_subset"),
            unparsed_atoms_lines=unparsed_audit if unparsed_audit else None,
            unparsed_atoms_lines_truncated=unparsed_audit_truncated,
            hetero_atoms_modified=bool(hetero_selected_details),
            hetero_atoms_modified_count=len(hetero_selected_details),
            hetero_atoms_modified_details=hetero_selected_details if hetero_selected_details else None,
        )
    
    def _write_patched_itp(
        self,
        original_itp: Path,
        delta_per_atom: float,
        safe_indices: List[int],
        output_dir: Path,
        run_id: str,
        mol_name: str,
        target_count: int,
        original_system_q: float,
        effective_method: str,
        method_desc: str,
    ) -> Path:
        """
        Write a patched ITP file with adjusted charges using token-based replacement.
        
        Task D: Uses token-based replacement instead of fragile regex.
        Task H: Includes comprehensive audit header and run-scoped filename.
        
        Args:
            original_itp: Path to original ITP
            delta_per_atom: Charge adjustment per atom
            safe_indices: List of 0-based atom indices to adjust
            output_dir: Output directory
            run_id: Run identifier
            mol_name: Target molecule name (only this moleculetype is patched)
            target_count: Number of target molecules in the system
            original_system_q: Original total system charge before correction
            effective_method: Effective correction method used
            method_desc: Method description/audit detail
            
        Returns:
            Path to patched ITP
        """
        output_dir.mkdir(parents=True, exist_ok=True)
        content = original_itp.read_text(encoding="utf-8")
        lines = content.splitlines()
        result = []
        
        target_key = _normalize_molname(mol_name)
        current_section = ""
        current_moleculetype: Optional[str] = None
        current_moleculetype_is_target = False
        in_atoms = False
        atom_idx = 0  # 0-based index within current moleculetype's atoms section
        safe_set = set(safe_indices)
        atoms_patched = 0
        target_atoms_sections = 0
        
        for line_no, line in enumerate(lines, 1):
            # Check for section header
            match = ChargeParser.SECTION_RE.match(line)
            if match:
                section = match.group(1).lower()
                current_section = section
                if section == 'moleculetype':
                    # Reset for new moleculetype
                    in_atoms = False
                    atom_idx = 0
                    current_moleculetype = None
                    current_moleculetype_is_target = False
                elif section == 'atoms':
                    in_atoms = current_moleculetype_is_target
                    atom_idx = 0
                    if in_atoms:
                        target_atoms_sections += 1
                else:
                    in_atoms = False
                result.append(line)
                continue
            
            # Parse moleculetype name
            if current_section == "moleculetype" and current_moleculetype is None and not in_atoms:
                stripped = line.strip()
                if stripped and not stripped.startswith(';') and not stripped.startswith('#'):
                    if not stripped.startswith('['):
                        data, _ = ChargeParser._split_data_comment(line)
                        if data:
                            parts = data.split()
                            if parts:
                                current_moleculetype = parts[0]
                                current_moleculetype_is_target = (
                                    _normalize_molname(current_moleculetype) == target_key
                                )
            
            # Only patch atoms in target moleculetype (Task D)
            if not in_atoms:
                result.append(line)
                continue
            
            # Skip comments/empty/preprocessor lines; only data lines are parse candidates.
            if not ChargeParser._is_atoms_data_line(line):
                result.append(line)
                continue
            
            # Task D: Use robust parser for token-based replacement
            parsed = ChargeParser._parse_atoms_line(
                line,
                line_no,
                strict=self.strict_mode,
            )
            if parsed is None:
                if self.allow_unparsed_atoms_lines:
                    result.append(line)
                    continue
                raise ChargeNeutralityError(
                    "Unparseable [ atoms ] data line while patching target moleculetype "
                    f"'{mol_name}' in {original_itp} at line {line_no}: "
                    f"{ChargeParser._line_excerpt(line)}"
                )
            
            # Only adjust if this atom index is in safe set
            if atom_idx in safe_set:
                old_charge = parsed['charge']
                new_charge = old_charge + delta_per_atom
                charge_col_idx = parsed['charge_col_idx']
                tokens = parsed['tokens']
                comment = parsed['comment']
                
                # Task D: Token-based replacement - replace the charge token
                # Preserve original line indentation
                leading_ws = line[:len(line) - len(line.lstrip())]
                
                # Format new charge with high precision (12 decimals)
                tokens[charge_col_idx] = f"{new_charge:.12e}"
                
                # Reconstruct line with normalized spacing + original comment
                new_data = "    ".join(tokens)  # Use reasonable spacing
                if comment:
                    new_line = f"{leading_ws}{new_data} {comment}"
                else:
                    new_line = f"{leading_ws}{new_data}"
                
                result.append(new_line)
                atoms_patched += 1
            else:
                result.append(line)
            
            atom_idx += 1

        if target_atoms_sections == 0:
            raise ChargeNeutralityError(
                f"Target moleculetype '{mol_name}' not found in {original_itp} while writing patched ITP."
            )
        if target_atoms_sections > 1:
            raise ChargeNeutralityError(
                f"Ambiguous target moleculetype '{mol_name}' in {original_itp}: "
                f"found {target_atoms_sections} [ atoms ] sections matching case-insensitively."
            )
        
        # Task D: Verify we patched the expected number of atoms
        if atoms_patched != len(safe_indices):
            raise ChargeNeutralityError(
                f"Patching mismatch: expected to patch {len(safe_indices)} atoms, "
                f"actually patched {atoms_patched}. This may indicate ITP parsing issues.\n"
                f"Target molecule: {mol_name}, ITP: {original_itp}\n"
                f"Safe indices: {sorted(safe_indices)[:20]}{'...' if len(safe_indices) > 20 else ''}"
            )
        
        # Task H: Comprehensive audit header with run-scoped metadata
        per_molecule_delta = delta_per_atom * atoms_patched
        system_delta = per_molecule_delta * target_count
        expected_system_delta = -original_system_q
        consistency_tol = 1e-9
        if not math.isclose(
            system_delta,
            expected_system_delta,
            rel_tol=0.0,
            abs_tol=consistency_tol,
        ):
            print(
                "    [WARN] Audit consistency check: system correction delta "
                f"{system_delta:+.12e} does not match -Q {expected_system_delta:+.12e} "
                f"(tol={consistency_tol:.0e})"
            )

        header = [
            "; ========================================",
            "; CHARGE CORRECTION AUDIT",
            "; ========================================",
            f"; Run ID: {run_id}",
            f"; Original file: {original_itp}",
            f"; Molecule: {mol_name}",
            f"; Method: {effective_method}",
            f"; Method detail: {method_desc}",
            f"; ",
            f"; Total system charge correction: {system_delta:+.12e}",
            f"; Delta per atom: {delta_per_atom:+.12e}",
            f"; Atoms patched per molecule: {atoms_patched}",
            f"; Target molecule count: {target_count}",
            f"; Per-molecule charge correction: {per_molecule_delta:+.12e}",
            f"; Expected system correction (-Q): {expected_system_delta:+.12e}",
            f"; Atom indices: {sorted(safe_indices)[:15]}{'...' if len(safe_indices) > 15 else ''}",
            f"; ",
            f"; Thresholds used:",
            f";   warn: {self.warn_threshold:.0e}",
            f";   correct: {self.threshold:.0e}",
            f";   max_delta_per_atom: {self.max_delta_per_atom:.0e}",
            f";   max_total_delta: {self.max_total_delta:.0e}",
            f";   neutrality_tol: {self.neutrality_tol:.0e}",
            f"; ",
            f"; Target allowlist: {sorted(list(self.target_allowlist)[:8])}{'...' if len(self.target_allowlist) > 8 else ''}",
            f"; Disallowed atomtypes: {list(self.disallow_atomtypes_tokens[:8])}{'...' if len(self.disallow_atomtypes_tokens) > 8 else ''}",
            "; ========================================",
            "",
        ]
        if self.allowed_atomnames is not None:
            header.insert(
                -2,
                f"; Allowed atomnames (override): "
                f"{sorted(list(self.allowed_atomnames)[:8])}"
                f"{'...' if len(self.allowed_atomnames) > 8 else ''}",
            )
        else:
            header.insert(
                -2,
                f"; Allowed H prefixes: "
                f"{sorted(list(self.allowed_h_prefixes)[:8])}"
                f"{'...' if len(self.allowed_h_prefixes) > 8 else ''}",
            )
            header.insert(
                -2,
                f"; Allowed C tokens: "
                f"{sorted(list(self.allowed_c_prefixes)[:8])}"
                f"{'...' if len(self.allowed_c_prefixes) > 8 else ''}",
            )
            header.insert(
                -2,
                f"; Excluded atom names: "
                f"{sorted(list(self.excluded_atom_names)[:8])}"
                f"{'...' if len(self.excluded_atom_names) > 8 else ''}",
            )
        
        output_content = '\n'.join(header + result) + '\n'
        
        # Task H: Run-scoped filename includes run_id
        output_path = output_dir / f"{mol_name}_corrected_{run_id}.itp"

        # Atomic write: temp file in same directory, fsync, then replace.
        temp_path: Optional[Path] = None
        try:
            fd, tmp_name = tempfile.mkstemp(
                dir=str(output_dir),
                prefix=f".{output_path.name}.tmp.",
                suffix=".itp",
            )
            temp_path = Path(tmp_name)
            with os.fdopen(fd, "w", encoding="utf-8", newline="\n") as fh:
                fh.write(output_content)
                fh.flush()
                os.fsync(fh.fileno())
            temp_path.replace(output_path)
        finally:
            if temp_path is not None and temp_path.exists():
                temp_path.unlink()
        
        return output_path
    
    def run(
        self,
        molecule_itp_paths: Dict[str, Path],
        molecule_counts: Dict[str, int],
        output_dir: Path,
        run_id: str,
        target_molecule: Optional[str] = None,
    ) -> CorrectionResult:
        """
        Run full charge neutrality check and correction.
        
        Args:
            molecule_itp_paths: Dict of molecule name to ITP path
            molecule_counts: Dict of molecule name to count
            output_dir: Output directory for patched ITPs
            run_id: Run identifier
            target_molecule: Optional explicit target molecule
            
        Returns:
            CorrectionResult
            
        Raises:
            ChargeNeutralityError: If |Q| > threshold and cannot auto-correct
        """
        print("  Checking charge neutrality...")
        
        # Load charges (with atom details for safe-subset)
        self.load_molecule_charges(molecule_itp_paths)

        per_moleculetype_report, rounding_only_passed, dominant_contributors = (
            self._build_moleculetype_rounding_report(molecule_counts)
        )
        self._print_moleculetype_rounding_report(
            per_moleculetype_report,
            dominant_contributors,
        )
        for row in per_moleculetype_report:
            if not row["rounding_pass"] and int(row["count"]) > 0:
                print(
                    f"    [WARN] {row['moleculetype']}: "
                    f"d=|Q-round(Q)|={float(row['distance_to_integer']):.3e} "
                    f"> tol={self.rounding_moleculetype_tol:.3e}"
                )

        def _attach_rounding_meta(result: CorrectionResult) -> CorrectionResult:
            result.rounding_only_passed = rounding_only_passed
            result.rounding_moleculetype_tol = self.rounding_moleculetype_tol
            result.per_moleculetype_report = per_moleculetype_report
            result.dominant_charge_contributors = dominant_contributors
            if self.unparsed_atoms_audit:
                capped = self.unparsed_atoms_audit[: self.MAX_UNPARSED_ATOMS_AUDIT]
                result.unparsed_atoms_lines = capped
                result.unparsed_atoms_lines_truncated = max(
                    0,
                    len(self.unparsed_atoms_audit) - len(capped),
                )
            return result

        # Compute total
        is_neutral, q = self.check_neutrality(molecule_counts)
        print(f"  Total system charge: Q = {q:+.8f}")
        
        if is_neutral:
            if not rounding_only_passed:
                print(
                    "  [WARN] System total charge is neutral but per-moleculetype integer drift "
                    "exceeds rounding tolerance."
                )
            print("  [OK] System is charge-neutral")
            return _attach_rounding_meta(CorrectionResult(
                original_q=q,
                corrected=False,
                correction_method="none (already neutral)",
                correction_delta=0.0,
                total_delta=0.0,
                atoms_adjusted=0,
                run_id=run_id,
                thresholds_used=self._thresholds_used(),
                neutrality_status="neutral",
            ))
        
        # Warn even if below correction threshold (Issue E)
        if abs(q) > self.warn_threshold:
            print(
                f"  [WARN] System charge |Q|={abs(q):.6f} "
                f"exceeds warning threshold {self.warn_threshold}"
            )

        if not rounding_only_passed:
            offenders = [
                row for row in per_moleculetype_report
                if int(row["count"]) > 0 and not row["rounding_pass"]
            ]
            offender_txt = ", ".join(
                f"{r['moleculetype']}(d={float(r['distance_to_integer']):.3e})"
                for r in offenders[:8]
            )
            msg = (
                "Detected non-rounding per-moleculetype charge drift. "
                f"Offenders: {offender_txt}. "
                "This indicates a topology/forcefield stitching issue, not rounding noise."
            )
            if not self.allow_non_rounding_correction:
                raise ChargeNeutralityError(
                    msg
                    + " Refusing auto-correction by default. "
                    "Use --charge-fix-allow-non-rounding to override explicitly."
                )
            print(
                "  [WARN] "
                + msg
                + " Proceeding due to explicit override (--charge-fix-allow-non-rounding)."
            )

        polymer_policy = self._polymer_policy_decision(q, molecule_counts, target_molecule)
        decision = str(polymer_policy.get("decision"))
        if decision == "acceptable_drift_skip":
            target = str(polymer_policy.get("target"))
            net_q = float(polymer_policy.get("net_q") or 0.0)
            tol = float(polymer_policy.get("tol") or self.polymer_net_charge_tol)
            print(
                f"  [INFO] Protected polymer drift accepted: target={target}, "
                f"net_q={net_q:+.6e}, tol={tol:.6e}. No correction applied."
            )
            return _attach_rounding_meta(CorrectionResult(
                original_q=q,
                corrected=False,
                correction_method="SKIPPED: protected polymer acceptable drift",
                correction_delta=0.0,
                total_delta=0.0,
                atoms_adjusted=0,
                run_id=run_id,
                thresholds_used=self._thresholds_used(),
                neutrality_status="polymer_acceptable_drift",
                target_molecule=target,
                polymer_policy=polymer_policy,
            ))
        if decision == "skip_non_strict":
            target = str(polymer_policy.get("target"))
            net_q = float(polymer_policy.get("net_q") or 0.0)
            print(
                f"  [WARN] Protected polymer '{target}' has net_q={net_q:+.6e} "
                f"> tol={self.polymer_net_charge_tol:.6e} and is not in "
                "--charge-fix-target-allowlist; skipping correction in non-strict mode."
            )
            return _attach_rounding_meta(CorrectionResult(
                original_q=q,
                corrected=False,
                correction_method="SKIPPED: protected polymer not allowlisted (non-strict)",
                correction_delta=0.0,
                total_delta=0.0,
                atoms_adjusted=0,
                run_id=run_id,
                thresholds_used=self._thresholds_used(),
                neutrality_status="polymer_skipped_non_strict",
                correction_refused=True,
                refuse_reason=(
                    f"Protected polymer '{target}' requires explicit --charge-fix-target-allowlist "
                    "for correction."
                ),
                target_molecule=target,
                polymer_policy=polymer_policy,
            ))
        if decision == "strict_fail":
            target = str(polymer_policy.get("target"))
            net_q = float(polymer_policy.get("net_q") or 0.0)
            raise ChargeNeutralityError(
                f"Protected polymer '{target}' has net_q={net_q:+.6e} "
                f"> polymer_net_charge_tol={self.polymer_net_charge_tol:.6e} and is not explicitly "
                "allowlisted. Add it to --charge-fix-target-allowlist and use "
                "--charge-fix-polymer-method spread_safe to permit correction."
            )
        if decision == "allow_spread_safe":
            resolved = self._resolve_target_name(str(polymer_policy.get("target")))
            if resolved is not None:
                target_molecule = resolved
                print(
                    f"  [INFO] Protected polymer correction enabled for '{target_molecule}' "
                    f"with spread_safe (allowlisted, exclusion_bonds={self.polymer_exclusion_bonds})."
                )
        
        if self.can_auto_correct(q):
            requested_target = self._resolve_target_name(target_molecule) if target_molecule else None
            if requested_target and self._is_protected_polymer(requested_target):
                method_label = "spread_safe"
            else:
                method_label = self.correction_method
            print(f"  [WARN] Small charge imbalance detected (|Q| <= {self.threshold})")
            print(f"  Applying auto-correction using method: {method_label}...")
            
            result = self.apply_correction(
                q, molecule_counts, output_dir, run_id, target_molecule
            )
            
            if result.corrected:
                print(f"  [OK] Correction applied: {result.correction_method}")
                print(f"  Patched ITP: {result.patched_itp_path}")
            elif result.correction_refused:
                print(f"  [WARN] Correction refused: {result.refuse_reason}")
                allowed_refusal_statuses = {"polymer_acceptable_drift", "polymer_skipped_non_strict"}
                if self.strict_mode or result.neutrality_status not in allowed_refusal_statuses:
                    raise ChargeNeutralityError(f"Charge correction refused: {result.refuse_reason}")
            else:
                print(f"  [FAIL] Correction failed: {result.correction_method}")

            if result.polymer_policy is None:
                result.polymer_policy = polymer_policy
            
            return _attach_rounding_meta(result)
        else:
            # Fail fast
            msg = (
                f"System charge neutrality check FAILED.\n"
                f"  Total charge Q = {q:+.8f}\n"
                f"  Threshold for auto-correction: |Q| <= {self.threshold}\n"
                f"  Please fix the charge imbalance manually."
            )
            raise ChargeNeutralityError(msg)
