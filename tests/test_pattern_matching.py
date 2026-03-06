#!/usr/bin/env python3
"""Test the _matches_protected_pattern function directly without importing pipeline."""
import sys
import re
from pathlib import Path
from typing import Set

# Read the function directly from the source file
sanitizer_path = Path(__file__).parent.parent / 'pipeline' / 'stages' / 'sanitizer.py'
content = sanitizer_path.read_text(encoding='utf-8')

# Extract the function using exec
exec_globals = {}
fn_start = content.find('def _matches_protected_pattern')
fn_end = content.find('\n@dataclass', fn_start)
fn_code = content[fn_start:fn_end]
exec(fn_code, exec_globals)
_matches_protected_pattern = exec_globals['_matches_protected_pattern']

# Read safe-atom matcher directly from charge_neutrality.py
charge_neutrality_path = Path(__file__).parent.parent / 'pipeline' / 'charge_neutrality.py'
cn_content = charge_neutrality_path.read_text(encoding='utf-8')
cn_start = cn_content.find('ALLOWED_H_PREFIXES')
cn_end = cn_content.find('def _matches_custom_allowlist', cn_start)
cn_code = cn_content[cn_start:cn_end]
cn_globals = {'re': re, 'Set': Set}
exec(cn_code, cn_globals)
_normalize_atomname = cn_globals['_normalize_atomname']
_matches_default_safe_atomname = cn_globals['_matches_default_safe_atomname']

# Test cases
def test_pattern_matching():
    print("Testing pattern matching...")
    
    # Short patterns (<=3 chars) require exact match after normalization
    assert _matches_protected_pattern('PC', frozenset({'PC'})) == True, 'PC exact match should work'
    assert _matches_protected_pattern('POLYMER_CHAIN', frozenset({'PC'})) == False, 'PC should NOT match POLYMER_CHAIN'
    assert _matches_protected_pattern('LIPID_GPC', frozenset({'PC'})) == False, 'PC should NOT match LIPID_GPC'
    assert _matches_protected_pattern('EC', frozenset({'EC'})) == True, 'EC exact match should work'
    assert _matches_protected_pattern('ELECTROCHEMICAL', frozenset({'EC'})) == False, 'EC should NOT match ELECTROCHEMICAL'
    print("  [OK] Short patterns require exact match")
    
    # Standard patterns - exact match
    assert _matches_protected_pattern('SOL', frozenset({'SOL', 'WATER'})) == True, 'SOL exact match should work'
    assert _matches_protected_pattern('SOLVENT', frozenset({'SOL'})) == False, 'SOL should NOT match SOLVENT (short pattern)'
    print("  [OK] SOL pattern matching works correctly")
    
    # Longer patterns - require word boundaries for substring match
    assert _matches_protected_pattern('WATER', frozenset({'WATER'})) == True, 'WATER exact match should work'
    # WATERBOX: 'WATER' is at boundary start but NOT at boundary end (BOX follows), so should NOT match
    assert _matches_protected_pattern('WATERBOX', frozenset({'WATER'})) == False, 'WATER should NOT match WATERBOX (no boundary after)'
    print("  [OK] Longer patterns require word boundaries")
    
    # Ion patterns
    assert _matches_protected_pattern('NA', frozenset({'NA'})) == True, 'NA exact match should work'
    assert _matches_protected_pattern('SODIUM', frozenset({'NA'})) == False, 'NA should NOT match SODIUM'
    assert _matches_protected_pattern('CL', frozenset({'CL'})) == True, 'CL exact match should work'
    assert _matches_protected_pattern('CHLORIDE', frozenset({'CL'})) == False, 'CL should NOT match CHLORIDE'
    print("  [OK] Ion patterns require exact match")

    # Requested protected-pattern edge cases
    assert _matches_protected_pattern('LI', frozenset({'LI'})) == True, 'LI exact match should work'
    assert _matches_protected_pattern('LI+', frozenset({'LI'})) == True, 'LI+ should normalize to LI'
    assert _matches_protected_pattern('TFSI', frozenset({'TFSI'})) == True, 'TFSI exact match should work'
    assert _matches_protected_pattern('TFSI-', frozenset({'TFSI'})) == True, 'TFSI- should match TFSI boundary'
    assert _matches_protected_pattern('TFSI_1', frozenset({'TFSI'})) == True, 'TFSI_1 should match TFSI boundary'
    assert _matches_protected_pattern('TFSI-alpha', frozenset({'TFSI'})) == True, 'TFSI-alpha should match TFSI boundary'
    assert _matches_protected_pattern('G4_ALT', frozenset({'G4_ALT'})) == True, 'Exact long pattern should survive separator normalization'
    assert _matches_protected_pattern('PEG4', frozenset({'PC'})) == False, 'PC should NOT match PEG4'
    assert _matches_protected_pattern('POLYMER_CHAIN', frozenset({'PC'})) == False, 'PC should NOT match POLYMER_CHAIN'
    assert _matches_protected_pattern('LIPID_GPC', frozenset({'PC'})) == False, 'PC should NOT match LIPID_GPC'
    print("  [OK] Edge-case protected-pattern matching stays tight")

    # Underscore/dash normalized patterns
    assert _matches_protected_pattern('P_C', frozenset({'PC'})) == True, 'P_C should match PC after normalization'
    assert _matches_protected_pattern('P-C', frozenset({'PC'})) == True, 'P-C should match PC after normalization'
    print("  [OK] Normalization removes underscores and dashes")
    
    print('All pattern matching tests PASSED!')

def test_safe_atom_name_matching():
    print("Testing safe atom name matching...")
    assert _matches_default_safe_atomname(_normalize_atomname('H1')) is True, 'H1 should be safe'
    assert _matches_default_safe_atomname(_normalize_atomname('HC')) is True, 'HC should be safe'
    assert _matches_default_safe_atomname(_normalize_atomname('HGA')) is True, 'HGA should be safe'
    assert _matches_default_safe_atomname(_normalize_atomname('C1')) is True, 'C1 should be safe'
    assert _matches_default_safe_atomname(_normalize_atomname('CT')) is True, 'CT should be safe'
    assert _matches_default_safe_atomname(_normalize_atomname('CL')) is False, 'CL should be excluded'
    assert _matches_default_safe_atomname(_normalize_atomname('CA')) is False, 'CA should be excluded'
    assert _matches_default_safe_atomname(_normalize_atomname('CU')) is False, 'CU should be excluded'
    assert _matches_default_safe_atomname(_normalize_atomname('NA')) is False, 'NA should be excluded'
    print("  [OK] Conservative safe atom matcher works")

if __name__ == '__main__':
    test_pattern_matching()
    test_safe_atom_name_matching()
