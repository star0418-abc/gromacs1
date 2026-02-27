"""
Shared utility functions for the pipeline.

This module provides centralized helpers to ensure consistent behavior
across all pipeline stages.
"""

from typing import Any


def parse_bool(value: Any, default: bool = False) -> bool:
    """
    Parse a boolean from various input types.
    
    This is the canonical boolean parser for CLI/config flags. It prevents
    the common bug where bool("false") returns True in Python.
    
    Accepts:
    - bool: returned as-is
    - str: true/false/yes/no/on/off/1/0 (case-insensitive)
    - int: 0 = False, non-zero = True
    
    Returns default for None or unrecognized values.
    
    Examples:
        >>> parse_bool("false")
        False
        >>> parse_bool("true")
        True
        >>> parse_bool("FALSE")
        False
        >>> parse_bool(0)
        False
        >>> parse_bool(1)
        True
        >>> parse_bool(None, False)
        False
        >>> parse_bool("invalid", True)
        True
    """
    if value is None:
        return default
    if isinstance(value, bool):
        return value
    if isinstance(value, int):
        return value != 0
    if isinstance(value, str):
        lower = value.strip().lower()
        if lower in ("true", "yes", "on", "1"):
            return True
        if lower in ("false", "no", "off", "0"):
            return False
    return default
