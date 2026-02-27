"""
PathManager: Enforces strict path isolation for pipeline I/O.

- All input reads must come from IN/ directory only
- All output writes must go to OUT_GMX/<RUN_ID>/ only
- Provides validation and path construction utilities
"""

import os
import hashlib
from pathlib import Path
from typing import Optional, Union


class PathViolationError(Exception):
    """Raised when a path access violates IN/OUT isolation rules."""
    pass


class PathManager:
    """
    Manages and validates paths for the pipeline.
    
    Enforces:
    - Reads only from project_root/IN/
    - Writes only to project_root/OUT_GMX/<RUN_ID>/
    """
    
    def __init__(self, project_root: Union[str, Path]):
        """
        Initialize PathManager with project root directory.
        
        Args:
            project_root: Absolute path to the project root directory
        """
        self.project_root = Path(project_root).resolve()
        self.in_dir = self.project_root / "IN"
        self.out_dir = self.project_root / "OUT_GMX"
        
        # Validate project structure exists
        if not self.in_dir.exists():
            raise PathViolationError(f"IN directory does not exist: {self.in_dir}")
    
    def get_input_path(self, *parts: str) -> Path:
        """
        Construct and validate an input path under IN/.
        
        Args:
            *parts: Path components relative to IN/
            
        Returns:
            Absolute path under IN/
            
        Raises:
            PathViolationError: If path escapes IN/ directory
        """
        path = self.in_dir.joinpath(*parts).resolve()
        if not self._is_under(path, self.in_dir):
            raise PathViolationError(
                f"Input path escapes IN/ directory: {path}"
            )
        return path
    
    def get_output_path(self, run_id: str, *parts: str) -> Path:
        """
        Construct an output path under OUT_GMX/<RUN_ID>/.
        
        Args:
            run_id: Run identifier
            *parts: Path components relative to OUT_GMX/<RUN_ID>/
            
        Returns:
            Absolute path under OUT_GMX/<RUN_ID>/
            
        Raises:
            PathViolationError: If path escapes OUT_GMX/ directory
        """
        path = self.out_dir.joinpath(run_id, *parts).resolve()
        if not self._is_under(path, self.out_dir):
            raise PathViolationError(
                f"Output path escapes OUT_GMX/ directory: {path}"
            )
        return path
    
    def get_run_dir(self, run_id: str) -> Path:
        """Get the root directory for a specific run."""
        return self.get_output_path(run_id)
    
    def validate_input_path(self, path: Union[str, Path]) -> Path:
        """
        Validate that a path is under IN/ directory.
        
        Args:
            path: Path to validate
            
        Returns:
            Resolved absolute path
            
        Raises:
            PathViolationError: If path is not under IN/
        """
        resolved = Path(path).resolve()
        if not self._is_under(resolved, self.in_dir):
            raise PathViolationError(
                f"Path is not under IN/ directory: {resolved}"
            )
        return resolved
    
    def validate_output_path(self, path: Union[str, Path]) -> Path:
        """
        Validate that a path is under OUT_GMX/ directory.
        
        Args:
            path: Path to validate
            
        Returns:
            Resolved absolute path
            
        Raises:
            PathViolationError: If path is not under OUT_GMX/
        """
        resolved = Path(path).resolve()
        if not self._is_under(resolved, self.out_dir):
            raise PathViolationError(
                f"Path is not under OUT_GMX/ directory: {resolved}"
            )
        return resolved
    
    def ensure_output_dir(self, run_id: str, *parts: str) -> Path:
        """
        Ensure an output directory exists, creating it if necessary.
        
        Args:
            run_id: Run identifier
            *parts: Path components relative to OUT_GMX/<RUN_ID>/
            
        Returns:
            Path to the created/existing directory
        """
        path = self.get_output_path(run_id, *parts)
        path.mkdir(parents=True, exist_ok=True)
        return path
    
    def hash_file(self, path: Union[str, Path]) -> str:
        """
        Compute SHA256 hash of a file.
        
        Args:
            path: Path to file
            
        Returns:
            Hex-encoded SHA256 hash
        """
        path = Path(path)
        sha256 = hashlib.sha256()
        with open(path, "rb") as f:
            for chunk in iter(lambda: f.read(8192), b""):
                sha256.update(chunk)
        return sha256.hexdigest()
    
    def _is_under(self, path: Path, parent: Path) -> bool:
        """Check if path is under parent directory."""
        try:
            path.relative_to(parent)
            return True
        except ValueError:
            return False
    
    # === Convenience methods for common paths ===
    
    def get_molecule_dir(self, molecule_name: str) -> Path:
        """Get path to a molecule's directory under IN/molecules/."""
        return self.get_input_path("molecules", molecule_name)
    
    def get_system_dir(self, system_id: str) -> Path:
        """Get path to a system's directory under IN/systems/."""
        return self.get_input_path("systems", system_id)
    
    def get_forcefield_dir(self, ff_name: str) -> Path:
        """Get path to a forcefield directory under IN/forcefield/."""
        return self.get_input_path("forcefield", ff_name.lower().replace("-", ""))
    
    def get_stage_dir(self, run_id: str, stage_name: str) -> Path:
        """Get path to a stage's output directory."""
        stage_dirs = {
            "packmol": "01_packmol",
            "htpolynet": "02_htpolynet",
            "gmx_em": "03_gromacs/em",
            "gmx_eq": "03_gromacs/nvt",  # NVT equilibration
            "gmx_npt": "03_gromacs/npt",
            "gmx_prod": "03_gromacs/md",
            "analysis": "analysis",
        }
        subdir = stage_dirs.get(stage_name, stage_name)
        return self.get_output_path(run_id, subdir)
    
    def get_logs_dir(self, run_id: str) -> Path:
        """Get path to logs directory for a run."""
        return self.ensure_output_dir(run_id, "logs")
