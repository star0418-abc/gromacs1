"""
BaseStage: Abstract base class for all pipeline stages.
"""

from abc import ABC, abstractmethod
from pathlib import Path
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from ..context import PipelineContext


class BaseStage(ABC):
    """
    Abstract base class for pipeline stages.
    
    Each stage must implement:
    - name: Stage identifier
    - output_subdir: Directory under OUT_GMX/<RUN_ID>/ for outputs
    - run(): Execute the stage logic
    
    Provides:
    - Sentinel-based completion tracking (done.ok)
    - Skip logic for resume/force
    """
    
    @property
    @abstractmethod
    def name(self) -> str:
        """Stage identifier (e.g., 'packmol', 'gmx_em')."""
        pass
    
    @property
    @abstractmethod
    def output_subdir(self) -> str:
        """Directory name under OUT_GMX/<RUN_ID>/ for this stage's outputs."""
        pass
    
    @abstractmethod
    def run(self, ctx: "PipelineContext") -> bool:
        """
        Execute the stage.
        
        Args:
            ctx: Pipeline context with configuration and resources
            
        Returns:
            True if stage completed successfully, False otherwise
        """
        pass
    
    def get_output_dir(self, ctx: "PipelineContext") -> Path:
        """Get the output directory for this stage."""
        return ctx.ensure_output_dir(self.output_subdir)
    
    def get_sentinel_path(self, ctx: "PipelineContext") -> Path:
        """Get path to the done.ok sentinel file."""
        return self.get_output_dir(ctx) / "done.ok"
    
    def is_complete(self, ctx: "PipelineContext") -> bool:
        """Check if this stage has already completed (sentinel exists)."""
        return self.get_sentinel_path(ctx).exists()
    
    def mark_complete(self, ctx: "PipelineContext") -> None:
        """Mark the stage as complete by writing the sentinel file."""
        sentinel = self.get_sentinel_path(ctx)
        sentinel.parent.mkdir(parents=True, exist_ok=True)
        sentinel.write_text(f"Stage {self.name} completed successfully.\n")
        
        # Update manifest
        if ctx.manifest:
            ctx.manifest.set_stage_status(self.name, "completed")
    
    def mark_failed(self, ctx: "PipelineContext", error: str) -> None:
        """Mark the stage as failed in the manifest."""
        if ctx.manifest:
            ctx.manifest.set_stage_status(self.name, "failed", error)
    
    def should_skip(self, ctx: "PipelineContext") -> bool:
        """
        Determine if this stage should be skipped.
        
        Returns True if:
        - Stage is complete AND resume=True AND force=False
        """
        if ctx.force:
            return False
        if ctx.resume and self.is_complete(ctx):
            return True
        return False
    
    def execute(self, ctx: "PipelineContext") -> bool:
        """
        Execute the stage with skip/resume logic.
        
        This is the main entry point called by the dispatcher.
        
        Returns:
            True if stage completed (or was skipped), False if failed
        """
        print(f"\n{'='*60}")
        print(f"STAGE: {self.name}")
        print(f"{'='*60}")
        
        if self.should_skip(ctx):
            print(f"  [SKIP] Stage already complete (resume=True)")
            if ctx.manifest:
                ctx.manifest.set_stage_status(self.name, "skipped", "Already complete")
            return True
        
        if ctx.manifest:
            ctx.manifest.set_stage_status(self.name, "running")
        
        try:
            success = self.run(ctx)
            if success:
                # CRITICAL: dry-run must NOT create done.ok (read-only semantics)
                if not ctx.dry_run:
                    self.mark_complete(ctx)
                    print(f"  [OK] Stage {self.name} completed successfully")
                else:
                    print(f"  [DRY-RUN] Stage {self.name} would complete (no done.ok written)")
            else:
                self.mark_failed(ctx, "Stage returned False")
                print(f"  [FAIL] Stage {self.name} did not complete successfully")
            return success
        except Exception as e:
            self.mark_failed(ctx, str(e))
            print(f"  [ERROR] Stage {self.name} failed: {e}")
            raise
