"""
HTPolyNet configuration generator.

Generates system_config.yaml for HTPolyNet using Jinja2 templates.
Supports dynamic injection of composition, box size, and forcefield settings.
"""

from dataclasses import dataclass, field
import json
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Tuple

# Try to use Jinja2, fall back to pure-Python rendering if not available.
try:
    from jinja2 import Template

    JINJA2_AVAILABLE = True
except ImportError:
    Template = None
    JINJA2_AVAILABLE = False


# =============================================================================
# Constants
# =============================================================================

_POLICY_VALUES = {"off", "warn", "cap", "error"}
HTPOLYNET_SUPPORTS_SEARCH_RADIUS_SCHEDULE = False

_SOLVENT_BP_CATALOG: Tuple[Tuple[str, float, Tuple[str, ...]], ...] = (
    ("DMC", 363.15, ("DMC", "DIMETHYL CARBONATE")),
    ("DME", 358.15, ("DME", "DIMETHOXYETHANE", "1,2-DIMETHOXYETHANE")),
)


# =============================================================================
# Data Classes
# =============================================================================


@dataclass
class HTPolyNetConfig:
    """
    Configuration for HTPolyNet run.

    Attributes:
        project_name: Name for the HTPolyNet project
        forcefield: Forcefield (GAFF2 or OPLS-AA)
        initial_gro: Path to initial structure from PACKMOL
        workdir: Sandbox directory for HTPolyNet intermediates
        output_dir: Directory for final outputs
        rules_file: Path to reaction rules YAML

        molecule_counts: Dict of molecule counts from composition
        box_length_nm: Box size in nm
        box_vectors_nm: Box vectors (Lx, Ly, Lz) in nm

        max_conversion: Target crosslink conversion (0-1)
        search_radius_nm: Radius for finding reactive pairs
        initial_temperature_K: Starting temperature
        final_temperature_K: Annealing target temperature
        cure_temperature_K: Curing temperature target

        cure_temp_policy: Guard rail policy for cure temperature
        cure_temp_safety_margin_K: Safety margin below solvent boiling point
        allow_cure_above_solvent_bp: Explicit override for boiling-point risk
        allow_cure_above_solvent_bp_reason: Required reason when override is enabled

        conversion_policy: Guard rail policy for high max_conversion values
        recommended_max_conversion: Recommended upper conversion bound
        allow_high_conversion: Explicit override for high conversion
        allow_high_conversion_reason: Required reason when override is enabled

        search_radius_schedule: Optional staged schedule (only emitted when supported)

        reactive_site_counts: Per-species reactive site counts
        crosslinkers: Species with reactive_site_count >= 2
        chain_stoppers: Species with reactive_site_count == 1

        curing_steps: Optional list of step dicts
        gmx_maxwarn: Max warnings for grompp
        threads: Number of threads for GROMACS
    """

    project_name: str
    forcefield: str
    initial_gro: Path
    workdir: Path
    output_dir: Path
    rules_file: Optional[Path] = None

    # Composition
    molecule_counts: Dict[str, int] = field(default_factory=dict)
    box_length_nm: float = 5.0
    box_vectors_nm: Optional[Tuple[float, float, float]] = None

    # Crosslinking
    max_conversion: float = 0.85
    search_radius_nm: float = 0.35
    initial_temperature_K: float = 300.0
    final_temperature_K: float = 300.0
    cure_temperature_K: float = 350.0
    curing_steps: Optional[List[Dict[str, Any]]] = None
    search_radius_schedule: Optional[List[Dict[str, float]]] = None

    # Safety policy controls
    cure_temp_policy: str = "warn"
    cure_temp_safety_margin_K: float = 10.0
    allow_cure_above_solvent_bp: bool = False
    allow_cure_above_solvent_bp_reason: Optional[str] = None

    conversion_policy: str = "warn"
    recommended_max_conversion: float = 0.85
    allow_high_conversion: bool = False
    allow_high_conversion_reason: Optional[str] = None

    # Reactive-site metadata
    reactive_site_counts: Dict[str, int] = field(default_factory=dict)
    crosslinkers: Dict[str, Dict[str, Any]] = field(default_factory=dict)
    chain_stoppers: Dict[str, Dict[str, Any]] = field(default_factory=dict)

    # Resources
    gmx_maxwarn: int = 10
    threads: int = 4

    def __post_init__(self) -> None:
        if not (0.0 < self.max_conversion <= 1.0):
            raise ValueError(f"max_conversion must satisfy 0 < max_conversion <= 1, got {self.max_conversion}")
        if self.search_radius_nm <= 0.0:
            raise ValueError(f"search_radius_nm must be > 0, got {self.search_radius_nm}")

        for name in ("initial_temperature_K", "final_temperature_K", "cure_temperature_K"):
            value = getattr(self, name)
            if value <= 0.0:
                raise ValueError(f"{name} must be > 0, got {value}")

        if not isinstance(self.gmx_maxwarn, int) or isinstance(self.gmx_maxwarn, bool):
            raise ValueError(f"gmx_maxwarn must be an integer >= 0, got {self.gmx_maxwarn!r}")
        if self.gmx_maxwarn < 0:
            raise ValueError(f"gmx_maxwarn must be >= 0, got {self.gmx_maxwarn}")

        if not isinstance(self.threads, int) or isinstance(self.threads, bool):
            raise ValueError(f"threads must be an integer > 0, got {self.threads!r}")
        if self.threads <= 0:
            raise ValueError(f"threads must be > 0, got {self.threads}")

        if self.box_length_nm <= 0.0:
            raise ValueError(f"box_length_nm must be > 0, got {self.box_length_nm}")

        if self.box_vectors_nm is not None:
            if len(self.box_vectors_nm) != 3:
                raise ValueError(
                    f"box_vectors_nm must contain exactly 3 values (Lx, Ly, Lz), got {len(self.box_vectors_nm)}"
                )
            if any(v <= 0.0 for v in self.box_vectors_nm):
                raise ValueError(f"all box_vectors_nm entries must be > 0, got {self.box_vectors_nm}")

        for molecule_name, count in self.molecule_counts.items():
            if not isinstance(molecule_name, str) or not molecule_name.strip():
                raise ValueError(f"molecule_counts keys must be non-empty strings, got {molecule_name!r}")
            if not isinstance(count, int) or isinstance(count, bool):
                raise ValueError(
                    f"molecule_counts['{molecule_name}'] must be a non-negative integer, got {count!r}"
                )
            if count < 0:
                raise ValueError(
                    f"molecule_counts['{molecule_name}'] must be a non-negative integer, got {count}"
                )

        if self.cure_temp_policy not in _POLICY_VALUES:
            raise ValueError(
                f"cure_temp_policy must be one of {sorted(_POLICY_VALUES)}, got {self.cure_temp_policy!r}"
            )
        if self.conversion_policy not in _POLICY_VALUES:
            raise ValueError(
                f"conversion_policy must be one of {sorted(_POLICY_VALUES)}, got {self.conversion_policy!r}"
            )

        if self.cure_temp_safety_margin_K < 0.0:
            raise ValueError(
                f"cure_temp_safety_margin_K must be >= 0, got {self.cure_temp_safety_margin_K}"
            )
        if not (0.0 < self.recommended_max_conversion <= 1.0):
            raise ValueError(
                "recommended_max_conversion must satisfy 0 < recommended_max_conversion <= 1, "
                f"got {self.recommended_max_conversion}"
            )

        if self.search_radius_schedule is not None:
            if not isinstance(self.search_radius_schedule, list):
                raise ValueError("search_radius_schedule must be a list of dict entries")
            normalized_schedule: List[Dict[str, float]] = []
            for idx, entry in enumerate(self.search_radius_schedule):
                if not isinstance(entry, dict):
                    raise ValueError(
                        f"search_radius_schedule[{idx}] must be a dict with keys "
                        "'until_conversion' and 'search_radius_nm'"
                    )
                if "until_conversion" not in entry or "search_radius_nm" not in entry:
                    raise ValueError(
                        f"search_radius_schedule[{idx}] must contain keys "
                        "'until_conversion' and 'search_radius_nm'"
                    )
                until_conversion = float(entry["until_conversion"])
                search_radius = float(entry["search_radius_nm"])
                if not (0.0 < until_conversion <= 1.0):
                    raise ValueError(
                        f"search_radius_schedule[{idx}]['until_conversion'] must satisfy "
                        f"0 < until_conversion <= 1, got {until_conversion}"
                    )
                if search_radius <= 0.0:
                    raise ValueError(
                        f"search_radius_schedule[{idx}]['search_radius_nm'] must be > 0, got {search_radius}"
                    )
                normalized_schedule.append(
                    {
                        "until_conversion": until_conversion,
                        "search_radius_nm": search_radius,
                    }
                )
            normalized_schedule.sort(key=lambda item: item["until_conversion"])
            for idx in range(1, len(normalized_schedule)):
                if normalized_schedule[idx]["until_conversion"] <= normalized_schedule[idx - 1]["until_conversion"]:
                    raise ValueError(
                        "search_radius_schedule entries must have strictly increasing 'until_conversion' values"
                    )
            self.search_radius_schedule = normalized_schedule


# =============================================================================
# Helpers
# =============================================================================


def _yaml_quote(value: str) -> str:
    """Return a deterministic double-quoted YAML scalar."""
    return json.dumps(value)


def _normalize_for_substring_match(value: str) -> str:
    return "".join(ch for ch in value.upper() if ch.isalnum())


def _ordered_unique(values: Iterable[str]) -> List[str]:
    seen: Dict[str, None] = {}
    for value in values:
        if value not in seen:
            seen[value] = None
    return list(seen.keys())


def _sorted_molecule_counts(molecule_counts: Dict[str, int]) -> Dict[str, int]:
    return dict(sorted(molecule_counts.items(), key=lambda item: item[0]))


def _molecule_count_items(molecule_counts: Dict[str, int]) -> List[Dict[str, Any]]:
    return [{"name": name, "count": count} for name, count in molecule_counts.items()]


def _reactive_site_items(reactive_site_counts: Dict[str, int]) -> List[Dict[str, Any]]:
    return [
        {"name_q": _yaml_quote(name), "count": count}
        for name, count in sorted(reactive_site_counts.items(), key=lambda item: item[0])
    ]


def _crosslinker_items(crosslinkers: Dict[str, Dict[str, Any]]) -> List[Dict[str, Any]]:
    items: List[Dict[str, Any]] = []
    for name, spec in sorted(crosslinkers.items(), key=lambda item: item[0]):
        classification = spec.get("classification", "crosslinker")
        items.append(
            {
                "name_q": _yaml_quote(name),
                "count": spec.get("count"),
                "reactive_site_count": spec.get("reactive_site_count"),
                "classification_q": _yaml_quote(str(classification)),
            }
        )
    return items


def _chain_stopper_items(chain_stoppers: Dict[str, Dict[str, Any]]) -> List[Dict[str, Any]]:
    items: List[Dict[str, Any]] = []
    for name, spec in sorted(chain_stoppers.items(), key=lambda item: item[0]):
        items.append(
            {
                "name_q": _yaml_quote(name),
                "count": spec.get("count"),
                "reactive_site_count": spec.get("reactive_site_count"),
            }
        )
    return items


def _normalize_curing_steps(curing_steps: Optional[List[Dict[str, Any]]]) -> Optional[List[Dict[str, Any]]]:
    if not curing_steps:
        return None
    normalized: List[Dict[str, Any]] = []
    for idx, step in enumerate(curing_steps):
        if not isinstance(step, dict):
            raise ValueError(f"curing_steps[{idx}] must be a dict, got {type(step).__name__}")
        if "temperature_K" not in step or "duration_ps" not in step:
            raise ValueError(
                f"curing_steps[{idx}] must include keys 'temperature_K' and 'duration_ps', got {sorted(step.keys())}"
            )
        normalized_step: Dict[str, Any] = {
            "temperature_K": step["temperature_K"],
            "duration_ps": step["duration_ps"],
            "ramp_ps": step.get("ramp_ps"),
        }
        normalized.append(normalized_step)
    return normalized


def _detect_known_solvents(molecule_counts: Dict[str, int]) -> List[Dict[str, Any]]:
    detected: List[Dict[str, Any]] = []
    normalized_names = {
        name: _normalize_for_substring_match(name)
        for name, count in molecule_counts.items()
        if count > 0
    }
    for solvent_name, bp_k, aliases in _SOLVENT_BP_CATALOG:
        alias_patterns = [_normalize_for_substring_match(alias) for alias in aliases]
        matched_names = sorted(
            name
            for name, norm_name in normalized_names.items()
            if any(alias in norm_name for alias in alias_patterns)
        )
        if matched_names:
            detected.append(
                {
                    "name": solvent_name,
                    "bp_K": bp_k,
                    "matched_names": matched_names,
                }
            )
    return detected


def _solvent_bp_summary(detected_solvents: List[Dict[str, Any]]) -> str:
    parts = [f"{item['name']} ({item['bp_K']:.0f} K)" for item in detected_solvents]
    return " / ".join(parts)


def _resolve_cure_temperature_guardrail(
    config: HTPolyNetConfig,
    detected_solvents: List[Dict[str, Any]],
) -> Tuple[float, List[str], Optional[str]]:
    effective_temp = config.cure_temperature_K
    warning_comments: List[str] = []
    adjustment_comment: Optional[str] = None

    if not detected_solvents:
        return effective_temp, warning_comments, adjustment_comment

    min_bp_k = min(item["bp_K"] for item in detected_solvents)
    threshold_k = min_bp_k - config.cure_temp_safety_margin_K
    if threshold_k <= 0.0:
        raise ValueError(
            "cure_temp_safety_margin_K is too large for detected solvents; "
            f"bp_min={min_bp_k:.2f} K, margin={config.cure_temp_safety_margin_K:.2f} K."
        )

    if config.cure_temperature_K <= threshold_k:
        return effective_temp, warning_comments, adjustment_comment

    solvent_summary = _solvent_bp_summary(detected_solvents)
    risk_comment = (
        f"cure_temperature {config.cure_temperature_K:.1f} K exceeds (bp - margin={config.cure_temp_safety_margin_K:.1f} K) "
        f"for {solvent_summary} at ~1 bar; may cavitate."
    )
    warning_comments.append(risk_comment)

    if config.allow_cure_above_solvent_bp:
        reason = (config.allow_cure_above_solvent_bp_reason or "").strip()
        if not reason:
            raise ValueError(
                "allow_cure_above_solvent_bp=True requires a non-empty "
                "allow_cure_above_solvent_bp_reason when boiling-risk is detected."
            )
        warning_comments.append(
            f"Explicit override active (allow_cure_above_solvent_bp=true). Reason: {reason}"
        )
        return effective_temp, warning_comments, adjustment_comment

    policy = config.cure_temp_policy
    if policy == "error":
        raise ValueError(
            f"cure_temperature_K={config.cure_temperature_K:.1f} K exceeds safe threshold {threshold_k:.1f} K "
            f"(min solvent bp {min_bp_k:.1f} K - margin {config.cure_temp_safety_margin_K:.1f} K; detected {solvent_summary}). "
            "Lower cure_temperature_K, set cure_temp_policy='cap'/'warn', or set "
            "allow_cure_above_solvent_bp=true with allow_cure_above_solvent_bp_reason."
        )
    if policy == "cap":
        effective_temp = min(config.cure_temperature_K, threshold_k)
        adjustment_comment = (
            f"Requested cure_temperature_K={config.cure_temperature_K:.1f} K capped to {effective_temp:.1f} K "
            "by cure_temp_policy=cap."
        )
        warning_comments.append(f"Applied cure_temp_policy=cap at threshold {threshold_k:.1f} K.")
    elif policy == "warn":
        warning_comments.append("Policy cure_temp_policy=warn keeps requested cure_temperature_K.")
    elif policy == "off":
        warning_comments.append("Policy cure_temp_policy=off disables automatic cap/error actions.")

    return effective_temp, warning_comments, adjustment_comment


def _resolve_conversion_guardrail(
    config: HTPolyNetConfig,
) -> Tuple[float, List[str], Optional[str]]:
    effective_max_conversion = config.max_conversion
    warning_comments: List[str] = []
    adjustment_comment: Optional[str] = None

    if config.max_conversion <= config.recommended_max_conversion:
        return effective_max_conversion, warning_comments, adjustment_comment

    risk_comment = (
        f"max_conversion {config.max_conversion:.2f} exceeds recommended {config.recommended_max_conversion:.2f}; "
        "late-stage linking can inject stress and trigger GROMACS LINCS instability."
    )
    warning_comments.append(risk_comment)

    if config.allow_high_conversion:
        reason = (config.allow_high_conversion_reason or "").strip()
        if not reason:
            raise ValueError(
                "allow_high_conversion=True requires a non-empty allow_high_conversion_reason when "
                "max_conversion exceeds recommended_max_conversion."
            )
        warning_comments.append(f"Explicit override active (allow_high_conversion=true). Reason: {reason}")
        return effective_max_conversion, warning_comments, adjustment_comment

    policy = config.conversion_policy
    if policy == "error":
        raise ValueError(
            f"max_conversion={config.max_conversion:.2f} exceeds recommended_max_conversion="
            f"{config.recommended_max_conversion:.2f}. "
            "Lower max_conversion, set conversion_policy='cap'/'warn', or set "
            "allow_high_conversion=true with allow_high_conversion_reason."
        )
    if policy == "cap":
        effective_max_conversion = config.recommended_max_conversion
        adjustment_comment = (
            f"Requested max_conversion={config.max_conversion:.2f} capped to "
            f"{effective_max_conversion:.2f} by conversion_policy=cap."
        )
        warning_comments.append(
            f"Applied conversion_policy=cap at recommended_max_conversion={config.recommended_max_conversion:.2f}."
        )
    elif policy == "warn":
        warning_comments.append("Policy conversion_policy=warn keeps requested max_conversion.")
    elif policy == "off":
        warning_comments.append("Policy conversion_policy=off disables automatic cap/error actions.")

    return effective_max_conversion, warning_comments, adjustment_comment


def _search_radius_comments(config: HTPolyNetConfig) -> Tuple[List[str], List[str]]:
    warning_comments: List[str] = []
    recipe_comments: List[str] = []

    if HTPOLYNET_SUPPORTS_SEARCH_RADIUS_SCHEDULE:
        return warning_comments, recipe_comments

    recipe_comments = [
        "Staged strategy is recommended for dense GPE systems (run sequentially):",
        "stage 1: until conversion 0.30 with search-radius 0.35 nm",
        "stage 2: until conversion 0.60 with search-radius 0.45 nm",
        "stage 3: until conversion 0.80 with search-radius 0.60 nm",
    ]

    if config.search_radius_schedule:
        warning_comments.append(
            "search_radius_schedule was provided but HTPolyNet config schema support is not confirmed in this repo; "
            "only single search-radius is emitted."
        )

    return warning_comments, recipe_comments


# =============================================================================
# Jinja2 Template
# =============================================================================

HTPOLYNET_CONFIG_TEMPLATE = """# HTPolyNet System Configuration
# Generated by MD Simulation Pipeline
# Project: {{ project_name }}
{% if warning_comments %}
# ========================================
# Safety Warnings
# ========================================
{% for warning in warning_comments %}
# WARNING: {{ warning }}
{% endfor %}
{% endif %}
Title: {{ project_name_q }}

# Forcefield and charge settings
forcefield: {{ forcefield_q }}

# ========================================
# System Composition
# ========================================
# Total molecules: {{ total_molecules }}
# Box size: {{ box_length_nm }} nm
{% for item in molecule_counts_items %}
# {{ item.name }}: {{ item.count }}
{% endfor %}

# ========================================
# Initial Structure
# ========================================
initial-structure: {{ initial_gro_q }}

# ========================================
# Box Specification (nm)
# ========================================
box:
  type: "orthorhombic"
  vectors: [{{ box_vectors_nm[0] }}, {{ box_vectors_nm[1] }}, {{ box_vectors_nm[2] }}]

# ========================================
# Reactive-Site Classification
# ========================================
{% if reactive_site_count_items %}
reactive-site-counts:
{% for item in reactive_site_count_items %}
  {{ item.name_q }}: {{ item.count }}
{% endfor %}
{% endif %}
{% if crosslinker_items %}
crosslinkers:
{% for item in crosslinker_items %}
  {{ item.name_q }}:
    count: {{ item.count }}
    reactive_site_count: {{ item.reactive_site_count }}
    classification: {{ item.classification_q }}
{% endfor %}
{% endif %}
{% if chain_stopper_items %}
chain-stoppers:
{% for item in chain_stopper_items %}
  {{ item.name_q }}:
    count: {{ item.count }}
    reactive_site_count: {{ item.reactive_site_count }}
    classification: "chain_stopper"
{% endfor %}
{% endif %}

# ========================================
# Reaction Rules
# ========================================
{% if rules_file_q %}
reactions-file: {{ rules_file_q }}
{% else %}
# WARNING: No rules file provided - using default radical polymerization
reactions:
  - name: "radical_addition"
    type: "C=C + C=C"
    probability: 1.0
{% endif %}

# ========================================
# Crosslinking Parameters
# ========================================
crosslink:
  max-conversion: {{ max_conversion_effective }}
{% if max_conversion_adjustment_comment %}
  # {{ max_conversion_adjustment_comment }}
{% endif %}
  search-radius: {{ search_radius_nm }}  # nm
{% if search_radius_recipe_comments %}
{% for line in search_radius_recipe_comments %}
  # {{ line }}
{% endfor %}
{% endif %}

# Annealing/curing
cure:
  initial-temperature: {{ initial_temperature_K }}  # K
  cure-temperature: {{ cure_temperature_K_effective }}  # K
{% if cure_temperature_adjustment_comment %}
  # {{ cure_temperature_adjustment_comment }}
{% endif %}
  final-temperature: {{ final_temperature_K }}  # K
{% if curing_steps %}

# Optional multi-step curing protocol
curing-steps:
{% for step in curing_steps %}
  - temperature_K: {{ step["temperature_K"] }}
    duration_ps: {{ step["duration_ps"] }}
{% if step["ramp_ps"] is not none %}
    ramp_ps: {{ step["ramp_ps"] }}
{% endif %}
{% endfor %}
{% endif %}

# ========================================
# GROMACS Settings
# ========================================
gromacs:
  maxwarn: {{ gmx_maxwarn }}

# ========================================
# Resource Settings
# ========================================
resources:
  threads: {{ threads }}

# ========================================
# Output Settings
# ========================================
output:
  prefix: "system"
  formats:
    - "gro"
    - "itp"
    - "top"
"""


# =============================================================================
# Config Generation
# =============================================================================


def generate_system_config(config: HTPolyNetConfig) -> str:
    """
    Generate HTPolyNet system_config.yaml content.

    Uses Jinja2 if available, otherwise falls back to semantically equivalent
    pure-Python rendering.

    Args:
        config: HTPolyNetConfig with all settings

    Returns:
        YAML configuration string
    """
    if config.box_vectors_nm is None:
        box_vectors_nm = (config.box_length_nm, config.box_length_nm, config.box_length_nm)
    else:
        box_vectors_nm = config.box_vectors_nm

    sorted_counts = _sorted_molecule_counts(config.molecule_counts)
    detected_solvents = _detect_known_solvents(sorted_counts)
    cure_temp_effective, cure_warnings, cure_adjustment_comment = _resolve_cure_temperature_guardrail(
        config, detected_solvents
    )
    max_conversion_effective, conversion_warnings, conversion_adjustment_comment = _resolve_conversion_guardrail(
        config
    )
    schedule_warnings, schedule_recipe_comments = _search_radius_comments(config)

    warning_comments = cure_warnings + conversion_warnings + schedule_warnings

    template_vars = {
        "project_name": config.project_name,
        "project_name_q": _yaml_quote(config.project_name),
        "forcefield_q": _yaml_quote(config.forcefield),
        "initial_gro_q": _yaml_quote(str(config.initial_gro)),
        "molecule_counts_items": _molecule_count_items(sorted_counts),
        "total_molecules": sum(sorted_counts.values()),
        "box_length_nm": config.box_length_nm,
        "box_vectors_nm": box_vectors_nm,
        "rules_file_q": _yaml_quote(str(config.rules_file)) if config.rules_file else None,
        "max_conversion_effective": max_conversion_effective,
        "max_conversion_adjustment_comment": conversion_adjustment_comment,
        "search_radius_nm": config.search_radius_nm,
        "initial_temperature_K": config.initial_temperature_K,
        "final_temperature_K": config.final_temperature_K,
        "cure_temperature_K_effective": cure_temp_effective,
        "cure_temperature_adjustment_comment": cure_adjustment_comment,
        "curing_steps": _normalize_curing_steps(config.curing_steps),
        "reactive_site_count_items": _reactive_site_items(config.reactive_site_counts),
        "crosslinker_items": _crosslinker_items(config.crosslinkers),
        "chain_stopper_items": _chain_stopper_items(config.chain_stoppers),
        "gmx_maxwarn": config.gmx_maxwarn,
        "threads": config.threads,
        "warning_comments": warning_comments,
        "search_radius_recipe_comments": schedule_recipe_comments,
    }

    if JINJA2_AVAILABLE and Template is not None:
        template = Template(HTPOLYNET_CONFIG_TEMPLATE)
        return template.render(**template_vars)
    return _simple_config_generator(template_vars)


def _simple_config_generator(vars: Dict[str, Any]) -> str:
    """Fallback config generator without Jinja2 (semantically equivalent sections)."""
    lines = [
        "# HTPolyNet System Configuration",
        "# Generated by MD Simulation Pipeline",
        f"# Project: {vars['project_name']}",
    ]

    if vars["warning_comments"]:
        lines.extend(
            [
                "# ========================================",
                "# Safety Warnings",
                "# ========================================",
            ]
        )
        for warning in vars["warning_comments"]:
            lines.append(f"# WARNING: {warning}")

    lines.extend(
        [
            "",
            f"Title: {vars['project_name_q']}",
            "",
            "# Forcefield and charge settings",
            f"forcefield: {vars['forcefield_q']}",
            "",
            "# ========================================",
            "# System Composition",
            "# ========================================",
            f"# Total molecules: {vars['total_molecules']}",
            f"# Box size: {vars['box_length_nm']} nm",
        ]
    )

    for item in vars["molecule_counts_items"]:
        lines.append(f"# {item['name']}: {item['count']}")

    lines.extend(
        [
            "",
            "# ========================================",
            "# Initial Structure",
            "# ========================================",
            f"initial-structure: {vars['initial_gro_q']}",
            "",
            "# ========================================",
            "# Box Specification (nm)",
            "# ========================================",
            "box:",
            "  type: \"orthorhombic\"",
            (
                f"  vectors: [{vars['box_vectors_nm'][0]}, {vars['box_vectors_nm'][1]}, "
                f"{vars['box_vectors_nm'][2]}]"
            ),
            "",
            "# ========================================",
            "# Reactive-Site Classification",
            "# ========================================",
        ]
    )

    if vars["reactive_site_count_items"]:
        lines.append("reactive-site-counts:")
        for item in vars["reactive_site_count_items"]:
            lines.append(f"  {item['name_q']}: {item['count']}")

    if vars["crosslinker_items"]:
        lines.append("crosslinkers:")
        for item in vars["crosslinker_items"]:
            lines.append(f"  {item['name_q']}:")
            lines.append(f"    count: {item['count']}")
            lines.append(f"    reactive_site_count: {item['reactive_site_count']}")
            lines.append(f"    classification: {item['classification_q']}")

    if vars["chain_stopper_items"]:
        lines.append("chain-stoppers:")
        for item in vars["chain_stopper_items"]:
            lines.append(f"  {item['name_q']}:")
            lines.append(f"    count: {item['count']}")
            lines.append(f"    reactive_site_count: {item['reactive_site_count']}")
            lines.append("    classification: \"chain_stopper\"")

    lines.extend(
        [
            "",
            "# ========================================",
            "# Reaction Rules",
            "# ========================================",
        ]
    )
    if vars["rules_file_q"]:
        lines.append(f"reactions-file: {vars['rules_file_q']}")
    else:
        lines.extend(
            [
                "# WARNING: No rules file provided - using default radical polymerization",
                "reactions:",
                "  - name: \"radical_addition\"",
                "    type: \"C=C + C=C\"",
                "    probability: 1.0",
            ]
        )

    lines.extend(
        [
            "",
            "# ========================================",
            "# Crosslinking Parameters",
            "# ========================================",
            "crosslink:",
            f"  max-conversion: {vars['max_conversion_effective']}",
        ]
    )
    if vars["max_conversion_adjustment_comment"]:
        lines.append(f"  # {vars['max_conversion_adjustment_comment']}")
    lines.append(f"  search-radius: {vars['search_radius_nm']}  # nm")
    for line in vars["search_radius_recipe_comments"]:
        lines.append(f"  # {line}")

    lines.extend(
        [
            "",
            "# Annealing/curing",
            "cure:",
            f"  initial-temperature: {vars['initial_temperature_K']}  # K",
            f"  cure-temperature: {vars['cure_temperature_K_effective']}  # K",
        ]
    )
    if vars["cure_temperature_adjustment_comment"]:
        lines.append(f"  # {vars['cure_temperature_adjustment_comment']}")
    lines.append(f"  final-temperature: {vars['final_temperature_K']}  # K")

    if vars["curing_steps"]:
        lines.extend(
            [
                "",
                "# Optional multi-step curing protocol",
                "curing-steps:",
            ]
        )
        for step in vars["curing_steps"]:
            lines.append(f"  - temperature_K: {step['temperature_K']}")
            lines.append(f"    duration_ps: {step['duration_ps']}")
            if step.get("ramp_ps") is not None:
                lines.append(f"    ramp_ps: {step['ramp_ps']}")

    lines.extend(
        [
            "",
            "# ========================================",
            "# GROMACS Settings",
            "# ========================================",
            "gromacs:",
            f"  maxwarn: {vars['gmx_maxwarn']}",
            "",
            "# ========================================",
            "# Resource Settings",
            "# ========================================",
            "resources:",
            f"  threads: {vars['threads']}",
            "",
            "# ========================================",
            "# Output Settings",
            "# ========================================",
            "output:",
            "  prefix: \"system\"",
            "  formats:",
            "    - \"gro\"",
            "    - \"itp\"",
            "    - \"top\"",
        ]
    )

    return "\n".join(lines)


def write_system_config(config: HTPolyNetConfig, output_path: Path) -> Path:
    """
    Generate and write system_config.yaml.

    Args:
        config: HTPolyNetConfig with all settings
        output_path: Path to write config file

    Returns:
        Path to written config file
    """
    content = generate_system_config(config)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(content)
    return output_path


# =============================================================================
# Rules Resolution
# =============================================================================


def _normalize_ff_dir_name(forcefield: str) -> str:
    ff_dir_name = forcefield.replace("-", "").upper()
    if ff_dir_name == "OPLSAA":
        return "OPLSAA"
    if ff_dir_name == "GAFF2":
        return "GAFF2"
    return ff_dir_name


def _resolve_rules_file_in_dir(rules_dir: Path) -> Path:
    candidates = [
        "reactions.yaml",
        "reactions.yml",
        "rules.yaml",
        "rules.yml",
    ]
    for candidate in candidates:
        rules_file = rules_dir / candidate
        if rules_file.exists():
            return rules_file
    yamls = list(rules_dir.glob("*.yaml")) + list(rules_dir.glob("*.yml"))
    if len(yamls) == 1:
        return yamls[0]
    if len(yamls) > 1:
        raise FileNotFoundError(
            f"Ambiguous rules: found {len(yamls)} YAML files in {rules_dir}. "
            "Please ensure exactly one rules file or name it 'reactions.yaml'."
        )
    raise FileNotFoundError(
        f"No rules YAML found in {rules_dir}. Expected reactions.yaml (or rules.yaml)."
    )


def _validate_rules_dir(
    rules_dir: Path,
    system_id: str,
    forcefield: str,
) -> None:
    if not rules_dir.exists():
        raise FileNotFoundError(
            f"HTPolyNet rules directory missing.\n"
            f"Expected: {rules_dir}\n"
            f"System: {system_id}, forcefield: {forcefield}\n"
            f"Create the directory and add a rules YAML (e.g., reactions.yaml) and "
            f"a forcefield marker (e.g., forcefield.itp)."
        )
    entries = list(rules_dir.glob("*"))
    if not entries:
        raise FileNotFoundError(
            f"HTPolyNet rules directory is empty: {rules_dir}\n"
            f"System: {system_id}, forcefield: {forcefield}\n"
            f"Add a rules YAML (e.g., reactions.yaml) and forcefield marker files."
        )
    ff_markers = [
        rules_dir / "forcefield.itp",
        rules_dir / "forcefield.top",
        rules_dir / "ff" / "forcefield.itp",
        rules_dir / "ff" / "forcefield.top",
    ]
    if not any(p.exists() for p in ff_markers):
        raise FileNotFoundError(
            f"Missing forcefield markers in {rules_dir}.\n"
            f"System: {system_id}, forcefield: {forcefield}\n"
            f"Expected one of: forcefield.itp, forcefield.top (or under ff/)."
        )


def resolve_rules_bundle(
    system_id: str,
    forcefield: str,
    in_dir: Path,
) -> Tuple[Path, Path, str]:
    """
    Resolve the rules directory and rules file for a given forcefield.

    Search order (first match wins):
      1) IN/systems/<SYSTEM_ID>/htpolynet/rules/<FF>/
      2) IN/htpolynet/rules/<FF>/
      3) IN/forcefield/<FF>/htpolynet/rules/<FF>/
      4) IN/forcefield/<FF>/rules/<FF>/

    For robustness, each tier checks canonical + lowercase variants for the FF
    trailing rules subdirectory, and forcefield/<ff> base directory supports
    both lowercase and uppercase variants.
    """
    ff_dir_name = _normalize_ff_dir_name(forcefield)
    rules_tail_variants = _ordered_unique([ff_dir_name, ff_dir_name.lower(), ff_dir_name.upper()])
    ff_base_variants = _ordered_unique([ff_dir_name.lower(), ff_dir_name.upper(), ff_dir_name])

    search_dirs: List[Tuple[Path, str]] = []
    seen_paths: set = set()

    def _add_candidate(path: Path, source: str) -> None:
        path_key = path.as_posix()
        if path_key in seen_paths:
            return
        seen_paths.add(path_key)
        search_dirs.append((path, source))

    # Tier 1: system-local
    for tail in rules_tail_variants:
        _add_candidate(in_dir / "systems" / system_id / "htpolynet" / "rules" / tail, "system-local")

    # Tier 2: shared
    for tail in rules_tail_variants:
        _add_candidate(in_dir / "htpolynet" / "rules" / tail, "shared")

    # Tier 3: forcefield-htpolynet
    for ff_base in ff_base_variants:
        for tail in rules_tail_variants:
            _add_candidate(
                in_dir / "forcefield" / ff_base / "htpolynet" / "rules" / tail,
                "forcefield-htpolynet",
            )

    # Tier 4: forcefield-shared
    for ff_base in ff_base_variants:
        for tail in rules_tail_variants:
            _add_candidate(
                in_dir / "forcefield" / ff_base / "rules" / tail,
                "forcefield-shared",
            )

    for rules_dir, source in search_dirs:
        if rules_dir.exists():
            _validate_rules_dir(rules_dir, system_id, forcefield)
            rules_file = _resolve_rules_file_in_dir(rules_dir)
            return rules_dir, rules_file, source

    expected = "\n  - ".join(str(p[0]) for p in search_dirs)
    raise FileNotFoundError(
        f"No HTPolyNet rules directory found for forcefield {forcefield}.\n"
        f"System: {system_id}\n"
        f"Searched:\n  - {expected}\n"
        f"Create one of these directories and add reactions.yaml + forcefield.itp."
    )
