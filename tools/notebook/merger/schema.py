"""
Config loading, validation, and defaults merging for the merger tool.

Usage
-----
    from schema import load_config

    config = load_config("config.yaml")
    # config["tasks"] -> list of validated task dicts
    # config["defaults"] -> merged global defaults
"""

from __future__ import annotations

import copy
from pathlib import Path
from typing import Any, Dict, List, Optional

# Try PyYAML first, fall back to standard json if yaml not available
try:
    import yaml

    _YAML_AVAILABLE = True
except ImportError:
    import json as _json

    _YAML_AVAILABLE = False


# ═══════════════════════════════════════════════════════════════════════════════
# defaults
# ═══════════════════════════════════════════════════════════════════════════════

DEFAULT_CONFIG: Dict[str, Any] = {
    "input_type": "png",
    "sort": "natsort",
    "exclude": "merged",
    "dpi": 400,
    "canvas": {
        "width": 1920,
        "height": 1080,
    },
    "margins": 2.0,
    "cols": None,
}


# ═══════════════════════════════════════════════════════════════════════════════
# deep merge
# ═══════════════════════════════════════════════════════════════════════════════

def _deep_merge(base: dict, override: dict) -> dict:
    """Recursively merge *override* into *base*.  Lists are replaced, not merged."""
    result = copy.deepcopy(base)
    for key, val in override.items():
        if key in result and isinstance(result[key], dict) and isinstance(val, dict):
            result[key] = _deep_merge(result[key], val)
        else:
            result[key] = copy.deepcopy(val)
    return result


# ═══════════════════════════════════════════════════════════════════════════════
# load / parse
# ═══════════════════════════════════════════════════════════════════════════════

def load_config(path: str) -> dict:
    """Load and validate a YAML config file.

    Returns
    -------
    dict
        ``{"tasks": [...], "defaults": {...}}``
    """
    config_path = Path(path).expanduser().resolve()
    if not config_path.exists():
        raise FileNotFoundError(f"Config file not found: {config_path}")

    raw = _parse_file(config_path)

    return _process(raw, config_path)


def _parse_file(path: Path) -> dict:
    """Parse YAML or JSON config file."""
    content = path.read_text(encoding="utf-8")

    if _YAML_AVAILABLE:
        try:
            data = yaml.safe_load(content)
        except yaml.YAMLError as e:
            raise ValueError(f"Failed to parse YAML config {path}: {e}")
    else:
        try:
            data = _json.loads(content)
        except _json.JSONDecodeError as e:
            raise ValueError(
                f"Failed to parse config {path}. "
                f"Install PyYAML (`pip install pyyaml`) for YAML support, "
                f"or use valid JSON.  Error: {e}"
            )

    if data is None:
        raise ValueError(f"Config file is empty: {path}")

    if not isinstance(data, dict):
        raise ValueError(f"Config must be a mapping, got {type(data).__name__}: {path}")

    return data


def _process(raw: dict, config_path: Path) -> dict:
    """Merge defaults, normalise tasks, validate."""
    # ── defaults ──
    user_defaults = raw.get("defaults", {})
    global_defaults = _deep_merge(DEFAULT_CONFIG, user_defaults)

    # ── tasks ──
    tasks_raw = raw.get("tasks")

    if tasks_raw is None:
        raise ValueError(f"Config must have a 'tasks' key: {config_path}")

    if isinstance(tasks_raw, dict):
        # Named tasks: {name: config} → list with name embedded
        tasks = []
        for name, task_cfg in tasks_raw.items():
            t = copy.deepcopy(task_cfg)
            t["name"] = name
            tasks.append(t)
    elif isinstance(tasks_raw, list):
        tasks = copy.deepcopy(tasks_raw)
        # Ensure each task has a name
        for i, t in enumerate(tasks):
            if "name" not in t:
                t["name"] = f"task_{i}"
    else:
        raise ValueError(
            f"'tasks' must be a dict or list, got {type(tasks_raw).__name__}: {config_path}"
        )

    if not tasks:
        raise ValueError(f"No tasks defined in config: {config_path}")

    # ── validate each task ──
    validated = []
    for t in tasks:
        _validate_task(t, global_defaults, config_path)
        validated.append(t)

    return {
        "tasks": validated,
        "defaults": global_defaults,
        "_config_path": str(config_path),
    }


# ═══════════════════════════════════════════════════════════════════════════════
# validation
# ═══════════════════════════════════════════════════════════════════════════════

_VALID_TYPES = {"simple", "chunked", "page_group", "composite", "nested"}


def _validate_task(task: dict, defaults: dict, config_path: Path) -> None:
    """Validate one task dict in-place, raising ValueError on issues."""
    name = task.get("name", "?")

    # ── type ──
    task_type = task.get("type")
    if not task_type:
        raise ValueError(f"[{name}] Missing required key 'type': {config_path}")
    if task_type not in _VALID_TYPES:
        raise ValueError(
            f"[{name}] Unknown type '{task_type}'. Valid: {sorted(_VALID_TYPES)}"
        )

    # ── type-specific validation ──
    _validators = {
        "simple": _validate_simple,
        "chunked": _validate_chunked,
        "page_group": _validate_page_group,
        "composite": _validate_composite,
        "nested": _validate_nested,
    }
    _validators[task_type](task, defaults, config_path)


def _validate_input_section(task: dict, defaults: dict, config_path: Path) -> None:
    """Validate the 'input' section common to simple / chunked / composite groups."""
    name = task.get("name", "?")
    inp = task.get("input")

    if inp is None:
        # page_group uses 'input.pdfs'; composite uses 'groups[].input'
        return

    if not isinstance(inp, dict):
        raise ValueError(f"[{name}] 'input' must be a dict: {config_path}")

    globs = inp.get("globs")
    if not globs:
        raise ValueError(f"[{name}] 'input.globs' is required: {config_path}")
    if not isinstance(globs, list):
        raise ValueError(f"[{name}] 'input.globs' must be a list of glob patterns: {config_path}")


def _validate_simple(task: dict, defaults: dict, config_path: Path) -> None:
    name = task.get("name", "?")
    _validate_input_section(task, defaults, config_path)
    out = task.get("output")
    if not out or not out.get("path"):
        raise ValueError(f"[{name}] (simple) 'output.path' is required: {config_path}")


def _validate_chunked(task: dict, defaults: dict, config_path: Path) -> None:
    name = task.get("name", "?")
    _validate_input_section(task, defaults, config_path)
    if not task.get("chunk_size"):
        raise ValueError(f"[{name}] (chunked) 'chunk_size' is required: {config_path}")
    out = task.get("output")
    if not out or not out.get("dir"):
        raise ValueError(f"[{name}] (chunked) 'output.dir' is required: {config_path}")


def _validate_page_group(task: dict, defaults: dict, config_path: Path) -> None:
    name = task.get("name", "?")
    inp = task.get("input")
    if not inp or not inp.get("pdfs"):
        raise ValueError(f"[{name}] (page_group) 'input.pdfs' is required: {config_path}")
    out = task.get("output")
    if not out or not out.get("dir"):
        raise ValueError(f"[{name}] (page_group) 'output.dir' is required: {config_path}")


def _validate_composite(task: dict, defaults: dict, config_path: Path) -> None:
    name = task.get("name", "?")
    groups = task.get("groups")
    if not groups:
        raise ValueError(f"[{name}] (composite) 'groups' is required: {config_path}")
    if not isinstance(groups, list) or len(groups) < 2:
        raise ValueError(f"[{name}] (composite) 'groups' must be a list with ≥2 groups: {config_path}")

    for g in groups:
        if not isinstance(g, dict):
            raise ValueError(f"[{name}] (composite) each group must be a dict: {config_path}")
        if "name" not in g:
            raise ValueError(f"[{name}] (composite) each group needs a 'name': {config_path}")
        inp = g.get("input")
        if not inp or not inp.get("globs"):
            raise ValueError(
                f"[{name}] (composite) group '{g['name']}' missing 'input.globs': {config_path}"
            )

    out = task.get("output")
    if not out or not out.get("dir"):
        raise ValueError(f"[{name}] (composite) 'output.dir' is required: {config_path}")


def _validate_nested(task: dict, defaults: dict, config_path: Path) -> None:
    name = task.get("name", "?")
    groups = task.get("groups")
    if not groups:
        raise ValueError(f"[{name}] (nested) 'groups' is required: {config_path}")
    if not isinstance(groups, list) or len(groups) < 2:
        raise ValueError(f"[{name}] (nested) 'groups' must be a list with ≥2 groups: {config_path}")

    for g in groups:
        if not isinstance(g, dict):
            raise ValueError(f"[{name}] (nested) each group must be a dict: {config_path}")
        if "name" not in g:
            raise ValueError(f"[{name}] (nested) each group needs a 'name': {config_path}")
        inp = g.get("input")
        if not inp or not inp.get("globs"):
            raise ValueError(
                f"[{name}] (nested) group '{g['name']}' missing 'input.globs': {config_path}"
            )

    out = task.get("output")
    if not out or not out.get("dir"):
        raise ValueError(f"[{name}] (nested) 'output.dir' is required: {config_path}")


# ═══════════════════════════════════════════════════════════════════════════════
# dry-run: list what would happen
# ═══════════════════════════════════════════════════════════════════════════════

def dry_run(config: dict) -> str:
    """Return a human-readable summary of what the config will do (no file I/O)."""
    from core import collect_files  # lazy import to avoid circular

    lines: list[str] = []
    defaults = config.get("defaults", {})

    lines.append("=" * 60)
    lines.append("DRY RUN — 以下是将执行的操作（不实际生成文件）")
    lines.append("=" * 60)

    for task in config.get("tasks", []):
        name = task.get("name", "?")
        task_type = task["type"]
        desc = task.get("desc", "")

        lines.append(f"\n── Task: {name} (type={task_type})")
        if desc:
            lines.append(f"   {desc}")

        try:
            if task_type == "simple":
                inp = task["input"]
                files = collect_files(
                    globs=inp["globs"],
                    exclude=inp.get("exclude", defaults.get("exclude")),
                    sort=inp.get("sort", defaults.get("sort", "natsort")),
                    input_type=inp.get("type", defaults.get("input_type", "png")),
                )
                out_path = task["output"]["path"]
                lines.append(f"   输入: {len(files)} 个文件")
                lines.append(f"   输出: {out_path}")
                if files:
                    lines.append(f"   示例: {files[0].name}")
                else:
                    lines.append(f"   ⚠ 警告: 没有匹配到任何文件!")

            elif task_type == "chunked":
                inp = task["input"]
                files = collect_files(
                    globs=inp["globs"],
                    exclude=inp.get("exclude", defaults.get("exclude")),
                    sort=inp.get("sort", defaults.get("sort", "natsort")),
                    input_type=inp.get("type", defaults.get("input_type", "png")),
                )
                chunk_size = task["chunk_size"]
                n_chunks = (len(files) + chunk_size - 1) // chunk_size
                lines.append(f"   输入: {len(files)} 个文件")
                lines.append(f"   分块: 每 {chunk_size} 张一组 → {n_chunks} 个输出")
                lines.append(f"   输出目录: {task['output']['dir']}")

            elif task_type == "page_group":
                pdfs = task["input"]["pdfs"]
                lines.append(f"   输入 PDF: {len(pdfs)} 个")
                for p in pdfs:
                    pp = Path(p)
                    lines.append(f"     - {pp.name}" + (" (文件存在)" if pp.exists() else " ⚠ 文件不存在!"))
                lines.append(f"   输出目录: {task['output']['dir']}")

            elif task_type == "composite":
                groups = task["groups"]
                align = task.get("align", "strict")
                lines.append(f"   对齐方式: {align}")
                lens = []
                for g in groups:
                    inp = g.get("input", {})
                    files = collect_files(
                        globs=inp.get("globs", []),
                        exclude=inp.get("exclude", defaults.get("exclude")),
                        sort=inp.get("sort", defaults.get("sort", "natsort")),
                        input_type=inp.get("type", defaults.get("input_type", "png")),
                    )
                    lens.append(len(files))
                    lines.append(f"   组 '{g['name']}': {len(files)} 个文件")
                if align == "shortest":
                    n_rows = min(lens)
                    lines.append(f"   有效行数: {n_rows} (按最短截断)")
                else:
                    n_rows = lens[0] if lens else 0
                    lines.append(f"   行数: {n_rows}")
                lines.append(f"   输出目录: {task['output']['dir']}")

            elif task_type == "nested":
                groups = task["groups"]
                lines.append(f"   子组数: {len(groups)}")
                for g in groups:
                    inp = g.get("input", {})
                    files = collect_files(
                        globs=inp.get("globs", []),
                        exclude=inp.get("exclude", defaults.get("exclude")),
                        sort=inp.get("sort", defaults.get("sort", "natsort")),
                        input_type=inp.get("type", defaults.get("input_type", "png")),
                    )
                    chunk = g.get("chunk_size", None)
                    if chunk:
                        n = (len(files) + chunk - 1) // chunk
                        lines.append(f"   子组 '{g['name']}': {len(files)} 个文件 → {n} 个中间图")
                    else:
                        lines.append(f"   子组 '{g['name']}': {len(files)} 个文件 → 1 个中间图")
                lines.append(f"   最终输出目录: {task['output']['dir']}")

        except Exception as e:
            lines.append(f"   ⚠ dry-run 出错: {e}")

    lines.append("\n" + "=" * 60)
    return "\n".join(lines)
