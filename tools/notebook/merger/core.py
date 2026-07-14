"""
Core merging engine — wraps the existing collage_png_pages_to_single /
collage_pdf_pages_to_single from merge_pdf with config-driven orchestration.

Task types
----------
simple      Collect files → one collage image.
chunked     Collect files → chunk by N → one collage per chunk.
page_group  Multiple multi-page PDFs → group pages by index → one collage per group.
composite   N groups of files, aligned row-by-row → one collage per aligned row.
nested      Each subgroup collages independently → final collage of subgroup outputs.
"""

from __future__ import annotations

import re
import sys
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence

import natsort
from PIL import Image

# ── import existing merge functions ──────────────────────────────────────────
_TEMPLATE_DIR = Path(__file__).resolve().parent.parent.parent / "template"
sys.path.insert(0, str(_TEMPLATE_DIR))
from merge_pdf import (                                # noqa: E402
    collage_pdf_pages_to_single,
    collage_png_pages_to_single,
)


# ═══════════════════════════════════════════════════════════════════════════════
# helpers
# ═══════════════════════════════════════════════════════════════════════════════

def _resolve_sort(method: str, files: List[Path]) -> List[Path]:
    if method == "natsort" or method == "natural":
        return natsort.natsorted(files, key=lambda p: p.name)
    if method == "name":
        return sorted(files, key=lambda p: p.name)
    if method == "mtime":
        return sorted(files, key=lambda p: p.stat().st_mtime)
    if method == "none":
        return files
    raise ValueError(f"Unknown sort method: {method}")


def _render_template(template: str, **kwargs) -> str:
    """Render a filename template with variables like {idx}, {start}, {end} etc.

    Supports Python format-spec:  {idx:03d}, {idx:04d}, etc.
    """
    # Build a custom dict that handles format specs
    class _Fmt:
        def __init__(self, value):
            self.value = value

        def __format__(self, fmt):
            if fmt:
                return format(self.value, fmt)
            return str(self.value)

    fmt_kwargs = {k: _Fmt(v) for k, v in kwargs.items()}
    return template.format(**fmt_kwargs)


def _ensure_dir(p: Path) -> Path:
    p.parent.mkdir(parents=True, exist_ok=True)
    return p


# ═══════════════════════════════════════════════════════════════════════════════
# file collection
# ═══════════════════════════════════════════════════════════════════════════════

def collect_files(
    globs: List[str],
    exclude: Optional[str] = None,
    sort: str = "natsort",
    input_type: str = "png",
) -> List[Path]:
    """Resolve glob patterns, optionally exclude, sort, deduplicate.

    Parameters
    ----------
    globs : list[str]
        Glob patterns (e.g. ``["/path/**/*.png"]``).
    exclude : str | None
        Substring to exclude from results (e.g. ``"merged"``).
    sort : str
        ``natsort`` | ``name`` | ``mtime`` | ``none``.
    input_type : str
        File extension hint; ``png`` or ``pdf``.  Only used for validation
        feedback — the glob pattern itself controls which files are collected.

    Returns
    -------
    list[Path]
        Sorted, deduplicated list of absolute paths.
    """
    seen: set[str] = set()
    files: list[Path] = []

    for pattern in globs:
        for p in sorted(Path().glob(pattern) if not pattern.startswith("/")
                        else Path("/").glob(pattern.lstrip("/"))):
            # Path.glob with absolute path needs special handling
            pass

    # Use a two-step approach: expand user, then glob
    for pattern in globs:
        expanded = Path(pattern).expanduser()
        # If pattern is absolute, use Path(pattern).glob doesn't work directly.
        # Use Path(expanded).parent.glob(expanded.name) for exact file patterns,
        # or use the recursive glob from root.
        if expanded.is_absolute():
            # Use pathlib's rglob from root
            root = Path(expanded.anchor)
            rel = str(expanded).lstrip(str(root))
            for p in root.glob(rel):
                if p.is_file():
                    if p.as_posix() not in seen:
                        seen.add(p.as_posix())
                        files.append(p)
        else:
            for p in Path().glob(str(expanded)):
                if p.is_file():
                    if p.as_posix() not in seen:
                        seen.add(p.as_posix())
                        files.append(p)

    # Filter by type extension
    if input_type == "png":
        files = [f for f in files if f.suffix.lower() in (".png",)]
    elif input_type == "pdf":
        files = [f for f in files if f.suffix.lower() in (".pdf",)]

    # Exclude
    if exclude:
        files = [f for f in files if exclude not in f.name]

    # Sort
    files = _resolve_sort(sort, files)

    return files


# ═══════════════════════════════════════════════════════════════════════════════
# collage parameters
# ═══════════════════════════════════════════════════════════════════════════════

@dataclass
class CollageParams:
    page_width: float = 1920
    page_height: float = 1080
    cols: Optional[int] = None
    margin: float = 2.0


def _extract_canvas(cfg: dict, key: str = "canvas") -> CollageParams:
    """Pull CollageParams from a config dict section.

    ``cols`` and ``margins`` can appear at task level OR inside the
    ``canvas`` dict (task-level wins).
    """
    c = cfg.get(key, {})
    return CollageParams(
        page_width=c.get("width", 1920),
        page_height=c.get("height", 1080),
        cols=cfg.get("cols", c.get("cols", None)),
        margin=cfg.get("margins", c.get("margin", 2.0)),
    )


# ═══════════════════════════════════════════════════════════════════════════════
# task runners
# ═══════════════════════════════════════════════════════════════════════════════

def run_simple(task: dict, defaults: dict, work_dir: Optional[Path] = None) -> List[Path]:
    """Simple: collect → one collage.

    Config keys
    -----------
    input.globs, input.exclude, input.sort, input.type
    output.path
    canvas.width, canvas.height
    cols, margins
    """
    cfg = task
    inp = cfg["input"]
    out = cfg["output"]
    cp = _extract_canvas(cfg)

    files = collect_files(
        globs=inp["globs"],
        exclude=inp.get("exclude", defaults.get("exclude")),
        sort=inp.get("sort", defaults.get("sort", "natsort")),
        input_type=inp.get("type", defaults.get("input_type", "png")),
    )

    if not files:
        raise FileNotFoundError(f"[{cfg.get('name', '?')}] No files matched input globs: {inp['globs']}")

    output_path = Path(out["path"]).expanduser()
    if work_dir and not output_path.is_absolute():
        output_path = work_dir / output_path
    _ensure_dir(output_path)

    print(f"  [{cfg.get('name', 'simple')}] {len(files)} files → {output_path}")
    collage_png_pages_to_single(
        input_png=[str(f) for f in files],
        output_png=str(output_path),
        page_width=cp.page_width,
        page_height=cp.page_height,
        cols=cp.cols,
        margin=cp.margin,
    )
    print(f"  ✓ {output_path}")
    return [output_path]


def run_chunked(task: dict, defaults: dict, work_dir: Optional[Path] = None) -> List[Path]:
    """Chunked: collect → split into chunks → one collage per chunk.

    Config keys
    -----------
    input.globs, input.exclude, input.sort, input.type
    chunk_size : int
    output.dir, output.name_template
    canvas, cols, margins
    """
    cfg = task
    inp = cfg["input"]
    out = cfg["output"]
    cp = _extract_canvas(cfg)
    chunk_size = cfg["chunk_size"]

    files = collect_files(
        globs=inp["globs"],
        exclude=inp.get("exclude", defaults.get("exclude")),
        sort=inp.get("sort", defaults.get("sort", "natsort")),
        input_type=inp.get("type", defaults.get("input_type", "png")),
    )

    if not files:
        raise FileNotFoundError(f"[{cfg.get('name', '?')}] No files matched input globs: {inp['globs']}")

    output_dir = Path(out["dir"]).expanduser()
    if work_dir and not output_dir.is_absolute():
        output_dir = work_dir / output_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    name_template = out.get("name_template", "collage_{idx:03d}.png")
    output_paths: list[Path] = []

    total = (len(files) + chunk_size - 1) // chunk_size
    for i in range(0, len(files), chunk_size):
        chunk = files[i : i + chunk_size]
        idx = i // chunk_size
        start_name = chunk[0].stem
        end_name = chunk[-1].stem
        fname = _render_template(
            name_template,
            idx=idx,
            start=start_name,
            end=end_name,
            total=total,
        )
        output_path = output_dir / fname
        _ensure_dir(output_path)

        print(f"  [{cfg.get('name', 'chunked')}] chunk {idx+1}/{total} ({len(chunk)} files) → {output_path.name}")
        collage_png_pages_to_single(
            input_png=[str(f) for f in chunk],
            output_png=str(output_path),
            page_width=cp.page_width,
            page_height=cp.page_height,
            cols=cp.cols,
            margin=cp.margin,
        )
        print(f"  ✓ {output_path}")
        output_paths.append(output_path)

    return output_paths


def run_page_group(task: dict, defaults: dict, work_dir: Optional[Path] = None) -> List[Path]:
    """Page-group: multiple multi-page PDFs → group pages by index → collage each.

    Config keys
    -----------
    input.pdfs : list[str]
    input.dpi : int
    output.dir, output.name_template
    canvas, cols, margins
    keep_temp : bool
    """
    cfg = task
    inp = cfg["input"]
    out = cfg["output"]
    cp = _extract_canvas(cfg)
    dpi = inp.get("dpi", defaults.get("dpi", 400))
    keep_temp = cfg.get("keep_temp", False)

    pdf_paths = [Path(p).expanduser() for p in inp["pdfs"]]

    # Validate
    for p in pdf_paths:
        if not p.exists():
            raise FileNotFoundError(f"[{cfg.get('name', '?')}] PDF not found: {p}")

    output_dir = Path(out["dir"]).expanduser()
    if work_dir and not output_dir.is_absolute():
        output_dir = work_dir / output_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    name_template = out.get("name_template", "page_{idx:03d}.png")
    temp_dir = output_dir / "_temp_pages"

    # ── extract all pages ──
    import fitz  # PyMuPDF

    all_pages: list[list[Image.Image]] = []
    max_pages = 0

    for pdf_path in pdf_paths:
        print(f"  读取: {pdf_path.name}")
        doc = fitz.open(str(pdf_path))
        pages: list[Image.Image] = []
        for page_idx in range(len(doc)):
            mat = fitz.Matrix(dpi / 72, dpi / 72)
            pix = doc[page_idx].get_pixmap(matrix=mat, alpha=False)
            img = Image.frombytes("RGB", (pix.width, pix.height), pix.samples)
            pages.append(img)
        all_pages.append(pages)
        max_pages = max(max_pages, len(pages))
        doc.close()
        print(f"    → {len(pages)} 页")

    # ── group by page index & collage ──
    output_paths: list[Path] = []

    for page_idx in range(max_pages):
        group_pngs: list[str] = []

        for pdf_idx, pages in enumerate(all_pages):
            if page_idx < len(pages):
                temp_path = temp_dir / f"pdf{pdf_idx:02d}_pg{page_idx:03d}.png"
                temp_path.parent.mkdir(parents=True, exist_ok=True)
                pages[page_idx].save(temp_path)
                group_pngs.append(str(temp_path))

        if not group_pngs:
            print(f"  页码 {page_idx+1}/{max_pages}: 无页面，跳过")
            continue

        output_path = output_dir / _render_template(name_template, idx=page_idx, total=max_pages)
        _ensure_dir(output_path)

        print(f"  页码 {page_idx+1}/{max_pages}: {len(group_pngs)} 张图 → {output_path.name}")
        collage_png_pages_to_single(
            input_png=group_pngs,
            output_png=str(output_path),
            page_width=cp.page_width,
            page_height=cp.page_height,
            cols=cp.cols,
            margin=cp.margin,
        )
        print(f"  ✓ {output_path}")
        output_paths.append(output_path)

        # clean up
        if not keep_temp:
            for p in group_pngs:
                Path(p).unlink()

    # ── clean temp dir ──
    if not keep_temp and temp_dir.exists():
        try:
            temp_dir.rmdir()
        except OSError:
            pass

    return output_paths


def run_composite(task: dict, defaults: dict, work_dir: Optional[Path] = None) -> List[Path]:
    """Composite: N groups of files, aligned row-by-row → collage each row.

    Groups are sorted independently; files at the same index across groups
    are collaged into one output image.

    Config keys
    -----------
    groups : list[dict]
        Each group: {name, input: {globs, exclude, sort}}
    output.dir, output.name_template
    canvas.width, canvas.height
    cols : int  (should equal number of groups for side-by-side)
    margins
    align : "strict" (error if groups differ in length)
            | "shortest" (truncate to shortest group)
    """
    cfg = task
    groups_cfg: list[dict] = cfg["groups"]
    out = cfg["output"]
    cp = _extract_canvas(cfg)
    align = cfg.get("align", "strict")

    output_dir = Path(out["dir"]).expanduser()
    if work_dir and not output_dir.is_absolute():
        output_dir = work_dir / output_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    name_template = out.get("name_template", "composite_{idx:03d}.png")

    # ── collect files for each group ──
    all_group_files: list[list[Path]] = []
    for g in groups_cfg:
        inp = g.get("input", {})
        files = collect_files(
            globs=inp.get("globs", []),
            exclude=inp.get("exclude", defaults.get("exclude")),
            sort=inp.get("sort", defaults.get("sort", "natsort")),
            input_type=inp.get("type", defaults.get("input_type", "png")),
        )
        if not files:
            raise FileNotFoundError(
                f"[{cfg.get('name', '?')}] group '{g['name']}': no files matched"
            )
        all_group_files.append(files)
        print(f"  组 '{g['name']}': {len(files)} 个文件")

    # ── check lengths ──
    lengths = [len(f) for f in all_group_files]
    if align == "strict":
        if len(set(lengths)) != 1:
            msg = " / ".join(
                f"{g['name']}={len(f)}"
                for g, f in zip(groups_cfg, all_group_files)
            )
            raise ValueError(
                f"[{cfg.get('name', '?')}] groups have different file counts ({msg}). "
                f"Set align='shortest' to truncate."
            )
        n_rows = lengths[0]
    elif align == "shortest":
        n_rows = min(lengths)
        all_group_files = [f[:n_rows] for f in all_group_files]
    else:
        raise ValueError(f"Unknown align mode: {align}")

    # ── collage each row ──
    output_paths: list[Path] = []
    for row_idx in range(n_rows):
        row_files = [str(group[row_idx]) for group in all_group_files]
        output_path = output_dir / _render_template(name_template, idx=row_idx, total=n_rows)
        _ensure_dir(output_path)

        print(f"  行 {row_idx+1}/{n_rows}: {len(row_files)} 张图 → {output_path.name}")
        collage_png_pages_to_single(
            input_png=row_files,
            output_png=str(output_path),
            page_width=cp.page_width,
            page_height=cp.page_height,
            cols=cp.cols,
            margin=cp.margin,
        )
        print(f"  ✓ {output_path}")
        output_paths.append(output_path)

    return output_paths


def run_nested(task: dict, defaults: dict, work_dir: Optional[Path] = None) -> List[Path]:
    """Nested: each subgroup collages independently → final collage of results.

    Config keys
    -----------
    groups : list[dict]
        Each group is a mini-task spec: {name, input, chunk_size, canvas, cols, margins}
    output.dir, keep_intermediate : bool
    canvas.final, final_cols, final_margins
    """
    cfg = task
    groups_cfg: list[dict] = cfg["groups"]
    out = cfg["output"]

    output_dir = Path(out["dir"]).expanduser()
    if work_dir and not output_dir.is_absolute():
        output_dir = work_dir / output_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    keep_intermediate = out.get("keep_intermediate", False)
    intermediate_dir = output_dir / "_intermediate" if not keep_intermediate else output_dir
    intermediate_dir.mkdir(parents=True, exist_ok=True)

    # ── extract canvas params ──
    final_canvas = cfg.get("canvas", {})
    final_cp = CollageParams(
        page_width=final_canvas.get("final", {}).get("width", 1920),
        page_height=final_canvas.get("final", {}).get("height", 1080),
        cols=cfg.get("final_cols", None),
        margin=cfg.get("margins", {}).get("final", cfg.get("margins", 2.0)),
    )

    # ── run each subgroup ──
    subgroup_outputs: list[list[Path]] = []
    total_subgroups = len(groups_cfg)

    for g_idx, g in enumerate(groups_cfg):
        print(f"\n  ── 子组 {g_idx+1}/{total_subgroups}: '{g['name']}' ──")

        # Build a mini chunked task
        g_inp = g.get("input", {})
        g_chunk_size = g.get("chunk_size", None)
        g_canvas = g.get("canvas", {})
        g_cp = CollageParams(
            page_width=g_canvas.get("width", 1920),
            page_height=g_canvas.get("height", 1080),
            cols=g.get("cols", None),
            margin=g.get("margins", cfg.get("margins", {}).get("group", 2.0)),
        )

        files = collect_files(
            globs=g_inp.get("globs", []),
            exclude=g_inp.get("exclude", defaults.get("exclude")),
            sort=g_inp.get("sort", defaults.get("sort", "natsort")),
            input_type=g_inp.get("type", defaults.get("input_type", "png")),
        )

        if not files:
            raise FileNotFoundError(
                f"[{cfg.get('name', '?')}] subgroup '{g['name']}': no files"
            )

        if g_chunk_size:
            # Chunked mode within subgroup
            name_tmpl = g.get("name_template", f"{g['name']}_{{:03d}}.png")
            subgroup_outputs.append([])
            for i in range(0, len(files), g_chunk_size):
                chunk = files[i : i + g_chunk_size]
                chunk_idx = i // g_chunk_size
                op = intermediate_dir / _render_template(name_tmpl, idx=chunk_idx)
                _ensure_dir(op)

                collage_png_pages_to_single(
                    input_png=[str(f) for f in chunk],
                    output_png=str(op),
                    page_width=g_cp.page_width,
                    page_height=g_cp.page_height,
                    cols=g_cp.cols,
                    margin=g_cp.margin,
                )
                subgroup_outputs[g_idx].append(op)
        else:
            # Simple mode within subgroup — all files in one collage
            op = intermediate_dir / f"{g['name']}.png"
            _ensure_dir(op)
            collage_png_pages_to_single(
                input_png=[str(f) for f in files],
                output_png=str(op),
                page_width=g_cp.page_width,
                page_height=g_cp.page_height,
                cols=g_cp.cols,
                margin=g_cp.margin,
            )
            subgroup_outputs.append([op])

        print(f"    → {len(subgroup_outputs[g_idx])} 输出")

    # ── all subgroups must produce same number of outputs ──
    sub_counts = [len(s) for s in subgroup_outputs]
    if len(set(sub_counts)) != 1:
        msg = " / ".join(
            f"{g['name']}={len(o)}" for g, o in zip(groups_cfg, subgroup_outputs)
        )
        raise ValueError(
            f"[{cfg.get('name', '?')}] subgroups produced different output counts ({msg}). "
            f"Ensure all subgroups have the same number of output images."
        )

    n_final = sub_counts[0] if sub_counts else 0
    name_template = out.get("name_template", "final_{idx:03d}.png")

    # ── final collage: for each index, take that output from each subgroup ──
    final_outputs: list[Path] = []
    for idx in range(n_final):
        row_files = [str(subgroup_outputs[s][idx]) for s in range(total_subgroups)]
        final_path = output_dir / _render_template(name_template, idx=idx, total=n_final)
        _ensure_dir(final_path)

        print(f"\n  最终合并 {idx+1}/{n_final}: {len(row_files)} 子组图 → {final_path.name}")
        collage_png_pages_to_single(
            input_png=row_files,
            output_png=str(final_path),
            page_width=final_cp.page_width,
            page_height=final_cp.page_height,
            cols=final_cp.cols,
            margin=final_cp.margin,
        )
        print(f"  ✓ {final_path}")
        final_outputs.append(final_path)

    # ── cleanup intermediate ──
    if not keep_intermediate and intermediate_dir != output_dir:
        import shutil
        shutil.rmtree(intermediate_dir, ignore_errors=True)

    return final_outputs


# ═══════════════════════════════════════════════════════════════════════════════
# runner dispatch
# ═══════════════════════════════════════════════════════════════════════════════

_RUNNERS = {
    "simple": run_simple,
    "chunked": run_chunked,
    "page_group": run_page_group,
    "composite": run_composite,
    "nested": run_nested,
}


def run_task(task: dict, defaults: dict, work_dir: Optional[Path] = None) -> List[Path]:
    """Execute a single task config.

    Returns list of output Paths produced.
    """
    task_type = task["type"]
    name = task.get("name", task_type)

    if task_type not in _RUNNERS:
        raise ValueError(
            f"Unknown task type '{task_type}' for task '{name}'. "
            f"Valid types: {list(_RUNNERS)}"
        )

    print(f"\n{'=' * 60}")
    print(f"Task: {name}  (type={task_type})")
    if task.get("desc"):
        print(f"  {task['desc']}")
    print(f"{'=' * 60}")

    return _RUNNERS[task_type](task, defaults, work_dir=work_dir)
