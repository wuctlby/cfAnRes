from __future__ import annotations

import os
import re
from dataclasses import dataclass, field

import yaml


# ── common physics aliases for THnSparse axes ─────────────────────────────────
_COMMON_AXIS_ALIASES: dict[str, list[str]] = {
    'inv. mass':
        ['mass', 'm', 'invmass', 'invariant mass', 'inv_mass', 'minv'],

    'centrality':
        ['cent', 'c', 'central'],

    'pt':
        ['pt', 'p_t', 'transverse momentum', 'mom', 'p_t'],

    'sign':
        ['charge', 'q'],

    'bkg score':
        ['bkg', 'background', 'bkgscore'],

    'fd score':
        ['fd', 'fdscore', 'fd_score'],

    'eta':
        ['pseudorapidity'],

    'mean pt a':
        ['mpt_a', 'meanpt_a', 'ptbar_a', '<pt>_a'],

    'mean pt b':
        ['mpt_b', 'meanpt_b', 'ptbar_b', '<pt>_b'],

    'pt product':
        ['product', 'cb_product', 'ptcand_mpt'],

    'mean pt product':
        ['meanpt_product', 'mm_product'],

    'n tracks a':
        ['ntrk_a', 'ntracks_a', 'mult_a', 'n_a', 'num_a'],

    'n tracks b':
        ['ntrk_b', 'ntracks_b', 'mult_b', 'n_b', 'num_b'],
}

def _norm(s: str) -> str:
    """Normalise punctuation for fuzzy matching.
    """
    s = s.lower()
    s = re.sub(r'[_\-./()\[\],]', ' ', s)
    s = re.sub(r'\s+', ' ', s)
    return s.strip()

_ALIAS_LOOKUP: list[tuple[str, list[str]]] = [
    (_norm(key), aliases) for key, aliases in _COMMON_AXIS_ALIASES.items()
]


@dataclass
class THnSparseInfo:
    """Metadata container for a THnSparse histogram.

    Attributes:
        name:         Logical identifier (e.g. ``'CharmBulk'``).
        aliname:      Alias name for this method.
        cand:         Particle species (e.g. ``'D0'``), empty for generic.
        thname:       Actual THnSparse histogram name in the ROOT file.
        thpath:       Path to the ROOT file containing the THnSparse.
        thurl:        Full URL / task path inside the ROOT file.
        datatype:     Type of data (e.g. ``'DATA'``, ``'MC'``).
        axisnames:    Labels of each axis (in dimension order).
        axisids:      IDs of each axis (as defined in META.yaml, in dimension order).
                      After ``pre=True``, these become sequential positions ``[0, 1, 2, …]``.
                      Use :meth:`pre_axis_id` to get the original IDs.
        axisids_kept: Original axis IDs kept after pre-filtering.
        pre_thname:   THnSparse name after pre-preprocessing.
        pre_thpath:   THnSparse path after pre-preprocessing.
        ali:          Per-axis extra aliases, e.g. ``{2: ['pt', 'p_t']}``.
    """
    name: str
    aliname: str
    cand: list[str]
    thname: str
    thpath: str
    thurl: str
    datatype: str
    axisnames: list[str] = field(default_factory=list)
    axisids: list[int] = field(default_factory=list)
    axisids_kept: list[int] = field(default_factory=list)
    pre_thname: str = field(default_factory=str)
    pre_thpath: str = field(default_factory=str)
    ali: dict[int, list[str]] = field(default_factory=dict)

    def __post_init__(self) -> None:
        """Build the alias to ID lookup table and initialise pre‑axes mapping."""
        self._alias_to_id: dict[str, int] = {}
        claimed: set[str] = set()

        for i in self.axisids:
            axname = self.axisnames[i]

            # 1) Lowercased canonical name
            key_lower = axname.lower()
            self._alias_to_id[key_lower] = i

            # 2) Normalised form
            key_norm = _norm(axname)
            self._alias_to_id[key_norm] = i

            # 3) Common physics aliases
            for alikey_norm, aliases in _ALIAS_LOOKUP:
                if alikey_norm in claimed:
                    continue
                if key_norm.startswith(alikey_norm):
                    claimed.add(alikey_norm)
                    for alias in aliases:
                        self._alias_to_id[alias.lower()] = i
                        self._alias_to_id[_norm(alias)] = i

            # 4) Per-instance custom aliases
            if i in self.ali:
                for alias in self.ali[i]:
                    self._alias_to_id[alias.lower()] = i
                    self._alias_to_id[_norm(alias)] = i

        # ── pre‑axes mapping (default: identity) ──
        self._pre_axisnames_map: dict[int, int] = {i: i for i in self.axisids}
        self._full_axisnames: list[str] = self.axisnames
        self._full_axisids: list[int] = self.axisids

    def axis_id(self, name: str) -> int | None:
        """Return the axis index (position in current axes list) for *name* (or alias).

        Returns *None* if the name is not found.
        """
        key = name.lower()
        if key in self._alias_to_id:
            return self._alias_to_id[key]
        norm = _norm(name)
        if norm in self._alias_to_id:
            return self._alias_to_id[norm]
        return None

    def axis_name(self, idx: int) -> str:
        """Return the canonical axis label for the given 0-based *idx*."""
        return self.axisnames[idx]

    # ── pre‑axes mapping ────────────────────────────────────────────────────

    def pre_axis_id(self, name: str) -> int:
        """Return the **original** axis index from the full axis set."""
        seq_idx = self.axis_id(name)
        if seq_idx is None:
            raise KeyError(f"Axis '{name}' not found")
        return self._pre_axisnames_map.get(seq_idx, seq_idx)

    def pre_axis_name(self, idx: int) -> str:
        """Return the **original** axis name from the full axis set."""
        orig_id = self._pre_axisnames_map.get(idx, idx)
        if orig_id < len(self._full_axisnames):
            return self._full_axisnames[orig_id]
        return self.axisnames[idx]

    # ── name → id dictionary ───────────────────────────────────────────────

    @property
    def axis_name_id_map(self) -> dict[str, int]:
        """``{axis_name: axis_id}`` for the current axes list.

        If ``pre=True`` was used, the IDs are the new sequential positions.
        Use :meth:`pre_axis_id` to get the original pre‑pre IDs.
        """
        return dict(zip(self.axisnames, self.axisids))


# ── META.yaml loading ─────────────────────────────────────────────────────────

def _meta_path() -> str:
    """Return the path to ``META.yaml``, co-located with this file."""
    return os.path.join(os.path.dirname(os.path.abspath(__file__)), 'META.yaml')


def _load_meta(meta_file: str | None = None) -> dict:
    """Load and return the META.yaml dictionary."""
    path = meta_file or _meta_path()
    with open(path) as f:
        return yaml.safe_load(f)


def _build_thinfo(method_name: str, thn_def: dict, pre: bool = False) -> THnSparseInfo:
    """Construct a :class:`THnSparseInfo` from a single META.yaml method entry."""
    axesdef = thn_def['axes']
    thpath = thn_def.get('thpath', '')
    thname = thn_def['thname']
    datatype = thn_def.get('datatype', 'DATA')
    aliname = thn_def.get('aliname', '')
    cand = _ensure_cand_list(thn_def.get('cand', []))

    full_axisnames = [ax['name'] for ax in axesdef]
    full_axisids = [ax['id'] for ax in axesdef]

    pre_axisids = thn_def.get('pre_axisids', None)
    pre_thname = thn_def.get('pre_thname', None)
    pre_thpath = thn_def.get('pre_thpath', None)

    if pre and pre_thname is not None:
        thname = pre_thname

    if pre and pre_thpath is not None:
        thpath = pre_thpath

    if pre and pre_axisids is not None:
        ax_by_id = {ax['id']: ax for ax in axesdef}
        ordered_axes = [ax_by_id[aid] for aid in pre_axisids]
        axisnames = [ax['name'] for ax in ordered_axes]
        axisids = list(range(len(ordered_axes)))
        _pre_axisnames_map: dict[int, int] = {
            new_idx: orig_id for new_idx, orig_id in enumerate(pre_axisids)
        }
    else:
        axisnames = full_axisnames
        axisids = full_axisids
        _pre_axisnames_map = {i: i for i in range(len(axisnames))}

    info = THnSparseInfo(
        name=thn_def.get('name', method_name),
        aliname=aliname,
        cand=cand,
        thname=thname,
        thpath=thpath,
        thurl=f'{thpath}/{thname}',
        axisnames=axisnames,
        axisids=axisids,
        datatype=datatype,
        axisids_kept=pre_axisids,
        pre_thname=pre_thname or '',
        pre_thpath=pre_thpath or '',
    )
    info._pre_axisnames_map = _pre_axisnames_map
    info._full_axisnames = full_axisnames
    info._full_axisids = full_axisids
    return info


# ── species alias normalisation ──────────────────────────────────────────────

_CAND_ALIASES: dict[str, str] = {
    # D⁰  — PDG 421
    'd0': 'D0',
    'dzero': 'D0',
    '421': 'D0',
    # D⁺  — PDG 411
    # 'dplus': 'D⁺',
    # 'd+': 'D⁺',
    # '411': 'D⁺',
    # Dₛ⁺ — PDG 431
    # 'ds': 'Dₛ⁺',
    # 'ds+': 'Dₛ⁺',
    # '431': 'Dₛ⁺',
}

def _norm_cand(s: str) -> str:
    """Normalise a species (cand) string for comparison.

    Handles case‑insensitive matching and PDG MC particle IDs::

        'D0', 'd0', 'Dzero', '421'  →  all match ``cand='D0'``
    """
    return _CAND_ALIASES.get(s.strip().lower(), s.strip().lower())


def _ensure_cand_list(raw: str | list[str]) -> list[str]:
    """Normalise the ``cand`` YAML field to a list of strings.

    Handles backward compatibility with the old single‑string format::

        ``cand: "D0"``  →  ``["D0"]``
        ``cand: []``    →  ``[]``
        ``cand: ["D0", "D+"]``  →  unchanged
    """
    if isinstance(raw, str):
        return [raw] if raw else []
    return list(raw)


class GetTHnInfo:
    """THnSparseInfo registry backed by :file:`META.yaml`.

    Usage::

        info = GetTHnInfo.default()
        info.axis_id('mass')   # 0

        info = GetTHnInfo.thn('bulk')
        info.axis_id('Centrality')   # 0

        print(GetTHnInfo.list_thns())   # ['charm_bulk', 'bulk', ...]
    """

    _meta: dict | None = None
    _meta_file: str | None = None
    _name_to_key: dict[str, str] | None = None

    @classmethod
    def _ensure_loaded(cls) -> None:
        """Lazy-load META.yaml if not already cached."""
        if cls._meta is None:
            cls._meta = _load_meta(cls._meta_file)
            alias_map: dict[str, str] = {}
            for key, thn_def in cls._meta.get('thns', {}).items():
                alias_map[key] = key
                for alias in thn_def.get('aliname', []):
                    alias_map[alias] = key
            cls._name_to_key = alias_map

    @classmethod
    def configure(cls, meta_file: str | None = None) -> None:
        """Point to a different META.yaml file (optional)."""
        cls._meta_file = meta_file
        cls._meta = None

    @classmethod
    def thn(cls, name: str, pre: bool = False, cand: str | None = None) -> THnSparseInfo:
        """Build a :class:`THnSparseInfo` for the named *method*.

        Parameters
        ----------
        name :
            Method key in the YAML ``thns:`` section, or any entry in the
            method's ``aliname`` list.
        pre :
            If *True*, only axes listed in ``pre_axisids`` are kept and
            re‑indexed sequentially.
        cand :
            If given (non‑empty), only return the method if one of its
            ``cand`` list entries matches (after normalisation).
            Supports various names for the same species:
            ``'D0'``, ``'d0'``, ``'Dzero'``, ``'421'`` are all equivalent.

        Raises ``KeyError`` if *name* cannot be resolved, or if *cand*
        is given but the method's ``cand`` doesn't match.
        """
        cls._ensure_loaded()
        assert cls._meta is not None
        assert cls._name_to_key is not None

        key = cls._name_to_key.get(name)
        if key is None:
            raise KeyError(
                f"Unknown THnSparse method {name!r}. "
                f"Available: {list(cls._name_to_key)}"
            )
        thn_def = cls._meta['thns'][key]

        if cand:
            method_cand_list = _ensure_cand_list(thn_def.get('cand', []))
            normed = _norm_cand(cand)
            if not any(_norm_cand(c) == normed for c in method_cand_list):
                raise KeyError(
                    f"THnSparse method {name!r} (→ '{key}') has cand={method_cand_list!r}, "
                    f"does not match requested cand={cand!r}"
                )

        return _build_thinfo(key, thn_def, pre=pre)

    @classmethod
    def list_thns(cls) -> list[str]:
        """Return the list of available method names."""
        cls._ensure_loaded()
        assert cls._meta is not None
        return list(cls._meta.get('thns', {}).keys())

    @classmethod
    def default(cls) -> THnSparseInfo:
        """Return the default THnSparseInfo."""
        cls._ensure_loaded()
        assert cls._meta is not None
        default_name = cls._meta.get('default', '')
        if not default_name:
            thns = cls.list_thns()
            if not thns:
                raise RuntimeError("META.yaml defines no methods.")
            default_name = thns[0]
        return cls.thn(default_name)

    @classmethod
    def charm_bulk(cls) -> THnSparseInfo:
        return cls.thn('charm_bulk')

    @classmethod
    def CharmBulk(cls) -> THnSparseInfo:
        return cls.thn('charm_bulk')

    @classmethod
    def bulk(cls) -> THnSparseInfo:
        return cls.thn('bulk')

    @classmethod
    def Bulk(cls) -> THnSparseInfo:
        return cls.thn('bulk')
