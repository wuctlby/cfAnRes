"""Canvas factory in the plot_syst.py / plot_2pc_sp.py style.

get_canvas_temp(template, cname, ...) returns

  'main'       -> TCanvas                     (pad: canvas.GetPad(1))
  'ratio_lr'   -> TCanvas                     (pads: canvas.GetPad(1) = main,
  'ratio_tb'                                   canvas.GetPad(2) = ratio)
  'main_ratio' -> (main_canvas, ratio_canvas) (two separate figures)

The subpads are real TCanvas subpads, so the standard ROOT addressing works:
canvas.cd(1) / canvas.cd(2) / canvas.GetPad(1) / canvas.GetPad(2).

Notes:
- ratio_lr's right pad is a full-height pad and uses main-pad-like fonts;
  the ratio_tb strip is short and needs larger relative text (see _style_axis).
"""

import ROOT

# ---- global style, reused from plot_syst.py / plot_2pc_sp.py ----
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)

# Canvas sizes (px) from the two style sources.
_DEFAULT_CANVAS_SIZE = (1920, 1536)   # main canvas (plot_syst.py)
_RATIO_TB_CANVAS_SIZE = (1200, 1400)  # stacked main+ratio canvas (plot_2pc_sp.py)
_RATIO_CANVAS_SIZE = (1200, 1000)     # standalone ratio canvas (plot_2pc_sp.py c2)

# (left, right, bottom, top) margins relative to each pad.
_MARGIN_MAIN = (0.12, 0.03, 0.12, 0.05)       # plot_syst.py main canvas
_MARGIN_TB_MAIN = (0.15, 0.05, 0.0, 0.05)     # plot_2pc_sp.py pad1 (top, joined)
_MARGIN_TB_RATIO = (0.15, 0.05, 0.3, 0.0)     # plot_2pc_sp.py pad2 (bottom, joined)
_MARGIN_SINGLE = (0.15, 0.05, 0.12, 0.05)     # plot_2pc_sp.py "no ratio" canvas
_MARGIN_LR_MAIN = (0.15, 0.02, 0.12, 0.05)    # ratio_lr left pad
_MARGIN_LR_RATIO = (0.16, 0.045, 0.12, 0.05)  # ratio_lr right pad (full height)

# Pad positions (NDC) for the split layouts: (main, ratio).
_LR_NDC = ((0.0, 0.0, 0.49, 1.0), (0.5, 0.0, 1.0, 1.0))    # main | ratio
_TB_NDC = ((0.0, 0.35, 1.0, 1.0), (0.0, 0.0, 1.0, 0.35))   # main over ratio

# Color palette shared with plot_syst.py / plot_2pc_sp.py.
COLOR_CENTRAL = ROOT.TColor.GetColor("#D55E00")  # Vermillion
COLOR_LM = ROOT.TColor.GetColor("#56B4E9")       # Sky Blue
COLOR_FIT = ROOT.TColor.GetColor("#E69F00")      # Soft orange
COLOR_FD = ROOT.TColor.GetColor("#009E73")       # Bluish green
COLOR_OUTLINE = ROOT.kGray + 1
ALPHA_BAND = 0.3

# Standard ALICE axis labels (caller passes these into x_title / y_title).
X_TITLE_PT = "#it{p}_{T}^{D} (GeV/#it{c})"
Y_TITLE_V2 = "#it{v}_{2}^{prompt}"

# Per-call counter so repeated calls (even with the same cname) never produce
# colliding ROOT object names; the caller's Python variables are irrelevant.
_name_counter = 0


def _validate_range(value, name):
    """Range must be (min, max) with min < max, or None."""
    if value is None:
        return
    if (not isinstance(value, (tuple, list)) or len(value) != 2
            or not (value[0] < value[1])):
        raise ValueError(f"{name} must be (min, max) with min < max, got {value!r}")


def _style_pad(pad, margins):
    """Apply the notebook's compact margins to a pad."""
    left, right, bottom, top = margins
    pad.SetLeftMargin(left)
    pad.SetRightMargin(right)
    pad.SetBottomMargin(bottom)
    pad.SetTopMargin(top)


def _style_axis(frame, kind):
    """Set axis text sizes.

    kind:
      'main'       - large main pad
      'ratio_full' - full-size ratio pad (ratio_lr right pad, standalone ratio
                     canvas): same scale as the main pad
      'ratio_tb'   - short ratio strip of the stacked layout: larger relative
                     text so it stays readable
    """
    if kind == "main":
        frame.GetXaxis().SetTitleSize(0.045)
        frame.GetXaxis().SetLabelSize(0.04)
        frame.GetXaxis().SetTitleOffset(1.0)
        frame.GetYaxis().SetTitleSize(0.06)
        frame.GetYaxis().SetLabelSize(0.05)
        frame.GetYaxis().SetTitleOffset(1.0)
    elif kind == "ratio_full":
        frame.GetXaxis().SetTitleSize(0.045)
        frame.GetXaxis().SetLabelSize(0.04)
        frame.GetXaxis().SetTitleOffset(1.0)
        frame.GetYaxis().SetTitleSize(0.05)
        frame.GetYaxis().SetLabelSize(0.048)
        frame.GetYaxis().SetTitleOffset(1.2)
        frame.GetYaxis().SetNdivisions(505)
    else:  # ratio_tb
        for ax in (frame.GetXaxis(), frame.GetYaxis()):
            ax.SetTitleSize(0.12)
            ax.SetLabelSize(0.10)
        frame.GetYaxis().SetTitleOffset(0.55)
        frame.GetXaxis().SetTitleOffset(1.0)
        frame.GetYaxis().SetNdivisions(505)


def _make_frame(pad, x_title, y_title, x_range, y_range, kind="main"):
    """Draw an empty frame that owns the axes/titles/ranges.

    The caller later draws data with `Draw("same")` on the same pad.  A None
    range falls back to a (0,1) placeholder so the axis titles still land
    somewhere; the caller's histogram then re-defines the real range.
    """
    xmin, xmax = x_range if x_range is not None else (0.0, 1.0)
    ymin, ymax = y_range if y_range is not None else (0.0, 1.0)
    pad.cd()
    frame = pad.DrawFrame(xmin, ymin, xmax, ymax)
    frame.GetXaxis().SetTitle(x_title)
    frame.GetYaxis().SetTitle(y_title)
    _style_axis(frame, kind)
    return frame


def _new_split_canvas(cname, canvas_size, ndc_pair, margin_main, margin_ratio,
                      x_title, y_title, x_range, y_range,
                      ratio_x_title, ratio_y_title, ratio_x_range, ratio_y_range,
                      main_kind, ratio_kind):
    """Canvas with two subpads ('main' + 'ratio') addressable via cd/GetPad."""
    canvas = ROOT.TCanvas(cname, cname, canvas_size[0], canvas_size[1])
    canvas.Divide(2, 1)
    main_pad, ratio_pad = canvas.GetPad(1), canvas.GetPad(2)
    main_pad.SetPad(*ndc_pair[0])
    ratio_pad.SetPad(*ndc_pair[1])
    main_pad.SetName(f"{cname}_main")
    ratio_pad.SetName(f"{cname}_ratio")
    _style_pad(main_pad, margin_main)
    _style_pad(ratio_pad, margin_ratio)
    _make_frame(main_pad, x_title, y_title, x_range, y_range, kind=main_kind)
    _make_frame(ratio_pad, ratio_x_title, ratio_y_title, ratio_x_range, ratio_y_range, kind=ratio_kind)
    return canvas


def get_canvas_temp(
    template="main",
    cname=None,
    x_title="",
    y_title="",
    x_range=None,
    y_range=None,
    ratio_x_title="",
    ratio_y_title="",
    ratio_x_range=None,
    ratio_y_range=None,
    canvas_size=None,
):
    """Create a canvas in the plot_syst.py / plot_2pc_sp.py style.

    Returns
    -------
    'main', 'ratio_lr', 'ratio_tb' : the TCanvas.
        The content pads are addressed the standard ROOT way:
        ``canvas.cd(1)`` / ``canvas.GetPad(1)`` = main pad,
        ``canvas.cd(2)`` / ``canvas.GetPad(2)`` = ratio pad (split templates).
    'main_ratio' : tuple ``(main_canvas, ratio_canvas)`` (two separate figures);
        each canvas may be used directly (``c.cd()``) or via ``c.GetPad(1)``.

    Parameters
    ----------
    template : {"main", "main_ratio", "ratio_lr", "ratio_tb"}
        ``main`` is a single pad; ``main_ratio`` uses two independent canvases;
        ``ratio_lr`` / ``ratio_tb`` use one canvas with the main and ratio pads
        side-by-side or stacked.
    cname : str, optional
        Base name used to build unique ROOT object names.
    x_title, y_title : str
        Main pad axis titles.
    x_range, y_range : (min, max), optional
        Main pad axis ranges; None lets the later histogram decide.
    ratio_x_title, ratio_y_title : str
        Ratio pad axis titles.
    ratio_x_range, ratio_y_range : (min, max), optional
        Ratio pad ranges; ``ratio_x_range`` inherits ``x_range`` when unset.
    canvas_size : (width, height), optional
        Canvas pixel size; defaults to the template's standard size.
    """
    allowed = ("main", "main_ratio", "ratio_lr", "ratio_tb")
    if template not in allowed:
        raise ValueError(f"unknown template '{template}'; allowed: {', '.join(allowed)}")

    _validate_range(x_range, "x_range")
    _validate_range(y_range, "y_range")
    _validate_range(ratio_x_range, "ratio_x_range")
    _validate_range(ratio_y_range, "ratio_y_range")

    if canvas_size is None:
        canvas_size = _RATIO_TB_CANVAS_SIZE if template == "ratio_tb" else _DEFAULT_CANVAS_SIZE

    # ratio inherits the main x-range only when it sets none of its own
    if ratio_x_range is None:
        ratio_x_range = x_range

    # unique ROOT names: base name + per-call counter (never rely on Python refs)
    global _name_counter
    _name_counter += 1
    base = cname if cname else "canvas"
    uid = _name_counter
    main_canvas_name = f"{base}_{uid}"

    # ---- template dispatch ----
    if template == "main":
        canvas = ROOT.TCanvas(main_canvas_name, main_canvas_name, canvas_size[0], canvas_size[1])
        canvas.Divide(1, 1)
        pad = canvas.GetPad(1)
        pad.SetName(f"{main_canvas_name}_main")
        _style_pad(pad, _MARGIN_MAIN)
        _make_frame(pad, x_title, y_title, x_range, y_range, kind="main")
        return canvas

    if template == "ratio_lr":
        return _new_split_canvas(
            main_canvas_name, canvas_size, _LR_NDC,
            _MARGIN_LR_MAIN, _MARGIN_LR_RATIO,
            x_title, y_title, x_range, y_range,
            ratio_x_title, ratio_y_title, ratio_x_range, ratio_y_range,
            "main", "ratio_full")

    if template == "ratio_tb":
        return _new_split_canvas(
            main_canvas_name, canvas_size, _TB_NDC,
            _MARGIN_TB_MAIN, _MARGIN_TB_RATIO,
            x_title, y_title, x_range, y_range,
            ratio_x_title, ratio_y_title, ratio_x_range, ratio_y_range,
            "main", "ratio_tb")

    # template == "main_ratio": two separate canvases
    main_canvas = ROOT.TCanvas(main_canvas_name, main_canvas_name, canvas_size[0], canvas_size[1])
    main_canvas.Divide(1, 1)
    main_pad = main_canvas.GetPad(1)
    main_pad.SetName(f"{main_canvas_name}_main")
    _style_pad(main_pad, _MARGIN_MAIN)
    _make_frame(main_pad, x_title, y_title, x_range, y_range, kind="main")

    rw, rh = _RATIO_CANVAS_SIZE if canvas_size == _DEFAULT_CANVAS_SIZE else (canvas_size[0], canvas_size[1])
    ratio_canvas = ROOT.TCanvas(f"{base}_{uid}_ratio", f"{base}_{uid}_ratio", rw, rh)
    ratio_canvas.Divide(1, 1)
    ratio_pad = ratio_canvas.GetPad(1)
    ratio_pad.SetName(f"{base}_{uid}_ratio_main")
    _style_pad(ratio_pad, _MARGIN_SINGLE)
    _make_frame(ratio_pad, ratio_x_title, ratio_y_title, ratio_x_range, ratio_y_range, kind="ratio_full")

    return main_canvas, ratio_canvas
