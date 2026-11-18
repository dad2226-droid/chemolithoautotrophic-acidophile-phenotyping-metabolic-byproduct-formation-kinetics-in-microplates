"""Ferric-indicator stability and indicator response in acidophile media."""

import string

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.collections import LineCollection
from matplotlib.colors import Normalize
from matplotlib.lines import Line2D
from matplotlib.ticker import AutoMinorLocator, FixedLocator, FormatStrFormatter, LogLocator, NullFormatter, NullLocator
from scipy.stats import gaussian_kde


IRON = ("#D55E00", "#CC79A7", "#E69F00", "#882255")
MEDIA = {
    "AFM4": ("AFM4", "#0072B2", "o"),
    "AFM5": ("AFM5", "#009E73", "s"),
    "F2S": ("F2S", "#332288", "^"),
    "Starkey_brown": ("St-BS", "#56B4E9", "v"),
    "Starkey_colloidal": ("St-YS", "#882255", "D"),
}
DYES = ("no_thymol_blue", "thymol_blue")
SIGNALS = ("A380", "A380_corrected")
PLATES = ("UV_transparent", "non_UV_transparent")
ACIDS = (1., .1, .01, .001)
IRON_AXIS = r"measured $[\mathrm{Fe}]_{\mathrm{total}}$ (M)"
IRON_SIGNAL = r"$[\mathrm{Fe}^{3+}]_{\mathrm{app.}}$ indicator (AU)"
ABSORBANCE_LABELS = (r"$A_{380}$", r"$A^{*}_{380}$")
APPARENT_IRON = r"$\mathrm{Fe}^{3+}_{\mathrm{app.}}$"
CORRECTED_IRON = r"$\mathrm{Fe}^{3+}_{\mathrm{app.,\,corr.}}$"
BOTH_IRON = r"$\mathrm{Fe}^{3+}_{\mathrm{app.,\,(corr.)}}$"
FERRIC_CONCENTRATION = r"[Fe$^{3+}$]"
INCUBATION_TEMPERATURE = r"$T_{\mathrm{incubate}}$"
ACID_CONCENTRATION = r"[H$_2$SO$_4$]"


# Publication layout dimensions are in physical PDF points.
PANEL_PT = 97.2
SHORT_DASH = (0, (1.4, 1.1))


def _grid(rows):
    if rows == 3:
        width, height, left = 517.62, 492.12, 24.516
        tops = (12.379, 178.232, 358.413)
    else:
        width, height, left = 510.12, 127.68, 23.855
        tops = (12.379,)
    fig = plt.figure(figsize=(width / 72, height / 72))
    axes = []
    for top in tops:
        for column in range(4):
            axes.append(fig.add_axes(((left + 129.6 * column) / width,
                                      1 - (top + PANEL_PT) / height,
                                      PANEL_PT / width, PANEL_PT / height)))
    return fig, axes


def _panel(ax, letter, title):
    # Keep every title and panel letter aligned, including multiline handling titles.
    bottom = {"i": 11.7, "j": 11.7, "k": 11.7, "l": 2.5}.get(letter, 5.4)
    ax.set_title(title, fontweight="normal", fontsize=7, y=1, pad=bottom)
    ax.text(-24.1 / PANEL_PT, 1 + 2.8 / PANEL_PT, letter, transform=ax.transAxes,
            ha="left", va="baseline", fontsize=8, fontweight="bold")
    ax.tick_params(which="both", pad=2, direction="in", top=True, right=True, width=.5)
    ax.tick_params(which="major", labelsize=6, length=3)
    ax.tick_params(which="minor", length=1.6)
    ax.xaxis.labelpad = ax.yaxis.labelpad = 1
    ax.xaxis.label.set_size(7)
    ax.yaxis.label.set_size(7)
    for spine in ax.spines.values():
        spine.set_linewidth(.5)


def _key_text(ax, x, y, text, *, ha="left", rotation=0, va="center"):
    ax.text(x / PANEL_PT, 1 - y / PANEL_PT, text, transform=ax.transAxes,
            fontsize=5, ha=ha, va=va, rotation=rotation, zorder=12)


def _key_sample(ax, x, y, line, span=6.4):
    """Place a line and marker inside the key in print coordinates."""
    y = 1 - y / PANEL_PT
    opts = dict(transform=ax.transAxes, zorder=11, clip_on=False)
    path = Line2D([(x-span/2)/PANEL_PT, (x+span/2)/PANEL_PT], [y, y],
                  color=line.get_color(), lw=.8, ls=line.get_linestyle(), **opts)
    if line.is_dashed():
        path.set_linestyle(SHORT_DASH)
    ax.add_artist(path)
    ax.add_artist(Line2D([x/PANEL_PT], [y], color=line.get_color(),
                         marker=line.get_marker(), markerfacecolor=line.get_markerfacecolor(),
                         markeredgewidth=line.get_markeredgewidth(),
                         markersize=line.get_markersize(), ls="none", **opts))


def _matrix(ax, rows, columns, handles, *, x, y, header=None, side=None,
            side_x=None, vertical=False):
    """Draw a factor key at fixed column and row pitches in the plot frame."""
    centers = [x + 8.6 * j for j in range(len(columns))]
    column_y = y - (9.7 if vertical else 5.6)
    if header:
        _key_text(ax, np.mean(centers), column_y - 5.6, header, ha="center")
    for xc, label in zip(centers, columns):
        _key_text(ax, xc, column_y, label, ha="center", rotation=90 if vertical else 0)
    for i, (label, samples) in enumerate(zip(rows, handles)):
        yc = y + 5.6 * i
        _key_text(ax, centers[-1] + 4.6, yc, label)
        for xc, sample in zip(centers, samples):
            _key_sample(ax, xc, yc, sample)
    if side:
        _key_text(ax, side_x, (column_y + y + 5.6 * (len(rows)-1)) / 2,
                  side, ha="center", rotation=90)


def _symbol(label, color, marker="o", style="-", open_symbol=False):
    return Line2D([], [], label=label, color=color, lw=.8, marker=marker, ms=3,
                  ls=SHORT_DASH if style == "--" else style,
                  markerfacecolor="none" if open_symbol else color,
                  markeredgewidth=.45)


def _curve(ax, x, center, sd, color, marker="o", style="-", open_symbol=False, *, series, point_x, point_y, plates, sd_kind="sample SD"):
    x, center, sd = (np.asarray(v, dtype=float) for v in (x, center, sd))
    order = np.argsort(x)
    x, center, sd = x[order], center[order], sd[order]
    band = ax.fill_between(x, center - sd, center + sd, color=color, alpha=.18, linewidth=0)
    band.set_gid(f"{series} | mean ± {sd_kind}")
    line, = ax.plot(x, center, color=color, ls=SHORT_DASH if style == "--" else style,
                    lw=.8)
    points = ax.scatter(point_x, point_y, s=9, marker=marker, edgecolors=color, linewidths=.45,
                        facecolors="none" if open_symbol else color, zorder=3,
                        gid=series + " | physical plate")
    points.plotted_value_labels = list(plates)
    line.set_gid(series + " | mean")


def _plate_curve(ax, data, color, marker="o", style="-", open_symbol=False, *, series):
    data = data[data.predictor.gt(0)]
    group = data.groupby("level_id", as_index=False).agg(
        x=("predictor", "mean"), center=("response_mean", "mean"), sd=("response_mean", "std"))
    group = group[group.x.gt(0)]
    _curve(ax, group.x / 1000, group.center, group.sd, color, marker, style, open_symbol, series=series,
           point_x=data.predictor / 1000, point_y=data.response_mean, plates=data.physical_plate_id)


def _heatmap(ax, values, row_labels, col_labels, norm, rotation=0):
    palette = plt.get_cmap("BrBG").with_extremes(bad="#969696")
    x = np.arange(len(col_labels) + 1) - .5
    y = np.arange(len(row_labels) + 1) - .5
    image = ax.pcolormesh(x, y, np.ma.masked_invalid(values), cmap=palette, norm=norm,
                          edgecolors="black", linewidth=.25, antialiased=False)
    ax.set_xlim(x[0], x[-1])
    ax.set_ylim(y[-1], y[0])
    ax.set_yticks(np.arange(len(row_labels)), row_labels)
    ax.set_xticks(np.arange(len(col_labels)), col_labels, rotation=rotation)
    ax.tick_params(axis="both", labelsize=5, length=1.5, pad=1.5)
    ax.tick_params(axis="x", labelsize=6, pad=1)
    ax.xaxis.set_minor_locator(NullLocator())
    ax.yaxis.set_minor_locator(NullLocator())
    return image


def _fold_key(fig, ax, image, label, *, storage):
    # Space the full-width key below the rotated acquisition labels or factor labels.
    gap = 20.554 if storage else 16.934
    height = 1.4 * 72 / 25.4
    key = ax.inset_axes([0, -(gap+height) / PANEL_PT, 1, height/PANEL_PT], label="colorbar")
    spacing = 1 if storage else .5
    ticks = np.arange(-3, 4) * spacing
    colorbar = fig.colorbar(image, cax=key, orientation="horizontal", ticks=ticks)
    key.xaxis.set_major_formatter(FormatStrFormatter("%d" if storage else "%.1f"))
    key.xaxis.set_minor_locator(AutoMinorLocator(2))
    key.tick_params(which="both", direction="in", bottom=True, top=True, pad=1, labelsize=6)
    key.tick_params(which="major", length=1.5)
    key.tick_params(which="minor", length=1)
    colorbar.set_label(label, fontsize=7, labelpad=.5)
    colorbar.outline.set_linewidth(.5)


def _iron_level(value):
    if value == 0:
        return "0.000"
    decimals = max(0, 1 - int(np.floor(np.log10(value))))
    return f"{value:.{decimals}f}"


def _storage_map(fig, ax, data, dye, signal, norm):
    g = data[data.condition.eq(dye) & data.signal.eq(signal)]
    levels = g.groupby("column").predictor.mean().sort_values()
    days = g.groupby("read_index").days.mean().sort_index()
    cells = g.pivot_table(index="column", columns="read_index", values="log2_fold")
    values = cells.reindex(index=levels.index, columns=days.index).to_numpy()
    image = _heatmap(ax, values, [_iron_level(v) for v in levels],
                     [f"{v:.1f}" for v in days], norm, rotation=90)
    ax.set(xlabel="incubation time (d)", ylabel=r"initial $[\mathrm{Fe}]_{\mathrm{total}}$ (mM)")
    image.set_gid(f"{signal} | {dye} | benchtop non-UV microplates | mean log2 fold from time zero")
    formula = (r"$\log_2(A_{380,t}/A_{380,0})$" if signal == "A380"
               else r"$\log_2(A^{*}_{380,t}/A^{*}_{380,0})$")
    _fold_key(fig, ax, image, formula, storage=True)


def _edges(ax, tables, storage):
    keys = ["physical_plate_id", "condition", "signal"]
    spans = storage[keys + ["initial_response_range"]].drop_duplicates()
    source = tables["edge_replicates"].query("assay == 'iron'").merge(
        spans, on=keys, how="left", validate="many_to_one")
    source["percent"] = 100 * source.edge_minus_interior / source.initial_response_range.where(
        source.initial_response_range.gt(0))
    plates = source.groupby(keys + ["read_index"], as_index=False).agg(
        percent=("percent", "mean"), days=("days", "mean"))
    handles = []
    for j, signal in enumerate(SIGNALS):
        for k, dye in enumerate(DYES):
            subset = plates[plates.condition.eq(dye) & plates.signal.eq(signal)]
            means = subset.groupby("read_index").agg(t=("days", "mean"),
                                                      avg=("percent", "mean"), sd=("percent", "std"))
            color, marker = IRON[j], ("o", "s")[k]
            style = ("-", "--")[k]
            _curve(ax, means.t, means.avg, means.sd, color, marker, style, j == 1,
                   series=f"{signal} | {dye} | edge–interior drift (% of initial span)",
                   point_x=subset.days, point_y=subset.percent, plates=subset.physical_plate_id)
            handles.append(_symbol(f"{ABSORBANCE_LABELS[j]}, {10*k} µM", color, marker,
                                   style, open_symbol=j == 1))
    ax.axhline(0, color=".4", ls=":", lw=.5)
    ax.set(xlabel="time (d)",
           ylabel=r"$\frac{(A_{\mathrm{e},t}-A_{\mathrm{e},0})-(A_{\mathrm{i},t}-A_{\mathrm{i},0})}{\operatorname{span}(A_0)}$ (%)",
           xlim=(0, 16), xticks=[0, 4, 8, 12, 16], ylim=(-1.816, 1.1132),
           yticks=np.arange(-1.5, 1.1, .5))
    ax.xaxis.set_minor_locator(AutoMinorLocator(2))
    ax.yaxis.set_minor_locator(AutoMinorLocator(2))
    ax.yaxis.set_major_formatter(FormatStrFormatter("%.1f"))
    _matrix(ax, ABSORBANCE_LABELS, ("0", "10"), [handles[:2], handles[2:]],
            x=69.8, y=19.2, header="[dye] (µM)")


def _reversed_acquisitions(levels):
    result = set()
    measured = levels[levels.signal.eq("A380") & levels.step_order.between(1, 6)
                      & levels.technical_n.ge(2)].dropna(subset=["predictor", "response_mean"])
    group_keys = ["condition", "plate_type", "acid_M", "step_order", "physical_plate_id"]
    for (dye, material, acid, step, _), curve in measured.groupby(group_keys):
        if curve.predictor.nunique() > 1 and curve.response_mean.nunique() > 1:
            if curve.predictor.rank().corr(curve.response_mean.rank()) < 0:
                result.add((dye, material, acid, int(step)))
    return result


def _handling_map(fig, ax, data, dye, signal, norm, reversed_acquisitions):
    g = data[data.condition.eq(dye) & data.signal.eq(signal)]
    plate_acid = pd.MultiIndex.from_product([PLATES, ACIDS], names=["plate_type", "acid_M"])
    cells = g.pivot_table(index=["plate_type", "acid_M"], columns="step_order", values="log2_fold")
    values = cells.reindex(index=plate_acid, columns=range(1, 7)).to_numpy()
    rows = [rf"$10^{{{int(np.log10(a))}}}$" for _, a in plate_acid]
    image = _heatmap(ax, values, rows, ["23", "4", "30", "75", "−", "+"], norm)
    ax.tick_params(axis="x", pad=2)
    image.set_gid(f"{signal} | {dye} | mean log2 fold vs UV 23 °C reference")
    image._plot_data_y_labels = [f"{material}, {acid:g} M H2SO4" for material, acid in plate_acid]
    image._plot_data_x_labels = ["23 °C", "4 °C", "30 °C", "75 °C", "no membrane", "membrane"]
    ax.axvline(3.5, color="black", lw=1.5)
    ax.axhline(3.5, color="black", lw=1.5)
    ax.set_ylabel(r"$[\mathrm{H_2SO_4}]$ (M)")
    for ypos, label in ((.75, "UV +"), (.25, "UV −")):
        ax.text(1.02, ypos, label, fontsize=7, rotation=90, va="center", ha="left",
                transform=ax.transAxes)
    ax.text(1/3, -.09, INCUBATION_TEMPERATURE + " (°C)", transform=ax.transAxes, ha="center", va="top", fontsize=7)
    ax.text(5/6, -.09, "membrane", transform=ax.transAxes, ha="center", va="top", fontsize=7)
    if signal == "A380":
        for row, (material, acid) in enumerate(plate_acid):
            for step in range(1, 7):
                if (dye, material, acid, step) in reversed_acquisitions:
                    ax.text(step-1, row, "†", fontsize=7, fontweight="bold", ha="center", va="center")
        if dye == "no_thymol_blue":
            ax.text(0, 0, "*", fontsize=7, fontweight="bold", ha="center", va="center")
    formula = (r"$\log_2(A_{380}/A_{380,\mathrm{UV,ref}})$" if signal == "A380"
               else r"$\log_2(A^{*}_{380}/A_{380,\mathrm{UV,ref}})$")
    _fold_key(fig, ax, image, formula, storage=False)


def _log_iron_axis(ax, xlim, upper):
    ax.set_xscale("log")
    ax.tick_params(axis="x", pad=.5)
    ax.set(xlabel=IRON_AXIS, ylabel=IRON_SIGNAL, xlim=xlim, ylim=(0, upper),
           yticks=np.arange(0, int(upper)+1))
    ax.xaxis.set_major_locator(LogLocator(base=10, numticks=10))
    ax.xaxis.set_minor_locator(LogLocator(base=10, subs=10**(np.arange(1, 10)/10), numticks=100))
    ax.xaxis.set_minor_formatter(NullFormatter())
    ax.yaxis.set_major_formatter(FormatStrFormatter("%.1f"))
    ax.yaxis.set_minor_locator(AutoMinorLocator(2))


def figure_3(tables):
    """Render optical stability, spatial drift and handling comparisons."""
    fig, grid = _grid(3)
    axes = list(grid)
    storage_factors = f" = f({FERRIC_CONCENTRATION}, $t$)"
    handling_factors = f"f(UV, {ACID_CONCENTRATION}, {INCUBATION_TEMPERATURE}, membrane)"
    titles = [
        f"{BOTH_IRON} = f(lid, {FERRIC_CONCENTRATION})",
        APPARENT_IRON + storage_factors,
        CORRECTED_IRON + storage_factors,
        "+ dye " + APPARENT_IRON + storage_factors,
        "+ dye " + CORRECTED_IRON + storage_factors,
        f"{BOTH_IRON} = f(well position, $t$)",
        f"{APPARENT_IRON} = f(UV, {INCUBATION_TEMPERATURE}, {FERRIC_CONCENTRATION})",
        f"{APPARENT_IRON} = f(UV, {ACID_CONCENTRATION}, {FERRIC_CONCENTRATION})",
        f"{APPARENT_IRON} detectability =\n{handling_factors}",
        f"{CORRECTED_IRON} detectability =\n{handling_factors}",
        f"+ dye {APPARENT_IRON} detectability =\n{handling_factors}",
        f"+ dye {CORRECTED_IRON}\ndetectability = f(UV, {ACID_CONCENTRATION},\n{INCUBATION_TEMPERATURE}, membrane)",
    ]
    for letter, ax, label in zip(string.ascii_lowercase, axes, titles):
        _panel(ax, letter, label)

    levels = tables["pub_robustness_levels"].query("assay == 'iron'")
    lid = levels[levels.family.eq("lid_comparison") & levels.condition.eq("no_thymol_blue")]
    handles = []
    for j, signal in enumerate(SIGNALS):
        for k, state in enumerate(("lid_absent", "lid_present")):
            g = lid[lid.signal.eq(signal) & lid.reading_label.eq(state)]
            style, marker = ("-", "--")[k], ("o", "s")[k]
            _plate_curve(axes[0], g, IRON[j], marker, style, j == 1,
                         series=f"{signal} | {state} | UV microplate | 0 µM dye")
            handles.append(_symbol(f"{ABSORBANCE_LABELS[j]}, lid {('−','+')[k]}",
                                   IRON[j], marker, style, j == 1))
    _log_iron_axis(axes[0], (7.478e-5, .12439), 4)
    _matrix(axes[0], ABSORBANCE_LABELS, ("−", "+"), [handles[:2], handles[2:]],
            x=7.6, y=16.9, header="lid")

    storage = tables["storage_folds"].query("assay == 'iron'")
    handling = tables["iron_handling_folds"]
    storage_means = storage.groupby(["condition", "signal", "row_group", "column", "read_index"]).log2_fold.mean()
    handling_means = handling.groupby(["condition", "signal", "plate_type", "acid_M", "step_order"]).log2_fold.mean()
    storage_limit = storage_means.abs().max()
    handling_limit = handling_means.abs().max()
    storage_norm = Normalize(-storage_limit, storage_limit)
    handling_norm = Normalize(-handling_limit, handling_limit)
    reversed_cells = _reversed_acquisitions(levels)
    for n, (dye, signal) in enumerate((dye, signal) for dye in DYES for signal in SIGNALS):
        _storage_map(fig, axes[n+1], storage, dye, signal, storage_norm)
        _handling_map(fig, axes[n+8], handling, dye, signal, handling_norm, reversed_cells)
    _edges(axes[5], tables, storage)

    raw = levels[levels.signal.eq("A380") & levels.condition.eq("no_thymol_blue")]
    for ax, data, factor, factors, labels in (
        (axes[6], raw[raw.acid_M.eq(1) & raw.step_order.between(1, 4)], "step_order",
         (1, 2, 3, 4), ("23 °C", "4 °C", "30 °C", "75 °C")),
        (axes[7], raw[raw.step_order.eq(1)], "acid_M", tuple(reversed(ACIDS)),
         tuple(rf"$10^{{{int(np.log10(v))}}}$ M" for v in reversed(ACIDS)))):
        handles = []
        for k, (value, label) in enumerate(zip(factors, labels)):
            for j, material in enumerate(PLATES):
                g = data[data[factor].eq(value) & data.plate_type.eq(material)]
                marker, style = ("o", "s")[j], ("-", "--")[j]
                _plate_curve(ax, g, IRON[k], marker, style,
                             series=f"A380 | {material} | {label} | 0 µM dye")
                handles.append(_symbol(f"{label}, {('UV+','UV−')[j]}", IRON[k], marker, style))
        is_temp = factor == "step_order"
        _log_iron_axis(ax, (.00022491, 1.3560) if is_temp else (.00022257, 1.3565),
                       3.4041 if is_temp else 3.42)
        rows = [s.removesuffix(" °C").removesuffix(" M") for s in labels]
        # The UV-negative square precedes the UV-positive circle in each row.
        samples = [[handles[2*k+1], handles[2*k]] for k in range(len(rows))]
        _matrix(ax, rows, ("−", "+"), samples, x=6.25, y=16.5 if is_temp else 16.6,
                header="UV", side=(INCUBATION_TEMPERATURE+" (°C)" if is_temp
                                   else ACID_CONCENTRATION+" (M)"), side_x=28.7 if is_temp else 33.5)
    return fig


def _technical_standard_curves(data):
    """Combine technical-well counts, means and sample SDs across equivalent standards."""
    group = data.copy()
    group["weighted"] = group.technical_n * group.response_mean.fillna(0)
    pooled = group.groupby("predictor").agg(n=("technical_n", "sum"), total=("weighted", "sum"))
    pooled["center"] = pooled.total / pooled.n.where(pooled.n.gt(0))
    centers = group.predictor.map(pooled.center)
    within = (group.technical_n - 1).clip(lower=0) * group.technical_sd.fillna(0)**2
    group["squared_deviations"] = within + group.technical_n * (group.response_mean - centers).fillna(0)**2
    pooled["sd"] = np.sqrt(group.groupby("predictor").squared_deviations.sum() /
                           (pooled.n-1).where(pooled.n.gt(1)))
    return pooled


def _medium_curves(ax, data, corrected):
    signal = "ratio_corrected" if corrected else "ratio_raw"
    box = [.10, .30, .64, .49] if corrected else [.10, .43, .48, .40]
    inset = ax.inset_axes(box)
    for medium, (label, color, marker) in MEDIA.items():
        q = data[data.medium_id.eq(medium) & data.signal.eq(signal) & data.condition.eq("with_dye")]
        g = _technical_standard_curves(q)
        for axis in (ax, inset):
            _curve(axis, g.index, g.center, g.sd, color, marker, open_symbol=corrected,
                   series=f"{medium} | {signal} | with dye", point_x=q.predictor, point_y=q.response_mean,
                   plates=q.physical_plate_id, sd_kind="pooled technical SD")
    # Make the complete uncertainty envelope visible while retaining the labeled 0–60 scale.
    lower = np.floor((ax.dataLim.ymin - .3) * 10) / 10 if corrected and ax.dataLim.ymin < 0 else 0
    ax.set(xlabel="measured pH", ylabel=r"$R^{*}$" if corrected else r"$R$",
           xticks=[1, 2, 3, 4], xlim=(4, .5),
           ylim=(lower, 60) if corrected else (.5, 3.),
           yticks=np.arange(0, 61, 10) if corrected else np.arange(.5, 3.01, .5))
    ax.xaxis.set_major_formatter(FormatStrFormatter("%.1f"))
    inset.set(xticks=[1, 2, 3] if corrected else [1, 2, 3, 4],
              xlim=(3, .5) if corrected else (4, .5),
              yticks=[0, 1, 2] if corrected else [.6, 1., 1.4],
              ylim=(0, 2) if corrected else (.5, 1.4))
    inset.tick_params(which="both", labelsize=5, pad=.6, direction="in", top=True, right=True)
    inset.tick_params(which="major", length=1.5)
    inset.tick_params(which="minor", length=.8)
    for axis in (ax, inset):
        axis.xaxis.set_minor_locator(AutoMinorLocator(2))
        axis.yaxis.set_minor_locator(AutoMinorLocator(2))
    if corrected:
        ax.yaxis.set_minor_locator(FixedLocator(np.arange(5, 60, 10)))
    # A single row places each medium's name next to its own sample.
    for (name, color, marker), label_x in zip(MEDIA.values(), (9.6, 28.0, 46.5, 60.8, 79.0)):
        line = _symbol(name, color, marker, open_symbol=corrected)
        _key_sample(ax, label_x-2.6, 7., line, span=3.9)
        _key_text(ax, label_x, 7., name)


def _quality(ax, data):
    media = ["AFM1", "AFM3", "F2S", "Starkey_brown", "Starkey_colloidal", "AFM4", "AFM5"]
    for readout, marker in (("A380", "o"), ("ratio_raw", "s")):
        subset = data[data.reference_signal.eq(readout)]
        for j, medium in enumerate(media):
            group = subset[subset.medium_id.eq(medium)].sort_values("physical_plate_id")
            values = group.contrast_log2_Q_ratio.to_numpy(float)
            valid = np.isfinite(values)
            paired = medium == "F2S"
            center = j + ((-.18 if readout == "A380" else .18) if paired else 0)
            if not valid.any():
                continue
            color = IRON[j] if readout == "A380" else MEDIA[medium][1]
            offsets = np.resize([-.09, .09, 0], len(values)) * (.6 if paired else 1)
            dots = ax.scatter((center+offsets)[valid], values[valid], s=9,
                              color=color, marker=marker, linewidths=0, zorder=4)
            dots.set_gid(f"{medium} | {readout} | physical-microplate log2(Q*/Q)")
            dots.plotted_value_labels = group.physical_plate_id.to_numpy()[valid]
            values = values[valid]
            if len(values) < 2 or np.ptp(values) == 0:
                continue
            density = gaussian_kde(values, bw_method="scott")
            support = np.linspace(values.min(), values.max(), 256)
            scale = (.15 if paired else .32) / density(support).max()
            halfwidth = scale * density(support)
            shape = ax.fill_betweenx(support, center-halfwidth, center+halfwidth,
                                     color=color, linewidth=.55, alpha=.18, zorder=1)
            shape.set_gid(f"{medium} | {readout} | distribution of log2(Q*/Q)")
            for side in (-1, 1):
                ax.add_collection(LineCollection([np.column_stack((center+side*halfwidth, support))],
                                                  colors=color, linewidths=.55, zorder=2))
            for statistic, probabilities, width, style in (
                ("quartiles", [.25, .75], .55, ":"), ("median", [.5], .9, "-")):
                y = np.quantile(values, probabilities)
                widths = scale * density(y)
                marks = ax.hlines(y, center-widths, center+widths, color=color,
                                  linewidth=width, linestyles=style, zorder=3)
                marks.set_gid(f"{medium} | {readout} | {statistic} log2(Q*/Q)")
    for cut in (1.5, 2.5, 4.5, 5.5):
        ax.axvline(cut, color=".65", lw=.55)
    ax.axhline(0, color=".55", lw=.5, ls=":", zorder=0)
    ax.set_xticks(range(len(media)), [MEDIA.get(m,(m,))[0] for m in media], rotation=90)
    ax.set(xlim=(-.65, 6.65), ylim=(-1.277, 1.2341), yticks=np.arange(-1., 1.01, .5),
           ylabel=r"$\Delta\log_{2}(Q)$")
    for index, (label, marker) in enumerate(((ABSORBANCE_LABELS[0], "o"), (r"$R$", "s"))):
        y = 7.2 + 5.6 * index
        _key_sample(ax, 6.5, y, _symbol(label, ".2", marker, style="none"))
        _key_text(ax, 8.7, y, label)


def figure_4(tables):
    """Render medium-specific calibration and same-well correction estimates."""
    fig, grid = _grid(1)
    axes = list(grid)
    titles = (
        r"$\mathrm{pH}_{\mathrm{app.}}$ = f(media, pH)",
        r"$\mathrm{pH}_{\mathrm{app.,\,corr.}}$ = f(media, pH)",
        f"{BOTH_IRON} = f(media, {FERRIC_CONCENTRATION})",
        "detectability = f(media)",
    )
    for ax, letter, label in zip(axes, "abcd", titles):
        _panel(ax, letter, label)
    data = tables["pub_media_levels"]
    _medium_curves(axes[0], data, False)
    _medium_curves(axes[1], data, True)
    handles = []
    iron = data[data.assay.eq("iron") & data.condition.eq("no_dye") & data.predictor.gt(0)]
    for j, medium in enumerate(("AFM1", "AFM3", "F2S")):
        for i, signal in enumerate(SIGNALS):
            q = iron[iron.medium_id.eq(medium) & iron.signal.eq(signal)]
            curve = _technical_standard_curves(q)
            marker, style = ("o", "s")[i], ("-", "--")[i]
            _curve(axes[2], curve.index / 1000, curve.center, curve.sd, IRON[j], marker, style, i==1,
                   series=f"{medium} | {signal} | no dye", point_x=q.predictor / 1000, point_y=q.response_mean,
                   plates=q.physical_plate_id, sd_kind="pooled technical SD")
            handles.append(_symbol(f"{medium}, {ABSORBANCE_LABELS[i]}", IRON[j], marker, style, i==1))
    _log_iron_axis(axes[2], (.0013287, .12890), 3.2503)
    axes[2].set_xlabel(r"added $[\mathrm{Fe}]_{\mathrm{total}}$ (M)")
    _matrix(axes[2], ABSORBANCE_LABELS, ("AFM1", "AFM3", "F2S"),
            [handles[::2], handles[1::2]], x=9.4, y=24.1, vertical=True)
    _quality(axes[3], tables["media_correction_effects__physical_plate_support"])
    return fig
