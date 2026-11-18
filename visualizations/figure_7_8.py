"""Plate–flask comparisons and soluble-sulfur validation (Figures 7 and 8)."""
from __future__ import annotations

import json

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.offsetbox import AnnotationBbox, DrawingArea
from matplotlib.text import Text
from matplotlib.ticker import (
    AutoMinorLocator, FixedLocator, FormatStrFormatter, FuncFormatter, NullFormatter, NullLocator,
)
import numpy as np
import pandas as pd


METAL = {"HH": "#D55E00", "In": "#CC79A7", "Ga": "#E69F00"}
TREATMENT = {"no_added_metal": "HH", "In_10mM": "In", "Ga_100mM": "Ga"}
CELLS = {"Ab": "#666666", "AF": "#0072B2", "AT": "#009E73"}
ADDITIVES = {
    ("Al", 0): (r"Al$^{3+}$", "#CC79A7", "o"),
    ("Cu", 0): (r"Cu$^{2+}$", "#009E73", "s"),
    ("Ga", 0): (r"Ga$^{3+}$", "#0072B2", "^"),
    ("NaCl", 0): ("NaCl", "#E69F00", "v"),
    ("VO", 0): (r"VO$^{2+}$", "#777777", "D"),
    ("citric_acid", 0): ("citric acid", "#D55E00", "P"),
    ("citric_acid", 1): ("citric acid + Ga", "#332288", "P"),
}



def _grid(rows, columns, titles):
    """Position square axes and their letter baselines in physical points."""
    size = 97.2
    if rows == 3:
        width, height, left, column_pitch = 385.32, 398.76, 23.855, 131.4
        row_tops = (12.38, 141.98, 285.50)
        title_pads = (5.18, 5.18, 5.18, 5.53, 5.52, 5.52, 2.84, 2.84, 3.84)
    else:
        width, height, left, column_pitch = 512.88, 272.04, 23.855, 129.6
        row_tops = (20.47, 158.79)
        title_pads = (4.27, 14.47, 6.35, 4.27, 2.84, 3.84, 2.84, 3.84)
    fig = plt.figure(figsize=(width / 72, height / 72))
    axes = []
    for index, title in enumerate(titles):
        row, col = divmod(index, columns)
        ax = fig.add_axes(((left + col * column_pitch) / width,
                           1 - (row_tops[row] + size) / height, size / width, size / height))
        ax.text(-24.09 / size, 1 + 2.84 / size, chr(97 + index), transform=ax.transAxes,
                weight="bold", fontsize=8, ha="left", va="baseline")
        ax.set_title(title, fontsize=7, pad=title_pads[index])
        if rows == 2 and index == 6:
            ax.title.set_linespacing(1.4)
        ax.tick_params(which="both", direction="in", top=True, right=True, width=.5)
        ax.tick_params(which="major", labelsize=6, length=3, pad=2.)
        ax.tick_params(which="minor", length=1.6)
        ax.xaxis.labelpad = ax.yaxis.labelpad = 2.
        ax.xaxis.set_minor_locator(AutoMinorLocator(2))
        ax.yaxis.set_minor_locator(AutoMinorLocator(2))
        axes.append(ax)
    return fig, axes


class _PanelKey:
    """Place marker, line and text samples in a transparent, point-sized legend."""
    def __init__(self, ax):
        self.area = DrawingArea(97.2, 97.2, clip=False)
        box = AnnotationBbox(self.area, (0, 1), xycoords='axes fraction',
                             box_alignment=(0, 1), pad=0, frameon=False)
        box.set_zorder(8)
        ax.add_artist(box)

    def label(self, x, baseline, text, *, ha='left', rotation=0):
        self.area.add_artist(Text(x, 97.2 - baseline, text, fontsize=5, ha=ha,
                                  va='baseline', rotation=rotation, rotation_mode='anchor'))

    def sample(self, x, y, color, marker, ls='-', *, hollow=False, width=5.6):
        height = 97.2 - y
        if ls != 'none':
            self.area.add_artist(Line2D([x - width/2, x + width/2], [height, height],
                                        color=color, linewidth=.65, linestyle=ls))
        self.area.add_artist(Line2D([x], [height], color=color, marker=marker,
                                    markersize=3, markeredgewidth=.4, linestyle='none',
                                    markerfacecolor='white' if hollow else color))


def _cell_matrix(ax, headers, rows, row_title, samples, *, x=0, y=7.69):
    key = _PanelKey(ax)
    for x_col, heading in zip((7.55, 16.15), headers):
        key.label(x + x_col, y, heading, ha='center')
    key.label(x + 20.8, y, row_title)
    for i, (label, row_samples) in enumerate(zip(rows, samples), 1):
        baseline = y + i * 5.6
        key.label(x + 20.8, baseline, label)
        for x_col, (color, marker, style, hollow) in zip((7.55, 16.15), row_samples):
            key.sample(x + x_col, baseline - 1.6, color, marker, style, hollow=hollow)


def _vertical_cell_key(ax, labels, colors, first_baseline, *, hollow=False):
    key = _PanelKey(ax)
    for i, (label, color) in enumerate(zip(labels, colors)):
        baseline = first_baseline + i * 5.6
        key.label(12.27, baseline, label)
        key.sample(7.55, baseline - 1.6, color, 's' if label=='Ab' else 'o',
                   '--' if label=='Ab' else '-', hollow=hollow)


def _citrate_key(ax, doses, dx, dy):
    key = _PanelKey(ax)
    key.label(4.16 + dx, 8.71 + dy, r'$[\mathrm{Ga}^{3+}]$ (mM)')
    key.label(6.95 + dx, 48.74 + dy, '[citric acid] (mM)', rotation=90)
    for col, ga in enumerate(('0', '1')):
        key.label(13.86 + dx + 8.6 * col, 13.68 + dy, ga, ha='center')
    for i, dose in enumerate(doses):
        baseline = 16.78 + dy + 5.6 * i
        key.label(26.16 + dx, baseline, f'{dose:g}')
        for j, (marker, style) in enumerate((('o', '-'), ('s', '--'))):
            key.sample(13.86 + dx + 8.6 * j, baseline - 1.7, _citrate_color(dose), marker, style)


def _spread(ax, data, x, y, step, color, series, *, marker="o", ls="-", open_marker=False,
            alpha=.10, size=9):
    """Plot independent-unit points at recorded times and acquisition means ± sample SD."""
    grouped = data.groupby(step, sort=True).agg(t=(x, "mean"), avg=(y, "mean"), sd=(y, "std"))
    # Keep acquisition rows with an unavailable ordinate so lines break at missing values.
    ax.fill_between(grouped.t, grouped.avg - grouped.sd, grouped.avg + grouped.sd,
                    facecolor=color, linewidth=0, alpha=alpha, zorder=1,
                    gid=series + " | mean ± sample SD")
    ax.plot(grouped.t, grouped.avg, color=color, linestyle=ls, linewidth=.8, zorder=2,
            gid=series + " | mean")
    points = ax.scatter(data[x], data[y], marker=marker, s=size, edgecolors=color,
                        facecolors="none" if open_marker else color, linewidths=.32,
                        alpha=1, zorder=3, gid=series + " | observations")
    identifier = next(col for col in ("plate_id", "flask_ID", "tube_ID") if col in data)
    points.plotted_value_labels = list(data[identifier].astype(str))





def _citrate_color(dose):
    levels = [0, 1, 3, 6, 12, 18, 24, 30, 36, 50]
    return plt.get_cmap("inferno")(np.interp(dose, levels, np.linspace(.1, .85, len(levels))))


def _citrate_trajectories(t, axes):
    flask = t["validation_flask_points"]
    batches = flask.loc[flask.nominal_citrate_mM.gt(0), "batch"].unique()
    other_additions = ["nominal_In_mM", "nominal_Al_mM", "nominal_V_mM",
                       "nominal_Cu_mM", "nominal_NaCl_mM"]
    flask = flask[flask.batch.isin(batches) & flask.kind.eq("biotic")
                  & flask[other_additions].eq(0).all(axis=1)].copy()
    plate = t["citrate_wavelength_pairs"]
    flask = flask.rename(columns={"nominal_Ga_mM": "Ga_mM", "nominal_citrate_mM": "citrate_mM"})
    keys = ["citrate_mM", "Ga_mM"]
    matched = flask[keys].drop_duplicates().merge(plate[keys].drop_duplicates(), on=keys)
    flask = flask.merge(matched, on=keys, validate="many_to_one")
    plate = plate.merge(matched, on=keys, validate="many_to_one")
    for (dose, ga, batch), g in flask.groupby(keys + ["batch"]):
        _spread(axes[0], g, "days_mean", "biotic_change", "sampling_event", _citrate_color(dose),
                f"AF flask apparent Fe(III) | {dose:g} mM citric acid | {ga:g} mM Ga(III) | {batch}",
                marker="s" if ga else "o", ls="--" if ga else "-", alpha=.09)
    for (dose, ga), g in plate.groupby(keys):
        for ax, measure in zip(axes[1:3], ("net_change_380", "abiotic_change_380")):
            _spread(ax, g, "days_mean_380", measure, "ordinal", _citrate_color(dose),
                    f"{measure} | {dose:g} mM citric acid | {ga:g} mM Ga(III)", marker="s" if ga else "o", ls="--" if ga else "-", alpha=.09)
    doses = sorted(matched.citrate_mM.unique())
    for i, (ax, ylabel) in enumerate(zip(axes, (r"$[\mathrm{Fe}^{3+}]_{\mathrm{app.}}$ (mM)",
                            r"$\Delta A_{380,\mathrm{net}}$ (AU)", r"$\Delta A_{380,\mathrm{Ab}}$ (AU)"))):
        ax.set(xlabel="time (d)", ylabel=ylabel, xlim=(0,7), xticks=range(8))
        if i == 0:
            # Flasks include a later acquisition and signed changes from baseline.
            stop = max(7, int(np.ceil(ax.dataLim.xmax)))
            low = min(0, 5 * np.floor(ax.dataLim.ymin / 5))
            high = max(80, 20 * np.ceil(ax.dataLim.ymax / 20))
            ax.set(xlim=(0,stop), xticks=range(stop+1), ylim=(low,high),
                   yticks=np.arange(20 * np.ceil(low/20), high+1, 20))
        else:
            ax.set(ylim=(-.2344,3.507), yticks=range(4))
            ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
        _citrate_key(ax, doses, (.24,0,0)[i], (.2,.75,0)[i])


def _paired_wavelengths(ax, pairs):
    keys = ["species", "experiment_id", "citrate_mM", "Ga_mM", "batch_id", "plate_id"]
    final = pairs.ordinal.eq(pairs.groupby(keys, dropna=False).ordinal.transform("max"))
    end = pairs[final & pairs.pair_available].drop_duplicates(keys)
    key = _PanelKey(ax)
    key.label(55.04, 8.91, r'$[\mathrm{Ga}^{3+}]$ (mM)')
    key.label(77.94, 16.4, r'$\lambda$ (nm)')
    colors = {380: ('#0072B2','#D55E00'), 428: ('#009E73','#CC79A7')}
    for col in (0,1):
        key.label(64.74 + col*8.6, 17.04, str(col), ha='center')
    for wave, marker in ((380, "o"), (428, "s")):
        field = f"net_change_{wave}"
        for ga, g in end.groupby("Ga_mM"):
            batch = g.groupby(["citrate_mM", "batch_id"])[field].mean()
            summary = batch.groupby("citrate_mM").agg(["mean", "std"])
            color = colors[wave][int(ga)]
            style = "--" if ga else "-"
            series = f"net endpoint change | {wave} nm | {ga:g} mM Ga(III)"
            ax.fill_between(summary.index, summary["mean"] - summary["std"],
                            summary["mean"] + summary["std"], color=color, alpha=.12, lw=0,
                            gid=series + " | equal-batch mean ± between-batch SD")
            ax.plot(summary.index, summary["mean"], color=color, ls=style, lw=.35,
                    gid=series + " | equal-batch mean")
            points = ax.scatter(g.citrate_mM, g[field], edgecolors=color,
                                facecolors="white" if ga else color, marker=marker,
                                s=9, linewidths=.45, alpha=1, gid=series + " | paired microplates")
            points.plotted_value_labels = list(g.plate_id.astype(str) + " | " + g.batch_id.astype(str))
            batch_points = ax.scatter(batch.index.get_level_values("citrate_mM"), batch.values,
                                      color=color, marker="_", s=18, linewidths=.65, zorder=4,
                                      gid=series + " | batch means")
            batch_points.plotted_value_labels = list(batch.index.get_level_values("batch_id"))
        baseline = 22.64 if wave == 380 else 28.24
        key.label(77.94, baseline, str(wave))
        for ga, color in enumerate(colors[wave]):
            key.sample(64.74+ga*8.6, baseline-1.6, color, marker, '--' if ga else '-', hollow=bool(ga))
    ax.set(xlabel="[citric acid] (mM)", ylabel=r"$\Delta A_{\lambda,\mathrm{net}}$ (AU)",
           xlim=(-1.8,37.8), xticks=[0,10,20,30], ylim=(.3693,3.43), yticks=[.8,1.6,2.4,3.2])
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))


def _format_comparisons(ax, folds, *, timing):
    field = "time_log2_fold" if timing else "response_log2_fold"
    data = folds[folds.readout.eq("raw_A380") & folds.plate_dose_mM.gt(0)
                 & np.isfinite(folds[field])]
    for (additive, ga), (label, color, marker) in ADDITIVES.items():
        g = data[data.additive.eq(additive)]
        if additive == "citric_acid":
            g = g[g.Ga_mM.eq(ga)]
        if g.empty:
            continue
        for exact, subset in g.groupby(g.match_status.eq("equal_reported_dose")):
            series = f"{field} | {additive}"
            if additive == "citric_acid":
                series += f" | {ga} mM Ga(III)"
            series += " | equal concentrations" if exact else " | neighboring concentrations"
            points = ax.scatter(subset.plate_dose_mM / 1000., subset[field], s=9, marker=marker,
                                facecolors=color if exact else "none", edgecolors=color, linewidths=.5,
                                zorder=3, gid=series)
            points.plotted_value_labels = list(subset.comparison_id.astype(str))
    if timing:
        last = data[data.final_read_maximum_present]
        points = ax.scatter(last.plate_dose_mM / 1000., last[field], marker="x", color="#303030",
                            linewidths=.45, s=5, zorder=4, gid="timing ratio | final-read maximum")
        points.plotted_value_labels = list(last.additive + " | " + last.comparison_id.astype(str))
    ax.axhline(0, lw=.5, ls=":", color=".4", gid="equal microplate and flask response")
    ax.set(xscale="log", xlabel="[plate additive] (M)",
           ylabel=r"$\log_2(t_{\max,\mathrm{p}}/t_{\max,\mathrm{f}})$" if timing else
                  r"$\log_2(\mathrm{AUC}_{\mathrm{rel,p}}/\mathrm{AUC}_{\mathrm{rel,f}})$")
    ax.set_xlim(2.255e-7,.293)
    majors = 10.0 ** np.arange(-6,0)
    ax.xaxis.set_major_locator(FixedLocator(majors))
    ax.xaxis.set_major_formatter(FuncFormatter(lambda v, pos: rf"$10^{{{int(np.rint(np.log10(v)))}}}$"
                                               if int(np.rint(np.log10(v))) % 2 == 0 else ""))
    ax.xaxis.set_minor_locator(FixedLocator(10.0 ** (np.arange(-66,-5) / 10)))
    ax.xaxis.set_minor_formatter(NullFormatter())
    if timing:
        ax.set(ylim=(-1.501,2.618), yticks=[-1,0,1,2])
    else:
        ax.set(ylim=(-2.7742,6.6116), yticks=[-2.5,0,2.5,5])
        ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    key = _PanelKey(ax)
    for ((additive, _), (label,color,marker)), baseline in zip(ADDITIVES.items(), (8.78,15.58,22.38,28.32,34.56)):
        key.label(7.93, baseline, label)
        key.sample(4.6, baseline-1.6, color, marker, 'none')
    extras = [('citric acid, 0 mM Ga', ADDITIVES['citric_acid',0][1], 'P', False),
              ('citric acid, 1 mM Ga', ADDITIVES['citric_acid',1][1], 'P', False),
              ('equal dose', '#444444', 'o', False), ('nearby doses','#444444','o',True)]
    if timing:
        extras.append(('final-read maximum', '#444444', 'x', False))
    for i,(label,color,marker,hollow) in enumerate(extras):
        baseline = 9.35 + 5.375*i
        key.label(29.93, baseline, label)
        key.sample(26.4, baseline-1.6, color, marker, 'none', hollow=hollow)


def _lf_microplates(ax, t):
    raw = t["validation_lf_plate_trajectories"]
    raw = raw[raw.signal.eq("raw_A380")].copy()
    composition = raw.component_final_mM_json.map(json.loads)
    benchmark = t["lf_plate_benchmark"]
    for row in benchmark.itertuples():
        label = row.LF_plate_context
        if label == "HH":
            keep = composition.map(lambda q: all(v == 0 for v in q.values()))
        else:
            # The formulation metadata record salt molarity. Match the cation
            # molarity in the comparison, including both cations per sulfate salt.
            compound = "gallium_sulfate" if label == "Ga" else "indium_sulfate"
            ion_mM = row.LF_Ga_mM if label == "Ga" else row.LF_In_mM
            keep = composition.map(lambda q: np.isclose(2 * q.get(compound, 0), ion_mM))
        subset = raw[keep]
        if subset.plate_id.nunique() != row.LF_plate_n:
            raise ValueError(f"LF {label} trajectory count differs from its endpoint comparison")
        for component, marker, ls in (("biotic", "o", "-"), ("abiotic", "s", "--")):
            series = f"{'LF' if component == 'biotic' else 'Ab'} | {label} | change in A380"
            _spread(ax, subset, "days_mean", f"{component}_change", "ordinal", METAL[label], series,
                    marker=marker, ls=ls)
    row = benchmark.iloc[0]
    batches = t["lf_plate_benchmark__batch_estimates"]
    points = ax.scatter(batches.target_day_d, batches["mean"], s=9, marker="D",
                        edgecolors="#882255", facecolors="none", linewidths=.45, zorder=3,
                        gid="AF | interpolated A380 change | batch estimates")
    points.plotted_value_labels = list(batches.batch_id)
    line, caps, bars = ax.errorbar(row.AF_reference_day_d, row.AF_reference_mean_AU,
                                   yerr=row.AF_reference_batch_sd_AU, fmt="D", ms=2.7,
                                   color="#882255", elinewidth=.65, capsize=2)
    line.set_gid("AF | interpolated A380 change | equal-batch mean")
    for artist in (*caps, *bars):
        artist.set_gid("AF | interpolated A380 change | between-batch SD")
    ax.annotate("AF at LF endpoint", xy=(row.AF_reference_day_d, row.AF_reference_mean_AU),
                xytext=(2.92 / 97.2, 1 - 20.06 / 97.2), textcoords="axes fraction",
                ha="left", va="baseline", fontsize=6, color="#882255",
                arrowprops=dict(arrowstyle="->", color="#882255", linewidth=.45))
    _metal_key(ax, plate=True)
    ax.set(xlabel="time (d)", ylabel=r"$\Delta A_{380}$ (AU)",
           xlim=(0,5.5), xticks=range(6), ylim=(0,3.865), yticks=range(4))
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))


def _metal_key(ax, *, plate=False, heading=True):
    prefix = '+ ~' if plate else '+ '
    rows = ['HH', prefix+r'$10^{-2}$ M In$^{3+}$', prefix+r'$10^{-1}$ M Ga$^{3+}$']
    symbols = [[(color, 's', '--', False),(color, 'o', '-', False)] for color in METAL.values()]
    _cell_matrix(ax, ['Ab','LF'], rows, 'medium' if heading else '', symbols,
                 x=38.48 if plate else 0, y=55.45 if plate else 7.68)


def _tube_time(ax, frame, *, orp, inset=False):
    measure = "orp_mV" if orp else "pH"
    shared = frame[frame.is_shared_initial_reference].drop_duplicates("initial_reference_ID")
    measured = frame[frame.is_independent_measurement]
    for (org, treatment), group in measured.groupby(["organism_code", "treatment_ID"]):
        color = METAL[TREATMENT[treatment]] if orp else CELLS[org]
        marker = "s" if org == "Ab" else "o"
        ls = "--" if org == "Ab" else "-"
        series = f"{org} | {treatment} | tube {measure}"
        _spread(ax, group, "time_average_d", measure, "sampling_event", color, series, marker=marker, ls=ls)
        initial = shared[shared.treatment_ID.eq(treatment)]
        first = group[group.sampling_event.eq(group.sampling_event.min())]
        x = [initial.time_average_d.mean(), first.time_average_d.mean()]
        avg = np.array([initial[measure].mean(), first[measure].mean()])
        sd = np.array([initial[measure].std(), first[measure].std()])
        ax.plot(x, avg, color=color, ls=ls, lw=.8, gid=series + " | shared start to first sample mean")
        ax.fill_between(x, avg - sd, avg + sd, color=color, lw=0, alpha=.10,
                        gid=series + " | shared start to first sample mean ± sample SD")
    for treatment, base in shared.groupby("treatment_ID"):
        color = METAL[TREATMENT[treatment]] if orp else ".35"
        points = ax.scatter(base.time_average_d, base[measure], marker="D", s=9,
                            color=color, linewidths=.3, zorder=4,
                            gid=f"{treatment} | shared starting-medium {measure}")
        points.plotted_value_labels = list(base.initial_reference_ID.astype(str))
    if inset:
        ax.set(xlim=(0,7), ylim=(390,515), xticks=[0,3,6], yticks=[400,450,500])
        ax.tick_params(which='major', length=2, labelsize=5, pad=.6, direction='in', top=True, right=True)
        ax.tick_params(which='minor', length=1, direction='in', top=True, right=True)
        ax.xaxis.set_minor_locator(AutoMinorLocator(2))
        ax.yaxis.set_minor_locator(AutoMinorLocator(2))
        return
    ax.set_xlabel("time (d)")
    if orp:
        ax.set(ylabel="ORP (mV)", xlim=(0,7), xticks=range(8),
               ylim=(383.12,728.6), yticks=range(400,701,50))
        _metal_key(ax)
        child = ax.inset_axes([.1355,.25,.45,.45])
        _tube_time(child, frame, orp=True, inset=True)
    else:
        ax.set(ylabel="pH", xlim=(-1.05,22.05), xticks=range(0,21,5),
               ylim=(1.297,2.551), yticks=np.arange(1.4,2.41,.2))
        _vertical_cell_key(ax, list(CELLS), list(CELLS.values()), 81.40)


def _violin(ax, group, position, color, marker, ls, *, lf):
    endpoints = group.dropna(subset=["apparent_OD600"]).sort_values("apparent_OD600")
    values = endpoints.apparent_OD600.to_numpy()
    if not len(values):
        return
    series = f"{endpoints.organism_code.iloc[0]} | {endpoints.treatment_ID.iloc[0]} | endpoint apparent OD600"
    if len(values) >= 2 and np.ptp(values) > 0:
        violin = ax.violinplot(values, positions=[position], widths=.27 if lf else .65,
                               showmeans=False, showmedians=False, showextrema=False)
        for body in violin["bodies"]:
            body.set(facecolor=color, edgecolor=color, alpha=.15, linewidth=.45,
                     gid=series + " | violin density")
    median_half_width = .10 if lf else .24
    ax.hlines(np.median(values), position - median_half_width, position + median_half_width,
              color=color, lw=.65, ls=ls,
              gid=series + " | median")
    half_range = .055 if lf else .13
    offsets = np.linspace(-half_range, half_range, len(values)) if len(values) > 1 else [0.]
    points = ax.scatter(position + np.array(offsets), values, color=color, marker=marker,
                        s=9, linewidths=.3, zorder=4, gid=series + " | individual tubes")
    points.plotted_value_labels = list(endpoints.tube_ID.astype(str))


def _tube_endpoints(ax, frame, *, lf):
    end = frame[frame.record_role.eq("endpoint")]
    if end.tube_ID.duplicated().any():
        raise ValueError("A tube has multiple terminal observations")
    if lf:
        for i, (treatment, label) in enumerate(TREATMENT.items(), 1):
            for org, offset, marker, ls in (("Ab", -.17, "s", "--"), ("LF", .17, "o", "-")):
                part = end[end.treatment_ID.eq(treatment) & end.organism_code.eq(org)]
                _violin(ax, part, i + offset, METAL[label], marker, ls, lf=True)
        ax.set_xticks([1,2,3], ['HH\n(baseline)',r'$+\ 10^{-2}$ M'+'\n'+r'In$^{3+}$',
                                r'$+\ 10^{-1}$ M'+'\n'+r'Ga$^{3+}$'])
        ax.set(xlim=(.564,3.436), ylim=(0,.006), yticks=[0,.002,.004,.006])
        ax.yaxis.set_major_formatter(FormatStrFormatter('%.3f'))
        _metal_key(ax, heading=False)
    else:
        for i, (org, color) in enumerate(CELLS.items(), 1):
            _violin(ax, end[end.organism_code.eq(org)], i, color, {"Ab":"s","AF":"o","AT":"^"}[org], "-", lf=False)
        ax.set_xticks([1, 2, 3], list(CELLS))
        ax.set(xlim=(.5,3.5), ylim=(0,end.apparent_OD600.max()*1.08), yticks=[0,.02,.04,.06,.08])
        ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    ax.set_ylabel(r"planktonic $\mathrm{OD}_{600,\mathrm{app.}}$")
    ax.xaxis.set_minor_locator(NullLocator())
    ax.yaxis.set_minor_locator(NullLocator())


def figure_7(tables):
    """Render the citrate, format and L. ferriphilum panels."""
    fig, axes = _grid(3, 3, [
        r"$\mathrm{Fe}^{3+}_{\mathrm{app.},\,\mathrm{AF}}$ in flasks",
        r"$\mathrm{Fe}^{3+}_{\mathrm{app.},\,\mathrm{AF}}$ in plates",
        r"$\mathrm{Fe}^{3+}_{\mathrm{app.},\,\mathrm{Ab}}$ in plates",
        r"$\mathrm{Fe}^{3+}_{\mathrm{app.},\,\Delta\mathrm{AF}}$ in plates = f($A_{\lambda}$)",
        r"plate-to-flask $\mathrm{AUC}_{\mathrm{Fe}^{3+}_{\mathrm{app.},\,\mathrm{AF}}}$",
        r"plate-to-flask $t_{\mathrm{max},\,\mathrm{Fe}^{3+}_{\mathrm{app.},\,\mathrm{AF}}}$",
        r"$\mathrm{Fe}^{3+}_{\mathrm{app.}}$ in plates =" "\n" "f(cells, [additive], t)",
        r"$\mathrm{Fe}^{3+}_{\mathrm{app.}}$ in tubes =" "\n" "f(cells, [additive], t)",
        r"$\mathrm{OD}_{600,\,\mathrm{app.},\,\mathrm{endpoint}}$" "\n" "in tubes = f(cells, [additive])",
    ])
    _citrate_trajectories(tables, axes[:3])
    _paired_wavelengths(axes[3], tables["citrate_wavelength_pairs"])
    for ax, timing in zip(axes[4:6], (False, True)):
        _format_comparisons(ax, tables["format_folds"], timing=timing)
    _lf_microplates(axes[6], tables)
    tubes = tables["tube_od600"]
    lf = tubes[tubes.medium.eq("HH") & tubes.incubation_temperature_C.eq(37)
               & tubes.organism_code.isin(["Ab", "LF"]) & tubes.treatment_ID.isin(TREATMENT)]
    _tube_time(axes[7], lf, orp=True)
    _tube_endpoints(axes[8], lf, lf=True)
    return fig


def _sulfur_panel(ax, t, group, corrected, *, inset=False):
    registry = t["sulfur_baseline_registry"]
    refs = registry[registry.source_panel.eq(group)].copy()
    if group == "b":
        refs = refs[refs.reference_id.eq("af_tt_5")]
    refs = refs.sort_values("initial_OD600_200_uL", kind="stable")
    readout = "corrected_R545_485" if corrected else "raw_R545_485"
    curves = t["sulfur_baseline_trajectories"]
    organism = "AT" if group == "c" else "AF"
    labels, legend_samples = [], []
    for i, (label, group_refs) in enumerate(refs.groupby("legend_label", sort=False)):
        selected = curves[curves.reference_id.isin(group_refs.reference_id) & curves.signal.eq(readout)]
        selected = selected.drop_duplicates(["batch_id", "plate_id", "ordinal"])
        marker = ('o','s','^')[i]
        for batch_id, batch in selected.groupby("batch_id"):
            for component, org, ls in (("biotic", organism, "-"), ("abiotic", "Ab", "--")):
                series = f"{org} | {label} | {batch_id} | {readout}"
                cell_marker = ('s' if org == 'Ab' else 'o') if group == 'b' else marker
                _spread(ax, batch, "days_mean", f"{component}_change", "ordinal", CELLS[org], series,
                        marker=cell_marker, ls=ls, open_marker=corrected, size=5 if inset else 9)
        if group == "c":
            labels.append(r"$\mathrm{S_2O_3^{2-}}$" if label.endswith("TS") else r"$\mathrm{S_4O_6^{2-}}$")
        elif group == "a":
            labels.append(rf"$10^{{{int(round(np.log10(group_refs.initial_OD600_200_uL.iloc[0])))}}}$")
        legend_samples.append([(CELLS['Ab'],marker,'--',corrected),(CELLS[organism],marker,'-',corrected)])
    if inset:
        ax.set(xlim=(0,10), ylim=(-16.1,3.3), xticks=[0,5,10], yticks=[-15,0])
        ax.tick_params(labelsize=5, pad=.6, length=1.5, direction="in", top=True, right=True)
        ax.xaxis.set_minor_locator(NullLocator())
        ax.yaxis.set_minor_locator(NullLocator())
        return
    x_end, x_step = {'a':(10,2),'b':(7,1),'c':(9,3)}[group]
    ax.set(xlabel="time (d)", ylabel=r"$R^*_t-R^*_0$" if corrected else r"$R_t-R_0$",
           xlim=(0,x_end),xticks=range(0,x_end+1,x_step))
    ordinate = {
        ('a',False): ((-.13,.30),[-.1,0,.1,.2,.3], '%.1f'),
        ('b',False): ((-.0177,.0955),[0,.02,.04,.06,.08], '%.2f'),
        ('c',False): ((-.0397,.2203),[0,.1,.2], '%.1f'),
        ('a',True): ((-2,1.6),[-2,-1,0,1], '%d'),
        ('b',True): ((-.1248,.2519),[-.1,0,.1,.2], '%.1f'),
        ('c',True): ((-.1418,.85),[0,.2,.4,.6], '%.1f'),
    }
    bounds, ticks, precision = ordinate[group,corrected]
    ax.set(ylim=bounds, yticks=ticks)
    ax.yaxis.set_major_formatter(FormatStrFormatter(precision))
    if group == 'c' and corrected:
        ax.yaxis.set_minor_locator(FixedLocator([-.1,.1,.3,.5,.7,.8]))
    if group == "b":
        _vertical_cell_key(ax, ['Ab','AF'], [CELLS['Ab'],CELLS['AF']],6.92,hollow=corrected)
    else:
        dx = -2.1 if group=='c' and corrected else 0
        y = 8.48 if group=='a' else 7.09 if corrected else 7.69
        _cell_matrix(ax, ['Ab',organism], labels, r'$\mathrm{OD}_{600}$' if group=='a' else '20 mM',
                     legend_samples, x=dx, y=y)
    if group == "a" and corrected:
        child = ax.inset_axes([.52,.09,.44,.30])
        _sulfur_panel(child, t, group, True, inset=True)


def figure_8(tables):
    """Render measured and scatter-corrected sulfur indicators and tube validation."""
    fig, axes = _grid(2, 4, [
        r"$\mathrm{pH}_{\mathrm{app.},\,\mathrm{S}_{2}\mathrm{O}_{3}^{2-}}$ in plates ="
        "\n" r"f(cells, $\mathrm{OD}_{600}$, t)",
        r"$\mathrm{pH}_{\mathrm{app.},\,\mathrm{S}_{4}\mathrm{O}_{6}^{2-}}$ in plates = f(cells, t)",
        r"$\mathrm{pH}_{\mathrm{app.}}$ in plates =" "\n" "f(cells, energy source, t)",
        r"$\mathrm{pH}_{\mathrm{app.},\,\mathrm{corr.},\,\mathrm{S}_{2}\mathrm{O}_{3}^{2-}}$ in plates ="
        "\n" r"f(cells, $\mathrm{OD}_{600}$, t)",
        r"$\mathrm{pH}_{\mathrm{app.},\,\mathrm{corr.},\,\mathrm{S}_{4}\mathrm{O}_{6}^{2-}}$"
        "\n" "in plates = f(cells, t)",
        r"$\mathrm{pH}_{\mathrm{app.},\,\mathrm{corr.}}$ in plates =" "\n" "f(cells, energy source, t)",
        r"pH, $\mathrm{S}_{4}\mathrm{O}_{6}^{2-}$" "\n" "in tubes = f(cells, t)",
        r"$\mathrm{OD}_{600,\,\mathrm{app.},\,\mathrm{endpoint}}$" "\n" "in tubes = f(cells)",
    ])
    for index, ax in enumerate(axes[:6]):
        _sulfur_panel(ax, tables, "abc"[index % 3], index >= 3)
    tubes = tables["tube_od600"]
    sulfur = tubes[tubes.medium.eq("SM4_tetrathionate") & tubes.incubation_temperature_C.eq(30)
                  & tubes.treatment_ID.eq("no_added_metal") & tubes.organism_code.isin(CELLS)]
    _tube_time(axes[6], sulfur, orp=False)
    _tube_endpoints(axes[7], sulfur, lf=False)
    return fig
