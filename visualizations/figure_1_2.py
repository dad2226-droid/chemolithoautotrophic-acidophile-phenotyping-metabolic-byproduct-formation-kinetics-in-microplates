"""Chemical calibration and pH stability panels for Figures 1 and 2."""
from __future__ import annotations

import numpy as np
import pandas as pd
from matplotlib import colormaps, pyplot as plt
from matplotlib.cm import ScalarMappable
from matplotlib.colors import LogNorm, Normalize
from matplotlib.lines import Line2D
from matplotlib.ticker import AutoMinorLocator, FormatStrFormatter, FuncFormatter, LogLocator, MultipleLocator, NullFormatter, NullLocator
from matplotlib.transforms import Affine2D, ScaledTranslation


PH = ("#0072B2", "#009E73", "#332288", "#56B4E9")
IRON = ("#D55E00", "#CC79A7", "#E69F00", "#882255")
DYES = {0: ("#555555", "^"), 1: (PH[0], "o"), 10: (PH[1], "s")}
R = r"$R$"
RS = r"$R^{*}$"
TOTAL = r"measured $[\mathrm{Fe}]_{\mathrm{total}}$ (M)"
SIGNALS = ("ratio_raw", "ratio_corrected")
PLATES = ("UV_transparent", "non_UV_transparent")


def _canvas(titles, *, width, height, row_pitch, title_pads):
    """Position square axes from the printed panel dimensions, in points."""
    size, column_pitch = 97.2, 129.6
    fig = plt.figure(figsize=(width/72, height/72))
    axes = []
    for index, title in enumerate(titles):
        row, column = divmod(index, 4)
        left = 23.855 + column * column_pitch
        top = 10.575 + row * row_pitch
        ax = fig.add_axes((left/width, (height-top-size)/height, size/width, size/height))
        ax.set_title(title, fontsize=7, pad=title_pads[index])
        ax.text(-.248, 1.029, chr(97 + index), transform=ax.transAxes, fontsize=8,
                fontweight="bold", ha="left", va="baseline", clip_on=False)
        ax.tick_params(which="major", direction="in", top=True, right=True, length=3,
                       labelsize=6, pad=2)
        ax.tick_params(which="minor", direction="in", top=True, right=True, length=1.6)
        for axis in (ax.xaxis, ax.yaxis):
            axis.set_minor_locator(AutoMinorLocator(2))
            axis.labelpad = .5
        axes.append(ax)
    return fig, axes


def _points(ax, x, y, color, marker, *, series, ids, opened=False, edge=.3, alpha=1):
    points = ax.scatter(x, y, s=9, marker=marker, edgecolors=color,
               facecolors="none" if opened else color, linewidths=edge, alpha=alpha, zorder=3, gid=series)
    points.plotted_value_labels = list(ids)


def _handle(color, marker, label="", linestyle="-", opened=False, edge=.3):
    return Line2D([], [], color=color, marker=marker, ls=linestyle, lw=.8, ms=3,
                  mew=edge, mec=color, mfc="none" if opened else color, label=label)


def _key_transform(ax):
    """Legend coordinates in physical points, measured from the upper-left corner."""
    return (Affine2D().scale(1/72, -1/72) + ax.figure.dpi_scale_trans
            + ScaledTranslation(0, 1, ax.transAxes))


def _key_text(ax, x, y, text, *, align="left", rotation=0):
    ax.text(x, y, text, transform=_key_transform(ax), fontsize=5, ha=align,
            va="center", rotation=rotation, zorder=6)


def _key_sample(ax, x, y, handle, width):
    sample = Line2D([x-width/2, x, x+width/2], [y, y, y],
                    color=handle.get_color(), lw=handle.get_linewidth(),
                    linestyle=handle.get_linestyle(), marker=handle.get_marker(),
                    markersize=3, markeredgewidth=handle.get_markeredgewidth(),
                    markerfacecolor=handle.get_markerfacecolor(), markevery=[1],
                    transform=_key_transform(ax), zorder=6)
    if handle.get_linestyle() == "--":
        sample.set_dashes((1.4, 1.1))
    ax.add_artist(sample)


def _list_key(ax, handles, *, x, y, text_x, title=None, width=7):
    if title:
        tx, ty, heading = title
        _key_text(ax, tx, ty, heading)
    for index, handle in enumerate(handles):
        row = y + 5.6 * index
        _key_sample(ax, x, row, handle, width)
        _key_text(ax, text_x, row, handle.get_label())


def _matrix_key(ax, row_labels, column_labels, rows, *, xs, row_y, label_x,
                column_y, heading, side_heading=None):
    """Place every condition directly under its column and beside its row label."""
    tx, ty, name, align = heading
    _key_text(ax, tx, ty, name, align=align)
    if side_heading:
        tx, ty, name = side_heading
        _key_text(ax, tx, ty, name, align="center", rotation=90)
    for x, label in zip(xs, column_labels):
        _key_text(ax, x, column_y, label, align="center")
    for i, (label, handles) in enumerate(zip(row_labels, rows)):
        y = row_y + 5.6 * i
        _key_text(ax, label_x, y, label)
        for x, handle in zip(xs, handles):
            if handle is not None:
                _key_sample(ax, x, y, handle, width=6.4)


def _band(ax, frame, group, x, y, color, marker, ls="-", opened=False, *, series, edge=.3, point_x=None):
    """Draw each physical plate and summarize plates at each standard or acquisition."""
    summary = frame.groupby(group, dropna=False).agg(
        coord=(x, "mean"), avg=(y, "mean"), sd=(y, "std")).sort_values("coord")
    line, = ax.plot(summary.coord, summary.avg, color=color, ls=ls, lw=.8, gid=series + " | mean")
    if ls == "--":
        line.set_dashes((1.4, 1.1))
    ax.fill_between(summary.coord, summary.avg-summary.sd, summary.avg+summary.sd,
                    color=color, alpha=.18, lw=0, gid=series + " | mean ± SD")
    _points(ax, frame[point_x or x], frame[y], color, marker, series=series + " | physical plate",
            ids=frame.physical_plate_id, opened=opened, edge=edge)


def _inset(ax, bounds):
    child = ax.inset_axes(bounds)
    child.tick_params(axis="both", which="major", labelsize=5, direction="in",
                      top=True, right=True, length=1.5, pad=.5)
    child.tick_params(axis="both", which="minor", direction="in", top=True, right=True, length=.8)
    for axis in (child.xaxis, child.yaxis):
        axis.label.set_size(5)
        axis.labelpad = .5
        axis.set_minor_locator(AutoMinorLocator(2))
    return child


def _color_key(ax, norm, cmap, bounds, title, ticks):
    cax = ax.inset_axes(bounds)
    bar = ax.figure.colorbar(ScalarMappable(norm=norm, cmap=cmap), cax=cax, ticks=ticks)
    cax.tick_params(axis="y", which="major", labelsize=5, length=2, pad=1,
                    direction="in", left=True, right=True)
    cax.tick_params(axis="y", which="minor", length=1, direction="in", left=True, right=True)
    if isinstance(norm, LogNorm):
        cax.yaxis.set_minor_locator(LogLocator(subs=10.**(np.arange(1,10)/10),numticks=100))
        cax.yaxis.set_minor_formatter(NullFormatter())
    else:
        cax.yaxis.set_minor_locator(MultipleLocator(.5))
    if title:
        cax.set_title(title, fontsize=5, pad=1)
    return bar


def _fit(ax, frame, color, scale=1., ls="-", *, series):
    frame = frame.sort_values("x")
    ax.plot(frame.x * scale, frame.fitted, color=color, ls=ls, lw=.9, gid=series + " | fitted mean")
    ax.fill_between(frame.x * scale, frame.ci_low, frame.ci_high, color=color, alpha=.17, lw=0,
                    gid=series + " | fitted 95% confidence interval")


def _spectral_ph(ax, frame, replicates):
    norm = Normalize(frame.mean_measured_pH.min(), frame.mean_measured_pH.max())
    cmap = colormaps["viridis"]
    zoom = _inset(ax, [.25, .39, .50, .55])
    for (block, column), series in frame.groupby(["row_block", "plate_column"]):
        series = series.sort_values("wavelength_nm")
        ph = series.mean_measured_pH.iloc[0]
        group = f"standard {block}/{column} | measured pH {ph:.6g}"
        color = cmap(norm(ph))
        points = replicates[replicates.row_block.eq(block) & replicates.plate_column.eq(column)]
        for panel in (ax, zoom):
            visible = points if panel is ax else points[points.wavelength_nm.between(350, 650)]
            _points(panel, visible.wavelength_nm, visible.mean_absorbance, color, "o", edge=0, alpha=.27,
                    ids=visible.physical_plate_id, series=group + " | physical plate")
            panel.plot(series.wavelength_nm, series["mean"], color=color, lw=.65, gid=group + " | mean")
            panel.fill_between(series.wavelength_nm, series["mean"]-series.sd,
                               series["mean"]+series.sd, color=color, alpha=.10, lw=0, gid=group + " | mean ± SD")
    window = frame[frame.wavelength_nm.between(350, 650)]
    zoom.set(xlim=(350, 650), xticks=[400, 500, 600], yticks=[0, .1],
             ylim=(min(0,window["mean"].min()), 1.06*(window["mean"]+window.sd).max()))
    zoom.yaxis.set_major_formatter(FormatStrFormatter("%.1f"))
    zoom.xaxis.set_minor_locator(MultipleLocator(50))
    zoom.yaxis.set_minor_locator(MultipleLocator(.05))
    for wave in (485, 545):
        for panel in (ax, zoom):
            panel.axvline(wave, color=".4", ls=":", lw=.55, gid=f"reference wavelength | {wave} nm")
        zoom.text(wave, .03, str(wave), rotation=90, transform=zoom.get_xaxis_transform(),
                  ha="right", va="bottom", fontsize=5)
    ax.set(xlabel="wavelength (nm)", ylabel=r"$A$ (AU)", xlim=(200,800), ylim=(0,2),
           xticks=[200,400,600,800], yticks=np.arange(0,2.1,.5))
    ax.yaxis.set_major_formatter(FormatStrFormatter("%.1f"))
    _color_key(ax,norm,cmap,[.80,.39,.025,.48],"pH",range(1,6))


def _ph_calibration(axes, tables):
    raw, corrected, precision = axes
    levels = tables["calibration_level_means"]
    fits = tables["response_curves"]
    dye = levels[levels.experiment_id.isin(["pH_configuration_23_C_UV",
                                            "pH_UV_20uM_stock_initial_and_spectra"])
                 & levels.condition.eq("thymol_blue")].copy()
    blank = levels[levels.experiment_id.eq("final_indicator_acid")
                   & levels.signal.eq("ratio_raw") & levels.condition.eq("no_thymol_blue")].copy()
    coord = ["row_block", "plate_column"]
    blank["ph"] = blank.groupby(coord).reference_pH.transform("mean")
    _band(raw, blank, coord, "ph", "mean", *DYES[0], series="R | 0 µM dye", point_x="reference_pH")
    dye["ph"] = dye.groupby(["signal","experiment_id",*coord]).reference_pH.transform("mean")
    for ax, sig, name in zip((raw,corrected),SIGNALS,(R,RS)):
        open_symbol = sig==SIGNALS[1]
        for dose, sub in dye[dye.signal.eq(sig)].groupby("nominal_dye_uM"):
            color, marker = DYES[dose]
            series = f"{sig} | {dose:g} µM dye"
            _points(ax, sub.reference_pH, sub["mean"], color, marker, opened=open_symbol, edge=1,
                    series=series + " | physical plate", ids=sub.physical_plate_id)
            _fit(ax,fits[fits.signal.eq(sig) & fits.nominal_dye_uM.eq(dose)],color,series=series)
        handles = [_handle(c, m, str(d), "none" if open_symbol and d<10 else "-", open_symbol, edge=.5 if d<10 and open_symbol else (1 if d else .3))
                   for d, (c, m) in DYES.items()]
        _list_key(ax, handles, x=14.17, y=13.82, text_x=18.27, title=(7.33, 8.8, "dye (µM)"))
        ax.set(xlabel="measured pH",ylabel=name)
        ax.yaxis.set_major_formatter(FormatStrFormatter("%.1f"))
    for sig,ls in zip(SIGNALS,("-","--")):
        sub = dye[dye.nominal_dye_uM.eq(10) & dye.signal.eq(sig)]
        _band(precision,sub,coord,"ph","technical_sd",*DYES[10],ls,opened=sig==SIGNALS[1],
              series=f"{sig} | 10 µM dye | within-standard SD", edge=.5, point_x="reference_pH")
    _list_key(precision,[_handle(*DYES[10],R,edge=.5), _handle(*DYES[10],RS,"--",True,edge=.5)],
              x=10.83,y=11.42,text_x=14.93)
    precision.set(xlabel="measured pH",ylabel="technical-replicate SD")
    precision.yaxis.set_major_formatter(FormatStrFormatter("%.1f"))
    raw.set_ylim(.34413,3.30790)
    corrected.set_ylim(-.55649,8.36790)
    precision.set_ylim(-.03564,.64280)
    for ax in axes:
        ax.set_xlim(7,0)
        ax.xaxis.set_major_locator(MultipleLocator(1))
        ax.xaxis.set_major_formatter(FormatStrFormatter("%.1f"))


def _iron_spectra(ax,tables):
    from acidophile_data.indicator_and_tube_summaries import iron_spectra

    spectrum=tables["iron_spectral_replicates"]
    positive = spectrum.total_Fe_M[spectrum.total_Fe_M.gt(0)]
    norm = LogNorm(positive.min(),positive.max())
    cmap = colormaps["inferno"]
    for column, group in spectrum.groupby("plate_column"):
        c = group.total_Fe_M.mean()
        color = cmap(norm(c)) if c>0 else "#777777"
        label = f"standard column {column} | measured total Fe {c:.6g} M"
        curve=group.groupby("wavelength_nm").agg(avg=("difference","mean"),
                    ferric_sd=("difference","std"),ferrous_sd=("ferrous_reference_sd","mean"))
        sd=np.hypot(curve.ferric_sd,curve.ferrous_sd)
        _points(ax, group.wavelength_nm, group.difference, color, "o", edge=0, alpha=.27,
                series=label + " | ferric physical plate minus ferrous reference", ids=group.physical_plate_id)
        ax.plot(curve.index,curve.avg,color=color,lw=.8,gid=label + " | mean ferric minus ferrous")
        ax.fill_between(curve.index,curve.avg-sd,curve.avg+sd,color=color,alpha=.15,lw=0,
                        gid=label + " | mean ± propagated SD")
    for wave in (380,428):
        ax.axvline(wave,color=".4",ls=":",lw=.55,gid=f"reference wavelength | {wave} nm")
        ax.text(wave,.99,str(wave),rotation=90,ha="center",va="top",
                transform=ax.get_xaxis_transform(),fontsize=5)
    # Use the exact endpoint acquisitions; 428 nm is absent from the native scan grid.
    levels=tables["calibration_level_means"]
    rows=levels[levels.experiment_id.isin(["final_indicator_ferric","final_indicator_ferrous_1M"])
                & levels.condition.eq("no_thymol_blue") & levels.signal.isin(["a380","a428"])].copy()
    rows["wavelength_nm"]=rows.signal.map({"a380":380,"a428":428})
    endpoints=iron_spectra({"calibration_spectral_summary":rows.rename(columns={"mean":"mean_absorbance"})})
    zoom=_inset(ax,[.45086,.496,.34,.38])
    handles=[]
    for wave,marker in ((380,"o"),(428,"^")):
        data=endpoints[endpoints.wavelength_nm.eq(wave)&endpoints.total_Fe_M.gt(0)]
        points = zoom.scatter(data.total_Fe_M,data.difference,c=data.total_Fe_M,cmap=cmap,norm=norm,
                     s=9,marker=marker,linewidths=0,gid=f"ferric minus ferrous | {wave} nm | physical plate")
        points.plotted_value_labels = list(data.physical_plate_id)
        handles.append(_handle(".2",marker,rf"$A_{{{wave}}}$","none",edge=0))
    zoom.set(xscale="log", xlabel=r"$[\mathrm{Fe}]_{\mathrm{total}}$ (M)", ylabel=r"$A$ (AU)",
             xlim=(1e-5,.1), ylim=(-.14371,3), yticks=[0,1,2,3])
    zoom.xaxis.set_major_locator(LogLocator(numticks=12))
    zoom.xaxis.set_major_formatter(FuncFormatter(
        lambda x, _: rf"$10^{{{int(round(np.log10(x)))}}}$" if int(round(np.log10(x))) % 2 == 0 else ""))
    zoom.xaxis.set_minor_locator(LogLocator(subs=[2,4,6,8],numticks=100))
    zoom.xaxis.set_minor_formatter(NullFormatter())
    _list_key(zoom, handles, x=3.62,y=6.53,text_x=5.44)
    ax.set(xlabel="wavelength (nm)",ylabel=r"$A_{\mathrm{Fe}^{3+}}-A_{\mathrm{Fe}^{2+}}$ (AU)",
           xlim=(200,800),xticks=[200,400,600,800])
    ax.yaxis.set_major_formatter(FormatStrFormatter("%.1f"))
    decades=10.**np.arange(np.ceil(np.log10(positive.min())),np.floor(np.log10(positive.max()))+1)
    _color_key(ax,norm,cmap,[.83586,.496,.025,.37],None,decades)


def _iron_calibration(axes,tables):
    acid_ax, reference_ax, apparent_ax = axes
    levels=tables["calibration_level_means"]
    iron=levels[levels.experiment_id.str.startswith("iron_configuration_")
                & levels.condition.eq("no_thymol_blue") & levels.signal.eq("a380") & levels.axis_value.gt(0)]
    curves=tables["response_curves"]
    handles=[]
    for (acid,sub),color,marker in zip(iron.groupby("acid_M"),IRON,("o","s","^","v")):
        series=f"no dye | H₂SO₄ {acid:g} M"
        _points(acid_ax, sub.axis_value/1000, sub["mean"], color, marker,
                series=series + " | physical plate", ids=sub.physical_plate_id, edge=1)
        fit=curves[curves.scope.eq("iron_acid_specific") & curves.acid_M.eq(acid)&curves.x.gt(0)]
        _fit(acid_ax,fit,color,.001,series=series)
        handles.append(_handle(color,marker,rf"$10^{{{int(round(np.log10(acid)))}}}$", edge=1))
    acid_ax.set(xscale="log",xlabel=TOTAL,ylabel=r"$A_{380}$ (AU)",ylim=(0,3.41664),
                yticks=[0,1,2,3], xlim=(.000197,1.5326))
    acid_ax.yaxis.set_major_formatter(FormatStrFormatter("%.1f"))
    acid_ax.xaxis.set_major_locator(LogLocator(numticks=20))
    acid_ax.xaxis.set_minor_locator(LogLocator(subs=10.**(np.arange(1,10)/10),numticks=100))
    acid_ax.xaxis.set_minor_formatter(NullFormatter())
    _list_key(acid_ax, handles, x=13.6,y=14.75,text_x=17.7,
              title=(4.9,6.99,r"$[\mathrm{H_2SO_4}]$ (M)"))
    ref=tables["iron_chemical_reference_summary"]
    lines=tables["feii_regression_curves"]
    keys={}
    for name,sub in ref.groupby("series_ID"):
        acid=sub.nominal_H2SO4_concentration_M.iloc[0]
        ferric=name=="ferric"
        marker,color,ls=("s",IRON[1],"-") if ferric else ("o",IRON[0],"-" if acid==1 else "--")
        _points(reference_ax,sub.total_Fe_M,sub.FeII_M,color,marker,series=f"{name} | mean Fe(II) by cerimetry", ids=sub.reference_sample_ID, edge=1)
        _fit(reference_ax,lines[lines.series_ID.eq(name)],color,ls=ls,series=str(name))
        keys[(acid,ferric)]=_handle(color,marker," ",ls, edge=1)
    maximum=ref.total_Fe_M.max()*1.04
    reference_ax.plot([0,maximum],[0,maximum],color=".4",ls="--",lw=.65,zorder=0,gid="identity reference | Fe(II) equals total Fe")
    reference_ax.annotate(r"$[\mathrm{Fe}^{2+}]=[\mathrm{Fe}]_{\mathrm{total}}$",
                          xy=(.071,.073), xytext=(.0426,.0191), fontsize=5,
                          arrowprops=dict(arrowstyle="-", lw=.4, color=".4"))
    _matrix_key(reference_ax,[r"$10^{0}$",r"$10^{-3}$"],[r"$\mathrm{Fe}^{2+}$",r"$\mathrm{Fe}^{3+}$"],
                [[keys.get((acid,f)) for f in (False,True)] for acid in [1,.001]],
                xs=(11.89,24.89),row_y=18.12,label_x=33,column_y=6.3,
                heading=(33,5.5,r"$[\mathrm{H_2SO_4}]$ (M)","left"))
    reference_ax.tick_params(axis="y",which="both",pad=.25)
    reference_ax.set(xlabel=TOTAL,ylabel=r"measured $[\mathrm{Fe}^{2+}]$ (M)",ylim=(-.005,.100),
                     yticks=np.arange(0,.101,.02))
    chemical=tables["iron_optical_reference_levels"]
    chemical=chemical[chemical.series_ID.eq("ferric")].copy()
    chemical["fe_molar"]=chemical.FeIII_mM/1000
    handles=[]
    for sig,color,opened in (("a380",IRON[0],False),("a380_corrected",IRON[1],True)):
        row=[]
        for cond,marker,ls in (("no_thymol_blue","o","-"),("thymol_blue","s","--")):
            sub=chemical[chemical.signal.eq(sig)&chemical.condition.eq(cond)]
            _band(apparent_ax,sub,"plate_column","fe_molar","mean",color,marker,ls,opened,
                  series=f"{sig} | {cond.replace('_', ' ')}")
            row.append(_handle(color,marker," ",ls,opened))
        handles.append(row)
    _matrix_key(apparent_ax,[r"$A_{380}$",r"$A^{*}_{380}$"],["0","10"],handles,
                xs=(66.04,74.64),row_y=79.64,label_x=79.25,column_y=73.77,
                heading=(70.33,68.22,"[dye] (µM)","center"))
    apparent_ax.set(xlabel=r"$[\mathrm{Fe}^{3+}]_{\mathrm{app.}}=[\mathrm{Fe}]_{\mathrm{total}}-[\mathrm{Fe}^{2+}]$",
                    ylabel=r"$[\mathrm{Fe}^{3+}]_{\mathrm{app.}}$ indicator (AU)",
                    ylim=(-.15175,3.07170))
    for ax in (reference_ax,apparent_ax):
        ax.set_xticks([0,.1],["0",r"$10^{-1}$"])
        ax.xaxis.set_minor_locator(MultipleLocator(.01))


def figure_1(tables):
    fig,axes=_canvas([
        r"$\mathrm{pH}_{\mathrm{app.}}$ = f(pH, wavelength)",
        r"$\mathrm{pH}_{\mathrm{app.}}$ = f(pH)",
        r"$\mathrm{pH}_{\mathrm{app.,corr.}}$ = f(pH)",
        r"within-plate $\mathrm{pH}_{\mathrm{app.,(corr.)}}$ noise",
        r"$A_{\mathrm{Fe}^{3+}}-\mathrm{Fe}^{2+}$ = f([Fe], wavelength)",
        r"$\mathrm{Fe}^{3+}_{\mathrm{app.}}$ = f([$\mathrm{H}_{2}\mathrm{SO}_{4}$], [$\mathrm{Fe}^{3+}$])",
        "Fe speciation validation",
        r"$\mathrm{Fe}^{3+}_{\mathrm{app.,(corr.)}}$ = f([dye], [$\mathrm{Fe}^{3+}$])",
    ],width=513.72,height=262.5,row_pitch=137.028,
        title_pads=(4.57,4.57,4.57,4.57,5.38,5.38,6.36,5.38))
    _spectral_ph(axes[0], tables["spectral_means"], tables["spectral_replicates"])
    _ph_calibration(axes[1:4],tables)
    _iron_spectra(axes[4],tables)
    _iron_calibration(axes[5:],tables)
    return fig


def _lid_panel(ax, levels):
    lid=levels[levels.family.eq("lid_comparison") & levels.assay.eq("pH")]
    zoom=_inset(ax,[.10944,.28448,.67552,.67552])
    handles=[]
    for sig,color,opened in zip(SIGNALS,PH,(False,True)):
        row=[]
        for state,marker,ls in (("lid_absent","o","-"),("lid_present","s","--")):
            sub=lid[lid.signal.eq(sig) & lid.reading_label.eq(state)]
            series=f"{sig} | {state.replace('_',' ')}"
            _band(ax,sub,"level_id","predictor","response_mean",color,marker,ls,opened,series=series)
            if not opened:
                _band(zoom,sub,"level_id","predictor","response_mean",color,marker,ls,opened,series=series)
            row.append(_handle(color,marker," ",ls,opened))
        handles.append(row)
    ax.set(xlabel="measured pH",ylabel=r"ratiometric $\mathrm{pH}_{\mathrm{app.}}$ indicator",xlim=(5,.4),xticks=range(1,6),ylim=(-1.820,38.697),yticks=range(0,31,10))
    ax.xaxis.set_major_formatter(FormatStrFormatter("%.1f"))
    zoom.set(xlim=(5,.5),ylim=(.4,2.1),xticks=range(1,6),yticks=[1,2])
    zoom.tick_params(axis="x", pad=1)
    zoom.get_xticklabels()[0].set_ha("right")
    _matrix_key(zoom,[R,RS],["−","+"],handles,xs=(6.3,14.9),row_y=16.52,label_x=19.5,
                column_y=10.66,heading=(10.6,5.1,"lid","center"))


def _heatmap(ax, values, row_labels, column_labels, norm, series):
    cmap=colormaps["PuOr"].copy()
    cmap.set_bad("#969696")
    data=np.ma.masked_invalid(values)
    h,w=data.shape
    ax.pcolormesh(np.arange(w+1)-.5,np.arange(h+1)-.5,data,cmap=cmap,norm=norm,
                  edgecolors="black",linewidth=.25,antialiased=False,gid=series)
    ax.set(xlim=(-.5,w-.5),ylim=(h-.5,-.5),xticks=range(w),xticklabels=column_labels,
           yticks=range(h),yticklabels=row_labels)
    ax.tick_params(length=1.5,labelsize=5,top=True,right=True,pad=1)
    ax.xaxis.set_minor_locator(NullLocator())
    ax.yaxis.set_minor_locator(NullLocator())
    return cmap


def _horizontal_key(ax, norm, cmap, text, bounds):
    key=ax.inset_axes(bounds)
    bar=ax.figure.colorbar(ScalarMappable(norm=norm,cmap=cmap),cax=key,orientation="horizontal")
    step = 2. if norm.vmax>2 else .5
    bar.set_ticks(np.arange(-2,3)*step)
    key.xaxis.set_major_formatter(FormatStrFormatter("%d" if step>1 else "%.1f"))
    key.xaxis.set_minor_locator(MultipleLocator(step/2))
    key.tick_params(which="major",direction="in",labelsize=6,pad=1.75,length=1.5,top=True,bottom=True)
    key.tick_params(which="minor",direction="in",length=1,top=True,bottom=True)
    key.set_xlabel(text,fontsize=6,labelpad=.5)


def _storage_maps(axes, storage):
    ph=storage[storage.assay.eq("pH") & storage.condition.eq("thymol_blue")]
    level_cols=["row_group","column"]
    ph_levels=ph.groupby(level_cols).predictor.mean().sort_values()
    days=ph.groupby("read_index").days.mean().sort_index()
    matrix=ph.groupby(["signal",*level_cols,"read_index"]).log2_fold.mean()
    limit=matrix.abs().max()
    norm=Normalize(-limit,limit)
    for ax,sig,expr in zip(axes,SIGNALS,(r"$\log_2(R_t/R_0)$",r"$\log_2(R_t^{*}/R_0^{*})$")):
        values=matrix.loc[sig].unstack("read_index").reindex(index=ph_levels.index,columns=days.index)
        cmap=_heatmap(ax,values.to_numpy(),[f"{p:.1f}" for p in ph_levels],
                      [f"{d:.1f}" for d in days],norm,
                      f"{sig} | mean within-plate log2 ratio to initial")
        ax.tick_params(axis="x",labelrotation=90,labelsize=6)
        ax.tick_params(axis="y",pad=2)
        ax.set(xlabel="incubation time (d)",ylabel="initial pH")
        _horizontal_key(ax,norm,cmap,expr,[0,-.25238,1,.040833])


def _edge_panel(ax,tables):
    edge=tables["edge_replicates"]
    keys=["assay","physical_plate_id","condition","signal"]
    span=tables["storage_folds"][keys+["initial_response_range"]].drop_duplicates()
    ph=edge[edge.assay.eq("pH") & edge.condition.eq("thymol_blue")].merge(span,on=keys,validate="many_to_one")
    ph["fraction"]=100*ph.edge_minus_interior/ph.initial_response_range.where(ph.initial_response_range.gt(0))
    unit=ph.groupby(["signal","physical_plate_id","read_index"],as_index=False).agg(
        drift=("fraction","mean"),days=("days","mean"))
    handles=[]
    for sig,name,color,marker,opened in zip(SIGNALS,(R,RS),PH,("o","s"),(False,True)):
        sub=unit[unit.signal.eq(sig)]
        _band(ax,sub,"read_index","days","drift",color,marker,opened=opened,
              series=f"{sig} | matched edge minus interior, % of initial span")
        handles.append(_handle(color,marker,name,opened=opened))
    ax.axhline(0,color=".5",ls=":",lw=.5,zorder=0,gid="zero drift reference")
    ax.set(xlabel="time (d)",xlim=(0,16),xticks=np.arange(0,17,4),ylim=(-400,400),
           ylabel=r"$\frac{(R_{\mathrm{e},t}-R_{\mathrm{e},0})-(R_{\mathrm{i},t}-R_{\mathrm{i},0})}{\operatorname{span}(R_0)}$ (%)",
           yticks=np.arange(-400,401,200))
    ax.yaxis.set_major_formatter(FuncFormatter(lambda value, _: f"{value/100:.0f}"))
    ax.text(.02,.985,r"$\times10^{2}$", transform=ax.transAxes, ha="left",va="top",fontsize=5)
    _list_key(ax,handles,x=8.16,y=85.25,text_x=12.26)


def _temperature_panels(axes,levels):
    ph=levels[levels.assay.eq("pH") & levels.step_order.between(1,4) & levels.condition.eq("thymol_blue")].copy()
    ph["ph"]=ph.groupby(["plate_type","step_order","signal","level_id"]).predictor.transform("mean")
    temperatures=(23,4,30,75)
    for ax,sig,name in zip(axes,SIGNALS,(R,RS)):
        rows=[]
        for step,color in enumerate(PH,1):
            row=[]
            for plate,ls,marker in zip(PLATES,("-","--"),("o","s")):
                sub=ph[ph.signal.eq(sig) & ph.step_order.eq(step) & ph.plate_type.eq(plate)]
                _band(ax,sub,"level_id","ph","response_mean",color,marker,ls,opened=sig==SIGNALS[1],
                      series=f"{sig} | {temperatures[step-1]} °C | {plate.replace('_',' ')}", point_x="predictor")
                row.append(_handle(color,marker," ",ls,opened=sig==SIGNALS[1]))
            rows.append(row)
        _matrix_key(ax,[str(v) for v in temperatures],["−","+"],[row[::-1] for row in rows],
                    xs=(6.5,15.1),row_y=16.76,label_x=19.7,column_y=10.9,
                    heading=(10.8,5.15,"UV","center"),
                    side_heading=(29.02,22.4,r"$T_{\mathrm{incubate}}$ (°C)"))
        ax.set(xlabel="measured pH",ylabel=name)
        ax.set_ylim((.33325,2.88745) if sig == "ratio_raw" else (-.57123,13.28184))
        ax.yaxis.set_major_locator(MultipleLocator(.8 if sig == "ratio_raw" else 4))
        ax.set_xlim(7,0)
        ax.xaxis.set_major_locator(MultipleLocator(1))
        ax.xaxis.set_major_formatter(FormatStrFormatter("%.1f"))


def _handling_maps(axes,folds):
    # Average standard-specific contrasts within each plate, then across independent plates.
    plates=folds.groupby(["signal","plate_type","step_order","physical_plate_id"]).log2_fold.mean()
    means=plates.groupby(level=["signal","plate_type","step_order"]).mean()
    limit=folds.groupby(["signal","plate_type","step_order"]).log2_fold.mean().abs().max()
    norm=Normalize(-limit,limit)
    for ax,sig in zip(axes,SIGNALS):
        matrix=means.loc[sig].unstack("step_order").reindex(index=PLATES,columns=range(1,7))
        cmap=_heatmap(ax,matrix.to_numpy(),["+","−"],["23","4","30","75","−","+"],norm,
                      f"{sig} | mean log2 ratio to measured UV comparator")
        ax.set_ylabel("UV transparency")
        ax.tick_params(axis="x",labelsize=6,pad=2)
        ax.tick_params(axis="y",labelsize=7,pad=2)
        ax.axvline(3.5,color="black",lw=1.5)
        ax.axhline(.5,color="black",lw=1.5)
        ax.text(1.5,-.09,r"$T_{\mathrm{incubate}}$ (°C)",transform=ax.get_xaxis_transform(),
                ha="center",va="top",fontsize=7)
        ax.text(4.5,-.09,"membrane",transform=ax.get_xaxis_transform(),
                ha="center",va="top",fontsize=7)
        if sig=="ratio_raw":
            for x in (0,4):
                ax.annotate("*",(x,0),xytext=(0,-1.5),textcoords="offset points",
                            ha="center",va="center",fontsize=8,fontweight="bold")
        ratio=r"$\log_2(R/R_{\mathrm{UV,ref}})$" if sig=="ratio_raw" else r"$\log_2(R^{*}/R_{\mathrm{UV,ref}})$"
        _horizontal_key(ax,norm,cmap,ratio,[0,-.215,1,.040833])


def figure_2(tables):
    fig,axes=_canvas([
        r"$\mathrm{pH}_{\mathrm{app.,(corr.)}}$ = f(lid, pH)",
        r"$\mathrm{pH}_{\mathrm{app.}}$ = f(pH, $t$)",
        r"$\mathrm{pH}_{\mathrm{app.,corr.}}$ = f(pH, $t$)",
        r"$\mathrm{pH}_{\mathrm{app.,(corr.)}}$ = f(well position, $t$)",
        r"$\mathrm{pH}_{\mathrm{app.}}$ = f(UV, pH)",
        r"$\mathrm{pH}_{\mathrm{app.,corr.}}$ = f(UV, pH)",
        r"$\mathrm{pH}_{\mathrm{app.}}$ detectability =" + "\n" + r"f(UV, $T_{\mathrm{incubate}}$, membrane)",
        r"$\mathrm{pH}_{\mathrm{app.,corr.}}$ detectability =" + "\n" + r"f(UV, $T_{\mathrm{incubate}}$, membrane)",
    ],width=512.94,height=312.48,row_pitch=168.248,
        title_pads=(4.57,4.57,4.57,4.57,13.47,13.47,4.27,4.27))
    levels=tables["pub_robustness_levels"]
    _lid_panel(axes[0],levels)
    _storage_maps(axes[1:3],tables["storage_folds"])
    _edge_panel(axes[3],tables)
    _temperature_panels(axes[4:6],levels)
    _handling_maps(axes[6:],tables["ph_handling_folds"])
    return fig
