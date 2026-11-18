"""Microplate and shake-flask iron responses and apparent cell density."""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import AutoMinorLocator, FixedLocator, FormatStrFormatter, FuncFormatter, LogFormatterMathtext, NullFormatter
from matplotlib.transforms import Affine2D
from scipy.stats import gaussian_kde


IRON = ('#D55E00', '#CC79A7', '#E69F00', '#882255',
        '#AA4499', '#332288', '#999933', '#000000')
OD = ('#AA4499', '#332288', '#CC79A7', '#882255', '#44AA99')
MARKERS = ('o', 's', '^', 'v', 'D', 'P', 'X', '>')
OD_READOUTS = (('raw_A600', r'plate $A_{600}$'),
               ('A600_minus_A750', r'plate $A^*_{600}$'),
               ('flask_OD600_apparent', r'flask $\mathrm{OD}_{600,\mathrm{app.}}$'))
PANEL_SIZE = 97.2  # PDF points


def _panels(titles, title_offsets, height, top):
    """Four 1.35-inch square data frames separated by 0.45-inch gutters."""
    width = 510.12
    fig = plt.figure(figsize=(width / 72, height / 72))
    axes = []
    for i, (letter, title, offset) in enumerate(zip('abcd', titles, title_offsets)):
        left = 23.8545 + i * 129.6
        ax = fig.add_axes([left / width, (height - top - PANEL_SIZE) / height,
                           PANEL_SIZE / width, PANEL_SIZE / height])
        ax.set_title(title, loc='center', fontsize=7, pad=offset, y=1,
                     va='baseline', linespacing=1.15)
        ax.annotate(letter, xy=(0, 1), xycoords='axes fraction',
                    xytext=(-24.09, 2.83), textcoords='offset points',
                    va='baseline', ha='left', fontsize=8, fontweight='bold')
        ax.tick_params(which='both', direction='in', top=True, right=True, pad=2)
        ax.tick_params(which='major', length=3)
        ax.tick_params(which='minor', length=1.6)
        ax.xaxis.labelpad = 2
        ax.yaxis.labelpad = 1
        axes.append(ax)
    return fig, axes


def _key_coordinates(ax):
    """Locate legend text and samples in points from the upper left of the data frame."""
    return Affine2D().scale(1 / PANEL_SIZE, -1 / PANEL_SIZE).translate(0, 1) + ax.transAxes


def _sample(ax, x, y, color, marker, *, width=7, style='-', opened=False, edge_width=1):
    # Separate line and marker so dashes do not depend on the label or marker width.
    transform = _key_coordinates(ax)
    line = Line2D([x, x + width], [y, y], color=color, lw=.8, ls=style,
                  transform=transform, zorder=10)
    ax.add_artist(line)
    ax.add_artist(Line2D([x + width / 2], [y], color=color, marker=marker,
                         markerfacecolor='none' if opened else color, ms=3,
                         markeredgewidth=edge_width, linestyle='none',
                         transform=transform, zorder=10))


def _key_text(ax, x, baseline, label, **kwargs):
    return ax.text(x, baseline, label, transform=_key_coordinates(ax),
                   fontsize=5, va='baseline', zorder=10, **kwargs)


def _time_axis(ax):
    # Retain the common daily tick interval through the last observed vessel.
    last = max(7., ax.dataLim.xmax)
    ax.set(xlim=(0, last + .04 if last > 7 else 7), xticks=np.arange(np.floor(last)+1), xlabel='time (d)')
    ax.xaxis.set_minor_locator(AutoMinorLocator(2))
    ax.yaxis.set_minor_locator(AutoMinorLocator(2))


def _points(ax, frame, x, y, color, marker, series, *, opened=False, edge_width=.3):
    valid = frame.loc[frame[x].notna() & frame[y].notna()]
    if valid.empty:
        return
    points = ax.scatter(valid[x], valid[y], s=9, marker=marker, edgecolor=color,
                        facecolor='none' if opened else color, linewidth=edge_width,
                        zorder=3, gid=series + ' | vessel observations')
    id_column = next(col for col in ('unit_id', 'flask_ID', 'plate_id') if col in valid)
    points.plotted_value_labels = valid[id_column].astype(str).tolist()


def _vessel_series(ax, frame, x, y, group, color, marker, series, style='-', opened=False):
    """One point per vessel and one mean and sample SD per recorded acquisition."""
    if not frame[y].notna().any():
        return
    _points(ax, frame, x, y, color, marker, series, opened=opened)
    observations = frame.assign(plot_x=frame[x].where(frame[y].notna()))
    summary = observations.groupby(group, sort=True).agg(
        plot_x=('plot_x', 'mean'), average=(y, 'mean'), sd=(y, 'std'))
    ax.plot(summary.plot_x, summary.average, color=color, ls=style, lw=.8,
            gid=series + ' | within-batch mean')
    ax.fill_between(summary.plot_x, summary.average-summary.sd, summary.average+summary.sd,
                    color=color, lw=0, alpha=.18, gid=series + ' | within-batch mean ± SD')


def _batch_trajectories(ax, frame, batch, value, event, palette, prefix, key_position):
    key_x, first_y = key_position
    for i, (batch_id, sample) in enumerate(frame.groupby(batch, sort=True)):
        color, marker = palette[i % len(palette)], MARKERS[i % len(MARKERS)]
        label = f'{prefix}{i+1}'
        _vessel_series(ax, sample, 'days_mean', value, event, color, marker, f'{label} ({batch_id})')
        y = first_y + 5.6 * i
        _sample(ax, key_x, y, color, marker)
        _key_text(ax, key_x + 7.6, y + 1.75, label)
    _time_axis(ax)


def _control_flasks(tables):
    f = tables['validation_flask_points']
    return f.loc[f.no_addition & f.kind.eq('biotic') & f.organism.eq('AF') & f.medium.eq('AFM1')]


def _repeatability(ax, tables, metric):
    baseline = tables['bio_baseline_repeatability'].set_index('metric').loc[metric]
    ga = tables['bio_flask_repeatability']
    ga = ga.loc[ga.metric.eq(metric) & ga.numerator_Ga_mM.eq(0)].set_index('denominator_Ga_mM')
    conditions = (baseline.plate_condition_id, baseline.flask_condition_id,
                  ga.loc[1., 'denominator_condition_id'], ga.loc[100., 'denominator_condition_id'])
    batch = tables['condition_repeatability_batches']
    categories = ('AFM1 microplates, no additive', 'AFM1 flasks, no additive',
                  'AFM1 flasks, 1 mM Ga³⁺', 'AFM1 flasks, 100 mM Ga³⁺')
    for pos, condition in enumerate(conditions):
        q = batch.loc[batch.metric.eq(metric) & batch.condition_id.eq(condition)
                      & batch.eligible & batch.deviation_percent.notna()].sort_values('batch_id')
        values = q.deviation_percent.to_numpy()
        if not len(values):
            continue
        color = IRON[0] if pos == 0 else '#444444'
        density = None
        if len(values) > 1 and np.ptp(values) > 0:
            density = gaussian_kde(values, bw_method='scott')
            span = np.linspace(values.min(), values.max(), 256)
            heights = density(span)
            scale = .36 / heights.max()
            widths = heights * scale
            ax.fill_betweenx(span, pos-widths, pos+widths, fc=color, ec=color,
                             lw=.55, alpha=.18, zorder=1)
            for side in (-1, 1):
                ax.plot(pos + side * widths, span, color=color, lw=.55, zorder=2,
                        gid=categories[pos] + ' | density outline')
        quartiles = np.quantile(values, [.25, .5, .75])
        for y, role, style, lw in zip(quartiles, ('lower quartile', 'median', 'upper quartile'),
                                      (':', '-', ':'), (.55, .9, .55)):
            width = float(density(y)[0] * scale) if density is not None else .13
            ax.hlines(y, pos-width, pos+width, color=color, lw=lw, ls=style, zorder=3,
                      gid=categories[pos] + ' | ' + role)
        has_iodate = q.source_context_ids.str.contains(r'potassium_iodate(?:_|/)', na=False)
        offsets = np.resize([-.15, .05, .15, -.05, -.10, .10, -.20, .20], len(q))
        for is_iodate, marker in ((False, 'o'), (True, 'D')):
            keep = has_iodate.to_numpy() == is_iodate
            if keep.any():
                points = ax.scatter(pos+offsets[keep], values[keep], s=9, c=color,
                                    marker=marker, lw=0, zorder=4,
                                    gid=categories[pos] + ' | batch deviations')
                points.plotted_value_labels = q.loc[keep, 'batch_id'].tolist()
    ax.axhline(0, ls=':', lw=.5, color='#777777', zorder=0)
    ax.axvline(1.5, lw=.7, color='black')
    ax.set_xlim(-.5, 3.5)
    ax.set_xlabel('[Ga] (M)', fontsize=6)
    ax.set_xticks(range(4), ['plates\n0', 'flasks\n0', 'flasks\n$10^{-3}$', 'flasks\n$10^{-1}$'])
    ax.tick_params(axis='x', labelsize=5)
    ax.xaxis.labelpad = 3
    for tick in ax.get_xticklabels():
        tick.set_linespacing(1.0)
    if metric == 'AUC_to_own_max':
        ax.set(ylim=(-74.883, 127.983), yticks=[-50, 0, 50, 100], ylabel='batch AUC deviation (%)')
    else:
        ax.set(ylim=(-24.221, 34.421), yticks=range(-20, 31, 10), ylabel='batch timing deviation (%)')


def figure_5(tables):
    ferric = r'\mathrm{Fe}^{3+}_{\mathrm{app.,AF}}'
    titles = (rf'${ferric}$ in plates = f(batch, t)',
              rf'${ferric}$ in flasks = (batch, t)',
              'plate-to-flask\n' + rf'$\mathrm{{AUC}}_{{{ferric}}}$ repeatability',
              'plate-to-flask\n' + rf'$t_{{\mathrm{{max,}}{ferric}}}$ repeatability')
    fig, (a, b, c, d) = _panels(titles, (14.31, 14.31, 5.72, 5.72), height=139.8, top=21.3162)
    baseline = tables['bio_baseline_repeatability'].iloc[0]
    plates = tables['af_iron_control_plate_trajectories']
    plates = plates.loc[plates.signal.eq('raw_A380') & plates.plate_id.isin(baseline.plate_unit_ids.split(';'))]
    if plates.duplicated(['plate_id', 'ordinal']).any():
        raise ValueError('No-additive iron control contains duplicate measurements per plate and acquisition')
    _batch_trajectories(a, plates, 'batch_id', 'biotic_change', 'ordinal', IRON, 'P', (4.90, 6.73))
    a.set(ylim=(-.7655, max(3.6942, a.dataLim.ymax + .08)), yticks=np.arange(-.5, 4., .5), ylabel=r'$\Delta A_{380,\mathrm{AF}}$ (AU)')
    a.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    flasks = _control_flasks(tables)
    _batch_trajectories(b, flasks.assign(ferric_M=flasks.biotic_change / 1000), 'batch', 'ferric_M',
                        'sampling_event', IRON, 'F', (80.06, 68.45))
    b.set(ylim=(min(0, b.dataLim.ymin - .001), .08), yticks=np.arange(0, .081, .02),
          ylabel=r'$[\mathrm{Fe}^{3+}]_{\mathrm{app.}}$ (M)')
    b.yaxis.set_major_formatter(FuncFormatter(lambda value, _: f'{100 * value:.0f}'))
    _key_text(b, 1.94, 6.38, r'$\times 10^{-2}$')
    _repeatability(c, tables, 'AUC_to_own_max')
    _repeatability(d, tables, 'time_to_observed_max')
    return fig


def _paired_controls(ax, tables):
    controls = tables['od_plate_trajectories']
    controls = controls.loc[controls.Ga_mM.eq(0)]
    inset = ax.inset_axes([.49, .1229, .47, .47])
    _key_text(ax, 3.91, 7.05, 'cell condition')
    for column, cell in enumerate(('Ab', 'AF')):
        _key_text(ax, 14.06 + 8.6 * column, 12.61, cell, ha='center')
    for i, (signal, _) in enumerate(OD_READOUTS[:2]):
        q = controls.loc[controls.signal.eq(signal)]
        if q.duplicated(['batch_id', 'plate_id', 'ordinal']).any():
            raise ValueError('Paired A600 control contains duplicate records for a physical plate')
        y = 16.92 + 5.6 * i
        _key_text(ax, 27.26, y + 1.1, r'$A_{600}$' if i == 0 else r'$A^*_{600}$')
        for j, (cell, value, style, marker) in enumerate([('Ab', 'abiotic_signal', '--', 's'),
                                                           ('AF', 'biotic_signal', '-', 'o')]):
            _sample(ax, 10.86 + 8.6*j, y, OD[i], marker, width=6.4,
                    style=style, opened=i == 1, edge_width=.3)
            for batch_id, cohort in q.groupby('batch_id', sort=True):
                series = f'{cell}, {signal} ({batch_id})'
                _vessel_series(ax, cohort, 'days_mean', value, 'ordinal', OD[i], marker, series, style, i == 1)
                if i == 1 or cell == 'Ab':
                    _vessel_series(inset, cohort, 'days_mean', value, 'ordinal',
                                   OD[i], marker, series, style, i == 1)
    _time_axis(ax)
    ax.set(ylim=(-.024, 1.4), yticks=np.arange(0., 1.41, .2), ylabel=r'$A_{600}$ (AU)')
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    inset.set(xlim=(0, 7), ylim=(-.02, .08), xticks=[0, 7], yticks=[0, .08])
    inset.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    inset.xaxis.set_minor_locator(AutoMinorLocator(2))
    inset.yaxis.set_minor_locator(AutoMinorLocator(2))
    inset.tick_params(which='both', labelsize=5, pad=.5, direction='in', top=True, right=True)
    inset.tick_params(which='major', length=1.5)
    inset.tick_params(which='minor', length=.8)


def _od_treatment(ax, tables, value):
    units = tables['od_dose_units']
    selected = units.loc[~units.is_reference & units.Ga_mM.gt(0)].copy()
    selected['ga_M'] = selected.Ga_mM / 1000
    for i, (signal, label) in enumerate(OD_READOUTS):
        data = selected.loc[selected.signal.eq(signal)].dropna(subset=[value])
        color, marker = OD[i], MARKERS[i]
        for batch_id, batch in data.groupby('batch_id', sort=True):
            series = f'{signal} ({batch_id})'
            _points(ax, batch, 'ga_M', value, color, marker, series, opened=i == 1, edge_width=.4)
            summary = batch.groupby('ga_M', sort=True)[value].agg(['mean', 'std'])
            mean, sd = summary['mean'], summary['std']
            ax.plot(mean.index, mean, color=color, lw=.8, gid=series + ' | within-batch mean')
            ax.fill_between(mean.index, mean-sd, mean+sd, lw=0, color=color, alpha=.12,
                            gid=series + ' | within-batch mean ± SD')
        y = 8.47 + 5.6 * i
        _sample(ax, 4.66, y, color, marker, opened=i == 1, edge_width=.4)
        _key_text(ax, 12.26, y + 1.5, label)
    timing = value == 'delay_d'
    ax.axhline(0 if timing else 100, color='#777777', lw=.5, ls=':', zorder=0)
    ax.set_xscale('log')
    ax.margins(x=.05)
    ax.set_xlabel('[Ga] (M)', fontsize=6)
    xmin, xmax = ax.get_xlim()
    decades = np.arange(np.ceil(np.log10(xmin)), 1 + np.floor(np.log10(xmax)))
    ax.set_xticks(10. ** decades)
    ax.xaxis.set_major_formatter(LogFormatterMathtext())
    # Nine intermediate positions evenly spaced in logarithmic display coordinates.
    minor = np.arange(np.ceil(10*np.log10(xmin)), 1 + np.floor(10*np.log10(xmax)))
    ax.xaxis.set_minor_locator(FixedLocator(10. ** (minor / 10)))
    ax.xaxis.set_minor_formatter(NullFormatter())
    ax.yaxis.set_minor_locator(AutoMinorLocator(2))
    if timing:
        ymin, ymax = ax.get_ylim()
        ax.set_yticks(np.arange(np.ceil(ymin), 1 + np.floor(ymax)))
        ax.set_ylabel('delay to maximum (d)')
    else:
        ax.set(ylim=(-18.526, max(460.439, ax.dataLim.ymax + 5)), yticks=range(0, 401, 100),
               ylabel='integrated response\n(% of control)')


def figure_6(tables):
    density = r'planktonic cell density ($\mathrm{OD}_{600}$)'
    titles = ('$A_{600}$ and $A^{*}_{600}$\nin plates = f(cells, t)',
              r'$\mathrm{OD}_{\mathrm{600,app.,AF}}$ in flasks = f(t)',
              'relative AUC: $A_{600}$, $A^{*}_{600}$, apparent\n' + density,
              r'$\Delta t_{\mathrm{max}}$: $A_{600}$, $A^{*}_{600}$, apparent' + '\n' + density)
    fig, (a, b, c, d) = _panels(titles, (2.83, 11.94, 1.92, 1.92), height=131.82, top=17.9472)
    _paired_controls(a, tables)
    _batch_trajectories(b, _control_flasks(tables), 'batch', 'calculated_planktonic_OD600_equivalent',
                        'sampling_event', OD, 'F', (4.75, 5.17))
    b.set(ylim=(min(0, b.dataLim.ymin - .001), .12), yticks=np.arange(0, .121, .03),
          ylabel=r'planktonic $\mathrm{OD}_{600,\mathrm{app.}}$')
    b.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    _od_treatment(c, tables, 'response_percent')
    _od_treatment(d, tables, 'delay_d')
    return fig
