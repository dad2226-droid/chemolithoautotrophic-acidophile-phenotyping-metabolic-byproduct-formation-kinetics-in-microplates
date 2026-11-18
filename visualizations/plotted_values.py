"""Write data coordinates, series identities and uncertainty intervals for each panel."""

import csv
from string import ascii_lowercase

import numpy as np
from matplotlib.collections import FillBetweenPolyCollection, LineCollection, PathCollection, QuadMesh
from matplotlib.container import ErrorbarContainer
from matplotlib.lines import Line2D


FIELDS = ('figure', 'panel', 'view', 'series', 'mark', 'group', 'index', 'observation',
          'x_axis', 'y_axis', 'x', 'y', 'lower_y', 'upper_y', 'x_category', 'y_category', 'value',
          'in_view')
CATEGORICAL = {(4, 'd'), (5, 'c'), (5, 'd'), (7, 'i'), (8, 'h')}


def _array(values):
    return np.ma.asarray(values, dtype=float).filled(np.nan)


def _series(artist):
    name = artist.get_gid() or artist.get_label()
    return name if isinstance(name, str) and not name.startswith('_') else ''


def _text(value):
    return str(value).replace('\n', ' ').strip()


def _tick_labels(axis):
    return {float(pos): _text(label.get_text())
            for pos, label in zip(axis.get_ticklocs(), axis.get_ticklabels())}


def _colorbar(ax):
    return ax.get_label() in ('colorbar', '<colorbar>') or getattr(ax, '_colorbar', None) is not None


def _line_rows(ax, line, categorical, heatmap):
    xs, ys = _array(line.get_xdata(orig=False)), _array(line.get_ydata(orig=False))
    if not len(xs):
        return
    transform = line.get_transform()
    if transform == ax.get_yaxis_transform():
        if not heatmap:
            yield dict(mark='horizontal reference', y=ys[0])
    elif transform == ax.get_xaxis_transform():
        if not (heatmap or categorical):
            yield dict(mark='vertical reference', x=xs[0])
    elif transform == ax.transData:
        role = 'identity reference' if len(xs) == 2 and np.array_equal(xs, ys) else 'line'
        for idx, (x, y) in enumerate(zip(xs, ys)):
            if np.isfinite(x) or np.isfinite(y):
                yield dict(mark=role, index=idx + 1, x=x, y=y)


def _point_rows(artist):
    values = artist.get_array()
    scalar = _array(values).reshape(-1) if values is not None else None
    offsets = _array(artist.get_offsets())
    labels = getattr(artist, 'plotted_value_labels', None)
    if labels is not None and len(labels) != len(offsets):
        raise ValueError(f'Unequal observation IDs and points for {_series(artist)}')
    for idx, (x, y) in enumerate(offsets):
        if np.isfinite(x) and np.isfinite(y):
            yield dict(mark='point', index=idx + 1, x=x, y=y,
                       observation=labels[idx] if labels is not None else '',
                       value=scalar[idx] if scalar is not None else '')


def _band_rows(artist):
    # Vertical violins describe a density. Their underlying observations and
    # quantile marks are exported directly from their point and line artists.
    if artist.t_direction != 'x':
        return
    for vertices in artist.get_paths():
        points = vertices.vertices
        # Matplotlib closes each contiguous interval polygon with its starting
        # point. Its two data boundaries are traversed in opposite directions.
        n = (len(points) - 3) // 2
        first = points[1:1 + n]
        second = points[n + 2:2 * n + 2][::-1]
        if len(points) != 2 * n + 3 or not np.array_equal(first[:, 0], second[:, 0]):
            raise ValueError(f'Unsupported interval polygon for {_series(artist)}')
        for i, ((x, low), (_, high)) in enumerate(zip(first, second), 1):
            yield dict(mark='interval band', index=i, x=x, lower_y=min(low, high),
                       upper_y=max(low, high))


def _mesh_rows(ax, mesh):
    data = _array(mesh.get_array())
    corners = mesh.get_coordinates()
    height, width = np.subtract(corners.shape[:2], 1)
    data = data.reshape(height, width)
    x_labels = getattr(mesh, '_plot_data_x_labels', None)
    y_labels = getattr(mesh, '_plot_data_y_labels', None)
    tick_x, tick_y = _tick_labels(ax.xaxis), _tick_labels(ax.yaxis)
    for r in range(height):
        for c in range(width):
            x, y = corners[r:r+2, c:c+2].mean(axis=(0, 1))
            yield dict(mark='heatmap cell', index=r * width + c + 1, x=x, y=y,
                       x_category=_text(x_labels[c]) if x_labels is not None else tick_x.get(x, ''),
                       y_category=_text(y_labels[r]) if y_labels is not None else tick_y.get(y, ''),
                       value=data[r, c])


def _summary_marks(lines):
    segments = lines.get_segments()
    for i, segment in enumerate(segments):
        if len(segment) != 2:
            continue
        (x0, y0), (x1, y1) = segment
        if y0 == y1:
            if len(segments) == 3:
                role = ('lower quartile', 'median', 'upper quartile')[i]
            elif len(segments) == 2 and 'quartile' in _series(lines).lower():
                role = ('lower quartile', 'upper quartile')[i]
            else:
                role = 'median' if 'median' in _series(lines).lower() else 'horizontal summary'
            yield dict(mark=role, index=i+1, x=(x0+x1)/2, y=y0)


def _uncertainty(ax):
    """Combine a point, error bars and caps into a single interval record."""
    excluded = set()
    rows = []
    for container in ax.containers:
        if not isinstance(container, ErrorbarContainer):
            continue
        point, caps, bars = container.lines
        excluded.update([point, *caps, *bars])
        xs, ys = point.get_data(orig=False)
        bar = bars[0]
        for i, (x, y, limits) in enumerate(zip(xs, ys, bar.get_segments()), 1):
            if len(limits) and limits[0][0] == limits[1][0]:
                rows.append(dict(mark='point with interval', series=_series(point) or _series(bar),
                                 index=i, x=x, y=y, lower_y=limits[:, 1].min(),
                                 upper_y=limits[:, 1].max()))
    return excluded, rows


def _axis_rows(ax, number, panel, main, view):
    categorical = (number, panel) in CATEGORICAL
    category_labels = _tick_labels(ax.xaxis) if categorical else {}
    mesh_present = any(isinstance(a, QuadMesh) for a in ax.collections)
    base = dict(figure=number, panel=panel, view=view, x_axis=_text(ax.get_xlabel() or main.get_xlabel()),
                y_axis=_text(ax.get_ylabel() or main.get_ylabel()))
    excluded, error_rows = _uncertainty(ax)
    for i, row in enumerate(error_rows, 1):
        yield {**base, 'group': f'errorbar-{i}', **row}
    count = 0
    for artist in ax.get_children():
        if artist in excluded or not artist.get_visible():
            continue
        if isinstance(artist, Line2D):
            rows = _line_rows(ax, artist, categorical, mesh_present)
        elif isinstance(artist, PathCollection):
            rows = _point_rows(artist)
        elif isinstance(artist, FillBetweenPolyCollection):
            rows = _band_rows(artist)
        elif isinstance(artist, QuadMesh):
            rows = _mesh_rows(ax, artist)
        elif isinstance(artist, LineCollection):
            rows = _summary_marks(artist)
        else:
            continue
        count += 1
        for row in rows:
            if categorical and category_labels and 'x' in row and np.isfinite(row['x']):
                level = min(category_labels, key=lambda v: abs(v-row['x']))
                if abs(level - row['x']) < .5:
                    row['x_category'] = category_labels[level]
            yield {**base, 'series': _series(artist), 'group': count, **row}


def write_plotted_values(figure, figure_number, destination):
    """Write one UTF-8 CSV containing the visible data marks and numerical guides."""
    figure.canvas.draw()
    parents = [ax for ax in figure.axes if ax.get_visible() and not _colorbar(ax)]
    with destination.open('w', encoding='utf-8', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=FIELDS, lineterminator='\n')
        writer.writeheader()
        for panel, main in zip(ascii_lowercase, parents):
            children = [ax for ax in main.child_axes if ax.get_visible() and not _colorbar(ax)]
            for i, ax in enumerate([main, *children]):
                view = f'inset {i}' if i else 'main'
                xmin, xmax = sorted(ax.get_xlim())
                ymin, ymax = sorted(ax.get_ylim())
                for row in _axis_rows(ax, figure_number, panel, main, view):
                    # Lines and intervals use full data coordinates even when an
                    # inset shows a smaller window. Record which sampled values
                    # fall inside that window, including partially visible bands.
                    x_in = 'x' not in row or xmin <= row['x'] <= xmax
                    y_in = 'y' not in row or ymin <= row['y'] <= ymax
                    if 'lower_y' in row:
                        y_in = row['upper_y'] >= ymin and row['lower_y'] <= ymax
                    row['in_view'] = bool(x_in and y_in)
                    for key, value in row.items():
                        if isinstance(value, (float, np.floating)) and not np.isfinite(value):
                            row[key] = ''
                    writer.writerow(row)
