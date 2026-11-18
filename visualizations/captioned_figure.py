"""Write a one-page vector PDF containing the plot and its manuscript caption."""

from functools import lru_cache
import re

from matplotlib import rcParams
from matplotlib.backends.backend_agg import RendererAgg
from matplotlib.font_manager import FontProperties
from matplotlib.transforms import Bbox


INLINE = re.compile(r"(\$[^$]+\$|\*\*[^*]+\*\*|\*[^*]+\*)")


def _caption_words(text):
    """Keep inline formulae and species names intact when wrapping lines."""
    word = []
    for part in INLINE.split(text):
        if part.startswith("$"):
            word.append(("math", part))
        elif part.startswith("**"):
            word.append(("bold", part[2:-2]))
        elif part.startswith("*"):
            word.append(("italic", part[1:-1]))
        else:
            for piece in re.split(r"(\s+)", part):
                if piece.isspace():
                    if word:
                        yield word
                        word = []
                elif piece:
                    word.append(("text", piece))
    if word:
        yield word


def _caption_text(path):
    paragraphs = path.read_text(encoding="utf-8").strip().split("\n\n")
    heading = paragraphs[0].removeprefix("# ")
    caption = []
    for paragraph in paragraphs[1:]:
        if paragraph.startswith("[Figure (PDF)]"):
            break
        if not paragraph.startswith("!["):
            caption.append(paragraph.strip())
    if not heading.startswith("Figure ") or not any(re.match(r"\*\*\(a(?:[–,][a-l])*\)\*\*", p) for p in caption):
        raise ValueError(f"Missing manuscript title or first panel caption in {path}")
    return heading, " ".join(caption)


def artwork_bounds(figure):
    """Keep the plot scale and allow room for all axis labels and scripts."""
    figure.canvas.draw()
    labels = figure.get_tightbbox(figure.canvas.get_renderer())
    return Bbox.union([figure.bbox_inches, labels]).padded(1 / 72)


def save_captioned_pdf(figure, caption_path, destination, bounds=None):
    """Print the complete caption below the lowest part of the artwork."""
    heading, caption = _caption_text(caption_path)
    # Work in physical points so line widths and margins do not depend on the
    # raster resolution of figure previews or on individual figure dimensions.
    width, height = figure.get_size_inches()
    if bounds is None:
        bounds = artwork_bounds(figure)
    margin = 14.
    available = bounds.width * 72 - 2 * margin
    renderer = RendererAgg(1, 1, 720)
    caption_artists = []

    @lru_cache(maxsize=None)
    def font(style, is_title):
        return FontProperties(family=[*rcParams["font.family"], "DejaVu Sans"],
                              size=8.5 if is_title else 7.8,
                              weight="bold" if is_title or style == "bold" else "normal",
                              style="italic" if style == "italic" else "normal")

    @lru_cache(maxsize=None)
    def metrics(style, text, is_title):
        raw = renderer.get_text_width_height_descent(text, font(style, is_title), style == "math")
        return tuple(value / 10 for value in raw)

    def wrap(text, is_title):
        words, used = [], 0.
        space = metrics("text", " ", is_title)[0]
        for word in _caption_words(text):
            length = sum(metrics(style, value, is_title)[0] for style, value in word)
            if words and used + space + length > available:
                yield words
                words, used = [], 0.
            used += (space if words else 0) + length
            words.append(word)
        if words:
            yield words

    depth = 4. - 72 * bounds.y0
    try:
        for text, is_title in ((heading, True), (caption, False)):
            for words in wrap(text, is_title):
                runs = []
                for i, word in enumerate(words):
                    for style, content in ([("text", " ")] if i else []) + word:
                        if runs and runs[-1][0] == style and style != "math":
                            runs[-1][1] += content
                        else:
                            runs.append([style, content])
                extents = [metrics(style, content, is_title) for style, content in runs]
                ascent = max(h - desc for _, h, desc in extents)
                descent = max(desc for _, _, desc in extents)
                baseline = depth + ascent
                offset = 72 * bounds.x0 + margin
                for (style, content), (w, _, _) in zip(runs, extents):
                    artist = figure.text(offset / (72*width), -baseline / (72*height), content,
                                         fontproperties=font(style, is_title), va="baseline",
                                         ha="left", parse_math=style == "math")
                    caption_artists.append(artist)
                    offset += w
                depth = baseline + descent + 2.9
            if is_title:
                depth += 1.8
        page = Bbox.from_extents(bounds.x0, -(depth + 7) / 72, bounds.x1, bounds.y1)
        figure.savefig(destination, format="pdf", bbox_inches=page, pad_inches=0,
                       metadata={"CreationDate": None, "ModDate": None})
    finally:
        for artist in caption_artists:
            artist.remove()
