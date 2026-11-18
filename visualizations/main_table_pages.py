"""Render manuscript overview tables and their component-level source links as Markdown."""
from html import escape
import re
from urllib.parse import unquote, urlsplit


SOURCE_FIELDS = {
    "Experimental design URL": "Experimental design",
    "Readout definition URL": "Readout",
    "Preparation definition URL": "Preparation",
    "Standard and reference URLs": "Standards and references",
    "Measurement URLs": "Measurements",
    "Well assignments URL": "Well assignments",
    "Condition time-course records": "Condition time-course records",
    "Acquisition records": "Acquisition registry",
    "Raw measurement directories": "Raw measurement directories",
    "Example raw acquisition URL": "Example raw acquisition",
    "Measurement definitions URL": "Measurement definitions",
    "Measurement methods URL": "Measurement methods",
    "Medium recipe URL": "Medium recipe",
    "Recorded observation URL": "Recorded observation",
}


def text_markup(text):
    """Format scientific names and subscripts, for the table cells."""
    result = escape(str(text), quote=False).replace("|", "\\|")
    result = re.sub(r"\b([nN])_(prep|well|microplate)(?==)", r"\1<sub>\2</sub>", result)
    result = re.sub(r"\b([AL]\. (?:ferrooxidans|ferridurans|ferriphilum|thiooxidans|caldus))\b",
                    r"*\1*", result)
    return result.replace("Aλ", "A<sub>λ</sub>").replace("R*", "R<sup>*</sup>")


def link_label(url):
    part = urlsplit(url)
    name = unquote(part.path.rsplit("/", 1)[-1])
    if not part.fragment:
        return name
    lines = part.fragment.removeprefix("L").replace("-L", "–")
    return f"{name}, {lines}"


def group_anchor(key):
    return "sources-for-" + key.replace(":", "-").lower()


def group_name(key):
    if ":" in key:
        group, temperature = key.split(":")
        return f"{group} ({temperature} °C)"
    return key


def render_page(number, definition, headings, components, destination, data_href):
    out = [
        f"# Table {number}. {text_markup(definition['title'])}",
        "",
        text_markup(definition["note"]),
        "",
        "| " + " | ".join(headings) + " |",
        "| " + " | ".join("---" for _ in headings) + " |",
    ]
    for group, cells in definition["rows"].items():
        label = f"[{cells[0]}](#{group_anchor(group)})"
        contents = [text_markup(cell).replace("\n", "<br><br>") for cell in cells[1:]]
        out.append("| " + " | ".join([label, *contents]) + " |")
    out.extend(["", "[Main-text tables](README.md)", ""])
    for group in definition["rows"]:
        out.extend([f'<a name="{group_anchor(group)}"></a>', "",
                    f"## Sources for {group_name(group)}", "",
                    "| Component | Acquisition or condition | Field | Records or values |",
                    "|---|---|---|---|"])
        for row in components[group]:
            component = text_markup(row["Component"])
            condition = text_markup(row["Acquisition or condition"])
            for field, value in row.items():
                if not value or field in {"Component", "Acquisition or condition"}:
                    continue
                if field in SOURCE_FIELDS:
                    records = [f"[{link_label(ref)}]({data_href}/{ref})"
                               for ref in value.splitlines()]
                else:
                    records = [text_markup(v) for v in value.splitlines()]
                label = SOURCE_FIELDS.get(field, field)
                out.append(f"| {component} | {condition} | {label} | {'<br>'.join(records)} |")
        out.append("")
    (destination / f"Table_{number}.md").write_text("\n".join(out), encoding="utf-8")
