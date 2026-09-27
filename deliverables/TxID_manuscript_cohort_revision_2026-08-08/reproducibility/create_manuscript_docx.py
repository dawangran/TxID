#!/usr/bin/env python3
"""Build a bilingual-ready TxID v1.1 manuscript as a self-contained DOCX."""

from __future__ import annotations

import argparse
import html
import re
import sys
import zipfile
from dataclasses import dataclass
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "manuscript" / "txid-gigascience-technical-note-v1.1.md"
OUTPUT = ROOT / "manuscript" / "TxID_GigaScience_Technical_Note_v1.1.docx"
ZIP_TIME = (2026, 8, 8, 0, 0, 0)
EMU_PER_INCH = 914_400
PAGE_TEXT_WIDTH_EMU = int(6.55 * EMU_PER_INCH)
W_NS = "http://schemas.openxmlformats.org/wordprocessingml/2006/main"


@dataclass(frozen=True)
class ImageAsset:
    source: Path
    archive_name: str
    relationship_id: str
    width_px: int
    height_px: int
    alt_text: str


def xml_text(value: str, *, attribute: bool = False) -> str:
    return html.escape(value, quote=attribute)


def contains_cjk(value: str) -> bool:
    return any("\u3400" <= char <= "\u9fff" for char in value)


def run(
    value: str,
    *,
    bold: bool = False,
    italic: bool = False,
    code: bool = False,
    size: int | None = None,
    color: str | None = None,
) -> str:
    properties: list[str] = []
    if bold:
        properties.append("<w:b/>")
    if italic:
        properties.append("<w:i/>")
    if code:
        properties.append(
            '<w:rFonts w:ascii="Consolas" w:hAnsi="Consolas" w:eastAsia="Microsoft YaHei"/>'
        )
        properties.append('<w:shd w:val="clear" w:color="auto" w:fill="F3F3F3"/>')
    elif contains_cjk(value):
        properties.append(
            '<w:rFonts w:ascii="Arial" w:hAnsi="Arial" w:eastAsia="Microsoft YaHei"/>'
        )
        properties.append('<w:lang w:val="en-US" w:eastAsia="zh-CN"/>')
    if size is not None:
        properties.append(f'<w:sz w:val="{size}"/><w:szCs w:val="{size}"/>')
    if color is not None:
        properties.append(f'<w:color w:val="{color}"/>')
    rpr = f"<w:rPr>{''.join(properties)}</w:rPr>" if properties else ""
    preserve = ' xml:space="preserve"' if value[:1].isspace() or value[-1:].isspace() else ""
    return f"<w:r>{rpr}<w:t{preserve}>{xml_text(value)}</w:t></w:r>"


INLINE_PATTERN = re.compile(r"(\*\*.+?\*\*|`.+?`|\*[^*]+?\*)")


def inline_runs(value: str, *, default_color: str | None = None) -> str:
    output: list[str] = []
    position = 0
    for match in INLINE_PATTERN.finditer(value):
        if match.start() > position:
            output.append(run(value[position : match.start()], color=default_color))
        token = match.group(0)
        if token.startswith("**"):
            output.append(run(token[2:-2], bold=True, color=default_color))
        elif token.startswith("`"):
            output.append(run(token[1:-1], code=True, color=default_color))
        else:
            output.append(run(token[1:-1], italic=True, color=default_color))
        position = match.end()
    if position < len(value):
        output.append(run(value[position:], color=default_color))
    return "".join(output)


def paragraph(
    content: str,
    *,
    style: str | None = None,
    align: str | None = None,
    before: int = 0,
    after: int = 120,
    line: int = 276,
    keep_next: bool = False,
    page_break_before: bool = False,
    left_indent: int | None = None,
    first_line: int | None = None,
) -> str:
    properties: list[str] = [
        f'<w:spacing w:before="{before}" w:after="{after}" '
        f'w:line="{line}" w:lineRule="auto"/>'
    ]
    if style:
        properties.append(f'<w:pStyle w:val="{style}"/>')
    if align:
        properties.append(f'<w:jc w:val="{align}"/>')
    if keep_next:
        properties.append("<w:keepNext/>")
    if page_break_before:
        properties.append("<w:pageBreakBefore/>")
    if left_indent is not None or first_line is not None:
        attrs = []
        if left_indent is not None:
            attrs.append(f'w:left="{left_indent}"')
        if first_line is not None:
            attrs.append(f'w:firstLine="{first_line}"')
        properties.append(f"<w:ind {' '.join(attrs)}/>")
    return f"<w:p><w:pPr>{''.join(properties)}</w:pPr>{content}</w:p>"


def text_paragraph(value: str, *, chinese: bool | None = None, **kwargs: object) -> str:
    is_chinese = contains_cjk(value) if chinese is None else chinese
    color = "3F3F3F" if is_chinese else None
    return paragraph(inline_runs(value, default_color=color), **kwargs)


def png_dimensions(payload: bytes, path: Path) -> tuple[int, int]:
    if payload[:8] != b"\x89PNG\r\n\x1a\n" or payload[12:16] != b"IHDR":
        raise ValueError(f"not a supported PNG: {path}")
    return int.from_bytes(payload[16:20], "big"), int.from_bytes(payload[20:24], "big")


def image_paragraph(asset: ImageAsset) -> str:
    max_height_emu = int(7.15 * EMU_PER_INCH)
    width_emu = PAGE_TEXT_WIDTH_EMU
    height_emu = round(width_emu * asset.height_px / asset.width_px)
    if height_emu > max_height_emu:
        height_emu = max_height_emu
        width_emu = round(height_emu * asset.width_px / asset.height_px)
    alt = xml_text(asset.alt_text, attribute=True)
    name = xml_text(asset.archive_name, attribute=True)
    drawing = f"""
<w:r><w:drawing><wp:inline distT="0" distB="0" distL="0" distR="0">
  <wp:extent cx="{width_emu}" cy="{height_emu}"/>
  <wp:effectExtent l="0" t="0" r="0" b="0"/>
  <wp:docPr id="{int(asset.relationship_id[3:])}" name="{name}" descr="{alt}"/>
  <wp:cNvGraphicFramePr><a:graphicFrameLocks noChangeAspect="1"/></wp:cNvGraphicFramePr>
  <a:graphic><a:graphicData uri="http://schemas.openxmlformats.org/drawingml/2006/picture">
    <pic:pic>
      <pic:nvPicPr><pic:cNvPr id="0" name="{name}"/><pic:cNvPicPr/></pic:nvPicPr>
      <pic:blipFill><a:blip r:embed="{asset.relationship_id}"/><a:stretch><a:fillRect/></a:stretch></pic:blipFill>
      <pic:spPr><a:xfrm><a:off x="0" y="0"/><a:ext cx="{width_emu}" cy="{height_emu}"/></a:xfrm>
        <a:prstGeom prst="rect"><a:avLst/></a:prstGeom>
      </pic:spPr>
    </pic:pic>
  </a:graphicData></a:graphic>
</wp:inline></w:drawing></w:r>"""
    return paragraph(drawing, align="center", before=80, after=160, line=240)


def table_xml(rows: list[list[str]]) -> str:
    if not rows:
        return ""
    columns = max(len(row) for row in rows)
    cell_width = max(700, 9_200 // columns)
    grid = "".join(f'<w:gridCol w:w="{cell_width}"/>' for _ in range(columns))
    table_rows: list[str] = []
    for row_index, values in enumerate(rows):
        cells: list[str] = []
        for column_index in range(columns):
            value = values[column_index] if column_index < len(values) else ""
            fill = "F2F2F2" if row_index == 0 else "FFFFFF"
            cell_content = inline_runs(value)
            if row_index == 0:
                cell_content = run(value, bold=True, size=14)
            elif contains_cjk(value):
                cell_content = inline_runs(value, default_color="3F3F3F")
            cells.append(
                f'<w:tc><w:tcPr><w:tcW w:w="{cell_width}" w:type="dxa"/>'
                f'<w:shd w:val="clear" w:color="auto" w:fill="{fill}"/>'
                '<w:tcMar><w:top w:w="70" w:type="dxa"/><w:left w:w="70" w:type="dxa"/>'
                '<w:bottom w:w="70" w:type="dxa"/><w:right w:w="70" w:type="dxa"/></w:tcMar>'
                f'</w:tcPr>{paragraph(cell_content, after=0, line=210)}</w:tc>'
            )
        table_rows.append(f"<w:tr>{''.join(cells)}</w:tr>")
    return (
        '<w:tbl><w:tblPr><w:tblW w:w="0" w:type="auto"/>'
        '<w:tblLayout w:type="fixed"/><w:tblBorders>'
        '<w:top w:val="single" w:sz="4" w:color="BFBFBF"/>'
        '<w:left w:val="single" w:sz="4" w:color="BFBFBF"/>'
        '<w:bottom w:val="single" w:sz="4" w:color="BFBFBF"/>'
        '<w:right w:val="single" w:sz="4" w:color="BFBFBF"/>'
        '<w:insideH w:val="single" w:sz="4" w:color="D9D9D9"/>'
        '<w:insideV w:val="single" w:sz="4" w:color="D9D9D9"/>'
        f'</w:tblBorders></w:tblPr><w:tblGrid>{grid}</w:tblGrid>{"".join(table_rows)}</w:tbl>'
        + paragraph("", after=100)
    )


def collect_assets(
    lines: list[str],
    *,
    source_path: Path = SOURCE,
    figure_dir: Path | None = None,
) -> tuple[dict[str, ImageAsset], dict[str, bytes]]:
    image_pattern = re.compile(r"^!\[(.*?)\]\((.*?)\)$")
    assets: dict[str, ImageAsset] = {}
    payloads: dict[str, bytes] = {}
    for line in lines:
        match = image_pattern.match(line.strip())
        if not match:
            continue
        relative = match.group(2)
        if relative in assets:
            continue
        source = (
            (figure_dir / Path(relative).name)
            if figure_dir is not None
            else (source_path.parent / relative)
        ).resolve()
        if not source.is_file():
            raise FileNotFoundError(source)
        payload = source.read_bytes()
        width, height = png_dimensions(payload, source)
        index = len(assets) + 3
        archive_name = f"image-{index - 2}-{source.name}"
        assets[relative] = ImageAsset(
            source=source,
            archive_name=archive_name,
            relationship_id=f"rId{index}",
            width_px=width,
            height_px=height,
            alt_text=match.group(1),
        )
        payloads[archive_name] = payload
    return assets, payloads


def markdown_blocks(lines: list[str], assets: dict[str, ImageAsset]) -> list[str]:
    blocks: list[str] = []
    paragraph_lines: list[str] = []
    table_rows: list[list[str]] = []
    first_title = True
    figure_page_break = False

    def flush_paragraph() -> None:
        nonlocal paragraph_lines
        if not paragraph_lines:
            return
        value = " ".join(line.strip().rstrip("  ") for line in paragraph_lines).strip()
        is_chinese = contains_cjk(value)
        blocks.append(
            text_paragraph(
                value,
                chinese=is_chinese,
                after=155 if is_chinese else 55,
                line=270,
                keep_next=value.startswith("**Figure") or value.startswith("**图 "),
                page_break_before=figure_page_break and value.startswith(("**Figure", "**Supplementary")),
            )
        )
        paragraph_lines = []

    def flush_table() -> None:
        nonlocal table_rows
        if table_rows:
            blocks.append(table_xml(table_rows))
            table_rows = []

    image_pattern = re.compile(r"^!\[(.*?)\]\((.*?)\)$")
    for raw_line in lines:
        line = raw_line.rstrip()
        stripped = line.strip()
        if not stripped:
            flush_paragraph()
            flush_table()
            continue
        if stripped.startswith("|") and stripped.endswith("|"):
            flush_paragraph()
            values = [cell.strip() for cell in stripped[1:-1].split("|")]
            if all(re.fullmatch(r":?-{3,}:?", value) for value in values):
                continue
            table_rows.append(values)
            continue
        flush_table()
        image_match = image_pattern.match(stripped)
        if image_match:
            flush_paragraph()
            blocks.append(image_paragraph(assets[image_match.group(2)]))
            continue
        if stripped.startswith("#"):
            flush_paragraph()
            level = len(stripped) - len(stripped.lstrip("#"))
            title = stripped[level:].strip()
            if level == 1:
                blocks.append(
                    paragraph(
                        inline_runs(title),
                        style="Title" if first_title else "Subtitle",
                        align="center",
                        after=80 if first_title else 220,
                        keep_next=True,
                    )
                )
                first_title = False
            elif level == 2:
                page_break = title.startswith(("References", "Figure legends"))
                if title.startswith("Figure legends"):
                    figure_page_break = True
                blocks.append(
                    paragraph(
                        inline_runs(title),
                        style="Heading1",
                        before=180,
                        after=100,
                        keep_next=True,
                        page_break_before=page_break,
                    )
                )
            else:
                blocks.append(
                    paragraph(
                        inline_runs(title),
                        style="Heading2",
                        before=140,
                        after=80,
                        keep_next=True,
                    )
                )
            continue
        if stripped.startswith("- "):
            flush_paragraph()
            value = stripped[2:]
            blocks.append(
                paragraph(
                    inline_runs(value, default_color="3F3F3F" if contains_cjk(value) else None),
                    after=55,
                    line=250,
                    left_indent=360,
                    first_line=-240,
                )
            )
            continue
        if re.match(r"^\d+\. ", stripped):
            flush_paragraph()
            blocks.append(
                paragraph(
                    inline_runs(stripped),
                    after=75,
                    line=250,
                    left_indent=360,
                    first_line=-300,
                )
            )
            continue
        paragraph_lines.append(stripped)
    flush_paragraph()
    flush_table()
    return blocks


def document_xml(lines: list[str], assets: dict[str, ImageAsset]) -> str:
    body = "".join(markdown_blocks(lines, assets))
    section = """
<w:sectPr>
  <w:pgSz w:w="11906" w:h="16838"/>
  <w:pgMar w:top="850" w:right="850" w:bottom="850" w:left="850" w:header="500" w:footer="500" w:gutter="0"/>
  <w:cols w:space="708"/>
  <w:docGrid w:linePitch="360"/>
</w:sectPr>"""
    return f"""<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<w:document xmlns:w="{W_NS}"
 xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships"
 xmlns:wp="http://schemas.openxmlformats.org/drawingml/2006/wordprocessingDrawing"
 xmlns:a="http://schemas.openxmlformats.org/drawingml/2006/main"
 xmlns:pic="http://schemas.openxmlformats.org/drawingml/2006/picture">
<w:body>{body}{section}</w:body></w:document>"""


CONTENT_TYPES = """<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types">
  <Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/>
  <Default Extension="xml" ContentType="application/xml"/>
  <Default Extension="png" ContentType="image/png"/>
  <Override PartName="/word/document.xml" ContentType="application/vnd.openxmlformats-officedocument.wordprocessingml.document.main+xml"/>
  <Override PartName="/word/styles.xml" ContentType="application/vnd.openxmlformats-officedocument.wordprocessingml.styles+xml"/>
  <Override PartName="/word/settings.xml" ContentType="application/vnd.openxmlformats-officedocument.wordprocessingml.settings+xml"/>
  <Override PartName="/docProps/core.xml" ContentType="application/vnd.openxmlformats-package.core-properties+xml"/>
  <Override PartName="/docProps/app.xml" ContentType="application/vnd.openxmlformats-officedocument.extended-properties+xml"/>
</Types>"""

ROOT_RELS = """<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">
  <Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/officeDocument" Target="word/document.xml"/>
  <Relationship Id="rId2" Type="http://schemas.openxmlformats.org/package/2006/relationships/metadata/core-properties" Target="docProps/core.xml"/>
  <Relationship Id="rId3" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/extended-properties" Target="docProps/app.xml"/>
</Relationships>"""

STYLES = """<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<w:styles xmlns:w="http://schemas.openxmlformats.org/wordprocessingml/2006/main">
  <w:docDefaults>
    <w:rPrDefault><w:rPr><w:rFonts w:ascii="Arial" w:hAnsi="Arial" w:eastAsia="Microsoft YaHei"/><w:sz w:val="20"/><w:szCs w:val="20"/><w:lang w:val="en-US" w:eastAsia="zh-CN"/></w:rPr></w:rPrDefault>
    <w:pPrDefault><w:pPr><w:spacing w:line="270" w:lineRule="auto"/></w:pPr></w:pPrDefault>
  </w:docDefaults>
  <w:style w:type="paragraph" w:default="1" w:styleId="Normal"><w:name w:val="Normal"/><w:qFormat/></w:style>
  <w:style w:type="paragraph" w:styleId="Title"><w:name w:val="Title"/><w:qFormat/><w:rPr><w:b/><w:sz w:val="32"/><w:szCs w:val="32"/><w:color w:val="1F1F1F"/></w:rPr></w:style>
  <w:style w:type="paragraph" w:styleId="Subtitle"><w:name w:val="Subtitle"/><w:qFormat/><w:rPr><w:b/><w:sz w:val="28"/><w:szCs w:val="28"/><w:color w:val="3F3F3F"/></w:rPr></w:style>
  <w:style w:type="paragraph" w:styleId="Heading1"><w:name w:val="heading 1"/><w:qFormat/><w:rPr><w:b/><w:sz w:val="27"/><w:szCs w:val="27"/><w:color w:val="1F1F1F"/></w:rPr></w:style>
  <w:style w:type="paragraph" w:styleId="Heading2"><w:name w:val="heading 2"/><w:qFormat/><w:rPr><w:b/><w:sz w:val="23"/><w:szCs w:val="23"/><w:color w:val="2B2B2B"/></w:rPr></w:style>
</w:styles>"""

SETTINGS = """<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<w:settings xmlns:w="http://schemas.openxmlformats.org/wordprocessingml/2006/main">
  <w:zoom w:percent="100"/><w:defaultTabStop w:val="720"/><w:compat/>
</w:settings>"""

CORE = """<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<cp:coreProperties xmlns:cp="http://schemas.openxmlformats.org/package/2006/metadata/core-properties"
 xmlns:dc="http://purl.org/dc/elements/1.1/" xmlns:dcterms="http://purl.org/dc/terms/"
 xmlns:dcmitype="http://purl.org/dc/dcmitype/" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance">
  <dc:title>TxID GigaScience Technical Note — cohort-scale identity revision</dc:title>
  <dc:creator>TxID contributors</dc:creator>
  <dc:subject>Technical Note submission manuscript</dc:subject>
  <dc:description>Bilingual-ready Technical Note manuscript with three code-generated scientific figures and three result tables.</dc:description>
  <cp:revision>1</cp:revision>
  <dcterms:created xsi:type="dcterms:W3CDTF">2026-08-08T00:00:00Z</dcterms:created>
  <dcterms:modified xsi:type="dcterms:W3CDTF">2026-08-08T00:00:00Z</dcterms:modified>
</cp:coreProperties>"""

APP = """<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<Properties xmlns="http://schemas.openxmlformats.org/officeDocument/2006/extended-properties"
 xmlns:vt="http://schemas.openxmlformats.org/officeDocument/2006/docPropsVTypes">
  <Application>TxID publication workflow</Application><AppVersion>1.0</AppVersion>
</Properties>"""


def document_relationships(assets: dict[str, ImageAsset]) -> str:
    relations = [
        '<Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/styles" Target="styles.xml"/>',
        '<Relationship Id="rId2" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/settings" Target="settings.xml"/>',
    ]
    for asset in assets.values():
        relations.append(
            f'<Relationship Id="{asset.relationship_id}" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/image" Target="media/{xml_text(asset.archive_name, attribute=True)}"/>'
        )
    return (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">'
        + "".join(relations)
        + "</Relationships>"
    )


def write_member(archive: zipfile.ZipFile, name: str, payload: str | bytes) -> None:
    info = zipfile.ZipInfo(name, ZIP_TIME)
    info.compress_type = zipfile.ZIP_DEFLATED
    info.external_attr = 0o644 << 16
    archive.writestr(info, payload.encode("utf-8") if isinstance(payload, str) else payload)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, default=SOURCE)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument(
        "--figure-dir",
        type=Path,
        help="Use PNG files with the manuscript figure basenames from this directory.",
    )
    args = parser.parse_args(argv)
    source_path = args.source.resolve()
    output_path = args.output.resolve()
    figure_dir = args.figure_dir.resolve() if args.figure_dir is not None else None
    lines = source_path.read_text(encoding="utf-8").splitlines()
    assets, payloads = collect_assets(
        lines,
        source_path=source_path,
        figure_dir=figure_dir,
    )
    if len(assets) != 3:
        raise ValueError(f"expected 3 manuscript figures, found {len(assets)}")
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(output_path, "w") as archive:
        write_member(archive, "[Content_Types].xml", CONTENT_TYPES)
        write_member(archive, "_rels/.rels", ROOT_RELS)
        write_member(archive, "docProps/core.xml", CORE)
        write_member(archive, "docProps/app.xml", APP)
        write_member(archive, "word/document.xml", document_xml(lines, assets))
        write_member(archive, "word/styles.xml", STYLES)
        write_member(archive, "word/settings.xml", SETTINGS)
        write_member(archive, "word/_rels/document.xml.rels", document_relationships(assets))
        for archive_name, payload in payloads.items():
            write_member(archive, f"word/media/{archive_name}", payload)
    print(output_path)
    return 0


if __name__ == "__main__":
    sys.exit(main())
