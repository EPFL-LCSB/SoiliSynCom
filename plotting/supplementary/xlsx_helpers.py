"""Minimal Strict OOXML worksheet reader used by the supplementary plotter."""
import re
import xml.etree.ElementTree as ET
from pathlib import Path, PurePosixPath
from zipfile import ZipFile


def _clean(value):
    return "" if value is None else str(value).strip()


def strict_xlsx_sheets(path: Path):
    with ZipFile(path) as archive:
        shared = []
        if "xl/sharedStrings.xml" in archive.namelist():
            root = ET.fromstring(archive.read("xl/sharedStrings.xml"))
            for item in root.findall(".//{*}si"):
                shared.append("".join(node.text or "" for node in item.findall(".//{*}t")))
        workbook = ET.fromstring(archive.read("xl/workbook.xml"))
        rels = ET.fromstring(archive.read("xl/_rels/workbook.xml.rels"))
        targets = {
            rel.attrib["Id"]: rel.attrib["Target"]
            for rel in rels.findall(".//{*}Relationship")
        }
        result = {}
        for sheet in workbook.findall(".//{*}sheet"):
            name = sheet.attrib["name"]
            rel_id = next(
                value for key, value in sheet.attrib.items() if key.endswith("}id")
            )
            target = targets[rel_id].lstrip("/")
            if not target.startswith("xl/"):
                target = str(PurePosixPath("xl") / target)
            xml = ET.fromstring(archive.read(target))
            parsed_rows = []
            for row in xml.findall(".//{*}sheetData/{*}row"):
                cells = {}
                for cell in row.findall("{*}c"):
                    ref = cell.attrib.get("r", "A1")
                    match = re.match(r"[A-Z]+", ref)
                    if match is None:
                        continue
                    column = 0
                    for letter in match.group(0):
                        column = column * 26 + ord(letter) - 64
                    value_node = cell.find("{*}v")
                    value = "" if value_node is None else (value_node.text or "")
                    if cell.attrib.get("t") == "s" and value:
                        value = shared[int(value)]
                    elif cell.attrib.get("t") == "inlineStr":
                        value = "".join(
                            node.text or "" for node in cell.findall(".//{*}t")
                        )
                    cells[column - 1] = _clean(value)
                if cells:
                    parsed_rows.append(
                        [cells.get(i, "") for i in range(max(cells) + 1)]
                    )
            result[name] = parsed_rows
        return result
