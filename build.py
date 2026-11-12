#!/usr/bin/env python3
"""Build index.html by inlining src/physics.js into src/template.html (the page is one self-contained file)."""
from pathlib import Path

root = Path(__file__).parent
# Read and write as UTF-8 explicitly: the sources contain characters such as θ and ₁, and without this Python on
# Windows (before 3.15) falls back to the cp1252 code page and fails.
template = (root / "src" / "template.html").read_text(encoding="utf-8")
physics = (root / "src" / "physics.js").read_text(encoding="utf-8")
assert "/*PHYSICS*/" in template, "template is missing the /*PHYSICS*/ placeholder"
(root / "index.html").write_text(template.replace("/*PHYSICS*/", physics), encoding="utf-8", newline="\n")
print("Wrote index.html")
