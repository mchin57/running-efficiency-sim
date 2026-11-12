#!/usr/bin/env python3
"""Build index.html by inlining src/physics.js into src/template.html (the page is one self-contained file)."""
from pathlib import Path

root = Path(__file__).parent
template = (root / "src" / "template.html").read_text()
physics = (root / "src" / "physics.js").read_text()
assert "/*PHYSICS*/" in template, "template is missing the /*PHYSICS*/ placeholder"
(root / "index.html").write_text(template.replace("/*PHYSICS*/", physics))
print("Wrote index.html")
