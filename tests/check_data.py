#!/usr/bin/env python3
"""核对随库科学数据的大小与 SHA-256。"""
from pathlib import Path
import hashlib
import json
root = Path(__file__).resolve().parents[1] / 'data'
for row in json.loads((root / 'manifest.json').read_text()):
    path = root / row['file']
    assert path.is_file() and not path.is_symlink(), path
    assert path.stat().st_size == row['bytes'], path
    assert hashlib.sha256(path.read_bytes()).hexdigest() == row['sha256'], path
print('PASS: 4 scientific data files')
