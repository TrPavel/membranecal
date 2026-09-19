"""Prepare a separate source-reconstruction workspace from RC1; no downloads or scores."""
import argparse
import hashlib
import json
import shutil
from pathlib import Path


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--release-root', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    args = p.parse_args()
    root, out = args.release_root.resolve(), args.output.resolve()
    if out.exists() or root == out or root in out.parents:
        raise ValueError('Choose a new output directory outside the release.')
    manifest = json.loads((root / 'manifest.json').read_text())['files']
    for name, digest in manifest.items():
        path = (root / name).resolve()
        if not path.is_relative_to(root) or hashlib.sha256(path.read_bytes()).hexdigest() != digest:
            raise ValueError(name)
    shutil.copytree(Path(__file__).parent / 'scripts', out / 'scripts')
    work = out / 'workspace/v02'
    work.mkdir(parents=True)
    development = json.loads((root / 'data/development_sequences.json').read_text())
    old = out / 'data/frozen/membranecal_real_pilot_v0.1'
    old.mkdir(parents=True)
    # Full original development record fields already present in RC1; includes the
    # accession/entry/construct exclusion keys used by the original acquisition code.
    (old / 'proteins.json').write_text(json.dumps({'records': development['records']}, indent=2))
    shutil.copyfile(root / 'data/development_sequences.json', work / 'development_sequences.json')
    shutil.copyfile(root / 'data/proteins.json', out / 'validation_records.json')
    shutil.copyfile(root / 'data/sources.json', out / 'expected_source_hashes.json')
    print('Scaffold created. Source downloads require explicit subsequent acquisition commands; no new study freeze was created.')


if __name__ == '__main__':
    main()
