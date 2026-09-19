# Reproduction from a clean checkout or archive

Use Python 3.12 in a new virtual environment. On Windows PowerShell activate it
with `.venv/Scripts/Activate.ps1`; on POSIX use `source .venv/bin/activate`.
Unset any existing `PYTHONPATH` before replay (`Remove-Item Env:PYTHONPATH` in
PowerShell if set; `unset PYTHONPATH` in POSIX). Assertions must remain enabled.

```bash
python -m pip install --require-hashes -r publication/v0.2/code/numerical/environment/requirements.lock
python tools/verify.py
python tools/replay.py --mode compact --output ../compact-output
python tools/replay.py --mode full --output ../full-output
```

Outputs must be new or empty and outside the checkout. Dependency installation
uses the network; replay does not. Compact mode verifies all 554 compact-coordinate
score records against the frozen scores, without creating new scientific estimates.
Full mode runs the existing standalone pipeline: all 36 scenarios, independent
primary bootstrap verification (20,000 draws each), nine tables and five figures.
It verifies outputs against existing frozen results, not a new analysis protocol.

The full figure acceptance environment is Windows, Python 3.12 and the hash-locked
dependencies, with Arial regular/bold fonts described in
[tested_environment.json](../publication/v0.2/code/numerical/environment/tested_environment.json).
Fonts are not redistributed. Full mode fails if figure pixels differ; numerical-only
success is never presented as full figure acceptance. Linux compact replay is
checked separately. A fresh renderer cache is created inside the output directory.

For a stronger Python-level access check, run from a fresh environment:

```bash
python -I tools/isolate.py . ../isolated-output --mode full
```

The guard rejects Python file opens outside the checkout, output, interpreter,
environment and system fonts, and rejects socket connections and subprocess
fallbacks. This is not an OS security sandbox; hosted CI additionally starts on a
fresh runner with no research repository checkout or private cache.

To build and test the exact release ZIP:

```bash
python tools/package.py --output ../membranecal-assets
```

Extract the ZIP into a new directory, install its included dependency lock into
a fresh environment, and run the same verify/replay commands there. The package
builder reports archive and manifest hashes. Git is not needed for extracted replay.
Historical source acquisition is described separately in the payload's upstream
README and is not required for the offline replay above.
