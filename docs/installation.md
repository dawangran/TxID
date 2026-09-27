# Installation and containers

TxID **0.1.3 is a research alpha**. The package requires Python **3.10 or later**
and uses only the Python standard library at runtime. The commands below assume
Bash or another POSIX shell.

## Install from source

```bash
git clone https://github.com/dawangran/TxID.git
cd TxID
python -m venv .venv
source .venv/bin/activate
python -m pip install .
txid --version
```

The version command prints `txid 0.1.3`. On Windows, activate the environment
with `.venv\Scripts\Activate.ps1` in PowerShell. A local
[Conda recipe](../conda-recipe/meta.yaml) is included for building the package.

For development, use `python -m pip install -e .` in the virtual environment.
Run the source tests from the repository root:

```bash
python -m unittest discover -s tests -v
```

Continue with the [tutorial](tutorial.md) for a registry built from your reference
FASTA and caller annotations. The [minimal example](../examples/minimal/README.md)
uses a small synthetic fixture and requires no genome download.

## Build a local container

The [Dockerfile](../Dockerfile) builds TxID with JupyterLab on a base image pinned
by digest. Run these commands from the repository root:

```bash
docker build -t txid:0.1.3-jupyter .
docker run --rm txid:0.1.3-jupyter txid --version
docker run --rm -p 8888:8888 -v "$PWD:/workspace" txid:0.1.3-jupyter
```

The final command starts JupyterLab and prints its access token. The container
runs as user `txid` with UID 1000; mounted directories must be writable by that
user. Record the final image digest for an archival run. The base-image digest
alone does not identify the complete built image.

## Workflow integration

| Workflow | Use |
| --- | --- |
| [`txid_multi_add.wdl`](../workflows/txid_multi_add.wdl) | Same-tool imports; [example inputs](../workflows/txid_multi_add.inputs.example.json) and [中文 guide](txid-multi-add-wdl-guide.zh-CN.md). |
| [`txid_batch.wdl`](../workflows/txid_batch.wdl) | Cohort imports with per-input provenance; [example inputs](../workflows/txid_batch.inputs.example.json) and [中文 guide](txid-batch-wdl-guide.zh-CN.md). |
| [Benchmark workflow](../workflows/publication_benchmark/README.md) | Snakemake evaluation with separate caller and comparator environments. |

WDL imports execute within one task because isolated scatter tasks cannot safely
mutate a shared SQLite registry. Set the workflow's `docker_image` input to a
published final image digest for archival execution.

The real-read benchmark runner, `benchmarks/run_real_interop_callers.py`, requires
Python **3.11 or later** and separately installed alignment and caller tools.
These benchmark prerequisites do not change the TxID package's Python 3.10
minimum.

Record the source commit, software version, reference FASTA, annotation release
and command options for reproducible analyses. The registry records input
checksums and import provenance; see the [reproducibility checklist](reproducibility.md).
