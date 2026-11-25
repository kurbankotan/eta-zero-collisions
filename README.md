# eta-zero-collisions

Companion repository for a manuscript by Kurban Kotan (under review) on zero collisions in a family of Fourier transforms built from the Dedekind eta function. This page only explains how to rerun the computations. The mathematics, the statements that are proved and their limits are in the manuscript.

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/kurbankotan/eta-zero-collisions/blob/main/eta_edge_zeros_reproduce.ipynb)

## Run it

**Google Colab.** Click the badge, then choose *Runtime → Run all*. A CPU session is enough and the run takes one to two minutes. When it ends, the browser downloads `eta_edge_zeros_results.zip`.

**Terminal** (Python 3.11 or newer):

```sh
python -m venv .venv
. .venv/bin/activate              # Windows PowerShell: .\.venv\Scripts\Activate.ps1
python -m pip install -r requirements.txt jupyter
jupyter nbconvert --to notebook --execute eta_edge_zeros_reproduce.ipynb --output executed.ipynb
```

**Older mpmath script and its checker:**

```sh
python check_certificate.py data/sonuclar.json
python replay_certificates.py --output replay.json
python check_certificate.py replay.json --compare data/sonuclar.json
```

Adding `--include-prime-scan` to the second command also runs an optional prime experiment.

## What is where

| Path | What it is |
| --- | --- |
| `eta_edge_zeros_reproduce.ipynb` | Main notebook (python-flint, Arb). Writes `results.json`, `versions.json`, `run_log.txt`, two PNG figures and the zip file. |
| `results/colab/`, `results/dev/` | Output of the main notebook from two machines |
| `results/archived_replay/` | Output of `replay_certificates.py` |
| `figures/` | PNG figures |
| `replay_certificates.py` | Older mpmath script |
| `check_certificate.py` | Checker for the JSON files of that script (standard library only) |
| `data/sonuclar.json` | JSON file from the original run of the mpmath script |
| `Kubik_L_arastirmasi.ipynb` | Original notebook (Turkish) |
| `requirements.txt`, `CITATION.cff`, `LICENSE` | Dependencies, citation metadata, MIT License |

Package versions of every stored run are listed in its `versions.json`.

## Notes

- Large language models were used in this project. The manuscript describes how.
- Title, journal and DOI will be added here once the article is published.
