# Result files

| Folder | Produced by | Files |
| --- | --- | --- |
| `colab/` | `eta_edge_zeros_reproduce.ipynb` on a Google Colab CPU runtime (Python 3.13.18, python-flint 0.9.0, NumPy 2.4.6, Matplotlib 3.10.0) | `results.json`, `versions.json`, `run_log.txt` |
| `dev/` | The same notebook in the development environment (Python 3.13.18, python-flint 0.9.0, NumPy 2.5.5, Matplotlib 3.11.2) | `results.json`, `versions.json`, `run_log.txt` |
| `archived_replay/` | `python replay_certificates.py --output replay.json` (Python 3.13.18, mpmath 1.3.0) | `replay.json`, `replay_log.txt` |

`results.json` contains every number reported in the Results section of the manuscript. Interval endpoints are rounded outward. The two notebook runs agree in all 262 recorded fields except `lattice_check_s2[3].rel_diff` (4.3e-15 versus 4.0e-15) and `scan.float64_max_abs_error_vs_arb` (1.0e-17 versus 1.6e-17), both float64 rounding-level quantities.

`archived_replay/replay.json` can be compared with the archived file `data/sonuclar.json`:

```sh
python check_certificate.py results/archived_replay/replay.json --compare data/sonuclar.json
```

Of the 53 fields common to the two files, 52 are identical; the remaining one is the recorded Python version.
