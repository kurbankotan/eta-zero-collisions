# Stored outputs

| Folder | Made by | Files |
| --- | --- | --- |
| `colab/` | main notebook, Google Colab CPU session | `results.json`, `versions.json`, `run_log.txt` |
| `dev/` | main notebook, second machine | same three files |
| `archived_replay/` | `python replay_certificates.py --output replay.json` | `replay.json`, `replay_log.txt` |

Package versions are in each `versions.json`.

Compare the replay with the original JSON file:

```sh
python check_certificate.py results/archived_replay/replay.json --compare data/sonuclar.json
```
