# Certified Edge Zeros in a Dedekind Eta Scale Family

Interval certificates, an independent ball-arithmetic replay, and an AI-assisted verification workflow.

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/kurbankotan/eta-zero-collisions/blob/main/eta_edge_zeros_reproduce.ipynb)

This repository contains the code, notebooks and result files for the manuscript *Certified Edge Zeros in a Dedekind Eta Scale Family* by Kurban Kotan (submitted). The manuscript studies the one-parameter family

$$
H(a,z)=\int_0^\infty e^{u/2}\eta\left(ie^{u+a}\right)\eta\left(ie^{u-a}\right)\cos(zu)\mathrm{d}u ,
$$

where $\eta$ is the Dedekind eta function. For real $a$, the slice $H(a,\cdot)$ is an even entire function that is real on the real axis. At the seven Artin points $a=\frac{1}{2}\log(B/A)$ with $A\in\lbrace 1,2,3,4,6,8,12\rbrace$ and $A+B=24$, it is the completed $L$-function of the weight-one eta product $\eta(A\tau)\eta(B\tau)$. At $a_0=\frac{1}{2}\log 3$ this is the newform 108.1.c.a, whose $L$-function is the Artin $L$-function attached to $T^3-2$.

## Certified results

Every statement in the table is proved by interval or ball arithmetic combined with explicit analytic error bounds.

| Certificate | Statement |
| --- | --- |
| Right edge $a_R$ | A box test of radius $10^{-19}$ about $(a_c,x_c)=(1.45200626937306544414\ldots, 14.83804011307100864779\ldots)$ shows that $H=H_x=0$ has exactly one solution $(a_R,x_R)$ in the box, so $x_R$ is a double zero of $H(a_R,\cdot)$. On the box $H_a<0<H_{xx}$, so near this point a conjugate pair of nonreal zeros returns to the real axis as $a$ increases through $a_R$. |
| Left edge $a_L$ | The same test with radius $10^{-30}$ about $(1.34898938954386827787\ldots, 13.82918665572891503246\ldots)$ certifies a second double zero. On the box $H_a>0$ and $H_{xx}>0$, so near this point two real zeros meet and leave the real axis as $a$ increases through $a_L$. |
| Nonreal zeros at $a=23/16$ | The Laguerre expression $H_x^2-HH_{xx}$ at $x=74/5$ lies in $[-1.74472143,-1.74472142]\times10^{-20}$, and Rouché disks of radius $10^{-12}$ and $10^{-19}$ about $14.67872714035986218218+0.31107587891741700810i$ each contain exactly one zero. |
| Nonreal zeros at seven parameter values | A negative Laguerre value is certified at one parameter value inside the tracked interval of each of the seven excursions (A to G) found by the exploratory scan of $0\le a\le1.75$, $0<x<20$. The test proves that nonreal zeros exist at these parameter values; it does not locate them. |
| LMFDB zeros at $a_0$ | Arb sign changes show that $H(a_0,\cdot)$ has a real zero within $10^{-25}$ of each of the ten smallest zero ordinates of 108.1.c.a listed by the LMFDB. |

**Finite checks.** The multiplicativity test of the twelve products $\eta(A\tau)\eta(B\tau)$ up to $n\le3000$ (exactly the seven products at the Artin points pass), Gauss's cubic-residue criterion for $\eta(6\tau)\eta(18\tau)$ at 428 primes, and the lattice-sum identity at $s=2$ for the seven Artin points (float64 lattice sum, relative difference below $10^{-14}$).

**Not claimed.** The certificates are local. They establish two opposite transitions, but they do not provide a certified continuation joining $a_L$ and $a_R$, so they do not prove that the pair leaving the real axis at $a_L$ is the pair returning at $a_R$; that identification rests on float64 tracking. The computations also do not test the generalized Riemann hypothesis for the seven $L$-functions and do not give a global zero count. The zero map ([figures/fig_zero_landscape.png](figures/fig_zero_landscape.png)) is an exploratory float64 computation, not a certificate.

## Repository contents

| Path | Contents |
| --- | --- |
| [eta_edge_zeros_reproduce.ipynb](eta_edge_zeros_reproduce.ipynb) | Reproduction notebook (Arb through python-flint). It computes every number reported in the Results section of the manuscript, writes `results.json`, `versions.json`, `run_log.txt` and the figures, and packs them into `eta_edge_zeros_results.zip`, which is downloaded automatically in Google Colab. |
| [results/colab/](results/colab) | Output of the reproduction notebook on a Google Colab CPU runtime |
| [results/dev/](results/dev) | Output of the same notebook in the development environment |
| [results/archived_replay/](results/archived_replay) | Replay of the archived mpmath certificates with `replay_certificates.py` |
| [figures/](figures) | Figures of the manuscript; `workflow.dot` is the Graphviz source of the workflow diagram |
| [replay_certificates.py](replay_certificates.py) | Archived mpmath interval script (Laguerre, Rouché and right-edge certificates) |
| [check_certificate.py](check_certificate.py) | Checker that re-verifies the stored inequalities in exact rational arithmetic |
| [data/sonuclar.json](data/sonuclar.json) | Archived result file with all interval endpoints (Python 3.14.6, mpmath 1.3.0) |
| [Kubik_L_arastirmasi.ipynb](Kubik_L_arastirmasi.ipynb) | Original Colab notebook (in Turkish); it also contains a supplementary prime scan that is not used in the manuscript |
| [requirements.txt](requirements.txt) | Dependencies |
| [CITATION.cff](CITATION.cff) | Citation metadata |
| [LICENSE](LICENSE) | MIT License |

## Reproduce in Google Colab

1. Click the **Open in Colab** badge above.
2. Keep the default CPU runtime and choose **Runtime → Run all**.
3. The notebook installs python-flint with pip and finishes in one to two minutes. Its last cell downloads `eta_edge_zeros_results.zip` with `results.json`, `versions.json`, `run_log.txt` and the two figures.

The Colab run and the development run agree in every certificate and finite check. Of the 262 recorded result fields, only two differ, and both are float64 rounding-level quantities: one relative difference in the lattice-sum check (`lattice_check_s2`) and the accuracy bound of the float64 evaluator (`float64_max_abs_error_vs_arb`).

| Run | Python | python-flint | NumPy | Matplotlib |
| --- | --- | --- | --- | --- |
| Google Colab, CPU | 3.13.18 | 0.9.0 | 2.4.6 | 3.10.0 |
| Development environment | 3.13.18 | 0.9.0 | 2.5.5 | 3.11.2 |

The figures in [figures/](figures) were produced by the development run. The Colab run produced the same underlying data.

## Reproduce from a terminal

Use Python 3.11 or newer.

```sh
python -m venv .venv
```

Activate the environment with `. .venv/bin/activate` on macOS or Linux, or `.\.venv\Scripts\Activate.ps1` in Windows PowerShell. Then run the following commands from the repository root.

```sh
python -m pip install -r requirements.txt jupyter

# reproduction notebook (Arb)
jupyter nbconvert --to notebook --execute eta_edge_zeros_reproduce.ipynb --output executed.ipynb

# archived mpmath certificates
python check_certificate.py data/sonuclar.json
python replay_certificates.py --output replay.json
python check_certificate.py replay.json --compare data/sonuclar.json
```

The replay recomputes the three archived interval certificates on the CPU in about half a minute. Add `--include-prime-scan` to also regenerate the supplementary prime experiment. The checker uses exact Python fractions; with `--compare` it also checks that the replayed enclosures are compatible with the archived ones.

## Method in brief

Each integral is evaluated over $[0,6]$; the tail beyond $u=6$ is below $10^{-80}$. The Arb computations use rigorous numerical integration at 320 bits, and $\eta$ is evaluated with Arb's modular-form routine. A negative Laguerre expression $H_x^2-HH_{xx}$ at a real point certifies nonreal zeros of the slice. A Rouché test isolates a single nonreal zero. A Banach/Krawczyk-type box test, with an analytic bound on the variation of the Jacobian, certifies a unique solution of $H=H_x=0$, that is, a double zero. The signs of $H_a$ and $H_{xx}$ at that point decide whether the two real zeros meet and leave the real axis or return to it as $a$ increases. The archived certificates use mpmath interval arithmetic with Lagrange quadrature weights computed in exact rational arithmetic and Cauchy remainder bounds.

## AI-assisted workflow

The archived certificates, the mpmath code and the first draft of the manuscript were developed with ChatGPT and Codex (OpenAI). Claude (Anthropic) was then used to review that work, replay the certificates in Arb, certify the left edge and the excursion points, compute the exploratory zero map, verify the references and expand the manuscript. A numerical statement was accepted only if its inequality held in interval or ball arithmetic. The conversations with the language models are not part of this repository.

## Citation

If you use this code or data, please cite the manuscript. The metadata in [CITATION.cff](CITATION.cff) will be updated when the article is published.

## License

MIT License; see [LICENSE](LICENSE).
