# A Certified Double Zero in a Dedekind Eta Scale Family

Certified nonreal zeros and a double-zero collision in a Dedekind eta Fourier family, with reproducible interval-arithmetic code and data.

This repository reproduces a computer-assisted result for the entire functions

$$
H_a(z)=\int_0^\infty e^{u/2}\eta(ie^{u+a})\eta(ie^{u-a})\cos(zu)\,du.
$$

The double-zero system has a unique solution within a box of coordinate radius $10^{-19}$ centered at

$$
(a_c,x_c)=(1.45200626937306544414442672465215400624761,\;
14.8380401130710086477956510630327906219061).
$$

The verified derivative signs $\partial_a H<0<\partial_x^2 H$ imply a local transition. Just below the enclosed parameter there is a conjugate pair of nonreal zeros near the collision. Just above it there are two real zeros near the collision.

At $a=23/16$, a separate computation proves a negative Laguerre expression at $x=74/5$ and certifies one simple zero in each of two conjugate disks of radius $10^{-12}$. The upper disk is centered at $14.67872714035986218218+0.31107587891741700810i$.

The modular Artin specialization occurs at $a=\tfrac12\log 3$. Its zero distribution remains open in this study. The notebook and replay also provide an optional prime count for $t^3\equiv 2\bmod p$ through $10^7$ as a finite computational experiment.

## Reproduce from the terminal

Use Python 3.11 or newer. The bundled report and terminal replay were checked with Python 3.14.6 and mpmath 1.3.0.

```sh
python -m venv .venv
```

Activate the environment with `. .venv/bin/activate` on macOS or Linux, or `.\.venv\Scripts\Activate.ps1` in Windows PowerShell. Then run the following commands from the repository root.

```sh
python -m pip install -r requirements.txt
python check_certificate.py data/sonuclar.json
python replay_certificates.py --output replay.json
python check_certificate.py replay.json --compare data/sonuclar.json
```

The replay recomputes all three interval certificates on the CPU. Add `--include-prime-scan` to also regenerate the prime experiment through $10^7$. The checker uses exact Python fractions to check the exported bounds and certificate inequalities. Its optional comparison checks matching certificate inputs and compatible interval enclosures. Use the replay to recompute the enclosures from the eta product.

## Reproduce in Google Colab

1. Open [Google Colab](https://colab.research.google.com/) and upload [Kubik_L_arastirmasi.ipynb](Kubik_L_arastirmasi.ipynb).
2. Choose a CPU runtime. If necessary, run `%pip install mpmath==1.3.0 plotly==6.9.0` in a setup cell.
3. Run all cells in order. The last cell downloads `kubik_L_sonuclari.zip`, containing the generated `sonuclar.json`.

## Files

| File | Contents |
| --- | --- |
| [replay_certificates.py](replay_certificates.py) | Terminal implementation of the notebook calculations and JSON exporter |
| [check_certificate.py](check_certificate.py) | Standard-library checks for stored inequalities and replay comparisons |
| [Kubik_L_arastirmasi.ipynb](Kubik_L_arastirmasi.ipynb) | Original Turkish notebook with the calculations, analytic bounds, prime scan, and export cell |
| [data/sonuclar.json](data/sonuclar.json) | Unmodified JSON from the available local notebook replay, including exact interval endpoints and the prime samples |
| [requirements.txt](requirements.txt) | Pinned interval arithmetic and notebook plotting dependencies |

## Certification method

Each enclosure uses the Dedekind eta inversion formula and product, exact rational quadrature weights, interval evaluation at quadrature nodes, and added analytic bounds for truncation and integration errors. A fixed-point contraction certifies the unique solution of $H=\partial_xH=0$. The Laguerre inequality and Rouché’s theorem provide two nonreal-zero certificates at the fixed scale. The notebook records the analytic estimates and the local splitting argument. [NIST DLMF, eta inversion](https://dlmf.nist.gov/23.18.E5) · [Csordas and Escassut, 2005](https://doi.org/10.5802/ambp.210)
