"""Recompute the eta interval certificates on the CPU and export their bounds."""

import argparse
import json
from fractions import Fraction
from math import isqrt, log, log10, sqrt
from pathlib import Path
from platform import python_version
from time import perf_counter

import mpmath
from mpmath import iv


scale = Fraction(23, 16)
frequency = Fraction(74, 5)
resolution = 32
degree = 24
end = 6


def box(rational):
    rational = Fraction(rational)
    return iv.mpf(rational.numerator) / rational.denominator


def weights_for_interpolation(degree):
    weights = []
    for j in range(degree + 1):
        coefficients = [1]
        denominator = 1
        for k in range(degree + 1):
            if j == k:
                continue
            expanded = [0] * (len(coefficients) + 1)
            for power, coefficient in enumerate(coefficients):
                expanded[power] -= k * coefficient
                expanded[power + 1] += coefficient
            coefficients = expanded
            denominator *= j - k
        weights.append(sum(Fraction(c * degree**p, (p + 1) * denominator)
                           for p, c in enumerate(coefficients)))
    assert sum(weights) == 1
    return weights


def eta_on_axis(log_y):
    y = iv.exp(box(abs(log_y)))
    q = iv.exp(-2 * iv.pi * y)
    product = iv.mpf(1)
    for n in range(1, 15):
        product *= 1 - q**n
    omitted = q**15 / (1 - q)
    product *= 1 - iv.mpf([0, 1]) * omitted
    return product * iv.exp(-iv.pi * y / 12 + box(max(-log_y, 0) / 2))


def certify_laguerre():
    iv.dps = 45
    weights = weights_for_interpolation(degree)
    segments = end * resolution
    node_count = segments * degree + 1
    node_weights = [Fraction(0)] * node_count
    for start in range(0, node_count - 1, degree):
        for j, weight in enumerate(weights):
            node_weights[start + j] += weight / resolution

    values = [iv.mpf(0), iv.mpf(0), iv.mpf(0)]
    for i, weight in enumerate(node_weights):
        u = Fraction(i, resolution * degree)
        u_iv = box(u)
        w = box(weight) * iv.exp(u_iv / 2) * eta_on_axis(u + scale) * eta_on_axis(u - scale)
        angle = box(frequency * u)
        c = iv.cos(angle)
        values[0] += w * c
        values[1] -= w * u_iv * iv.sin(angle)
        values[2] -= w * u_iv**2 * c

    # Analytic error bounds for each integral, independent of the interval arithmetic.
    radius = Fraction(1, 4)
    sup_bound = 10**6
    error = end * sup_bound * (1 / (resolution * radius))**(degree + 1) + Fraction(1, 10**80)
    values = [value + iv.mpf([-1, 1]) * box(error) for value in values]
    h, dh, ddh = values
    laguerre = dh**2 - h * ddh
    assert laguerre.b < 0
    return values, laguerre, error


def certify_nonreal_zero():
    """İlk doğrulama hücresindeki eta ve rasyonel ağırlık yardımcılarını kullanır."""
    iv.dps = 60
    scale = Fraction(23, 16)
    real_text = '14.67872714035986218218'
    imag_text = '0.31107587891741700810'
    real_center, imag_center = Fraction(real_text), Fraction(imag_text)
    radius = Fraction(1, 10**12)
    resolution, degree, end = 32, 40, 6
    weights = weights_for_interpolation(degree)
    count = end * resolution * degree + 1
    node_weights = [Fraction(0)] * count
    for start in range(0, count - 1, degree):
        for i, weight in enumerate(weights):
            node_weights[start + i] += weight / resolution

    real_value, imag_value, real_slope = [iv.mpf(0) for _ in range(3)]
    for i, weight in enumerate(node_weights):
        u = Fraction(i, resolution * degree)
        u_iv = box(u)
        w = box(weight) * iv.exp(u_iv / 2) * eta_on_axis(u + scale) * eta_on_axis(u - scale)
        angle = box(real_center * u)
        exp_beta = iv.exp(box(imag_center * u))
        cosh = (exp_beta + 1 / exp_beta) / 2
        sinh = (exp_beta - 1 / exp_beta) / 2
        real_value += w * iv.cos(angle) * cosh
        imag_value -= w * iv.sin(angle) * sinh
        real_slope -= w * u_iv * iv.sin(angle) * cosh

    # Complex-u Cauchy bound for these integrands and explicit u>6 tail.
    error = end * 10**8 * Fraction(1, 8)**(degree + 1) + Fraction(1, 10**80)
    uncertainty = iv.mpf([-1, 1]) * box(error)
    real_value += uncertainty
    imag_value += uncertainty
    real_slope += uncertainty
    assert real_slope.a > 0
    # On the disk |H''| < 20. Rouché compares H with H'(center)(z-center).
    remainder = abs(real_value) + abs(imag_value) + box(10 * radius**2)
    linear = real_slope * box(radius)
    assert remainder.b < linear.a
    assert imag_center > radius
    return {
        "scale": scale,
        "positive_imaginary_center": {"real": real_text, "imaginary": imag_text},
        "radius": radius,
        "zeros_per_disk_counting_multiplicity": 1,
        "H_at_center": {"real": real_value, "imaginary": imag_value},
        "real_H_prime_at_center": real_slope,
        "rouche_remainder_upper_bound": remainder.b,
        "rouche_linear_lower_bound": linear.a,
        "analytic_error_bound_per_integral": error,
        "interval_precision_digits": iv.dps,
    }


def certify_collision(scale_text, frequency_text, inverse_text):
    """Verilen rasyonel merkezin yakınında H=H_x=0 çözümünü doğrular."""
    iv.dps = 70
    scale = Fraction(scale_text)
    frequency = Fraction(frequency_text)
    resolution = 32
    degree = 48
    end = 6
    terms = 24
    radius = Fraction(1, 10**19)
    assert Fraction('1.4') <= scale - radius <= scale + radius <= Fraction('1.5')
    assert Fraction('14.5') <= frequency - radius <= frequency + radius <= Fraction('15.5')

    def box(number):
        number = Fraction(number)
        return iv.mpf(number.numerator) / number.denominator

    def norm(vector):
        return max(abs(item).b for item in vector)

    def matrix_norm(matrix):
        return max(sum(abs(item) for item in row).b for row in matrix)

    def interpolation_weights():
        weights = []
        for j in range(degree + 1):
            coefficients = [1]
            divisor = 1
            for k in range(degree + 1):
                if j != k:
                    expanded = [0] * (len(coefficients) + 1)
                    for power, coefficient in enumerate(coefficients):
                        expanded[power] -= k * coefficient
                        expanded[power + 1] += coefficient
                    coefficients = expanded
                    divisor *= j - k
            weights.append(sum(Fraction(c * degree**p, (p + 1) * divisor)
                               for p, c in enumerate(coefficients)))
        # Check polynomial exactness over [0,1] in exact rational arithmetic.
        for power in range(degree + 1):
            assert sum(w * Fraction(j, degree)**power for j, w in enumerate(weights)) == Fraction(1, power + 1)
        return weights

    unit_interval = iv.mpf([0, 1])
    symmetric_interval = iv.mpf([-1, 1])

    def eta_and_log_derivative(log_y):
        y = iv.exp(box(abs(log_y)))
        q = iv.exp(-2 * iv.pi * y)
        product = iv.mpf(1)
        series = iv.mpf(0)
        power = q
        for n in range(1, terms + 1):
            factor = 1 - power
            product *= factor
            series += n * power / factor
            power *= q
        product *= 1 - unit_interval * power / (1 - q)
        series += unit_interval * power * (terms + 1 - terms * q) / (1 - q)**3
        value = product * iv.exp(-iv.pi * y / 12)
        slope = iv.pi * y * (2 * series - iv.mpf(1) / 12)
        if log_y < 0:
            value *= iv.exp(box(-log_y / 2))
            slope = -iv.mpf(1) / 2 - slope
        return value, slope

    weights = interpolation_weights()
    node_count = end * resolution * degree + 1
    node_weights = [Fraction(0)] * node_count
    for start in range(0, node_count - 1, degree):
        for index, weight in enumerate(weights):
            node_weights[start + index] += weight / resolution

    values = [iv.mpf(0) for _ in range(5)]
    for i, weight in enumerate(node_weights):
        u = Fraction(i, resolution * degree)
        u_iv = box(u)
        plus, plus_slope = eta_and_log_derivative(u + scale)
        minus, minus_slope = eta_and_log_derivative(u - scale)
        w = box(weight) * iv.exp(u_iv / 2) * plus * minus
        angle = box(frequency * u)
        cosine = iv.cos(angle)
        sine = iv.sin(angle)
        wa = w * (plus_slope - minus_slope)
        values[0] += w * cosine            # H
        values[1] -= w * u_iv * sine       # H_x
        values[2] += wa * cosine           # H_a
        values[3] -= wa * u_iv * sine      # H_ax
        values[4] -= w * u_iv**2 * cosine  # H_xx

    analytic_error = end * 10**7 * Fraction(1, 8)**(degree + 1) + Fraction(1, 10**80)
    h, hx, ha, hax, hxx = [v + symmetric_interval * box(analytic_error) for v in values]
    jacobian = [(ha, hx), (hax, hxx)]
    rational_inverse = [[Fraction(entry) for entry in row] for row in inverse_text]
    assert (rational_inverse[0][0] * rational_inverse[1][1]
            - rational_inverse[0][1] * rational_inverse[1][0]) != 0
    inverse = [[box(entry) for entry in row] for row in rational_inverse]

    residual = norm([sum(coefficient * v for coefficient, v in zip(row, (h, hx)))
                     for row in inverse])
    inverse_error = matrix_norm([
        [int(i == j) - sum(inverse[i][k] * jacobian[k][j] for k in range(2))
         for j in range(2)]
        for i in range(2)
    ])
    derivative_change = box(50_000 * radius)
    contraction = inverse_error + matrix_norm(inverse) * derivative_change

    # Banach's theorem on the closed box, with all comparisons in intervals.
    assert contraction.b < 1
    assert (residual + contraction * box(radius)).b < box(radius).a
    ha_box = ha + symmetric_interval * derivative_change
    hxx_box = hxx + symmetric_interval * derivative_change
    assert ha_box.b < 0 and hxx_box.a > 0
    scale_interval = box(scale) + symmetric_interval * box(radius)
    frequency_interval = box(frequency) + symmetric_interval * box(radius)
    return {
        "center": {"scale": scale_text, "frequency": frequency_text},
        "box_radius": radius,
        "scale_interval": scale_interval,
        "frequency_interval": frequency_interval,
        "inverse_jacobian_approximation": rational_inverse,
        "center_displacement_upper_bound": residual,
        "contraction_factor_upper_bound": contraction.b,
        "H_a_on_box": ha_box,
        "H_xx_on_box": hxx_box,
        "analytic_error_bound_per_integral": analytic_error,
        "interval_precision_digits": iv.dps,
    }


def base_odd_primes(limit):
    """Tek sayı 2*i+1, i indisinde saklanır."""
    flags = bytearray(b"\x01") * ((limit + 1) // 2)
    if flags:
        flags[0] = 0
    for prime in range(3, isqrt(limit) + 1, 2):
        if flags[prime // 2]:
            start = prime * prime // 2
            count = (len(flags) - 1 - start) // prime + 1
            flags[start::prime] = b"\x00" * count
    return [2 * i + 1 for i, is_prime in enumerate(flags) if is_prime]


def primes_after_three(limit, odds_per_segment):
    """5'ten başlayıp üst sınırı da kapsayan bölümlü Eratosthenes eleği."""
    base_primes = base_odd_primes(isqrt(limit))
    last_odd = limit if limit % 2 else limit - 1
    for lo in range(5, last_odd + 1, 2 * odds_per_segment):
        hi = min(last_odd, lo + 2 * (odds_per_segment - 1))
        size = (hi - lo) // 2 + 1
        flags = bytearray(b"\x01") * size
        for prime in base_primes:
            if prime * prime > hi:
                break
            first = max(prime * prime, ((lo + prime - 1) // prime) * prime)
            if first % 2 == 0:
                first += prime
            if first <= hi:
                start = (first - lo) // 2
                count = (size - 1 - start) // prime + 1
                flags[start::prime] = b"\x00" * count

        index = flags.find(b"\x01")
        while index != -1:
            yield lo + 2 * index
            index = flags.find(b"\x01", index + 1)


def root_count(prime):
    """Asal prime > 3 için t³ = 2 (mod prime) denkleminin kök sayısı."""
    if prime % 3 == 2:
        return 1
    return 3 if pow(2, (prime - 1) // 3, prime) == 1 else 0


def check_small_primes():
    trial = [p for p in range(5, 201)
             if all(p % d for d in range(2, isqrt(p) + 1))]
    # Segment sınırlarını da sınamak için testte küçük segmentler kullanılır.
    assert list(primes_after_three(200, 7)) == trial
    assert list(primes_after_three(199, 7)) == trial
    for prime in trial:
        direct = sum((t**3 - 2) % prime == 0 for t in range(prime))
        assert root_count(prime) == direct, (prime, direct)
    assert root_count(7) == 0 and root_count(31) == 3
    return len(trial)


def compute_discrepancy(limit, odds_per_segment, points_per_decade):
    if limit < 10 or odds_per_segment < 1 or points_per_decade < 1:
        raise ValueError("LIMIT ≥ 10, segment ve nokta sayısı pozitif olmalı.")

    powers_of_ten = {10**k for k in range(1, int(log10(limit)) + 1) if 10**k <= limit}
    log_steps = (
        round(10 ** (k / points_per_decade))
        for k in range(points_per_decade, int(log10(limit) * points_per_decade) + 1)
    )
    checkpoints = sorted(powers_of_ten | {x for x in log_steps if 10 <= x <= limit} | {limit})
    rows = []
    next_point = 0
    counts = [0, 0, 0, 0]
    total = 0.0
    correction = 0.0

    for prime in primes_after_three(limit, odds_per_segment):
        while next_point < len(checkpoints) and checkpoints[next_point] < prime:
            x = checkpoints[next_point]
            rows.append((x, counts[0], counts[1], counts[3], total, total / sqrt(x)))
            next_point += 1

        roots = root_count(prime)
        counts[roots] += 1
        if roots != 1:
            # Kahan toplamı, kayan nokta yuvarlama hatasını azaltır.
            summand = (roots - 1) * log(prime) - correction
            updated = total + summand
            correction = (updated - total) - summand
            total = updated

    for x in checkpoints[next_point:]:
        rows.append((x, counts[0], counts[1], counts[3], total, total / sqrt(x)))

    return rows, powers_of_ten


def decimal_endpoint(endpoint):
    # mpmath bir aralığın uçlarını işaret, mantis ve 2'nin üssüyle saklar.
    # Bu bilgi tam ondalık yazıma çevrilir, kullanıcıya gösterilen yuvarlak çıktıya değil.
    sign, mantissa, exponent, bit_count = endpoint
    if bit_count < 0:
        raise ValueError("Aralık ucu sonlu olmalı.")
    if not mantissa:
        return "0"
    if exponent >= 0:
        digits = str(mantissa << exponent)
    else:
        places = -exponent
        digits = str(mantissa * 5**places).zfill(places + 1)
        digits = digits[:-places] + "." + digits[-places:]
    return "-" + digits if sign else digits


def encode_exact(value):
    if isinstance(value, Fraction):
        return str(value)
    if isinstance(value, type(iv.mpf(0))):
        lower, upper = value._mpi_
        return {"lower": decimal_endpoint(lower), "upper": decimal_endpoint(upper)}
    raise TypeError(f"JSON'a aktarılamayan tür {type(value).__name__}")


def prime_report():
    limit, odds_per_segment, points_per_decade = 10_000_000, 1_000_000, 30
    small_primes_checked = check_small_primes()
    start = perf_counter()
    rows, _ = compute_discrepancy(limit, odds_per_segment, points_per_decade)
    elapsed = perf_counter() - start
    return {
        "limit": limit,
        "odd_integers_per_segment": odds_per_segment,
        "points_per_decade": points_per_decade,
        "excluded_primes": [2, 3],
        "small_primes_checked_to_200": small_primes_checked,
        "prime_count": sum(rows[-1][1:4]),
        "elapsed_seconds": elapsed,
        "samples": [
            {
                "x": x,
                "zero_roots": roots_0,
                "one_root": roots_1,
                "three_roots": roots_3,
                "D": total,
                "D_over_sqrt_x": scaled,
            }
            for x, roots_0, roots_1, roots_3, total, scaled in rows
        ],
    }


def replay(include_prime_scan):
    if not __debug__:
        raise RuntimeError("Run without Python optimization to enable certificate assertions.")

    print("Recomputing the Laguerre certificate...", flush=True)
    values, laguerre, error = certify_laguerre()
    laguerre_report = {
        "scale": scale,
        "frequency": frequency,
        "H": values[0],
        "H_x": values[1],
        "H_xx": values[2],
        "laguerre_expression": laguerre,
        "analytic_error_bound_per_integral": error,
        "interval_precision_digits": iv.dps,
    }

    print("Recomputing the Rouche certificate...", flush=True)
    nonreal_report = certify_nonreal_zero()

    print("Recomputing the double-zero certificate...", flush=True)
    collision_report = certify_collision(
        "1.45200626937306544414442672465215400624761",
        "14.8380401130710086477956510630327906219061",
        [[-387250867, 0], [-948213121, 1602446079]],
    )

    report = {
        "schema_version": 1,
        "software": {"python": python_version(), "mpmath": mpmath.__version__},
        "laguerre_counterexample": laguerre_report,
        "nonreal_zero_pair": nonreal_report,
        "double_zero": collision_report,
    }
    if include_prime_scan:
        print("Recomputing the finite prime scan...", flush=True)
        report["prime_scan"] = prime_report()
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=Path("replay.json"))
    parser.add_argument("--include-prime-scan", action="store_true")
    args = parser.parse_args()
    report = replay(args.include_prime_scan)
    serialized = json.dumps(report, allow_nan=False, indent=2, default=encode_exact)
    args.output.write_text(serialized + "\n", encoding="utf-8")
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
