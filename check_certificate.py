"""Check the exported certificate inequalities using exact rational arithmetic."""

import argparse
import json
from fractions import Fraction
from math import isclose, isfinite, sqrt
from pathlib import Path


def require(condition, explanation):
    if not condition:
        raise ValueError(explanation)


def bounds(record):
    require(isinstance(record["lower"], str) and isinstance(record["upper"], str),
            "Interval endpoints must be strings encoding exact numbers.")
    lower = Fraction(record["lower"])
    upper = Fraction(record["upper"])
    require(lower <= upper, "Interval endpoints are reversed.")
    return lower, upper


def product_interval(first, second):
    endpoints = [a * b for a in first for b in second]
    return min(endpoints), max(endpoints)


def check_error(record, digits, degree, sup_bound):
    expected = 6 * sup_bound * Fraction(1, 8) ** (degree + 1) + Fraction(1, 10**80)
    require(record["interval_precision_digits"] == digits, "Unexpected interval precision.")
    require(Fraction(record["analytic_error_bound_per_integral"]) == expected,
            "The reported analytic error does not match its quadrature and tail formula.")


def check_laguerre(record):
    require(Fraction(record["scale"]) == Fraction(23, 16), "Unexpected Laguerre scale.")
    require(Fraction(record["frequency"]) == Fraction(74, 5), "Unexpected Laguerre frequency.")
    check_error(record, digits=45, degree=24, sup_bound=10**6)
    h_product = product_interval(bounds(record["H"]), bounds(record["H_xx"]))
    hx_lower, hx_upper = bounds(record["H_x"])
    hx_square_lower = 0 if hx_lower <= 0 <= hx_upper else min(hx_lower**2, hx_upper**2)
    hx_square_upper = max(hx_lower**2, hx_upper**2)
    lower, upper = bounds(record["laguerre_expression"])
    require(lower <= hx_square_lower - h_product[1], "Laguerre lower endpoint is too high.")
    require(hx_square_upper - h_product[0] <= upper < 0,
            "Laguerre upper endpoint does not certify the strict negative sign.")


def check_rouche(record):
    require(Fraction(record["scale"]) == Fraction(23, 16), "Unexpected Rouché scale.")
    check_error(record, digits=60, degree=40, sup_bound=10**8)
    center = record["positive_imaginary_center"]
    require(Fraction(center["real"]) == Fraction("14.67872714035986218218"),
            "Unexpected Rouché disk center.")
    imag_center = Fraction(center["imaginary"])
    require(imag_center == Fraction("0.31107587891741700810"), "Unexpected Rouché disk center.")
    radius = Fraction(record["radius"])
    require(radius == Fraction(1, 10**12) and radius < imag_center < 1 - radius,
            "The disk must lie strictly inside the upper half of the verified strip.")

    center_abs_upper = sum(
        max(abs(endpoint) for endpoint in bounds(record["H_at_center"][component]))
        for component in ("real", "imaginary")
    )
    slope_lower, _ = bounds(record["real_H_prime_at_center"])
    remainder_lower, remainder = bounds(record["rouche_remainder_upper_bound"])
    linear, linear_upper = bounds(record["rouche_linear_lower_bound"])
    require(remainder_lower == remainder and linear == linear_upper,
            "Expected scalar Rouché bounds.")
    require(center_abs_upper + 10 * radius**2 <= remainder,
            "The reported Rouché remainder bound is too small.")
    require(0 < linear <= slope_lower * radius,
            "The reported Rouché linear bound is too large or nonpositive.")
    require(remainder < linear, "Rouché's strict inequality did not hold.")
    require(record["zeros_per_disk_counting_multiplicity"] == 1,
            "The certified count is one zero per disk, counting multiplicity.")


def check_collision(record):
    check_error(record, digits=70, degree=48, sup_bound=10**7)
    radius = Fraction(record["box_radius"])
    require(radius == Fraction(1, 10**19), "Unexpected collision box radius.")
    expected_centers = {
        "scale": "1.45200626937306544414442672465215400624761",
        "frequency": "14.8380401130710086477956510630327906219061",
    }
    for coordinate, center_text in expected_centers.items():
        center = Fraction(record["center"][coordinate])
        require(center == Fraction(center_text), f"Unexpected collision {coordinate} center.")
        lower, upper = bounds(record[f"{coordinate}_interval"])
        require(lower <= center - radius and center + radius <= upper,
                f"The stored {coordinate} interval does not contain the rational box.")
    inverse = [[Fraction(value) for value in row] for row in record["inverse_jacobian_approximation"]]
    require(inverse == [[-387250867, 0], [-948213121, 1602446079]],
            "Unexpected collision preconditioner.")
    require(inverse[0][0] * inverse[1][1] - inverse[0][1] * inverse[1][0] != 0,
            "The collision preconditioner is singular.")

    delta_lower, delta = bounds(record["center_displacement_upper_bound"])
    k_lower, k = bounds(record["contraction_factor_upper_bound"])
    require(delta_lower == delta and k_lower == k, "Expected scalar upper bounds.")
    require(delta >= 0 and 0 <= k < 1, "Invalid nonnegative bound or contraction factor.")
    require(delta + k * radius < radius, "The contraction does not map the box inside itself.")
    require(bounds(record["H_a_on_box"])[1] < 0 < bounds(record["H_xx_on_box"])[0],
            "The local splitting derivative signs were not certified.")


def check_prime_scan(record):
    require(record["excluded_primes"] == [2, 3], "Unexpected excluded primes.")
    require(record["small_primes_checked_to_200"] == 44, "Unexpected small-prime test count.")
    require(record["odd_integers_per_segment"] > 0 and record["points_per_decade"] > 0,
            "Prime scan settings must be positive.")
    require(isfinite(record["elapsed_seconds"]) and record["elapsed_seconds"] >= 0,
            "The recorded scan time is invalid.")
    last_x = 0
    last_counts = (0, 0, 0)
    for row in record["samples"]:
        x = row["x"]
        counts = tuple(row[key] for key in ("zero_roots", "one_root", "three_roots"))
        require(type(x) is int and last_x < x <= record["limit"], "Invalid sample ordering.")
        require(all(type(count) is int and count >= previous
                    for count, previous in zip(counts, last_counts)),
                "Prime class counts must be nondecreasing integers.")
        require(isfinite(row["D"]) and isfinite(row["D_over_sqrt_x"])
                and isclose(row["D_over_sqrt_x"], row["D"] / sqrt(x), rel_tol=1e-12, abs_tol=1e-12),
                "A stored normalized prime sum is inconsistent.")
        last_x, last_counts = x, counts
    require(last_x == record["limit"] and sum(last_counts) == record["prime_count"],
            "Final sample and prime count disagree.")
    if record["limit"] == 10_000_000:
        require(last_counts == (221607, 332383, 110587),
                "Unexpected terminal class counts through 10,000,000.")


CHECKS = {
    "laguerre_counterexample": check_laguerre,
    "nonreal_zero_pair": check_rouche,
    "double_zero": check_collision,
}


def check_report(report):
    require(report["schema_version"] == 1, "Unsupported report schema.")
    require(all(isinstance(report["software"][name], str) for name in ("python", "mpmath")),
            "Incomplete software provenance.")
    for key, check in CHECKS.items():
        check(report[key])
    if "prime_scan" in report:
        check_prime_scan(report["prime_scan"])


def overlaps(first, second):
    a_lower, a_upper = bounds(first)
    b_lower, b_upper = bounds(second)
    return max(a_lower, b_lower) <= min(a_upper, b_upper)


def compare_reports(report, reference):
    interval_keys = {
        "laguerre_counterexample": ("H", "H_x", "H_xx", "laguerre_expression"),
        "nonreal_zero_pair": ("real_H_prime_at_center",),
        "double_zero": ("scale_interval", "frequency_interval", "H_a_on_box", "H_xx_on_box"),
    }
    for section, keys in interval_keys.items():
        for key in keys:
            require(overlaps(report[section][key], reference[section][key]),
                    f"Disjoint enclosures for {section}.{key}.")
    for component in ("real", "imaginary"):
        key = "nonreal_zero_pair"
        require(overlaps(report[key]["H_at_center"][component],
                         reference[key]["H_at_center"][component]),
                f"Disjoint enclosures of the {component} part at the disk center.")
    if "prime_scan" in report and "prime_scan" in reference:
        actual, expected = report["prime_scan"], reference["prime_scan"]
        for key in ("limit", "odd_integers_per_segment", "points_per_decade", "prime_count"):
            require(actual[key] == expected[key], f"Prime scan mismatch for {key}.")
        require(len(actual["samples"]) == len(expected["samples"]), "Sample counts differ.")
        for row, previous in zip(actual["samples"], expected["samples"]):
            for key in ("x", "zero_roots", "one_root", "three_roots"):
                require(row[key] == previous[key], f"Prime sample mismatch for {key}.")
            require(isclose(row["D"], previous["D"], rel_tol=1e-12, abs_tol=1e-12),
                    "The replayed weighted prime sum differs from the reference.")
    return all(report[key] == reference[key] for key in CHECKS)


def read_checked(path):
    report = json.loads(path.read_text(encoding="utf-8"))
    check_report(report)
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("report", type=Path)
    parser.add_argument("--compare", type=Path, metavar="REFERENCE")
    args = parser.parse_args()
    try:
        report = read_checked(args.report)
        if args.compare:
            exact = compare_reports(report, read_checked(args.compare))
    except (OSError, ValueError, TypeError, KeyError, ZeroDivisionError) as error:
        parser.exit(1, f"Certificate check failed: {error}\n")

    print("Passed: Laguerre sign, Rouche inequality, and double-zero contraction and signs.")
    if "prime_scan" in report:
        print("Passed: finite prime scan consistency checks.")
    if args.compare:
        print("All stored interval certificate fields exactly match the reference."
              if exact else "The replayed and reference enclosures overlap.")


if __name__ == "__main__":
    main()
