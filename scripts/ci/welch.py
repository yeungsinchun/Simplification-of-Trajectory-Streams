import math
import sys


T95 = {
    1: 6.314, 2: 2.920, 3: 2.353, 4: 2.132, 5: 2.015, 6: 1.943,
    7: 1.895, 8: 1.860, 9: 1.833, 10: 1.812, 11: 1.796, 12: 1.782,
    13: 1.771, 14: 1.761, 15: 1.753, 16: 1.746, 17: 1.740,
    18: 1.734, 19: 1.729, 20: 1.725, 21: 1.721, 22: 1.717,
    23: 1.714, 24: 1.711, 25: 1.708, 26: 1.706, 27: 1.703,
    28: 1.701, 29: 1.699, 30: 1.697, 40: 1.684, 60: 1.671,
    120: 1.658, 1000: 1.645,
}


def t95(df):
    if df <= 1:
        return T95[1]
    if df >= 1000:
        return T95[1000]
    keys = sorted(T95)
    for lower, upper in zip(keys, keys[1:]):
        if lower <= df <= upper:
            return T95[lower] + (df - lower) * (T95[upper] - T95[lower]) / (upper - lower)
    raise ValueError("invalid degrees of freedom")


def verdict(n, limit, orig_means, orig_stds, new_means, new_stds):
    m = len(orig_means)
    if n < 2 or m == 0 or any(len(values) != m for values in (orig_stds, new_means, new_stds)):
        raise ValueError("expected equal nonempty samples with at least two runs")
    if not all(math.isfinite(value) for value in (limit, *orig_means, *orig_stds, *new_means, *new_stds)):
        raise ValueError("expected finite benchmark statistics")
    if limit <= 0 or any(std < 0 for std in (*orig_stds, *new_stds)):
        raise ValueError("expected positive limit and nonnegative standard deviations")

    difference = (sum(new_means) - limit * sum(orig_means)) / m
    variances = [std * std / (n * m * m) for std in new_stds]
    variances += [(limit * std) ** 2 / (n * m * m) for std in orig_stds]
    variance = sum(variances)
    if variance == 0:
        return difference > 0
    df = variance * variance / sum(component * component / (n - 1) for component in variances)
    return difference - t95(df) * math.sqrt(variance) > 0


if __name__ == "__main__":
    if len(sys.argv) != 7:
        sys.exit("usage: welch.py runs limit orig_means orig_stds new_means new_stds")
    try:
        runs = int(sys.argv[1])
        limit = float(sys.argv[2])
        groups = [[float(value) for value in arg.split()] for arg in sys.argv[3:]]
        print(int(verdict(runs, limit, *groups)))
    except ValueError as error:
        sys.exit(str(error))
