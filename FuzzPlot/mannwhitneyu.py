#!/usr/bin/env python3

import argparse
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu, describe
from a12 import a12s


def read_fuzz_log(prefix: str, key: str) -> list[int]:
    path = f"{prefix}{key}.log"
    values = []
    with open(path) as fp:
        for line in fp:
            values.extend(int(x) for x in line.split())
    return values


def basic_stats(name: str, values: list[int]) -> None:
    print(f"Basic stats for {name}")
    print(describe(values))


def compare_mwu(name1: str, name2: str, d: dict[str, list[int]]) -> None:
    print(f"MW U {name1} vs. {name2}")
    u_stat, p_value = mannwhitneyu(
        d[name1],
        d[name2],
        alternative="two-sided",
        method="auto",
    )
    print(f"U = {u_stat}, p = {p_value}\n")


def compute_a12(data: dict[str, list[int]]) -> None:
    rxs = [[name, *vals] for name, vals in data.items()]
    for thr in (0.56, 0.64, 0.71):
        print(f"A12 over {thr}")
        for rx in a12s(rxs, enough=thr):
            print(rx)
        print()


def plot_boxplots(data: dict[str, list[int]]) -> None:
    labels = list(data.keys())
    values = [data[k] for k in labels]

    plt.figure()
    plt.boxplot(values)
    plt.xticks(range(1, len(labels) + 1), labels, rotation=15)
    plt.ylabel("Coverage (Num. PCs)")
    plt.tight_layout()
    plt.savefig("coverage_boxplots.png", dpi=150)
    plt.close()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--prefix",
        default="fuzz-",
        help="Prefix for log files (default: %(default)s)",
    )
    parser.add_argument(
        "--experiments",
        nargs="+",
        default=["baseline", "exp1", "exp2", "exp3"],
        help="Experiment names to load",
    )
    args = parser.parse_args()

    data: dict[str, list[int]] = {}
    for key in args.experiments:
        data[key] = read_fuzz_log(args.prefix, key)
        basic_stats(key, data[key])

    # pairwise MWU against baseline
    if "baseline" in data:
        for k in args.experiments:
            if k != "baseline":
                compare_mwu("baseline", k, data)

    compute_a12(data)
    plot_boxplots(data)


if __name__ == "__main__":
    main()
