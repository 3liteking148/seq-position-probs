import sys
import os
import argparse
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
from scipy import stats

DISTS = {
    "gumbel": {"dist": stats.gumbel_r,   "names": ("loc", "scale")},
    "gev":    {"dist": stats.genextreme, "names": ("shape", "loc", "scale")},
}


def read_numbers(path):
    if path is None:
        stream = sys.stdin
    else:
        stream = open(path, "r", encoding="utf-8")
    vals = []
    with stream as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            line = line.replace(",", " ")
            for tok in line.split():
                vals.append(float(tok))
    x = np.asarray(vals, dtype=float)
    return x[np.isfinite(x)]


def fit_distribution(x, name):
    dist_class, names = DISTS[name]["dist"], DISTS[name]["names"]
    params = dist_class.fit(x)
    frozen = dist_class(*params)
    n = x.size
    param_str = ", ".join(f"{nm}={v:.6g}" for nm, v in zip(names, params))
    xs = np.linspace(x.min(), x.max(), 600)
    x_sorted = np.sort(x)
    surv_emp = 1 - (np.arange(1, n + 1) - 1) / n
    return {
        "frozen": frozen,
        "param_str": param_str,
        "xs": xs,
        "x_sorted": x_sorted,
        "surv_emp": surv_emp,
        "n": n,
    }


def plot_survival(ax, data, dist_name):
    frozen = data["frozen"]
    ax.plot(data["x_sorted"], data["surv_emp"], marker=".", linestyle="none", label="observed")
    ax.plot(data["xs"], frozen.sf(data["xs"]), linewidth=2, label="expected")
    ax.set_xlabel("similarity score")
    ax.set_ylabel("cumulative frequency")
    ax.yaxis.set_major_formatter(FuncFormatter(lambda v, pos: f"{v * data['n']:.0f}"))
    ax.legend()

    # --- hidden for now: histogram + fitted PDF ---
    # pdf = frozen.pdf(data["xs"])
    # ax.hist(x, bins="fd", density=True, alpha=0.6)
    # ax.plot(xs, pdf, linewidth=2)
    # ax.set_title("Histogram + fitted PDF")

    # --- hidden for now: Q-Q plot ---
    # p = (np.arange(1, n + 1) - 0.5) / n
    # theo_q = frozen.ppf(p)
    # ax.scatter(theo_q, x_sorted, s=14)
    # q1x, q3x = np.quantile(theo_q, [0.25, 0.75])
    # q1y, q3y = np.quantile(x_sorted, [0.25, 0.75])
    # slope = (q3y - q1y) / (q3x - q1x)
    # intercept = q1y - slope * q1x
    # line_x = np.array([theo_q.min(), theo_q.max()])
    # ax.plot(line_x, intercept + slope * line_x, linewidth=2)
    # ax.set_title("Q-Q plot")


def main():
    parser = argparse.ArgumentParser(
        description="Fit Gumbel/GEV to scores; one score file per family (label = basename). "
                    "With no files, reads numbers from stdin.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("score_files", nargs="*", help="Score files, one per family (up to 3 families per figure)")
    parser.add_argument("--title", type=str, default="Score distribution", help="Plot title")
    parser.add_argument("--fit", type=str, choices=sorted(DISTS), default="gumbel", help="Distribution to fit")
    parser.add_argument("--output", type=str, default=None, help="Output PNG path")
    args = parser.parse_args()

    if args.score_files:
        families = []
        for path in args.score_files:
            x = read_numbers(path)
            if x.size < 5:
                print(f"# Skipping {path}: only {x.size} finite data points")
                continue
            label = os.path.splitext(os.path.basename(path))[0]
            families.append((label, x))
        if not families:
            sys.exit(1)
    else:
        x = read_numbers(None)
        families = [("stdin", x)]

    fitted = [(label, fit_distribution(x, args.fit)) for label, x in families]

    n_fams = len(families)
    fig, axes = plt.subplots(n_fams, 1, sharex=True, figsize=(5, 5 * n_fams), squeeze=False)
    for ax, (label, data) in zip(axes[:, 0], fitted):
        plot_survival(ax, data, args.fit)
        ax.set_title(f"{label}\n{data['param_str']}")
        print(f"Fitted {args.fit} ({label}): {data['param_str']}")

    fig.suptitle(f"{args.title} ({args.fit} fit)")
    fig.tight_layout()

    if args.output:
        plt.savefig(args.output, dpi=150, bbox_inches="tight")
        print(f"Saved to {args.output}")
    else:
        plt.show()
    plt.close(fig)


if __name__ == "__main__":
    main()
