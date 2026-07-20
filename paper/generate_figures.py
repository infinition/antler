"""Generate the figures and numerical tables used by the ANTLER paper."""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
PAPER = ROOT / "paper"
FIGURES = PAPER / "figures"


def load(relative: str) -> dict:
    return json.loads((ROOT / relative).read_text(encoding="utf-8"))


def configure() -> None:
    plt.rcParams.update(
        {
            "font.size": 10,
            "axes.labelsize": 10,
            "axes.titlesize": 11,
            "legend.fontsize": 8,
            "figure.dpi": 150,
            "savefig.dpi": 300,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )


def save(fig: plt.Figure, name: str) -> None:
    fig.tight_layout()
    fig.savefig(FIGURES / f"{name}.pdf", bbox_inches="tight")
    fig.savefig(FIGURES / f"{name}.png", bbox_inches="tight")
    plt.close(fig)


def protocol_schematic() -> None:
    fig, ax = plt.subplots(figsize=(7.1, 2.55))
    length = 14
    upper, lower = 1.10, 0.38
    blue = "#1769aa"
    red = "#c62828"
    for rail, y in (("upper rail", upper), ("lower rail", lower)):
        ax.plot([0, length - 1], [y, y], color="0.35", lw=1.4)
        ax.scatter(np.arange(length), np.full(length, y), s=15, color="0.25", zorder=3)
        ax.text(-0.35, y, rail, ha="right", va="center", fontsize=8)
    ax.scatter([0, 0], [upper, lower], s=66, facecolors="none", edgecolors=blue, lw=1.6, zorder=5)
    ax.scatter([length - 1, length - 1], [upper, lower], s=66, facecolors="none", edgecolors=red, lw=1.6, zorder=5)

    arrow = {"arrowstyle": "-|>", "color": blue, "lw": 2.0, "mutation_scale": 13}
    ax.annotate("", xy=(3.86, upper), xytext=(0.27, upper), arrowprops=arrow, zorder=6)
    ax.annotate("", xy=(4.00, lower + 0.06), xytext=(4.00, upper - 0.06), arrowprops=arrow, zorder=6)
    ax.annotate("", xy=(0.27, lower), xytext=(3.73, lower), arrowprops=arrow, zorder=6)
    ax.annotate("", xy=(0.00, upper - 0.06), xytext=(0.00, lower + 0.06), arrowprops=arrow, zorder=6)

    ax.text(1.85, 1.48, "1  outbound", color=blue, ha="center", va="bottom")
    ax.text(4.20, 0.77, "2  first handoff", color=blue, ha="left", va="center")
    ax.text(2.05, -0.02, "3  return", color=blue, ha="center", va="top")
    ax.text(0.25, 0.77, "4  second handoff", color=blue, ha="left", va="center")
    ax.text(4.00, 1.27, "turning rung", ha="center", va="bottom", fontsize=8)
    ax.text(0.00, 1.82, "mobile logical branch", color=blue, ha="left")
    ax.text(length - 1, 1.82, "spectator logical branch", color=red, ha="right")
    ax.set_xlim(-0.9, length - 0.1)
    ax.set_ylim(-0.20, 1.98)
    ax.set_axis_off()
    save(fig, "protocol_schematic")


def deep_limit_figure(closure: dict) -> None:
    rows = closure["rows"]
    depth = np.asarray([row["depth"] for row in rows], dtype=float)
    slope = np.asarray([row["odd_slope"] for row in rows], dtype=float)
    fit = closure["fit"]
    grid = np.linspace(depth.min(), depth.max(), 400)
    predicted = -1.0 + fit["a"] * grid ** (-fit["p"])
    fig, (left, right) = plt.subplots(1, 2, figsize=(7.1, 3.0), gridspec_kw={"width_ratios": [1.08, 1.0]})
    left.loglog(depth, 1.0 + slope, "o", color="#1f77b4", label="ED data")
    left.loglog(grid, 1.0 + predicted, "-", color="#d62728", label="fixed-limit fit")
    left.set_xlabel("trap depth $D$")
    left.set_ylabel("$1+\\Delta\\phi/\\theta$")
    left.set_title("Finite-depth correction")
    left.legend(frameon=False, loc="upper right")
    left.grid(alpha=0.25, which="both")
    right.plot(depth, slope, "o", color="#1f77b4", label="ED data")
    right.plot(grid, predicted, "-", color="#d62728", label="fit")
    right.axhline(-1.0, color="0.35", lw=1.0, ls="--", label="digital limit")
    right.set_xlabel("trap depth $D$")
    right.set_ylabel("$\\Delta\\phi/\\theta$")
    right.set_title("Approach to the digital limit")
    right.set_ylim(-1.004, -0.967)
    right.grid(alpha=0.25)
    save(fig, "deep_limit")


def path_figure(path: dict) -> None:
    rows = [row for row in path["rows"] if row["variant"] != "baseline"]
    names = [row["variant"].replace("_", " ") for row in rows]
    phase = 1.0e5 * np.asarray([row["delta_phase_to_baseline"] for row in rows])
    distance = 1.0e5 * np.asarray([row["unitary_distance_to_baseline"] for row in rows])
    fig, (left, right) = plt.subplots(1, 2, figsize=(7.1, 3.45), sharey=True)
    y = np.arange(len(rows))
    left.axvline(0.0, color="0.35", lw=0.8)
    left.plot(phase, y, "o", color="#1f77b4")
    left.set_xlabel("phase shift from baseline ($10^{-5}$ rad)")
    left.set_yticks(y, names)
    left.set_title("Registered deformations")
    left.grid(axis="x", alpha=0.25)
    right.plot(distance, y, "o", color="#d62728")
    right.set_xlabel("relative logical-unitary distance ($10^{-5}$)")
    right.set_title("Change in reconstructed gate")
    right.grid(axis="x", alpha=0.25)
    save(fig, "path_invariance")


def composition_figure(composition: dict) -> None:
    rows = composition["rows"]
    cycles = np.asarray([row["n"] for row in rows])
    phase_error = np.asarray([row["phase_additivity_error"] for row in rows])
    leakage = np.asarray([row["leak_worst"] for row in rows])
    coherent = np.asarray([row["coherent_distance_to_power_of_one_cycle"] for row in rows])
    fig, (left, right) = plt.subplots(1, 2, figsize=(7.1, 3.0))
    finite = cycles > 1
    left.semilogy(cycles[finite], phase_error[finite], "o-", color="#1f77b4", label="phase additivity error")
    left.semilogy(cycles[finite], coherent[finite], "s--", color="#d62728", label="distance to one-cycle power")
    left.axvline(1, color="0.55", lw=0.8, ls=":")
    left.text(1.12, 1.6e-6, "$n=1$: zero by definition", fontsize=8, va="bottom")
    left.set_xlabel("number of cycles $n$")
    left.set_ylabel("error")
    left.set_xticks(cycles)
    left.set_ylim(1e-6, 5e-2)
    left.set_title("Composition")
    left.grid(alpha=0.25, which="both")
    left.legend(frameon=False)
    right.plot(cycles, leakage, "o-", color="#2ca02c")
    right.set_xlabel("number of cycles $n$")
    right.set_ylabel("worst-case leakage")
    right.set_xticks(cycles)
    right.set_title("Leakage under composition")
    right.grid(alpha=0.25)
    save(fig, "composition")


def latex_tables(closure: dict, composition: dict, path: dict) -> None:
    deep_rows = []
    for row in closure["rows"]:
        deep_rows.append(
            f"{row['depth']} & {row['T']:.0f} & {row['odd_slope']:.6f} & "
            f"{row['leak_worst']:.2e} & {row['sigma_min']:.6f} & "
            f"{row['minimum_handoff_isolation_gap']:.2e} \\\\"
        )
    time_rows = []
    for depth in (6, 8):
        for dt in (0.5, 0.25, 0.125):
            directory = "deep_limit" if dt == 0.25 else "deep_limit/dt_convergence"
            data = load(f"results/phase4_7/{directory}/D{depth}_dt{str(dt).replace('.', 'p')}.json")
            metrics = data["metrics"]
            time_rows.append(
                f"{depth} & {dt:.3f} & {metrics['odd_slope']:.6f} & "
                f"{metrics['leak_worst']:.2e} & {metrics['sigma_min']:.6f} \\\\"
            )
    composition_rows = []
    for row in composition["rows"]:
        composition_rows.append(
            f"{row['n']} & {row['phase_additivity_error']:.3e} & "
            f"{row['leak_worst']:.2e} & {row['coherent_distance_to_power_of_one_cycle']:.2e} & "
            f"{row['axis_drift_from_z']:.2e} \\\\"
        )
    path_rows = []
    for key, value in (
        ("largest phase shift", path["max_abs_phase_shift_from_baseline"]),
        ("largest logical-unitary distance", path["max_unitary_distance_from_baseline"]),
        ("largest error to $-\\theta$", path["max_abs_phase_error_to_minus_theta"]),
    ):
        path_rows.append(f"{key} & {value:.3e} \\\\"
        )
    content = "\n".join(
        [
            "% Generated by paper/generate_figures.py. Do not edit by hand.",
            "\\begin{table}[t]",
            "\\centering",
            "\\caption{Fine-step deep-limit data at $\\theta=0.3$ and $\\Delta t=0.125$. The protocol duration follows $T(D)=20000(D/4)^2$.}",
            "\\label{tab:deep-limit}",
            "\\begin{tabular}{rrrrrr}",
            "\\toprule",
            "$D$ & $T$ & $\\Delta\\phi/\\theta$ & leakage & $\\sigma_{\\min}$ & handoff gap \\\\",
            "\\midrule",
            *deep_rows,
            "\\bottomrule",
            "\\end{tabular}",
            "\\end{table}",
            "",
            "\\begin{table}[t]",
            "\\centering",
            "\\caption{Independent time-step check. The $\\Delta t=0.5$ rows are retained as a resolved coarse-step failure, rather than being used in the closure fit.}",
            "\\label{tab:timestep}",
            "\\begin{tabular}{rrrrr}",
            "\\toprule",
            "$D$ & $\\Delta t$ & $\\Delta\\phi/\\theta$ & leakage & $\\sigma_{\\min}$ \\\\",
            "\\midrule",
            *time_rows,
            "\\bottomrule",
            "\\end{tabular}",
            "\\end{table}",
            "",
            "\\begin{table}[t]",
            "\\centering",
            "\\caption{Repeated-cycle composition at $D=6$, $\\theta=0.3$ and $\\Delta t=0.25$.}",
            "\\label{tab:composition}",
            "\\begin{tabular}{rrrrr}",
            "\\toprule",
            "$n$ & phase error (rad) & leakage & one-cycle-power distance & $Z$-axis drift \\\\",
            "\\midrule",
            *composition_rows,
            "\\bottomrule",
            "\\end{tabular}",
            "\\end{table}",
            "",
            "\\begin{table}[t]",
            "\\centering",
            "\\caption{Summary of the eleven finite-depth protocol deformations at $D=8$, $\\theta=0.3$ and $\\Delta t=0.25$.}",
            "\\label{tab:path}",
            "\\begin{tabular}{lr}",
            "\\toprule",
            "quantity & maximum value \\\\",
            "\\midrule",
            *path_rows,
            "\\bottomrule",
            "\\end{tabular}",
            "\\end{table}",
            "",
        ]
    )
    (PAPER / "generated_results.tex").write_text(content, encoding="utf-8")


def main() -> None:
    FIGURES.mkdir(parents=True, exist_ok=True)
    configure()
    closure = load("results/phase4_7/publication_closure.json")
    path = load("results/phase4_7/path_invariance/summary.json")
    composition = load("results/phase4_7/composition/summary.json")
    protocol_schematic()
    deep_limit_figure(closure)
    path_figure(path)
    composition_figure(composition)
    latex_tables(closure, composition, path)


if __name__ == "__main__":
    main()
