"""Persistent monitoring for the constrained Phase 7C search loop.

The monitor is solver-agnostic: a classical optimizer or an external search
human can append the same audited local-candidate record.  It writes JSONL for
lossless replay, a compact CSV for quick inspection, a one-line console status
and a Matplotlib dashboard with loss, gap/capture and parameter trajectories.
"""
from __future__ import annotations

from dataclasses import dataclass, field
import csv
import json
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np

from .phase7_microscopic_optimizer import candidate_from_payload, interaction_connectivity


METRIC_FIELDS = (
    "minimum_monomer_overlap_singular_value",
    "low_to_high_gap",
    "gap_over_mott_u",
    "maximum_low_state_mott_violation",
    "parity_commutator_frobenius_max",
    "observed_target_coefficient",
    "target_coefficient_ratio",
    "target_alignment",
    "fixed_scale_spectral_algebraic_residual",
    "best_scale_operator_residual",
    "unwanted_pauli_norm_over_target",
)


def _flat_parameters(candidate_payload: dict) -> dict[str, float]:
    candidate = candidate_from_payload(candidate_payload)
    out = {"mott_u": candidate.mott_u}
    for edge, bias in enumerate(candidate.rail_biases):
        out[f"bias_{edge}"] = bias
    for channel in candidate.channels:
        stem = f"channel.{channel.name}"
        out[f"{stem}.detuning"] = channel.detuning
        out[f"{stem}.coupling"] = channel.coupling
        out[f"{stem}.coupling_over_detuning"] = channel.coupling / channel.detuning
        out[f"{stem}.phase"] = channel.phase
    for position, hop in enumerate(candidate.rail_hops):
        out[f"hop.{position}.amplitude"] = hop.amplitude
        out[f"hop.{position}.phase"] = hop.phase
    for position, coupling in enumerate(candidate.zz_couplings):
        out[f"zz.{position}.strength"] = coupling.strength
    return out


def _record(iteration: int, audit: dict, tag: str) -> dict[str, Any]:
    state = audit["state_vector"]
    candidate = audit["candidate"]
    connectivity = interaction_connectivity(candidate_from_payload(candidate))
    return {
        "iteration": int(iteration),
        "tag": str(tag),
        "reward": float(audit["reward"]),
        "loss": float(state["fixed_scale_spectral_algebraic_residual"]),
        "target": audit["target"]["label"],
        "hard_failure_count": len(audit["hard_failures"]),
        "connected_patch": int(connectivity["patch_connected"]),
        "metric": {field: float(state[field]) for field in METRIC_FIELDS},
        "parameters": _flat_parameters(candidate),
    }


@dataclass
class OptimizationMonitor:
    """Append audited candidates and render a replayable phase-7C dashboard."""

    output_dir: Path
    reset: bool = False
    rows: list[dict[str, Any]] = field(default_factory=list, init=False)

    def __post_init__(self) -> None:
        self.output_dir = Path(self.output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.jsonl_path = self.output_dir / "history.jsonl"
        self.csv_path = self.output_dir / "history.csv"
        self.png_path = self.output_dir / "dashboard.png"
        if self.reset:
            for path in (self.jsonl_path, self.csv_path, self.png_path):
                path.unlink(missing_ok=True)
        elif self.jsonl_path.exists():
            self.rows = [json.loads(line) for line in self.jsonl_path.read_text(encoding="utf-8").splitlines() if line]

    def record(self, iteration: int, audit: dict, tag: str = "candidate") -> dict[str, Any]:
        row = _record(iteration, audit, tag)
        self.rows.append(row)
        with self.jsonl_path.open("a", encoding="utf-8") as stream:
            stream.write(json.dumps(row, sort_keys=True) + "\n")
        self._write_csv()
        self.render()
        print(
            f"iter={row['iteration']:04d} target={row['target']} loss={row['loss']:.6g} "
            f"reward={row['reward']:.6g} gap/U={row['metric']['gap_over_mott_u']:.6g} "
            f"capture={row['metric']['minimum_monomer_overlap_singular_value']:.6g} "
            f"connected={row['connected_patch']} failures={row['hard_failure_count']}"
        )
        return row

    def _write_csv(self) -> None:
        parameter_names = sorted({name for row in self.rows for name in row["parameters"]})
        fieldnames = ["iteration", "tag", "target", "loss", "reward", "hard_failure_count", "connected_patch"]
        fieldnames += list(METRIC_FIELDS) + parameter_names
        with self.csv_path.open("w", newline="", encoding="utf-8") as stream:
            writer = csv.DictWriter(stream, fieldnames=fieldnames)
            writer.writeheader()
            for row in self.rows:
                flat = {key: row[key] for key in ("iteration", "tag", "target", "loss", "reward", "hard_failure_count", "connected_patch")}
                flat.update(row["metric"])
                flat.update(row["parameters"])
                writer.writerow(flat)

    def render(self) -> None:
        if not self.rows:
            return
        iterations = np.asarray([row["iteration"] for row in self.rows])
        loss = np.asarray([row["loss"] for row in self.rows])
        reward = np.asarray([row["reward"] for row in self.rows])
        gaps = np.asarray([row["metric"]["gap_over_mott_u"] for row in self.rows])
        capture = np.asarray([row["metric"]["minimum_monomer_overlap_singular_value"] for row in self.rows])
        target_ratio = np.asarray([row["metric"]["target_coefficient_ratio"] for row in self.rows])
        unwanted = np.asarray([row["metric"]["unwanted_pauli_norm_over_target"] for row in self.rows])
        all_parameter_names = {name for row in self.rows for name in row["parameters"]}
        priority = [
            "mott_u",
            *sorted(name for name in all_parameter_names if name.endswith("coupling_over_detuning")),
            *sorted(name for name in all_parameter_names if name.endswith("detuning")),
            *sorted(name for name in all_parameter_names if name.endswith("coupling")),
            *sorted(name for name in all_parameter_names if name.endswith("phase")),
            *sorted(name for name in all_parameter_names if name.startswith("hop.") or name.startswith("zz.")),
            *sorted(name for name in all_parameter_names if name.startswith("bias_") and any(abs(row["parameters"].get(name, 0.0)) > 1e-15 for row in self.rows)),
        ]
        parameter_names = list(dict.fromkeys(priority))[:7]
        figure, axes = plt.subplots(2, 2, figsize=(12, 7), constrained_layout=True)
        axes[0, 0].plot(iterations, loss, marker="o", label="loss")
        axes[0, 0].plot(iterations, -reward, marker="s", label="-reward")
        axes[0, 0].set_title("Loss / reward")
        axes[0, 0].set_xlabel("iteration")
        axes[0, 0].legend()
        axes[0, 1].plot(iterations, gaps, marker="o", label="gap / U")
        axes[0, 1].plot(iterations, capture, marker="s", label="minimum capture")
        axes[0, 1].set_title("Monomer manifold")
        axes[0, 1].set_xlabel("iteration")
        axes[0, 1].legend()
        axes[1, 0].plot(iterations, target_ratio, marker="o", label="target coefficient / (-J/2)")
        axes[1, 0].plot(iterations, unwanted, marker="s", label="unwanted / target")
        axes[1, 0].set_title("Operator selectivity")
        axes[1, 0].set_xlabel("iteration")
        axes[1, 0].legend()
        for name in parameter_names:
            values = [row["parameters"].get(name, np.nan) for row in self.rows]
            axes[1, 1].plot(iterations, values, marker=".", label=name)
        axes[1, 1].set_title("Registered parameter trajectories")
        axes[1, 1].set_xlabel("iteration")
        if parameter_names:
            axes[1, 1].legend(fontsize=7, ncol=2)
        figure.savefig(self.png_path, dpi=160)
        plt.close(figure)
