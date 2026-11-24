"""Compare bp_fem_solver S-parameter outputs with the HFSS TFE reference.

The HFSS reference file `wg_bp_filter_S_parameters.csv` carries:
    "Freq [GHz]","mag(S(1,1)) []","mag(S(2,1)) []"

The solver writes `s_parameters.csv` with:
    freq_Hz,S11_real,S11_imag,S11_dB,S21_real,S21_imag,S21_dB

The script accepts any number of solver CSVs (each with a label) and overlays
them against the HFSS reference. It reports:
- max / mean / RMS |Δ| in dB on each curve
- max / mean / RMS |Δ| on linear magnitude (less sensitive to deep-null
  dB blowup, often the more honest comparison metric for bandpass filters)

Usage:
    python compare_with_hfss.py <hfss_csv> <out_dir>
        --solver <label>=<csv> [--solver <label>=<csv> ...]
        [--baseline <label>=<old_csv>]

    python compare_with_hfss.py <hfss_csv> <out_dir> --batch-root result
        # scans immediate child directories containing s_parameters.csv
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import sys
import tempfile
import uuid
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


def _load_solver(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, float_precision="round_trip").copy()
    df["freq_GHz"] = df["freq_Hz"].astype(float) / 1e9
    df["S11_mag"] = np.hypot(df["S11_real"], df["S11_imag"])
    df["S21_mag"] = np.hypot(df["S21_real"], df["S21_imag"])
    return df.sort_values("freq_GHz").reset_index(drop=True)


def _validate_run(path, data, allow_legacy):
    if path.parent.name.startswith(".bpfem-pending-") or os.path.lexists(path.parent / "run.json.pending"):
        raise ValueError(f"{path}: unfinished or failed run")
    frequencies = data["freq_Hz"].to_numpy(dtype=float)
    if (not len(frequencies) or not np.all(np.isfinite(data[["freq_Hz", "S11_mag", "S21_mag"]]))
            or frequencies[0] <= 0 or np.any(np.diff(frequencies) <= 0)):
        raise ValueError(f"{path}: expected finite S parameters at distinct positive frequencies")
    sidecar = path.parent / "run.json"
    report = json.loads(sidecar.read_text(encoding="utf-8")) if sidecar.exists() else None
    if report is None:
        if allow_legacy:
            return {"csv": str(path.resolve()), "verification": "legacy-csv-only"}
        raise ValueError(f"{path}: missing run.json; use --allow-legacy for historical CSVs")
    if report.get("status") != "ok":
        raise ValueError(f"{path}: solver run did not succeed")
    samples = report.get("sweep")
    modern = report.get("schema", 1) >= 2 or "inputs" in report
    if modern or samples is not None:
        if not isinstance(samples, list) or not np.array_equal(
                frequencies, sorted(point["freq_hz"] for point in samples)):
            raise ValueError(f"{path}: frequency samples differ between CSV and run.json")
    options = report["options"]
    port = options["portMethod"]
    physical = {"basis_order": int(options["basisOrder"]), "port_method": port,
                "tfe_modes_per_port": int(options["tfeModesPerPort"]) if port == "tfe" else None}
    verification = {"csv": str(path.resolve()), "physical": physical, "paths": {
        key: str(Path(report["cwd"]) / options[key]) for key in ("aedt", "mesh")
    }}
    if modern:
        fingerprints = report.get("inputs") or {}
        for key in ("aedt_sha256", "mesh_sha256"):
            if not re.fullmatch("[0-9a-f]{64}", str(fingerprints.get(key, "")).lower()):
                raise ValueError(f"{path}: missing or invalid recorded input SHA-256")
        verification.update(verification="sha256", inputs=fingerprints)
    elif allow_legacy:
        verification["verification"] = "legacy-path-only"
    else:
        raise ValueError(f"{path}: report has no input hashes; use --allow-legacy for historical runs")
    return verification


def _check_compatibility(records, grids, allow_physical_variation, interpolate):
    with_reports = [item for item in records.values() if "physical" in item]
    if with_reports:
        reference = next((item for item in with_reports if "inputs" in item), with_reports[0])
        for run in with_reports:
            if "inputs" in reference and "inputs" in run:
                matched = reference["inputs"] == run["inputs"]
            else:
                matched = all(os.path.normcase(os.path.normpath(reference["paths"][key])) ==
                              os.path.normcase(os.path.normpath(run["paths"][key])) for key in ("aedt", "mesh"))
            if not matched:
                raise ValueError(f"{run['csv']}: input identity differs from {reference['csv']}")
            if run["physical"] != reference["physical"] and not allow_physical_variation:
                raise ValueError(f"{run['csv']}: physical settings differ; use --allow-physical-variation deliberately")
    if not interpolate and any(not np.array_equal(grids[0], grid) for grid in grids[1:]):
        raise ValueError("Solver frequency grids differ; use --interpolate for an interpolated comparison")


def _load_hfss(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    df.columns = [c.strip().strip('"') for c in df.columns]
    rename = {}
    for c in df.columns:
        cl = c.lower()
        if "freq" in cl:
            rename[c] = "freq_GHz"
        elif "s(1,1)" in cl:
            rename[c] = "S11_mag"
        elif "s(2,1)" in cl:
            rename[c] = "S21_mag"
    df = df.rename(columns=rename)[["freq_GHz", "S11_mag", "S21_mag"]].copy()
    return df.sort_values("freq_GHz").reset_index(drop=True)


def _to_db(mag: np.ndarray) -> np.ndarray:
    return 20.0 * np.log10(np.maximum(mag, 1e-300))


def _interp(grid: np.ndarray, df: pd.DataFrame, col: str) -> np.ndarray:
    return np.interp(grid, df["freq_GHz"].to_numpy(), df[col].to_numpy())


def _interp_complex(grid: np.ndarray, df: pd.DataFrame, prefix: str) -> np.ndarray:
    """Interpolate a complex solver trace component-wise on the comparison grid."""
    real = np.interp(grid, df["freq_GHz"].to_numpy(), df[f"{prefix}_real"].to_numpy())
    imag = np.interp(grid, df["freq_GHz"].to_numpy(), df[f"{prefix}_imag"].to_numpy())
    return real + 1j * imag


def _stats(name: str, diff: np.ndarray, unit: str) -> str:
    return (
        f"  {name}: max|Δ|={np.max(np.abs(diff)):.4f} {unit}, "
        f"mean|Δ|={np.mean(np.abs(diff)):.4f} {unit}, "
        f"RMS|Δ|={np.sqrt(np.mean(diff**2)):.4f} {unit}"
    )


def main(argv: list[str]) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("hfss_csv", type=Path, help="HFSS reference CSV")
    parser.add_argument("out_dir", type=Path, help="Output directory")
    parser.add_argument(
        "--solver",
        action="append",
        default=[],
        help="Solver run as label=path/to/csv. Repeatable.",
    )
    parser.add_argument(
        "--baseline",
        action="append",
        default=[],
        help="Optimized-vs-baseline complex comparison as label=path/to/s_parameters.csv.",
    )
    parser.add_argument(
        "--batch-root",
        type=Path,
        help="Scan immediate child directories for s_parameters.csv and compare all runs.",
    )
    parser.add_argument("--allow-legacy", action="store_true", help="Allow historical runs without input hashes.")
    parser.add_argument("--allow-physical-variation", action="store_true",
                        help="Allow deliberate changes to basis order or port formulation.")
    parser.add_argument("--interpolate", action="store_true",
                        help="Allow different solver frequency grids and interpolate on the HFSS overlap.")
    args = parser.parse_args(argv[1:])

    hfss = _load_hfss(args.hfss_csv)
    if hfss.empty or not np.all(np.isfinite(hfss)) or np.any(np.diff(hfss["freq_GHz"]) <= 0):
        parser.error("HFSS reference must have finite data at distinct frequencies")
    source_records = {}

    def load_run(label, filename):
        if label in source_records:
            parser.error(f"Repeated run label: {label}")
        path = Path(filename)
        data = _load_solver(path)
        source_records[label] = _validate_run(path, data, args.allow_legacy)
        return data

    baselines: dict[str, pd.DataFrame] = {}
    for spec in args.baseline:
        if "=" not in spec:
            parser.error(f"--baseline must be label=path, got {spec}")
        label, path = spec.split("=", 1)
        baselines[label] = load_run(f"baseline:{label}", path)
    solver_runs = []
    solver_specs = list(args.solver)
    if args.batch_root is not None:
        for child in sorted(args.batch_root.iterdir()):
            csv_path = child / "s_parameters.csv"
            if child.is_dir() and not child.name.startswith(".") and csv_path.exists():
                label = child.name
                if label.startswith("result_"):
                    label = label[len("result_") :]
                elif label.startswith("results_"):
                    label = label[len("results_") :]
                solver_specs.append(f"{label}={csv_path}")
    if not solver_specs:
        parser.error("provide at least one --solver label=csv or use --batch-root")
    for spec in solver_specs:
        if "=" not in spec:
            parser.error(f"--solver must be label=path, got {spec}")
        label, path = spec.split("=", 1)
        solver_runs.append((label, load_run(label, path)))

    all_runs = [run for _, run in solver_runs] + list(baselines.values())
    _check_compatibility(source_records, [run["freq_Hz"].to_numpy() for run in all_runs],
                         args.allow_physical_variation, args.interpolate)
    f_lo = max(hfss["freq_GHz"].min(), max(s["freq_GHz"].min() for s in all_runs))
    f_hi = min(hfss["freq_GHz"].max(), min(s["freq_GHz"].max() for s in all_runs))
    grid = hfss[(hfss["freq_GHz"] >= f_lo) & (hfss["freq_GHz"] <= f_hi)]["freq_GHz"].to_numpy()
    if not len(grid):
        parser.error("HFSS and solver runs have no common frequency region")
    points_match = all(np.isin(grid, run["freq_GHz"].to_numpy()).all() for run in all_runs)
    if not points_match and not args.interpolate:
        parser.error("HFSS has unmatched points in the common region; use --interpolate to compare them")

    preferred = args.out_dir.resolve()
    preferred.parent.mkdir(parents=True, exist_ok=True)
    out_dir = Path(tempfile.mkdtemp(prefix=".bpfem-compare-pending-", dir=preferred.parent))

    s11_ref = _interp(grid, hfss, "S11_mag")
    s21_ref = _interp(grid, hfss, "S21_mag")
    s11_ref_db = _to_db(s11_ref)
    s21_ref_db = _to_db(s21_ref)

    summary_lines = [
        f"HFSS reference: {args.hfss_csv}",
        f"Compared band: {f_lo:.4f} - {f_hi:.4f} GHz, {len(grid)} points",
        f"Frequency comparison: {'exact points' if points_match else 'interpolated overlap'}",
        "",
    ]

    def _three_db_band(f_arr: np.ndarray, s21_db: np.ndarray):
        peak = int(np.argmax(s21_db))
        thr = s21_db[peak] - 3.0
        inside = np.where(s21_db >= thr)[0]
        if inside.size == 0:
            return float(f_arr[peak]), np.nan, np.nan
        return float(f_arr[peak]), float(f_arr[inside.min()]), float(f_arr[inside.max()])

    pf, lo, hi = _three_db_band(grid, s21_ref_db)
    summary_lines.append(
        f"HFSS |S21| peak @ {pf:.4f} GHz; 3-dB band {lo:.4f} - {hi:.4f} GHz (center {(lo+hi)/2:.4f} GHz)"
    )
    summary_lines.append("")

    out_table = pd.DataFrame({"freq_GHz": grid, "S11_hfss": s11_ref, "S21_hfss": s21_ref})
    summary_records = []

    # Plot 1: |S11| dB and |S21| dB overlay.
    fig, axes = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
    axes[0].plot(grid, s11_ref_db, label="HFSS TFE (ref)", linewidth=1.8, color="black")
    axes[1].plot(grid, s21_ref_db, label="HFSS TFE (ref)", linewidth=1.8, color="black")

    cycle = ["#d62728", "#1f77b4", "#2ca02c", "#9467bd", "#ff7f0e"]
    for idx, (label, df) in enumerate(solver_runs):
        s11 = _interp(grid, df, "S11_mag")
        s21 = _interp(grid, df, "S21_mag")
        s11_db = _to_db(s11)
        s21_db = _to_db(s21)
        color = cycle[idx % len(cycle)]
        axes[0].plot(grid, s11_db, label=label, linewidth=1.1, linestyle="--", color=color)
        axes[1].plot(grid, s21_db, label=label, linewidth=1.1, linestyle="--", color=color)

        out_table[f"S11_{label}"] = s11
        out_table[f"S21_{label}"] = s21
        out_table[f"S11_diff_dB_{label}"] = s11_db - s11_ref_db
        out_table[f"S21_diff_dB_{label}"] = s21_db - s21_ref_db
        out_table[f"S11_diff_lin_{label}"] = s11 - s11_ref
        out_table[f"S21_diff_lin_{label}"] = s21 - s21_ref

        summary_lines.append(f"[{label}]")
        pf_s, lo_s, hi_s = _three_db_band(grid, s21_db)
        summary_lines.append(
            f"  |S21| peak @ {pf_s:.4f} GHz; 3-dB band {lo_s:.4f} - {hi_s:.4f} GHz "
            f"(center {(lo_s+hi_s)/2:.4f} GHz, shift vs HFSS = {(lo_s+hi_s)/2 - (lo+hi)/2:+.4f} GHz)"
        )
        summary_lines.append(_stats("|S11| dB", s11_db - s11_ref_db, "dB"))
        summary_lines.append(_stats("|S21| dB", s21_db - s21_ref_db, "dB"))
        summary_lines.append(_stats("|S11| lin", s11 - s11_ref, ""))
        summary_lines.append(_stats("|S21| lin", s21 - s21_ref, ""))
        complex_metrics = {
            "baseline_max_abs_complex_s11": np.nan,
            "baseline_relative_l2_complex_s11": np.nan,
            "baseline_max_abs_complex_s21": np.nan,
            "baseline_relative_l2_complex_s21": np.nan,
        }
        if label in baselines:
            optimized_s11 = _interp_complex(grid, df, "S11")
            optimized_s21 = _interp_complex(grid, df, "S21")
            baseline_s11 = _interp_complex(grid, baselines[label], "S11")
            baseline_s21 = _interp_complex(grid, baselines[label], "S21")
            delta_s11 = optimized_s11 - baseline_s11
            delta_s21 = optimized_s21 - baseline_s21
            complex_metrics = {
                "baseline_max_abs_complex_s11": float(np.max(np.abs(delta_s11))),
                "baseline_relative_l2_complex_s11": float(
                    np.linalg.norm(delta_s11) / max(np.linalg.norm(baseline_s11), 1e-300)
                ),
                "baseline_max_abs_complex_s21": float(np.max(np.abs(delta_s21))),
                "baseline_relative_l2_complex_s21": float(
                    np.linalg.norm(delta_s21) / max(np.linalg.norm(baseline_s21), 1e-300)
                ),
            }
            summary_lines.append(
                "  optimized vs baseline complex: "
                f"S11 max={complex_metrics['baseline_max_abs_complex_s11']:.3e}, "
                f"relL2={complex_metrics['baseline_relative_l2_complex_s11']:.3e}; "
                f"S21 max={complex_metrics['baseline_max_abs_complex_s21']:.3e}, "
                f"relL2={complex_metrics['baseline_relative_l2_complex_s21']:.3e}"
            )
        summary_lines.append("")
        summary_records.append(
            {
                "label": label,
                "points": int(len(grid)),
                "freq_min_GHz": float(f_lo),
                "freq_max_GHz": float(f_hi),
                "s21_peak_GHz": float(pf_s),
                "s21_3db_lo_GHz": float(lo_s),
                "s21_3db_hi_GHz": float(hi_s),
                "s21_3db_center_shift_GHz": float((lo_s + hi_s) / 2 - (lo + hi) / 2),
                "max_abs_db_s11": float(np.max(np.abs(s11_db - s11_ref_db))),
                "mean_abs_db_s11": float(np.mean(np.abs(s11_db - s11_ref_db))),
                "rms_db_s11": float(np.sqrt(np.mean((s11_db - s11_ref_db) ** 2))),
                "max_abs_db_s21": float(np.max(np.abs(s21_db - s21_ref_db))),
                "mean_abs_db_s21": float(np.mean(np.abs(s21_db - s21_ref_db))),
                "rms_db_s21": float(np.sqrt(np.mean((s21_db - s21_ref_db) ** 2))),
                "max_abs_lin_s11": float(np.max(np.abs(s11 - s11_ref))),
                "mean_abs_lin_s11": float(np.mean(np.abs(s11 - s11_ref))),
                "rms_lin_s11": float(np.sqrt(np.mean((s11 - s11_ref) ** 2))),
                "max_abs_lin_s21": float(np.max(np.abs(s21 - s21_ref))),
                "mean_abs_lin_s21": float(np.mean(np.abs(s21 - s21_ref))),
                "rms_lin_s21": float(np.sqrt(np.mean((s21 - s21_ref) ** 2))),
                **complex_metrics,
            }
        )

    axes[0].set_ylabel("|S11| [dB]")
    axes[0].grid(True, alpha=0.4)
    axes[0].legend(loc="best", fontsize=9)
    axes[0].set_title("BP filter: bp_fem_solver vs HFSS TFE reference")
    axes[1].set_ylabel("|S21| [dB]")
    axes[1].set_xlabel("Frequency [GHz]")
    axes[1].grid(True, alpha=0.4)
    axes[1].legend(loc="best", fontsize=9)
    fig.tight_layout()
    fig.savefig(out_dir / "tfe_compare.png", dpi=160)
    plt.close(fig)

    # Plot 2: linear magnitude residual (in 0..1 units), more honest near nulls.
    fig2, ax = plt.subplots(figsize=(10, 5))
    for idx, (label, df) in enumerate(solver_runs):
        s11 = _interp(grid, df, "S11_mag")
        s21 = _interp(grid, df, "S21_mag")
        color = cycle[idx % len(cycle)]
        ax.plot(grid, s11 - s11_ref, label=f"{label} ΔS11", linewidth=1.0, color=color)
        ax.plot(grid, s21 - s21_ref, label=f"{label} ΔS21", linewidth=1.0, linestyle=":", color=color)
    ax.axhline(0.0, color="black", linewidth=0.6)
    ax.set_xlabel("Frequency [GHz]")
    ax.set_ylabel("solver - HFSS  (linear |S|)")
    ax.grid(True, alpha=0.4)
    ax.legend(loc="best", fontsize=8)
    ax.set_title("BP filter linear-magnitude residual vs HFSS TFE reference")
    fig2.tight_layout()
    fig2.savefig(out_dir / "tfe_compare_residual.png", dpi=160)
    plt.close(fig2)

    out_table.to_csv(out_dir / "tfe_compare.csv", index=False, float_format="%.10g")
    pd.DataFrame(summary_records).to_csv(
        out_dir / "benchmark_summary.csv", index=False, float_format="%.10g"
    )
    summary = "\n".join(summary_lines)
    (out_dir / "tfe_compare.txt").write_text(summary, encoding="utf-8")
    metadata = {
        "hfss_reference": str(args.hfss_csv.resolve()),
        "hfss_csv_sha256": hashlib.sha256(args.hfss_csv.read_bytes()).hexdigest(),
        "grid": "exact-points" if points_match else "interpolated-overlap",
        "solver_runs": source_records,
    }
    (out_dir / "comparison_inputs.json").write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
    if os.path.lexists(preferred):
        target = preferred.with_name(f"{preferred.name}-{uuid.uuid4().hex}")
        out_dir.rename(target)
    else:
        target = preferred
        try:
            out_dir.rename(target)
        except OSError:
            if not os.path.lexists(target):
                raise
            target = preferred.with_name(f"{preferred.name}-{uuid.uuid4().hex}")
            out_dir.rename(target)
    print(summary)
    print(f"Output directory: {target}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
