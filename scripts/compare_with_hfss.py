"""Compare bp_fem_solver S-parameter outputs with the HFSS TFE reference.

The HFSS reference file `S Parameter Plot 1.csv` carries:
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

    python compare_with_hfss.py <hfss_csv> <out_dir> --batch-root result
        # scans immediate child directories containing s_parameters.csv
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


def _load_solver(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path).copy()
    df["freq_GHz"] = df["freq_Hz"].astype(float) / 1e9
    df["S11_mag"] = np.hypot(df["S11_real"], df["S11_imag"])
    df["S21_mag"] = np.hypot(df["S21_real"], df["S21_imag"])
    return df.sort_values("freq_GHz").reset_index(drop=True)


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
        "--batch-root",
        type=Path,
        help="Scan immediate child directories for s_parameters.csv and compare all runs.",
    )
    args = parser.parse_args(argv[1:])

    out_dir = args.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    hfss = _load_hfss(args.hfss_csv)
    solver_runs = []
    solver_specs = list(args.solver)
    if args.batch_root is not None:
        for child in sorted(args.batch_root.iterdir()):
            csv_path = child / "s_parameters.csv"
            if child.is_dir() and csv_path.exists():
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
        solver_runs.append((label, _load_solver(Path(path))))

    f_lo = max(hfss["freq_GHz"].min(), max(s["freq_GHz"].min() for _, s in solver_runs))
    f_hi = min(hfss["freq_GHz"].max(), min(s["freq_GHz"].max() for _, s in solver_runs))
    grid = hfss[(hfss["freq_GHz"] >= f_lo) & (hfss["freq_GHz"] <= f_hi)]["freq_GHz"].to_numpy()

    s11_ref = _interp(grid, hfss, "S11_mag")
    s21_ref = _interp(grid, hfss, "S21_mag")
    s11_ref_db = _to_db(s11_ref)
    s21_ref_db = _to_db(s21_ref)

    summary_lines = [
        f"HFSS reference: {args.hfss_csv}",
        f"Compared band: {f_lo:.4f} - {f_hi:.4f} GHz, {len(grid)} points",
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
    print(summary)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
