"""Archive a completed schema-2 BP-FEM run in an existing benchmark case."""

import argparse
import csv
import hashlib
import json
import math
import os
import shutil
from pathlib import Path
from uuid import uuid4


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as source:
        for block in iter(lambda: source.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def is_link(path):
    info = path.lstat()
    return path.is_symlink() or bool(getattr(info, "st_file_attributes", 0) & 0x400)


def inventory(root):
    paths = []
    for directory, subdirs, files in os.walk(root):
        for name in (*subdirs, *files):
            path = Path(directory) / name
            if is_link(path):
                raise ValueError(f"Symbolic link or reparse point in run: {path}")
            if path.is_file():
                paths.append(path.relative_to(root))
    return {p: sha256(root / p) for p in sorted(paths)}


def archive(run, case, label, aedt, mesh):
    if not label or label in (".", "..") or Path(label).name != label or label.startswith("."):
        raise ValueError("Archive label must be a single visible directory name")
    for folder in (run, case, case / "project", case / "results"):
        if not folder.is_dir() or is_link(folder):
            raise ValueError(f"Expected a regular directory: {folder}")
    if run.name.startswith(".bpfem-pending-") or os.path.lexists(run / "run.json.pending"):
        raise ValueError("The run has not been published")
    report = json.loads((run / "run.json").read_text(encoding="utf-8"))
    if report.get("status") != "ok" or report.get("schema", 0) < 2:
        raise ValueError("Expected a successful run with recorded input SHA-256")
    algorithm = report.get("sweep_strategy")
    required_outputs = ["run.log", "s_parameters.csv", "timing.json", "diagnostics.json"]
    if algorithm == "wcawe":
        required_outputs.append("basis_condition.csv")
    for filename in required_outputs:
        output = run / filename
        if not output.is_file() or output.stat().st_size == 0:
            raise ValueError(f"Required output is missing or empty: {filename}")
    timing = json.loads((run / "timing.json").read_text(encoding="utf-8"))
    diagnostics = json.loads((run / "diagnostics.json").read_text(encoding="utf-8"))
    if (not algorithm or timing.get("status") != "ok" or diagnostics.get("algorithm") != algorithm
            or diagnostics.get("reduced_solve_succeeded") is not True):
        raise ValueError("Timing or solver diagnostics disagree with the successful run")
    required = ("freq_Hz", "S11_real", "S11_imag", "S21_real", "S21_imag")
    with (run / "s_parameters.csv").open(encoding="utf-8", newline="") as stream:
        reader = csv.DictReader(stream)
        headers = reader.fieldnames or []
        if len(headers) != len(set(headers)) or not set(required).issubset(headers):
            raise ValueError("S-parameter CSV is missing required complex columns or has duplicate headers")
        checked = (*required, *(key for key in ("S11_dB", "S21_dB") if key in headers))
        frequencies = []
        for number, row in enumerate(reader, 2):
            try:
                values = [float(row[key]) for key in checked]
            except (TypeError, ValueError) as error:
                raise ValueError(f"Incomplete or invalid S-parameter values at CSV line {number}") from error
            if None in row or not all(math.isfinite(value) for value in values):
                raise ValueError(f"Invalid or non-finite S-parameter values at CSV line {number}")
            frequency = values[0]
            if frequency <= 0 or (frequencies and frequency <= frequencies[-1]):
                raise ValueError("S-parameter frequencies must be positive and strictly increasing")
            frequencies.append(frequency)
    if not frequencies or frequencies != [point["freq_hz"] for point in report["sweep"]]:
        raise ValueError("S-parameter frequency samples differ from the successful report")
    for name, input_file in (("aedt", aedt), ("mesh", mesh)):
        if not input_file.is_file() or is_link(input_file) or report.get("inputs", {}).get(f"{name}_sha256") != sha256(input_file):
            raise ValueError(f"Archived project {name} differs from the run input: {input_file}")
    dest = case / "results" / label
    if os.path.lexists(dest):
        raise FileExistsError(f"Archive already exists: {dest}")
    manifest = inventory(run)
    if Path("SHA256SUMS.txt") in manifest:
        raise ValueError("Run already contains an archive manifest")
    # Keep the rename on the same volume and inherit the benchmark directory's access policy.
    staged = dest.parent / f".bpfem-archive-pending-{uuid4().hex}"
    staged.mkdir()
    for relative in manifest:
        target = staged / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(run / relative, target)
    if inventory(staged) != manifest or inventory(run) != manifest:
        raise ValueError(f"Input changed during archival. Inspect {staged}")
    checksum = "\n".join(f"{digest}  {relative.as_posix()}" for relative, digest in manifest.items()) + "\n"
    (staged / "SHA256SUMS.txt").write_text(checksum, encoding="ascii", newline="\n")
    staged.rename(dest)
    return dest


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run", type=Path, help="Published run directory reported by the solver")
    parser.add_argument("case", type=Path, help="Existing benchmark case with project/ and results/")
    parser.add_argument("--label", required=True, help="New directory name under results/")
    parser.add_argument("--aedt", required=True, type=Path, help="AEDT filename relative to the case project/")
    parser.add_argument("--mesh", default=Path("current.ngmesh"), type=Path, help="Mesh filename relative to the case project/")
    args = parser.parse_args()
    chosen = (args.aedt, args.mesh)
    if any(path.name != str(path) for path in chosen):
        parser.error("AEDT and mesh must be direct children of the benchmark project directory")
    project = (args.case / "project").absolute()
    paths = [project / path for path in chosen]
    print(archive(args.run.absolute(), args.case.absolute(), args.label, *paths))
