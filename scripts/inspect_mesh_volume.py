"""Sanity-check the bp_filter NGMesh: total volume of Box1, vs the analytic
volume the AEDT design parameters describe.

If the iris boolean subtractions are missing from the mesh, Box1 mesh volume
will equal the bare a*b*L rectangular box. If the irises ARE in the mesh, the
volume will be smaller (waveguide minus iris solids).
"""
from __future__ import annotations
import re
from pathlib import Path

import numpy as np


def main() -> None:
    p = Path(r"C:\Users\82061\Desktop\TFE-Direct-2\example\bp_filter\current.ngmesh")
    text = p.read_text(encoding="utf-8", errors="ignore")

    # Coordinates are stored in METERS in the file body (the file declares
    # user_unit_name=mm and user_units_per_one_meter=1000, but the actual
    # numeric coords are SI-meters - we'll keep meters and multiply at the end).
    pts = {}
    for m in re.finditer(r"^pid\s+(\d+)\s+coords\s+(\S+)\s+(\S+)\s+(\S+)", text, re.MULTILINE):
        pid = int(m.group(1))
        pts[pid] = np.array([float(m.group(2)), float(m.group(3)), float(m.group(4))])

    # Tetrahedra: "veid <id> ... body_id <bid> ... vert_ids <a> <b> <c> <d>"
    tet_pat = re.compile(
        r"^veid\s+\d+.*?body_id\s+(\d+).*?vert_ids\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)",
        re.MULTILINE,
    )
    vols_m3 = {}
    counts = {}
    for m in tet_pat.finditer(text):
        bid = int(m.group(1))
        v = [pts[int(m.group(2 + i))] for i in range(4)]
        a = v[1] - v[0]
        b = v[2] - v[0]
        c = v[3] - v[0]
        vol = abs(np.dot(a, np.cross(b, c))) / 6.0  # m^3
        vols_m3[bid] = vols_m3.get(bid, 0.0) + vol
        counts[bid] = counts.get(bid, 0) + 1

    print(f"Total tets        : {sum(counts.values())}")
    print(f"Per-body breakdown (volumes in mm^3):")
    for bid in sorted(vols_m3):
        print(
            f"  body_id={bid:3d}  tets={counts[bid]:6d}"
            f"  volume={vols_m3[bid] * 1e9:.4f} mm^3"
        )

    # Box1 nominal: a=7.112, b=3.556, L=2*l1+2*l2+2*l3+5*t = 31.2174 mm
    a, b = 7.112, 3.556
    L = 2 * 6 + 2 * 3.3327 + 2 * 3.776 + 5 * 1
    nominal_box1 = a * b * L  # mm^3
    print()
    print(f"Box1 nominal hollow rectangular volume  : a*b*L = {a}*{b}*{L:.4f}")
    print(f"                                       = {nominal_box1:.4f} mm^3")
    if 6 in vols_m3:
        actual_mm3 = vols_m3[6] * 1e9
        ratio = actual_mm3 / nominal_box1
        print()
        print(f"Box1 mesh volume = {actual_mm3:.4f} mm^3")
        print(f"Box1 mesh volume / nominal hollow volume = {ratio:.6f}")
        if ratio < 0.999:
            print(f"  => irises DISPLACE  {(1.0-ratio)*100:.2f} %  of the bare-WG volume")
        else:
            print("  => mesh appears to be empty rectangular waveguide; no irises")


if __name__ == "__main__":
    main()
