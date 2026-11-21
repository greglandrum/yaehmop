#!/usr/bin/env python3
"""Compare the YAeHMOP bands with the WanTiBEXOS bands of the converted H(R), S(R).

    python3 compare_bands.py

For each case (yaehmop_soc/wtb_soc and yaehmop_nosoc/wtb_nosoc) it reads the
YAeHMOP .band file and the wtb.x out/bands.dat. Both are taken relative to their
valence band maximum (wtb.x writes bands.dat that way). The two programs use the same
k points (40 intervals per line in YAeHMOP = 41 points per segment in
kpoints.dat), so the comparison is point by point.

It also diagonalizes, in double precision (numpy/scipy), the H(k), S(k) rebuilt
from the converted file (huckel-*.dat) on the same k points. This checks the
conversion apart from wtb.x, which works in single precision. Writes validation.txt
and bands_yaehmop_vs_wtb.png.
"""
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
CASES = [("nosoc", "without SOC", 18 // 2), ("soc", "with SOC", 18)]   # MoS2: 18 electrons


def read_yaehmop_band(path):
    lines = Path(path).read_text().splitlines()
    nst = int(next(l for l in lines if "orbitals in the unit cell" in l).split()[0])
    nsp = int(next(l for l in lines if "Special points were done with" in l).split()[0])
    j = next(n for n, l in enumerate(lines) if "orbitals in the unit cell" in l)
    labels = [l.split()[0] for l in lines[j + 1:j + 1 + nsp]]
    ev, cur = [], []
    for l in lines[next(n for n, l in enumerate(lines) if "Begin band data" in l) + 1:]:
        if l.startswith("#END"):
            break
        if l.startswith("; K point"):
            if cur:
                ev.append(cur)
            cur = []
        elif l.strip():
            cur.append(float(l))
    ev.append(cur)
    return np.array(ev), labels, nst


def read_wtb_bands(path):
    """bands.dat: '#band i' blocks of 'k-path energy sx sy sz orbital'."""
    bands, cur = [], None
    for l in Path(path).read_text().splitlines():
        s = l.strip()
        if s.startswith("#band"):
            if cur:
                bands.append(cur)
            cur = []
        elif s and not s.startswith("#") and cur is not None:
            p = s.split()
            cur.append((float(p[0]), float(p[1])))
    if cur:
        bands.append(cur)
    a = np.array(bands)                     # (nbands, nk, 2)
    return a[0, :, 0], a[:, :, 1].T         # kdist, energies (nk, nbands)


def read_tb(path):
    """WanTiBEXOS TB file -> lattice vectors, size, [(R [A], H(R), S(R))]."""
    L = Path(path).read_text().splitlines()
    cell = np.array([[float(x) for x in L[i].split()] for i in (3, 4, 5)])
    n, nv = int(L[6]), int(L[7])
    d = np.loadtxt(path, skiprows=9).reshape(nv, n * n, 8)
    return cell, n, [(b[0, :3], (b[:, 5] + 1j * b[:, 6]).reshape(n, n), b[:, 7].reshape(n, n)) for b in d]


def exact_bands(tbfile, kfrac):
    """Double-precision eigenvalues of the converted H(R), S(R) at fractional k."""
    import scipy.linalg as sl
    cell, n, cells = read_tb(tbfile)
    B = 2 * np.pi * np.linalg.inv(cell).T
    out = []
    for k in kfrac:
        kc = k @ B
        H = sum(np.exp(1j * (kc @ R)) * h for R, h, _ in cells)
        S = sum(np.exp(1j * (kc @ R)) * s for R, _, s in cells)
        out.append(sl.eigh(H, S, eigvals_only=True))
    return np.array(out)


def main():
    report, fig_data = [], []
    for tag, title, nocc in CASES:
        ey, labels, nst = read_yaehmop_band(HERE / f"yaehmop_{tag}" / "MoS2.bind.band")
        kd, ew = read_wtb_bands(HERE / f"wtb_{tag}" / "out" / "bands.dat")
        nseg = len(labels) - 1
        npts = len(kd) // nseg
        ppl = (len(ey) - 1) // nseg
        assert npts == ppl + 1, "the k paths do not match"
        # wtb point i of segment s  <->  YAeHMOP point s*ppl + i
        idx = np.array([s * ppl + i for s in range(nseg) for i in range(npts)])
        ey_rel = ey - ey[:, :nocc].max()
        ey_on_w = ey_rel[idx]
        assert ew.shape[1] == nst
        diff = np.abs(ew - ey_on_w)
        win = (ey_on_w > -6) & (ey_on_w < 4)
        gy = ey_rel[:, nocc:].min()
        gw = ew[:, nocc:].min() - ew[:, :nocc].max()
        # double-precision bands of the converted file on the same k points
        lines = (HERE / f"wtb_{tag}" / "kpoints.dat").read_text().split()
        kp = [np.array([float(x) for x in lines[2 + 4 * j:5 + 4 * j]]) for j in range(2 * nseg)]
        kfrac = [kp[2 * s] + (kp[2 * s + 1] - kp[2 * s]) * i / (npts - 1)
                 for s in range(nseg) for i in range(npts)]
        tbfile = next((HERE / f"wtb_{tag}").glob("huckel-*.dat"))
        ex = exact_bands(tbfile, kfrac)
        ex -= ex[:, :nocc].max()
        dex = np.abs(ex - ey_on_w)
        report += [f"MoS2 {title}: {nst} bands, {len(kd)} k points (WanTiBEXOS) on the YAeHMOP path",
                   f"  converted file ({tbfile.name}) diagonalized in double precision vs YAeHMOP:",
                   f"      all bands {dex.max():.2e} eV;  bands within VBM-6 .. VBM+4 eV {dex[win].max():.2e} eV",
                   f"  WanTiBEXOS (wtb.x, single precision) vs YAeHMOP:",
                   f"      all bands {diff.max():.2e} eV;  bands within VBM-6 .. VBM+4 eV {diff[win].max():.2e} eV",
                   f"  band gap: YAeHMOP {gy:.6f} eV, WanTiBEXOS {gw:.6f} eV"]
        if "SOC" in title and "without" not in title:
            iK = labels.index("K") * ppl
            iKw = int(np.where(idx == iK)[0][0])
            sy = ey_rel[iK, nocc - 1] - ey_rel[iK, nocc - 2]
            sw = ew[iKw, nocc - 1] - ew[iKw, nocc - 2]
            report.append(f"  SOC splitting of the top valence band at K: YAeHMOP {sy:.6f} eV, "
                          f"WanTiBEXOS {sw:.6f} eV")
        report.append("")
        fig_data.append((title, kd, ew, ey_on_w, labels, nseg, npts))
    (HERE / "validation.txt").write_text("\n".join(report))
    print("\n".join(report))

    fig, axs = plt.subplots(1, 2, figsize=(11, 5), sharey=True)
    for ax, (title, kd, ew, ey, labels, nseg, npts) in zip(axs, fig_data):
        for n in range(ew.shape[1]):
            ax.plot(kd, ey[:, n], color="0.2", lw=1.3, label="YAeHMOP" if n == 0 else None)
            ax.plot(kd[::2], ew[::2, n], "o", ms=2.6, mfc="none", color="#c8372d",
                    label="WanTiBEXOS (converted H, S)" if n == 0 else None)
        ticks = [kd[s * npts] for s in range(nseg)] + [kd[-1]]
        for t in ticks:
            ax.axvline(t, color="0.8", lw=0.6, zorder=0)
        ax.axhline(0, color="0.5", lw=0.6, ls=":")
        ax.set_xticks(ticks)
        ax.set_xticklabels(["Γ" if l.upper() in ("G", "GAMMA") else l for l in labels])
        ax.set_xlim(kd[0], kd[-1]); ax.set_ylim(-6, 4)
        ax.set_title(f"MoS$_2$ {title}", fontsize=11)
    axs[0].set_ylabel("E − E$_{VBM}$ (eV)")
    h, l = axs[1].get_legend_handles_labels()
    fig.legend(h, l, loc="upper right", fontsize=8, frameon=False, ncol=2)
    fig.suptitle("YAeHMOP → yaehmop2wtb.py → WanTiBEXOS (wtb.x)", fontsize=11, x=0.35)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(HERE / "bands_yaehmop_vs_wtb.png", dpi=200)
    print("wrote", HERE / "bands_yaehmop_vs_wtb.png")


if __name__ == "__main__":
    main()
