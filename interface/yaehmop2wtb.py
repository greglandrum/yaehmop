#!/usr/bin/env python3
"""YAeHMOP H and S  ->  WanTiBEXOS non-orthogonal tight-binding files.

    python3 yaehmop2wtb.py RUN.bind [-o OUTDIR] [--zero auto|fermi|midgap|<eV>]

RUN.bind is the input file of a YAeHMOP run made with 'Dump Hamil' and
'Dump Overlap' and a 'K points' block holding a full Gamma-centred mesh. The
outputs of that run are found next to it and the kind of run is detected:

  spinless            RUN.bind.HAM/.OV/.out                    -> huckel-NP.dat,  basis_set-NP
  SOC (--soc)         RUN.bind.HAM/.OV (2n x 2n), .out         -> huckel-SOC.dat, basis_set-SOC
  spin polarized      RUN_up.bind.*, RUN_down.bind.*           -> huckel-SP.dat,  basis_set-SP

What it does
  1. reads the structure, the lattice, the overlap radius r and the k points
     from RUN.bind, and the orbitals of each element, the electron count and the
     Fermi level from the .out;
  2. reads H(k), S(k) from the binary dumps (int nk, int n, then per k an n x n
     double array; element (i,j), i<j, is m[i][j] - i*m[j][i]);
  3. gets the real-space matrices by an inverse Fourier transform on the mesh,
         H(R) = (1/N) sum_k exp(-2 pi i k.R) H(k),   |R_i| <= r_i,
     which is exact when the mesh has N_i >= 2 r_i + 1 points along each
     periodic direction (YAeHMOP sums the lattice up to r cells), and checks it
     by rebuilding H(k), S(k) on the mesh;
  4. writes the WanTiBEXOS file (same layout as siesta2wtb.py, read with
     DFT= S): type, scissors, Fermi level, 3 lattice vectors [A], basis size,
     number of R vectors, then 'Rx Ry Rz [A]  i  j  ReH  ImH  S';
  5. writes the basis file (index species x y z l m spin), an input file and
     a k-path file for wtb.x (the band path of RUN.bind).

Energy zero: wtb.x diagonalizes (H(k) - E_F*I, S(k)), which is not a rigid shift
when S != I. The file therefore holds H(R) - E0*S(R) with Fermi level 0, so wtb.x
gives exactly the YAeHMOP eigenvalues minus E0. --zero auto uses the midgap of
the mesh eigenvalues for a gapped system and the Fermi level of the .out
otherwise. Spin polarized: with a common Fermi level (bind --common-fermi, read
from RUN.bind.spin) both channels get that one zero; in a fixed-moment run each
channel gets its own Fermi level, so that both channels have E_F at 0, as when
the channels were fitted to a DFT run with that moment. SOC with two spin
tables (one 2n x 2n run) is converted like any SOC run.

Needs numpy (and scipy for the eigenvalues).
"""
import argparse, re, sys
from pathlib import Path
import numpy as np

ORBITALS = {"s": ["s"], "p": ["px", "py", "pz"],
            "d": ["dx2-y2", "dz2", "dxy", "dxz", "dyz"],
            "f": ["fz3", "fxz2", "fyz2", "fxyz", "fz(x2-y2)", "fx(x2-3y2)", "fy(3x2-y2)"]}
# (l, m) of the real orbitals, as in the basis files of WanTiBEXOS/sisl
ORBITAL_LM = {"s": (0, 0), "px": (1, 1), "py": (1, -1), "pz": (1, 0),
              "dx2-y2": (2, 2), "dz2": (2, 0), "dxy": (2, -2), "dxz": (2, 1), "dyz": (2, -1),
              "fz3": (3, 0), "fxz2": (3, 1), "fyz2": (3, -1), "fxyz": (3, -2),
              "fz(x2-y2)": (3, 2), "fx(x2-3y2)": (3, 3), "fy(3x2-y2)": (3, -3)}


# ------------------------------------------------------------------ input file
def read_bind_input(path):
    """Structure, lattice, overlap radius, k points, band path of a .bind file."""
    raw = Path(path).read_text().splitlines()
    lines = [l.split(";")[0].rstrip() for l in raw]
    lines = [l for l in lines if l.strip()]
    up = [l.strip().upper() for l in lines]
    info = dict(kpoints=None, band=None, electrons=None, charge=0.0, crystal=None)
    i = 0
    while i < len(lines):
        u = up[i]
        if u.startswith("GEOMETRY"):
            info["crystallographic"] = "CRYSTAL" in u
            n = int(lines[i + 1].split()[0])
            atoms = []
            for l in lines[i + 2:i + 2 + n]:
                p = l.split()
                atoms.append((p[1], np.array([float(x) for x in p[2:5]])))
            info["atoms"] = atoms
            i += 2 + n
            continue
        if u == "LATTICE":
            dim = int(lines[i + 1].split()[0])
            ov = [int(x) for x in lines[i + 2].split()[:3]]
            vec = [tuple(int(x) for x in lines[i + 3 + j].split()[:2]) for j in range(dim)]
            info.update(dim=dim, overlaps=ov, lattice_atoms=vec)
            i += 3 + dim
            continue
        if u.startswith("CRYSTAL SPEC"):
            info["crystal"] = ([float(x) for x in lines[i + 1].split()[:3]],
                               [float(x) for x in lines[i + 2].split()[:3]])
            i += 3
            continue
        if u.startswith("K POINTS") or u == "KPOINTS":
            n = int(lines[i + 1].split()[0])
            info["kpoints"] = np.array([[float(x) for x in l.split()[:3]] for l in lines[i + 2:i + 2 + n]])
            i += 2 + n
            continue
        if u == "BAND" or u == "BANDS":
            ppl = int(lines[i + 1].split()[0]); nsp = int(lines[i + 2].split()[0])
            pts = []
            for l in lines[i + 3:i + 3 + nsp]:
                p = l.split()
                pts.append((p[0], np.array([float(x) for x in p[1:4]])))
            info["band"] = (ppl, pts)
            i += 3 + nsp
            continue
        if u.startswith("ELECTRONS"):
            p = lines[i].split()
            info["electrons"] = float(p[1]) if len(p) > 1 else float(lines[i + 1].split()[0])
            i += 1 if len(p) > 1 else 2
            continue
        if u.startswith("CHARGE"):
            p = lines[i].split()
            info["charge"] = float(p[1]) if len(p) > 1 else float(lines[i + 1].split()[0])
            i += 1 if len(p) > 1 else 2
            continue
        i += 1
    for key in ("atoms", "dim"):
        if key not in info:
            raise SystemExit(f"{path}: no {'Geometry' if key == 'atoms' else 'Lattice'} section found")
    return info


def cell_from_spec(lengths, angles):
    """Lattice vectors as YAeHMOP builds them: a along x, b in the xy plane."""
    a, b, c = lengths
    al, be, ga = np.radians(angles)
    A = np.zeros((3, 3))
    A[0] = [a, 0, 0]
    A[1] = [b * np.cos(ga), b * np.sin(ga), 0]
    cx = c * np.cos(be)
    cy = c * (np.cos(al) - np.cos(be) * np.cos(ga)) / np.sin(ga)
    A[2] = [cx, cy, np.sqrt(max(c * c - cx * cx - cy * cy, 0.0))]
    return A


def structure(info):
    """Real atoms (symbol, cartesian position) and the 3 lattice vectors [A]."""
    atoms, dim = info["atoms"], info["dim"]
    nreal = len(atoms) - dim            # the ends of the lattice vectors are the last atoms
    if info["crystallographic"]:
        if info["crystal"] is None:
            raise SystemExit("Crystallographic geometry without a 'Crystal Spec' section")
        cell = cell_from_spec(*info["crystal"])
        cart = [(s, f @ cell) for s, f in atoms]
    else:
        cart = list(atoms)
        cell = np.zeros((3, 3))
        for j, (b, e) in enumerate(info["lattice_atoms"]):
            cell[j] = cart[e - 1][1] - cart[b - 1][1]
        pos = np.array([p for _, p in cart[:nreal]])
        # vacuum directions of 1D/2D systems: unit vectors orthogonal to the
        # periodic ones, long enough to hold the atoms plus 15 A of vacuum
        for j in range(dim, 3):
            if j == 1:                      # 1D: any direction orthogonal to a
                t = np.array([0.0, 0.0, 1.0])
                if abs(cell[0] @ t) > 0.9 * np.linalg.norm(cell[0]):
                    t = np.array([1.0, 0.0, 0.0])
                v = np.cross(cell[0], t)
            else:                           # third vector: orthogonal to the first two
                v = np.cross(cell[0], cell[1])
            v = v / np.linalg.norm(v)
            ext = np.ptp(pos @ v) if len(pos) > 1 else 0.0
            cell[j] = v * max(20.0, ext + 15.0)
    real = [(s, p) for s, p in cart[:nreal] if not s.startswith("&")]
    if any(s.startswith("*") for s, _ in real):
        raise SystemExit("Atoms with parameters given in the input ('*') are not supported")
    return real, cell


# ------------------------------------------------------------------ .out file
def read_out(path):
    text = Path(path).read_text()
    shells, cur = {}, None
    for l in text.splitlines():
        m = re.match(r"ATOM:\s+(\S+)\s+Atomic number:.*# Valence Electrons:\s*(\d+)", l)
        if m:
            cur = m.group(1).upper(); shells[cur] = dict(shells=[], valence=int(m.group(2)))
            continue
        m = re.match(r"\s+\d+([SPDF]):", l)
        if cur and m:
            shells[cur]["shells"].append(m.group(1).lower())
            continue
        if cur and not l.startswith("\t"):
            cur = None
    def key(k, cast=float):
        v = re.findall(rf"#{k}:\s*(\S+)", text)
        return cast(v[-1]) if v else None
    return dict(shells=shells, fermi=key("Fermi_Energy"), soc_states=key("SOC_STATES", int),
                channel=(re.findall(r"#SPIN_CHANNEL:\s*(\S+)", text) or [None])[-1],
                # with --common-fermi the channel's electrons come from the shared Fermi level
                channel_electrons=key("ELECTRONS_THIS_CHANNEL_COMMON_FERMI")
                or key("ELECTRONS_THIS_CHANNEL"),
                lattice=[np.array([float(x) for x in l.split()[1:4]])
                         for l in re.findall(r"^\([abc]\)\s+.*$", text, re.M)])


def orbital_labels(real_atoms, shells):
    labels = []
    for ia, (sym, pos) in enumerate(real_atoms):
        el = sym.upper()
        if el not in shells:
            raise SystemExit(f"element {sym} not found in the parameters of the .out")
        for sh in "spdf":
            if sh in shells[el]["shells"]:
                labels += [(ia, sym, orb) for orb in ORBITALS[sh]]
    return labels


# ------------------------------------------------------------------ matrices
def read_dump(path):
    raw = Path(path).read_bytes()
    nk, n = (int(x) for x in np.frombuffer(raw[:8], dtype=np.int32))
    data = np.frombuffer(raw[8:8 + 8 * nk * n * n], dtype=np.float64)
    if data.size != nk * n * n:
        raise SystemExit(f"{path}: file is shorter than its header says ({nk} k points, n = {n})")
    data = data.reshape(nk, n, n)
    iu = np.triu_indices(n, 1)
    out = np.zeros((nk, n, n), complex)
    for q in range(nk):
        m = data[q]
        H = np.diag(np.diag(m)).astype(complex)
        H[iu] = m[iu] - 1j * m.T[iu]
        H[(iu[1], iu[0])] = np.conj(H[iu])
        out[q] = H
    return out


def check_mesh(kpts, dim, overlaps):
    """The k points must be a full Gamma-centred N1 x N2 x N3 mesh; returns N."""
    N = []
    for d in range(3):
        if d >= dim:
            if np.abs(kpts[:, d]).max() > 1e-6:
                raise SystemExit("k points have a component along a non-periodic direction")
            N.append(1); continue
        frac = np.mod(kpts[:, d] + 1e-9, 1.0)
        vals = np.unique(np.round(frac, 6))
        n = len(vals)
        if not np.allclose(vals, np.arange(n) / n, atol=2e-6):
            raise SystemExit(f"the k points are not a uniform Gamma-centred mesh along direction {d+1}")
        N.append(n)
    if len(kpts) != int(np.prod(N)):
        raise SystemExit(f"{len(kpts)} k points, but a {N[0]}x{N[1]}x{N[2]} mesh has {int(np.prod(N))}")
    alias = [d for d in range(dim) if N[d] < 2 * overlaps[d] + 1]
    # rounding of the k points in the input: YAeHMOP uses them as written, so the
    # Fourier transform is exact only if they are close to j/N
    dev = max(np.abs(kpts[:, d] - np.round(kpts[:, d] * N[d]) / N[d]).max() for d in range(dim))
    return N, alias, dev


def real_space(Hk, Sk, kpts, dim, overlaps):
    import itertools
    rng = [range(-overlaps[d], overlaps[d] + 1) if d < dim else (0,) for d in range(3)]
    cells = []
    for R in itertools.product(*rng):
        R = np.array(R)
        ph = np.exp(-2j * np.pi * (kpts @ R))[:, None, None]
        cells.append((R, (ph * Hk).sum(0) / len(kpts), (ph * Sk).sum(0) / len(kpts)))
    return cells


def rebuild_error(cells, Hk, Sk, kpts):
    eH = eS = 0.0
    for q, k in enumerate(kpts):
        H = sum(np.exp(2j * np.pi * (k @ R)) * h for R, h, _ in cells)
        S = sum(np.exp(2j * np.pi * (k @ R)) * s for R, _, s in cells)
        eH = max(eH, np.abs(H - Hk[q]).max()); eS = max(eS, np.abs(S - Sk[q]).max())
    return eH, eS


def mesh_gap(Hk, Sk, nocc):
    import scipy.linalg as sl
    if nocc is None or abs(nocc - round(nocc)) > 1e-6 or round(nocc) < 1:
        return None
    nocc = int(round(nocc))
    e = np.array([sl.eigh(h, s, eigvals_only=True) for h, s in zip(Hk, Sk)])
    if nocc >= e.shape[1]:
        return None
    vbm, cbm = e[:, nocc - 1].max(), e[:, nocc].min()
    return vbm, cbm


# ------------------------------------------------------------------ main
def channel_paths(run):
    run = Path(run)
    stem, ext = (run.stem, run.suffix) if run.suffix else (run.name, "")
    def base(tag):
        return run.with_name(f"{stem}{tag}{ext}")
    return base


def energy_zero(choice, out, Hk, Sk, nocc, label):
    if choice not in ("auto", "fermi", "midgap"):
        return float(choice), "value given with --zero"
    gap = mesh_gap(Hk, Sk, nocc) if choice in ("auto", "midgap") else None
    if gap and (choice == "midgap" or gap[1] - gap[0] > 1e-3):
        return 0.5 * (gap[0] + gap[1]), f"midgap of the mesh eigenvalues (VBM {gap[0]:.6f}, CBM {gap[1]:.6f} eV)"
    if choice == "midgap":
        raise SystemExit(f"{label}: no gap found on the mesh; use --zero fermi or a value")
    if out["fermi"] is None:
        raise SystemExit(f"{label}: no gap and no #Fermi_Energy in the .out "
                         "(add 'Average Properties' to the input or give --zero <eV>)")
    return out["fermi"], "Fermi level of the .out"


def common_fermi(run):
    """common Fermi level of a spin-polarized run (bind --common-fermi), from <run>.spin"""
    p = Path(str(run) + ".spin")
    if not p.exists():
        return None
    v = re.findall(r"^#Fermi_Energy:\s*(\S+)", p.read_text(), re.M)
    return float(v[-1]) if v else None


def write_tb(path, kind, cell, n2, blocks):
    with open(path, "w") as f:
        f.write(f"{kind}\n0.000000\n0.000000\n")
        for v in cell:
            f.write(f"{v[0]:.8f}  {v[1]:.8f}  {v[2]:.8f}\n")
        f.write(f"{n2}\n{len(blocks)}\n#rcell x   rcell y   rcell z   i   j   ReH   ImH   S\n")
        for R, H, S in blocks:
            r = R @ cell
            for i in range(n2):
                for j in range(n2):
                    f.write(f"{r[0]:.8f}   {r[1]:.8f}   {r[2]:.8f}   {i + 1}   {j + 1}   "
                            f"{H[i, j].real:.10f}   {H[i, j].imag:.10f}   {S[i, j].real:.10f}\n")


def write_basis(path, labels, spins, offset_by_spin):
    with open(path, "w") as f:
        f.write("#bindex aspecie ax ay az l m spin\n")
        b = 0
        for spin in spins:
            for ia, sym, orb in labels:
                b += 1
                l, m = ORBITAL_LM[orb]
                p = offset_by_spin[ia]
                f.write(f"{b}   {sym}   {p[0]:.8f}   {p[1]:.8f}   {p[2]:.8f}   {l}   {m}   {spin}\n")


def write_wtb_input(outdir, kind, dim, info, params_name):
    sysdim = {1: '"1D"', 2: '"2D"', 3: '"3D"'}[dim] if dim in (1, 2, 3) else '"3D"'
    (outdir / f"input-{kind}.dat").write_text(
        "NTHREADS= 1\n"
        f"SYSDIM= {sysdim}\n"
        "DFT= S\n\n"
        f'PARAMS_FILE= "{params_name}"\n'
        'KPATH_FILE= "kpoints.dat"\n'
        'OUTPUT= "./out/"\n'
        'CALC_DATA= "./out/"\n\n'
        "BANDS= T\n")
    if info["band"]:
        ppl, pts = info["band"]
        with open(outdir / "kpoints.dat", "w") as f:
            f.write(f"{2 * (len(pts) - 1)}\n{ppl + 1}\n")
            for (l1, k1), (l2, k2) in zip(pts[:-1], pts[1:]):
                f.write(f"{k1[0]:.8f} {k1[1]:.8f} {k1[2]:.8f} {l1}\n"
                        f"{k2[0]:.8f} {k2[1]:.8f} {k2[2]:.8f} {l2}\n\n")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("run", help="the .bind input file of the YAeHMOP run (outputs next to it)")
    ap.add_argument("-o", "--outdir", default=".", help="where to write the WanTiBEXOS files")
    ap.add_argument("--zero", default="auto", help="energy zero: auto, fermi, midgap or a value in eV")
    ap.add_argument("--name", default="huckel", help="prefix of the TB file (default huckel -> huckel-NP.dat)")
    ap.add_argument("--cutoff", type=float, default=1e-10,
                    help="R vectors with all |H|,|S| below this are dropped (default 1e-10)")
    ap.add_argument("--allow-aliasing", action="store_true",
                    help="convert even if the mesh is smaller than 2r+1 (H(R) is then approximate)")
    a = ap.parse_args()

    run = Path(a.run)
    outdir = Path(a.outdir); outdir.mkdir(parents=True, exist_ok=True)
    info = read_bind_input(run)
    if info["kpoints"] is None:
        raise SystemExit(f"{run}: no 'K points' block (the H/S dumps are written on those k points)")
    dim, ov = info["dim"], info["overlaps"]
    kpts = info["kpoints"]
    N, alias, kdev = check_mesh(kpts, dim, ov)
    if kdev > 1e-9:
        print(f"WARNING: the k points of the input are rounded (max deviation from j/N {kdev:.1e}); "
              "write them with 12 or more decimals for an exact H(R)", file=sys.stderr)
    if alias and not a.allow_aliasing:
        need = "x".join(str(2 * ov[d] + 1) if d < dim else "1" for d in range(3))
        raise SystemExit(f"the mesh {N[0]}x{N[1]}x{N[2]} is too small for the overlap radius {ov}: "
                         f"use at least {need} (or --allow-aliasing)")
    real_atoms, cell = structure(info)

    base = channel_paths(run)
    if Path(str(run) + ".HAM").exists():
        runs = {"": run}
    elif Path(str(base("_up")) + ".HAM").exists() and Path(str(base("_down")) + ".HAM").exists():
        runs = {"up": base("_up"), "down": base("_down")}
    else:
        raise SystemExit(f"no {run}.HAM nor {base('_up')}.HAM/{base('_down')}.HAM: "
                         "run YAeHMOP with 'Dump Hamil' and 'Dump Overlap'")

    data = {}
    for tag, r in runs.items():
        out = read_out(str(r) + ".out")
        Hk, Sk = read_dump(str(r) + ".HAM"), read_dump(str(r) + ".OV")
        if len(Hk) != len(kpts):
            raise SystemExit(f"{r}.HAM has {len(Hk)} k points, the input has {len(kpts)}")
        data[tag] = dict(out=out, Hk=Hk, Sk=Sk)
    out0 = data[next(iter(data))]["out"]
    labels = orbital_labels(real_atoms, out0["shells"])
    n = len(labels)
    # lattice check against the .out (6 decimals)
    for j, v in enumerate(out0["lattice"][:dim]):
        if np.abs(v - cell[j]).max() > 2e-4:
            raise SystemExit(f"lattice vector {j+1} differs from the .out: {cell[j]} vs {v}")

    if "up" in data:
        kind = "SP"
    elif out0["soc_states"]:
        kind = "SOC"
    else:
        kind = "NP"
    nexp = 2 * n if kind == "SOC" else n
    for tag, d in data.items():
        if d["Hk"].shape[1] != nexp:
            raise SystemExit(f"matrix size {d['Hk'].shape[1]} but the basis has {nexp} functions")

    # electrons -> occupied states (for the midgap)
    if info["electrons"] is not None:
        nelec = info["electrons"]
    else:
        nelec = sum(out0["shells"][s.upper()]["valence"] for s, _ in real_atoms) - info["charge"]

    report = [f"YAeHMOP run: {run}", f"type: {kind}", f"basis: {n} orbitals"
              + (f" x 2 spins = {2*n}" if kind != "NP" else ""),
              f"mesh: {N[0]}x{N[1]}x{N[2]} ({len(kpts)} k points), overlap radius {ov}"
              + ("  (WARNING: smaller than 2r+1, H(R) approximate)" if alias else "")
              + (f"  (WARNING: k points rounded by up to {kdev:.1e})" if kdev > 1e-9 else "")]
    blocks_by_tag = {}
    # spin polarized with a common Fermi level: one energy zero for both channels,
    # so that the up and down bands keep their relative position
    ef_common = common_fermi(run) if kind == "SP" and a.zero in ("auto", "fermi") else None
    if kind == "SP" and ef_common is None and a.zero in ("auto", "fermi"):
        print("note: fixed-moment run (no common Fermi level in the .spin file): each spin "
              "channel gets its own Fermi level as energy zero; with bind --common-fermi "
              "both channels share one zero", file=sys.stderr)
    for tag, d in data.items():
        lab = f"channel {tag}" if tag else kind
        if kind == "SP":
            nocc = d["out"]["channel_electrons"]
        elif kind == "SOC":
            nocc = nelec
        else:
            nocc = nelec / 2
        if ef_common is not None:
            e0, why = ef_common, f"common Fermi level of {run}.spin (--common-fermi)"
        else:
            e0, why = energy_zero(a.zero, d["out"], d["Hk"], d["Sk"], nocc, lab)
        cells = real_space(d["Hk"], d["Sk"], kpts, dim, ov)
        eH, eS = rebuild_error(cells, d["Hk"], d["Sk"], kpts)
        imH = max(np.abs(h.imag).max() for _, h, _ in cells)
        imS = max(np.abs(s.imag).max() for _, _, s in cells)
        if kind != "SOC" and imH > 1e-6:
            raise SystemExit(f"{lab}: H(R) is not real (max imaginary part {imH:.2e})")
        if imS > 1e-6:
            raise SystemExit(f"{lab}: S(R) is not real (max imaginary part {imS:.2e})")
        kept = [(R, h - e0 * s, s.real) for R, h, s in cells
                if max(np.abs(h).max(), np.abs(s).max()) > a.cutoff]
        blocks_by_tag[tag] = kept
        report += [f"[{lab}] energy zero E0 = {e0:.6f} eV ({why}); file holds H(R) - E0*S(R), Fermi level 0",
                   f"[{lab}] rebuilt H(k), S(k) on the mesh: max error {eH:.2e} eV, {eS:.2e}; "
                   f"max |Im H(R)| {imH:.2e}, |Im S(R)| {imS:.2e}; {len(kept)} R vectors kept"]

    pos = {ia: p for ia, (s, p) in enumerate(real_atoms)}
    if kind == "SP":
        up, dn = blocks_by_tag["up"], blocks_by_tag["down"]
        keys = sorted({tuple(R) for R, _, _ in up} | {tuple(R) for R, _, _ in dn},
                      key=lambda R: (sum(abs(x) for x in R), R))
        du = {tuple(R): (h, s) for R, h, s in up}; dd = {tuple(R): (h, s) for R, h, s in dn}
        Z = np.zeros((n, n))
        blocks = []
        for R in keys:
            hu, su = du.get(R, (Z, Z)); hd, sd = dd.get(R, (Z, Z))
            blocks.append((np.array(R), np.block([[hu, Z], [Z, hd]]), np.block([[su, Z], [Z, sd]])))
        write_tb(outdir / f"{a.name}-SP.dat", "SP", cell, 2 * n, blocks)
        write_basis(outdir / "basis_set-SP", labels, (1, -1), pos)
    elif kind == "SOC":
        write_tb(outdir / f"{a.name}-SOC.dat", "SOC", cell, 2 * n, blocks_by_tag[""])
        write_basis(outdir / "basis_set-SOC", labels, (1, -1), pos)
    else:
        write_tb(outdir / f"{a.name}-NP.dat", "NP", cell, n, blocks_by_tag[""])
        write_basis(outdir / "basis_set-NP", labels, (0,), pos)
    write_wtb_input(outdir, kind, dim, info, f"{a.name}-{kind}.dat")
    report += [f"written: {a.name}-{kind}.dat, basis_set-{kind}, input-{kind}.dat"
               + (", kpoints.dat (band path of the input)" if info["band"] else ""),
               "run WanTiBEXOS with:  wtb.x < input-" + kind + ".dat"]
    (outdir / f"yaehmop2wtb-{kind}.log").write_text("\n".join(report) + "\n")
    print("\n".join(report))


if __name__ == "__main__":
    main()
