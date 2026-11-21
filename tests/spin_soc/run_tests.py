#!/usr/bin/env python3
"""Tests of the spin-polarized (--parm-up/--parm-down) and SOC (--soc) runs.

    python3 tests/spin_soc/run_tests.py --bind build/bind [--ref-bind build-ref/bind]

Needs numpy and scipy. Every run happens in a temporary directory.
--ref-bind (a build of the original, unmodified program) adds a regression
test: without the new options the outputs must be byte-identical.

The SOC reference is an independent numpy version of the on-site lambda L.S
model (the same one as parameter_fitter/soc.py in yaehmop-bandfit), applied to
the spinless H(k), S(k) dumps of the program.  The same model gives the
reference Mulliken charges of a SOC run (spinor eigenvectors, common filling).

The regression test compares the .out files without the projected DOS curves
(they changed with the fix of the k point weights in the projected DOS, #47),
and a separate test checks that the projections add up to the total DOS.
"""
import argparse, filecmp, re, shutil, subprocess, sys, tempfile
from pathlib import Path
import numpy as np
import scipy.linalg as sl

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
DATA = HERE / "data"
EXAMPLES = ROOT / "examples"
TABLE = ROOT / "tightbind" / "eht_parms.dat"
FAILS = []


def check(name, ok, detail=""):
    print(f"[{'ok' if ok else 'FAIL'}] {name}" + (f"  ({detail})" if detail else ""))
    if not ok:
        FAILS.append(name)


def run(bind, args, cwd, expect_ok=True):
    r = subprocess.run([str(bind)] + [str(a) for a in args], cwd=cwd,
                       capture_output=True, text=True)
    if expect_ok and r.returncode != 0:
        raise RuntimeError(f"bind {args} failed:\n{r.stdout}\n{r.stderr}")
    return r


def read_band(path):
    text = Path(path).read_text()
    nst = int([l for l in text.splitlines() if "orbitals in the unit cell" in l][0].split()[0])
    body = text.split("; Begin band data.")[1].split("#END")[0].strip().splitlines()
    kp, ev, cur = [], [], []
    for l in body:
        if l.startswith("; K point"):
            if cur:
                ev.append(cur)
            cur = []
            kp.append([float(x) for x in l.split(":")[1].split()])
        else:
            cur.append(float(l))
    ev.append(cur)
    ev = np.array(ev)
    assert ev.shape[1] == nst
    return np.array(kp), ev


def read_dump(path):
    """.HAM/.OV: int nk, int n, per k an n*n double array; (i,j), i<j, is m[i,j] - i m[j,i]."""
    raw = Path(path).read_bytes()
    nk, n = np.frombuffer(raw[:8], dtype=np.int32)
    data = np.frombuffer(raw[8:8 + 8 * nk * n * n], dtype=np.float64).reshape(nk, n, n)
    iu = np.triu_indices(n, 1)
    out = np.zeros((nk, n, n), complex)
    for q in range(nk):
        m = data[q]
        H = np.diag(np.diag(m)).astype(complex)
        H[iu] = m[iu] - 1j * m.T[iu]
        H[(iu[1], iu[0])] = np.conj(H[iu])
        out[q] = H
    return out


def grep(path, key):
    return [l for l in Path(path).read_text().splitlines() if key in l]


def fermi(path):
    return float(grep(path, "#Fermi_Energy")[-1].split()[-1])


def net_charges(path):
    lines = Path(path).read_text().splitlines()
    i = lines.index("# AVERAGE NET CHARGES")
    out = []
    for l in lines[i + 1:]:
        p = l.split()
        if len(p) != 3 or not p[0].isdigit():
            break
        out.append(float(p[2]))
    return np.array(out)


def curves(path):
    """all '#BEGIN CURVE' ... '#END CURVE' blocks of an .out: list of (weight, E) arrays"""
    out, cur = [], None
    for l in Path(path).read_text().splitlines():
        if l.startswith("#BEGIN CURVE"):
            cur = []
        elif l.startswith("#END CURVE"):
            out.append(np.array(cur).reshape(-1, 2)); cur = None
        elif cur is not None and l.strip() and not l.startswith(";"):
            cur.append([float(x) for x in l.split()[:2]])
    return out


def binned(c, de=0.05, emin=-60.0, emax=40.0):
    h, _ = np.histogram(c[:, 1], bins=np.arange(emin, emax, de), weights=c[:, 0])
    return h


def without_pdos(path):
    """the .out with the projected DOS curves left out"""
    keep, skip = [], False
    for l in Path(path).read_text().splitlines():
        if l.startswith("### PROJECTED DENSITY OF STATES"):
            skip = True
        if not skip:
            keep.append(l)
        if skip and l.startswith("#END CURVE"):
            skip = False
    return keep


def avg_coops(path):
    lines = Path(path).read_text().splitlines()
    i = max(k for k, l in enumerate(lines) if l.startswith("; Average Values of COOP"))
    out = []
    for l in lines[i + 1:]:
        if not l.startswith("("):
            break
        out.append(float(l.split(":")[-1]))
    return np.array(out)


def molecular_charges(path):
    lines = Path(path).read_text().splitlines()
    i = max(k for k, l in enumerate(lines) if "Net Atomic Charges" in l)
    out = []
    for l in lines[i + 1:]:
        p = l.split()
        if len(p) != 3 or not p[0].isdigit():
            break
        out.append(float(p[2]))
    return np.array(out)


def spin_summary(path):
    out = {}
    for l in Path(path).read_text().splitlines():
        if l.startswith("#") and ":" in l:
            k, v = l[1:].split(":", 1)
            try:
                out[k.strip()] = float(v.split()[0])
            except (ValueError, IndexError):
                out[k.strip()] = v.strip()
    return out


# ---------------- independent SOC model (numpy) ----------------
REAL = {"px": (1, "c"), "py": (1, "s"), "pz": (0, "0"),
        "dx2-y2": (2, "c"), "dz2": (0, "0"), "dxy": (2, "s"), "dxz": (1, "c"), "dyz": (1, "s"),
        "fz3": (0, "0"), "fxz2": (1, "c"), "fyz2": (1, "s"), "fxyz": (2, "s"),
        "fz(x2-y2)": (2, "c"), "fx(x2-3y2)": (3, "c"), "fy(3x2-y2)": (3, "s")}
SHELL = {"s": ["s"], "p": ["px", "py", "pz"], "d": ["dx2-y2", "dz2", "dxy", "dxz", "dyz"],
         "f": ["fz3", "fxz2", "fyz2", "fxyz", "fz(x2-y2)", "fx(x2-3y2)", "fy(3x2-y2)"]}


def angular_momentum(l, names):
    U = np.zeros((len(names), 2 * l + 1), complex)
    for a, nm in enumerate(names):
        m, kind = REAL[nm]
        if kind == "0":
            U[a, l] = 1
        elif kind == "c":
            U[a, l + m] = (-1) ** m / np.sqrt(2); U[a, l - m] = 1 / np.sqrt(2)
        else:
            U[a, l - m] = 1j / np.sqrt(2); U[a, l + m] = -1j * (-1) ** m / np.sqrt(2)
    m = np.arange(-l, l + 1)
    Lp = np.zeros((2 * l + 1, 2 * l + 1), complex)
    for i, mm in enumerate(m[:-1]):
        Lp[i + 1, i] = np.sqrt(l * (l + 1) - mm * (mm + 1))
    Ls = ((Lp + Lp.conj().T) / 2, (Lp - Lp.conj().T) / 2j, np.diag(m).astype(complex))
    return [U.conj() @ L @ U.T for L in Ls]


def soc_reference(symbols, shells_of, lambdas):
    """H_SOC (2n x 2n, spin up first); shells_of[el] = shells in YAeHMOP order."""
    blocks, n = [], 0
    for at, el in enumerate(symbols):
        for sh in shells_of[el]:
            blocks.append((el, sh, n)); n += len(SHELL[sh])
    H = np.zeros((2 * n, 2 * n), complex)
    for el, sh, s0 in blocks:
        lam = lambdas.get((el.upper(), sh), 0.0)
        if sh == "s" or lam == 0.0:
            continue
        Lx, Ly, Lz = angular_momentum("spdf".index(sh), SHELL[sh])
        u = np.arange(s0, s0 + len(SHELL[sh])); d = u + n
        H[np.ix_(u, u)] += lam / 2 * Lz;  H[np.ix_(d, d)] -= lam / 2 * Lz
        H[np.ix_(u, d)] += lam / 2 * (Lx - 1j * Ly); H[np.ix_(d, u)] += lam / 2 * (Lx + 1j * Ly)
    return H


def shells_from_table(path):
    out = {}
    for l in Path(path).read_text().splitlines():
        p = l.split()
        if len(p) >= 6 and not l.startswith(";") and p[5] in "spdf":
            out.setdefault(p[0].upper(), []).append(p[5])
    return {el: [s for s in "spdf" if s in v] for el, v in out.items()}


# ---------------- tests ----------------
def diamond_input(dst):
    """examples/spin_soc/diamond.bind (C from the table, H/S dumps, energy levels)."""
    shutil.copy(EXAMPLES / "spin_soc" / "diamond.bind", dst / "diamond.bind")
    for f in ("C_up.dat", "C_down.dat", "soc.dat"):
        shutil.copy(EXAMPLES / "spin_soc" / f, dst / f)


def test_regression(bind, ref_bind, tmp):
    for name in ("diamond", "XeF4", "Te2Br", "CoNb4Si"):
        out = {}
        for tag, exe in (("ref", ref_bind), ("new", bind)):
            d = tmp / f"reg_{tag}"; d.mkdir(exist_ok=True)
            shutil.copy(EXAMPLES / name, d / name)
            run(exe, [name, TABLE], d)
            out[tag] = d
        for ext in (".out", ".band"):
            a, b = out["ref"] / (name + ext), out["new"] / (name + ext)
            if not a.exists():
                continue
            if ext == ".out" and "PROJECTED DENSITY" in a.read_text():
                check(f"regression {name}{ext} identical to the original program "
                      "(projected DOS curves left out)", without_pdos(a) == without_pdos(b))
            else:
                check(f"regression {name}{ext} identical to the original program",
                      filecmp.cmp(a, b, shallow=False))


def test_spin_same_tables(bind, tmp):
    d = tmp / "sp_same"; d.mkdir()
    diamond_input(d)
    ref = d / "ref"; ref.mkdir()
    shutil.copy(d / "diamond.bind", ref)
    run(bind, ["diamond.bind", TABLE], ref)
    run(bind, ["--parm-up", TABLE, "--parm-down", TABLE, "diamond.bind"], d)
    for ch in ("up", "down"):
        for ext in (".band", ".HAM", ".OV"):
            check(f"spin, same tables: {ch} {ext} = spinless",
                  filecmp.cmp(d / f"diamond_{ch}.bind{ext}", ref / f"diamond.bind{ext}", shallow=False))
        check(f"spin, same tables: {ch} Fermi = spinless",
              fermi(d / f"diamond_{ch}.bind.out") == fermi(ref / "diamond.bind.out"))
        e = float(grep(d / f"diamond_{ch}.bind.out", "#ELECTRONS_THIS_CHANNEL")[0].split()[-1])
        check(f"spin: {ch} channel has half of the electrons", e == 4.0, f"{e}")


def test_spin_different_tables(bind, tmp):
    d = tmp / "sp_diff"; d.mkdir()
    up, dn = d / "up.dat", d / "dn.dat"
    shutil.copy(TABLE, up)
    txt = TABLE.read_text()
    old = "CO      27      9      2     3     d  -13.180"
    assert old in txt
    dn.write_text(txt.replace(old, "CO      27      9      2     3     d  -12.180"))
    src = (EXAMPLES / "CoNb4Si").read_text()
    src = src.replace("\nprint\n", "\nprint\nnet charges\n").replace("\ncharge \n", "\nDump Hamil\nDump Overlap\ncharge \n")
    (d / "CoNb4Si.bind").write_text(src)
    ref = d / "ref"; ref.mkdir()
    shutil.copy(d / "CoNb4Si.bind", ref)
    run(bind, ["CoNb4Si.bind", up], ref)
    run(bind, ["--parm-up", up, "--parm-down", dn, "CoNb4Si.bind"], d)
    check("spin, different tables: up .band = spinless (up table)",
          filecmp.cmp(d / "CoNb4Si_up.bind.band", ref / "CoNb4Si.bind.band", shallow=False))
    check("spin, different tables: up .HAM = spinless (up table)",
          filecmp.cmp(d / "CoNb4Si_up.bind.HAM", ref / "CoNb4Si.bind.HAM", shallow=False))
    check("spin, different tables: down .HAM differs",
          not filecmp.cmp(d / "CoNb4Si_down.bind.HAM", ref / "CoNb4Si.bind.HAM", shallow=False))
    check("spin, different tables: down .OV = spinless (same zetas)",
          filecmp.cmp(d / "CoNb4Si_down.bind.OV", ref / "CoNb4Si.bind.OV", shallow=False))
    q_up, q_ref = net_charges(d / "CoNb4Si_up.bind.out"), net_charges(ref / "CoNb4Si.bind.out")
    check("spin: net charges of the up channel = half of the spinless ones",
          np.allclose(2 * q_up, q_ref, atol=2e-6), f"max {np.abs(2*q_up-q_ref).max():.1e}")


def test_spin_const(bind, tmp):
    d = tmp / "sp_const"; d.mkdir()
    diamond_input(d)
    ref = d / "ref"; ref.mkdir()
    txt = (d / "diamond.bind").read_text()
    (ref / "diamond.bind").write_text(txt.replace("\nprint\n", "\nTHE CONST 1.60\nprint\n", 1))
    run(bind, ["diamond.bind", TABLE], ref)
    run(bind, ["--parm-up", TABLE, "--parm-down", TABLE, "--const-up", "1.60",
               "diamond.bind"], d)
    check("spin: --const-up 1.60 = spinless run with THE CONST 1.60",
          filecmp.cmp(d / "diamond_up.bind.band", ref / "diamond.bind.band", shallow=False))
    check("spin: down channel keeps the K of the input",
          not filecmp.cmp(d / "diamond_down.bind.band", ref / "diamond.bind.band", shallow=False))


def test_soc_lambda0(bind, tmp):
    d = tmp / "soc0"; d.mkdir()
    diamond_input(d)
    (d / "soc0.dat").write_text("C p 0.0\n")
    ref = d / "ref"; ref.mkdir()
    shutil.copy(d / "diamond.bind", ref)
    run(bind, ["diamond.bind", TABLE], ref)
    run(bind, ["--soc", "soc0.dat", "diamond.bind", TABLE], d)
    _, e0 = read_band(ref / "diamond.bind.band")
    _, e1 = read_band(d / "diamond.bind.band")
    diff = np.abs(np.sort(np.repeat(e0, 2, axis=1), 1) - e1).max()
    check("SOC lambda=0: 2n bands = spinless bands twice", diff == 0.0, f"{diff:.1e}")
    I2 = np.eye(2)
    for ext in (".HAM", ".OV"):
        A, B = read_dump(ref / f"diamond.bind{ext}"), read_dump(d / f"diamond.bind{ext}")
        diff = max(np.abs(np.kron(I2, a) - b).max() for a, b in zip(A, B))
        check(f"SOC lambda=0: spinor {ext} = I2 x spinless", diff == 0.0, f"{diff:.1e}")
    check("SOC lambda=0: Fermi level unchanged",
          fermi(d / "diamond.bind.out") == fermi(ref / "diamond.bind.out"))


def test_soc_mos2(bind, tmp):
    d = tmp / "mos2"; d.mkdir()
    params, socf = DATA / "MoS2_params.dat", DATA / "MoS2_soc.dat"
    shutil.copy(DATA / "MoS2_path.bind", d / "path.bind")
    run(bind, ["path.bind", params], d)
    _, e0 = read_band(d / "path.bind.band")
    run(bind, ["--soc", socf, "path.bind", params], d)
    kp, esoc = read_band(d / "path.bind.band")
    # the same k points as a K points list: spinless and spinor dumps
    head = (d / "path.bind").read_text().split("\nBand\n")[0]
    klist = "\n".join(f"{k[0]:.6f} {k[1]:.6f} {k[2]:.6f} 1.0" for k in kp)
    (d / "kl.bind").write_text(head + f"\n\nDump Hamil\nDump Overlap\n\nK points\n{len(kp)}\n{klist}\n")
    run(bind, ["kl.bind", params], d)
    H, S = read_dump(d / "kl.bind.HAM"), read_dump(d / "kl.bind.OV")
    run(bind, ["--soc", socf, "kl.bind", params], d)
    Hc, Sc = read_dump(d / "kl.bind.HAM"), read_dump(d / "kl.bind.OV")
    lam = {}
    for l in socf.read_text().splitlines():
        if l.strip() and not l.startswith("#"):
            el, sh, v = l.split(); lam[(el.upper(), sh)] = float(v)
    Hsoc = soc_reference(["MO", "S", "S"], shells_from_table(params), lam)
    I2 = np.eye(2)
    dH = max(np.abs(hc - (np.kron(I2, h) + Hsoc)).max() for hc, h in zip(Hc, H))
    dS = max(np.abs(sc - np.kron(I2, s)).max() for sc, s in zip(Sc, S))
    check("SOC MoS2: spinor H(k) = I2 x H(k) + H_SOC (numpy model)", dH < 1e-10, f"{dH:.1e}")
    check("SOC MoS2: spinor S(k) = I2 x S(k)", dS < 1e-12, f"{dS:.1e}")
    epy = np.array([sl.eigh(np.kron(I2, h) + Hsoc, np.kron(I2, s), eigvals_only=True)
                    for h, s in zip(H, S)])
    nocc = 18                       # Mo 6 + 2 x S 6 electrons, one per spinor state
    vbm = esoc[:, :nocc].max()
    win = (esoc > vbm - 6) & (esoc < vbm + 4)
    dE = np.abs(esoc - epy)[win].max()
    check("SOC MoS2: bands (VBM-6..VBM+4 eV) = numpy model", dE < 1e-4, f"{dE:.1e} eV")
    iK = np.argmin(np.linalg.norm(kp - np.array([1 / 3, 1 / 3, 0]), axis=1))
    split = esoc[iK, nocc - 1] - esoc[iK, nocc - 2]
    gap0 = e0[:, nocc // 2:].min() - e0[:, :nocc // 2].max()
    gap = esoc[:, nocc:].min() - esoc[:, :nocc].max()
    check("SOC MoS2: VB splitting at K = 0.165 eV (yaehmop-bandfit fit)", abs(split - 0.1650) < 1e-3, f"{split:.4f}")
    check("SOC MoS2: gap 2.313 -> 2.268 eV (yaehmop-bandfit fit)",
          abs(gap0 - 2.3130) < 1e-3 and abs(gap - 2.2681) < 1e-3, f"{gap0:.4f} -> {gap:.4f}")
    check("SOC MoS2: .SOC file written", (d / "kl.bind.SOC").exists())


def conb4si_input(d, extra=""):
    src = (EXAMPLES / "CoNb4Si").read_text()
    src = src.replace("\nprint\n", "\nprint\nnet charges\n" + extra)
    (d / "CoNb4Si.bind").write_text(src)


def test_pdos_sum_rule(bind, tmp):
    d = tmp / "pdos"; d.mkdir()
    shutil.copy(EXAMPLES / "CoNb4Si", d / "CoNb4Si")
    run(bind, ["CoNb4Si", TABLE], d)
    c = curves(d / "CoNb4Si.out")
    tot, proj = binned(c[0]), sum(binned(x) for x in c[1:4])   # 3 projections cover all atoms
    diff = np.abs(tot - proj).max()
    check("projected DOS: the atom projections add up to the total DOS", diff < 1e-4, f"{diff:.1e}")


def test_soc_props_lambda0(bind, tmp):
    """SOC with lambda = 0: properties equal the spinless ones (each level twice)"""
    d = tmp / "socprop0"; d.mkdir()
    conb4si_input(d)
    (d / "soc0.dat").write_text("Co d 0.0\nNb d 0.0\nSi p 0.0\n")
    ref = d / "ref"; ref.mkdir()
    shutil.copy(d / "CoNb4Si.bind", ref)
    run(bind, ["CoNb4Si.bind", TABLE], ref)
    run(bind, ["--soc", "soc0.dat", "CoNb4Si.bind", TABLE], d)
    a, b = ref / "CoNb4Si.bind.out", d / "CoNb4Si.bind.out"
    check("SOC lambda=0: Fermi level = spinless", fermi(a) == fermi(b))
    check("SOC lambda=0: average energy = spinless",
          grep(a, "#Average Energy") == grep(b, "#Average Energy"))
    dq = np.abs(net_charges(a) - net_charges(b)).max()
    check("SOC lambda=0: net charges = spinless", dq < 2e-6, f"{dq:.1e}")
    dc = np.abs(avg_coops(a) - avg_coops(b)).max()
    check("SOC lambda=0: average COOPs = spinless", dc < 2e-6, f"{dc:.1e}")
    ca, cb = curves(a), curves(b)
    dd = np.abs(2 * binned(ca[0]) - binned(cb[0])).max()
    check("SOC lambda=0: total DOS = 2 x spinless (twice the states)", dd < 1e-9, f"{dd:.1e}")
    dp = max(np.abs(2 * binned(x) - binned(y)).max() for x, y in zip(ca[1:4], cb[1:4]))
    check("SOC lambda=0: projected DOS = 2 x spinless", dp < 1e-4, f"{dp:.1e}")
    dco = max(np.abs(binned(x) - binned(y)).max() for x, y in zip(ca[4:], cb[4:]))
    check("SOC lambda=0: COOP curves = spinless", dco < 1e-5, f"{dco:.1e}")
    check("SOC: spin moments printed (zero without magnetism)",
          abs(float(grep(b, "#TOTAL_SPIN_MOMENT")[0].split()[-1])) < 1e-6)


def test_soc_charges_mos2(bind, tmp):
    """SOC Mulliken charges of MoS2 against the numpy spinor model"""
    d = tmp / "mos2q"; d.mkdir()
    params, socf = DATA / "MoS2_params.dat", DATA / "MoS2_soc.dat"
    head = (DATA / "MoS2_path.bind").read_text().split("\nBand\n")[0]
    g = np.linspace(0, 1, 6, endpoint=False)
    kl = [(x, y, 0.0) for x in g for y in g]
    klist = "\n".join(f"{k[0]:.6f} {k[1]:.6f} {k[2]:.6f} 1.0" for k in kl)
    (d / "q.bind").write_text(head + "\n\nprint\nnet charges\nend_print\nAverage Properties\n"
                              "Electrons\n18\n"
                              f"Dump Hamil\nDump Overlap\n\nK points\n{len(kl)}\n{klist}\n")
    run(bind, ["q.bind", params], d)
    H, S = read_dump(d / "q.bind.HAM"), read_dump(d / "q.bind.OV")
    run(bind, ["--soc", socf, "q.bind", params], d)
    q_c = net_charges(d / "q.bind.out")
    lam = {}
    for l in socf.read_text().splitlines():
        if l.strip() and not l.startswith("#"):
            el, sh, v = l.split(); lam[(el.upper(), sh)] = float(v)
    shells = shells_from_table(params)
    Hsoc = soc_reference(["MO", "S", "S"], shells, lam)
    n = H.shape[1]
    owner = np.concatenate([[a] * sum(len(SHELL[sh]) for sh in shells[el])
                            for a, el in enumerate(["MO", "S", "S"])])
    I2 = np.eye(2)
    E, W = [], []
    for h, s in zip(H, S):
        Ss = np.kron(I2, s)
        e, c = sl.eigh(np.kron(I2, h) + Hsoc, Ss)
        w = np.real(np.conj(c) * (Ss @ c))            # Mulliken weight of each spinor AO
        E.append(e); W.append(w[:n] + w[n:])           # add both spin parts of orbital mu
    E, W = np.array(E), np.array(W)
    nel = 18 * len(kl)
    order = np.argsort(E, axis=None, kind="stable")[:nel]
    pop = np.zeros(3)
    for idx in order:
        k, j = np.unravel_index(idx, E.shape)
        for a in range(3):
            pop[a] += W[k][owner == a, j].sum()
    q_py = np.array([6, 6, 6]) - pop / len(kl)
    dq = np.abs(q_py - q_c).max()
    check("SOC MoS2: net charges = numpy spinor model", dq < 2e-5, f"{dq:.1e}")
    check("SOC MoS2: net charges add up to zero", abs(q_c.sum()) < 1e-5, f"{q_c.sum():.1e}")


def test_soc_molecule(bind, tmp):
    d = tmp / "socmol"; d.mkdir()
    src = (EXAMPLES / "XeF4").read_text()
    src = re.sub(r"\nfmo\n(.*\n){9}", "\n", src)                   # no FMO block
    src = re.sub(r"mo print\n15\n(\d+\n){15}", "", src)             # no .MO file
    assert "fmo" not in [l.strip().lower() for l in src.splitlines()]
    (d / "x.bind").write_text(src)
    (d / "soc0.dat").write_text("Xe p 0.0\nF p 0.0\n")
    (d / "soc.dat").write_text("Xe p 0.8\nF p 0.03\n")
    ref = d / "ref"; ref.mkdir()
    shutil.copy(d / "x.bind", ref)
    run(bind, ["x.bind", TABLE], ref)
    run(bind, ["--soc", "soc0.dat", "x.bind", TABLE], d)
    a, b = ref / "x.bind.out", d / "x.bind.out"
    dq = np.abs(molecular_charges(a) - molecular_charges(b)).max()
    check("SOC molecule lambda=0: net charges = spinless", dq < 2e-6, f"{dq:.1e}")
    ta = float(grep(a, "Total_Energy")[0].split()[-1]); tb = float(grep(b, "Total_Energy")[0].split()[-1])
    check("SOC molecule lambda=0: total energy = spinless", abs(ta - tb) < 1e-3, f"{ta} {tb}")
    run(bind, ["--soc", "soc.dat", "x.bind", TABLE], d)
    lev = [float(l.split()[1]) for l in Path(b).read_text().splitlines() if ":--->" in l]
    lev = np.array(lev)
    check("SOC molecule: 2n levels in Kramers pairs", len(lev) == 40 and
          np.abs(lev[0::2] - lev[1::2]).max() < 1e-4)


def test_soc_thin_walsh(bind, tmp):
    d = tmp / "socthin"; d.mkdir()
    diamond_input(d)
    txt = (d / "diamond.bind").read_text().replace("average properties\n", "")
    (d / "fat.bind").write_text(txt)
    (d / "thin.bind").write_text(txt.replace("Dump Hamil\n", "thin\nDump Hamil\n"))
    run(bind, ["--soc", "soc.dat", "fat.bind", TABLE], d)
    run(bind, ["--soc", "soc.dat", "thin.bind", TABLE], d)
    _, ef = read_band(d / "fat.bind.band"); _, et = read_band(d / "thin.bind.band")
    check("SOC THIN mode: bands = FAT mode", np.abs(ef - et).max() < 1e-9)
    # Walsh diagram of a molecule with SOC (lambda = 0: same total energy and charge)
    w = (EXAMPLES / "ethane.b").read_text().rstrip() + (
        "\n\nprint\norbital energy walsh\norbital 7 1\nenergy walsh\norbital 1 1\n"
        "net charge walsh\natom 3 1\nend_print\n")
    (d / "eth").write_text(w)
    ref = d / "ref"; ref.mkdir(); shutil.copy(d / "eth", ref)
    (d / "c0.dat").write_text("C p 0.0\n")
    run(bind, ["eth", TABLE], ref)
    run(bind, ["--soc", "c0.dat", "eth", TABLE], d)
    A = np.loadtxt(ref / "eth.walsh"); B = np.loadtxt(d / "eth.walsh")
    # columns: variable, net charge, total energy, orbital energy
    check("SOC Walsh lambda=0: net charge and total energy = spinless",
          A.shape == B.shape and np.abs(A[:, :3] - B[:, :3]).max() < 1e-5)


def fe_inputs(d):
    shutil.copy(DATA / "Fe_bcc.bind", d / "Fe.bind")
    shutil.copy(DATA / "Fe_up.dat", d / "Fe_up.dat")
    shutil.copy(DATA / "Fe_down.dat", d / "Fe_down.dat")
    up = [l.split() for l in (DATA / "Fe_up.dat").read_text().splitlines() if l.startswith("FE")]
    dn = [l.split() for l in (DATA / "Fe_down.dat").read_text().splitlines() if l.startswith("FE")]
    rows = ["; Fe spin-down Hii with the spin-up zetas"]
    for u, w in zip(up, dn):
        u = list(u); u[6] = w[6]; rows.append("  ".join(u))
    (d / "Fe_dn_samez.dat").write_text("\n".join(rows + ["END"]) + "\n")
    return ["--const-up", "1.752806708980162", "--const-down", "1.7490786283756976"]


def test_spin_common_fermi(bind, tmp):
    d = tmp / "cfermi"; d.mkdir()
    K = fe_inputs(d)
    run(bind, ["--common-fermi"] + K + ["--parm-up", "Fe_up.dat", "--parm-down", "Fe_down.dat",
                                         "Fe.bind"], d)
    S = spin_summary(d / "Fe.bind.spin")
    nu, nd = S["ELECTRONS_UP"], S["ELECTRONS_DOWN"]
    check("spin --common-fermi: N_up + N_down = electrons of the input", abs(nu + nd - 8) < 1e-9,
          f"{nu} + {nd}")
    check("spin --common-fermi: common E_F = highest channel E_F",
          S["Fermi_Energy"] == max(S["Fermi_Energy_up"], S["Fermi_Energy_down"]))
    check("spin --common-fermi: moment = N_up - N_down, Fe ferromagnetic",
          abs(S["MAGNETIC_MOMENT"] - (nu - nd)) < 1e-9 and S["MAGNETIC_MOMENT"] > 1.0,
          f"{S['MAGNETIC_MOMENT']:.3f}")
    fixed = d / "fixed"; fixed.mkdir()
    for f in ("Fe.bind", "Fe_up.dat", "Fe_down.dat"):
        shutil.copy(d / f, fixed / f)
    run(bind, K + ["--nelec-up", repr(nu), "--nelec-down", repr(nd), "--parm-up", "Fe_up.dat",
                   "--parm-down", "Fe_down.dat", "Fe.bind"], fixed)
    for ch in ("up", "down"):
        check(f"spin --common-fermi: {ch} channel = run with --nelec-* fixed to its N",
              fermi(d / f"Fe_{ch}.bind.out") == fermi(fixed / f"Fe_{ch}.bind.out") and
              net_charges(d / f"Fe_{ch}.bind.out").tolist() ==
              net_charges(fixed / f"Fe_{ch}.bind.out").tolist())
    F = spin_summary(fixed / "Fe.bind.spin")
    check("spin summary without --common-fermi: COMMON_FERMI no", str(F["COMMON_FERMI"]).startswith("no"))


def test_spin_coop_scale(bind, tmp):
    """each channel holds half of the electrons: COOP curves and values are half"""
    d = tmp / "spcoop"; d.mkdir()
    conb4si_input(d)
    ref = d / "ref"; ref.mkdir(); shutil.copy(d / "CoNb4Si.bind", ref)
    run(bind, ["CoNb4Si.bind", TABLE], ref)
    run(bind, ["--parm-up", TABLE, "--parm-down", TABLE, "CoNb4Si.bind"], d)
    a, b = ref / "CoNb4Si.bind.out", d / "CoNb4Si_up.bind.out"
    dv = np.abs(2 * avg_coops(b) - avg_coops(a)).max()
    check("spin, same tables: average COOP of a channel = half of the spinless", dv < 2e-6, f"{dv:.1e}")
    dc = max(np.abs(2 * binned(x) - binned(y)).max() for x, y in zip(curves(b)[4:], curves(a)[4:]))
    check("spin, same tables: COOP curve of a channel = half of the spinless", dc < 1e-5, f"{dc:.1e}")


def test_soc_spin(bind, tmp):
    d = tmp / "socspin"; d.mkdir()
    K = fe_inputs(d)
    (d / "fe0.dat").write_text("Fe d 0.0\n")
    (d / "fe.dat").write_text("Fe d 0.06\n")
    run(bind, ["--common-fermi"] + K + ["--parm-up", "Fe_up.dat", "--parm-down", "Fe_dn_samez.dat",
                                         "Fe.bind"], d)
    S = spin_summary(d / "Fe.bind.spin")
    sd = d / "soc"; sd.mkdir()
    for f in ("Fe.bind", "Fe_up.dat", "Fe_dn_samez.dat", "fe0.dat", "fe.dat"):
        shutil.copy(d / f, sd / f)
    run(bind, ["--soc", "fe0.dat"] + K + ["--parm-up", "Fe_up.dat", "--parm-down",
                                           "Fe_dn_samez.dat", "Fe.bind"], sd)
    o = sd / "Fe.bind.out"
    _, eu = read_band(d / "Fe_up.bind.band"); _, ed = read_band(d / "Fe_down.bind.band")
    _, es = read_band(sd / "Fe.bind.band")
    db = np.abs(np.sort(np.hstack([eu, ed]), 1) - es).max()
    check("SOC + spin, lambda=0: bands = spin up and spin down bands together", db == 0.0, f"{db:.1e}")
    check("SOC + spin, lambda=0: Fermi level = common Fermi level of the spin run",
          fermi(o) == S["Fermi_Energy"])
    m = float(grep(o, "#TOTAL_SPIN_MOMENT")[0].split()[-1])
    check("SOC + spin, lambda=0: spin moment = N_up - N_down of the spin run",
          abs(m - S["MAGNETIC_MOMENT"]) < 1e-5, f"{m:.6f} vs {S['MAGNETIC_MOMENT']:.6f}")
    run(bind, ["--soc", "fe.dat"] + K + ["--parm-up", "Fe_up.dat", "--parm-down",
                                          "Fe_dn_samez.dat", "Fe.bind"], sd)
    m1 = float(grep(o, "#TOTAL_SPIN_MOMENT")[0].split()[-1])
    check("SOC + spin, lambda=0.06: moment changes a little", 0 < abs(m1 - m) < 0.1, f"{m1:.4f}")
    r = run(bind, ["--soc", "fe.dat", "--parm-up", "Fe_up.dat", "--parm-down", "Fe_down.dat",
                   "Fe.bind"], sd, expect_ok=False)
    check("SOC + spin with different zetas refused",
          any("zetas" in l for l in grep(sd / "Fe.bind.status", "FATAL")))


def test_errors(bind, tmp):
    d = tmp / "err"; d.mkdir()
    diamond_input(d)
    shutil.copy(EXAMPLES / "XeF4", d / "XeF4")
    r = run(bind, ["--parm-up", "C_up.dat", "diamond.bind"], d, expect_ok=False)
    check("error: --parm-up without --parm-down", r.returncode != 0 and "both" in r.stderr)
    r = run(bind, ["--soc", "soc.dat", "--parm-up", "C_up.dat", "--parm-down", "C_down.dat",
                   "--nelec-up", "4", "diamond.bind"], d, expect_ok=False)
    check("error: SOC + spin with --nelec-up refused", r.returncode != 0 and "don't apply" in r.stderr)
    r = run(bind, ["--common-fermi", "diamond.bind"], d, expect_ok=False)
    check("error: --common-fermi without spin tables", r.returncode != 0)
    run(bind, ["--soc", "soc.dat", "XeF4", TABLE], d, expect_ok=False)
    check("error: SOC with FMO refused (message names FMO)",
          any("FMO" in l for l in grep(d / "XeF4.status", "FATAL")))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bind", required=True, help="bind built from this branch")
    ap.add_argument("--ref-bind", help="bind of the original program (regression test)")
    ap.add_argument("--keep", action="store_true", help="keep the temporary directory")
    a = ap.parse_args()
    bind = Path(a.bind).resolve()
    tmp = Path(tempfile.mkdtemp(prefix="yaehmop_spin_soc_"))
    try:
        if a.ref_bind:
            test_regression(bind, Path(a.ref_bind).resolve(), tmp)
        test_spin_same_tables(bind, tmp)
        test_spin_different_tables(bind, tmp)
        test_spin_const(bind, tmp)
        test_soc_lambda0(bind, tmp)
        test_soc_mos2(bind, tmp)
        test_pdos_sum_rule(bind, tmp)
        test_soc_props_lambda0(bind, tmp)
        test_soc_charges_mos2(bind, tmp)
        test_soc_molecule(bind, tmp)
        test_soc_thin_walsh(bind, tmp)
        test_spin_common_fermi(bind, tmp)
        test_spin_coop_scale(bind, tmp)
        test_soc_spin(bind, tmp)
        test_errors(bind, tmp)
    finally:
        if a.keep:
            print("outputs in", tmp)
        else:
            shutil.rmtree(tmp, ignore_errors=True)
    print(f"\n{len(FAILS)} failure(s)" if FAILS else "\nall tests passed")
    sys.exit(1 if FAILS else 0)


if __name__ == "__main__":
    main()
