# Spin-polarized runs and spin-orbit coupling

Two additions to `bind`. Runs without the new options give the same output as before.

```
bind [options] <inputfile> [paramfile]

  --parm-up FILE --parm-down FILE   spin-polarized run, one parameter table per channel
  --nelec-up X   --nelec-down X     electrons in each channel (default: half of the
                                    electrons of the input in each one)
  --const-up K   --const-down K     Wolfsberg-Helmholz K of each channel (default:
                                    THE CONST of the input)
  --common-fermi                    one Fermi level for both channels (N_up, N_down
                                    come out of it)
  --soc FILE                        spin-orbit coupling (on-site lambda L.S); with
                                    --parm-up/--parm-down: SOC on top of the two
                                    spin channels
  -h, --help                        list the options
```

## Spin-polarized runs

```
bind --parm-up up.dat --parm-down down.dat MoS2.bind
```

Each spin channel is an independent extended Hückel calculation that uses its own
parameter table (same format as `eht_parms.dat`), so the channels can have different
Hii values, and different zetas too. The two channels run in two separate processes
at the same time.

**Outputs.** The suffix goes before the extension of the input name. Every output gets it:

| spinless              | spin up                  | spin down                  |
|-----------------------|--------------------------|----------------------------|
| `MoS2.bind.out`       | `MoS2_up.bind.out`       | `MoS2_down.bind.out`       |
| `MoS2.bind.band`      | `MoS2_up.bind.band`      | `MoS2_down.bind.band`      |
| `MoS2.bind.HAM`, `.OV`| `MoS2_up.bind.HAM`, `.OV`| `MoS2_down.bind.HAM`, `.OV`|
| `MoS2.bind.status`    | `MoS2_up.bind.status`    | `MoS2_down.bind.status`    |

(An input without an extension, `diamond`, gives `diamond_up.out`, `diamond_up.band`, ...)

**Occupations.** Each channel holds one electron per state. By default each channel gets
half of the electrons of the input (`Electrons` or `Charge` keyword). To fix the moment,
use `--nelec-up` and `--nelec-down`; each channel then has its own Fermi level
(`#Fermi_Energy` in its `.out`).

**Common Fermi level (`--common-fermi`).** The two channel processes send their energies
to the parent process, which finds the Fermi level shared by both channels and gives each
channel its number of electrons (`#ELECTRONS_THIS_CHANNEL_COMMON_FERMI` in its `.out`).
This needs an extended system with `Average Properties`, and excludes `--nelec-*`.

**Summary `<input>.spin`.** Every spin-polarized run writes it (e.g. `Fe.bind.spin`):

```
#COMMON_FERMI: yes
#Fermi_Energy:  -14.376729          (only with --common-fermi)
#Fermi_Energy_up:  -14.376729
#Fermi_Energy_down:  -14.400753
#ELECTRONS_UP:  4.982000
#ELECTRONS_DOWN:  3.018000
#MAGNETIC_MOMENT:  1.964000   ; N_up - N_down per cell (Bohr magnetons)
#Average_Energy:  -134.104417   ; up + down (eV)
; atom      pop_up     pop_down   net_charge   moment
1 FE     4.982000   3.018000  -0.000000   1.964000
```

**COOP.** Each state of a channel holds one electron, so the COOP curves of a channel are
for one electron per state (the original curve assumed two).

**K per channel.** Channel fits can come with different Wolfsberg-Helmholz constants.
`--const-up` / `--const-down` replace `THE CONST` of the input in each channel, and the
`.out` records the value used (`#THE_CONST`).

**Charges.** Net charges are given per channel: half of the atom's valence electrons
minus the channel population. Add up and down to get the total. Charge iteration is
refused in spin-polarized runs.

**H and S.** With `Dump Hamil` / `Dump Overlap` in the input, each channel writes its own
`.HAM`/`.OV` (n×n, same binary format as without spin). If both tables use the same zetas,
the two `.OV` files are identical.

The `.out` header of each channel records the channel (`#SPIN_CHANNEL: up`), its
parameter file and its electron count (`#ELECTRONS_THIS_CHANNEL`).

## Spin-orbit coupling

```
bind --soc soc.dat MoS2.bind params.dat
```

`soc.dat` has one line per element and shell, with lambda in eV. Lines starting with `#`
or `;` are comments. This is the same format as the `soc_parameters.txt` written by
yaehmop-bandfit:

```
# element shell lambda_eV
Mo  d  0.063352
S   p  0.197577
```

Model (the same as `parameter_fitter/soc.py` in yaehmop-bandfit). The coupling acts only
between orbitals of the same shell on the same atom:

    H_SOC = sum over atoms A and shells l >= 1 of  lambda_{A,l} L.S
    H_spinor(k) = I2 (x) H(k) + H_SOC,      S_spinor(k) = I2 (x) S(k)

The spinor basis has 2n states: the n orbitals with spin up (YAeHMOP order), then the same
n orbitals with spin down. A p shell splits into j=3/2 (+lambda/2) and j=1/2 (-lambda). A d
shell splits into j=5/2 (+lambda) and j=3/2 (-3lambda/2). Each spinor state holds one
electron.

**Outputs with SOC**

- `.band`: 2n bands per k point. The header says `2n orbitals in the unit cell.`, so
  existing readers work unchanged.
- `.HAM`, `.OV` (with `Dump Hamil` / `Dump Overlap` and a `K points` block): spinor
  H(k) and S(k), 2n×2n, in the usual binary format. That format is: int nk, int dim, then
  per k a dim×dim double array, where element (i,j) with i<j is `m[i][j] - i*m[j][i]`.
  The `read_matrix_dump` reader of yaehmop-bandfit reads these files directly.
- `.SOC`: text copy of the on-site H_SOC (nonzero elements `i j Re Im`). With it, the
  spinor H(R) can be rebuilt from a spinless run: `I2 (x) H(R) + delta_{R,0} H_SOC`.
- `.out`: the lambdas used for each atom, `#SOC_STATES: 2n`, and all the usual property
  sections: Fermi level, total and projected DOS, COOP/COHP, AO occupations, net charges,
  average overlap populations (with `Average Properties`); energies and occupations,
  wave functions (spin-up and spin-down parts), Mulliken populations, net charges, charge
  matrices, spinor H and S (with the `print` options). Spin moments
  m_z = N_up - N_down per atom: `# AVERAGE SPIN MOMENTS` and `#TOTAL_SPIN_MOMENT`.

**Properties with SOC.** S_spinor = I2 (x) S does not couple the spins, so the Mulliken
populations, projected DOS and COOP of a spinor state are exactly the sums of those of its
spin-up and spin-down parts. Orbital and atom numbers in `Projected DOS` and `COOP` keep
their spinless meaning (an orbital projection adds both spins). Energy levels, charge
matrix rows and Walsh orbital numbers refer to the 2n spinor states. COHP uses the
spinless EHT matrix elements (the on-site SOC term is not in it). This works for crystals,
molecules (levels in Kramers pairs, populations, charges), charge iteration, Walsh
diagrams and THIN mode (bands).

**Not implemented with SOC yet.** The program stops with a message naming the option:
FMO/FCO, `MO print` (.MO files for viewkel), zeta variation, modified Mulliken analysis,
sparse dumps, the molecular electrostatic energy, netCDF, and the Walsh items that follow
single MO coefficients. SOC needs the build with LAPACK (`-DUSE_BLAS_LAPACK=ON`).

## SOC on top of the two spin channels

```
bind --soc soc.dat --parm-up up.dat --parm-down down.dat [--const-up K1 --const-down K2] input.bind
```

One calculation (one process, outputs without suffix) with
H = diag(H_up(k), H_down(k)) + H_SOC and S = I2 (x) S(k). The two tables may differ in
Hii (and K), not in the basis or the zetas (the program stops if they do). The electrons
are those of the input, with one Fermi level; `--nelec-*` and `--common-fermi` do not
apply. The spin moments in the `.out` give the magnetization. With lambda = 0 the result
equals a spin-polarized run with `--common-fermi`.

## Tests

The tests are automated in `tests/spin_soc/run_tests.py` and need only numpy and scipy:

```
python3 tests/spin_soc/run_tests.py --bind build/bind [--ref-bind <bind of master>]
```

The SOC reference inside the script is an independent numpy version of the lambda L.S
model. The MoS2 inputs (`tests/spin_soc/data/`) come from a fit of yaehmop-bandfit
(https://github.com/tromer-unb/yaehmop-bandfit) to HSE06+SOC bands; the bcc Fe tables
(one per spin channel) from a fit to spin-polarized PBE bands.

- Without the new options: `.out` and `.band` of the examples diamond, XeF4, Te2Br and
  CoNb4Si are identical to master (same LAPACK).
- Spin polarization with the same table in both channels: each channel's `.band` and
  `.HAM` equal the spinless run; energies, occupations, net charges and COOPs are exactly
  half. With a different Co d level in the down channel (CoNb4Si), only that channel
  changes.
- `--common-fermi` (bcc Fe): N_up + N_down = N, and each channel equals a run with
  `--nelec-*` fixed to its N.
- SOC with lambda = 0: each spinless band appears twice, the spinor `.HAM`/`.OV` are
  exactly I2 (x) H(k), I2 (x) S(k); Fermi level, average energy, net charges and COOPs
  equal the spinless run, and the total and projected DOS are twice the spinless ones.
  Molecules (XeF4), THIN mode and Walsh diagrams are checked the same way.
- SOC on MoS2: the spinor H(k), S(k) equal the numpy model (difference 0.0); bands agree
  to 1e-5 eV; the VB splitting at K is 0.1650 eV and the gap 2.2681 eV, as in the fit; the
  average net charges equal a numpy calculation from the spinor eigenvectors (4e-7).
- SOC + spin with lambda = 0: bands, Fermi level and moment equal the spin-polarized run
  with `--common-fermi`.

## Example

`examples/spin_soc/`: diamond with a carbon table per channel (`C_up.dat`,
`C_down.dat`) and `soc.dat`:

```
bind diamond.bind C_up.dat                              # spinless
bind --soc soc.dat diamond.bind C_up.dat                # SOC
bind --parm-up C_up.dat --parm-down C_down.dat diamond.bind   # spin polarized
bind --common-fermi --parm-up C_up.dat --parm-down C_down.dat diamond.bind
bind --soc soc.dat --parm-up C_up.dat --parm-down C_down.dat diamond.bind   # SOC + spin
```
