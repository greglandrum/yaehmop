/*******************************************************

Copyright (C) 2026 Alexandre C. Dias
All rights reserved

This file is part of yaehmop.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

1. Redistributions of source code must retain the above copyright
notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright
notice, this list of conditions and the following disclaimer in the
documentation and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

********************************************************************/

/*******************************************************
 * Spin-orbit coupling for the tightbinding program.
 *
 * On-site approximation, as in DFTB+ and in the yaehmop-bandfit (soc.py):
 *
 *     H_SOC = sum over atoms A and shells l >= 1 of  lambda_{A,l} L.S
 *
 * acting only between orbitals of the same shell of the same atom.  The
 * spinor basis has 2n functions: the n orbitals with spin up, then the same
 * n orbitals with spin down.  With S = sigma/2,
 *
 *     up-up = +(lambda/2) Lz            dn-dn = -(lambda/2) Lz
 *     up-dn =  (lambda/2) (Lx - i Ly)   dn-up =  (lambda/2) (Lx + i Ly)
 *
 *     H(k) -> I2 (x) H(k) + H_SOC,      S(k) -> I2 (x) S(k)
 *
 * (orbitals of one atom don't overlap each other, so H_SOC needs no overlap
 * correction).  p levels split into j = 3/2 (+lambda/2) and j = 1/2 (-lambda),
 * d levels into j = 5/2 (+lambda) and j = 3/2 (-3 lambda/2).
 *
 * The real orbitals are the YAeHMOP ones (px py pz; dx2-y2 dz2 dxy dxz dyz;
 * fz3 fxz2 fyz2 fxyz fz(x2-y2) fx(x2-3y2) fy(3x2-y2)), built from complex
 * spherical harmonics with the Condon-Shortley phase.
 *
 * Complex hermitian matrices use the packed layout of the program: element
 * (i,j), i<j, of the physical matrix is mat[i*n+j] - I*mat[j*n+i] (this is
 * what the LAPACK copy in kpoints.c/bands.c and the readers of .HAM/.OV in
 * the yaehmop-bandfit assume), the diagonal is real.
 *
 * The spinor states (2n per k point) go through the usual average
 * properties code with num_spin_comp = 2: the spin-up part of state J is
 * kept in avg_prop_info.orbs, the spin-down part in orbs_dn, and the charge
 * matrix chg_mat[J][mu] (2n x n) adds the two spin parts of orbital mu.
 * Since S_spinor = I2 (x) S does not couple the spins, Mulliken charges,
 * projected DOS and COOPs of a spinor state are exactly the sums of those of
 * its spin-up and spin-down parts.
 *
 * Input file (--soc FILE), one line per element and shell, ';' or '#' start
 * comments:
 *      Mo  d  0.080
 *      S   p  0.050
 *******************************************************/
#include "bind.h"

#define SOC_MAX_ENTRIES 512

typedef struct {
  char symb[ATOM_SYMB_LEN];
  int l;
  real lambda;
  int used;
} soc_entry_type;

typedef struct { double re, im; } zc_type;

static int soc_dim = 0;            /* 2*num_orbs */
static real *socR = 0, *socI = 0;  /* dense H_SOC (2n x 2n) */
static real *Hs = 0, *Ss = 0;      /* packed spinor H(k), S(k) */
static complex *zH = 0, *zS = 0, *zwork = 0;
static real *soc_vals = 0, *soc_rwork = 0;
static int soc_ham_file = -1, soc_ov_file = -1;
/* eigenvectors (state J, orbital mu): spin-up and spin-down parts, n2 x n */
static real *vupR = 0, *vupI = 0, *vdnR = 0, *vdnI = 0;
static real *soc_chg = 0;          /* charge matrix, n2 x n */
static real *soc_occ = 0;          /* occupations of the n2 states (one k point) */
static real *SdR = 0, *SdI = 0;    /* dense S(k), n x n */
static real *soc_tmp = 0;          /* work, n2 x n2 */
static real *soc_spin = 0;         /* spin-up minus spin-down weight of each AO, n2 x n */
/* SOC + spin polarization: second spin channel */
static int soc_two_tables = 0;
static real (*coul_dn)[4] = 0;     /* Hii (s,p,d,f) of every atom, spin down table */
static real const_dn = -1.0;       /* K of the spin-down channel */
static hermetian_matrix_type HdnR = {0,0}, HdnK = {0,0};

int soc_active(void) { return soc_file_name != 0; }
int soc_num_states(void) { return soc_dim; }

/********
  Lx, Ly, Lz (hbar = 1) in the real YAeHMOP orbitals of shell l.
  real_a = sum_m U[a][m] Y_{l,m}, m = -l..l, and L_real = conj(U) L U^T.
********/
static void real_orbital_U(int l, zc_type U[7][7])
{
  /* (m, kind) of each real orbital, kind 0: m=0, 1: cos-like, 2: sin-like */
  static const int p_m[3] = {1,1,0}, p_k[3] = {1,2,0};
  static const int d_m[5] = {2,0,2,1,1}, d_k[5] = {1,0,2,1,2};
  static const int f_m[7] = {0,1,1,2,2,3,3}, f_k[7] = {0,1,2,2,1,1,2};
  const int *mm, *kk;
  int a, m, sgn, n = 2*l+1;
  real r2 = 1.0/sqrt(2.0);

  if( l == 1 ){ mm = p_m; kk = p_k; }
  else if( l == 2 ){ mm = d_m; kk = d_k; }
  else { mm = f_m; kk = f_k; }
  for(a=0;a<7;a++) for(m=0;m<7;m++){ U[a][m].re = 0.0; U[a][m].im = 0.0; }
  for(a=0;a<n;a++){
    m = mm[a];
    sgn = (m % 2) ? -1 : 1;
    if( kk[a] == 0 ){
      U[a][l].re = 1.0;
    } else if( kk[a] == 1 ){   /* (1/sqrt2) [ (-1)^m Y_{l,m} + Y_{l,-m} ] */
      U[a][l+m].re = sgn*r2;
      U[a][l-m].re = r2;
    } else {                   /* (i/sqrt2) [ Y_{l,-m} - (-1)^m Y_{l,m} ] */
      U[a][l-m].im = r2;
      U[a][l+m].im = -sgn*r2;
    }
  }
}

static void angular_momentum(int l, zc_type Lx[7][7], zc_type Ly[7][7], zc_type Lz[7][7])
{
  zc_type U[7][7], Lm[3][7][7];
  int n = 2*l+1, a, b, i, j, c;
  real cp;

  real_orbital_U(l,U);
  /* Lx, Ly, Lz in |l,m>, m = -l..l (index i = m+l) */
  for(c=0;c<3;c++) for(i=0;i<7;i++) for(j=0;j<7;j++){
    Lm[c][i][j].re = 0.0; Lm[c][i][j].im = 0.0;
  }
  for(i=0;i<n;i++) Lm[2][i][i].re = (real)(i-l);
  for(i=0;i<n-1;i++){
    real mm = (real)(i-l);
    cp = sqrt(l*(l+1) - mm*(mm+1));      /* <m+1| L+ |m> */
    /* Lx = (L+ + L-)/2 ; Ly = (L+ - L-)/(2i) */
    Lm[0][i+1][i].re += cp/2;  Lm[0][i][i+1].re += cp/2;
    Lm[1][i+1][i].im += -cp/2; Lm[1][i][i+1].im += cp/2;
  }
  /* L_real[a][b] = sum_ij conj(U[a][i]) L[i][j] U[b][j] */
  for(c=0;c<3;c++){
    zc_type (*out)[7] = (c == 0) ? Lx : (c == 1) ? Ly : Lz;
    for(a=0;a<n;a++) for(b=0;b<n;b++){
      double sr = 0.0, si = 0.0;
      for(i=0;i<n;i++) for(j=0;j<n;j++){
        double ur = U[a][i].re, ui = -U[a][i].im;          /* conj(U[a][i]) */
        double lr = Lm[c][i][j].re, li = Lm[c][i][j].im;
        double vr = U[b][j].re, vi = U[b][j].im;
        double tr = ur*lr - ui*li, ti = ur*li + ui*lr;
        sr += tr*vr - ti*vi;
        si += tr*vi + ti*vr;
      }
      out[a][b].re = sr; out[a][b].im = si;
    }
  }
}

static int same_symbol(const char *a, const char *b)
{
  while( *a == ' ' ) a++;
  while( *b == ' ' ) b++;
  while( *a && *b && *a != ' ' && *b != ' ' ){
    char ca = *a, cb = *b;
    if( ca >= 'a' && ca <= 'z' ) ca += 'A'-'a';
    if( cb >= 'a' && cb <= 'z' ) cb += 'A'-'a';
    if( ca != cb ) return 0;
    a++; b++;
  }
  return (!*a || *a == ' ') && (!*b || *b == ' ');
}

static int read_soc_file(const char *name, soc_entry_type *entries)
{
  FILE *f = fopen(name,"r");
  char line[MAX_STR_LEN], symb[80], shell[80], err[MAX_STR_LEN+100];
  int n = 0, lineno = 0;
  double lam;

  if( !f ){
    sprintf(err,"Can't open the SOC file %s",name);
    fatal(err);
  }
  while( fgets(line,MAX_STR_LEN,f) ){
    char *p = line;
    lineno++;
    while( *p == ' ' || *p == '\t' ) p++;
    if( !*p || *p == '\n' || *p == ';' || *p == '#' || *p == '!' ) continue;
    if( sscanf(p,"%79s %79s %lf",symb,shell,&lam) != 3 ){
      sprintf(err,"SOC file %s, line %d: expected 'Element shell lambda'",name,lineno);
      fatal(err);
    }
    if( n >= SOC_MAX_ENTRIES ) fatal("Too many lines in the SOC file.");
    if( strlen(symb) >= ATOM_SYMB_LEN ) fatal("Element symbol too long in the SOC file.");
    strcpy(entries[n].symb,symb);
    switch( shell[0] ){
    case 'p': case 'P': entries[n].l = 1; break;
    case 'd': case 'D': entries[n].l = 2; break;
    case 'f': case 'F': entries[n].l = 3; break;
    case 's': case 'S':
      fprintf(stderr,"SOC file: s shell of %s ignored (L = 0)\n",symb);
      continue;
    default:
      sprintf(err,"SOC file %s, line %d: shell must be p, d or f",name,lineno);
      fatal(err);
    }
    entries[n].lambda = lam;
    entries[n].used = 0;
    n++;
  }
  fclose(f);
  return n;
}

/****************************************************************************
 *
 *                   Procedure soc_setup_two_tables
 *
 *   SOC on top of a spin-polarized calculation: H_spinor = diag(H_up(k),
 *   H_down(k)) + H_SOC.  The main parameter table (already in cell->atoms)
 *   is the spin-up one; this reads the Hii of the spin-down table
 *   (soc_parm_down).  Both tables must have the same basis and zetas, so that
 *   S_spinor = I2 (x) S still holds.  Special atoms ('*', parameters in the
 *   input file) are the same in both channels.
 *
 ****************************************************************************/
static void soc_setup_two_tables(cell_type *cell, detail_type *details, int num_orbs)
{
  atom_type *dn;
  int i, na = cell->num_atoms;
  char err[MAX_STR_LEN];
  real input_const = details->the_const;

  if( details->Execution_Mode == THIN )
    fatal("SOC with two spin tables is not implemented in THIN mode.");
  if( details->do_chg_it || details->do_muller_it )
    fatal("Charge iteration is not implemented for SOC with two spin tables.");

  dn = (atom_type *)calloc(na,sizeof(atom_type));
  coul_dn = calloc(na,sizeof(*coul_dn));
  if( !dn || !coul_dn ) fatal("Can't allocate the spin-down parameters.");
  memcpy(dn,cell->atoms,na*sizeof(atom_type));
  /* special atoms keep the parameters of the input in both channels */
  for(i=0;i<na;i++) if( dn[i].symb[0] == '*' ) dn[i].symb[0] = '&';
  fill_atomic_parms(dn,na,(FILE *)0,soc_parm_down);
  for(i=0;i<na;i++){
    atom_type *u = &(cell->atoms[i]), *d = &(dn[i]);
    if( u->symb[0] == '&' ) continue;
    if( u->symb[0] == '*' ){
      coul_dn[i][0] = u->coul_s; coul_dn[i][1] = u->coul_p;
      coul_dn[i][2] = u->coul_d; coul_dn[i][3] = u->coul_f;
      continue;
    }
    if( u->ns != d->ns || u->np != d->np || u->nd != d->nd || u->nf != d->nf ||
        fabs(u->exp_s - d->exp_s) > 1e-8 || fabs(u->exp_p - d->exp_p) > 1e-8 ||
        fabs(u->exp_d - d->exp_d) > 1e-8 || fabs(u->exp_d2 - d->exp_d2) > 1e-8 ||
        fabs(u->coeff_d1 - d->coeff_d1) > 1e-8 || fabs(u->coeff_d2 - d->coeff_d2) > 1e-8 ||
        fabs(u->exp_f - d->exp_f) > 1e-8 || fabs(u->exp_f2 - d->exp_f2) > 1e-8 ||
        fabs(u->coeff_f1 - d->coeff_f1) > 1e-8 || fabs(u->coeff_f2 - d->coeff_f2) > 1e-8 ){
      sprintf(err,"SOC with two spin tables: %s has a different basis or zetas in the "
              "spin-down table.\n\tBoth channels must share the orbitals (only Hii may differ).",
              u->symb);
      fatal(err);
    }
    coul_dn[i][0] = d->coul_s; coul_dn[i][1] = d->coul_p;
    coul_dn[i][2] = d->coul_d; coul_dn[i][3] = d->coul_f;
  }
  free(dn);

  if( soc_const_up > 0.0 ) details->the_const = soc_const_up;
  const_dn = soc_const_down > 0.0 ? soc_const_down : input_const;

  HdnR.dim = HdnK.dim = num_orbs;
  HdnR.mat = (real *)calloc(num_orbs*num_orbs,sizeof(real));
  HdnK.mat = (real *)calloc(num_orbs*num_orbs,sizeof(real));
  if( !HdnR.mat || !HdnK.mat ) fatal("Can't allocate the spin-down hamiltonian.");
  soc_two_tables = 1;

  fprintf(output_file,"\n; SOC on top of two spin channels: H = diag(H_up, H_down) + H_SOC\n");
  fprintf(output_file,"#SOC_SPIN_TABLES: up (main table)  down %s\n",soc_parm_down);
  fprintf(output_file,"#THE_CONST_UP: %lf\n#THE_CONST_DOWN: %lf\n",details->the_const,const_dn);
  fprintf(output_file,"; atom   Hii up (s p d f)   /   Hii down (s p d f)\n");
  for(i=0;i<na;i++){
    atom_type *u = &(cell->atoms[i]);
    if( u->symb[0] == '&' ) continue;
    fprintf(output_file,"; %s%-3d %8.4lf %8.4lf %8.4lf %8.4lf / %8.4lf %8.4lf %8.4lf %8.4lf\n",
            u->symb,i+1,u->coul_s,u->coul_p,u->coul_d,u->coul_f,
            coul_dn[i][0],coul_dn[i][1],coul_dn[i][2],coul_dn[i][3]);
  }
}

/****************************************************************************
 *
 *                   Procedure soc_build_down_hamiltonian
 *
 *   SOC with two spin tables: the R space hamiltonian of the spin-down
 *   channel (same call as for spin up, with the spin-down Hii and K).
 *   Called from driver.c right after the spin-up one is built.
 *
 ****************************************************************************/
void soc_build_down_hamiltonian(cell_type *cell, detail_type *details,
                                hermetian_matrix_type overlapR, int num_orbs,
                                int *orbital_lookup_table)
{
  int i, na;
  real save_const, (*save)[4];
  if( !soc_two_tables ) return;
  na = cell->num_atoms;
  save = calloc(na,sizeof(*save));
  if( !save ) fatal("Can't build the spin-down hamiltonian.");
  for(i=0;i<na;i++){
    atom_type *a = &(cell->atoms[i]);
    save[i][0] = a->coul_s; save[i][1] = a->coul_p; save[i][2] = a->coul_d; save[i][3] = a->coul_f;
    if( a->symb[0] == '&' ) continue;
    a->coul_s = coul_dn[i][0]; a->coul_p = coul_dn[i][1];
    a->coul_d = coul_dn[i][2]; a->coul_f = coul_dn[i][3];
  }
  save_const = details->the_const;
  details->the_const = const_dn;
  full_R_space_Hamiltonian(cell,details,overlapR,HdnR,num_orbs,orbital_lookup_table,0);
  details->the_const = save_const;
  for(i=0;i<na;i++){
    atom_type *a = &(cell->atoms[i]);
    a->coul_s = save[i][0]; a->coul_p = save[i][1]; a->coul_d = save[i][2]; a->coul_f = save[i][3];
  }
  free(save);
}

/****************************************************************************
 *
 *                   Procedure soc_setup
 *
 *   Reads the SOC file, builds the dense on-site H_SOC (2n x 2n), allocates
 *   the spinor work space and refuses the options not done with SOC yet.
 *   Called once, after read_inputfile.
 *
 ****************************************************************************/
void soc_setup(cell_type *cell, detail_type *details, int num_orbs,
               int *orbital_lookup_table)
{
  soc_entry_type entries[SOC_MAX_ENTRIES];
  int num_entries, i, e, a, b, l, n2, off, lwork;
  zc_type Lx[7][7], Ly[7][7], Lz[7][7];
  char err[MAX_STR_LEN];

  if( !soc_active() ) return;
#ifndef USE_LAPACK
  fatal("SOC needs the program built with LAPACK (-DUSE_BLAS_LAPACK=ON).");
#endif

  /******** what is not done with SOC yet (see docs/SPIN_SOC.md) ********/
#define SOC_REFUSE(cond,what) if( cond ){ \
    sprintf(err,"%s is not implemented with SOC yet.",what); fatal(err); }
  SOC_REFUSE(details->num_FMO_frags || details->num_FCO_frags,"FMO/FCO");
  {
    /* Walsh diagrams work; the items that need single MO coefficients don't */
    printing_info_type *p_info = details->step_print_options;
    while( p_info ){
      SOC_REFUSE(p_info->which_to_print == PRT_WAVE_FUNC ||
                 p_info->which_to_print == PRT_ORB_COEFF,"Walsh printing of MO coefficients");
      SOC_REFUSE(p_info->which_to_print == PRT_RCHG_MAT,"Walsh printing of the reduced charge matrix");
      p_info = p_info->next;
    }
  }
  SOC_REFUSE(details->num_MOs_to_print,"MO printing for viewkel (Print MOs)");
  SOC_REFUSE(details->vary_zeta,"Zeta variation");
  SOC_REFUSE(details->mod_OP_mat_PRT || details->mod_ROP_mat_PRT || details->mod_net_chg_PRT,
             "Modified Mulliken analysis");
  SOC_REFUSE(details->dump_sparse_mats,"Sparse matrix dump");
  SOC_REFUSE(details->Execution_Mode == MOLECULAR && details->eval_electrostat,
             "Electrostatic energy");
#ifdef INCLUDE_NETCDF_SUPPORT
  SOC_REFUSE(details->do_netCDF,"netCDF output");
#endif
#undef SOC_REFUSE

  if( soc_parm_down ) soc_setup_two_tables(cell,details,num_orbs);

  num_entries = read_soc_file(soc_file_name,entries);
  max_occup = 1.0;
  num_spin_comp = 2;
  soc_dim = n2 = 2*num_orbs;

  socR = (real *)calloc(n2*n2,sizeof(real));
  socI = (real *)calloc(n2*n2,sizeof(real));
  Hs = (real *)calloc(n2*n2,sizeof(real));
  Ss = (real *)calloc(n2*n2,sizeof(real));
  zH = (complex *)calloc(n2*n2,sizeof(complex));
  zS = (complex *)calloc(n2*n2,sizeof(complex));
  lwork = n2*n2;
  zwork = (complex *)calloc(lwork,sizeof(complex));
  soc_vals = (real *)calloc(n2,sizeof(real));
  soc_rwork = (real *)calloc(3*n2,sizeof(real));
  vupR = (real *)calloc(n2*num_orbs,sizeof(real));
  vupI = (real *)calloc(n2*num_orbs,sizeof(real));
  vdnR = (real *)calloc(n2*num_orbs,sizeof(real));
  vdnI = (real *)calloc(n2*num_orbs,sizeof(real));
  soc_chg = (real *)calloc(n2*num_orbs,sizeof(real));
  soc_occ = (real *)calloc(n2,sizeof(real));
  SdR = (real *)calloc(num_orbs*num_orbs,sizeof(real));
  SdI = (real *)calloc(num_orbs*num_orbs,sizeof(real));
  soc_tmp = (real *)calloc(n2*n2,sizeof(real));
  soc_spin = (real *)calloc(n2*num_orbs,sizeof(real));
  if( !socR || !socI || !Hs || !Ss || !zH || !zS || !zwork || !soc_vals || !soc_rwork ||
      !vupR || !vupI || !vdnR || !vdnI || !soc_chg || !soc_occ || !SdR || !SdI || !soc_tmp || !soc_spin )
    fatal("Can't allocate memory for the SOC matrices.");

  fprintf(output_file,"\n; Spin-orbit coupling (on-site lambda L.S), file %s\n",soc_file_name);
  fprintf(output_file,"#SOC_STATES: %d\n",n2);
  fprintf(output_file,"; spinor basis: orbitals 1..%d spin up, %d..%d the same orbitals spin down\n",
          num_orbs,num_orbs+1,n2);
  fprintf(output_file,"; atom   shell  lambda (eV)\n");

  for(i=0;i<cell->num_atoms;i++){
    atom_type *atom = &(cell->atoms[i]);
    int begin = orbital_lookup_table[i], shell_start[4];
    if( begin < 0 ) continue;              /* dummy atom */
    off = begin;
    shell_start[0] = shell_start[1] = shell_start[2] = shell_start[3] = -1;
    if( atom->ns ) off += 1;
    if( atom->np ){ shell_start[1] = off; off += 3; }
    if( atom->nd ){ shell_start[2] = off; off += 5; }
    if( atom->nf ){ shell_start[3] = off; off += 7; }
    for(e=0;e<num_entries;e++){
      if( !same_symbol(entries[e].symb,atom->symb) ) continue;
      l = entries[e].l;
      if( shell_start[l] < 0 ){
        sprintf(err,"SOC file: %s has no %c shell in this basis; ignored.",
                atom->symb,"spdf"[l]);
        if( !entries[e].used ) error(err);
        entries[e].used = 1;
        continue;
      }
      entries[e].used = 1;
      if( entries[e].lambda == 0.0 ) continue;
      fprintf(output_file,"%4s%-4d  %c   %10.6lf\n",atom->symb,i+1,"spdf"[l],
              entries[e].lambda);
      angular_momentum(l,Lx,Ly,Lz);
      {
        real h = 0.5*entries[e].lambda;
        int s0 = shell_start[l], nl = 2*l+1;
        for(a=0;a<nl;a++) for(b=0;b<nl;b++){
          int up_a = s0+a, up_b = s0+b;
          int dn_a = num_orbs+s0+a, dn_b = num_orbs+s0+b;
          /* up-up: +h Lz ; dn-dn: -h Lz */
          socR[up_a*n2+up_b] += h*Lz[a][b].re;  socI[up_a*n2+up_b] += h*Lz[a][b].im;
          socR[dn_a*n2+dn_b] -= h*Lz[a][b].re;  socI[dn_a*n2+dn_b] -= h*Lz[a][b].im;
          /* up-dn: h (Lx - i Ly) ; dn-up: h (Lx + i Ly) */
          socR[up_a*n2+dn_b] += h*(Lx[a][b].re + Ly[a][b].im);
          socI[up_a*n2+dn_b] += h*(Lx[a][b].im - Ly[a][b].re);
          socR[dn_a*n2+up_b] += h*(Lx[a][b].re - Ly[a][b].im);
          socI[dn_a*n2+up_b] += h*(Lx[a][b].im + Ly[a][b].re);
        }
      }
    }
  }
  for(e=0;e<num_entries;e++){
    if( !entries[e].used ){
      sprintf(err,"SOC file: element %s is not in the structure; ignored.",entries[e].symb);
      error(err);
    }
  }

  /* hermiticity check of H_SOC */
  for(a=0;a<n2;a++) for(b=0;b<n2;b++){
    if( fabs(socR[a*n2+b]-socR[b*n2+a]) > 1e-10 || fabs(socI[a*n2+b]+socI[b*n2+a]) > 1e-10 )
      FATAL_BUG("H_SOC is not hermitian.");
  }

  /* text copy of H_SOC, to rebuild the spinor H(R) outside the program */
  {
    char name[600];
    FILE *f;
    sprintf(name,"%s.SOC",details->filename);
    f = fopen(name,"w");
    if( !f ) fatal("Can't open the .SOC file.");
    fprintf(f,"# on-site spin-orbit matrix H_SOC (eV) from %s\n",soc_file_name);
    fprintf(f,"# spinor basis: 1..%d spin up, %d..%d spin down (YAeHMOP orbital order)\n",
            num_orbs,num_orbs+1,n2);
    fprintf(f,"# H_spinor(k) = I2 (x) H(k) + H_SOC,  S_spinor(k) = I2 (x) S(k)\n");
    fprintf(f,"# nonzero elements: i j Re Im (1-based)\n");
    fprintf(f,"%d\n",n2);
    for(a=0;a<n2;a++) for(b=0;b<n2;b++){
      if( fabs(socR[a*n2+b]) > 1e-12 || fabs(socI[a*n2+b]) > 1e-12 )
        fprintf(f,"%d %d %.12f %.12f\n",a+1,b+1,socR[a*n2+b],socI[a*n2+b]);
    }
    fclose(f);
  }

  if( details->Execution_Mode == MOLECULAR && details->use_symmetry ){
    fprintf(output_file,"; SOC: the symmetry labels of the MOs are not determined (double group)\n");
  }
}

/* physical element (i,j) of a packed hermitian matrix of dimension n */
static void packed_get(const real *mat, int n, int i, int j, real *re, real *im)
{
  if( i == j ){ *re = mat[i*n+i]; *im = 0.0; }
  else if( i < j ){ *re = mat[i*n+j]; *im = -mat[j*n+i]; }
  else { *re = mat[j*n+i]; *im = mat[i*n+j]; }
}

/* packed spinor H(k), S(k) (dimension 2n) from the spinless ones */
static void build_spinor(hermetian_matrix_type hamilK, hermetian_matrix_type overlapK,
                         int num_orbs)
{
  int n2 = soc_dim, I, J;
  /* two spin tables: spin-down H(k) from its R space hamiltonian */
  if( soc_two_tables ){
    if( details->Execution_Mode == MOLECULAR )
      memcpy(HdnK.mat,HdnR.mat,num_orbs*num_orbs*sizeof(real));
    else
      build_k_hamil_FAT(unit_cell,HdnR,HdnK,overlapK,num_orbs);
  }
  for(I=0;I<n2;I++){
    for(J=I;J<n2;J++){
      int a = I/num_orbs, b = J/num_orbs, p = I%num_orbs, q = J%num_orbs;
      real hr = 0.0, hi = 0.0, sr = 0.0, si = 0.0;
      if( a == b ){
        if( a == 1 && soc_two_tables ) packed_get(HdnK.mat,num_orbs,p,q,&hr,&hi);
        else packed_get(hamilK.mat,num_orbs,p,q,&hr,&hi);
        packed_get(overlapK.mat,num_orbs,p,q,&sr,&si);
      }
      hr += socR[I*n2+J];
      hi += socI[I*n2+J];
      if( I == J ){
        Hs[I*n2+I] = hr;
        Ss[I*n2+I] = sr;
      } else {
        Hs[I*n2+J] = hr;  Hs[J*n2+I] = -hi;
        Ss[I*n2+J] = sr;  Ss[J*n2+I] = -si;
      }
    }
  }
}

/* eigenvalues of the packed spinor pair (Hs, Ss) into soc_vals; with
   want_vectors the eigenvectors are left in zH (column J = state J) */
static void diag_spinor(detail_type *details, int want_vectors)
{
#ifdef USE_LAPACK
  int n2 = soc_dim, j, k, itype = 1, info = 0, lwork = n2*n2;
  char jobz = want_vectors ? 'V' : 'N', uplo = 'L';
  for(j=0;j<n2;j++){
    int jtab = j*n2;
    for(k=j+1;k<n2;k++){
      int ktab = k*n2;
      zH[jtab+k].r = Hs[jtab+k];  zH[jtab+k].i = Hs[ktab+j];
      zS[jtab+k].r = Ss[jtab+k];  zS[jtab+k].i = Ss[ktab+j];
      zH[ktab+j].r = zH[ktab+j].i = 0.0;
      zS[ktab+j].r = zS[ktab+j].i = 0.0;
    }
    zH[jtab+j].r = Hs[jtab+j];  zH[jtab+j].i = 0.0;
    zS[jtab+j].r = Ss[jtab+j];  zS[jtab+j].i = 0.0;
  }
  if( !details->diag_wo_overlap ){
    zhegv((integer *)&itype,&jobz,&uplo,(integer *)&n2,zH,(integer *)&n2,zS,
          (integer *)&n2,soc_vals,zwork,(integer *)&lwork,soc_rwork,(integer *)&info);
  } else {
    zheev(&jobz,&uplo,(integer *)&n2,zH,(integer *)&n2,soc_vals,zwork,
          (integer *)&lwork,soc_rwork,(integer *)&info);
  }
  fprintf(status_file,"Error value from SOC Diagonalization (0 is good): %d\n",info);
  if( info != 0 ) error("Problems in the SOC diagonalization, try more overlaps.");
#endif
}

/* split the eigenvectors in zH into spin-up and spin-down parts */
static void extract_vectors(int num_orbs)
{
  int n2 = soc_dim, J, mu;
  for(J=0;J<n2;J++){
    complex *col = &zH[J*n2];
    for(mu=0;mu<num_orbs;mu++){
      vupR[J*num_orbs+mu] = col[mu].r;
      vupI[J*num_orbs+mu] = col[mu].i;
      vdnR[J*num_orbs+mu] = col[num_orbs+mu].r;
      vdnI[J*num_orbs+mu] = col[num_orbs+mu].i;
    }
  }
}

/* charge matrix of the spinor states:
     chg[J][mu] = sum over spins of Re( conj(c_mu) (S c)_mu )
   (S_spinor = I2 (x) S, so the spins don't mix); each row adds up to 1. */
static void spinor_charge_matrix(hermetian_matrix_type overlapK, int num_orbs)
{
  int n2 = soc_dim, J, mu, nu, sp;
  for(mu=0;mu<num_orbs;mu++) for(nu=0;nu<num_orbs;nu++)
    packed_get(overlapK.mat,num_orbs,mu,nu,&SdR[mu*num_orbs+nu],&SdI[mu*num_orbs+nu]);
  for(J=0;J<n2;J++){
    for(mu=0;mu<num_orbs;mu++) soc_chg[J*num_orbs+mu] = 0.0;
    for(sp=0;sp<2;sp++){
      real *vr = sp ? &vdnR[J*num_orbs] : &vupR[J*num_orbs];
      real *vi = sp ? &vdnI[J*num_orbs] : &vupI[J*num_orbs];
      for(mu=0;mu<num_orbs;mu++){
        real scr = 0.0, sci = 0.0;
        for(nu=0;nu<num_orbs;nu++){
          real sr = SdR[mu*num_orbs+nu], si = SdI[mu*num_orbs+nu];
          scr += sr*vr[nu] - si*vi[nu];
          sci += sr*vi[nu] + si*vr[nu];
        }
        soc_chg[J*num_orbs+mu] += vr[mu]*scr + vi[mu]*sci;
        if( sp == 0 ) soc_spin[J*num_orbs+mu] = vr[mu]*scr + vi[mu]*sci;
        else soc_spin[J*num_orbs+mu] -= vr[mu]*scr + vi[mu]*sci;
      }
    }
  }
}

/* dense real and imaginary parts of a packed matrix of dimension n */
static void print_spinor_matrix(const char *what, real *packed, int n, detail_type *details)
{
  int i, j;
  real re, im;
  for(i=0;i<n;i++) for(j=0;j<n;j++){
    packed_get(packed,n,i,j,&re,&im);
    soc_tmp[i*n+j] = re;
  }
  fprintf(output_file,";\t\t --- %s, real part (spinor basis: 1..%d up, %d..%d down) ---\n",
          what,n/2,n/2+1,n);
  printmat(soc_tmp,n,n,output_file,1e-4,0,details->line_width);
  for(i=0;i<n;i++) for(j=0;j<n;j++){
    packed_get(packed,n,i,j,&re,&im);
    soc_tmp[i*n+j] = im;
  }
  fprintf(output_file,";\t\t --- %s, imaginary part ---\n",what);
  printmat(soc_tmp,n,n,output_file,1e-4,0,details->line_width);
}

/****************************************************************************
 *
 *                   Procedure soc_postprocess
 *
 *   What postprocess_results does at one k point (or for a molecule), for
 *   the spinor states: wave functions, energies and occupations, Mulliken
 *   overlap populations and net charges, charge matrix.
 *
 ****************************************************************************/
static void soc_postprocess(cell_type *cell, detail_type *details,
                            hermetian_matrix_type overlapK, int num_orbs)
{
  int n2 = soc_dim, i, j, J, a, begin, end;
  eigenset_type tmpset;
  real tot, accumR;

  if( details->wave_fn_PRT ){
    char transp = details->wave_fn_PRT & PRT_TRANSPOSE_FLAG;
    fprintf(output_file,";\t*** Spinor wavefunctions *** (states in rows, AO's in columns, "
            "energies INCREASE down)\n");
    fprintf(output_file,";\t***> SPIN UP, REAL:\n");
    print_labelled_mat(vupR,n2,num_orbs,output_file,1e-4,cell->atoms,cell->num_atoms,
                       orbital_lookup_table,num_orbs,transp,LABEL_COLS,details->line_width);
    fprintf(output_file,";\t***> SPIN UP, IMAGINARY:\n");
    print_labelled_mat(vupI,n2,num_orbs,output_file,1e-4,cell->atoms,cell->num_atoms,
                       orbital_lookup_table,num_orbs,transp,LABEL_COLS,details->line_width);
    fprintf(output_file,";\t***> SPIN DOWN, REAL:\n");
    print_labelled_mat(vdnR,n2,num_orbs,output_file,1e-4,cell->atoms,cell->num_atoms,
                       orbital_lookup_table,num_orbs,transp,LABEL_COLS,details->line_width);
    fprintf(output_file,";\t***> SPIN DOWN, IMAGINARY:\n");
    print_labelled_mat(vdnI,n2,num_orbs,output_file,1e-4,cell->atoms,cell->num_atoms,
                       orbital_lookup_table,num_orbs,transp,LABEL_COLS,details->line_width);
  }

  /* occupations of the 2n states (one electron each) */
  tmpset.dim = n2;
  tmpset.val = soc_vals;
  tmpset.vectR = tmpset.vectI = 0;
  for(J=0;J<n2;J++) soc_occ[J] = 0.0;
  calc_occupations(details,cell->num_electrons,n2,soc_occ,tmpset);

  if( details->levels_PRT || details->Execution_Mode == MOLECULAR ){
    real total = 0.0;
    fprintf(output_file,"\n#\t******* Energies (in eV)  and Occupation Numbers *******\n");
    fprintf(output_file,"; SOC: %d spinor states, one electron each\n",n2);
    for(J=0;J<n2;J++){
      fprintf(output_file,"%d:--->  %8.6lg  [%4.3lf Electrons]\n",J+1,soc_vals[J],soc_occ[J]);
      total += soc_occ[J]*soc_vals[J];
    }
    fprintf(output_file,"Total_Energy: %8.6lg\n",total);
    properties.total_E = total;
  }

  /* Mulliken analysis (same formulas as eval_mulliken, adding both spins) */
  if( details->OP_mat_PRT || details->ROP_mat_PRT || details->net_chg_PRT ){
    real *OP = properties.OP_mat, *acc = work1;
    for(i=0;i<num_orbs;i++) acc[i] = 0.0;
    for(i=0;i<num_orbs;i++){
      for(j=0;j<num_orbs;j++){
        real sr, si;
        accumR = 0.0;
        for(J=0;J<n2;J++){
          if( soc_occ[J] == 0.0 ) continue;
          accumR += soc_occ[J]*(vupR[J*num_orbs+i]*vupR[J*num_orbs+j] +
                                vupI[J*num_orbs+i]*vupI[J*num_orbs+j] +
                                vdnR[J*num_orbs+i]*vdnR[J*num_orbs+j] +
                                vdnI[J*num_orbs+i]*vdnI[J*num_orbs+j]);
        }
        packed_get(overlapK.mat,num_orbs,i,j,&sr,&si);
        OP[i*num_orbs+j] = (i != j ? 2.0 : 1.0)*accumR*sr;
        acc[i] += accumR*sr;
      }
    }
    for(a=0;a<cell->num_atoms;a++){
      find_atoms_orbs(num_orbs,cell->num_atoms,a,orbital_lookup_table,&begin,&end);
      if( begin >= 0 ){
        real net = VALENCE_FRACTION*cell->atoms[a].num_valence;
        for(i=begin;i<end;i++) net -= acc[i];
        properties.net_chgs[a] = net;
      }
    }
    if( details->ROP_mat_PRT ){
      reduced_mulliken(cell->num_atoms,num_orbs,orbital_lookup_table,
                       properties.OP_mat,properties.ROP_mat);
    }
  }
  if( details->OP_mat_PRT ){
    fprintf(output_file,
            "\n\n; \t\tq-q-q-q-q-q-q-q  Mulliken Analysis q-q-q-q-q-q-q-q\n");
    fprintf(output_file,"; SOC: populations add the spin-up and spin-down parts\n");
    fprintf(output_file,
            "\n;    Mulliken Overlap Population Matrix for %6.3lf electrons:\n",
            cell->num_electrons);
    print_labelled_mat(properties.OP_mat,num_orbs,num_orbs,output_file,1e-5,
                       cell->atoms,cell->num_atoms,orbital_lookup_table,
                       num_orbs,details->OP_mat_PRT & PRT_TRANSPOSE_FLAG,LABEL_BOTH,
                       details->line_width);
  }
  if( details->ROP_mat_PRT ){
    fprintf(output_file,
            "\n;   Reduced Mulliken Overlap Population Matrix for %6.3lf electrons:\n",
            cell->num_electrons);
    print_sym_mat(properties.ROP_mat,cell->num_atoms,cell->num_atoms,output_file,(char *)0,
                  (char *)0,details->line_width);
  }
  if( details->net_chg_PRT ){
    fprintf(output_file,"\n;    Net Atomic Charges for %6.3lf electrons:\n",
            cell->num_electrons);
    tot = 0.0;
    for(j=0;j<cell->num_atoms;j++){
      tot += properties.net_chgs[j];
      fprintf(output_file,"%d %s: %8.6lf\n", j+1, cell->atoms[j].symb,
              properties.net_chgs[j]);
    }
    fprintf(output_file,";      Total Charge is: %8.6lf\n",tot);

    /* spin moments (z): Mulliken spin-up minus spin-down populations */
    fprintf(output_file,"\n;    Spin moments m_z (N_up - N_down, Bohr magnetons) with SOC:\n");
    tot = 0.0;
    for(a=0;a<cell->num_atoms;a++){
      real m = 0.0;
      find_atoms_orbs(num_orbs,cell->num_atoms,a,orbital_lookup_table,&begin,&end);
      if( begin < 0 ) continue;
      for(J=0;J<n2;J++){
        if( soc_occ[J] == 0.0 ) continue;
        for(i=begin;i<end;i++) m += soc_occ[J]*soc_spin[J*num_orbs+i];
      }
      tot += m;
      fprintf(output_file,"%d %s: %8.6lf\n",a+1,cell->atoms[a].symb,m);
    }
    fprintf(output_file,";      Total moment is: %8.6lf\n",tot);
  }

  /* charge matrix: 2n spinor states x n orbitals */
  if( details->chg_mat_PRT || details->Rchg_mat_PRT ){
    fprintf(output_file,"\n\n; \t\tq-q-q-q-q-q-q-q  Charge Matrix q-q-q-q-q-q-q-q\n");
    fprintf(output_file,"; SOC: %d spinor states (rows), spin-up + spin-down weight of each AO\n",n2);
  }
  if( details->chg_mat_PRT ){
    fprintf(output_file,";  Complete Charge Matrix Independant of Occupation\n");
    print_labelled_mat(soc_chg,n2,num_orbs,output_file,1e-5,
                       cell->atoms,cell->num_atoms,orbital_lookup_table,
                       num_orbs,details->chg_mat_PRT & PRT_TRANSPOSE_FLAG,
                       LABEL_COLS,details->line_width);
  }
  if( details->Rchg_mat_PRT ){
    /* atoms in rows, states in columns, as reduced_charge_mat */
    for(a=0;a<cell->num_atoms;a++){
      find_atoms_orbs(num_orbs,cell->num_atoms,a,orbital_lookup_table,&begin,&end);
      for(J=0;J<n2;J++){
        real t = 0.0;
        if( begin >= 0 ) for(i=begin;i<end;i++) t += soc_chg[J*num_orbs+i];
        soc_tmp[a*n2+J] = t;
      }
    }
    fprintf(output_file,"\n;   Reduced Charge Matrix Independant of Occupation\n");
    printmat(soc_tmp,cell->num_atoms,n2,output_file,1e-5,(char)0,details->line_width);
  }
}

static int open_dump(const char *name, int nk, int n)
{
  int fd;
#ifndef _MSC_VER
  fd = open(name, O_RDWR | O_TRUNC | O_CREAT, S_IRUSR | S_IWUSR);
#else
  fd = open(name, O_RDWR | O_TRUNC | O_CREAT, _S_IREAD | _S_IWRITE);
#endif
  if( fd == -1 ) fatal("Can't open a SOC matrix dump file.");
  write(fd,(const char *)&nk,sizeof(int));
  write(fd,(const char *)&n,sizeof(int));
  return fd;
}

/****************************************************************************
 *
 *                   Procedure soc_kpoint
 *
 *   Work at one k point of loop_over_k_points (or for the molecule) when
 *   SOC is on: spinor matrices, printing and binary dumps (.HAM/.OV, same
 *   format as without SOC but with dimension 2n), diagonalization, the
 *   per-k properties (soc_postprocess) and the storage of the 2n states
 *   for the average properties (Fermi level, DOS, projected DOS, COOP,
 *   occupations and charges).
 *
 ****************************************************************************/
void soc_kpoint(cell_type *cell, detail_type *details, int ik, int num_KPOINTS,
                hermetian_matrix_type hamilK, hermetian_matrix_type overlapK,
                int num_orbs)
{
  char name[600];
  int j, mu, n2 = soc_dim, want_vectors;

  build_spinor(hamilK,overlapK,num_orbs);

  if( details->overlap_mat_PRT )
    print_spinor_matrix(details->Execution_Mode == MOLECULAR ? "Spinor overlap S(R)" :
                        "Spinor overlap S(K)",Ss,n2,details);
  if( details->hamil_PRT )
    print_spinor_matrix(details->Execution_Mode == MOLECULAR ? "Spinor hamiltonian H(R)" :
                        "Spinor hamiltonian H(K)",Hs,n2,details);

  if( details->dump_overlap ){
    if( ik == 0 ){
      sprintf(name,"%s.OV",details->filename);
      soc_ov_file = open_dump(name,num_KPOINTS,soc_dim);
    }
    dump_hermetian_mat(soc_ov_file,Ss,soc_dim);
    if( ik == num_KPOINTS-1 ){ close(soc_ov_file); soc_ov_file = -1; }
  }
  if( details->dump_hamil ){
    if( ik == 0 ){
      sprintf(name,"%s.HAM",details->filename);
      soc_ham_file = open_dump(name,num_KPOINTS,soc_dim);
    }
    dump_hermetian_mat(soc_ham_file,Hs,soc_dim);
    if( ik == num_KPOINTS-1 ){ close(soc_ham_file); soc_ham_file = -1; }
  }
  if( details->just_matrices ) return;

  want_vectors = !details->just_avgE;
  diag_spinor(details,want_vectors);

  if( want_vectors ){
    extract_vectors(num_orbs);
    spinor_charge_matrix(overlapK,num_orbs);
    soc_postprocess(cell,details,overlapK,num_orbs);
  }

  /* keep the states for the average properties */
  if( details->avg_props && avg_prop_info && details->Execution_Mode != MOLECULAR ){
    avg_prop_info_type *api = &(avg_prop_info[ik]);
    for(j=0;j<n2;j++) api->energies[j] = (float)soc_vals[j];
    if( want_vectors ){
      for(j=0;j<n2*num_orbs;j++){
        api->orbs[j] = (float)vupR[j];
        api->orbsI[j] = (float)vupI[j];
        api->orbs_dn[j] = (float)vdnR[j];
        api->orbsI_dn[j] = (float)vdnI[j];
        api->chg_mat[j] = (float)soc_chg[j];
        api->chg_spin[j] = (float)soc_spin[j];
      }
    }
  }
  (void)mu;
}

/* SOC eigenvalues of one band structure k point (from the spinless H(k),S(k)) */
real *soc_band_energies(detail_type *details, hermetian_matrix_type hamilK,
                        hermetian_matrix_type overlapK, int num_orbs)
{
  build_spinor(hamilK,overlapK,num_orbs);
  diag_spinor(details,0);
  return soc_vals;
}

/* energy of spinor state J and charge matrix element (state J, orbital mu)
   at the last k point diagonalized (or the molecule), for the Walsh output */
real soc_state_energy(int J)
{
  if( J < 0 || J >= soc_dim ) fatal("Walsh: SOC state number out of range.");
  return soc_vals[J];
}
real soc_chg_element(int J, int mu)
{
  int n = soc_dim/2;
  if( J < 0 || J >= soc_dim || mu < 0 || mu >= n )
    fatal("Walsh: SOC charge matrix element out of range.");
  return soc_chg[J*n+mu];
}

/****************************************************************************
 *
 *                   Procedure soc_print_spin_moments
 *
 *   Average spin moments m_z per atom with SOC (Mulliken populations of
 *   the spin-up minus the spin-down parts of the occupied states), with the
 *   same k point averaging as calc_avg_occups.
 *
 ****************************************************************************/
void soc_print_spin_moments(detail_type *details, cell_type *cell, int num_orbs,
                            K_orb_ptr_type *orbital_ordering, avg_prop_info_type *api)
{
  real *mom, tot_num_K = 0.0, num_occup_bands = 0.0, tot = 0.0, m;
  int i, j, a, begin, end;
  if( !soc_active() || details->just_avgE ) return;
  mom = (real *)calloc(num_orbs,sizeof(real));
  if( !mom ) fatal("Can't allocate the spin moments.");
  i = 0;
  while( orbital_ordering[i].occup > .0001 ){
    int k = orbital_ordering[i].Kpoint, MO = orbital_ordering[i].MO;
    float *w = &(api[k].chg_spin[MO*num_orbs]);
    for(j=0;j<num_orbs;j++)
      mom[j] += (real)w[j]*orbital_ordering[i].occup*details->K_POINTS[k].weight;
    i++;
  }
  for(i=0;i<details->num_KPOINTS;i++){
    tot_num_K += details->K_POINTS[i].weight*details->K_POINTS[i].num_filled_bands;
    num_occup_bands += details->K_POINTS[i].num_filled_bands;
  }
  fprintf(output_file,"\n# AVERAGE SPIN MOMENTS (SOC, m_z = N_up - N_down, Bohr magnetons)\n");
  for(a=0;a<cell->num_atoms;a++){
    find_atoms_orbs(num_orbs,cell->num_atoms,a,orbital_lookup_table,&begin,&end);
    if( begin < 0 ) continue;
    m = 0.0;
    for(j=begin;j<end;j++) m += mom[j];
    m *= num_occup_bands/(tot_num_K*details->num_KPOINTS);
    tot += m;
    fprintf(output_file,"%d %s %lf\n",a+1,cell->atoms[a].symb,m);
  }
  fprintf(output_file,"#TOTAL_SPIN_MOMENT: %lf\n",tot);
  free(mom);
}
