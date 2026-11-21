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
 * Spin-polarized runs: communication between the two channel
 * processes and the parent process.
 *
 * Each channel (spin up, spin down) runs in its own process (main.c).  The
 * parent keeps a pipe to and from each child:
 *
 *  - Fermi level shared by the channels (--common-fermi): after the energies
 *    of all k points are known and sorted, each child sends them to the
 *    parent (message SPIN_MSG_ENERGIES).  The parent merges the two lists,
 *    fills the lowest N*nk states (one electron each, N = electrons of the
 *    input) and answers each child with its number of electrons per cell.
 *    The child then fills its own states with that number, which gives the
 *    same occupations as a common Fermi level.
 *
 *  - Summary: at the end each child sends its Fermi level, electrons,
 *    average energy and net charges (SPIN_MSG_SUMMARY); the parent writes
 *    <input>.spin with the populations and magnetic moments.
 *
 * Messages are binary (the processes are on the same machine).
 *******************************************************/
#include "bind.h"
#ifndef _MSC_VER
#include <sys/wait.h>
#include <signal.h>
#endif

#define SPIN_MSG_ENERGIES 1
#define SPIN_MSG_SUMMARY  2

int spin_pipe_out = -1;      /* child: write end to the parent */
int spin_pipe_in = -1;       /* child: read end from the parent */
int spin_common_fermi = 0;   /* --common-fermi */
real spin_input_electrons = 0.0;  /* electrons of the input (both channels) */

static int write_all(int fd, const void *buf, size_t n)
{
  const char *p = (const char *)buf;
  while( n > 0 ){
    ssize_t w = write(fd,p,n);
    if( w <= 0 ) return -1;
    p += w; n -= (size_t)w;
  }
  return 0;
}

static int read_all(int fd, void *buf, size_t n)
{
  char *p = (char *)buf;
  while( n > 0 ){
    ssize_t r = read(fd,p,n);
    if( r <= 0 ) return -1;
    p += r; n -= (size_t)r;
  }
  return 0;
}

/****************************************************************************
 *
 *                   Procedure spin_exchange_fermi   (child)
 *
 *   Called after sort_avg_prop_info.  Sends the sorted energies to the
 *   parent and sets cell->num_electrons to the number of electrons of this
 *   channel for the common Fermi level.
 *
 ****************************************************************************/
void spin_exchange_fermi(detail_type *details, cell_type *cell, int num_orbs,
                         K_orb_ptr_type *orbital_ordering)
{
  int tag = SPIN_MSG_ENERGIES, nk = details->num_KPOINTS;
  int nst = num_orbs*num_spin_comp, tot = nst*nk, i;
  real answer[2];
  float *E;

  if( !spin_channel || !spin_common_fermi || spin_pipe_out < 0 ) return;
  E = (float *)malloc(tot*sizeof(float));
  if( !E ) fatal("Can't allocate the energies for the common Fermi level.");
  for(i=0;i<tot;i++) E[i] = *(orbital_ordering[i].energy);
  if( write_all(spin_pipe_out,&tag,sizeof(int)) || write_all(spin_pipe_out,&nk,sizeof(int)) ||
      write_all(spin_pipe_out,&nst,sizeof(int)) ||
      write_all(spin_pipe_out,&spin_input_electrons,sizeof(real)) ||
      write_all(spin_pipe_out,E,tot*sizeof(float)) )
    fatal("Common Fermi level: can't send the energies to the parent process.");
  free(E);
  if( read_all(spin_pipe_in,answer,2*sizeof(real)) )
    fatal("Common Fermi level: no answer from the parent process.");
  cell->num_electrons = answer[0];
  fprintf(output_file,"\n; Fermi level shared by both spin channels (--common-fermi):\n");
  fprintf(output_file,";   common Fermi level %lf eV, %lf electrons in this channel\n",
          answer[1],answer[0]);
  fprintf(output_file,"#ELECTRONS_THIS_CHANNEL_COMMON_FERMI: %lf\n",answer[0]);
}

/****************************************************************************
 *
 *                   Procedure spin_send_summary   (child)
 *
 ****************************************************************************/
void spin_send_summary(cell_type *cell, detail_type *details, prop_type *props)
{
  int tag = SPIN_MSG_SUMMARY, i, have = 0, na = cell->num_atoms;
  real vals[3];
  if( !spin_channel || spin_pipe_out < 0 ) return;
  if( details->avg_props && !details->just_avgE && props->net_chgs ) have = 1;
  vals[0] = props->Fermi_E;
  vals[1] = cell->num_electrons;
  vals[2] = props->total_E;
  if( write_all(spin_pipe_out,&tag,sizeof(int)) || write_all(spin_pipe_out,&have,sizeof(int)) ||
      write_all(spin_pipe_out,vals,3*sizeof(real)) || write_all(spin_pipe_out,&na,sizeof(int)) )
    goto fail;
  for(i=0;i<na;i++){
    char symb[ATOM_SYMB_LEN+1];
    real av[2];
    memset(symb,0,sizeof(symb));
    strncpy(symb,cell->atoms[i].symb,ATOM_SYMB_LEN);
    av[0] = cell->atoms[i].num_valence;
    av[1] = (have && cell->atoms[i].at_number >= 0) ? props->net_chgs[i] : 0.0;
    if( cell->atoms[i].symb[0] == '&' ) av[0] = -1.0;   /* dummy atom */
    if( write_all(spin_pipe_out,symb,sizeof(symb)) || write_all(spin_pipe_out,av,2*sizeof(real)) )
      goto fail;
  }
  return;
 fail:
  fprintf(stderr,"Spin channel: can't send the summary to the parent process.\n");
}

#ifndef _MSC_VER
typedef struct {
  int got, have, na;
  real fermi, nelec, energy;
  char (*symb)[ATOM_SYMB_LEN+1];
  real *val, *net;
} channel_summary_type;

typedef struct { float e; int ch; } spin_level_type;

static int cmp_level(const void *a, const void *b)
{
  float d = ((const spin_level_type *)a)->e - ((const spin_level_type *)b)->e;
  return d > 0 ? 1 : (d < 0 ? -1 : 0);
}

/* common Fermi level from the energies of both channels: answers ne[c] */
static void common_fermi(float *E[2], int tot[2], int nk, real nelec, real ne[2], real *ef)
{
  int n = tot[0]+tot[1], i, last, b, e, cnt[2] = {0,0}, gcnt[2] = {0,0};
  real left = nelec*(real)nk, last_part = 0.0, grp_el;
  spin_level_type *L = (spin_level_type *)malloc(n*sizeof(spin_level_type));
  if( !L ){ fprintf(stderr,"common Fermi: out of memory\n"); exit(-1); }
  for(i=0;i<tot[0];i++){ L[i].e = E[0][i]; L[i].ch = 0; }
  for(i=0;i<tot[1];i++){ L[tot[0]+i].e = E[1][i]; L[tot[0]+i].ch = 1; }
  qsort(L,n,sizeof(spin_level_type),cmp_level);
  /* fill one electron per state */
  for(i=0;i<n && left > 0.0;i++){
    if( left >= 1.0 ){ left -= 1.0; last_part = 1.0; }
    else { last_part = left; left = 0.0; }
  }
  last = i-1;
  if( last < 0 ){ ne[0] = ne[1] = 0.0; *ef = L[0].e; free(L); return; }
  /* states degenerate with the last occupied one share its electrons */
  b = last; while( b > 0 && fabs(L[b-1].e - L[last].e) < DEGEN_TOL ) b--;
  e = last; while( e < n-1 && fabs(L[e+1].e - L[last].e) < DEGEN_TOL ) e++;
  for(i=0;i<b;i++) cnt[L[i].ch]++;
  for(i=b;i<=e;i++) gcnt[L[i].ch]++;
  grp_el = (real)(last-b) + last_part;     /* electrons in the degenerate group */
  for(i=0;i<2;i++)
    ne[i] = ((real)cnt[i] + grp_el*(real)gcnt[i]/(real)(e-b+1))/(real)nk;
  *ef = L[last].e;
  free(L);
}

/****************************************************************************
 *
 *                   Procedure spin_parent_loop   (parent)
 *
 *   Serves the two channel processes until both finish, then writes the
 *   summary file.  fd_in[c]: from child c, fd_out[c]: to child c.
 *   Returns 0 if both channels sent their summary.
 *
 ****************************************************************************/
int spin_parent_loop(int fd_in[2], int fd_out[2], pid_t pids[2],
                     const char *summary_name, char out_names[2][520],
                     int common, real kconst[2])
{
  channel_summary_type S[2];
  int c, ok = 0, i;
  real last_ef = 0.0;
  memset(S,0,sizeof(S));

  while( 1 ){
    int tag[2];
    for(c=0;c<2;c++){
      if( read_all(fd_in[c],&tag[c],sizeof(int)) ) tag[c] = -1;
    }
    if( tag[0] < 0 || tag[1] < 0 ){
      /* a channel finished or died: if the other one waits for us, stop it */
      for(c=0;c<2;c++) if( tag[c] >= 0 && tag[c] == SPIN_MSG_ENERGIES ) kill(pids[c],SIGTERM);
      if( tag[0] != tag[1] ){
        for(c=0;c<2;c++) if( tag[c] < 0 ) kill(pids[1-c],SIGTERM);
      }
      break;
    }
    if( tag[0] != tag[1] ){
      fprintf(stderr,"Spin channels out of step (messages %d and %d).\n",tag[0],tag[1]);
      kill(pids[0],SIGTERM); kill(pids[1],SIGTERM);
      break;
    }
    if( tag[0] == SPIN_MSG_ENERGIES ){
      int nk[2], nst[2], tot[2];
      real nel[2], ne[2], ef, ans[2];
      float *E[2];
      for(c=0;c<2;c++){
        if( read_all(fd_in[c],&nk[c],sizeof(int)) || read_all(fd_in[c],&nst[c],sizeof(int)) ||
            read_all(fd_in[c],&nel[c],sizeof(real)) ) goto broken;
        tot[c] = nk[c]*nst[c];
        E[c] = (float *)malloc(tot[c]*sizeof(float));
        if( !E[c] || read_all(fd_in[c],E[c],tot[c]*sizeof(float)) ) goto broken;
      }
      if( nk[0] != nk[1] ) goto broken;
      common_fermi(E,tot,nk[0],nel[0],ne,&ef);
      last_ef = ef;
      for(c=0;c<2;c++){
        ans[0] = ne[c]; ans[1] = ef;
        if( write_all(fd_out[c],ans,2*sizeof(real)) ) goto broken;
        free(E[c]);
      }
    } else if( tag[0] == SPIN_MSG_SUMMARY ){
      for(c=0;c<2;c++){
        real vals[3];
        if( read_all(fd_in[c],&S[c].have,sizeof(int)) || read_all(fd_in[c],vals,3*sizeof(real)) ||
            read_all(fd_in[c],&S[c].na,sizeof(int)) ) goto broken;
        S[c].fermi = vals[0]; S[c].nelec = vals[1]; S[c].energy = vals[2];
        S[c].symb = calloc(S[c].na,sizeof(*S[c].symb));
        S[c].val = (real *)calloc(S[c].na,sizeof(real));
        S[c].net = (real *)calloc(S[c].na,sizeof(real));
        if( !S[c].symb || !S[c].val || !S[c].net ) goto broken;
        for(i=0;i<S[c].na;i++){
          real av[2];
          if( read_all(fd_in[c],S[c].symb[i],ATOM_SYMB_LEN+1) ||
              read_all(fd_in[c],av,2*sizeof(real)) ) goto broken;
          S[c].val[i] = av[0]; S[c].net[i] = av[1];
        }
        S[c].got = 1;
      }
    } else goto broken;
  }

  if( S[0].got && S[1].got ){
    FILE *f = fopen(summary_name,"w");
    real mom = S[0].nelec - S[1].nelec;
    ok = 1;
    if( !f ){ fprintf(stderr,"Can't open %s\n",summary_name); return 0; }
    fprintf(f,"# Spin-polarized run: summary of the two channels\n");
    fprintf(f,"# channels: %s.out (up), %s.out (down)\n",out_names[0],out_names[1]);
    fprintf(f,"#COMMON_FERMI: %s\n",common ? "yes" : "no (fixed electrons per channel)");
    if( common ) fprintf(f,"#Fermi_Energy:  %lf\n",last_ef);
    fprintf(f,"#Fermi_Energy_up:  %lf\n",S[0].fermi);
    fprintf(f,"#Fermi_Energy_down:  %lf\n",S[1].fermi);
    fprintf(f,"#ELECTRONS_UP:  %lf\n",S[0].nelec);
    fprintf(f,"#ELECTRONS_DOWN:  %lf\n",S[1].nelec);
    fprintf(f,"#MAGNETIC_MOMENT:  %lf   ; N_up - N_down per cell (Bohr magnetons)\n",mom);
    if( kconst[0] > 0 || kconst[1] > 0 )
      fprintf(f,"; K per channel: up %lf, down %lf (<= 0: THE CONST of the input)\n",
              kconst[0],kconst[1]);
    if( S[0].have && S[1].have && S[0].na == S[1].na ){
      real tu = 0, td = 0;
      fprintf(f,"#Average_Energy:  %lf   ; up + down (eV)\n",S[0].energy+S[1].energy);
      fprintf(f,"\n# ATOMIC POPULATIONS AND MOMENTS\n");
      fprintf(f,"; atom      pop_up     pop_down   net_charge   moment\n");
      for(i=0;i<S[0].na;i++){
        real pu, pd;
        if( S[0].val[i] < 0 ) continue;
        pu = 0.5*S[0].val[i] - S[0].net[i];
        pd = 0.5*S[1].val[i] - S[1].net[i];
        tu += pu; td += pd;
        fprintf(f,"%d %-4s %10.6lf %10.6lf %10.6lf %10.6lf\n",i+1,S[0].symb[i],pu,pd,
                S[0].net[i]+S[1].net[i],pu-pd);
      }
      fprintf(f,"; total     %10.6lf %10.6lf              %10.6lf\n",tu,td,tu-td);
    } else {
      fprintf(f,"; no average properties: populations not available\n");
    }
    fclose(f);
  }
  for(c=0;c<2;c++){ free(S[c].symb); free(S[c].val); free(S[c].net); }
  return ok ? 0 : -1;

 broken:
  fprintf(stderr,"Spin-polarized run: communication with the channels broke.\n");
  kill(pids[0],SIGTERM); kill(pids[1],SIGTERM);
  for(c=0;c<2;c++){ free(S[c].symb); free(S[c].val); free(S[c].net); }
  return -1;
}
#endif
