/*******************************************************

Copyright (C) 1995 Greg Landrum
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

/****************************************************************************
*
*     this is the main procedure for the program bind
*
*  read the file README.bind to find out more
*
*  created:  greg landrum  August 1993
*
*****************************************************************************/
/***
  Edit History:

  March '98: WG
    - print FMO-FMO COOP's to .out

  07.05.98 gL:
    added netCDF support.
   14.02.2004 gL:
     added a "10th anniversary" message :-)

***/
#include "bind.h"
#ifndef _MSC_VER
#include <sys/wait.h>
#endif

// For booleans
#include "stdbool.h"

FILE *COHP_file;
const char greetings[]="Welcome to the 10th Anniversary edition of YAeHMOP!\n";

/* outputs of a spin channel: the suffix goes before the extension of the
   input name (MoS2.bind -> MoS2_up.bind, so MoS2_up.bind.out, .band, ...) */
static void channel_name(const char *in, const char *suffix, char *out, int len)
{
  const char *slash = strrchr(in,'/');
  const char *base = slash ? slash+1 : in;
  const char *dot = strrchr(base,'.');
  if( dot && dot != base ){
    int n = (int)(dot - in);
    snprintf(out,len,"%.*s%s%s",n,in,suffix,dot);
  } else {
    snprintf(out,len,"%s%s",in,suffix);
  }
}

static void usage(void)
{
  fprintf(stderr,
    "Usage: bind [options] <inputfile> [paramfile]\n"
    "       bind --use_stdin_stdout [paramfile]\n"
    "Options:\n"
    "  --parm-up FILE --parm-down FILE  spin-polarized run: one parameter table\n"
    "                    per spin channel; outputs get _up / _down before the\n"
    "                    input extension (one electron per state in each channel)\n"
    "  --nelec-up X --nelec-down X      electrons in each channel (default: half\n"
    "                    of the electrons of the input in each one)\n"
    "  --const-up K --const-down K      Wolfsberg-Helmholz K of each channel\n"
    "                    (default: THE CONST of the input)\n"
    "  --common-fermi    spin-polarized run with one Fermi level for both\n"
    "                    channels (N_up and N_down come out of it); summary\n"
    "                    of both channels in <inputfile>.spin\n"
    "  --soc FILE        spin-orbit coupling, lines 'Element shell lambda(eV)';\n"
    "                    with --parm-up/--parm-down: SOC on top of the two spin\n"
    "                    channels, H = diag(H_up, H_down) + H_SOC (same zetas)\n"
    "  -v                print the version\n");
}

int main(int argc, char **argv)
{
  char file_name[500];
  bool use_stdin_stdout = false;
  char *positional[2] = {0,0};
  int num_positional = 0;
  char *parm_up = 0, *parm_down = 0;
  real nelec_up = -1.0, nelec_down = -1.0;
  real const_up = -1.0, const_down = -1.0;
  int i;

  if( argc == 2 && strcmp(argv[1], "-v") == 0){
    fprintf(stdout, "version: %s\n", VERSION_STRING);
    exit(0);
  }

  for(i=1;i<argc;i++){
    if( !strcmp(argv[i],"--use_stdin_stdout") ){
      use_stdin_stdout = true;
    } else if( !strcmp(argv[i],"--parm-up") && i+1 < argc ){
      parm_up = argv[++i];
    } else if( !strcmp(argv[i],"--parm-down") && i+1 < argc ){
      parm_down = argv[++i];
    } else if( !strcmp(argv[i],"--nelec-up") && i+1 < argc ){
      nelec_up = atof(argv[++i]);
    } else if( !strcmp(argv[i],"--nelec-down") && i+1 < argc ){
      nelec_down = atof(argv[++i]);
    } else if( !strcmp(argv[i],"--const-up") && i+1 < argc ){
      const_up = atof(argv[++i]);
    } else if( !strcmp(argv[i],"--const-down") && i+1 < argc ){
      const_down = atof(argv[++i]);
    } else if( !strcmp(argv[i],"--common-fermi") ){
      spin_common_fermi = 1;
    } else if( !strcmp(argv[i],"--soc") && i+1 < argc ){
      soc_file_name = argv[++i];
    } else if( !strcmp(argv[i],"-h") || !strcmp(argv[i],"--help") ){
      usage();
      exit(0);
    } else if( argv[i][0] == '-' && argv[i][1] == '-' ){
      fprintf(stderr,"Unknown or incomplete option: %s\n",argv[i]);
      usage();
      exit(-1);
    } else if( num_positional < 2 ){
      positional[num_positional++] = argv[i];
    } else {
      usage();
      exit(-1);
    }
  }

  /* with stdin the only positional argument is the parameter file */
  char* parm_file_name = NULL;
  if( use_stdin_stdout ){
    strcpy(file_name,"stdin");
    if( num_positional > 1 ){ usage(); exit(-1); }
    if( num_positional == 1 ) parm_file_name = positional[0];
  } else {
    if( num_positional < 1 ){ usage(); exit(-1); }
    strncpy(file_name,positional[0],499);
    file_name[499] = 0;
    if( num_positional == 2 ) parm_file_name = positional[1];
  }

  /* install the sig_int handler */
  signal(SIGINT,handle_sigint);

  fprintf(stderr,greetings);

  if( !parm_up && !parm_down ){
    if( nelec_up >= 0 || nelec_down >= 0 || const_up > 0 || const_down > 0 ||
        spin_common_fermi ){
      fprintf(stderr,"--nelec-*, --const-* and --common-fermi need --parm-up/--parm-down.\n");
      exit(-1);
    }
    run_bind(file_name, use_stdin_stdout, parm_file_name);
    exit(0);
  }

  /********
    spin-polarized run: each channel is an independent EHT calculation with
    its own parameter table, done in its own process (the program keeps a
    lot of global and static state, so the channels must not share one).
  ********/
  if( !parm_up || !parm_down ){
    fprintf(stderr,"Spin-polarized runs need both --parm-up and --parm-down.\n");
    exit(-1);
  }
  if( use_stdin_stdout ){
    fprintf(stderr,"Spin-polarized runs can't use --use_stdin_stdout.\n");
    exit(-1);
  }
  if( parm_file_name ){
    fprintf(stderr,"Give the tables with --parm-up/--parm-down, not as a second argument.\n");
    exit(-1);
  }
  if( spin_common_fermi && (nelec_up >= 0 || nelec_down >= 0) ){
    fprintf(stderr,"--common-fermi and --nelec-up/--nelec-down exclude each other.\n");
    exit(-1);
  }
  if( soc_file_name ){
    /********
      SOC on top of the two spin channels: one calculation in one process
      with the spinor hamiltonian diag(H_up, H_down) + H_SOC (soc.c).  The
      up table is the main parameter table; soc.c reads the down one.
    ********/
    if( nelec_up >= 0 || nelec_down >= 0 || spin_common_fermi ){
      fprintf(stderr,"With --soc the spin channels are coupled: one Fermi level and\n"
              "the electrons of the input; --nelec-* and --common-fermi don't apply.\n");
      exit(-1);
    }
    soc_parm_down = parm_down;
    soc_const_up = const_up;
    soc_const_down = const_down;
    run_bind(file_name, use_stdin_stdout, parm_up);
    exit(0);
  }
#ifdef _MSC_VER
  fprintf(stderr,"Spin-polarized runs are not available on Windows.\n");
  exit(-1);
#else
  {
    const char *suffix[2] = {"_up","_down"};
    char *parms[2];
    real nelec[2], kconst[2];
    char out_name[2][520], summary_name[520];
    pid_t pids[2];
    int c, status, failed = 0;
    int to_parent[2][2], to_child[2][2], fd_in[2], fd_out[2];
    parms[0] = parm_up;   parms[1] = parm_down;
    nelec[0] = nelec_up;  nelec[1] = nelec_down;
    kconst[0] = const_up; kconst[1] = const_down;
    for(c=0;c<2;c++){
      channel_name(file_name,suffix[c],out_name[c],520);
      if( pipe(to_parent[c]) || pipe(to_child[c]) ){
        perror("pipe");
        exit(-1);
      }
    }
    snprintf(summary_name,520,"%s.spin",file_name);
    signal(SIGPIPE,SIG_IGN);
    fflush(stdout); fflush(stderr);
    for(c=0;c<2;c++){
      pids[c] = fork();
      if( pids[c] < 0 ){
        perror("fork");
        exit(-1);
      }
      if( pids[c] == 0 ){
        int o;
        /* keep only this channel's ends of the pipes */
        for(o=0;o<2;o++){
          close(to_parent[o][0]);
          close(to_child[o][1]);
          if( o != c ){ close(to_parent[o][1]); close(to_child[o][0]); }
        }
        spin_pipe_out = to_parent[c][1];
        spin_pipe_in = to_child[c][0];
        spin_channel = c+1;
        spin_nelec = nelec[c];
        spin_const = kconst[c];
        run_bind_ex(file_name, out_name[c], false, parms[c]);
        close(spin_pipe_out);
        exit(0);
      }
    }
    for(c=0;c<2;c++){
      close(to_parent[c][1]);
      close(to_child[c][0]);
      fd_in[c] = to_parent[c][0];
      fd_out[c] = to_child[c][1];
    }
    if( spin_parent_loop(fd_in,fd_out,pids,summary_name,out_name,
                         spin_common_fermi,kconst) ) failed = 1;
    for(c=0;c<2;c++){
      if( waitpid(pids[c],&status,0) < 0 || !WIFEXITED(status) ||
          WEXITSTATUS(status) != 0 ){
        fprintf(stderr,"Spin channel %s failed (see %s.status).\n",
                suffix[c]+1,out_name[c]);
        failed = 1;
      }
    }
    if( failed ) exit(-1);
    fprintf(stderr,"Spin-polarized run done: %s.* and %s.*, summary in %s\n",
            out_name[0],out_name[1],summary_name);
  }
#endif
  exit(0);
}
