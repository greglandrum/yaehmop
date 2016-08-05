/************************************************************************

  This program reads the values for a hermetian matrix out of a
   binary file, and then generates a .ps file which shows in a
   block way.

   Created by greg Landrum March 1994
************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <fcntl.h>

#include "fit_props.h"

/* the width and height in points */
#define BOX_WIDTH 432
#define BOX_HEIGHT 432

#define START_X 50
#define START_Y 100



#define HERMETIAN_R(matrix,i,j)  ( j>i ? matrix[i*num_orbs+j] : matrix[j*num_orbs+i] )
#define HERMETIAN_I(matrix,i,j)  ( i==j ? 0.0 : j>i ? matrix[j*num_orbs+i] : matrix[i*num_orbs+j] )


real tic_values[] = {3.0,2.5,2.0,1.5,1.0,0.5,0.25,0.1,0.05,-1.0};


#ifndef USING_THE_MAC
void main(argc, argv)
  int argc;
  char **argv;
#else
void main()
#endif
{
  int infile;
  FILE *psfile,*the_file;
  int num_orbs,num_mats;
  char draw_grid;
  real xp,yp;
  real xgap,ygap;
  int i,j,curr_mat;
  real *matrix;
  real realpart,imagpart;
  real mag,max_mag;
  real gray;
  int tot_cells, cells_written;
#ifdef USING_THE_MAC
  int argc;
  char argv[4][80];

        /* set up some stuff for Sioux */
        //SIOUXSettings.standalone = FALSE;
        SIOUXSettings.asktosaveonclose = FALSE;
        SIOUXSettings.autocloseonquit = FALSE;
        printf("Starting matrix_view.\n");

  the_file = choose_mac_file(argv[1],MAC_FOPEN_OPEN_CD);
  if( !the_file ) {
          fatal("User cancelled intial file open");
  }

  printf("Enter name of output (PS) file: ");
  scanf("%s",argv[2]);
  argc = 3;


  /* get the command line arguments */
//  argc = ccommand(&argv);

#endif
  if(argc < 3){
    fatal("Usage: matrix_view infile outfile");
  }

  if( argc == 4 && argv[3][0] == '1' ) draw_grid = 1;
  else draw_grid = 0;

  /* open the file */
#ifndef USING_THE_MAC
  infile = open(argv[1],O_RDONLY,"r");
#else
  infile = open(argv[1],O_RDONLY);
#endif
  if( infile == -1 ){
    error("Can't open file for binary I/O.");
    return;
  }

  /* read out the number of matrices */
  read(infile,(const char *)&num_mats,sizeof(int));
  /* now the number of orbitals */
  read(infile,(const char *)&num_orbs,sizeof(int));


  /* get space to store the matrices */
  matrix = (real *)calloc(num_orbs*num_orbs,sizeof(real));
  if(!matrix) fatal("Can't get space for the matrix\n");

  /* open the ps file */
  psfile = fopen(argv[2],"w+");
  if(!psfile) fatal("Can't open ps file\n");

  /* write out some header information */
  fprintf(psfile,"%%!PS-Adobe-2.0 EPSF-1.2\n");
  fprintf(psfile,"%%%%DocumentFonts: Times-Roman \n");
  fprintf(psfile,"%%%%Creator: matrix_view\n");
  fprintf(psfile,"%%%%Pages: %d\n",num_mats);
  fprintf(psfile,"%%%%BoundingBox: %d %d %d %d\n",START_X-20,START_Y-20,
          BOX_WIDTH+START_X+100,BOX_HEIGHT+START_Y+20);
  fprintf(psfile,"%%%%EndComments\n\n");

  fprintf(psfile,"/M {moveto} def\n");
  fprintf(psfile,"/L {lineto} def\n");
  fprintf(psfile,"/R {rmoveto} def\n");
  fprintf(psfile,"/RL {rlineto} def\n");
  fprintf(psfile,"/SG {setgray} def\n");

  xgap = (real)BOX_WIDTH / (real)num_orbs;
  ygap = (real)BOX_HEIGHT / (real)num_orbs;
  /* now read in the matrices one at a time */
  for(curr_mat=0;curr_mat<num_mats;curr_mat++){
    read(infile,(const char *)matrix,num_orbs*num_orbs*sizeof(real));

    /* find the maximum value */
    max_mag = 0.0;
    for(i=0;i<num_orbs;i++){
      for(j=0;j<num_orbs;j++){
        realpart=HERMETIAN_R(matrix,i,j);
        imagpart=HERMETIAN_I(matrix,i,j);
        mag = sqrt(realpart*realpart + imagpart*imagpart);
        if( mag > max_mag ) max_mag = mag;
      }
    }

    fprintf(psfile,"%%%%Page: %d %d\n",curr_mat+1, curr_mat+1);
    fprintf(psfile,"/Times-Roman findfont 12 scalefont setfont\n");
    fprintf(psfile,"%d %d M (Max: %6.4lf) show\n",
            START_X+BOX_WIDTH+10,START_Y+BOX_HEIGHT/2,max_mag);

    xp = START_X+BOX_WIDTH+20;
    yp = START_Y+BOX_HEIGHT/2-20;
    for( i=0; tic_values[i] != -1.0; i++){
      if( tic_values[i] <= max_mag ){
        fprintf(psfile,"%d %d M %d %d RL %d %d RL %d %d RL %d %d RL %4.2lf SG fill 0 SG stroke\n",
                (int)xp,(int)yp,10,0,0,-10,-10,0,0,10,
                fabs(max_mag-tic_values[i])/max_mag);
        fprintf(psfile,"%d %d M (%4.2lf) show\n",(int)xp+15,(int)yp-5,tic_values[i]);
        yp -= 20;
      }
    }

    tot_cells = 0;
    cells_written = 0;

    /* now fill in the blocks */
    for(i=0;i<num_orbs;i++){
      for(j=0;j<num_orbs;j++){
        realpart=HERMETIAN_R(matrix,i,j);
        imagpart=HERMETIAN_I(matrix,i,j);
        mag = sqrt(realpart*realpart + imagpart*imagpart);

        gray = fabs(max_mag-mag)/max_mag;
        xp = (real)START_X + (real)i*xgap;
        yp = (real)START_Y + (real)(num_orbs-j)*ygap;

        tot_cells++;
        if( mag >= .05 ){
          cells_written++;
          fprintf(psfile,
                  "%4.2lf %4.2lf M %4.2lf %4.2lf L %4.2lf %4.2lf L %4.2lf %4.2lf L %4.2lf %4.2lf L %4.2lf SG fill\n",
                  xp,yp,xp+xgap,yp,xp+xgap,yp-ygap,xp,yp-ygap,xp,yp,gray);
        }
      }
    }
    if( draw_grid ){
      /* draw the grid */
      fprintf(psfile,"0 SG\n");
      for(i=0;i<=num_orbs;i++){
        fprintf(psfile,"%4.2lf %4.2lf M\n",((real)START_X+xgap*i),(real)START_Y);
        fprintf(psfile,"%4.2lf %4.2lf L\n",((real)START_X+xgap*i),
                (real)(START_Y+BOX_HEIGHT));
        fprintf(psfile,"%4.2lf %4.2lf M\n",(real)START_X,(real)(START_Y+ygap*i));
        fprintf(psfile,"%4.2lf %4.2lf L\n",(real)(START_X+BOX_WIDTH),
                (real)(START_Y+ygap*i));
        fprintf(psfile,"stroke\n");
      }
    }
    fprintf(psfile,"showpage\n");
    mag = .05;
    printf("%d cells (of %d) with magnitude > %lf were written.\n",cells_written,tot_cells,
           mag);

  }

  fclose(psfile);
  close(infile);

}



