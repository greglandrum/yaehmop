/*

  This is a program to generate a YAeHMOP input file of
  arbitrary size from a simpler format input file.

  created by greg Landrum December, 1995


  file format: (lines starting with semicolons and blank lines are ignored)
    dimensionality of the system (1,2, or 3)
    number of overlaps along each direction
    number of electrons per cell
    num_atoms
    the atoms:  num symbol x y z
    the geometry line (i.e. Geometry Crystallographic or Geometry)
    any additional lines are passed on verbatim to the bind input file.

    NOTE:  it is assumed that the lattice vectors are defined by the first
        and last ndim atoms.
*/
#include <stdio.h>

#define real double

typedef struct{
  real x,y,z;
} point_type;

typedef struct{
  char type[4];
  point_type loc;
} atom_type;

typedef struct{
  int num_dim;
  int num_raw_atoms,num_atoms;
  real num_raw_electrons,num_electrons;
  atom_type *raw_atoms;
  atom_type *atoms;
  point_type vects[3];
  int overlaps[3];
} molec_type;


#define ERROR 11
#define FATAL 22
#define IGNORE 33
#define MAX_STR_LEN 400


void fatal( errorstring )
  char *errorstring;
{
  fprintf( stderr, "FATAL ERROR: %s.\nExecution Terminated.\n",
      errorstring );
  exit(-1);
}

/* Procedure error
 * prints an error message
 *
 * in case the job was queued, the error message is echoed to the status file
 */
void error( errorstring )
  char *errorstring;
{
  fprintf( stderr, "ERROR: %s.\n", errorstring );
  return;
}


/****************************************************************************
 *
 *                   Procedure skipcomments
 *
 * Arguments: file : a pointer to file type
 *        instring : pointer to type char
 *          toggle : a char
 * Returns: an integer
 *
 * Action: Reads in lines from 'file' until one is hit that does not begin
 *     with a ; or a return. puts the first non-comment line into instring
 *     and then returns.
 *
 *    if 'toggle is set to FATAL then hitting EOF will result in
 *      program termination with a call to fatal.
 *    if 'toggle is set to ERROR then EOF results in a call to error then
 *        the function returns.
 *    if 'toggle is set to IGNORE then EOF is ignored.
 *
 *   in any case, if the function returns and EOF has been hit the return
 *    value is -1.
 *
 ****************************************************************************/
int skipcomments(FILE *file,char *string,char toggle)
{
  int i;

  /* use the first element of string to check for EOF */
  string[0] = 0;
  fgets(string,MAX_STR_LEN,file);

  /*******
    deal with the fact that the string may contain only spaces, which
    we will want to skip
  ********/
  i = 0;
  while(string[i] == ' ') i++;
  while( string[i] == '\n' || string[i] == ';'
        && string[i] != 0 ){
    string[0] = 0;
    fgets(string,MAX_STR_LEN,file);
    i = 0;
    while(string[i] == ' ') i++;
  }

  if( string[0] != 0 ) return(0);
  else{
    switch(toggle){
    case FATAL:
      fatal("End of File (EOF) hit in skipcomments.");
      break;
    case ERROR:
      error("End of File (EOF) hit in skipcomments, execution continuing.");
      break;
    case IGNORE:
      break;
    }
    return(-1);
  }
}



void read_from_file(FILE *infile, molec_type *molec)
{
  char instring[MAX_STR_LEN];
  int foo_int;
  int i;

  skipcomments(infile,instring,FATAL);
  sscanf(instring,"%d",&molec->num_dim);
  skipcomments(infile,instring,FATAL);
  sscanf(instring,"%d %d %d",&molec->overlaps[0],
         &molec->overlaps[1],&molec->overlaps[2]);

  skipcomments(infile,instring,FATAL);
  sscanf(instring,"%lf",&molec->num_raw_electrons);
  skipcomments(infile,instring,FATAL);
  sscanf(instring,"%d",&molec->num_raw_atoms);

  /* get memory to store the raw atoms */
  molec->raw_atoms = (atom_type *)malloc(molec->num_raw_atoms*molec->num_dim*sizeof(atom_type));
  if( !molec->raw_atoms ){
    fatal("Can't get space for raw_atoms");
  }

  for(i=0;i<molec->num_raw_atoms;i++){
    skipcomments(infile,instring,FATAL);
    sscanf(instring,"%d %s %lf %lf %lf",&foo_int,molec->raw_atoms[i].type,
           &(molec->raw_atoms[i].loc.x),&(molec->raw_atoms[i].loc.y),
           &(molec->raw_atoms[i].loc.z));
  }

  /* figure out the lattice vectors */
  molec->num_raw_atoms -= molec->num_dim;
  for( i=0; i<molec->num_dim; i++){
    molec->vects[i].x = molec->raw_atoms[molec->num_raw_atoms+i].loc.x -
      molec->raw_atoms[0].loc.x;
    molec->vects[i].y = molec->raw_atoms[molec->num_raw_atoms+i].loc.y -
      molec->raw_atoms[0].loc.y;
    molec->vects[i].z = molec->raw_atoms[molec->num_raw_atoms+i].loc.z -
      molec->raw_atoms[0].loc.z;
  }

  /* that's that */
}

void grow_it(molec_type *solid,int num_a,int num_b, int num_c)
{
  int i,j,k,l;
  atom_type *temp_atoms;
  point_type dist_a,dist_b,dist_c;
  int num_added;

  /* get the memory we'll need */
  solid->num_atoms = solid->num_raw_atoms*num_a*num_b*num_c;
  temp_atoms = (atom_type *)calloc(solid->num_atoms,sizeof(atom_type));
  if( !temp_atoms ) fatal("Memory allocation: can't get space for more atoms.");


    /* now grow the crystal */
  num_added = 0;
  for(i=0;i<num_a;i++){
    dist_a.x = i*solid->vects[0].x;
    dist_a.y = i*solid->vects[0].y;
    dist_a.z = i*solid->vects[0].z;

    for(j=0;j<num_b;j++){
      dist_b.x = j*solid->vects[1].x;
      dist_b.y = j*solid->vects[1].y;
      dist_b.z = j*solid->vects[1].z;

      for(k=0;k<num_c;k++){
        dist_c.x = k*solid->vects[2].x;
        dist_c.y = k*solid->vects[2].y;
        dist_c.z = k*solid->vects[2].z;

        /* first copy in the old atom data */
        bcopy((char *)solid->raw_atoms,(char *)&(temp_atoms[num_added]),
              solid->num_raw_atoms*sizeof(atom_type));

        /* now update the locations */
        for(l=0;l<solid->num_raw_atoms;l++){
          temp_atoms[num_added].loc.x =
            solid->raw_atoms[l].loc.x + dist_a.x + dist_b.x + dist_c.x;
          temp_atoms[num_added].loc.y =
            solid->raw_atoms[l].loc.y + dist_a.y + dist_b.y + dist_c.y;
          temp_atoms[num_added].loc.z =
            solid->raw_atoms[l].loc.z + dist_a.z + dist_b.z + dist_c.z;
          num_added++;
        }
      }
    }
  }

  solid->atoms = temp_atoms;

  solid->num_electrons = solid->num_raw_electrons * num_a * num_b * num_c;

  solid->vects[0].x *= num_a;solid->vects[0].y *= num_a;solid->vects[0].z *= num_a;
  solid->vects[1].x *= num_b;solid->vects[1].y *= num_b;solid->vects[1].z *= num_b;
  solid->vects[2].x *= num_c;solid->vects[2].y *= num_c;solid->vects[2].z *= num_c;
}


void main(int argc,char **argv)
{
  FILE *infile,*outfile;
  int i;
  char instring[MAX_STR_LEN];
  molec_type molec;
  int num_a,num_b,num_c;

  if( argc != 3 )
    fatal("usage: grow_xtal <infilename> <outfilename>\n");

  /* open the infile */
  infile = fopen(argv[1],"r");
  if( !infile ) fatal("Can't open input file.");



  /* read the data */
  read_from_file(infile,&molec);

  /* prompt for the size */
  printf("The file system has: %d atoms and is %d dimensional\n",
         molec.num_raw_atoms,molec.num_dim);

  printf("Please enter the number of cells along each lattice direction on separate lines.\n");
  printf("(a)  ");
  scanf("%d",&num_a);
  if( num_a < 1 ){
    error("Don't enter dumb values!");
    num_a = 1;
  }
  if( molec.num_dim > 1 ){
    printf("(b)  ");
    scanf("%d",&num_b);
    if( num_b < 1 ){
      error("Don't enter dumb values!");
      num_b = 1;
    }
  }
  else{
    num_b = num_c = 1;
  }
  if( molec.num_dim > 2 ){
    printf("(c)  ");
    scanf("%d",&num_c);
    if( num_c < 1 ){
      error("Don't enter dumb values!");
      num_c = 1;
    }
  }
  else{
    num_c = 1;
  }

  /* okey dokey, grow it! */
  grow_it(&molec,num_a,num_b,num_c);

  /* open the output file  */
  outfile = fopen(argv[2],"w+");
  if( !outfile ) fatal("Can't open outfile");

  /* do some headers */
  fprintf(outfile,"Not new3!\n");
  fprintf(outfile,"Automatically Generated input file\n");

  /* read the geometry line and spew it again */
  skipcomments(infile,instring,FATAL);
  fputs(instring,outfile);

  printf("Writing %d atoms\n",molec.num_atoms);
  /* now write the number of atoms and the atoms themselves */
  fprintf(outfile,"%d\n",molec.num_atoms+molec.num_dim);
  for(i=0;i<molec.num_atoms;i++){
    fprintf(outfile,"%d %s %lf %lf %lf\n",i+1,molec.atoms[i].type,
            molec.atoms[i].loc.x,molec.atoms[i].loc.y,
            molec.atoms[i].loc.z);
  }
  /* write the lattice vectors */
  for(i=0;i<molec.num_dim;i++){
    fprintf(outfile,"%d %s %lf %lf %lf\n",molec.num_atoms+i+1,
            molec.atoms[i].type,
            molec.raw_atoms[0].loc.x + molec.vects[i].x,
            molec.raw_atoms[0].loc.y + molec.vects[i].y,
            molec.raw_atoms[0].loc.z + molec.vects[i].z);
  }

  /* the number of electrons */
  fprintf(outfile,"Electrons\n %lf\n",molec.num_electrons);

  /* the lattice stuff */
  fprintf(outfile,"Lattice\n%d\n",molec.num_dim);
  for(i=0;i<molec.num_dim;i++){
    fprintf(outfile,"%d ",molec.overlaps[i]);
  }
  fprintf(outfile,"\n");
  for(i=0;i<molec.num_dim;i++){
    fprintf(outfile,"1 %d\n",molec.num_atoms+i+1);
  }

  /* now read and write back out the rest of the file */
  while(skipcomments(infile,instring,IGNORE)>=0){
    fputs(instring,outfile);
  }

  printf("Done. Don't forget to put *'s back in for atoms that need parms\n");

  fclose(outfile);
}



