#include <stdio.h>

main()
{
  char instring[80];

  strcpy(instring,"# NUMBER OF ATOMS: ");
 
  printf("instring[0] = %c\n",instring[0]);
  printf("strstr: %d\n",strstr(instring,"DENSITY"));
  while( instring[0] != '#' && !strstr(instring,"DENSITY") )
    printf("Doing okay!\n");

  printf("Whoops!\n");
}
