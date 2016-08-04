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

/*******

  mem_debug.c:
    This contains a set of replacements for the standard
    memory allocation functions malloc, calloc, realloc, and free
    in order to allow the tracking down of core leaks and various
    other bad memory allocation problems.

    Here's how I use it:

1)  somewhere in a header file included by everything else,
    add the following definitions:
#ifndef MEM_DEBUG
#define D_MALLOC(_a_) malloc(_a_)
#define D_CALLOC(_a_,_b_) calloc(_a_,_b_)
#define D_REALLOC(_a_,_b_) realloc(_a_,_b_)
#define D_FREE(_a_) free(_a_)
#else
#define D_MALLOC(_a_) d_malloc(_a_,__FILE__,__LINE__)
#define D_CALLOC(_a_,_b_) d_calloc(_a_,_b_,__FILE__,__LINE__)
#define D_REALLOC(_a_,_b_) d_realloc(_a_,_b_,__FILE__,__LINE__)
#define D_FREE(_a_) d_free(_a_,__FILE__,__LINE__)
#endif

2) replace all memory allocation calls with the corresponding D_ call.
3) recompile everything with -DMEM_DEBUG
4) add a call to d_check_core() at the end of the program
5) if a given leak/problem is really difficult to find, try
 recompiling mem_debug.c with -DGAG_ME_WITH_VERBOSITY
 this will print out a short record of every transaction.

*********/

#include <stdio.h>
#include <signal.h>

#ifdef MEM_DEBUG
/****
  this is the structure used to maintain the list of
  allocated memory chunks.
****/
typedef struct malloc_record_type_def malloc_record_type;
struct malloc_record_type_def{
  char *address;  /* the pointer itself */
  int size;       /* the amount of memory associated with this pointer */
  char loc_file[240]; /* which file was the allocation done in */
  int loc_line;   /* which line of that file */
  malloc_record_type *prev,*next; /* do a doubly-linked list */
};


malloc_record_type *records=0;
int num_allocated = 0;

/*****

  here's the malloc replacement

*****/
char *d_malloc(int size,char *filename,int linenum){
  char *temp_ptr;
  malloc_record_type *new_record;

  /* start with the memory */
  temp_ptr = (char *)malloc(size);
  num_allocated++;
  if( temp_ptr ){
    /***
      allocate space for the malloc record, fill in the
      requisite information and slap it at the head of
      the list.
    ***/
    new_record = (malloc_record_type *)malloc(sizeof(malloc_record_type));
    if(!new_record) fatal("can't allocate a malloc_record.");
    new_record->prev = 0;
    new_record->next = 0;
    new_record->address = temp_ptr;
    new_record->size = size;
    new_record->loc_line = linenum;
    strncpy(new_record->loc_file,filename,240);
    new_record->next = records;
    if( records ) records->prev = new_record;
    records = new_record;
#ifdef GAG_ME_WITH_VERBOSITY
    printf("\tmalloc from line\t %d of file\t %s of\t %d bytes ok.\t(%0xd)\n",
           linenum,filename,size,temp_ptr);
#endif
  }

  return(temp_ptr);
}

/*****

  here's the calloc replacement, it's the same as d_malloc

*****/
char *d_calloc(int num,int size,char *filename,int linenum){
  char *temp_ptr;
  malloc_record_type *new_record;

  temp_ptr = (char *)calloc(num,size);
  num_allocated++;
  if( temp_ptr ){
    new_record = (malloc_record_type *)malloc(sizeof(malloc_record_type));
    if(!new_record) fatal("can't allocate a malloc_record.");
    new_record->prev = 0;
    new_record->next = 0;

    new_record->address = temp_ptr;
    new_record->size = size*num;
    new_record->loc_line = linenum;
    strncpy(new_record->loc_file,filename,240);
    new_record->next = records;
    if( records ) records->prev = new_record;
    records = new_record;
#ifdef GAG_ME_WITH_VERBOSITY
    printf("\tcalloc from line\t %d of file\t %s of\t %d bytes ok.\t(%0xd)\n",
           linenum,filename,size*num,temp_ptr);
#endif

  }

  return(temp_ptr);
}

/*******

  print out the information for each allocated record
   this is mainly useful for debugging purposes.

*******/
void d_crawl()
{
  malloc_record_type *record;

  record = records;
  while(record){
    printf("\t%0xd: line %d, file %s, size %d\n",record->address,
           record->loc_line,record->loc_file,record->size);
    record = record->next;
  }
}


/*****

  here's the free replacement

*****/
void d_free(char *address,char *filename,int linenum){
  malloc_record_type *record;

  record = records;

  /* find the relevant record in the list */
  while( record && record->address != address ) record = record->next;
  if( !record ){
    fprintf(stderr,"free called on a pointer which was not allocated, (%0xd).\n",address);
    fprintf(stderr,"\tfrom line %d of file %s\n",linenum,filename);
    kill(getpid(),SIGSEGV);
  }

  /* pull it back out of the list */
  if(record==records){
    records = record->next;
    if(records) records->prev = 0;
  }
  if(record->prev) record->prev->next = record->next;
  if(record->next) record->next->prev = record->prev;

  /* free the memory */
  free(record->address);
#ifdef GAG_ME_WITH_VERBOSITY
    printf("\tfree from line\t %d of file\t %s of\t %d bytes ok.\t(%0xd)\n",
           record->loc_line,record->loc_file,record->size,address);
#endif
  free(record);
  num_allocated--;
}

/******

  replacement for realloc

*******/
char *d_realloc(char *address,int size,char *filename,int linenum){
  char *temp_ptr;
  malloc_record_type *record;

  record = records;
  while( record && record->address != address ) record = record->next;
  if( !record ){
    fprintf(stderr,"realloc called on a pointer which was not allocated.\n");
    fprintf(stderr,"\tfrom line %d of file %s\n",linenum,filename);
    kill(getpid(),SIGSEGV);
  }

  temp_ptr = (char *)realloc(address,size);
  if( temp_ptr ){
    record->address = temp_ptr;
    record->size = size;
    record->loc_line = linenum;
    strncpy(record->loc_file,filename,240);
#ifdef GAG_ME_WITH_VERBOSITY
    printf("\trealloc from line\t %d of file\t %s of\t %d bytes ok.\t(%0xd)\n",
           linenum,filename,size,temp_ptr);
#endif

  }
  return(temp_ptr);
}

/******

  check for and display any leaks

*******/
void d_check_core(){
  int num_leaks;
  malloc_record_type *record;

  record = records;
  num_leaks = 0;
  if( record ){
    printf("\td_check_core: %d leaks present. Argh!\n",num_allocated);
    while(record){
      printf("\t\t%d) Line %d of file %s.  Size: %d\n",num_leaks+1,
             record->loc_line,record->loc_file,record->size);
      num_leaks++;
      record = record->next;
    }
  }else{
    printf("\td_check_core: We're clean!\n");
  }
}
#endif
