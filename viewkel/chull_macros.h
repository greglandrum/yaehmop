/*====================================================================
    macros.h
 
 	macros used to access data structures and perform quick tests.

  ====================================================================*/

/* general-purpose macros */
#define LSWAP(t,x,y)	{ t = x; x = y; y = t; }

#ifdef NEED_BOGUS_MALLOC_PROTO
char *D_MALLOC();
#endif

#define NEW(p,type)	if ((p=(type *) D_MALLOC (sizeof(type))) == NULL) {\
				printf ("Out of Memory!\n");\
				exit(0);\
			}

#define FREE(p)		if (p) { free ((char *) p); p = NULL; }


#define ADD( head, p )  if ( head )  { \
				p->next = head->next; \
				p->prev = head; \
				head->next = p; \
				p->next->prev = p; \
			} \
			else { \
				head = p; \
				head->next = head->prev = p; \
			}

#define DELETE( head, p )   if ( head )  { \
				if ( head == head->next ) \
					head = NULL;  \
				else if ( p == head ) \
					head = head->next; \
				p->next->prev = p->prev;  \
				p->prev->next = p->next;  \
				FREE( p ); \
			} 

