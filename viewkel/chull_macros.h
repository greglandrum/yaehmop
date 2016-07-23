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
/*====================================================================
    macros.h

        macros used to access data structures and perform quick tests.

  ====================================================================*/

/* general-purpose macros */
#define LSWAP(t, x, y)                                                         \
  {                                                                            \
    t = x;                                                                     \
    x = y;                                                                     \
    y = t;                                                                     \
  }

#ifdef NEED_BOGUS_MALLOC_PROTO
char *D_MALLOC();
#endif

#define NEW(p, type)                                                           \
  if ((p = (type *)D_MALLOC(sizeof(type))) == NULL) {                          \
    printf("Out of Memory!\n");                                                \
    exit(0);                                                                   \
  }

#define FREE(p)                                                                \
  if (p) {                                                                     \
    free((char *)p);                                                           \
    p = NULL;                                                                  \
  }

#define ADD(head, p)                                                           \
  if (head) {                                                                  \
    p->next = head->next;                                                      \
    p->prev = head;                                                            \
    head->next = p;                                                            \
    p->next->prev = p;                                                         \
  } else {                                                                     \
    head = p;                                                                  \
    head->next = head->prev = p;                                               \
  }

#define DELETE(head, p)                                                        \
  if (head) {                                                                  \
    if (head == head->next)                                                    \
      head = NULL;                                                             \
    else if (p == head)                                                        \
      head = head->next;                                                       \
    p->next->prev = p->prev;                                                   \
    p->prev->next = p->next;                                                   \
    FREE(p);                                                                   \
  }
