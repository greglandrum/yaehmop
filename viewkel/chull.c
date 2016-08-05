/*
        chull.c
        Written by Joseph O'Rourke
                and John Kutcher, Catherine Schevon, Susan Weller.
        Last modified: 14 August 1995

        Modified by greg Landrum for inclusion in viewkel
         03.05.98

*/
#include "viewkel.h"
#include "chull_macros.h"

/***
  Recent Edit History
   26.09.98 gL:
     various modifications to remove warning when compiled with -Wall
     under gcc (yeah yeah... it's anal)
***/

/* Define Boolean type */
#ifdef USE_ENUM_BOOL
typedef        enum { FALSE, TRUE }        bool;
#else
typedef char bool;
#ifndef FALSE
#define FALSE 0
#define TRUE 1
#endif
#endif


/* Define vertex indices. */
#define LXD   0
#define LYD   1
#define LZD   2

/* Define structures for vertices, edges and faces */
typedef struct tVertexStructure tsVertex;
typedef tsVertex *tVertex;

typedef struct tEdgeStructure tsEdge;
typedef tsEdge *tEdge;

typedef struct tFaceStructure tsFace;
typedef tsFace *tFace;

struct tVertexStructure {
        float        v[3];
        int        vnum;
        tEdge        duplicate;        /* pointer to incident cone edge (or NULL) */
        bool         onhull;                /* T iff point on hull. */
        bool        mark;                /* T iff point already processed. */
        tVertex next, prev;
};

struct tEdgeStructure {
        tFace         adjface[2];
        tVertex endpts[2];
        tFace        newface;        /* pointer to incident cone face. */
        bool         delete;                /* T iff edge should be delete. */
        tEdge         next, prev;
};

struct tFaceStructure {
        tEdge         edge[3];
        tVertex vertex[3];
        bool        visible;        /* T iff face visible from new point. */
        tFace         next, prev;
};

/* Define flags */
#define ONHULL           TRUE
#define REMOVED          TRUE
#define VISIBLE          TRUE
#define PROCESSED        TRUE

/* Global variable definitions */
tVertex vertices        = NULL;
tEdge edges                    = NULL;
tFace faces                    = NULL;
bool debug = FALSE;

tVertex vertices;
tEdge         edges;
tFace         faces;
bool         debug;

/* Function declarations */
tVertex MakeVertex( void );
void    CopyVertices( atom_type *,int);
void    CopyTris( triangle_list_type **);
void    Print( void );
void    Tetrahedron( void );
void    ConstructHull( void );
bool        AddOne( tVertex p );
float     Volume6(tFace f, tVertex p);
tFace        MakeStructs( tEdge e, tVertex p );
void    MakeCcw( tFace f, tEdge e, tVertex p );
tEdge   MakeEdge( void );
tFace   MakeFace( void );
void    CleanUp( void );
void    CleanEdges( void );
void    CleanFaces( void );
void    CleanVertices( void );
bool        Collinear( tVertex a, tVertex b, tVertex c );
void    CheckEuler(int, int, int );
double         Volumed( tFace, tVertex );
void        PrintPoint( tVertex p );
void    Checks( void );
void        Consistency( void );
void        Convexity( void );
void        PrintOut( tVertex v );
void        PrintVertices( void );
void        PrintEdges( void );
void        PrintFaces( void );
void    BlowOutGlobals(void);


/*--------------------------------------------------------------------------*/
void gen_chull( atom_type *active_atoms, int num_atoms, triangle_list_type **tri_list )
{

  debug = 0;
  CopyVertices(active_atoms,num_atoms);
  Tetrahedron();
  ConstructHull();
  CopyTris(tri_list);
  BlowOutGlobals();
}

void BlowOutGlobals(void)
{
  tVertex vert=vertices,tv;
  tEdge edge=edges,te;
  tFace face=faces,tf;

  while(vert){
    tv = vert->next;
    D_FREE(vert);
    vert = tv;
    if(vert==vertices) vert=0;
  }
  while(edge){
    te = edge->next;
    D_FREE(edge);
    edge = te;
    if(edge==edges) edge=0;
  }
  while(face){
    tf = face->next;
    D_FREE(face);
    face = tf;
    if(face==faces) face=0;
  }

  vertices = 0;
  edges = 0;
  faces = 0;
}


/*------------------------------------------------------------------
MakeVertex: Makes a vertex, nulls out fields.
--------------------------------------------------------------------*/
tVertex        MakeVertex( void )
{
        tVertex        v;

        NEW( v, tsVertex );
        v->duplicate = NULL;
        v->onhull = !ONHULL;
        v->mark = !PROCESSED;
        ADD( vertices, v );

        return v;
}
/*------------------------------------------------------------------
CopyVertices: Reads in the vertices, and links them into a circular
list with MakeVertex.  There is no need for the # of vertices to be
the first line: the function looks for EOF instead.
--------------------------------------------------------------------*/
void        CopyVertices( atom_type *active_atoms, int num_atoms )
{
  tVertex v;
  int        vnum = 0;
  int i;

  for(i=0;i<num_atoms;i++){
    v = MakeVertex();
    v->v[LXD] = active_atoms[i].loc.x;
    v->v[LYD] = active_atoms[i].loc.y;
    v->v[LZD] = active_atoms[i].loc.z;
    v->vnum = vnum++;
  }
}

/*------------------------------------------------------------------
Print: Prints out the vertices and the faces.  Uses the vnum indices
corresponding to the order in which the vertices were input.
--------------------------------------------------------------------*/
void        CopyTris( triangle_list_type **tri_list )
{
  /* Pointers to vertices, edges, faces. */
  triangle_list_type *tris=0,*new_tri;
  tVertex        v;
  tEdge   e;
  tFace   f;
  /* Counters for Euler's formula. */
  int         V = 0, E = 0 , F = 0;
  /* Note: lowercase==pointer, uppercase==counter. */

  /* Vertices. */
  printf("\n");
  v = vertices;
  do {
    V++;
    v = v->next;
  } while ( v != vertices );
  printf("Vertices:\tV = %d\n", V );

  printf("index:\tx\ty\tz\n");
  do {
    printf("%5d:\t%f\t%f\t%f\n",
           v->vnum, v->v[LXD], v->v[LYD], v->v[LZD] );
    v = v->next;
  } while ( v != vertices );

  /* Faces. */
  f = faces;
  do {
    ++F;
    f  = f ->next;
  } while ( f  != faces );
  printf("Faces:\t\tF = %d\n", F );

  printf("\tv0\tv1\tv2\t(vertex indices)\n");
  do {
    new_tri = (triangle_list_type *)D_CALLOC(1,sizeof(triangle_list_type));
    if(!new_tri) fatal("can't allocate a triangle\n");
    printf("\t%d\t%d\t%d\n",
           f->vertex[0]->vnum,
           f->vertex[1]->vnum,
           f->vertex[2]->vnum );
    new_tri->tri.vertices[0] = f->vertex[0]->vnum;
    new_tri->tri.vertices[1] = f->vertex[1]->vnum;
    new_tri->tri.vertices[2] = f->vertex[2]->vnum;

    new_tri->next = tris;
    tris = new_tri;

    f = f->next;
  } while ( f != faces );

  /* Edges. */
  e = edges;
  do {
    E++;
    e = e->next;
  } while ( e != edges );
  printf("Edges:\tE = %d\n", E );
  /* Edges not printed out (but easily added). */
  debug = TRUE;
  CheckEuler( V, E, F );

  *tri_list = tris;

}

/*-----------------------------------------------------------------------
 Tetrahedron builds the initial tetrahedron.  It first finds 3 noncollinear
 points and makes a face out of them, and then finds a fourth point that
 is not coplanar with that face.  The vertices are stored in the face
 structure in counterclockwise order so that the volume between the face
 and the point is negative. Lastly, the 3 newfaces to the fourth point
 are constructed and the data structures are cleaned up.
 -----------------------------------------------------------------------*/
void        Tetrahedron( void )
{
        tVertex v1, v4, t;
        tFace         f;
        tEdge         e1, e2, e3, s;
        float         vol;


        /* Find 3 non-Collinear points. */
        v1 = vertices;
        while ( Collinear( v1, v1->next, v1->next->next ) )
                            if ( ( v1 = v1->next ) == vertices ) {
                      printf("All points are Collinear!\n");
                      exit(0);
                }

        /* Mark the vertices as processed. */
        v1->mark = PROCESSED;
        v1->next->mark = PROCESSED;
        v1->next->next->mark = PROCESSED;

        /* Create edges of the initial triangle. */
        e1 = MakeEdge();
        e2 = MakeEdge();
        e3 = MakeEdge();
        e1->endpts[0] = v1;              e1->endpts[1] = v1->next;
        e2->endpts[0] = v1->next;        e2->endpts[1] = v1->next->next;
        e3->endpts[0] = v1->next->next;  e3->endpts[1] = v1;

        /* Create face for triangle. */
        f = MakeFace();
        f->edge[0] = e1;   f->edge[1] = e2;   f->edge[2] = e3;
        f->vertex[0] = v1;  f->vertex[1] = v1->next;
        f->vertex[2] = v1->next->next;

        /* Link edges to face. */
        e1->adjface[0] = e2->adjface[0] = e3->adjface[0] = f;


        /* Find a fourth, non-coplanar point to form tetrahedron. */
        v4 = v1->next->next->next;
        vol = Volume6( f, v4 );
        while ( !vol )   {
                      if ( ( v4 = v4->next ) == v1 ) {
                        printf("All points are coplanar!\n");
                       return;
                }
                       vol = Volume6( f, v4 );
               }
        v4->mark = PROCESSED;

        /* Store vertices in ccw order. */
        if( vol < 0 ) {
            LSWAP( t, f->vertex[1], f->vertex[2] );
            LSWAP( s, f->edge[1], f->edge[2] );
          }


        /* Construct the faces and edges between the original
              triangle and the fourth point. */
        e1->adjface[1] = MakeStructs( e1, v4 );
        e2->adjface[1] = MakeStructs( e2, v4 );
        e3->adjface[1] = MakeStructs( e3, v4 );

        CleanUp();
}

/*-------------------------------------------------------------------------
ConstructHull adds the vertices to the hull one at a time.  The hull
vertices are those in the list marked as onhull.
  -------------------------------------------------------------------------*/
void        ConstructHull( void )
{
        tVertex        v, vnext;
        tFace         f;
        bool        changed;        /* T if addition changes hull; not used. */

        v = vertices;
        f = faces;
        do {
                vnext = v->next;
                if ( !v->mark ) {
                        v->mark = PROCESSED;
                        changed = AddOne( v );
                        CleanUp();
                        if ( debug ) {
                                fprintf(stderr,"after cleanup:\n");
                                PrintOut( v );
                                Checks();
                        }
                }
                v = vnext;
        } while ( v != vertices );
}
/*-------------------------------------------------------------------------
AddOne is passed a vertex.  It first determines all faces visible from
that point.  If none are visible then the point is marked as not onhull.
Next is a loop over edges.  If both faces adjacent to an edge are
visible, then the edge is marked for deletion.  If just one of the adjacent
faces is visible then a new face is constructed.
--------------------------------------------------------------------------*/
bool         AddOne( tVertex p )
{
        tFace         f;
        tEdge         e;
        float         vol;
        bool        vis = FALSE;

        /* Mark faces visible from p. */
        f = faces;
        do {
                vol = Volume6( f, p );
/*
                if (debug) fprintf(stderr,"faddr: %6x   paddr: %6x   Vol = %f\n",f,p,vol);
*/
                if ( vol < 0 ) {
                        f->visible = VISIBLE;
                        vis = TRUE;
                }
                f = f->next;
        } while ( f != faces );

        /* If no faces are visible from p, then p is inside the hull. */
        if ( !vis ) {
                p->onhull = !ONHULL;
                return FALSE;
        }

        /* Mark edges in interior of visible region for deletion.
           Erect a newface based on each border edge. */
        e = edges;
        do {
                tEdge temp;
                temp = e->next;
                if ( e->adjface[0]->visible && e->adjface[1]->visible )
                        /* e interior: mark for deletion. */
                        e->delete = REMOVED;
                else if ( e->adjface[0]->visible || e->adjface[1]->visible )
                        /* e border: make a new face. */
                        e->newface = MakeStructs( e, p );
                e = temp;
        } while ( e != edges );
        return TRUE;
}

/*-------------------------------------------------------------------------
Volume6 returns six times the volume of the tetrahedron determined by f
and p.  Volume6 is positive iff p is on the negative side of f,
where the positive side is determined by the rh-rule.  So the volume
is positive if the ccw normal to f points outside the tetrahedron.
--------------------------------------------------------------------------*/
float         Volume6( tFace f, tVertex p )
{
        float         vol;
        float         ax, ay, az, bx, by, bz, cx, cy, cz, dx, dy, dz;
        float        bxdx, bydy, bzdz, cxdx, cydy, czdz;
        double        vold;
        int        i;

        ax = f->vertex[0]->v[LXD];
        ay = f->vertex[0]->v[LYD];
        az = f->vertex[0]->v[LZD];
        bx = f->vertex[1]->v[LXD];
        by = f->vertex[1]->v[LYD];
        bz = f->vertex[1]->v[LZD];
        cx = f->vertex[2]->v[LXD];
        cy = f->vertex[2]->v[LYD];
        cz = f->vertex[2]->v[LZD];
        dx = p->v[LXD];
        dy = p->v[LYD];
        dz = p->v[LZD];

        /* This is the expression used in the text.  Now replaced.
        vol =          -az * by * cx + ay * bz * cx + az * bx * cy - ax * bz * cy
                - ay * bx * cz + ax * by * cz + az * by * dx - ay * bz * dx
                - az * cy * dx + bz * cy * dx + ay * cz * dx - by * cz * dx
                - az * bx * dy + ax * bz * dy + az * cx * dy - bz * cx * dy
                - ax * cz * dy + bx * cz * dy + ay * bx * dz - ax * by * dz
                - ay * cx * dz + by * cx * dz + ax * cy * dz - bx * cy * dz;
        */
        /* This expression is algebraically equivalent to the above, but
        uses fewer multiplications.
        vol =         -(az-dz) * (by-dy) * (cx-dx)
                + (ay-dy) * (bz-dz) * (cx-dx)
                + (az-dz) * (bx-dx) * (cy-dy)
                - (ax-dx) * (bz-dz) * (cy-dy)
                - (ay-dy) * (bx-dx) * (cz-dz)
                + (ax-dx) * (by-dy) * (cz-dz); */
        /* And this one has even fewer arithmetic operations:
           (thanks to Robert Fraczkiewicz): */
        bxdx=bx-dx;
        bydy=by-dy;
        bzdz=bz-dz;
        cxdx=cx-dx;
        cydy=cy-dy;
        czdz=cz-dz;
        vol =    (az-dz)*(bxdx*cydy-bydy*cxdx)
               + (ay-dy)*(bzdz*cxdx-bxdx*czdz)
               + (ax-dx)*(bydy*czdz-bzdz*cydy);

        /* Compare integer volume with double volume for saftey. */
        vold = Volumed( f, p );
        if (debug)
        {
                fprintf(stderr,"Face = %6x\tVertex = %d, vol(int) = %f\tvol(double) = %lf\n",
                        f,p->vnum,vol,vold);
        }
        if ( fabs( vol - vold ) >= 1.0 ) {
                printf("Likely integer overflow in volume calculation\n");
                printf("because coordinates are too large:\n");
                printf("\tvol(int) = %f;\tvol(double) = %lf\n");
                printf("\tfour points:\n");
                for ( i=0; i < 3; i++ )
                        PrintPoint( f->vertex[i] );
                PrintPoint( p );
                printf("Basing decision based on double volume.\n");
                /* Return based on double volume. */
                if ( vold >= 1.0 )
                        return        1;
                else if ( vold <= -1.0 )
                        return        -1;
                else        return        0;
        }
        return vol;
}
/*-------------------------------------------------------------------------
Volumed is the same as Volume6 but computed with doubles.  For protection
against overflow.
--------------------------------------------------------------------------*/
double         Volumed( tFace f, tVertex p )
{
        double        vol;
        double         ax, ay, az, bx, by, bz, cx, cy, cz, dx, dy, dz;
        double        bxdx, bydy, bzdz, cxdx, cydy, czdz;

        ax = f->vertex[0]->v[LXD];
        ay = f->vertex[0]->v[LYD];
        az = f->vertex[0]->v[LZD];
        bx = f->vertex[1]->v[LXD];
        by = f->vertex[1]->v[LYD];
        bz = f->vertex[1]->v[LZD];
        cx = f->vertex[2]->v[LXD];
        cy = f->vertex[2]->v[LYD];
        cz = f->vertex[2]->v[LZD];
        dx = p->v[LXD];
        dy = p->v[LYD];
        dz = p->v[LZD];

        bxdx=bx-dx;
        bydy=by-dy;
        bzdz=bz-dz;
        cxdx=cx-dx;
        cydy=cy-dy;
        czdz=cz-dz;
        vol =    (az-dz)*(bxdx*cydy-bydy*cxdx)
               + (ay-dy)*(bzdz*cxdx-bxdx*czdz)
               + (ax-dx)*(bydy*czdz-bzdz*cydy);


        return vol;
}

/*-------------------------------------------------------------------------*/
void        PrintPoint( tVertex p )
{
        int        i;

        for ( i = 0; i < 3; i++ )
                printf("\t%f", p->v[i]);
        putchar('\n');
}
/*----------------------------------------------------------------------
MakeStructs makes a new face and two new edges between the
edge and the point that are passed to it. It returns a pointer to
the new face.
----------------------------------------------------------------------*/
tFace        MakeStructs( tEdge e, tVertex p )
{
        tEdge         new_edge[2];
        tFace         new_face;
        int         i, j;

        /* Make two new edges (if don't already exist). */
        for ( i=0; i < 2; ++i )
                /* If the edge exists, copy it into new_edge. */
                if ( !( new_edge[i] = e->endpts[i]->duplicate) ) {
                        /* Otherwise (duplicate is NULL), MakeEdge. */
                        new_edge[i] = MakeEdge();
                        new_edge[i]->endpts[0] = e->endpts[i];
                        new_edge[i]->endpts[1] = p;
                        e->endpts[i]->duplicate = new_edge[i];
                }

        /* Make the new face. */
        new_face = MakeFace();
        new_face->edge[0] = e;
        new_face->edge[1] = new_edge[0];
        new_face->edge[2] = new_edge[1];
        MakeCcw( new_face, e, p );

        /* Set the adjacent face pointers. */
        for ( i=0; i < 2; ++i )
        for ( j=0; j < 2; ++j )
                /* Only one NULL link should be set to new_face. */
                if ( !new_edge[i]->adjface[j] ) {
                        new_edge[i]->adjface[j] = new_face;
                        break;
                }

        return new_face;
}

/*------------------------------------------------------------------------
MakeCcw puts the vertices in the face structure in counterclockwise order.
If there is no adjacent face[1] then we know that
we are working with the first face of the initial tetrahedron.  In this
case we want to store the vertices in the opposite order from the
initial face.  Otherwise, we want to store the vertices in the same order
as in the visible face.  The third vertex is always p.
------------------------------------------------------------------------*/
void        MakeCcw( tFace f, tEdge e, tVertex p )
{
        int         i;        /* Index */
        tFace         fi;        /* The invisible face adjacent to e */
        tEdge        s;        /* Temporary, for swapping */

        /* If this is the initial tetrahedron, then e has only one
           adjacent face, and use that for fi.  Otherwise, use the
           invisible face. */
        if ( !e->adjface[1] )
                fi = e->adjface[0];
        else {
                if  ( !e->adjface[0]->visible )
                        fi = e->adjface[0];
                else        fi = e->adjface[1];
        }

        /* Set vertex[0] & [1] of f to have the opposite orientation
           as do the corresponding vertices of fi. */
        /* Find the index i of e->endpoint[1] in fi. */
        for ( i=0; fi->vertex[i] != e->endpts[1]; ++i )
                ;
        /* Orient f opposite that of fi. */
        if ( fi->vertex[ (i+1) % 3 ] != e->endpts[0] ) {
                f->vertex[0] = e->endpts[1];
                f->vertex[1] = e->endpts[0];
        }
        else {
                f->vertex[0] = e->endpts[0];
                f->vertex[1] = e->endpts[1];
                LSWAP( s, f->edge[1], f->edge[2] );
        }

        f->vertex[2] = p;
}

/*---------------------------------------------------------------------
MakeEdge creates a new cell and initializes all pointers to NULL
and sets all flags to off.  It returns a pointer to the empty cell.
---------------------------------------------------------------------*/
tEdge         MakeEdge( void )
{
        tEdge         e;

        NEW( e, tsEdge );
        e->adjface[0] = e->adjface[1] = e->newface = NULL;
        e->endpts[0] = e->endpts[1] = NULL;
        e->delete = !REMOVED;
        ADD( edges, e );
        return e;
}

/*---------------------------------------------------------------------
MakeFace creates a new face structure and initializes all of its
flags to NULL and sets all the flags to off.  It returns a pointer
to the empty cell.
----------------------------------------------------------------------*/
tFace         MakeFace( void )
{
        tFace         f;
        int         i;

        NEW( f, tsFace);
        for ( i=0; i < 3; ++i ) {
                  f->edge[i] = NULL;
                  f->vertex[i] = NULL;
        }
        f->visible = !VISIBLE;
        ADD( faces, f );
        return f;
}

/*-------------------------------------------------------------------------
CleanUp goes through each data structure list and clears all
flags and NULLs out some pointers.  The order of processing
(edges, faces, vertices) is important.
------------------------------------------------------------------------*/
void        CleanUp( void )
{
        CleanEdges();
        CleanFaces();
        CleanVertices();
}

/*------------------------------------------------------------------------
CleanEdges runs through the edge list and cleans up the structure.
If there is a newface then it will put that face in place of the
visible face and NULL out newface. It also deletes so marked edges.
-----------------------------------------------------------------------*/
void        CleanEdges( void )
{
        tEdge         e;        /* Primary index into edge list. */
        tEdge         t;        /* Temporary edge pointer. */

        /* Integrate the newface's into the data structure. */
        /* Check every edge. */
        e = edges;
        do {
                if ( e->newface ) {
                        if ( e->adjface[0]->visible )
                                e->adjface[0] = e->newface;
                        else        e->adjface[1] = e->newface;
                        e->newface = NULL;
                }
                e = e->next;
        } while ( e != edges );

        /* Delete any edges marked for deletion. */
        while ( edges && edges->delete ) {
                e = edges;
                DELETE( edges, e );
        }
        e = edges->next;
        do {
                if ( e->delete ) {
                        t = e;
                        e = e->next;
                        DELETE( edges, t );
                }
                else e = e->next;
        } while ( e != edges );
}

/*------------------------------------------------------------------------
CleanFaces runs through the face list and deletes any face marked visible.
-----------------------------------------------------------------------*/
void        CleanFaces( void )
{
        tFace         f;        /* Primary pointer into face list. */
        tFace         t;        /* Temporary pointer, for deleting. */


        while ( faces && faces->visible ) {
                f = faces;
                DELETE( faces, f );
        }
        f = faces->next;
        do {
                if ( f->visible ) {
                        t = f;
                        f = f->next;
                        DELETE( faces, t );
                }
                else f = f->next;
        } while ( f != faces );
}

/*-------------------------------------------------------------------------
CleanVertices runs through the vertex list and deletes the
vertices that are marked as processed but are not incident to any undeleted
edges.
-------------------------------------------------------------------------*/
void        CleanVertices( void )
{
        tEdge         e;
        tVertex        v, t;

        /* Mark all vertices incident to some undeleted edge as on the hull. */
        e = edges;
        do {
                 e->endpts[0]->onhull = e->endpts[1]->onhull = ONHULL;
                 e = e->next;
        } while (e != edges);

        /* Delete all vertices that have been processed but
           are not on the hull. */
        while ( vertices && vertices->mark && !vertices->onhull ) {
                v = vertices;
                DELETE( vertices, v );
        }
        v = vertices->next;
        do {
                if ( v->mark && !v->onhull ) {
                        t = v;
                        v = v->next;
                        DELETE( vertices, t )
                }
                else v = v->next;
        } while ( v != vertices );

        /* Reset flags. */
        v = vertices;
        do {
                v->duplicate = NULL;
                v->onhull = !ONHULL;
                v = v->next;
        } while ( v != vertices );
}
/*------------------------------------------------------------------
Collinear checks to see if the three points given are collinear,
by checking to see if each element of the cross product is zero.
---------------------------------------------------------------------*/
bool        Collinear( tVertex a, tVertex b, tVertex c )
{
        return ( c->v[LZD] - a->v[LZD] ) * ( b->v[LYD] - a->v[LYD] ) -
               ( b->v[LZD] - a->v[LZD] ) * ( c->v[LYD] - a->v[LYD] ) == 0
            && ( b->v[LZD] - a->v[LZD] ) * ( c->v[LXD] - a->v[LXD] ) -
               ( b->v[LXD] - a->v[LXD] ) * ( c->v[LZD] - a->v[LZD] ) == 0
            && ( b->v[LXD] - a->v[LXD] ) * ( c->v[LYD] - a->v[LYD] ) -
               ( b->v[LYD] - a->v[LYD] ) * ( c->v[LXD] - a->v[LXD] ) == 0  ;
}


/*------------------------------------------------------------------------
Consistency runs through the edge list and checks that all
adjacent faces have their endpoints in opposite order.  This verifies
that the vertices are in counterclockwise order.
-----------------------------------------------------------------------*/
void        Consistency( void )
{
        register tEdge e;
        register int i, j;

        e = edges;

        do {
               /* find index of endpoint[0] in adjacent face[0] */
           for ( i = 0; e->adjface[0]->vertex[i] != e->endpts[0]; ++i )
                 ;

              /* find index of endpoint[0] in adjacent face[1] */
              for ( j = 0; e->adjface[1]->vertex[j] != e->endpts[0]; ++j )
                  ;

           /* check if the endpoints occur in opposite order */
           if ( !( e->adjface[0]->vertex[ (i+1) % 3 ] ==
                   e->adjface[1]->vertex[ (j+2) % 3 ] ||
                   e->adjface[0]->vertex[ (i+2) % 3 ] ==
                   e->adjface[1]->vertex[ (j+1) % 3 ] )  )
               break;
           e = e->next;

           } while ( e != edges );

        if ( e != edges )
             fprintf( stderr, "Checks: edges are NOT consistent.\n");
        else
             fprintf( stderr, "Checks: edges consistent.\n");

}

/*----------------------------------------------------------------------
Convexity checks that the volume between every face and every
point is negative.  This shows that each point is inside every face
and therefore the hull is convex.
---------------------------------------------------------------------*/
void        Convexity( void )
{
        register tFace f;
        register tVertex v;
        float vol;

        f = faces;

        do {
              v = vertices;
              do {
                   if ( v->mark ) {
                      vol = Volume6( f, v );
                      if ( vol < 0 )
                        break;
                   }
                   v = v->next;
                   } while ( v != vertices );

              f = f->next;

            } while ( f != faces );

        if ( f != faces )
           fprintf( stderr, "Checks: NOT convex.\n");
        else if ( debug )
           fprintf( stderr, "Checks: convex.\n");
}

/*----------------------------------------------------------------------
CheckEuler checks Euler's relation, as well as its implications when
all faces are known to be triangles.  Only prints positive information
when debug is true, but always prints negative information.
  ---------------------------------------------------------------------*/
void        CheckEuler( int V, int E, int F )
{
        if ( debug )
             fprintf( stderr, "Checks: V, E, F = %d %d %d:\t", V, E, F);

        if ( (V - E + F) != 2 )
             fprintf( stderr, "Checks: V-E+F != 2\n");
        else if ( debug )
             fprintf( stderr, "V-E+F = 2\t");


        if ( F != (2 * V - 4) )
             fprintf( stderr, "Checks: F=%d != 2V-4=%d; V=%d\n",
                F, 2*V-4, V);
        else if ( debug )
             fprintf( stderr, "F = 2V-4\t");

        if ( (2 * E) != (3 * F) )
             fprintf( stderr, "Checks: 2E=%d != 3F=%d; E=%d, F=%d\n",
                2*E, 3*F, E, F );
        else if ( debug )
             fprintf( stderr, "2E = 3F\n");
}
/*-----------------------------------------------------------------------*/
void        Checks( void )
{
        tVertex v;
        tEdge   e;
        tFace   f;
        int         V = 0, E = 0 , F = 0;

        Consistency();
        Convexity();
        if ( v = vertices )
             do {
                if (v->mark) V++;
                v = v->next;
                } while ( v != vertices );
        if ( e = edges )
             do {
                 E++;
                e = e->next;
             } while ( e != edges );
        if ( f = faces )
             do {
                F++;
                f  = f ->next;
                } while ( f  != faces );
        CheckEuler( V, E, F );
}
/*================================================================
These functions are used whenever the debug flag is set.
They print out the entire contents of each data structure.
Printing is to standard error.  To grab the output in a file in the csh,
use this:
        chull < i.file >&! o.file
==================================================================*/
/*-------------------------------------------------------------------*/
void        PrintOut( tVertex v )
{
        fprintf( stderr, "\nAdding vertex %6x :\n", v );
        PrintVertices();
        PrintEdges();
        PrintFaces();
}

/*-------------------------------------------------------------------------*/
void        PrintVertices( void )
{
        tVertex temp;

        temp = vertices;
        fprintf (stderr, "Vertex List\n");
        if (vertices) do {
            fprintf(stderr,"  addr %6x\t", vertices );
            fprintf(stderr,"  vnum %4d", vertices->vnum );
            fprintf(stderr,"   (%6d,%6d,%6d)",vertices->v[LXD],
                    vertices->v[LYD], vertices->v[LZD] );
            fprintf(stderr,"   active:%3d", vertices->onhull );
            fprintf(stderr,"   dup:%5x", vertices->duplicate );
            fprintf(stderr,"   mark:%2d\n", vertices->mark );
            vertices = vertices->next;
            } while ( vertices != temp );

}
/*-------------------------------------------------------------------------*/
void        PrintEdges( void )
{
        tEdge         temp;
        int         i;

        temp = edges;
        fprintf (stderr, "Edge List\n");
        if (edges) do {
            fprintf( stderr, "  addr: %6x\t", edges );
            fprintf( stderr, "adj: ");
            for (i=0; i<2; ++i)
                 fprintf( stderr, "%6x", edges->adjface[i] );
            fprintf( stderr, "  endpts:");
            for (i=0; i<2; ++i)
                 fprintf( stderr, "%4d", edges->endpts[i]->vnum);
            fprintf( stderr, "  del:%3d\n", edges->delete );
            edges = edges->next;
        } while (edges != temp );

}
/*-------------------------------------------------------------------------*/
void        PrintFaces( void )
{
        int         i;
        tFace         temp;

        temp = faces;
        fprintf (stderr, "Face List\n");
        if (faces) do {
            fprintf(stderr, "  addr: %6x\t", faces );
            fprintf(stderr, "  edges:");
            for( i=0; i<3; ++i )
                 fprintf(stderr, "%6x", faces->edge[i] );
            fprintf(stderr, "  vert:");
            for ( i=0; i<3; ++i)
                  fprintf(stderr, "%4d", faces->vertex[i]->vnum );
            fprintf(stderr, "  vis: %d\n", faces->visible );
            faces= faces->next;
        } while ( faces != temp );

}
