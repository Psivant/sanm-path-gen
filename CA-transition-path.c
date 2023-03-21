/*
Released under MIT License

Copyright (c) 2023 Psivant Therapeutics, LLC.

Copyright (c) 2023 Istvan B Kolossvary.

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

#include <float.h>      /* For DBL_EPSILON */
#include <math.h>       /* For sqrt, fabs, atan2 */
#include <stdio.h>      /* For fprintf, sprintf, sscanf, fopen, fclose, fgets, ferror, feof */
#include <stdlib.h>     /* For atoi, atof */
#include <string.h>     /* For memcpy, strcmp, strncmp, strcpy, strstr */
#include <assert.h>     /* For assert */

#include "common.h"

/* These macros are not robust: */
#define SQR(a)          ( (a)*(a) )
#define ABS(a)          ( (a)>=0  ? (a) :(-a) )
#define MAX(a,b)        ( (a)>(b) ? (a) : (b) )
#define MIN(a,b)        ( (a)<(b) ? (a) : (b) )

#define PI              3.141592653589793
#define ZERO            0.0
#define ONE             1.0
#define TRUE            1
#define FALSE           0
#define YES             1
#define NO              0
#define MAXBUF          1024    /* char array buffer */
#define MAXFRAMES       1000    /* max number of path frames */
/* M-LBFGS general params: */
#define M_LBFGS         3
#define ITER_MAX        20000
#define GRMS_MAX        5e-02
#define VERBOSITY       0
#define VISUAL_DEBUG    0
#define SAVE_PATH       0
/* Options for updating L-BFGS matrix: */
#define UNITMATRIX      1
#define SCALING         2
#define DIAGONAL        3
/* MCSRCH line minimization: */
#define GTOL            9e-01
#define FTOL            1e-04
#define MAXFEV          10
#define MAXCRDMOV       0.1
/* ANM model parameters: */
#define NPOT            6       /* Number of potentials used in the model */
#define ANM_CUTOFF      10.0    /* CA-CA distance cutoff in the ANM model */
#define MAX_EN          11      /* Dimesion of the energy[] array */

# define DOT(n,v1,v2,dot)           \
BEGIN_PRIMITIVE {                   \
  int _n=(n),_i;                    \
  double * __restrict _v1=(v1);     \
  double * __restrict _v2=(v2);     \
  double _dot=0.0;                  \
  for( _i=0; _i<_n; _i++ )          \
    _dot += _v1[_i] * _v2[_i];      \
  (dot) = _dot;                     \
} END_PRIMITIVE

# define NORM2(n,v,norm)            \
BEGIN_PRIMITIVE {                   \
  int _n=(n),_i;                    \
  double * __restrict _v=(v);       \
  double _norm=0.0;                 \
  for( _i=0; _i<_n; _i++ )          \
    _norm += _v[_i] * _v[_i];       \
  (norm) = sqrt( _norm );           \
} END_PRIMITIVE

# define NORMINF(n,v,norm)          \
BEGIN_PRIMITIVE {                   \
  int _n=(n),_i;                    \
  double * __restrict _v=(v);       \
  double _norm=0.0, _t;             \
  for( _i=0; _i<_n; _i++ )          \
    if( (_t=fabs(_v[_i])) >=_norm ) \
      _norm = _t;                   \
  (norm) = _norm;                   \
} END_PRIMITIVE

# define NORMINF_XYZ(n,v,norm)      \
BEGIN_PRIMITIVE {                   \
  int _n=(n),_i;                    \
  double * __restrict _v=(v);       \
  double _norm=0.0, _t;             \
  for( _i=0; _i<_n; _i++ )      {   \
    _t=(_v[_i*3  ]*_v[_i*3  ] +     \
        _v[_i*3+1]*_v[_i*3+1] +     \
        _v[_i*3+2]*_v[_i*3+2] );    \
    if( _t > _norm ) _norm = _t;}   \
  (norm) = sqrt(_norm);             \
} END_PRIMITIVE


typedef struct min_args_t {

  int iter_in;
  int iter_out;
  double grms_in;
  double grms_out;
  double rmsd_gap;
  int m_lbfgs;
  int verbosity;
  int visual_debug;
  int save_path;
  int which_path;
  int n_path_frame[2];
  double ** path;
} min_args_t;

typedef struct lbfgs_t {

  double rho;
  double gamma;
  double * s;
  double * y;
} lbfgs_t;
    
typedef struct func_args_t {

  double fc_str;
  double r0_str;
  double fbhw_str;
  double fc_ang;
  double t0_ang;
  double fbhw_ang;
  double k0_ang;
  double fc_dih;
  double p0_dih;
  double fbhw_dih;
  double k0_dih;
  double fc_bump;
  double dmin_bump;
  double fc_posre;
  double * posre_xyz;   /* allow for position restraints */
  double fc_sanm;
  double * min1_xyz;    /* start of path */
  double * min2_xyz;    /* end of path */
  double * r0_min1;
  double * r0_min2;
  char   * pdb_file_image;
} func_args_t;

typedef error_msg
(*eval_derivs_t)( int n,
                  double * x,
                  double * fx,  /* array of energy terms, fx[0] = total energy */
                  double * grad,
                  void * args );


/*****
        R M S D  F I T
*****/
error_msg
diag( int n,
      double aa[3][3],
      double * d,
      double x[3][3] ) {

  int i, j, k, l, np, n2, j1;
  double *e=NULL, *a=NULL, ha, ep, tl, h, g, f, s, h1, ze, b, el, r, p, dx, c, ei;

  MALLOC( e, n   );
  MALLOC( a, n*n );

  for (i = 0; i < n; i++)
     for (j = 0; j < n; j++)
        a[i * n + j] = aa[i][j];
  tl = 1.e-29;
  ha = .5;
  np = n - 1;
  n2 = np;
  ep = 1.e-08;
  ze = ZERO;
  for (i = 0; i <= np; i++)
     for (j = i + 1; j <= np; j++)
        a[i * n + j] = ZERO;
  if (np > 1)
     for (i = np; i >= 2; i--) {
        l = i - 2;
        h = 0;
        g = a[i * n + i - 1];
        for (k = 0; k <= l; k++) {
           f = a[i * n + k];
           h += f * f;
        }
        if (h > 0) {
           s = h + g * g;
           if (s >= tl) {
              l++;
              f = g;
              g = sqrt(s);
              g = (g + s / g) * ha;
              if (f > 0.)
                 g = -g;
              h = s - f * g;
              a[i * n + i - 1] = f - g;
              f = 0.;
              h1 = 1. / h;
              for (j = 0; j <= l; j++) {
                 a[j * n + i] = a[i * n + j] * h1;
                 s = 0.;
                 for (k = 0; k <= j; k++)
                    s += a[j * n + k] * a[i * n + k];
                 j1 = j + 1;
                 if (j1 <= l)
                    for (k = j1; k <= l; k++)
                       s += a[k * n + j] * a[i * n + k];
                 e[j] = s * h1;
                 f += s * a[j * n + i];
              }
              h1 *= f * ha;
              for (j = 0; j <= l; j++) {
                 f = a[i * n + j];
                 s = e[j];
                 s -= h1 * f;
                 e[j] = s;
                 for (k = 0; k <= j; k++)
                    a[j * n + k] -= f * e[k] + a[i * n + k] * s;
              }
           }
        }
        d[i] = h;
        e[i - 1] = g;
     }
  e[0] = a[n];
  d[0] = a[0];
  a[0] = 1.;
  for (i = 1; i <= np; i++) {
     l = i - 1;
     if (d[i] > 0.) {
        for (j = 0; j <= l; j++) {
           s = 0.0;
           for (k = 0; k <= l; k++)
              s += a[k * n + j] * a[i * n + k];
           for (k = 0; k <= l; k++)
              a[k * n + j] -= s * a[k * n + i];
        }
     }
     d[i] = a[i * (n + 1)];
     a[i * (n + 1)] = 1.;
     for (j = 0; j <= l; j++)
        a[i * n + j] = a[j * n + i] = ze;
  }
  b = 0.0;
  f = 0.0;
  e[np] = 0.0;
  for (l = 0; l <= np; l++) {
     el = e[l];
     r = fabs(el);
     h = ep * (fabs(d[l]) + r);
     if (h > b)
        b = h;
     if (r > b) {
        for (i = l + 1; i <= np; i++)
           if (fabs(e[i]) <= b) {
              j = i;
              i = np + 1;
           }
        do {
           p = (d[l + 1] - d[l]) * ha / el;
           dx = p * p + 1.;
           r = sqrt(dx);
           r = (r + dx / r) * ha;
           if (p < 0.0)
              p -= r;
           else
              p += r;
           h = d[l] - el / p;
           for (i = l; i <= np; i++)
              d[i] -= h;
           f += h;
           p = d[j];
           c = 1.;
           s = 0.0;
           j1 = j - 1;
           for (i = j1; i >= l; i--) {
              ei = e[i];
              g = c * ei;
              h = c * p;
              if (fabs(p) < fabs(ei)) {
                 c = p / ei;
                 dx = c * c + 1.;
                 r = sqrt(dx);
                 r = (r + dx / r) * ha;
                 e[i + 1] = s * ei * r;
                 s = 1. / r;
                 c /= r;
              } else {
                 c = ei / p;
                 dx = c * c + 1.;
                 r = sqrt(dx);
                 r = (r + dx / r) * ha;
                 e[i + 1] = s * p * r;
                 s = c / r;
                 c = 1. / r;
              }
              p = c * d[i] - s * g;
              d[i + 1] = h + s * (c * g + s * d[i]);
              for (k = 0; k <= np; k++) {
                 h = a[k * n + i + 1];
                 g = a[k * n + i];
                 a[k * n + i + 1] = g * s + h * c;
                 a[k * n + i] = g * c - h * s;
              }
           }
           el = s * p;
           e[l] = el;
           d[l] = c * p;
        } while (fabs(el) > b);
     }
     d[l] += f;
  }
  for (i = 0; i <= n2; i++) {
     k = i;
     p = d[i];
     for (j = i + 1; j <= np; j++)
        if (d[j] > p) {
           k = j;
           p = d[j];
        }
     if (k != i) {
        d[k] = d[i];
        d[i] = p;
        for (j = 0; j <= np; j++) {
           p = a[j * n + i];
           a[j * n + i] = a[j * n + k];
           a[j * n + k] = p;
        }
     }
  }
  for (i = 0; i < n; i++)
     for (j = 0; j < n; j++)
        x[i][j] = a[i * n + j];

  FREE(e);
  FREE(a);

  return SUCCESS;
}


error_msg
rmsfit( int n_atom,
        double * ref,
        double * xyz,
        double * rms,
        int inplace ) {
/*
   Based on  W. Kabsch; Acta Cryst. (1976). A32, 922-923,
             W. Kabsch; Acta Cryst. (1978). A34, 827-828.
*/

  int i, j, n, n3;
  double rmatrx[3][3], rtrmat[3][3], d[3], z[3][3], b[3][3], cntr_ref[3], cntr_xyz[3], umat[3][3];
  double fac, xsq, ysq, zsq, dsq, dmax, savex, savey, savez;
  double *save_ref=NULL, *save_xyz=NULL;

  if( ref==NULL ||
      xyz==NULL ) return ERROR_MSG( "Bad input arguments" );

  MALLOC( save_ref, 3*n_atom );
  MALLOC( save_xyz, 3*n_atom );

  COPY( save_ref, ref, 3*n_atom );
  COPY( save_xyz, xyz, 3*n_atom );

  for( i = 0; i < 3; i++ )
    cntr_ref[i] = cntr_xyz[i] = ZERO;
  for( n = 0; n < n_atom; n++) {
    n3 = 3*n;
    cntr_ref[0] += ref[n3  ];
    cntr_ref[1] += ref[n3+1];
    cntr_ref[2] += ref[n3+2];
    cntr_xyz[0] += xyz[n3  ];
    cntr_xyz[1] += xyz[n3+1];
    cntr_xyz[2] += xyz[n3+2];
  }
  cntr_ref[0] /= n_atom;
  cntr_ref[1] /= n_atom;
  cntr_ref[2] /= n_atom;
  cntr_xyz[0] /= n_atom;
  cntr_xyz[1] /= n_atom;
  cntr_xyz[2] /= n_atom;
  for (n = 0; n < n_atom; n++) {
     n3 = 3*n;
     ref[n3  ] -= cntr_ref[0];
     ref[n3+1] -= cntr_ref[1];
     ref[n3+2] -= cntr_ref[2];
     xyz[n3  ] -= cntr_xyz[0];
     xyz[n3+1] -= cntr_xyz[1];
     xyz[n3+2] -= cntr_xyz[2];
  }
  for( j = 0; j < 3; j++ )
    for( i = 0; i < 3; i++ )
      rmatrx[i][j] = ZERO;
  for( n = 0; n < n_atom; n++ ) {
    n3 = 3*n;
    rmatrx[0][0] += 0.1 * ref[n3  ] * xyz[n3  ];
    rmatrx[1][0] += 0.1 * ref[n3+1] * xyz[n3  ];
    rmatrx[2][0] += 0.1 * ref[n3+2] * xyz[n3  ];
    rmatrx[0][1] += 0.1 * ref[n3  ] * xyz[n3+1];
    rmatrx[1][1] += 0.1 * ref[n3+1] * xyz[n3+1];
    rmatrx[2][1] += 0.1 * ref[n3+2] * xyz[n3+1];
    rmatrx[0][2] += 0.1 * ref[n3  ] * xyz[n3+2];
    rmatrx[1][2] += 0.1 * ref[n3+1] * xyz[n3+2];
    rmatrx[2][2] += 0.1 * ref[n3+2] * xyz[n3+2];
  }
  for( j = 0; j < 3; j++ )
    for( i = 0; i < 3; i++ )
      rtrmat[j][i] = rmatrx[0][i] * rmatrx[0][j]
                   + rmatrx[1][i] * rmatrx[1][j]
                   + rmatrx[2][i] * rmatrx[2][j];
  TRY( diag( 3, rtrmat, d, z ) );
  z[0][2] = z[1][0] * z[2][1] - z[2][0] * z[1][1];
  z[1][2] = z[2][0] * z[0][1] - z[0][0] * z[2][1];
  z[2][2] = z[0][0] * z[1][1] - z[1][0] * z[0][1];
  for( j = 0; j < 2; j++ ) {
    for( i = 0; i < 3; i++ )
      b[i][j] = rmatrx[i][0] * z[0][j]
              + rmatrx[i][1] * z[1][j]
              + rmatrx[i][2] * z[2][j];
    fac = sqrt(b[0][j] * b[0][j] + b[1][j] * b[1][j] + b[2][j] * b[2][j]);
    for( i = 0; i < 3; i++ )
      b[i][j] /= MAX(fac, 1e-15);
  }
  b[0][2] = b[1][0] * b[2][1] - b[2][0] * b[1][1];
  b[1][2] = b[2][0] * b[0][1] - b[0][0] * b[2][1];
  b[2][2] = b[0][0] * b[1][1] - b[1][0] * b[0][1];
  for( j = 0; j < 3; j++ )
    for( i = 0; i < 3; i++ )
      umat[i][j] = b[i][0] * z[j][0]
                 + b[i][1] * z[j][1]
                 + b[i][2] * z[j][2];
  savex = cntr_xyz[0];
  savey = cntr_xyz[1];
  savez = cntr_xyz[2];
  for( n = 0, dsq = ZERO; n < n_atom; n++ ) {
    n3 = 3*n;
    savex = xyz[n3  ];
    savey = xyz[n3+1];
    savez = xyz[n3+2];
    xyz[n3  ] = umat[0][0] * savex + umat[0][1] * savey + umat[0][2] * savez;
    xsq = SQR( xyz[n3  ] - ref[n3  ] );
    dsq += xsq;
    xyz[n3+1] = umat[1][0] * savex + umat[1][1] * savey + umat[1][2] * savez;
    ysq = SQR( xyz[n3+1] - ref[n3+1]);
    dsq += ysq;
    xyz[n3+2] = umat[2][0] * savex + umat[2][1] * savey + umat[2][2] * savez;
    zsq = SQR( xyz[n3+2] - ref[n3+2]);
    dsq += zsq;
  }
  rms[0] = sqrt(dsq / n_atom);

  /* After optimal rotation translate back everything to the original reference frame */
  for (n = 0; n < n_atom; n++) {
     n3 = 3*n;
     ref[n3  ] += cntr_ref[0];
     ref[n3+1] += cntr_ref[1];
     ref[n3+2] += cntr_ref[2];
     xyz[n3  ] += cntr_ref[0];
     xyz[n3+1] += cntr_ref[1];
     xyz[n3+2] += cntr_ref[2];
  }

  if( inplace ) {  /* move ref and xyz structures back to their original position */
  
    COPY( ref, save_ref, 3*n_atom );
    COPY( xyz, save_xyz, 3*n_atom );
  }
  FREE( save_ref );
  FREE( save_xyz );

  return SUCCESS;
}


/*****
        S T R U C T U R E  I / O
*****/
error_msg
read_pdb_file( int n_atom,
               void * __restrict args,
               double * x,
               FILE * pdb_file ) {

  func_args_t * __restrict func_args = (func_args_t *) args;
  int i;
  char line[MAXBUF];
  error_msg err = SUCCESS;

  i = 0;
  while( fgets( line, MAXBUF, pdb_file ) != NULL ) {

    if( strstr( line, " CA") == NULL )
      continue;
    sscanf( line+30, "%lf %lf %lf", &x[3*i], &x[3*i+1], &x[3*i+2] );
    if( ++i >= n_atom )
      break;
  }

  if( ferror( pdb_file) ) {

    err = ERROR_MSG( "File read error" );
    goto cleanup;
  }

 cleanup:

  return err;
}


error_msg
write_pdb_file( int n,  /* n_atom */
                void * __restrict args,
                double * x,
                FILE * pdb_file,
                int format ) {  /* write pdb_file_image updated with current coordinates */

  func_args_t * __restrict func_args = (func_args_t *) args;
  int i;
  char xyz_field[25], line[MAXBUF];
  error_msg err = SUCCESS;

  if( format==1 ) fprintf( pdb_file, "MODEL\n" );  /* normal */
  for( i=0; i<n; i++ ) {

    sprintf( xyz_field, "%8.3f%8.3f%8.3f", x[i*3], x[i*3+1], x[i*3+2] );
    strcpy( line, &func_args->pdb_file_image[i*MAXBUF] );
    strncpy( &line[30], xyz_field, 24 );
    if( format==2 )
      strncpy( &line[54], "  1.00  1.00", 12 ); /* over write occupancy and b-factor fields with PLUMED RMSD weights */
    fprintf( pdb_file, "%s", line );
  }
  if( format==1 )
    fprintf( pdb_file, "ENDMDL\n" );  /* normal */
  else fprintf( pdb_file, "END\n" );  /* plumed */

 cleanup:

  fflush( pdb_file );

  return err;
}


error_msg
read_reference_structures( int n_atom,
                           void * __restrict args,
                           FILE * pdb_file ) {  /* read first and last frame from input path */

  func_args_t * __restrict func_args = (func_args_t *) args;
  int i, j, do_transform=0, in_place=1;
  double dist, rmsd[1];
  double * x;
  char line[MAXBUF];
  error_msg err = SUCCESS;

  x = func_args->min1_xyz;  /* read start structure */

  rewind( pdb_file ); /* just in case */
  i = 0;
  while( fgets( line, MAXBUF, pdb_file ) != NULL ) {

    if( strstr( line, " CA") == NULL )
      continue;

    sscanf( line+30, "%lf %lf %lf", &x[3*i], &x[3*i+1], &x[3*i+2] );        /* read coordinates */
    strcpy( &func_args->pdb_file_image[i*MAXBUF], line );                   /* create memory image of the PDB file */

    if( ++i >= n_atom )
      break;
  }

  if( ferror( pdb_file) ) {

    err = ERROR_MSG( "File read error" );
    goto cleanup;
  }

  for( i=1; i<n_atom; i++ )
  for( j=0; j<i;      j++ ) {

    dist = sqrt( SQR(x[3*i  ] - x[3*j  ]) + 
                 SQR(x[3*i+1] - x[3*j+1]) + 
                 SQR(x[3*i+2] - x[3*j+2]) );
    func_args->r0_min1[i*n_atom+j] = 
    func_args->r0_min1[j*n_atom+i] = dist;
  }

  x = func_args->min2_xyz;  /* read end structure */

  i = 0;
  while( fgets( line, MAXBUF, pdb_file ) != NULL ) {

    if( strstr( line, " CA") == NULL )
      continue;

    sscanf( line+30, "%lf %lf %lf", &x[3*i], &x[3*i+1], &x[3*i+2] );    /* read coordinates */

    if( ++i >= n_atom )
      i = 0;    /* reset counter */
  }

  if( ferror( pdb_file) ) {

    err = ERROR_MSG( "File read error" );
    goto cleanup;
  }

  /* Superimpose end structure to start structure: */
  TRY( rmsfit( n_atom, func_args->min1_xyz, x, rmsd, do_transform ) );

  for( i=1; i<n_atom; i++ )
  for( j=0; j<i;      j++ ) {

    dist = sqrt( SQR(x[3*i  ] - x[3*j  ]) + 
                 SQR(x[3*i+1] - x[3*j+1]) + 
                 SQR(x[3*i+2] - x[3*j+2]) );
    func_args->r0_min2[i*n_atom+j] = 
    func_args->r0_min2[j*n_atom+i] = dist;
  }

cleanup:

  return err;
}


/*****
        D E R I V S
*****/
error_msg
calc_derivs( eval_derivs_t * func_list,
             int n,  /* n_dim */
             double * x,
             double * fx,
             double * grad,
             void * func_args ) {

  int i;

  if( n<=0 ) return ERROR_MSG( "Bad input arguments" );
  if( x==NULL  || 
      fx==NULL ||
      func_list==NULL ) {
    return ERROR_MSG( "Bad input arguments" );
  }
  else {
    CLEAR( fx, MAX_EN );        /* clear function value and components */
  }
  if( grad==NULL ) {
    return ERROR_MSG( "Bad input arguments" );
  }
  else {
    CLEAR( grad, n );           /* clear gradient */
  }

  for( i=0; func_list[i] != NULL; i++ ) {  /* reduce grad for a list of functions */

    func_list[i]( n,
                  x,
                  fx,
                  grad,
                  func_args );
  }

  return SUCCESS;
}


error_msg
calc_sanm_derivs( int n,  /* n_dim */
                  double * __restrict x,
                  double * __restrict fx,
                  double * __restrict grad,
                  void * __restrict args ) {  /* This is the Smooth Anisotropic Network Model potential */

  static unsigned long int select=0;  /* two consecutive calls to this function will use min1 and min2, respectively */

  func_args_t * __restrict func_args = (func_args_t *) args;
  int N = n/3; /* n_atom */
  int i, j, i3, j3, k;  
  double en, eni, enij;
  double r2, one_r, r0, dr0, dr02, tmp, test1, test2;
  double dx, dy, dz, rix, riy, riz, gx, gy, gz, gix, giy, giz, g_r;
  double fc=func_args->fc_sanm;
  double * __restrict r0_min1 = func_args->r0_min1,
         * __restrict r0_min2 = func_args->r0_min2;  /* pairwise CA-CA distances in ref structures */

  ++select;  /* odd computes w/ min1 and even computes w/ min2 */

  en = 0;
  
  for( i=1; i<N; i++ ) {

    i3 = 3*i;
    
    rix = x[i3  ];
    riy = x[i3+1];
    riz = x[i3+2];

    eni = 0;
    
    gix  = 0;
    giy  = 0;
    giz  = 0;
    
    for( j=0; j<i; j++ ) {

      k = i*N+j;
      test1 = 0.5 * (r0_min1[k] + r0_min2[k]);

      if( test1 > ANM_CUTOFF )         continue;  /* residues outside the range of the ANM model */
      if( j == (i-1) && test1 <= 4.2 ) continue;  /* skip 1-2 interactions unless there is a gap */
      if( j == (i-2) ) {                          /* skip 1-3 interactions unless there is a gap */

        k = i*(N+1)-1;
        test1 = 0.5 * (r0_min1[k] + r0_min2[k]);
        k -= N+1;
        test2 = 0.5 * (r0_min1[k] + r0_min2[k]);
        if( test1 < 4.2 && test2 < 4.2 ) continue;
      }

      j3 = 3*j;
 
      dx  = rix - x[j3  ];
      dy  = riy - x[j3+1];
      dz  = riz - x[j3+2];
      
      r2 = dx*dx + (dy*dy + dz*dz);     /* current CA-CA distance */

      if( select %2 == 1 ) 
        r0 = r0_min1[i*N+j];            /* CA-CA distance in min1 */
      else
        r0 = r0_min2[i*N+j];            /* CA-CA distance in min2 */
           
      one_r  = 1.0/sqrt(r2);
      dr0    = r2*one_r - r0;
      dr02   = dr0*dr0;
      tmp    = 1.0/(1.0+dr02);
      enij   = 0.5*fc*dr02*tmp;

      eni += enij;
           
      g_r = fc*dr0*tmp*tmp*one_r;

      gx = dx*g_r;  // grad wrt 'x', dr/dx=dx/r, but 1/r is absorbed in g_r
      gy = dy*g_r;
      gz = dz*g_r;
                      
      gix        += gx;
      giy        += gy;
      giz        += gz;
      
      grad[j3  ] -= gx;
      grad[j3+1] -= gy;
      grad[j3+2] -= gz;
      
    }  /* min2 j */
    
    en += eni;

    grad[i3  ]      += gix;
    grad[i3+1]      += giy;
    grad[i3+2]      += giz;

  }  /* min2 i */
  
  fx[0] += en;  /* reduce energy */
  fx[6] += en;
  if( select %2 == 1 ) fx[7] += en; else fx[8] += en;

  return SUCCESS;
}


error_msg
calc_angle_derivs( int n,  /* n_dim */
                   double * __restrict x,
                   double * __restrict fx,
                   double * __restrict grad,
                   void * __restrict args ) {

  func_args_t * __restrict func_args = (func_args_t *) args;
  int N = n/3; /* n_atom */
  double * __restrict r0_min1 = func_args->r0_min1,
         * __restrict r0_min2 = func_args->r0_min2;  /* pairwise CA-CA distances in ref structures */

  const double tiny = 1e-15;
  int i,j,k,i3,j3,k3,l,m;                   // Particle ids
  double eni, en;                           // Interaction energy
  double eni0, en0;                         // Penalty term
  double fc    =func_args->fc_ang;          // Force constant
  double t0    =func_args->t0_ang;          // Flat-bottom well center
  double sigma =func_args->fbhw_ang;        // Flat-bottom half-width
  double k0    =func_args->k0_ang;          // Penalty force constant
  double fxi, fyi, fzi;                     // Force on the i particle
  double fxk, fyk, fzk;                     // Force on the k particle
  double v0, v1, v2, v3, v4;                // Computation temporaries
  int c0;

// Algorithm by Kevin J. Bowers
// The following algorithm reduces the forces and potential energy
// on particles associated with generic angle bmin2ing terms.
// Roughly, in numerically stabilized form:
//
// xij  = x[i]-y[j],        yij  = y[i]-y[j],        zij  = z[i]-z[j]
// rij = sqrt( xij^2 + yij^2 + zij^2 + tiny )
// xij /= rij,              yij /= rij,              zij /= rij
//
// xkj  = x[k]-y[j],        ykj  = y[k]-y[j],        zkj  = z[k]-z[j]
// rkj = sqrt( xkj^2 + ykj^2 + zkj^2 + tiny )
// xkj /= rkj,              ykj /= rkj,              zkj /= rkj
//
// ax   = xij + xkj,        ay   = yij + ykj,        az   = zij + zkj
// bx   = xij - xkj,        by   = yij - ykj,        bz   = zij - zkj
// if xij*xkj + yij*ykj + zij*zkj >= 0, 
//   c   = sqrt( ax^2 + ay^2 + az^2 + tiny )
//   ax /= c,    ay /= c,    az /= c
//   c  /= 2
//   d   = ax*bx + ay*by + az*bz
//   bx -= d*ax, by -= d*ay, bz -= d*az
//   s   = sqrt( bx^2 + by^2 + bz^2 + tiny )
//   bx /= s,    by /= s,    bz /= s
//   s  /= 2
// else
//   s   = sqrt( bx^2 + by^2 + bz^2 + tiny )
//   bx /= s,    by /= s,    bz /= s
//   s  /= 2
//   d   = ax*bx + ay*by + az*bz
//   ax -= d*bx, ay -= d*by, az -= d*bz
//   c   = sqrt( ax^2 + ay^2 + az^2 + tiny )
//   ax /= c,    ay /= c,    az /= c
//   c  /= 2
// min2if
// ax *= s, ay *= s, az *= s
// bx *= c, by *= c, bz *= c
// en    =  fc * c  (cos(phi))
// fr    = -fc * s  (sin(phi))
//
// fxi = fr(ax-bx)/|rij|, fyi = fr(ay-by)/|rij|, fzi = fr(az-bz)/|rij|
// fxk = fr(ax+bx)/|rkj|, fyk = fr(ay+by)/|rkj|, fzk = fr(az+bz)/|rkj|
//
// fx[i] += fxi,          fy[i] += fyi,          fz[i] += fzi
//
// fx[j] -= (fxi+fxk),    fy[j] -= (fyi+fyk),    fz[j] -= (fzi+fzk)
//
// fx[k] += fxk,          fy[k] += fyk,          fz[k] += fzk
//
// Angle bending terms have a fan in/out of three.
//
// Note that an additional flat-bottom penalty function will be
// applied to prevent bond angles to become too close to 180 deg
// which causes a physical divergence in the torsion angle term
// that cannot be helped numerically and would wreak havoc.
// Moreover, the flat bottom potential can be used to keep bond
// angles in a certain range.

  en = en0 = 0;
  
  for( i=0; i<N-2; i++ ) {  /* loop over all bond angles */

    j  = i+1;
    k  = i+2;
    l  = i*N+j;
    m  = j*N+k;
    if( (0.5 * (r0_min1[l] + r0_min2[l]) > 4.2) 
    ||  (0.5 * (r0_min1[m] + r0_min2[m]) > 4.2) ) continue;  /* one of them is a gap and not a bond */

    i3 = 3*i;
    j3 = 3*j;
    k3 = 3*k;
    
    // Swizzle in the particle data

    fxi = x[i3  ];
    fyi = x[i3+1];
    fzi = x[i3+2];

    v0  = x[j3  ];
    v1  = x[j3+1];
    v2  = x[j3+2];

    fxk = x[k3  ];
    fyk = x[k3+1];
    fzk = x[k3+2];
    
    // Compute bond vectors

    fxi -= v0;
    fyi -= v1;
    fzi -= v2;

    fxk -= v0;
    fyk -= v1;
    fzk -= v2;


    // Compute the force and the interaction energy

    v3   = 1/sqrt(fxi*fxi+fyi*fyi+fzi*fzi+tiny); // v3 = 1/|rij|
    fxi *= v3;
    fyi *= v3;
    fzi *= v3;                                   // fi = rij_hat

    v4   = 1/sqrt(fxk*fxk+fyk*fyk+fzk*fzk+tiny); // v4 = 1/|rkj|
    fxk *= v4;
    fyk *= v4;
    fzk *= v4;                                   // fk = rkj_hat

    c0   = (fxi*fxk+fyi*fyk+fzi*fzk)<0;          // true if a has cancellation

    v0   = fxi;
    v1   = fyi;
    v2   = fzi;
    fxi += fxk;
    fyi += fyk;
    fzi += fzk;                                  // fi = rij_hat+rkj_hat = a
    fxk  = v0-fxk;
    fyk  = v1-fyk;
    fzk  = v2-fzk;                               // fk = rij_hat-rkj_hat = b

    if(c0) {
      v0 = fxi, fxi = fxk, fxk = v0;
      v0 = fyi, fyi = fyk, fyk = v0;
      v0 = fzi, fzi = fzk, fzk = v0;             // fk is more inaccurate one
    }                                           

    v2   = fxi*fxi + fyi*fyi + fzi*fzi + tiny;   // v2 = |a|^2
    v0   = 1/sqrt(v2);                           // v0 = 1/|a|
    v2  *= v0/2;                                 // v2 = |a|/2 = cos(theta/2)
    fxi *= v0;
    fyi *= v0;
    fzi *= v0;                                   // fi = a_hat

    v0   = fxi*fxk + fyi*fyk + fzi*fzk;          // v0 = b.a_hat ~ 0

    fxk -= v0*fxi;
    fyk -= v0*fyi;
    fzk -= v0*fzi;                               // fk  = b - (b.a_hat)a_hat

    v1   = fxk*fxk + fyk*fyk + fzk*fzk + tiny;   // v1  = |b|^2
    v0   = 1/sqrt(v1);                           // v0  = 1/|b|
    v1  *= v0/2;                                 // v1  = |b|/2 = sin(theta/2)
    v0  *= v2;                                   // v0  = cos(theta/2)/|b|

    fxi *= v1;
    fyi *= v1;
    fzi *= v1;                                   // fi = a_hat sin(theta/2)

    fxk *= v0;
    fyk *= v0;
    fzk *= v0;                                   // fk = b_hat cos(theta/2)

    if(c0) {
      v0 = fxi, fxi = fxk, fxk = v0;
      v0 = fyi, fyi = fyk, fyk = v0;
      v0 = fzi, fzi = fzk, fzk = v0;             // get a_hat, b_hat in order
      v0 = v1,  v1  = v2,  v2  = v0;             // get sin/cos in order
    }

    eni  =  fc*(2*v2*v2-1);                      // en =  fc*cos(theta)
    v0   = -fc*2*v1*v2;                          // fr = -fc*sin(theta) = en'

    // Add the flat-bottom well penalty

    v1   = atan2( v1, v2 );                      // v1 = theta/2
    v1  += v1-t0;                                // v1 = theta-theta0

    // apply fbhw condition
    v1   = v1 < -sigma ? v1+sigma :
           v1 >  sigma ? v1-sigma :
           0;

    v2   = k0*v1;                                // v2 =  k0(theta-theta0)
    eni0 = v2*v1;                                // en =  k0(theta-theta0)^2
    v2  += v2;                                   // fr = 2k0(theta-theta0)

    v3  *= v0+v2;                                // v3 = fr/|rij|
    v4  *= v0+v2;                                // v4 = fr/|rkj|

    v0   = fxi;
    v1   = fyi;
    v2   = fzi;
    fxi  = v3*(fxi-fxk);
    fyi  = v3*(fyi-fyk);
    fzi  = v3*(fzi-fzk);                         // fi = fi
    fxk  = v4*(v0 +fxk);
    fyk  = v4*(v1 +fyk);
    fzk  = v4*(v2 +fzk);                         // fk = fk

    en  += eni;
    en0 += eni0;

    // Reduce the forces --- grad .eq. negative force

    grad[i3  ] -= fxi;
    grad[i3+1] -= fyi;
    grad[i3+2] -= fzi; 

    grad[j3  ] += fxi+fxk;
    grad[j3+1] += fyi+fyk;
    grad[j3+2] += fzi+fzk;

    grad[k3  ] -= fxk;
    grad[k3+1] -= fyk;
    grad[k3+2] -= fzk;

  }

  // Reduce the energy

  fx[0] += en+en0;
  fx[2] += en0;
  fx[3] += en;
    
  return SUCCESS;
}


error_msg
calc_bond_derivs( int n,  /* n_dim */
                  double * __restrict x,
                  double * __restrict fx,
                  double * __restrict grad,
                  void * __restrict args ) {

  func_args_t * __restrict func_args = (func_args_t *) args;
  int N = n/3; /* n_atom */
  double * __restrict r0_min1 = func_args->r0_min1,
         * __restrict r0_min2 = func_args->r0_min2;  /* pairwise CA-CA distances in ref structures */

  const double tiny = 1e-15;
  int i,j,i3,j3,k;                          // Particle ids
  double eni, en;                           // Interaction energy
  double fc=func_args->fc_str;              // Force constant
  double r0=func_args->r0_str;              // Bond length
  double sigma =func_args->fbhw_str;        // Flat-bottom half-width
  double fxi, fyi, fzi;                     // Force on the i particle
  double v0, v1, v2;                        // Computation temporaries

// Algorithm by Kevin J. Bowers
// The following algorithm reduces the forces and potential energy
// on particles associated with generic bond stretching terms.
// Roughly, in numerically stabilized form:
//
// xij = xi-xj
// yij = yi-yj
// zij = zi-zj
// rij = sqrt( xij^2 + yij^2 + zij^2 )
//
// fx = -2 fc ( rij - r0 ) xij / rij
// fy = -2 fc ( rij - r0 ) yij / rij
// fz = -2 fc ( rij - r0 ) zij / rij
// en =    fc ( rij - r0 )^2
//
// fxi += fx
// fyi += fy
// fzi += fz
//
// fxj -= fx
// fyj -= fy
// fzj -= fz
//
// Bond strecthing terms have a fan in/out of two.


  en = 0;
  
  for( i=0; i<N-1; i++ ) {  /* loop over all bonds */

    j  = i+1;
    k  = i*N+j;
    if( 0.5 * (r0_min1[k] + r0_min2[k]) > 4.2 ) continue;  /* this is a gap and not a bond */

    i3 = 3*i;
    j3 = 3*j;
    
    // Swizzle in the particle data

    fxi = x[i3  ];
    fyi = x[i3+1];
    fzi = x[i3+2];

    v0  = x[j3  ];
    v1  = x[j3+1];
    v2  = x[j3+2];

    // Compute bond vector

    fxi -= v0;
    fyi -= v1;
    fzi -= v2;

    // Compute the force and the interaction energy

    eni  = fxi*fxi + fyi*fyi + fzi*fzi + tiny; // en = rij^2
    v0   = 1/sqrt(eni);                        // v0 = 1/rij
    eni  = r0 - v0*eni;                        // en =    -( rij-r0 )

    // apply fbhw condition
    eni  = eni < -sigma ? eni+sigma :
           eni >  sigma ? eni-sigma :
           0;

    v1   = fc*eni;                             // v1 =  -fc( rij-r0 )
    eni *= v1;                                 // en =   fc( rij-r0 )^2
    v1  *= v0;                                 // v1 =  -fc( rij-r0 )/rij
    v1  += v1;                                 // v1 = -2fc( rij-r0 )
    fxi *= v1;                                 // fx = -2fc( rij-r0 )xij/rij
    fyi *= v1;                                 // fy = -2fc( rij-r0 )yij/rij
    fzi *= v1;                                 // fz = -2fc( rij-r0 )zij/rij

    en += eni;

    // Reduce the forces --- grad .eq. negative force

    grad[i3  ] -= fxi;
    grad[i3+1] -= fyi;
    grad[i3+2] -= fzi; 

    grad[j3  ] += fxi;
    grad[j3+1] += fyi;
    grad[j3+2] += fzi;

  }

  // Reduce the energy

  fx[0] += en;
  fx[1] += en;
    
  return SUCCESS;
}


error_msg
calc_posre_derivs( int n,  /* n_dim */
                   double * __restrict x,
                   double * __restrict fx,
                   double * __restrict grad,
                   void * __restrict args ) {  /* Pull one conformation to another conformation on the SANM surface */

  func_args_t * __restrict func_args = (func_args_t *) args;
  int N = n/3; /* n_atom */

  int i, i3;                                // Particle ids
  double eni, en;                           // Interaction energy
  double fc=func_args->fc_posre;            // Force constant
  double fxi, fyi, fzi;                     // Force on the i particle
  double v0, v1, v2;                        // Computation temporaries
  double *ref=func_args->posre_xyz;         // Reference conformation

  en = 0;
  
  for( i=0; i<N; i++ ) {  /* loop over all atoms */

    i3 = 3*i;
    
    // Swizzle in the particle data

    fxi = x[i3  ];
    fyi = x[i3+1];
    fzi = x[i3+2];

    v0  = ref[i3  ];
    v1  = ref[i3+1];
    v2  = ref[i3+2];

    // Compute distance vector

    fxi -= v0;
    fyi -= v1;
    fzi -= v2;

    // Compute the force and the interaction energy

    eni  = fxi*fxi + fyi*fyi + fzi*fzi;        // en = rij^2
    eni *= fc;                                 // en = fc*rij^2
    v0   = 2*fc;                               // v0 = 2fc
    fxi *= v0;                                 // fx = 2fc*rij*xij/rij
    fyi *= v0;                                 // fy = 2fc*rij*xij/rij
    fzi *= v0;                                 // fz = 2fc*rij*xij/rij

    en += eni;

    // Reduce the gradient

    grad[i3  ] += fxi;
    grad[i3+1] += fyi;
    grad[i3+2] += fzi; 

  }

  // Reduce the energy

  fx[0] += en;
  fx[9] += en;
    
  return SUCCESS;
}


error_msg
calc_bump_derivs( int n,  /* n_dim */
                  double * __restrict x,
                  double * __restrict fx,
                  double * __restrict grad,
                  void * __restrict args ) {  /* Avoid nonbonded contacts */

  func_args_t * __restrict func_args = (func_args_t *) args;
  int N = n/3; /* n_atom */
  int i, j, i3, j3;  
  double en, eni, enij;
  const double tiny=1e-15;
  double r2, one_r, dx, dy, dz, rix, riy, riz, gx, gy, gz, gix, giy, giz, g_r;
  double fc=func_args->fc_bump, dmin=func_args->dmin_bump, dmin2;

  dmin2 = dmin*dmin;
  en = 0;
  
  for( i=3; i<N; i++ ) {

    i3 = 3*i;
    
    rix = x[i3  ];
    riy = x[i3+1];
    riz = x[i3+2];

    eni = 0;
    
    gix  = 0;
    giy  = 0;
    giz  = 0;
    
    for( j=0; j<i-2; j++ ) {  /* 1-2 and 1-3 neighbors are skipped */

      j3 = 3*j;
 
      dx  = rix - x[j3  ];
      dy  = riy - x[j3+1];
      dz  = riz - x[j3+2];
      
      r2  = dx*dx + dy*dy + dz*dz + tiny;

      if( r2 > dmin2 ) continue;

      one_r = 1./sqrt(r2);
      enij  = one_r*r2 - dmin;  // enij = rij-rmin
      g_r   = fc*enij;          // g_r  = fc*(rij-rmin)
      enij *= g_r;              // enij = fc*(rij-rmin)^2
      
      eni += enij;
           
      g_r *= one_r;             // g_r  = g_r/rij
      g_r += g_r;               // g_r  = 2*fc*(rij-rmin)/rij

      gx = dx*g_r;              // grad wrt 'x', dr/dx=dx/r, but 1/r is absorbed in g_r
      gy = dy*g_r;
      gz = dz*g_r;
                      
      gix        += gx;
      giy        += gy;
      giz        += gz;
      
      grad[j3  ] -= gx;
      grad[j3+1] -= gy;
      grad[j3+2] -= gz;
      
    }  /* end j */
    
    en += eni;

    grad[i3  ]      += gix;
    grad[i3+1]      += giy;
    grad[i3+2]      += giz;

  }  /* end i */
  
  fx[0]  += en;  /* reduce energy */
  fx[10] += en;

  return SUCCESS;
}


double
calc_num_derivs( eval_derivs_t * obj_func,
                 int n,  /* ndim */
                 double * x,
                 double * fx,
                 double * grad,
                 void * func_args,
                 int i ) {  /* only for debugging */

  int k;
  double eps, dot, hold;
  double f1, f2;

  for( k=0, dot=ZERO; k<n; k++ ) dot += SQR( x[k] );
  eps = sqrt(DBL_EPSILON) * (ONE + sqrt(dot));

  hold = x[i];
  x[i] = hold - eps;
  TRY( calc_derivs( obj_func, n, x, fx, grad, func_args ) );
  f1   = fx[0];
  x[i] = hold + eps;
  TRY( calc_derivs( obj_func, n, x, fx, grad, func_args ) );
  f2   = fx[0];
  x[i] = hold;

  TRY( calc_derivs( obj_func, n, x, fx, grad, func_args ) );  /* restore grad */
  return (f2-f1) / (2.*eps);
}

error_msg
debug_derivs( eval_derivs_t * obj_func,
              int n,  /* ndim */
              double * x,
              double * fx,
              double * grad,
              void * func_args ) {

  int i;
  double numdrv, reldiff, absdiff;

  log_printf( "__________\n" );
  log_printf( "1st derivs\n" );
  log_printf( "----------\n" );
  log_printf( "            atom  xyz    anal.drv     num.drv    rel.diff    abs.diff\n" );
  log_printf( "_____________________________________________________________________\n" );

  for( i=0; i<n; i++ ) {  /* 1st derivs */

    absdiff = fabs( grad[i] - (numdrv = calc_num_derivs( obj_func, n, x, fx,
                                                         grad, func_args, i )) );
    if( absdiff != 0 ) {
      reldiff = absdiff / MAX( fabs(grad[i]), fabs(numdrv) );
    }
    else{
      reldiff = 0;
    }
    if( absdiff >= ZERO || reldiff >= ZERO ) {
      log_printf( "            %4d   %c   %10.3f  %10.3f  %10.3e  %10.3f\n",
                  i/3+1,'x'+i%3,grad[i],numdrv,(reldiff)<ONE?reldiff:absdiff,absdiff );
    }
  }

  return SUCCESS;
}


/*****
        L I N E S E A R C H
*****/
error_msg
update_mcsrch_step( double * stx,
                    double * fx,
                    double * dx,
                    double * sty,
                    double * fy,
                    double * dy,
                    double * stp,
                    double fp,
                    double dp,
                    int * brackt,
                    double stpmin,
                    double stpmax ) {
  /* MCSTEP by:
     Jorge J. More', David J. Thuente
     Argonne National Laboratory. MINPACK project. June 1983
     
     C version by IK.
     
     The purpose of MCSTEP is to compute a safeguarded step for
     a linesearch and to update an interval of uncertainty for
     a minimizer of the function.
     
     The parameter stx contains the step with the least function
     value. The parameter stp contains the current step. It is
     assumed that the derivative at stx is negative in the
     direction of the step. If brackt is set true then a
     minimizer has been bracketed in an interval of uncertainty
     with min2points stx and sty.
     
     stx, fx, and dx are variables which specify the step,
     the function, and the derivative at the best step obtained
     so far. The derivative must be negative in the direction
     of the step, that is, dx and stp-stx must have opposite
     signs. On output these parameters are updated appropriately.
     
     sty, fy, and dy are variables which specify the step,
     the function, and the derivative at the other min2point of
     the interval of uncertainty. On output these parameters are
     updated appropriately.
     
     stp, fp, and dp are variables which specify the step,
     the function, and the derivative at the current step.
     If brackt is set true then on input stp must be
     between stx and sty. On output stp is set to the new step.
     
     brackt is a logical variable which specifies if a minimizer
     has been bracketed. If the minimizer has not been bracketed
     then on input brackt must be set false. If the minimizer
     is bracketed then on output brackt is set true.
     
     stpmin and stpmax are input variables which specify lower
     and upper bounds for the step. */

  error_msg err = SUCCESS;
  double gamma, p, q, r, s, sgnd, stpc, stpf, stpq, theta;
  int bound;

  /* Check input arguments for inconsistencies: */

  if( (*brackt && (*stp <= MIN(*stx, *sty) || *stp >= MAX(*stx, *sty))) ||
      (*dx * (*stp - *stx) >= 0) || stpmax < stpmin ) {
    err = ERROR_MSG( "Inconsistent input arguments " );
    return err;
  }

  /* Determine if the derivatives have opposite sign: */
  sgnd = dp * ( *dx / fabs(*dx) );

  if( fp > *fx ) {

    /* First case. A higher function value.
       The minimum is bracketed. If the cubic step is closer
       to stx than the quadratic step, the cubic step is taken,
       else the average of the cubic and quadratic steps is taken. */
    
    bound = TRUE;
    theta = 3 * (*fx - fp) / (*stp - *stx) + *dx + dp;
    s = MAX( fabs(theta), MAX(fabs(*dx), fabs(dp)) );
    gamma = s * sqrt( SQR( theta/s ) - (*dx/s) * (dp/s) );
    if( *stp < *stx ) gamma = -gamma;
    p =  (gamma - *dx) + theta;
    q = ((gamma - *dx) + gamma) + dp;
    r = p / q;
    stpc = *stx + r * (*stp - *stx);
    stpq = *stx +
      ((*dx / ( (*fx - fp) / (*stp - *stx) + *dx )) / 2) * (*stp - *stx);
    if( fabs(stpc - *stx) < fabs(stpq - *stx) ) stpf = stpc;
    else                                        stpf = stpc + (stpq-stpc)/2;
    *brackt = TRUE;

  } else if( sgnd < 0 ) {

    /* Second case. A lower function value and derivatives of
       opposite sign. The minimum is bracketed. If the cubic
       step is closer to stx than the quadratic (secant) step,
       the cubic step is taken, else the quadratic step is taken. */

    bound = FALSE;
    theta = 3 * (*fx - fp) / (*stp - *stx) + *dx + dp;
    s = MAX( fabs(theta), MAX(fabs(*dx), fabs(dp)) );
    gamma = s * sqrt( SQR( theta/s ) - (*dx/s) * (dp/s) );
    if( *stp > *stx ) gamma = -gamma;
    p =  (gamma - dp) + theta;
    q = ((gamma - dp) + gamma) + *dx;
    r = p / q;
    stpc = *stp + r * (*stx - *stp);
    stpq = *stp + (dp / (dp - *dx)) * (*stx - *stp);
    if( fabs(stpc - *stp) > fabs(stpq - *stp) ) stpf = stpc;
    else                                        stpf = stpq;
    *brackt = TRUE;
  } else if( fabs(dp) < fabs(*dx) ) {

    /* Third case. A lower function value, derivatives of the
       same sign, and the magnitude of the derivative decreases.
       The cubic step is only used if the cubic tmin2s to infinity
       in the direction of the step or if the minimum of the cubic
       is beyond stp. Otherwise the cubic step is defined to be
       either stpmin or stpmax. The quadratic (secant) step is also
       computed and if the minimum is bracketed then the the step
       closest to stx is taken, else the step farthest away is taken. */

    bound = TRUE;
    theta = 3 * (*fx - fp) / (*stp - *stx) + *dx + dp;
    s = MAX( fabs(theta), MAX(fabs(*dx), fabs(dp)) );

    /* The case gamma=0 only arises if the cubic does not
       tmin2 to infinity in the direction of the step. */

    gamma = s * sqrt( MAX(ZERO, SQR( theta/s ) - (*dx/s) * (dp/s)) );
    if( *stp > *stx ) gamma = -gamma;
    p = (gamma - dp) + theta;
    q = (gamma + (*dx - dp)) + gamma;
    r = p / q;
    if( r < 0 && gamma != 0 ) stpc = *stp + r * (*stx - *stp);
    else if( *stp > *stx )    stpc = stpmax;
    else                      stpc = stpmin;
    stpq = *stp + (dp / (dp - *dx)) * (*stx - *stp);
    if( *brackt ) {
      if( fabs(*stp - stpc) < fabs(*stp - stpq) )
        stpf = stpc;
      else
        stpf = stpq;
    } else {
      if( fabs(*stp - stpc) > fabs(*stp - stpq) )
        stpf = stpc;
      else
        stpf = stpq;
    }
  } else {

    /* Fourth case. A lower function value, derivatives of the
       same sign, and the magnitude of the derivative does
       not decrease. If the minimum is not bracketed, the step
       is either stpmin or stpmax, else the cubic step is taken. */

    bound = FALSE;
    if( *brackt ) {
      theta = 3 * (fp - *fy) / (*sty - *stp) + *dy + dp;
      s = MAX( fabs(theta), MAX(fabs(*dy), fabs(dp)) );
      gamma = s * sqrt( SQR( theta/s ) - (*dy/s) * (dp/s) );
      if( *stp > *sty )
        gamma = -gamma;
      p =  (gamma - dp) + theta;
      q = ((gamma - dp) + gamma) + *dy;
      r = p / q;
      stpc = *stp + r * (*sty - *stp);
      stpf = stpc;
    } else if( *stp > *stx )
      stpf = stpmax;
    else
      stpf = stpmin;
  }

  /* Update the interval of uncertainty. This update does not
     depend on the new step or the case analysis above. */

  if( fp > *fx ) {
    *sty = *stp;
    *fy = fp;
    *dy = dp;
  } else {
    if( sgnd < 0 ) {
      *sty = *stx;
      *fy = *fx;
      *dy = *dx;
    }
    *stx = *stp;
    *fx = fp;
    *dx = dp;
  }

  /* Compute the new step and safeguard it. */

  stpf = MIN(stpmax, stpf);
  stpf = MAX(stpmin, stpf);
  *stp = stpf;
  if( *brackt && bound ) {
    if( *sty > *stx ) *stp = MIN( *stx + 0.66 * (*sty - *stx), *stp );
    else              *stp = MAX( *stx + 0.66 * (*sty - *stx), *stp );
  }

  return err;
}

error_msg
line_minimize_mcsrch( eval_derivs_t * obj_func,
                      int n,
                      double * s,
                      double * x,
                      double * f,
                      double * g,
                      void * func_args,
                      int verbosity ) {

  /* MCSRCH by:
     Jorge J. More', David J. Thuente
     Argonne National Laboratory. MINPACK project. June 1983
     
     Slightly modified and rewritten in C by IK.
     Still the best line minimizer I have encountered.
     
     The purpose of MCSRCH is to find a step which satisfies
     a sufficient decrease condition and a curvature condition.
     
     At each stage the subroutine updates an interval of
     uncertainty with min2points stx and sty. The interval of
     uncertainty is initially chosen so that it contains a
     minimizer of the modified function
     
       f(x+stp*s) - f(x) - ftol*stp*(gradf(x)'s).
     
     If a step is obtained for which the modified function
     has a nonpositive function value and nonnegative derivative,
     then the interval of uncertainty is chosen so that it
     contains a minimizer of f(x+stp*s).
     
     The algorithm is designed to find a step which satisfies
     the sufficient decrease condition
     
       f(x+stp*s) <= f(x) + ftol*stp*(gradf(x)'s),
     
     and the curvature condition
     
       abs(gradf(x+stp*s)'s)) <= gtol*abs(gradf(x)'s).
     
     If ftol is less than gtol and if, for example, the function
     is bounded below, then there is always a step which satisfies
     both conditions. If no step can be found which satisfies both
     conditions, then the algorithm usually stops when rounding
     errors prevent further progress. In this case stp only
     satisfies the sufficient decrease condition. */

  error_msg err = SUCCESS;
  int i, maxfev, nfev;
  double gtol, ftol, xtol, stp, stpmin, stpmax;
  double *x0=NULL;
  int brackt, stage1;
  double dg, dgm, dginit, dgtest, dgx, dgxm, dgy, dgym, finit, ftest1;
  double fm, fx, fxm, fy, fym;
  double stx, sty, stmin, stmax, width, width1;
  double snorm, smax, xnorm;

  maxfev = MAXFEV;
  gtol   = GTOL;
  ftol   = FTOL;
  xtol   = DBL_EPSILON;
  stpmax = 99;
  stpmin =  0;
  stp    = ONE;  /* initial step is set to  full step */
  MALLOC( x0, n );

  /* Check search direction: */

  for( i=0, dginit=ZERO; i<n; i++ ) dginit += g[i] * s[i];

  if( dginit ==ZERO ) { /* nothing to do, aleady at line minimum */
    goto cleanup;
  }

  /* check if initial step too agressive */

  for( i=0, smax=snorm=xnorm=ZERO; i<n; i++ ) {
    snorm += SQR( s[i] );
    //xnorm += SQR( x[i] );
    if( fabs( s[i] ) > smax )
        smax = fabs( s[i] );
  }
  snorm = sqrt( snorm );
  //xnorm = sqrt( xnorm );

  if( smax > MAXCRDMOV ) stp = MAXCRDMOV / smax;             /* limit max coord movement */

  /* Initialize local variables: */

  brackt = FALSE;
  stage1 = TRUE;
  nfev = 0;
  finit = *f;
  dgtest = ftol * dginit;
  width = stpmax - stpmin;
  width1 = 2 * width;
  for( i=0; i<n; i++ ) x0[i] = x[i];

  /* The variables stx, fx, dgx contain the values of the step,
     function, and directional derivative at the best step.
     The variables sty, fy, dgy contain the value of the step,
     function, and derivative at the other min2point of
     the interval of uncertainty.
     The variables stp, f, dg contain the values of the step,
     function, and derivative at the current step. */

  stx = ZERO;
  fx = finit;
  dgx = dginit;
  sty = ZERO;
  fy = finit;
  dgy = dginit;


  /* START ITERATION LOOP: */

  for( ;; ) {

    /* Set the minimum and maximum steps to correspond
       to the present interval of uncertainty: */

    if( brackt ) {
      stmin = MIN( stx, sty );
      stmax = MAX( stx, sty );
    } else {
      stmin = stx;
      stmax = stp + 4 * (stp - stx);
    }

    /* Force the step to be within the bounds stpmax and stpmin: */

    stp = MAX( stp, stpmin );
    stp = MIN( stp, stpmax );

    /* If an unusual termination is to occur then let stp be the lowest point
       obtained so far: */

    if( nfev >= maxfev-1 ||
        (brackt && (stp <= stmin || stp >= stmax)) ||
        (brackt && (stmax-stmin <= xtol*stmax)) )
      stp = stx;

    /* Evaluate the function and gradient at stp and compute the directional
       derivative: */

    for( i=0; i<n; i++ ) x[i] = x0[i] + stp * s[i];

    if( (err = calc_derivs( obj_func, n, x, f, g, func_args )) != SUCCESS ) {
      ERROR(( err ));
      goto cleanup;
    }

    nfev++;
    for( i=0, dg=ZERO; i<n; i++ ) dg += g[i] * s[i];
    ftest1 = finit + stp * dgtest;

    /* TEST FOR CONVERGENCE: */

    if( *f <= ftest1 && fabs(dg) <= gtol * (-dginit) ) {
                                                            /* line minimization succesful */
      break;
    }

    /* TEST FOR ERRORS: */

    if( (brackt && (stp <= stmin || stp >= stmax)) ) {
      err = ERROR_MSG( "Line minimizer aborted: rounding error" );
    }

    if( stp == stpmax && *f <= ftest1 && dg <= dgtest ) {
      err = ERROR_MSG( "Line minimizer aborted: step at upper bound" );
    }

    if( stp == stpmin && (*f > ftest1 || dg >= dgtest) ) {
      err = ERROR_MSG( "Line minimizer aborted: step at lower bound" );
    }

    if( nfev >= maxfev ) {
      err = ERROR_MSG( "Line minimizer aborted: too many steps" );
    }

    if ( brackt && (stmax - stmin <= xtol * stmax) ) {
      err = ERROR_MSG( "Line minimizer aborted: interval of uncertainty too small" );
    }

    if( err != SUCCESS ) goto cleanup;  /* err can be rewritten multiple times,
                                           only the final message is relevant */

    /* MORE WORK TO DO.
       In the first stage we seek a step for which the modified
       function has a nonpositive value and a nonnegative derivative: */

    if( stage1 && *f <= ftest1 && dg >= MIN( ftol, gtol ) * dginit )
      stage1 = FALSE;

    /* A modified function is used to predict the step only if
       we have not obtained a step for which the modified
       function has a nonpositive function value and nonnegative
       derivative, and if a lower function value has been
       obtained but the decrease is not sufficient: */

    if( stage1 && *f <= fx && *f > ftest1 ) {

      /* Define the modified function and derivative values: */

      fm = *f - stp * dgtest;
      fxm = fx - stx * dgtest;
      fym = fy - sty * dgtest;
      dgm = dg - dgtest;
      dgxm = dgx - dgtest;
      dgym = dgy - dgtest;

      /* Update the interval of uncertainty and compute the new step: */

      if( (err = update_mcsrch_step( &stx, &fxm, &dgxm, &sty, &fym, &dgym, &stp,
                                     fm, dgm, &brackt, stmin, stmax )) != SUCCESS ) {
        ERROR(( err ));
        goto cleanup;
      }

      /* Reset the function and gradient values for f: */

      fx = fxm + stx * dgtest;
      fy = fym + sty * dgtest;
      dgx = dgxm + dgtest;
      dgy = dgym + dgtest;

    } else if( (err = update_mcsrch_step( &stx, &fx, &dgx, &sty, &fy, &dgy, &stp,
                                          *f, dg, &brackt, stmin, stmax )) != SUCCESS ) {
             ERROR(( err ));
             goto cleanup;
           }

    /* Force a sufficient decrease in the size of the
       interval of uncertainty: */

    if( brackt ) {
      if( fabs(sty-stx) >= .66*width1 )
        stp = stx + .5 * (sty - stx);
      width1 = width;
      width = fabs(sty - stx);
    }
  }

  /* END ITERATION LOOP. */

 cleanup:
  
  if( verbosity == 2 ) {
    log_printf( "   LS: rel_step= %15.8g  abs_step=%15.8g  it=%2d\n", stp, stp*snorm, nfev );
  }

  FREE( x0 );

  return err;
}


/*****
        M I N I M I Z E R
*****/
error_msg
init_lbfgs_matrix( lbfgs_t * lbfgsmat,
                   int m,
                   int n ) {

  error_msg err = SUCCESS;
  int i,j;
  
  for( i=0; i<m; i++ ) {
    MALLOC( lbfgsmat[i].s, n );
    MALLOC( lbfgsmat[i].y, n );
  }

  return err;
}

error_msg
load_lbfgs_matrix( lbfgs_t * lbfgsmat,
                   int current,
                   int n,
                   double * s,
                   double * y ) {

  error_msg err = SUCCESS;
  int i;
  double rho, gamma;

  for( i=0, rho=gamma=ZERO; i<n; i++ ) {
    lbfgsmat[current].s[i] = s[i];
    lbfgsmat[current].y[i] = y[i];
    rho   += s[i] * y[i];
    gamma += y[i] * y[i];
  }

  if( rho==ZERO || gamma==ZERO ) {
    err = ERROR_MSG( "L-BFGS matrix singular: YS=0 or YY=0" );
    return err;
  }

  lbfgsmat[current].rho   = ONE / rho;
  lbfgsmat[current].gamma = fabs(rho / gamma);

  return err;
}

error_msg
nocedal( lbfgs_t * lbfgsmat,
         int m,
         int current,
         int n,
         double sign,
         int update,
         double * hdiag,
         double * vec_in,
         double * vec_out ) {

  error_msg err = SUCCESS;
  int i, j;
  double *alpha=NULL, beta, dot, *tmp_in=NULL;

  MALLOC( alpha, m );
  MALLOC( tmp_in, n );

  for( i=0; i<n; i++ ) tmp_in[i] = sign * vec_in[i];
  for( j=m-1; j>=0; j-- ) {
    for( i=0, dot=0; i<n; i++ )
      dot += lbfgsmat[j].s[i] * tmp_in[i];
    alpha[j] = lbfgsmat[j].rho * dot;
    for( i=0; i<n; i++ )
      tmp_in[i] -= alpha[j] * lbfgsmat[j].y[i];
  }
  switch( update ) {

    case DIAGONAL:
      for( i=0; i<n; i++ )
        vec_out[i] = tmp_in[i] / fabs( hdiag[i] );  /* diagonal - hdiag should be an estimate */
      break;                                        /* at the minimum, NOT at current iterate */

    case SCALING:
      for( i=0; i<n; i++ )
        vec_out[i] = tmp_in[i] * lbfgsmat[current].gamma;  /* scaling */
      break;
    default:
      memcpy( vec_out, tmp_in, n*sizeof(vec_out[0]) );  /* unit matrix */
      break;
  }
  for( j=0; j<m; j++ ) {
    for( i=0, dot=0; i<n; i++ )
      dot += lbfgsmat[j].y[i] * vec_out[i];
    beta = lbfgsmat[j].rho * dot;
    for( i=0; i<n; i++ )
      vec_out[i] += lbfgsmat[j].s[i] * (alpha[j] - beta);
  }
  
 cleanup:
  
  FREE( alpha );
  FREE( tmp_in );

  return err;
}

error_msg
lbfgs_minimize( eval_derivs_t * obj_func,
                void * func_args,
                int ndim,
                double * x,
                double * fx,
                double * grad,
                min_args_t * min_args ) {

  func_args_t * _func_args = (func_args_t *) func_args;
  error_msg err = SUCCESS;
  int i, min_iter, min_conv, do_transform=0, in_place=1;
  int n_frame, save_path=min_args->save_path, which_path=min_args->which_path;
  double grms, dgrad, rmsd[1], rmsd1[1], rmsd2[1];
  double *p=NULL, *g_old=NULL, *g_dif=NULL;
  double *step=NULL, *x_old=NULL, *x_save=NULL;
  double ** path=min_args->path;
  double rmsd_gap=min_args->rmsd_gap, rmsd_save=0, rmsd_frame, *frame_poi;
  lbfgs_t *lbfgs_matrix=NULL;
  FILE * debug_file=NULL;
  char en_term[MAX_EN][20]={ "total energy", "bond stretch", "angle restraint", "angle bending",
                             "torsion restraint", "torsion", "total SANM", "SANM 1", "SANM 2",
                             "POSRE", "nonbonded bump" };

  /* Check input arguments: */
  if( x==NULL    ||
      fx==NULL   ||
      grad==NULL ||
      min_args==NULL ) {
    err = ERROR_MSG( "Bad input arguments" );
    return err;
  }
  if( ndim<2                ||
      min_args->m_lbfgs<=0  ||
      min_args->iter_in<0   ||  /* =0 -> function evaluation only */
      min_args->grms_in<0   ||
      min_args->verbosity<0 ||
      min_args->verbosity>2 ) {
    err = ERROR_MSG( "Bad input arguments" );
    return err;
  }

  if( min_args->visual_debug ) {

    if( (debug_file = fopen("min-debug.pdb", "a")) == NULL ) {

      err = ERROR_MSG( "Cannot open min-debug.pdb" );
      return err;
    }
  }

  MALLOC( p,      ndim );
  MALLOC( g_old,  ndim );
  MALLOC( g_dif,  ndim );
  MALLOC( step,   ndim );
  MALLOC( x_old,  ndim );
  MALLOC( x_save, ndim );
  CLEAR( g_old, ndim );
  CLEAR( x_old, ndim );

  MALLOC( lbfgs_matrix, min_args->m_lbfgs );
  if( (err = init_lbfgs_matrix( lbfgs_matrix, min_args->m_lbfgs, ndim )) != SUCCESS ) {

    ERROR(( err ));
    goto cleanup;
  }

  /* START MAIN MINIMIZATION LOOP: */

  if( min_args->verbosity > 0 ) {
    log_printf( "________________________________________________________________\n" );
  }

  min_conv = min_iter = 0;

  do{

    if( min_args->visual_debug ) TRY( write_pdb_file( ndim/3, func_args, x, debug_file, 1 ) );

    /* Update derivs: */
    if( min_iter == 0 ) {

      if( (err = calc_derivs( obj_func, ndim, x, fx, grad, func_args )) != SUCCESS ) {

        ERROR(( err ));
        goto cleanup;
      }
      if( save_path ) {

        /* Initialize path: */

        n_frame = min_args->n_path_frame[min_args->which_path];

        if( n_frame==0 ) {
          log_printf( "\n--------------\n" );
          log_printf( " Path Frames\n" );
          log_printf( "--------------\n" );
          log_printf( "        Energy    RMSD  RMSD-A  RMSD-B\n" );
          log_printf( "______________________________________\n" );
          MALLOC( path[n_frame], ndim );
          COPY( path[n_frame], x, ndim );   /* first frame */
          TRY( rmsfit( ndim/3, _func_args->min1_xyz, x, rmsd1, in_place ) );
          TRY( rmsfit( ndim/3, _func_args->min2_xyz, x, rmsd2, in_place ) );
          log_printf( "%4d  %8.2f          %6.2f  %6.2f\n", n_frame+1, fx[6], rmsd1[0], rmsd2[0] );  /* print only sanm energy */
          ++n_frame;
        } else { log_printf( "--------------------------------------\n" ); }
      }
    }

    if( min_iter > 0 ) {
      /* store last m_lbfgs step[] and g_dif[] vectors in circular order */
      for( i=0; i<ndim; i++ ) {
        step[i] = x[i] - x_old[i];
        g_dif[i] = grad[i] - g_old[i];
      }
      if( (err = load_lbfgs_matrix( lbfgs_matrix, (min_iter-1)%min_args->m_lbfgs, ndim, step, g_dif )) != SUCCESS ) {
        ERROR(( err ));
        WARNING(( "L-BFGS minimizer aborted: saving current coordinates" ));
        break;
      }
    }

    /* Check for exit */

    for( i=0, grms=ZERO; i<ndim; i++ ) grms  += SQR( grad[i] );
    grms = sqrt( grms/ndim );
    if( min_args->verbosity > 0 ) {
      log_printf( " MIN:                          It= %4d  E= %10.2f (%7.3f)\n", min_iter, fx[0], grms );
    }
    if( grms <= min_args->grms_in || min_iter == min_args->iter_in ) {
      min_conv = YES;
      break;
    }

    /* Update L-BFGS search direction: */

    if( min_iter > 0 ) {

      /* Update L-BFGS search direction: */

      if( (err = nocedal( lbfgs_matrix, MIN(min_iter, min_args->m_lbfgs), (min_iter-1)%min_args->m_lbfgs,
                          ndim, -ONE, SCALING, NULL, grad, p )) != SUCCESS ) {
        ERROR(( err ));
        goto cleanup;
      }
      /*  make sure that search direction is downhill! */
      for (i=0, dgrad=ZERO; i<ndim; i++)
        dgrad += grad[i] * p[i];
      if (dgrad > ZERO) {  /* Uphill movement! Replace LBFGS step with -grad */
        for (i=0; i<ndim; i++)
          p[i] = -grad[i];
      }
    } else for( i=0; i<ndim; i++ ) p[i] = -grad[i];  /* first step always steepest descent */

    memcpy( x_old, x,    ndim*sizeof(x[0]) );     /* save old coord before LS move */
    memcpy( g_old, grad, ndim*sizeof(grad[0]) );  /* save old grad  before LS move */

    /* At this point L-BFGS search direction is given in p[].
       Find the line minimum: */

    if( (err = line_minimize_mcsrch( obj_func, ndim, p, x, fx, grad, 
                                     func_args, min_args->verbosity )) != SUCCESS ) {

      if( min_args->visual_debug ) TRY( write_pdb_file( ndim/3, func_args, x, debug_file, 1 ) );
      if( min_args->verbosity == 2 ) {
        ERROR(( err ));
      }
      /*goto cleanup;*/  /* ignore error, continue with minimization */
    }

    if( save_path ) {

      /* Add frame to path if rmsd to last frame >= rmsd gap: */
    
      TRY( rmsfit( ndim/3, path[n_frame -1], x, rmsd, do_transform ) );
      if( rmsd[0] >= rmsd_gap ) {
    
        MALLOC( path[n_frame], ndim );

        /* Check if structure in previous min_iter was actually closer to the rmsd gap value: */

        if( (rmsd[0] - rmsd_gap) < (rmsd_gap - rmsd_save) || rmsd[0] >= 2*rmsd_gap ) {

          frame_poi = x;
          rmsd_frame = rmsd[0];
          rmsd[0] = 0.0;  /* we are right at the new frame */
        } else {

          frame_poi = x_save;
          rmsd_frame = rmsd_save;
          TRY( rmsfit( ndim/3, x_save, x, rmsd, do_transform ) );  /* re-fit the latest x coords to the new frame x_save */
        }

        COPY( path[n_frame], frame_poi, ndim );   /* add frame */
        ++n_frame;
        TRY( rmsfit( ndim/3, _func_args->min1_xyz, frame_poi, rmsd1, in_place ) );
        TRY( rmsfit( ndim/3, _func_args->min2_xyz, frame_poi, rmsd2, in_place ) );
        log_printf( "%4d  %8.2f  %6.2f  %6.2f  %6.2f\n", n_frame, fx[6], rmsd_frame, rmsd1[0], rmsd2[0] );  /* print only sanm energy */
#       if 0
        for( i=0; i<MAX_EN; i++ ) log_printf( "%20s  %20.10g\n", en_term[i], fx[i] );
#       endif
      }
      COPY( x_save, x, ndim );
      rmsd_save = rmsd[0];
    }
  } while( ++min_iter <= min_args->iter_in );

  /* END MAIN MINIMIZATION LOOP. */

  /* Sanity check: */
  if( (err = calc_derivs( obj_func, ndim, x, fx, grad, func_args )) != SUCCESS ) {
    ERROR(( err ));
    goto cleanup;
  }

  for( i=0, grms=ZERO; i<ndim; i++ ) grms += SQR( grad[i] );
  grms = sqrt( grms/ndim );
  if( min_args->verbosity > 0 ) {
    log_printf( "----------------------------------------------------------------\n" );
    log_printf( " FIN:             %s                    E= %10.2f (%7.3f)\n\n",
                min_conv==YES?":-)":":-(", fx[0], grms );
  }

# if 0
  for( i=0; i<MAX_EN; i++ ) log_printf( "%20s  %20.10g\n", en_term[i], fx[i] );
# endif

 cleanup:

  FREE( p );
  FREE( g_old );
  FREE( g_dif );
  FREE( step );
  FREE( x_old );
  FREE( x_save );
  for( i=0; i<min_args->m_lbfgs; i++ ) {
    FREE( lbfgs_matrix[i].s );
    FREE( lbfgs_matrix[i].y );
  }
  FREE( lbfgs_matrix );
  if( debug_file !=NULL ) fclose( debug_file );

  min_args->iter_out = min_iter;
  min_args->grms_out = grms;
  min_args->n_path_frame[min_args->which_path] = n_frame;

  return err;
}


error_msg
minimize_reference_structures( int n_atom,
                               void * __restrict args,
                               FILE * pdb_file,
                               double ** ref_confs,
                               eval_derivs_t * obj_func,
                               int n_dim,
                               double * x,
                               double * fx,
                               double * grad,
                               min_args_t * min_args ) {

  func_args_t * __restrict func_args = (func_args_t *) args;
  int i=0, n=0, do_transform=0, in_place=1;
  double rmsd[1], rmsd1[1], rmsd2[1];
  char line[MAXBUF];
  error_msg err = SUCCESS;

  log_printf( " -----------------------------------\n" );
  log_printf( " SANM-Minimized Reference Structures\n" );
  log_printf( " -----------------------------------\n" );
  log_printf( "        Energy    RMSD  RMSD-A  RMSD-B\n" );
  log_printf( "______________________________________\n" );

  rewind( pdb_file );  /* just in case */
  while( fgets( line, MAXBUF, pdb_file ) != NULL ) {

    if( strstr( line, " CA") == NULL )
      continue;
    
    sscanf( line+30, "%lf %lf %lf", &x[3*i], &x[3*i+1], &x[3*i+2] );

    if( ++i < n_atom )
      continue; /* just read the coordinates */

    i = 0;      /* whole conformation read, reset counter */

    TRY( lbfgs_minimize( obj_func, func_args, n_dim, x, fx, grad, min_args ) );
    if( n>0 )
      TRY( rmsfit( n_atom, ref_confs[0], x, rmsd, do_transform ) );

    MALLOC( ref_confs[n], n_dim );
    COPY( ref_confs[n], x, n_dim );

    TRY( rmsfit( n_atom, func_args->min1_xyz, x, rmsd1, in_place ) );
    TRY( rmsfit( n_atom, func_args->min2_xyz, x, rmsd2, in_place ) );
    if( n==0 )
      log_printf( "   %c  %8.2f          %6.2f  %6.2f\n", 'A'+n, fx[6], rmsd1[0], rmsd2[0] );  /* print only sanm energy */
    else
      log_printf( "   %c  %8.2f  %6.2f  %6.2f  %6.2f\n", 'A'+n, fx[6], rmsd[0], rmsd1[0], rmsd2[0] );
    ++n;
  }

  if( ferror( pdb_file) ) {

    err = ERROR_MSG( "File read error" );
    goto cleanup;
  }

cleanup:

  return err;
}



/*****
        M A I N
*****/
int
main( int argc,
      char ** argv ) {

  int i, j, i_save, j_save, n_atom, n_dim, do_transform=0, in_place=1;
  int normal_style=1, plumed_style=2, n_path_frame[2], n_tail_frame[2];
  double rmsd[1], rmsd_gap, rmsd_avg, rmsd_min, rmsd_max;
  double fc_posre=0.1;
  double *xyz=NULL, *grad=NULL, energy[MAX_EN], *tmp=NULL;
  double *ref_confs_min[2], *path_frames[2*MAXFRAMES];
  double elapsed;
  FILE *input_file=NULL, *output_file=NULL, *plumed_file=NULL, *output_file2=NULL, *plumed_file2=NULL;
  char input_fname[MAXBUF], output_fname[MAXBUF], log_fname[MAXBUF], plumed_fname[MAXBUF];
  char output_fname2[MAXBUF], plumed_fname2[MAXBUF];
  min_args_t min_args[1];
  func_args_t func_args[1]; func_args->r0_min1=NULL; func_args->r0_min2=NULL;
  eval_derivs_t sanm_energy[NPOT+1] = { calc_bond_derivs,
                                        calc_angle_derivs,
                                        calc_bump_derivs,
                                        calc_sanm_derivs,     /* min1 anm */
                                        calc_sanm_derivs,     /* min2 anm */
                                        NULL
                                      },
                pull_energy[NPOT+1] = { calc_bond_derivs,
                                        calc_angle_derivs,
                                        calc_bump_derivs,
                                        calc_sanm_derivs,     /* min1 anm */
                                        calc_sanm_derivs,     /* min2 anm */
                                        calc_posre_derivs,
                                        NULL
                                      };

  if( argc != 5 && argc != 6 ) {
    ERROR(( "Usage:  %s  n_CA_atoms  endpoints_confs_pdb_fname  CA_transition_path_fname(no extension)  rmsd_gap_between_output_frames_Angs  [fc_pull]", argv[0] ));
    goto cleanup;
  }

  /* Read command line */

  n_atom   = atoi( argv[1] );       /* number of CA atoms in the ANM model */
  strcpy( input_fname,  argv[2] );
  strcpy( output_fname, argv[3] );
  strcpy( plumed_fname,  output_fname );
  strcpy( output_fname2, output_fname );
  strcpy( plumed_fname2, output_fname );
  rmsd_gap = atof( argv[4] );
  if( argc == 6 ) {

    fc_posre = atof( argv[5] );
  }

  for( i=0; i<argc; i++ ) {

    log_printf( "%s  ", argv[i] );
  }
  log_printf( "\n\n" );

  /* Allocate local arrays */

  n_dim = 3*n_atom;
  MALLOC( xyz,       n_dim );
  MALLOC( grad,      n_dim );
  MALLOC( func_args->min1_xyz, n_dim );
  MALLOC( func_args->min2_xyz, n_dim );
  MALLOC( func_args->posre_xyz, n_dim );
  MALLOC( func_args->r0_min1, n_atom*n_atom );  /* could use a linked list including only r0 <= ANM_CUTOFF */
  MALLOC( func_args->r0_min2, n_atom*n_atom );
  MALLOC( func_args->pdb_file_image, n_atom*MAXBUF );
  MALLOC( tmp,       n_dim );

  /* Initialize func args */

  func_args->fc_posre           =   fc_posre;
  func_args->fc_str             =   1e5;
  func_args->r0_str             =   3.8;            /* typical CA-CA dist ~ 3.8 Angs in the PDB */
  func_args->fbhw_str           =   0.2;
  func_args->fc_ang             =   0.0; //  1.0
  func_args->t0_ang             = 110.0*(PI/180);
  func_args->fbhw_ang           =  25.0*(PI/180);   /* 75 <= CA-CA-CA ang <= 145 Deg in the PDB */
  func_args->k0_ang             =   1e5;
  func_args->fc_dih             =   0.0; // -0.5
  func_args->p0_dih             =  -PI;
  func_args->fbhw_dih           = 175*(PI/180);     /* keep CA-CA-CA-CA dihedrals from colinear */
  func_args->k0_dih             =   0.0;
  func_args->fc_bump            =   1e5;
  func_args->dmin_bump          =   3.8;
  func_args->fc_sanm            =   1.0;

  input_file  = fopen( input_fname, "r" );
  if( input_file ==NULL ) goto cleanup;
 
  min_args->iter_in      = ITER_MAX;
  min_args->grms_in      = GRMS_MAX;
  min_args->m_lbfgs      = M_LBFGS;
  min_args->verbosity    = VERBOSITY;
  min_args->visual_debug = VISUAL_DEBUG;
  min_args->save_path    = SAVE_PATH;
  min_args->rmsd_gap     = rmsd_gap;

#if 0
  goto test;
#endif

  elapsed = wallclock();

  TRY( read_reference_structures( n_atom, func_args, input_file ) );  /* read both reference structures, copy them to func_args->min[1,2]_xyz */
  TRY( minimize_reference_structures( n_atom, func_args, input_file, 
                                      ref_confs_min, sanm_energy,  
                                      n_dim, xyz, energy, grad, min_args ) );

  /* Build a path from min1 to min2 and also a reverse path from min2 to min1: */

  n_path_frame[0] = 
  n_path_frame[1] = 
  min_args->n_path_frame[0] = 
  min_args->n_path_frame[1] = 0;
  COPY( xyz, func_args->min1_xyz, n_dim );
  log_printf( "\n\n F O R W A R D" );
  for( i=0; i<2; i++ ) {  /* forward direction and reverse */

    min_args->save_path = YES;
    min_args->path = &path_frames[i*MAXFRAMES];
    min_args->which_path = i;
  //min_args->visual_debug=1;
    TRY( lbfgs_minimize( sanm_energy, func_args, n_dim, xyz, energy, grad, min_args ) );
    n_tail_frame[i] = min_args->n_path_frame[i];

    COPY(                  xyz, path_frames[i*MAXFRAMES+n_tail_frame[i] -1], n_dim );
    COPY( func_args->posre_xyz, ref_confs_min[1-i], n_dim );
    TRY( lbfgs_minimize( pull_energy, func_args, n_dim, xyz, energy, grad, min_args ) );
    n_path_frame[i] = min_args->n_path_frame[i];

    COPY(                  xyz, func_args->min2_xyz, n_dim );
    if( i==0 ) log_printf( "\n\n R E V E R S E" );
  }

  /* Add opposite tails: */

  for( i=0; i<n_tail_frame[1]; i++ ) {  /* add back tail to forward path */

    MALLOC( path_frames[n_path_frame[0] +i], n_dim );
    COPY(   path_frames[n_path_frame[0] +i], path_frames[MAXFRAMES+n_tail_frame[1] -1 -i], n_dim );
    /* no rmsd fit necessary, it is done below */
  }

  for( i=0; i<n_tail_frame[0]; i++ ) {  /* add front tail to reverse path */

    MALLOC( path_frames[MAXFRAMES+n_path_frame[1] +i], n_dim );
    COPY(   path_frames[MAXFRAMES+n_path_frame[1] +i], path_frames[n_tail_frame[0] -1 -i], n_dim );
    /* no rmsd fit necessary, it is done below */
  }

  n_path_frame[0] += n_tail_frame[1];
  n_path_frame[1] += n_tail_frame[0];

  /* Save both forward and reverse paths in two output files, one multi-pdb and the other PLUMED style: */

  /* Forward: */
  output_file = fopen( strcat(output_fname, "-forward-models.pdb"), "w" );
  if( output_file==NULL ) goto cleanup;
  plumed_file = fopen( strcat(plumed_fname, "-forward-plumed.pdb"), "w" );
  if( plumed_file==NULL ) goto cleanup;

  rmsd_avg = ZERO;
  rmsd_min = 1e10;
  rmsd_max = ZERO;
  for( i=0; i <n_path_frame[0]; i++ ) {  /* forward path */

    if( i>0 ) { 
      TRY( rmsfit( n_atom, path_frames[i-1], path_frames[i], rmsd, do_transform ) );
      rmsd_avg += rmsd[0];
      if( rmsd[0] < rmsd_min ) rmsd_min = rmsd[0];
      if( rmsd[0] > rmsd_max ) rmsd_max = rmsd[0];
      //log_printf( "%6.4f, ", rmsd[0] );
    }
    write_pdb_file( n_atom, func_args, path_frames[i], output_file, normal_style );
    write_pdb_file( n_atom, func_args, path_frames[i], plumed_file, plumed_style );
  }
  rmsd_avg /= (n_path_frame[0] -1);

  /* Write PLUMED commands to run meta-eABF on this path CV: */

  log_printf( "\n\n------------------------------------\n" );
  log_printf( " PLUMED meta-eABF commands - Forward\n" );
  log_printf( "------------------------------------\n" );
  log_printf( "    plumed_force = \"\"\"\n" );
  log_printf( "#RESTART\n" );
  log_printf( "UNITS LENGTH=nm TIME=ps ENERGY=kj/mol\n" );
  log_printf( "TIME LABEL=tim\n" );
  log_printf( "PATHMSD REFERENCE=%s LAMBDA=%d NEIGH_STRIDE=4 NEIGH_SIZE=8 LABEL=path\n", plumed_fname, (int)(230./SQR(rmsd_avg)) );
  log_printf( "CUSTOM ARG=path.zzz FUNC=sqrt(x) PERIODIC=NO LABEL=path_z\n" );
  log_printf( "DRR ARG=path.sss FULLSAMPLES=2000 GRID_MIN=1 GRID_MAX=%d GRID_SPACING=0.05 TEMP=310 FRICTION=8.0 TAU=0.5 TEXTOUTPUT OUTPUTFREQ=500000 HISTORYFREQ=5000000 DRR_RFILE=drr LABEL=drr\n", n_path_frame[0] );
  log_printf( "METAD ARG=drr.path.sss_fict SIGMA=0.25 HEIGHT=1.8 PACE=500 GRID_MIN=0 GRID_MAX=%d GRID_SPACING=0.05 BIASFACTOR=20 TEMP=310 FILE=HILLS LABEL=metad\n", n_path_frame[0] +1 );
  log_printf( "CUSTOM ARG=tim FUNC=0.5*(erf((2000-x)*0.1)+1) PERIODIC=NO LABEL=stp1\n" );
  log_printf( "CUSTOM ARG=tim FUNC=0.5*(erf((x-2000)*0.1)+1) PERIODIC=NO LABEL=stp2\n" );
  log_printf( "CUSTOM ARG=tim FUNC=0.5*(erf((3000-x)*0.1)+1) PERIODIC=NO LABEL=stp3\n" );
  log_printf( "CUSTOM ARG=tim FUNC=0.5*(erf((x-3000)*0.1)+1) PERIODIC=NO LABEL=stp4\n" );
  log_printf( "CUSTOM ARG=tim FUNC=0.5*(erf((4000-x)*0.1)+1) PERIODIC=NO LABEL=stp5\n" );
  log_printf( "CUSTOM ARG=tim FUNC=0.5*(erf((x-4000)*0.1)+1) PERIODIC=NO LABEL=stp6\n" );
  log_printf( "CUSTOM ARG=tim FUNC=0.5*(erf((5000-x)*0.1)+1) PERIODIC=NO LABEL=stp7\n" );
  log_printf( "CUSTOM ARG=tim FUNC=0.5*(erf((x-5000)*0.1)+1) PERIODIC=NO LABEL=stp8\n" );
  log_printf( "CUSTOM ARG=stp1,stp2,stp3,stp4,stp5,stp6,stp7,stp8 VAR=x1,x2,x3,x4,x5,x6,x7,x8 FUNC=0.15*x1+0.14*x2*x3+0.13*x4*x5+0.12*x6*x7+0.11*x8 PERIODIC=NO LABEL=max_z\n" );
  log_printf( "CUSTOM ARG=path_z,max_z FUNC=(y-x)^-1 PERIODIC=NO LABEL=uwz_arg\n" );
  log_printf( "BIASVALUE ARG=uwz_arg LABEL=uw_path_z\n" );
  log_printf( "UPPER_WALLS ARG=path.sss AT=%d KAPPA=1000 LABEL=uw_s\n", n_path_frame[0] );
  log_printf( "LOWER_WALLS ARG=path.sss AT=1  KAPPA=1000 LABEL=lw_s\n" );
  log_printf( "UPPER_WALLS ARG=drr.path.sss_fict AT=%d KAPPA=1000 LABEL=uw_s_fict\n", n_path_frame[0] );
  log_printf( "LOWER_WALLS ARG=drr.path.sss_fict AT=1  KAPPA=1000 LABEL=lw_s_fict\n" );


  log_printf( "PRINT ARG=* STRIDE=1000 FILE=forward-plumed-force.log FMT=%c12.4f\n", 37 );
  log_printf( "FLUSH STRIDE=1000\n" );
  log_printf( "\"\"\"\n" );

  /* Check for kink in the path: */
  if( fabs( rmsd_avg - rmsd_min ) / rmsd_avg > 0.3 
   || fabs( rmsd_avg - rmsd_max ) / rmsd_avg > 0.3 ) {

    log_printf( "\n WARNING: kink found in the path, try slightly changing rmsd_gap\n\n" );
    log_printf( " rmsd_avg=%6.2f  rmsd_min=%6.2f  rmsd_max=%6.2f\n", rmsd_avg, rmsd_min, rmsd_max );
  }


  /* Reverse: */
  output_file2 = fopen( strcat(output_fname2, "-reverse-models.pdb"), "w" );
  if( output_file2==NULL ) goto cleanup;
  plumed_file2 = fopen( strcat(plumed_fname2, "-reverse-plumed.pdb"), "w" );
  if( plumed_file2==NULL ) goto cleanup;

  rmsd_avg = ZERO;
  rmsd_min = 1e10;
  rmsd_max = ZERO;
  for( i=0; i <n_path_frame[1]; i++ ) {  /* reverse path */

    if( i>0 ) {
      TRY( rmsfit( n_atom, path_frames[MAXFRAMES+i-1], path_frames[MAXFRAMES+i], rmsd, do_transform ) );
      rmsd_avg += rmsd[0];
      if( rmsd[0] < rmsd_min ) rmsd_min = rmsd[0];
      if( rmsd[0] > rmsd_max ) rmsd_max = rmsd[0];
      //log_printf( "%6.4f, ", rmsd[0] );
    }
    write_pdb_file( n_atom, func_args, path_frames[MAXFRAMES+i], output_file2, normal_style );
    write_pdb_file( n_atom, func_args, path_frames[MAXFRAMES+i], plumed_file2, plumed_style );
  }
  rmsd_avg /= (n_path_frame[1] -1);

  /* Write PLUMED commands to run meta-eABF on this path CV: */

  log_printf( "\n\n------------------------------------\n" );
  log_printf( " PLUMED meta-eABF commands - Reverse\n" );
  log_printf( "------------------------------------\n" );
  log_printf( "    plumed_force = \"\"\"\n" );
  log_printf( "#RESTART\n" );
  log_printf( "UNITS LENGTH=nm TIME=ps ENERGY=kj/mol\n" );
  log_printf( "TIME LABEL=tim\n" );
  log_printf( "PATHMSD REFERENCE=%s LAMBDA=%d NEIGH_STRIDE=4 NEIGH_SIZE=8 LABEL=path\n", plumed_fname2, (int)(230./SQR(rmsd_avg)) );
  log_printf( "CUSTOM ARG=path.zzz FUNC=sqrt(x) PERIODIC=NO LABEL=path_z\n" );
  log_printf( "DRR ARG=path.sss FULLSAMPLES=2000 GRID_MIN=1 GRID_MAX=%d GRID_SPACING=0.05 TEMP=310 FRICTION=8.0 TAU=0.5 TEXTOUTPUT OUTPUTFREQ=500000 HISTORYFREQ=5000000 DRR_RFILE=drr LABEL=drr\n", n_path_frame[1] );
  log_printf( "METAD ARG=drr.path.sss_fict SIGMA=0.25 HEIGHT=1.8 PACE=500 GRID_MIN=0 GRID_MAX=%d GRID_SPACING=0.05 BIASFACTOR=20 TEMP=310 FILE=HILLS LABEL=metad\n", n_path_frame[1] +1 );
  log_printf( "CUSTOM ARG=tim FUNC=0.5*(erf((2000-x)*0.1)+1) PERIODIC=NO LABEL=stp1\n" );
  log_printf( "CUSTOM ARG=tim FUNC=0.5*(erf((x-2000)*0.1)+1) PERIODIC=NO LABEL=stp2\n" );
  log_printf( "CUSTOM ARG=tim FUNC=0.5*(erf((3000-x)*0.1)+1) PERIODIC=NO LABEL=stp3\n" );
  log_printf( "CUSTOM ARG=tim FUNC=0.5*(erf((x-3000)*0.1)+1) PERIODIC=NO LABEL=stp4\n" );
  log_printf( "CUSTOM ARG=tim FUNC=0.5*(erf((4000-x)*0.1)+1) PERIODIC=NO LABEL=stp5\n" );
  log_printf( "CUSTOM ARG=tim FUNC=0.5*(erf((x-4000)*0.1)+1) PERIODIC=NO LABEL=stp6\n" );
  log_printf( "CUSTOM ARG=tim FUNC=0.5*(erf((5000-x)*0.1)+1) PERIODIC=NO LABEL=stp7\n" );
  log_printf( "CUSTOM ARG=tim FUNC=0.5*(erf((x-5000)*0.1)+1) PERIODIC=NO LABEL=stp8\n" );
  log_printf( "CUSTOM ARG=stp1,stp2,stp3,stp4,stp5,stp6,stp7,stp8 VAR=x1,x2,x3,x4,x5,x6,x7,x8 FUNC=0.15*x1+0.14*x2*x3+0.13*x4*x5+0.12*x6*x7+0.11*x8 PERIODIC=NO LABEL=max_z\n" );
  log_printf( "CUSTOM ARG=path_z,max_z FUNC=(y-x)^-1 PERIODIC=NO LABEL=uwz_arg\n" );
  log_printf( "BIASVALUE ARG=uwz_arg LABEL=uw_path_z\n" );
  log_printf( "UPPER_WALLS ARG=path.sss AT=%d KAPPA=1000 LABEL=uw_s\n", n_path_frame[1] );
  log_printf( "LOWER_WALLS ARG=path.sss AT=1  KAPPA=1000 LABEL=lw_s\n" );
  log_printf( "UPPER_WALLS ARG=drr.path.sss_fict AT=%d KAPPA=1000 LABEL=uw_s_fict\n", n_path_frame[1] );
  log_printf( "LOWER_WALLS ARG=drr.path.sss_fict AT=1  KAPPA=1000 LABEL=lw_s_fict\n" );
  log_printf( "PRINT ARG=* STRIDE=1000 FILE=reverse-plumed-force.log FMT=%c12.4f\n", 37 );
  log_printf( "FLUSH STRIDE=1000\n" );
  log_printf( "\"\"\"\n" );

  /* Check for kink in the path: */
  if( fabs( rmsd_avg - rmsd_min ) / rmsd_avg > 0.3 
   || fabs( rmsd_avg - rmsd_max ) / rmsd_avg > 0.3 ) {

    log_printf( "\n WARNING: kink found in the path, try slightly changing rmsd_gap\n\n" );
    log_printf( " rmsd_avg=%6.2f  rmsd_min=%6.2f  rmsd_max=%6.2f\n", rmsd_avg, rmsd_min, rmsd_max );
  }


  /* Done */

  elapsed = wallclock() - elapsed;

  log_printf( " \n Total execution time %.1f sec\n\n", elapsed );

#if 0
 test:

  TRY( read_reference_structures( n_atom, func_args, input_file ) );  /* read first and last frame, copy them to func_args->min[1,2]_xyz */
  COPY( xyz, func_args->min1_xyz, n_dim );
  TRY( calc_derivs( sanm_energy, n_dim, xyz, energy, grad, func_args ) );
  TRY( debug_derivs(sanm_energy, n_dim, xyz, energy, grad, func_args ) );
#endif

 cleanup:

  if( input_file   !=NULL ) fclose( input_file );
  if( output_file  !=NULL ) fclose( output_file );
  if( plumed_file  !=NULL ) fclose( plumed_file );
  if( output_file2 !=NULL ) fclose( output_file2 );
  if( plumed_file2 !=NULL ) fclose( plumed_file2 );
  FREE( xyz );
  FREE( grad );
  FREE( func_args->min1_xyz );
  FREE( func_args->min2_xyz );
  FREE( func_args->posre_xyz );
  FREE( func_args->r0_min1 );
  FREE( func_args->r0_min2 );
  FREE( func_args->pdb_file_image );
  for( i=0; i< 2; i++ ) FREE( ref_confs_min[i] );
  for( i=0; i< 2; i++ )
  for( j=0; j< n_path_frame[i]; j++ ) FREE( path_frames[i*MAXFRAMES+j] );

  return 0;
}

/*
    gcc -o CA-transition-path.exe -O3 CA-transition-path.c common.c -lm
*/
