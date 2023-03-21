/* Written by Michael Levitt, Andras Aszodi, and Istvan Kolossvary */


#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


#define SCALE 1.0
#define CBOND 1.5
#define OBOND 1.8
#define HBOND 1.0
#define NTOTRESTYP 33
#define CABOND_MIN_PDB     1.5  /* min CA-CA found in PDB 2.6 Angs */
#define CABOND_MAX_PDB     4.5  /* max CA-CA found in PDB 4.2 Angs */
#define MAXCHAIN 20
#define MAXBUF 150
#define ZERO 0.0
#define NO  0
#define YES 1


static double sqrarg;
#define SQR(a) ((sqrarg=(a)) == 0.0 ? 0.0 : sqrarg*sqrarg)


struct resnam{
    char *std_3name;
    char  std_1name;
    int   res_id;
} residue[] = {
    "CYS",'c',0,
    "MET",'m',1,
    "PHE",'f',2,
    "ILE",'i',3,
    "LEU",'l',4,
    "VAL",'v',5,
    "TRP",'w',6,
    "TYR",'y',7,
    "ALA",'a',8,
    "GLY",'g',9,
    "THR",'t',10,
    "SER",'s',11,
    "GLN",'q',12,
    "ASN",'n',13,
    "GLU",'e',14,
    "ASP",'d',15,
    "HIS",'h',16,
    "ARG",'r',17,
    "LYS",'k',18,
    "PRO",'p',19,
    "CSH",'c',20,
    "CYH",'c',21,
    "CSM",'c',22,
    "CYX",'c',23,
    "CYM",'c',24,
    "MSE",'m',25,
    "HID",'h',26,
    "HIE",'h',27,
    "HIP",'h',28,
    "TRY",'w',29,
    "CPR",'p',30,
    "PCA",'q',31,
    "SAC",'s',32
};


struct coord{
    double x;
    double y;
    double z;
};


int *sequence = NULL;
char *seq_string = NULL;
struct coord *pdb_chain = NULL;
int chain_endpoints[2*MAXCHAIN];


void
add_cg(struct coord *ca, struct coord *cb, struct coord *n, struct coord *cd, struct coord *cg)
{
    cg->x = (cb->x + cd->x) * 0.5 + (cb->x - ca->x + cd->x - n->x) * 0.2;
    cg->y = (cb->y + cd->y) * 0.5 + (cb->y - ca->y + cd->y - n->y) * 0.2;
    cg->z = (cb->z + cd->z) * 0.5 + (cb->z - ca->z + cd->z - n->z) * 0.2;
}


void
add_cb(struct coord *n, struct coord *ca, struct coord *c, struct coord *cb)
{
double cbond=CBOND,ang,su,sv;
struct coord nca,cca,u,v;
    nca.x = ca->x - n->x;
    nca.y = ca->y - n->y;
    nca.z = ca->z - n->z;
    cca.x = ca->x - c->x;
    cca.y = ca->y - c->y;
    cca.z = ca->z - c->z;
    u.x = nca.x + cca.x;
    u.y = nca.y + cca.y;
    u.z = nca.z + cca.z;
    v.x = nca.y * cca.z - nca.z * cca.y;
    v.y = nca.z * cca.x - nca.x * cca.z;
    v.z = nca.x * cca.y - nca.y * cca.x;
    ang = acos(-1.0) / 2.0 - asin(1.0 / sqrt(3.0));
    su = cbond * cos(ang) / sqrt( SQR(u.x) + SQR(u.y) + SQR(u.z) );
    sv = cbond * sin(ang) / sqrt( SQR(v.x) + SQR(v.y) + SQR(v.z) );
    cb->x = ca->x + u.x * su + v.x * sv;
    cb->y = ca->y + u.y * su + v.y * sv;
    cb->z = ca->z + u.z * su + v.z * sv;
}


void
add_mai_c(struct coord *nca, struct coord *ca, struct coord *cca, struct coord *c ,
          struct coord  *n , struct coord *o , struct coord  *h , struct coord *cb)
/*
    Add main chain atoms after CA.
*/
{
double scale=SCALE,fac,fac_cb;
struct coord cent,u,v,uxv;
    cent.x = (ca->x + cca->x) / 2.0;
    cent.y = (ca->y + cca->y) / 2.0;
    cent.z = (ca->z + cca->z) / 2.0;
    u.x = nca->x - ca->x; v.x = cca->x - ca->x;
    u.y = nca->y - ca->y; v.y = cca->y - ca->y;
    u.z = nca->z - ca->z; v.z = cca->z - ca->z;

    uxv.x = u.y * v.z - u.z * v.y;
    uxv.y = u.z * v.x - u.x * v.z;
    uxv.z = u.x * v.y - u.y * v.x;
    fac = scale / sqrt( SQR(uxv.x) + SQR(uxv.y) + SQR(uxv.z) );

    n->x = cent.x + v.x * 0.125;
    n->y = cent.y + v.y * 0.125;
    n->z = cent.z + v.z * 0.125;
    c->x = cent.x - v.x * 0.125 + uxv.x * 0.5 * fac;
    c->y = cent.y - v.y * 0.125 + uxv.y * 0.5 * fac;
    c->z = cent.z - v.z * 0.125 + uxv.z * 0.5 * fac;
    o->x = cent.x + uxv.x * 2.0 * fac;
    o->y = cent.y + uxv.y * 2.0 * fac;
    o->z = cent.z + uxv.z * 2.0 * fac;
    h->x = cent.x - uxv.x * fac;
    h->y = cent.y - uxv.y * fac;
    h->z = cent.z - uxv.z * fac;
}                             


void
add_mai_n(struct coord *nca, struct coord *ca, struct coord *cca, struct coord *c ,
          struct coord  *n , struct coord *o , struct coord  *h , struct coord *cb, int pro)
/*
    Add main chain atoms before CA.
*/
{
double cbond=CBOND,obond=OBOND,hbond=HBOND,scale=SCALE,fac;
struct coord cent,u,v,uxv;
    cent.x = (ca->x + nca->x) / 2.0;
    cent.y = (ca->y + nca->y) / 2.0;
    cent.z = (ca->z + nca->z) / 2.0;
    u.x = nca->x - ca->x; v.x = cca->x - ca->x;
    u.y = nca->y - ca->y; v.y = cca->y - ca->y;
    u.z = nca->z - ca->z; v.z = cca->z - ca->z;
    if( pro ){
        add_cb(nca,ca,cca,cb);
        uxv.x = cb->x - ca->x;
        uxv.y = cb->y - ca->y;
        uxv.z = cb->z - ca->z;
        fac = scale / cbond;
        hbond = obond;
    }
    else{
        uxv.x = u.y * v.z - u.z * v.y;
        uxv.y = u.z * v.x - u.x * v.z;
        uxv.z = u.x * v.y - u.y * v.x;
        fac = scale / sqrt( SQR(uxv.x) + SQR(uxv.y) + SQR(uxv.z) );
    }
    n->x = cent.x - u.x * 0.125 + uxv.x * 0.25 * fac;
    n->y = cent.y - u.y * 0.125 + uxv.y * 0.25 * fac;
    n->z = cent.z - u.z * 0.125 + uxv.z * 0.25 * fac;
    c->x = cent.x + u.x * 0.125 - uxv.x * 0.50 * fac;
    c->y = cent.y + u.y * 0.125 - uxv.y * 0.50 * fac;
    c->z = cent.z + u.z * 0.125 - uxv.z * 0.50 * fac;
    o->x = cent.x - uxv.x * obond * fac;
    o->y = cent.y - uxv.y * obond * fac;
    o->z = cent.z - uxv.z * obond * fac;
    h->x = cent.x + uxv.x * hbond * fac;
    h->y = cent.y + uxv.y * hbond * fac;
    h->z = cent.z + uxv.z * hbond * fac;
}                             


void
build_main_chain(struct coord *pdb_ca, struct coord *pdb_n, struct coord *pdb_c,
                 struct coord *pdb_o,  struct coord *pdb_h, struct coord *pdb_cb, int start, int end)
{
int i;
struct coord dummy;
    dummy.x = pdb_ca[start].x + pdb_ca[start+1].x - pdb_ca[start+2].x;
    dummy.y = pdb_ca[start].y + pdb_ca[start+1].y - pdb_ca[start+2].y;
    dummy.z = pdb_ca[start].z + pdb_ca[start+1].z - pdb_ca[start+2].z;
    add_mai_n(&dummy, &pdb_ca[start], &pdb_ca[start+1],
              &dummy, &pdb_n[ start],
              &dummy, &pdb_h[ start],
                      &pdb_cb[start],
              NO);
    for(i=start+2; i<=end; i++)
    {
        add_mai_n(&pdb_ca[i-2], &pdb_ca[i-1], &pdb_ca[i],
                  &pdb_c[ i-2], &pdb_n[ i-1],
                  &pdb_o[ i-2], &pdb_h[ i-1],
                                &pdb_cb[i-1],
                  (residue[sequence[i-1]].res_id == 19 || residue[sequence[i-1]].res_id == 30) ? YES : NO);
        add_cb(   &pdb_n[ i-2],
                  &pdb_ca[i-2],
                  &pdb_c[ i-2],
                  &pdb_cb[i-2]);
    }
    add_mai_c(&pdb_ca[end-2],  &pdb_ca[end-1], &pdb_ca[end],
                               &pdb_c[ end-1], &pdb_n[ end],
                               &pdb_o[ end-1], &pdb_h[ end],
                               &pdb_cb[end-1]);
    add_cb(   &pdb_n[ end-1],
              &pdb_ca[end-1],
              &pdb_c[ end-1],
              &pdb_cb[end-1]);
    dummy.x = pdb_ca[end].x + pdb_ca[end-1].x - pdb_ca[end-2].x;
    dummy.y = pdb_ca[end].y + pdb_ca[end-1].y - pdb_ca[end-2].y;
    dummy.z = pdb_ca[end].z + pdb_ca[end-1].z - pdb_ca[end-2].z;
    add_mai_c(&pdb_ca[end-1], &pdb_ca[end], &dummy,
                               &pdb_c[ end], &dummy,
                               &pdb_o[ end], &dummy,
                               &pdb_cb[end]);
    add_cb(   &pdb_n[ end],
              &pdb_ca[end],
              &pdb_c[ end],
              &pdb_cb[end]);
}


void
write_pdb_file(char *fname, struct coord *xyz, int nres, int nchain)
{
struct coord *pdb_n,*pdb_c,*pdb_o,*pdb_h,*pdb_cb,pdb_cg,pdb_oxt,unit,proj;
int i, j, k, chain_start, chain_end;
static int model=0;
double scale;
static FILE *pdb_file=NULL;
    pdb_n = (struct coord *)malloc(nres*sizeof(struct coord));
    pdb_c = (struct coord *)malloc(nres*sizeof(struct coord));
    pdb_o = (struct coord *)malloc(nres*sizeof(struct coord));
    pdb_h = (struct coord *)malloc(nres*sizeof(struct coord));
    pdb_cb = (struct coord *)malloc(nres*sizeof(struct coord));
/*
    Add main chain and CB to CA trace for visualization:
*/
    for(i=0; i<nchain; i++) 
    {
        chain_start = chain_endpoints[i*2];
        chain_end   = chain_endpoints[i*2+1];
        build_main_chain(xyz,pdb_n,pdb_c,pdb_o,pdb_h,pdb_cb,chain_start,chain_end);
    }
    if( pdb_file == NULL )
    {
        if((pdb_file = fopen(fname,"w")) == NULL)
        {
            fprintf(stderr, " ERROR: Cannot open %s for w - Aborting.\n",fname);
            exit(0);
        }
    }
    fprintf( pdb_file, "MODEL %8d\n", ++model );
    for( i=0, j=0; i<nres; i++ )
    {
        fprintf(pdb_file,"ATOM  %5d  N   %3s  %4d    %8.3f%8.3f%8.3f\n",
            ++j,residue[sequence[i]].std_3name,i+1,
            pdb_n[i].x,pdb_n[i].y,pdb_n[i].z);                                          /* N  */
        fprintf(pdb_file,"ATOM  %5d  CA  %3s  %4d    %8.3f%8.3f%8.3f\n",
            ++j,residue[sequence[i]].std_3name,i+1,
            xyz[i].x,xyz[i].y,xyz[i].z);                                                /* CA */
        fprintf(pdb_file,"ATOM  %5d  C   %3s  %4d    %8.3f%8.3f%8.3f\n",
            ++j,residue[sequence[i]].std_3name,i+1,
            pdb_c[i].x,pdb_c[i].y,pdb_c[i].z);                                          /* C  */
        fprintf(pdb_file,"ATOM  %5d  O   %3s  %4d    %8.3f%8.3f%8.3f\n",
            ++j,residue[sequence[i]].std_3name,i+1,
            pdb_o[i].x,pdb_o[i].y,pdb_o[i].z);                                          /* O  */
        if( residue[sequence[i]].res_id != 9 ){                             /* !GLY */
            fprintf(pdb_file,"ATOM  %5d  CB  %3s  %4d    %8.3f%8.3f%8.3f\n",
                ++j,residue[sequence[i]].std_3name,i+1,
                pdb_cb[i].x,pdb_cb[i].y,pdb_cb[i].z);                                   /* CB */
        }
#if 0
        if( residue[sequence[i]].res_id == 19 || residue[sequence[i]].res_id == 30 ){   /*  PRO */
            add_cg(&xyz[i], &pdb_cb[i], &pdb_n[i], &pdb_h[i], &pdb_cg);
            fprintf(pdb_file,"ATOM  %5d  CG  %3s  %4d    %8.3f%8.3f%8.3f\n",
                ++j,residue[sequence[i]].std_3name,i+1,
                pdb_cg.x,pdb_cg.y,pdb_cg.z);                                            /* CG */
            fprintf(pdb_file,"ATOM  %5d  CD  %3s  %4d    %8.3f%8.3f%8.3f\n",
                ++j,residue[sequence[i]].std_3name,i+1,
                pdb_h[i].x,pdb_h[i].y,pdb_h[i].z);                                      /* CD */
        }
#endif
        for( k=0; k<nchain; k++ ) {
            if( i == chain_endpoints[k*2+1] ) {  /* end of chain */
            unit.x = pdb_c[i].x - xyz[i].x;
            unit.y = pdb_c[i].y - xyz[i].y;
            unit.z = pdb_c[i].z - xyz[i].z;
                scale = unit.x * unit.x + unit.y * unit.y + unit.z * unit.z;
                scale = sqrt(scale);
                unit.x /= scale;
                unit.y /= scale;
                unit.z /= scale;
                scale = unit.x * (pdb_o[i].x - pdb_c[i].x)
                      + unit.y * (pdb_o[i].y - pdb_c[i].y)
                      + unit.z * (pdb_o[i].z - pdb_c[i].z);
                proj.x = unit.x * scale;
                proj.y = unit.y * scale;
                proj.z = unit.z * scale;
                pdb_oxt.x = 2. * (pdb_c[i].x + proj.x) - pdb_o[i].x;
                pdb_oxt.y = 2. * (pdb_c[i].y + proj.y) - pdb_o[i].y;
                pdb_oxt.z = 2. * (pdb_c[i].z + proj.z) - pdb_o[i].z;
                fprintf(pdb_file,"ATOM  %5d  OXT %3s  %4d    %8.3f%8.3f%8.3f\n",
                    ++j,residue[sequence[i]].std_3name,i+1,
                    pdb_oxt.x,pdb_oxt.y,pdb_oxt.z);                                         /* OXT */
                //fprintf( pdb_file, "TER\n" );
                break;
            }
        }
    }
    fprintf( pdb_file, "ENDMDL\n" );
    free(pdb_n); free(pdb_c); free(pdb_o); free(pdb_h); free(pdb_cb);
}


void
read_pdb_file(char *fname, int * nres_out, int * nchain_out, double user_defined_chain_gap, int * duplicate, int * eof)
{
int i, j, ires, error=NO;
static int nres=0, nchain=1, model=0, cnt=0;
char buf[MAXBUF],test1[7],test2[5],res[4];
double dist, x, y, z, x_save, y_save, z_save;
static FILE *pdb_file=NULL;
    chain_endpoints[0] = 0;  /* start */
    if( pdb_file == NULL )
    {
        if((pdb_file = fopen(fname,"r")) == NULL)
        {
            fprintf(stderr, " ERROR: Cannot open %s - Aborting.\n",fname);
            exit(0);
        }
        nres = 0;
        while(fgets(buf,MAXBUF,pdb_file) != NULL)
        {
            if( (strncmp( buf, "TER", 3 ) == 0) || (strncmp( buf, "MODEL", 5 ) == 0))  /* skip TER and MODEL cards */
                continue;
            if( strncmp( buf, "ENDMDL", 6 ) == 0 || strncmp( buf, "END", 3 ) == 0 )
                break;
            strncpy(test1,buf,6); test1[6]='\0';
            strncpy(test2,buf+12,4); test2[4]='\0';
            if( (!strcmp(test1,"ATOM  ")) &&                                        /* CA atom found */
                (!strcmp(test2,"CA  "  ) || !strcmp(test2," CA "  )  ||             /* PDB standard horrible */
                 !strcmp(test2,"  CA" )) && (buf[16] == ' ' || buf[16] == 'A') ){   /* Use altLoc A  */

                nres++;
                sscanf(buf+30,"%lf %lf %lf", &x, &y, &z);  /* read coords for chain gap check */
            }
            if( nres > 1 ){
                dist = ZERO;
                dist += SQR(x - x_save);
                dist += SQR(y - y_save);
                dist += SQR(z - z_save);
                dist = sqrt(dist);
                if( dist > user_defined_chain_gap ){  /* end of chain */

                    chain_endpoints[2*(nchain-1)+1] = nres-2;  /* end previous */
                    chain_endpoints[2*nchain]       = nres-1;  /* start next   */
                    nchain++;
                }
            }
            x_save = x;
            y_save = y;
            z_save = z;
        }
        chain_endpoints[2*(nchain-1)+1] = nres-1;  /* end */
        if(ferror(pdb_file))
        {
            fprintf(stderr, " ERROR: Input error in %s - Aborting.\n",fname);
            exit(0);
        }
        rewind(pdb_file);
        if( sequence   == NULL ) sequence   = (int *)malloc(nres*sizeof(int));
        if( seq_string == NULL ) seq_string = (char *)malloc((nres+1)*sizeof(char));
        if( pdb_chain  == NULL ) pdb_chain  = (struct coord *)malloc(nres*sizeof(struct coord));
        *nres_out   = nres;
        *nchain_out = nchain;
        return;
    }
    ires = 0;
    model++;
    while(fgets(buf,MAXBUF,pdb_file) != NULL)
    {
        if( (strncmp( buf, "TER", 3 ) == 0) || (strncmp( buf, "MODEL", 5 ) == 0))  /* skip TER and MODEL cards */
            continue;
        if( strncmp( buf, "ENDMDL", 6 ) == 0 || strncmp( buf, "END", 3 ) == 0 )
        {
            cnt++;
            if( cnt>1 ) {*duplicate = YES;} else {*duplicate = NO;}
            break;
        }
        cnt = 0;

        strncpy(test1,buf,6); test1[6]='\0';
        strncpy(test2,buf+12,4); test2[4]='\0';
        if( (!strcmp(test1,"ATOM  ")) &&                                        /* CA atom found */
            (!strcmp(test2,"CA  "  ) || !strcmp(test2," CA "  )  ||             /* PDB standard horrible */
             !strcmp(test2,"  CA" )) && (buf[16] == ' ' || buf[16] == 'A') )    /* Use altLoc A  */
        {
            strncpy(res,buf+17,3); res[3]='\0';
            for(i=0; i<NTOTRESTYP; i++)
            {
                if( !strcmp(res,residue[i].std_3name) )
                {
                    seq_string[ires] = residue[i].std_1name;
                    sequence[ires] = residue[i].res_id;
                    sscanf(buf+30,"%lf %lf %lf",&pdb_chain[ires].x,
                                                &pdb_chain[ires].y,
                                                &pdb_chain[ires].z);
                    if( ires >= 1 ){
                        dist = ZERO;
                        dist += SQR(pdb_chain[ires].x - pdb_chain[ires-1].x);
                        dist += SQR(pdb_chain[ires].y - pdb_chain[ires-1].y);
                        dist += SQR(pdb_chain[ires].z - pdb_chain[ires-1].z);
                        dist = sqrt(dist);
                        if( dist < CABOND_MIN_PDB || dist > CABOND_MAX_PDB ) {
                            error = YES;
                            for( j=0; j<2*nchain; j++ ) {
                                if( ires==chain_endpoints[j] ) { error = NO; break; }
                            }
                            if( error ) {
                                fprintf(stderr," ERROR: Distorted or missing residue in position %d, model %d,\n",ires+1,model);
                                fprintf(stderr,"        CA-CA dist= %5.2f Angs.\n",dist);
                            }
                        }
                    }
                    ires++; /* Read next residue. */
                    break;
                }
            }
            if( i == NTOTRESTYP ){ 
                fprintf(stderr," ERROR: Unknown residue %s - Aborting.\n", res);
                free(sequence); free(seq_string); free(pdb_chain);
                exit(0);
            }
        }
    }
    if(ferror(pdb_file))
    {
        fprintf(stderr, " ERROR: Input error in %s - Aborting.\n",fname);
        exit(0);
    }
    if( feof( pdb_file ) ) {

        *eof = YES;
    }
    seq_string[ires] = '\0';
    if( error ) exit(0);
}

int
main(int argc, char *argv[])
{
int nres, nchain, model, duplicate, eof=NO;
double chain_gap;

    if( argc != 4 ) {
    fprintf(stderr, "Usage:  %s  CA-pdb-fname  backbone-pdb-fname  chain_gap\n", argv[0] );
    return 1;
  }


    chain_gap = atof(argv[3]);
    read_pdb_file(argv[1], &nres, &nchain, chain_gap, &duplicate, &eof);
    for( ;; )
    {
        read_pdb_file(argv[1], &nres, &nchain, chain_gap, &duplicate, &eof);
        if( eof )
            break;
        if( !duplicate ) write_pdb_file(argv[2], pdb_chain, nres, nchain);
    }
    free(sequence); free(seq_string); free(pdb_chain);
    return 0;
}
