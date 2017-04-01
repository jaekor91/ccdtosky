/* spherematchDJE -- Daniel Eisenstein, June 2015

This code utilizes KD-tree infrastructure from the UW HPCC Smooth code
by Joachim Stadel.  The matching infrastructure written by Daniel Eisenstein.

This code will match two lists of objects, either in 2-d or 3-d.
Returns a list of id1, id2, dist12.

Includes a listing of parameters in lines starting with #.
Output can be sent to stdout or to a file (-o option).

If run angularly (-ang), then all distances are in units of arcseconds.
Can indicate (explicitly by the -auto flag) that the two lists are
the same and therefore limit the outputs to one listing per pair,
with ID1<ID2.

Input files are simple ASCII, with each object on its own line.
Either x,y,z or RA,Dec must be the first columns on the line.
Lines starting with '#' are skipped as comments.
Edit get_one_object() if you don't like this format.

The first input file is put into the KD-tree, so this should be the larger file.
This file is actually read twice, so it cannot be a pipe.

The second file is then streamed through to search for neighbors in the tree.
If the file name is -, then this input is taken from stdin.

*/



#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <ctype.h>
#include <assert.h>
#include "kd.h"
#include "smooth.h"

/* ================================================================ */

int global_angular;  // ==1 if this is RA/Dec and should be pushed to the unit sphere
int global_auto;   // ==1 if this is a auto-count and should require id1<id2

#define ARCSEC (3600.0*180.0/M_PI)

int get_one_object(FILE *fp, double pos[3]) {
    char line[200];
    while (fgets(line,200,fp)!=NULL) {
        if (line[0]=='#') continue;

	if (global_angular) {
	    double ra, dec;
	    int nread = sscanf(line,"%lf %lf", &ra, &dec);
	    assert(nread==2);
	    assert(dec<90.00000001&&dec>-90.0000001);
	    // We don't assert on RA because it might be useful to wrap the sphere.
	    ra *= M_PI/180.0;
	    dec *= M_PI/180.0;
	    pos[0] = ARCSEC*cos(ra)*cos(dec);
	    pos[1] = ARCSEC*sin(ra)*cos(dec);
	    pos[2] = ARCSEC*sin(dec);
	    // We use a sphere of radius 3600*180/PI so that the angular unit is arcsec
	} else {
	    int nread = sscanf(line,"%lf %lf %lf", pos, pos+1, pos+2);
	    assert(nread==3);
	}

	return 1;
    }
    return 0;
}

/* ============= Read the second input file to find matches ======= */

void FindMatches(SMX smx, char filename[], double match_radius, double min_radius) {
    FILE *fp;
    if (strcmp(filename,"-")) {
	printf("# Searching for neighbors from %s\n", filename);
	fp = fopen(filename,"r");
	assert(fp!=NULL);
    } else {
	printf("# Searching for neighbors from stdin\n");
        fp = stdin;
    }

    double pos[3], r2 = match_radius*match_radius;
    double min_radius2 = min_radius*min_radius;
    if (min_radius<0) min_radius2 = -1.0;   // This will never fail

    int id=0, npair=0;
    printf("# Using radius %lf, min %lf\n", match_radius, min_radius);
    if (global_auto) 
    	printf("# Auto-correlation: requiring ID1<ID2 (and no-self-matches)\n");

    while (get_one_object(fp,pos)) {
	int nSmooth = smBallGather(smx,r2,pos);
	for (int j=0;j<nSmooth;j++) {
	    if (smx->fList[j]<min_radius2) continue;
	    	// Omit very close pairs, as these are probably self-pairs
	    PARTICLE *p = smx->kd->p+smx->pList[j];
	    if (global_auto && p->iOrder>=id) continue;
	    	// Omit the duplicate copy
	    
	    npair++;
	    printf("%d %d %lf\n", p->iOrder, id, sqrt(smx->fList[j]));
	    // printf("  %lf %lf %lf   %lf %lf %lf\n", pos[0], pos[1], pos[2], p->r[0], p->r[1], p->r[2]);
	}
	id++;
    }
    printf("# Searched %d neighbors, found %d pairs\n", id, npair);
    fclose(fp);
}

/* ====================  Read the first input file and put it in KD ========= */

int kdReadFile(KD kd,char infile[])
{
	FILE *fp = fopen(infile,"r");
	assert(fp!=NULL);
	double pos[3];

	/* Count particles in input file */
	kd->nActive = 0;
	while (get_one_object(fp,pos)) {
	    kd->nActive++;
	}
	printf("# Expect %d objects from file %s\n", kd->nActive, infile);
	rewind(fp);
	fflush(NULL);

	/*
	 ** Allocate particles.
	 */
	kd->p = (PARTICLE *)malloc(kd->nActive*sizeof(PARTICLE));
	assert(kd->p != NULL);

	/*
	 ** Read Stuff!
	 */
	int nCnt = 0;
	while (get_one_object(fp,pos)) {
	    kd->p[nCnt].iOrder = nCnt;
	    for (int j=0;j<3;++j) kd->p[nCnt].r[j] = pos[j];
	    ++nCnt;
	}
	assert(nCnt==kd->nActive);
	printf("# Read %d objects from file %s\n", kd->nActive, infile);

	/*
	 ** Calculate Bounds.
	 */
	BND bnd;
	for (int j=0;j<3;++j) {
		bnd.fMin[j] = kd->p[0].r[j];
		bnd.fMax[j] = kd->p[0].r[j];
		}
	for (int i=1;i<kd->nActive;++i) {
		for (int j=0;j<3;++j) {
			if (bnd.fMin[j] > kd->p[i].r[j]) 
				bnd.fMin[j] = kd->p[i].r[j];
			else if (bnd.fMax[j] < kd->p[i].r[j])
				bnd.fMax[j] = kd->p[i].r[j];
			}
		}
	kd->bnd = bnd;
	return(kd->nActive);
}

/* ================================================================ */

void usage(char error[])
{
	fprintf(stderr,"\nUSAGE: \n");
	fprintf(stderr," spherematchDJE infile1 infile2 [-o <Output Name>] \n");
	fprintf(stderr,"       Using infile2 of '-' will stream from stdin\n");
	fprintf(stderr,"       Omitting -o will send output to stdout\n");
	fprintf(stderr,"       Outputs are (ID1,ID2,dist12), zero-indexed.\n");
	fprintf(stderr,"   [-ang] use the unit sphere.  RA/Dec in dec, radius in arcsec\n");
	fprintf(stderr,"   [-r <Match Radius>]\n");
	fprintf(stderr,"   [-auto] to require ID1<ID2 and avoid duplicate pairs\n");
	fprintf(stderr,"   [-rmin <Match Radius>] Minimum radius\n");
	fprintf(stderr,"\n Advanced options:\n");
	fprintf(stderr,"   [-max <maxMatch>] the maximum number of matches per object\n");
	fprintf(stderr,"     Code will crash if this is violated; not a silent failure\n");
	fprintf(stderr,"   [-b <nBucket>]\n");
	fprintf(stderr,"   [-p <xyzPeriod>] (will default to open boundary conditions)\n");
	fprintf(stderr,"ERROR: %s\n\n", error);
	exit(1);
	}

int main(int argc,char **argv)
{
	KD kd;
	SMX smx;
	FILE *fp;
	double fPeriod[3];
	for (int j=0;j<3;++j) fPeriod[j] = HUGE;
	char infile1[200], infile2[200], outfile[200];
	strcpy(infile1,"");
	strcpy(infile2,"");
	strcpy(outfile,"");
	double match_radius=0.0;
	double min_radius=-1.0;   // Off by default
	int infiles=0;
   	int nBucket = 16;
   	int maxMatch = 128;
	global_angular = 0;
	global_auto = 0;

	int i = 1;
	while (i < argc) {
		if (!strcmp(argv[i],"-b")) {
		    ++i; if (i >= argc) usage("Missing parameter after -b.");
		    nBucket = atoi(argv[i++]);
		}
		else if (!strcmp(argv[i],"-max")) {
		    ++i; if (i >= argc) usage("Missing parameter after -max.");
		    maxMatch = atoi(argv[i++]);
		}
		else if (!strcmp(argv[i],"-o")) {
		    ++i; if (i >= argc) usage("Missing parameter after -o.");
		    strcpy(outfile,argv[i++]);
		}
		else if (!strcmp(argv[i],"-r")) {
		    ++i; if (i >= argc) usage("Missing parameter after -r.");
		    match_radius = atof(argv[i++]);
		}
		else if (!strcmp(argv[i],"-rmin")) {
		    ++i; if (i >= argc) usage("Missing parameter after -rmin.");
		    min_radius = atof(argv[i++]);
		}
		else if (!strcmp(argv[i],"-ang")) {
		    global_angular = 1;
		    ++i; 
		}
		else if (!strcmp(argv[i],"-auto")) {
		    global_auto = 1;
		    ++i; 
		}
		else if (!strcmp(argv[i],"-p")) {
		    ++i; if (i >= argc) usage("Missing parameter after -p.");
		    fPeriod[0] = atof(argv[i]);
		    fPeriod[1] = atof(argv[i]);
		    fPeriod[2] = atof(argv[i]);
		    ++i;
		}
		else if (infiles<2) {
		    // printf("# File %d: %s\n", infiles, argv[i]);
		    if (infiles==0) 
		        strcpy(infile1,argv[i++]);
		    else
		        strcpy(infile2,argv[i++]);
		    infiles++;
		}
		else usage("Unknown parameter option.");
	}

	if (argc==1) usage("No parameters given.");
        if (!(match_radius>0)) 
	    usage("-r must be a positive number.");
        if (!(min_radius>=0||min_radius==-1.0)) 
	    usage("If used, -rmin must be a non-negative number.");
        if (!(strcmp(infile1,""))||!(strcmp(infile2,""))) 
	    usage("Must supply two input files names.");

	if (strcmp(outfile,"")) {
	    freopen(outfile,"w",stdout);
	}
	if (global_angular) {
	    printf("# Doing an angular match on the sphere in arcsec\n");
	} else { 
	    printf("# Doing a 3-d match in user-supplied units\n");
	}

	kdInit(&kd,nBucket);
	kdReadFile(kd,infile1);
	kdBuildTree(kd);
	smInit(&smx,kd,maxMatch,fPeriod);

	FindMatches(smx,infile2,match_radius,min_radius);

	// kdOrder(kd);
	smFinish(smx);
	kdFinish(kd);
	fclose(stdout);
	return 0;
}


