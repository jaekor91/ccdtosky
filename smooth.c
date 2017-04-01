#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "smooth.h"
#include "kd.h"


int smInit(SMX *psmx,KD kd,int nSmooth,double *fPeriod)
{
	// nSmooth will allocate the max number of neighbors per search
	SMX smx;
	int j;

	// assert(nSmooth <= kd->nActive);

	smx = (SMX)malloc(sizeof(struct smContext));
	assert(smx != NULL);
	smx->kd = kd;
	smx->nSmooth = nSmooth;

	smx->nListSize = smx->nSmooth;
	smx->fList = (double *)malloc(smx->nListSize*sizeof(double));
	assert(smx->fList != NULL);
	smx->pList = (int *)malloc(smx->nListSize*sizeof(int));
	assert(smx->pList != NULL);
	/*
	 ** Set for Periodic Boundary Conditions.
	 */
	for (j=0;j<3;++j) smx->fPeriod[j] = fPeriod[j];
	*psmx = smx;	
	return(1);
	}


void smFinish(SMX smx)
{
	free(smx);
}

int smBallGather(SMX smx,double fBall2,double *ri)
{
	KDN *c;
	PARTICLE *p;
	int pj,nCnt,cp,nSplit;
	double dx,dy,dz,x,y,z,lx,ly,lz,sx,sy,sz,fDist2;

	c = smx->kd->kdNodes;
	p = smx->kd->p;
	nSplit = smx->kd->nSplit;
	lx = smx->fPeriod[0];
	ly = smx->fPeriod[1];
	lz = smx->fPeriod[2];
	x = ri[0];
	y = ri[1];
	z = ri[2];
	nCnt = 0;
	cp = ROOT;
	while (1) {
		INTERSECT(c,cp,fBall2,lx,ly,lz,x,y,z,sx,sy,sz);
		/*
		 ** We have an intersection to test.
		 */
		if (cp < nSplit) {
			cp = LOWER(cp);
			continue;
			}
		else {
			for (pj=c[cp].pLower;pj<=c[cp].pUpper;++pj) {
				dx = sx - p[pj].r[0];
				dy = sy - p[pj].r[1];
				dz = sz - p[pj].r[2];
				fDist2 = dx*dx + dy*dy + dz*dz;
				if (fDist2 < fBall2) {
					smx->fList[nCnt] = fDist2;
					smx->pList[nCnt++] = pj;
					}
				}
			}
	GetNextCell:
		SETNEXT(cp);
		if (cp == ROOT) break;
		}
	assert(nCnt <= smx->nListSize);
	return(nCnt);
	}

