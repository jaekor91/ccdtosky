#ifndef SMOOTH_HINCLUDED
#define SMOOTH_HINCLUDED

#include "kd.h"

typedef struct smContext {
	KD kd;
	int nSmooth;
	double fPeriod[3];
	int nListSize;
	double *fList;
	int *pList;
	} * SMX;

int smInit(SMX *,KD,int,double *);
void smFinish(SMX);
int  smBallGather(SMX,double,double *);

#endif



