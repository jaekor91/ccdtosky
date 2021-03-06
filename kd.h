#ifndef KD_HINCLUDED
#define KD_HINCLUDED

#define ROOT		1
#define LOWER(i)	(i<<1)
#define UPPER(i)	((i<<1)+1)
#define PARENT(i)	(i>>1)
#define SIBLING(i) 	((i&1)?i-1:i+1)
#define SETNEXT(i)\
{\
	while (i&1) i=i>>1;\
	++i;\
	}

typedef struct Particle {
	double r[3];
	int iOrder;
	} PARTICLE;

typedef struct bndBound {
	double fMin[3];
	double fMax[3];
	} BND;

typedef struct kdNode {
	double fSplit;
	BND bnd;
	int iDim;
	int pLower;
	int pUpper;
	} KDN;

typedef struct kdContext {
	int nBucket;
	int nActive;
	double fTime;
	BND bnd;
	int nLevels;
	int nNodes;
	int nSplit;
	PARTICLE *p;
	KDN *kdNodes;
	int uSecond;
	int uMicro;
	} * KD;


#define INTERSECT(c,cp,fBall2,lx,ly,lz,x,y,z,sx,sy,sz)\
{\
	double dx,dy,dz,dx1,dy1,dz1,fDist2;\
	dx = c[cp].bnd.fMin[0]-x;\
	dx1 = x-c[cp].bnd.fMax[0];\
	dy = c[cp].bnd.fMin[1]-y;\
	dy1 = y-c[cp].bnd.fMax[1];\
	dz = c[cp].bnd.fMin[2]-z;\
	dz1 = z-c[cp].bnd.fMax[2];\
	if (dx > 0.0) {\
		dx1 += lx;\
		if (dx1 < dx) {\
			fDist2 = dx1*dx1;\
			sx = x+lx;\
			}\
		else {\
			fDist2 = dx*dx;\
			sx = x;\
			}\
		if (fDist2 > fBall2) goto GetNextCell;\
		}\
	else if (dx1 > 0.0) {\
		dx += lx;\
		if (dx < dx1) {\
			fDist2 = dx*dx;\
			sx = x-lx;\
			}\
		else {\
			fDist2 = dx1*dx1;\
			sx = x;\
			}\
		if (fDist2 > fBall2) goto GetNextCell;\
		}\
	else {\
		fDist2 = 0.0;\
		sx = x;\
		}\
	if (dy > 0.0) {\
		dy1 += ly;\
		if (dy1 < dy) {\
			fDist2 += dy1*dy1;\
			sy = y+ly;\
			}\
		else {\
			fDist2 += dy*dy;\
			sy = y;\
			}\
		if (fDist2 > fBall2) goto GetNextCell;\
		}\
	else if (dy1 > 0.0) {\
		dy += ly;\
		if (dy < dy1) {\
			fDist2 += dy*dy;\
			sy = y-ly;\
			}\
		else {\
			fDist2 += dy1*dy1;\
			sy = y;\
			}\
		if (fDist2 > fBall2) goto GetNextCell;\
		}\
	else {\
		sy = y;\
		}\
	if (dz > 0.0) {\
		dz1 += lz;\
		if (dz1 < dz) {\
			fDist2 += dz1*dz1;\
			sz = z+lz;\
			}\
		else {\
			fDist2 += dz*dz;\
			sz = z;\
			}\
		if (fDist2 > fBall2) goto GetNextCell;\
		}\
	else if (dz1 > 0.0) {\
		dz += lz;\
		if (dz < dz1) {\
			fDist2 += dz*dz;\
			sz = z-lz;\
			}\
		else {\
			fDist2 += dz1*dz1;\
			sz = z;\
			}\
		if (fDist2 > fBall2) goto GetNextCell;\
		}\
	else {\
		sz = z;\
		}\
	}


void kdTime(KD,int *,int *);
int kdInit(KD *,int);
int kdBuildTree(KD);
void kdOrder(KD);
void kdFinish(KD);

#endif











