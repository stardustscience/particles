#define NCIRC 21
#define FIRST 1
#define LAST 2
#define MID 3

#define SHALLOW 1
#define DEEP 2
#define LONG_TRACE 800
#define MAX_TRACE_LEN 1000
#define TRACE_OFFSET 100 //

#define NONE 0
#define START 1
#define END 2

#define MASK 1
#define ANTS 2

struct trSect{
	float tx1, ty1, tx2, ty2;  // in model coords 1/10th deg
	float flowVol; // the volume through the transect in Sv
	float shallowFlow, deepFlow;

};

class flowPath;

struct traceLine {
	int start, end;
	float traceX[MAX_TRACE_LEN], traceY[MAX_TRACE_LEN];
	int color;

};

class Layer
{
public:
	Layer(int nr, int nc, int L);
	void setVecs(float *Xvecs, float *Yvecs, float *spd) { vX_S = Xvecs; vY_S = Yvecs; speed_S = spd; };
	//void setVecs_D(float *Xvecs, float *Yvecs, float *spd) { vX_D = Xvecs; vY_D = Yvecs; speed_D = spd; };

	void setFlows(float *fx, float *fy, float *sp, float **N, float **E);
	void setTransect_S(float x1, float y1, float x2, float y2, int n);
	float computeTransect(float x1, float y1, float x2, float y2, int n);
	float computeTransectHR(float x1, float y1, float x2, float y2, float *layers, int n, int mid, int SD); // hi res calculation
	float computeTransectHR_2( int n); // hi res calculation

	float computePeakTransectHR(int n, float &maxp, float &pxMax, float &pyMax); // find the peak
	void measureTransctFlows();

	float computePeakNorth(int Lat, float *layer);
	float computePeakNorth_2(int Lat10, int Lon10, float *agN,int &c, bool north);
	void findCriticalPoint(int sLat,int sLon, int &fLat, int &fLon, bool north);
	void computeLRpath(int t1, int t2); // a path from t1 to t2 in low res.
	
	void drawTransectLine(int n);
	void drawOrthoFlowProfile(float *Layers, int n);

	bool traceToTransect(float xs, float ys, float *xt, float *yt, int t2, float &avgFlow, bool reverse);
	bool traceToTransect_2(float xs, float ys, traceLine *tr, int t2, float &avgFlow, bool reverse);
	void extendTrace_2(traceLine *tr, int p, int n, int randval, bool reverse);
	bool drawToTransect(float xs, float ys, float *xt, float *yt, int t2, bool reverse);
	void drawTrace(traceLine *tr, int rtype);
	void drawTrace_2(traceLine *tr);

	void smoothTrace(traceLine *tr);
	void intersectTransect(int ip, int t2);
	void getTransectFlow(int n, float &s, float &d);

	void computePathline(int tl, int t2, int tp, bool reverse,float flow); // compute a pathline from t1 to t2
	//void computeDrawPathLine(int t1, int t2, int t3, int tpe, float flow); // a two part path

	void setColors(int c,int scm,float &red,float &green,float &blue);
	
	void computePathLoop(int t1, int t2, int tpe, float flowVol);

	void computeDrawPathLine(int t1, int t2, float p, bool dir,int colors, float tr, int startEnd, int rtype, int scm);
	bool computeDrawPathLoop(int t1, int t2, float p, int colors, float tr, int rtype,int scm);
	void computeDrawPathABC(int t1, int t2, int t3, float p, int colors, float tr, int startEndExt, int rtype, int scm);

	void computePath2part(int t1, int t2, int t3, int tpe, float flow); // a two part path
	void getTransectAnchors(int i, float &x1, float &y1, float &x2, float &y2);
																		//void drawLongTrace();
	void drawPaths();
	//float shallowFlow, deepFlow;  // should be part of the struct	
	float startX[50], startY[50];  // the points used for path starts.
private:	
	int nRows, nCols;
	float *transOrth;
	float *smFlowX, *smFlowY, *smSpeed;  // the smoothed flows
	float **fN, **fE; // Northing and Easting the entire model
	float *vX_S, *vY_S, *speed_S;  // this is the unsmoothed flow

	float rn, gn, bn; // the node colors.


	float *smX, *smY; // temps for smoothing
	float VoX, VoY; // unit vector orthogonal to transect
	trSect transects_S[50];

	traceLine *trace_1, *trace_2;

	float netFlow[80]; // by layer

	float *cosCor;

	float circX[NCIRC], circY[NCIRC];

	/*--------------------------*/
	flowPath **flowPaths;
	float *trLongX, *trLongY;
	int longTraceLen;
	int pathCount;
};
