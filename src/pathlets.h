
#define MAXPATHS 30
#define LONG_TRACE 800

class flowPath;

class pathlets
{
public: 
	pathlets(int nr, int nc, int L, int colScm);
	void setVecs(float *Xvecs, float *Yvecs, float *spd) { vX_S = Xvecs; vY_S = Yvecs; speed_S = spd; };
	//void setVecs_D(float *Xvecs, float *Yvecs, float *spd) { vX_D = Xvecs; vY_D = Yvecs; speed_D = spd; };

	void makeTraces();
	void drawAnimatedTraces(float lineWid);
	void drawTracesStatic();
	void drawLongTrace();
	void drawPaths(bool top, bool bottom);
	void drawArrowHead(int si, float *xp, float *yp, float wid);
//	void loadPaths(char *fname);
	//void drawGlyphs();

	void trace(float *tx, float *ty,  float xs, float ys);
	void longTrace(float xs, float ys, bool reverse);
	void longTracer(float xs, float ys, float *xt, float *yt, bool reverse);
	void genPlets();

	void restart(int i);
	void backupStart(int i);
	void updateAllPlets();
	void drawPlets(float Lwid);
	void setPtColors(unsigned char *mask, bool first);
	int traceLen;
	int currentTrace;
	int colorScheme;
	float *trLongX, *trLongY;
private:
	float *vX_S,*vY_S, *speed_S;
	int thisLayer;
//	float *vX_D, *vY_D, *speed_D;
	//float *DvX, *DvY, *DvZ, *Dspeed; // deep vectors are they used?

	//float *tracesXv, *tracesYv; // temp for the trace field

	//flowPath **flowPaths;

	// for defining boundaries and main pathways
	//float  *pathsX[MAXPATHS], *pathsY[MAXPATHS];
	//int pathLens[MAXPATHS];
	//float pathStartX[MAXPATHS], pathStartY[MAXPATHS];
	//float pathWids[MAXPATHS];
	//float pathScale[MAXPATHS];
	int nPaths_S;
//	int nPaths_D;


	//------------------------------------------------

	float *cosCor;
	float nRows, nCols;

	void genPts(int sample, int L); 
	float **Xpaths; // a set of paths for background streaklets
	float **Ypaths;


	// for circular Q pathlets
	float **Xplets;
	float **Yplets;
	int *Heads;
	bool *isActive;
	//float **sideX;
	//float **sideY;
	float *Xstart, *Ystart;  // the set of start positions		
	unsigned char *traceRed, *traceGreen, *traceBlue;
	float *newXstart, *newYstart;	// backup in reverse directinos to flow.
	int *streakletAge;  // the ages of the particles
	int *streakletColorAge;
	int *Tails;
	float *RNDcols;

	int nStreakTraces; // should be calles nStreaklets

//	glyphArray *glyphs;

	float rgbLUT[1000];
	float redBlueLUT[1000][3];

	int betweens;
	int cycle;
	int nFrames;
	float traceScale, traceScale_2;

	float *fade;

};