
#define XLAYER 1
#define YLAYER 2
#define XSIZE 1000 //900
#define YSIZE 700 
//#define YSIZE 500 
#define IN_LAYERS 80
#define MAX_TRANS 30
#define MAX_CPTS 40

#define T_LAYERS 3 // transport layers

#define EAST 1
#define NORTH 2

#define SHALLOW 1
#define DEEP 2

#define MASK 1
#define ANTS 2
#define SCHEME_1 1
#define SCHEME_2 2

#define BETWEENS 2

class cmapTex;

class pathlets;
//class transect;
class Layer;

class flowBar;

struct defPt {
	float Lat, Lon;
};


class MPAS_vis
{
public:
	MPAS_vis();
	void openFile();
	void loadNetcdfMPAS(int yr, int mo);
	void readMPAS(int time);
	void genSaveCompact(int yr);
	void loadColorMaps();
	void loadTransects_2();
	void measureTransectFlows();
	void getTransectMeans();
	void printFlowGraph(float Lat);

	void computeAMOCflows(float Lat);
	void makeTraces(int L);

	void defineLand(int L);
	void makeSpeed(float *Ev, float *Nv, float *Sv);
	void makeTransportLayer(float **inVec, float *outVec, int top, int bottom);
	void makeSmoothVec_1(float *in, float *out, int top, int bottom);
	void makeSmoothVec_2(float *inVec, float *outVec, int top, int bottom);
	void makeSmoothRadius();
	void calcCellIO();

	void saveBinaryLayers(int modelNo);

	void GulfFix();
	void loadSmoothLayers(int time); //100,200,300 etc.
	void loadLayers(int time);
	void loadLayerSet(int time, float **E, float **N, float **HE, float **HN);

	//void getDTM(float **(&dd),int &R, int &C);
	void loadDepthLayers();

	void movieSetup();
	void interplolateLayers(float p, int L);

	void readMPVectors();

	void setPoint(float rx, float ry);
	void changeTraceLen(int ch);
	void incrementTransect();

	void drawVecField(int Layer);
	void drawTransects_Paths(); // draws the lines over the plan view
	void draw_Flow_Indicators_S(int ts);
	void draw_Flow_Indicators_D(int ts);
	void drawSlice(float Lat);

	void drawPathTraces_S(int rtype, float trans);
	void drawPathTraces_D(int rtype, float trans);

	void drawCurrentTransect();
	void drawBG();
	void drawMask(unsigned char *buf);
	void drawGainLoss();

	void setSamplePts(unsigned char *buff, int L, bool first);

	float minX, maxX, minY, maxY;
	float depthScale;

	double *doubleVec_1,*doubleVec_2;
	float *fLat, *fLon, *ssh;
	int nPts;
	float **Nvec, **Evec;  // full resolutoin data
	float *Land;
	float *kernel;
	float cellcount;
	int kIndexL[100], kIndexR[100];

	float *agF_E[T_LAYERS], *agF_N[T_LAYERS];  // unsmoothed flows (transport)
	float *agF_S[T_LAYERS]; // speed

	float *smoothF_E[T_LAYERS], *smoothF_N[T_LAYERS];
	float *smoothSpeed[T_LAYERS];

//-------------------MOVIE MODE----------------------------------
	float *smoothF_E_1[T_LAYERS], *smoothF_N_1[T_LAYERS];
	float *smoothF_E_2[T_LAYERS], *smoothF_N_2[T_LAYERS];

	//float *hires_E_1
	float *hires_E_1[T_LAYERS], *hires_N_1[T_LAYERS];
	float *hires_E_2[T_LAYERS], *hires_N_2[T_LAYERS];
//-------------------------------------------------------

	float *smoothRadius;

	float *tmp; // a working vector
	float *vecdat;
	
	int ntimes;
	int npoints;  // the number of locations
	int nTransects_S, nTransects_D;
	int nDefPts_S, nDefPts_D;
	int currentTransect_S, currentTransect_D;

	int inRows, inCols; // for a distorted grid
	int streakLayer;
	int ntriangles; // the number of triangles
/////////  The stuff below here is for the vectors	
	void readEastLayer(bool save);
	void readNorthLayer(int L,bool save);
	float depthLayers[IN_LAYERS]; // the thickness of the layers
	
	float layerDepths[IN_LAYERS]; // the actual depths

	int depthBoundaries[10];

	bool drawTraces;
	bool drawTransct;
	bool drawCurrents;
	bool drawFlowSpeed;
	bool drawLongTrace;

	bool topPaths;
	bool bottomPaths;

	int viewTransect;
	int colorScheme;

	int frameNumber;

private:
	defPt flowDefPts_S[MAX_CPTS];
	defPt flowDefPts_D[MAX_CPTS];

	float FV[30]; // flow volumes

	int ncmp;

	cmapTex *colormaps;
	pathlets *flowVis_S, *flowVis_D;

	flowBar *flowIndicator[20];

	float cosLat[YSIZE];

	//Layer *layer_1; // change name later to
	//Layer *layer_2; // change name later to
	Layer *layers[T_LAYERS];

	int ncid, eastid, northid,geomid;
	int evecid, nvedid, vertvecid,vecid;  // vectors
	int std;

	int nRows, nCols,skipFactor;
	int nTimes, nLayers, xDim, yDim;
	float binsize;


	short *vdat; // the file data
	float *vecs;  // the scaled data
	float *vals;  // used for reading the frame buffer

	float cellIOdelta[35][5050]; // a two deg grid of cell IO diffs.

};