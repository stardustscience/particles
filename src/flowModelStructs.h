#ifndef ModelStructs_H
#define ModelStructs_H

#define SALT 4
#define TEMP 5
#define PARTICLE_LAYER 1
#define STREAK_LAYER 2

class dragSquare;

struct Layer  // a single time layer
{
	float **uVelocityData;
	float **vVelocityData;
	float **Wvec;
	float **Scalar;  // can be salinity or temp
};

struct vecSet  // make this a timeSeries
{		
	Layer *VecGrids;  // a set of layers  one for each depth
};

struct globalParms
{
	float *sigma;
	int nRows;
	int nCols;	
	int nTimes;
	int nSigmaLayers;

	float binsize;

	float depthScale; // vertical exaggeration factor
	int sigmaSkip;
	int streakLayer;
	int layerType;
	int nLevels;  // this is the loaded layers
	int betweens;

	int scalar;
	float **bathy;

	dragSquare *Tab[200];  // the whole set of tabs
	int tabCounter;
	// should add bathy
};

#endif
