#include "netcdf.h"
#include <iostream>
#include <stdio.h>
#include <string.h> 
#include <math.h>
#include <GL/glut.h>
//#include "flowModelStructs.h"
#include "MPAS_vis.h"
#include "cmapTex.h"
#include "pathlets.h"
//#include "transect.h"
#include "Layer.h"
#include "flowBar.h"

#define LAT_TO_METERS 111.32

#define KSIZE 41 //must be odd was 41

#define TOP_DEP 0
//#define BOT_DEP 45

//#define TOP_DEP 45
//#define BOT_DEP 80
//#define MID_DEP 45
#define N_DEPTHS 80
#define MID_DEP 50
#define BOT_DEP 70


using namespace std;

MPAS_vis::MPAS_vis()
{
	std = 0;
	nTransects_S = 0;
	nTransects_D = 0;
	viewTransect = SHALLOW;
	colorScheme = SCHEME_2;

	drawTraces = true;
	drawTransct = false;
	drawCurrents = false;
	drawFlowSpeed = false;
	drawLongTrace = false;

	currentTransect_S = 0;
	currentTransect_D = 0;
	topPaths = true;
	bottomPaths = false;
	//longTrace = 1;

	// problem if the kernel is discrete it cannot be easily stretched
	int r, c, i, index;
	float mid, rad;
	int dx, dy, rad2;

	kernel = new float[KSIZE*KSIZE];
	for (r = 0; r < KSIZE*KSIZE; ++r) kernel[r] = 1.0f;
	mid = KSIZE / 2; 
	rad2 = mid*mid;
	for(r=0;r<KSIZE;++r)
		for(c = 0;c<KSIZE;++c)
		{ 
			dx = c - mid;
			dy = r - mid;
			index = r*KSIZE + c;
			if (rad2 < (dx*dx + dy*dy))kernel[index] = 0.0f;
		}
	cellcount = 0.0f;
	for (r = 0; r < KSIZE*KSIZE; ++r) cellcount = cellcount + kernel[r];
	cerr << "CELL COUNT " << cellcount << "\n\n";
	// find left and right edge indices
	for (r = 0; r < KSIZE; ++r)
	{
		c = 0;
		do
		{
			index = r*KSIZE + c;
			++c;
		} while (kernel[index] < 1.0);
		kIndexL[r] = c-1;
		index = KSIZE - c;// because of the c increment
		kIndexR[r] = index;

		cerr << r << " " << kIndexL[r] << " " << kIndexR[r] << "\n";
	}

	for (r = 0; r < 30; ++r) FV[r] = 1.0;

	for (i = 0; i < YSIZE; ++i)
		cosLat[i] = cos((double(i) / 10.0)*3.141592 / 180.0);

	depthBoundaries[0] = TOP_DEP;
	depthBoundaries[1] = MID_DEP;
	depthBoundaries[2] = BOT_DEP;
	depthBoundaries[3] = N_DEPTHS;

    flowIndicator[0] = new flowBar(120.0f, 360.0f, 1,BETWEENS, 55.0f);
	flowIndicator[1] = new flowBar(250.0f, 590.0f, 1, BETWEENS, 30.2f);
	flowIndicator[2] = new flowBar(1020.0f, 350.0f, 1, BETWEENS,30.31f);
	flowIndicator[3] = new flowBar(1020.0f, 580.0f, 1, BETWEENS, 4.55f);
	flowIndicator[4] = new flowBar(120.0f, 25.0f, 1, BETWEENS, 22.5);// 31.8f);

	flowIndicator[5] = new flowBar(70.0f, 260.0f, 2, BETWEENS, 22.44f);
	flowIndicator[6] = new flowBar(200.0f, 480.0f, 2, BETWEENS, 17.19f);
	flowIndicator[7] = new flowBar(850.0f, 50.0f, 2, BETWEENS, 25.5f);
	flowIndicator[8] = new flowBar(1040.0f, 220.0f, 2, BETWEENS,21.6f);

	for (i = 0; i < T_LAYERS; ++i)
	{
		smoothF_E_1[i] = new float[XSIZE*YSIZE];
		smoothF_N_1[i] = new float[XSIZE*YSIZE];

		smoothF_E_2[i] = new float[XSIZE*YSIZE];
		smoothF_N_2[i] = new float[XSIZE*YSIZE];

		hires_E_1[i] = new float[XSIZE*YSIZE];
		hires_N_1[i] = new float[XSIZE*YSIZE];

		hires_E_2[i] = new float[XSIZE*YSIZE];
		hires_N_2[i] = new float[XSIZE*YSIZE];
	}

	Nvec = new float*[80];
	Evec = new float*[80];
	for (i = 0; i < 80; ++i)
	{
		Evec[i] = new float[XSIZE*YSIZE];
		Nvec[i] = new float[XSIZE*YSIZE];
	}


}


void MPAS_vis::genSaveCompact(int yr)
{
	int  mo,yy;
	for(yy = yr; yy < (yr+5);++yy)
	for (mo = 1; mo <= 12; ++mo)
	{

		cerr << "YR  MO " << yy << " " << mo << "\n";
		loadNetcdfMPAS(yy, mo); //load the full 80 layers
		saveBinaryLayers(yy * 12 + mo);
	}

}
void MPAS_vis::loadNetcdfMPAS(int yr,int mo)
{

	int vecid, nattrib, ndim, i, L;
	int vecNorthId, vecEastId;
	long len;
	int dims[5], dimsz[5];
	long start[5], size[5];
	double max, min, v;

	char filename[120];

	nc_type ttype;
	doubleVec_1 = new double[XSIZE*YSIZE];

	//yr = 55; mo = 3;

	if(mo < 10)
		sprintf(filename, "MPAS_Files/UV.mpaso.hist.am.timeSeriesStatsMonthly.00%d-0%d-01.interp0.1x0.1_NorthAtl.nc", yr, mo);
		//sprintf(filename, "MPAS_Files/UV.mpasoMonthly.00%d-0%d-01.interp0.1x0.1_NorthAtl.nc", yr, mo);
	else sprintf(filename, "MPAS_Files/UV.mpaso.hist.am.timeSeriesStatsMonthly.00%d-%d-01.interp0.1x0.1_NorthAtl.nc", yr, mo);
		//sprintf(filename, "MPAS_Files/UV.mpasoMonthly.00%d-%d-01.interp0.1x0.1_NorthAtl.nc", yr, mo);

	ncmp = ncopen(filename, NC_NOWRITE); frameNumber = mo*100;
	
	//ncmp = ncopen("UV.mpasoMonthly.0055-01-01.interp0.1x0.1_NorthAtl.nc", NC_NOWRITE); frameNumber = 100;
	//ncmp = ncopen("UV.mpasoMonthly.0055-02-01.interp0.1x0.1_NorthAtl.nc", NC_NOWRITE); frameNumber = 200;
	//ncmp = ncopen("UV.mpasoMonthly.0055-03-01.interp0.1x0.1_NorthAtl.nc", NC_NOWRITE);frameNumber = 300;
	//ncmp = ncopen("UV.mpasoMonthly.0055-04-01.interp0.1x0.1_NorthAtl.nc", NC_NOWRITE); frameNumber = 400;
	//ncmp = ncopen("UV.mpasoMonthly.0055-05-01.interp0.1x0.1_NorthAtl.nc", NC_NOWRITE); frameNumber = 500;
	//ncmp = ncopen("UV.mpasoMonthly.0055-06-01.interp0.1x0.1_NorthAtl.nc", NC_NOWRITE);frameNumber = 600;
	//ncmp = ncopen("UV.mpasoMonthly.0055-07-01.interp0.1x0.1_NorthAtl.nc", NC_NOWRITE); frameNumber = 700;
	//ncmp = ncopen("UV.mpasoMonthly.0055-08-01.interp0.1x0.1_NorthAtl.nc", NC_NOWRITE); frameNumber = 800;
	//ncmp = ncopen("UV.mpasoMonthly.0055-09-01.interp0.1x0.1_NorthAtl.nc", NC_NOWRITE);frameNumber = 900;
	//ncmp = ncopen("UV.mpasoMonthly.0055-10-01.interp0.1x0.1_NorthAtl.nc", NC_NOWRITE);frameNumber = 1000;
	//ncmp = ncopen("UV.mpasoMonthly.0055-11-01.interp0.1x0.1_NorthAtl.nc", NC_NOWRITE); frameNumber = 1100;
	//ncmp = ncopen("UV.mpasoMonthly.0055-12-01.interp0.1x0.1_NorthAtl.nc", NC_NOWRITE); frameNumber = 1200;

	vecEastId = ncvarid(ncmp, "timeMonthly_avg_velocityZonal");
	vecNorthId = ncvarid(ncmp, "timeMonthly_avg_velocityMeridional");

	ncvarinq(ncmp, vecNorthId, 0, &ttype, &ndim, dims, &nattrib);

	cerr << "READ MPAS -----------------------------------\n";
	//cerr << "Dims " << ndim << ":  ";
	//cerr << "Type " << ttype << "\n";

	for (i = 0; i < ndim; ++i)
	{
		ncdiminq(ncmp, dims[i], 0, &len);
		size[i] = len;
		//cerr << i << " " << len << "\n ";
	}

	start[0] = 0;
	start[1] = 10;
	start[2] = 0;
	start[3] = 0;

	size[0] = 1; // get a single layer
	size[1] = 1;
	size[2] = YSIZE;
	size[3] = XSIZE;

	cerr << "\n READ LAYERS ";

	for (L = 0; L < IN_LAYERS; ++L)
	{
		start[1] = L;

		//cerr << "READ LAYER " << L << "\n";
		ncvarget(ncmp, vecNorthId, start, size, doubleVec_1);
		for (i = 0; i < YSIZE*XSIZE; ++i)
		{
			v = doubleVec_1[i];
			if (v < -100.0) v = 0.0;// -2.0;
			Nvec[L][i] = float(v);
		}

		ncvarget(ncmp, vecEastId, start, size, doubleVec_1);
		min = 10000.0;
		max = -1000.0;
		for (i = 0; i < YSIZE*XSIZE; ++i)
		{
			v = doubleVec_1[i];
			if (v < -100.0) v = 0.0;
			Evec[L][i] = float(v);
			if (v < min) min = doubleVec_1[i];
			if (v > max) max = doubleVec_1[i];
		}
		cerr << L << " ";
		if (L % 10 == 9) cerr << "\n";
	}
	cerr << "\n DATA LOADED\n";


	//defineLand(0);  // use layer 0 to define Land
					//cerr << "Speed made \n";

	for (i = 0; i < T_LAYERS; ++i)
	{
		makeTransportLayer(Evec, agF_E[i], depthBoundaries[i], depthBoundaries[i + 1]); // sum East and North VFs
		makeTransportLayer(Nvec, agF_N[i], depthBoundaries[i], depthBoundaries[i + 1]); // sum East and North VFs
		makeSpeed(agF_E[i], agF_N[i], agF_S[i]);
		//cerr << "Avg HRes Layer " << i << " \n\n";
		//} //SMOOTH_1or2
		makeSmoothVec_2(agF_E[i], smoothF_E[i], depthBoundaries[i], depthBoundaries[i + 1]); // Smoothed East and North VFs
		makeSmoothVec_2(agF_N[i], smoothF_N[i], depthBoundaries[i], depthBoundaries[i + 1]); // Smoothed East and North VFs
		//makeSmoothVec_1(agF_E[i], smoothF_E[i], depthBoundaries[i], depthBoundaries[i + 1]); // Smoothed East and North VFs
		//makeSmoothVec_1(agF_N[i], smoothF_N[i], depthBoundaries[i], depthBoundaries[i + 1]); // Smoothed East and North VFs
		makeSpeed(smoothF_E[i], smoothF_N[i], smoothSpeed[i]);
		//cerr << "Smooth Speed Layer " << i << " \n\n";
	}


}





void MPAS_vis::readMPAS(int time)
{
	int vecid, nattrib, ndim, i, L;
	int vecNorthId, vecEastId;
	long len;
	int dims[5], dimsz[5];
	long start[5], size[5];
	double max, min, v;

	nc_type ttype;
	doubleVec_1 = new double[XSIZE*YSIZE];

	

	for (i = 0; i < T_LAYERS; ++i)
	{
		agF_E[i] = new float[XSIZE*YSIZE];  // these are the unsmoothed depth layers
		agF_N[i] = new float[XSIZE*YSIZE];
		agF_S[i] = new float[XSIZE*YSIZE];
		smoothSpeed[i] = new float[XSIZE*YSIZE];
		smoothF_E[i] = new float[XSIZE*YSIZE];
		smoothF_N[i] = new float[XSIZE*YSIZE];
	}
	Land = new float[XSIZE*YSIZE];
	tmp = new float[XSIZE*YSIZE];

	frameNumber = 0;

	//

	loadDepthLayers();  // the distances between layers
	loadLayers(time);
	defineLand(0);
	makeSmoothRadius();	
	//loadNetcdfMPAS(55, 3); //load the full 80 layers
	//loadSmoothLayers(frameNumber + 100);

	//loadLayers(time);


	//loadSmoothLayers(frameNumber);
	//GulfFix();
	//for (i = 0; i < T_LAYERS; ++i)
	//	makeSpeed(smoothF_E[i], smoothF_N[i], smoothSpeed[i]);


	layers[0] = new Layer(YSIZE, XSIZE,0);
	layers[0]->setFlows(smoothF_E[0], smoothF_N[0], smoothSpeed[0], Nvec, Evec);
	layers[0]->setVecs(agF_E[0], agF_N[0], agF_S[0]);

	layers[1] = new Layer(YSIZE, XSIZE,1);
	layers[1]->setFlows(smoothF_E[1], smoothF_N[1], smoothSpeed[1], Nvec, Evec);
	layers[1]->setVecs(agF_E[1], agF_N[1], agF_S[1]);

	flowVis_S = new pathlets(YSIZE, XSIZE,0, colorScheme);
	flowVis_S->setVecs(smoothF_E[0], smoothF_N[0], smoothSpeed[0]);

	flowVis_S->genPlets(); //WORK

	flowVis_D = new pathlets(YSIZE, XSIZE, 1, colorScheme);
	flowVis_D->setVecs(smoothF_E[1], smoothF_N[1], smoothSpeed[1]);
	flowVis_D->genPlets(); //WORK


	flowVis_S->makeTraces();
	flowVis_D->makeTraces();

	//flowVis_S->loadPaths("TraceStarts.txt");

	cerr << "Traces made \n";
	flowVis_S->longTrace(500.0f, 500.0f , false);

	delete[] doubleVec_1;

}

void MPAS_vis::saveBinaryLayers(int modelNo)
{
	int L;
	char fnm[40];
	FILE *ofp;

	// the smoothed layers
	for (L = 0; L < T_LAYERS; ++L)  ///SMOOTH_1or2
	{
		sprintf(fnm, "BinLayers_2/TRL_S_%d_%d_%d.bin", modelNo, L, depthBoundaries[L + 1]);  // TRL = transport layer
		cerr << "SAVE BINARY " << fnm << "\n";
		ofp = fopen(fnm, "wb");
		//fwrite(agF_E[L], 4, XSIZE*YSIZE, ofp); // unsmoothed
		//fwrite(agF_N[L], 4, XSIZE*YSIZE, ofp);

		fwrite(smoothF_E[L], 4, XSIZE*YSIZE, ofp); // smoothed
		fwrite(smoothF_N[L], 4, XSIZE*YSIZE, ofp);
		fclose(ofp);
	}	
	//SMOOTH_1or2
	for (L = 0; L < T_LAYERS; ++L)
	{
		sprintf(fnm, "BinLayers_2/TRL_U_%d_%d_%d.bin", modelNo, L, depthBoundaries[L + 1]);  // TRL = transport layer
		cerr << "SAVE BINARY " << fnm << "\n";
		ofp = fopen(fnm, "wb");
		fwrite(agF_E[L], 4, XSIZE*YSIZE, ofp); // unsmoothed
		fwrite(agF_N[L], 4, XSIZE*YSIZE, ofp);
		fclose(ofp);
	}
}

void MPAS_vis::GulfFix() // create a bypass for gulf of mexico
{

	int r, c, jj;
	int rad2;

	float rad,dx,dy;
	float wid;
	int iwid;

	wid = 55.0f; // 5.5 deg height

	float vx;
	int cx, cy;
	cx = 25; cy = 215;
	//HACK
	for (r = 20; r < (YSIZE); ++r)
		for (c = 0; c < 30; ++c)
		{
			jj = r*XSIZE + c;

			dy = (r - cy);
			dx = (c - cx)*3; // 1/3 width

			vx = float(dy)*0.1 / wid;

			rad2 = dx*dx + dy*dy;

			if (rad2 < 3000) { // 5.5 deg half ht 3000

				smoothF_E[0][jj] = vx;
				smoothF_N[0][jj] = 0.10;
			}

		}

}


void MPAS_vis::loadLayers(int mn)
{
	int L;
	//loadLayerSet(time, smoothF_E, smoothF_N, agF_E, agF_N);
	loadLayerSet(mn, smoothF_E_1, smoothF_N_1, hires_E_1, hires_N_1);
	loadLayerSet(mn+1, smoothF_E_2, smoothF_N_2, hires_E_2, hires_N_2);
	for (L = 0; L < T_LAYERS; ++L)
	{
		//makeSpeed(smoothF_E[L], smoothF_N[L], smoothSpeed[L]);		
		interplolateLayers(0.0f, L);
		//makeSpeed(agF_E[L], agF_N[L], agF_S[L]);
	}

}

void MPAS_vis::loadLayerSet(int time, float **SE, float **SN, float **HE, float **HN)
{
	int L;
	char fnm[40];
	FILE *ifp;
	for (L = 0; L < T_LAYERS; ++L)
	{
		//sprintf(fnm, "Backup_SmoothLayers/TRL_S_%d_%d_%d.bin", time, L, depthBoundaries[L + 1]);  // TRL = transport layer
		sprintf(fnm, "BinLayers_1/TRL_S_%d_%d_%d.bin", time, L, depthBoundaries[L + 1]);  // TRL = transport layer
		cerr << "LOAD BINARY " << fnm << "\n";
		cerr << "LOAD BINARY " << fnm << "\n";
		ifp = fopen(fnm, "rb");
		fread(SE[L], 4, XSIZE*YSIZE, ifp); // smoothed
		fread(SN[L], 4, XSIZE*YSIZE, ifp);
		fclose(ifp);

	}
	for (L = 0; L < T_LAYERS; ++L)
	{
		sprintf(fnm, "BinLayers_1/TRL_U_%d_%d_%d.bin", time, L, depthBoundaries[L + 1]);  // TRL = transport layer
		cerr << "LOAD BINARY " << fnm << "\n";
		ifp = fopen(fnm, "rb");
		fread(HE[L], 4, XSIZE*YSIZE, ifp); // unsmoothed
		fread(HN[L], 4, XSIZE*YSIZE, ifp);
		fclose(ifp);
	}
}



void MPAS_vis::loadSmoothLayers(int time)
{
	int L;
	char fnm[40];
	FILE *ifp;
	for (L = 0; L < T_LAYERS; ++L)
	{
		sprintf(fnm, "BinLayers_1/TRL_S_%d_%d_%d.bin", time, L, depthBoundaries[L + 1]);  // TRL = transport layer
		cerr << "LOAD BINARY " << fnm << "\n";
		ifp = fopen(fnm, "rb");
		fread(smoothF_E[L], 4, XSIZE*YSIZE, ifp); // unsmoothed
		fread(smoothF_N[L], 4, XSIZE*YSIZE, ifp);
		fclose(ifp);
	}

}

void MPAS_vis::movieSetup()
{
	int i,j;

	cerr << "MOVIE SETUP \n";  // get two frames to interpolate.
	/*
	for (i = 0; i < T_LAYERS; ++i)
	{
		smoothF_E_1[i] = new float[XSIZE*YSIZE];
		smoothF_N_1[i] = new float[XSIZE*YSIZE];

		smoothF_E_2[i] = new float[XSIZE*YSIZE];
		smoothF_N_2[i] = new float[XSIZE*YSIZE];
	}
*/
	for (i = 0; i < T_LAYERS; ++i) // copy
		for (j = 0; j < XSIZE*YSIZE; ++j)
		{
			smoothF_E_1[i][j] = smoothF_E[i][j];
			smoothF_N_1[i][j] = smoothF_N[i][j];
		}

	loadSmoothLayers(frameNumber + 100);

	for (i = 0; i < T_LAYERS; ++i) // copy
		for (j = 0; j <XSIZE*YSIZE; ++j)
		{
			smoothF_E_2[i][j] = smoothF_E[i][j];
			smoothF_N_2[i][j] = smoothF_N[i][j];
		}
/*	
	for (i = 2000; i < 5000; i = i + 730)
	{
		cerr << "TEST " << i << " " << smoothF_E_1[0][i] << " " << smoothF_E_2[0][i] << "\n";
	}
*/
}

void MPAS_vis::interplolateLayers(float p,int L) 
{
	int j;
	//cerr << "MAKE frame " << p << "\n";

	for (j = 0; j <XSIZE*YSIZE; ++j)
	{

		smoothF_E[L][j] = smoothF_E_1[L][j] * (1.0 - p) + smoothF_E_2[L][j] * p;
		smoothF_N[L][j] = smoothF_N_1[L][j] * (1.0 - p) + smoothF_N_2[L][j] * p;
	}
	makeSpeed(smoothF_E[L], smoothF_N[L], smoothSpeed[L]);

	for (j = 0; j <XSIZE*YSIZE; ++j)
	{

		agF_E[L][j] = hires_E_1[L][j] * (1.0 - p) + hires_E_2[L][j] * p;
		agF_N[L][j] = hires_N_1[L][j] * (1.0 - p) + hires_N_2[L][j] * p;
	}
	makeSpeed(agF_E[L], agF_N[L], agF_S[L]);
	//flowVis_S->setVecs(smoothF_E[L], smoothF_N[L], smoothSpeed[L]);
	//makeTraces(1);

}

void MPAS_vis::setPoint(float rx, float ry)
{
	float Lat, Lon;
	float px, py; // in data space
	Lat = (698.0f - ry) / 10.0f;
	Lon = ((rx - 500.0f) / cosLat[int(700.0f -ry)] + 500.0f)/10.0f - 90.0f;

	px = (Lon + 90.0f) *10.0f;
	py = Lat * 10.0f;

	cerr << "LON LAT  x y " << Lon << " " << Lat << " " << rx << " " << 700.0f - ry << "\n";

	flowVis_S->longTrace(px, py, true);

}

void MPAS_vis::makeTraces(int L)
{

	// Smoothed vs not Smoothed
	flowVis_S->setVecs(smoothF_E[L], smoothF_N[L], smoothSpeed[L]);
	//flowVis_S->setVecs(agF_E[L], agF_N[L], agF_S[L]);
	flowVis_S->makeTraces();

	/*
	if( TB == 1){
		flowVis_S->setVecs(smoothF_E[0], smoothF_N[0], smoothSpeed[0]);
		flowVis_S->makeTraces();}
	if (TB == 2) {
		flowVis_S->setVecs(smoothF_E[1], smoothF_N[1], smoothSpeed[1]);
		flowVis_S->makeTraces();
	}
	if (TB == 3) {
		flowVis_S->setVecs(smoothF_E[2], smoothF_N[2], smoothSpeed[2]);
		flowVis_S->makeTraces();
	}
*/
}

void MPAS_vis::changeTraceLen(int ch)
{
	float elon, elat;
	flowVis_S->traceLen = flowVis_S->traceLen + ch;
	if (flowVis_S->traceLen >= LONG_TRACE)flowVis_S->traceLen = LONG_TRACE;
	cerr << "TRACE LEN " << flowVis_S->traceLen << "\n";
	cerr << "TRACE END " << flowVis_S->trLongX[flowVis_S->traceLen] << " " << flowVis_S->trLongY[flowVis_S->traceLen] * 0.1 << "\n";
	elat = flowVis_S->trLongY[flowVis_S->traceLen];
	elon = ((flowVis_S->trLongX[flowVis_S->traceLen] - 500.0f) / cosLat[int(elat)] + 500.0f) / 10.0f - 90.0f;

	cerr << "TRACE END Lat Lon " << elat << " " << elon << "\n";
}


void MPAS_vis::incrementTransect()
{
	float s, d;
	if (viewTransect == SHALLOW)
	{
		currentTransect_S = (currentTransect_S + 1) % nTransects_S;
		cerr << "CURRENT SHALLOW " << currentTransect_S << "\n";
		layers[0]->getTransectFlow(currentTransect_S, s, d);
	}
	if (viewTransect == DEEP)
	{
		currentTransect_D = (currentTransect_D + 1) % nTransects_D;
		cerr << "CURRENT DEEP " << currentTransect_D << "\n";
		layers[1]->getTransectFlow(currentTransect_D, s, d);
	}
	//cerr << "\nShallowFlow " << transects[currentTransect]->shallowFlow << "\n";
	//cerr << "DeepFlow " << transects[currentTransect]->deepFlow << "\n";
	//cerr << "FLOW VOLUME " << currentTransect << " " << transectFlow[currentTransect] << "\n";
}

void MPAS_vis::loadDepthLayers()
{
	int i, ii;
	float d, ad, th;
	FILE *fp;
	float dL[80];
	fp = fopen("Depths.txt", "r");

	for (i = 0; i < IN_LAYERS; ++i)
	{
		fscanf(fp, "%d %f %f", &ii, &d, &ad);
		depthLayers[i] = d;
		layerDepths[i] = ad;
	}

	// recald the detph layers considering samples in center of layer.
	dL[0] = depthLayers[0] + (depthLayers[0] + depthLayers[1]) / 2.0f;
	for (i = 1; i < IN_LAYERS-1; ++i)
	{
		dL[i] = (depthLayers[i] + depthLayers[i + 1]) / 2.0f;
		//cerr << "LAYER " << i << " " << depthLayers[i] << "\n";
	}
	dL[IN_LAYERS - 1] = depthLayers[IN_LAYERS - 1];  // leave the same thickness

	for (i = 0; i < IN_LAYERS; ++i)
		depthLayers[i] = dL[i];

}

void MPAS_vis::loadColorMaps()
{
	colormaps = new cmapTex();
	colormaps->load("Blue_Red.txt", false, 0);
	colormaps->load("Blue_Red.txt", false, 1);
	colormaps->load("Spiral.txt", false, 2);
	colormaps->load("Spiral.txt", false, 3);
	colormaps->mkTex2Dcmap();
}



void MPAS_vis::loadTransects_2()
{
	int i, p1, p2;
	FILE *fp;
	float lon1, lat1, lon2, lat2;
	int tf, d;
	bool yes;
	char Code[8];
	int cpLat, cpLon;


	cerr << "LOAD TRANSECTS 2\n";

	fp = fopen("FlowDef_S.txt", "r");
	fscanf(fp, "%s %d", Code, &nDefPts_S);
	for (i = 0; i < nDefPts_S; ++i)
	{
		fscanf(fp, "%d %f %f", &d, &lon1, &lat1);
		cerr << d << " " << lon1 << " " << lat1 << "\n";
		flowDefPts_S[i].Lat = lat1;
		flowDefPts_S[i].Lon = lon1;
	}
	fscanf(fp, "%s %d", Code, &nTransects_S);
	cerr << "LOAD N " << nTransects_S << "\n";

	/*
	layers[0]->findCriticalPoint(int(flowDefPts_S[0].Lat)*10, int(flowDefPts_S[0].Lon) * 10, cpLat, cpLon, true);
	flowDefPts_S[0].Lat = float(cpLat)/10.0f; flowDefPts_S[0].Lon = -90.0f + float(cpLon)/10.0f;
	cerr << "\n";
	layers[0]->findCriticalPoint(int(flowDefPts_S[1].Lat) * 10, int(flowDefPts_S[1].Lon) * 10,cpLat, cpLon, false);
	flowDefPts_S[1].Lat = float(cpLat) / 10.0f; flowDefPts_S[1].Lon = -90.0f + float(cpLon) / 10.0f;
	*/
	cerr << "\nSHALLOW TRANSECTS \n";
	for (i = 0; i < nTransects_S; ++i)
	{
		fscanf(fp, "%d %d %d ", &d, &p1, &p2);
		lon1 = flowDefPts_S[p1].Lon; lat1 = flowDefPts_S[p1].Lat;
		lon2 = flowDefPts_S[p2].Lon; lat2 = flowDefPts_S[p2].Lat;
		cerr << "TRANS " << i << " " << lon1 << " " << lat1 << "    " << lon2 << " " << lat2 << "   ";
		layers[0]->setTransect_S(lon1, lat1, lon2, lat2, i);
		//transectFlow[i] = layer->computeTransectHR(lon1, lat1, lon2, lat2, depthLayers, i, MID_DEP, SHALLOW);
		cerr << layers[0]->computeTransectHR_2(i) << "\n";
	}
	fclose(fp);

	float maxp, mlat, mlon;

	layers[0]->computePeakTransectHR(11, maxp,mlon,mlat);  //adjusts the peak along a line

	fp = fopen("FlowDef_D.txt", "r");
	fscanf(fp, "%s %d", Code, &nDefPts_D);
	for (i = 0; i < nDefPts_D; ++i)
	{
		fscanf(fp, "%d %f %f", &d, &lon1, &lat1);
		cerr << d << " " << lon1 << " " << lat1 << "\n";
		flowDefPts_D[i].Lat = lat1;
		flowDefPts_D[i].Lon = lon1;
	}

	cerr << "\nDEEP TRANSECTS \n";
	fscanf(fp, "%s %d", Code, &nTransects_D);
	for (i = 0; i < nTransects_D; ++i)
	{
		fscanf(fp, "%d %d %d ", &d, &p1, &p2);
		lon1 = flowDefPts_D[p1].Lon; lat1 = flowDefPts_D[p1].Lat; // generalize this
		lon2 = flowDefPts_D[p2].Lon; lat2 = flowDefPts_D[p2].Lat;
		layers[1]->setTransect_S(lon1, lat1, lon2, lat2, i);
		cerr << "TRANS DEEP " << i << " " << lon1 << " " << lat1 << "    " << lon2 << " " << lat2 << "\n";
		//if (i == 5)
		layers[1]->computeTransectHR_2(i);
		//	transectFlow[i] = layer->computeTransectHR(lon1, lat1, lon2, lat2, depthLayers, i, 40); //HACK
		//else
		//	transectFlow[i] = layer->computeTransectHR(lon1, lat1, lon2, lat2, depthLayers, i, MID_DEP);
	}
	fclose(fp);
	measureTransectFlows();
}

void MPAS_vis::getTransectMeans()
{
	int i;
	cerr << "TRANSECT MEANS \n ";
	for (i = 0; i < 9; ++i)
	{
		cerr << i << " " << flowIndicator[i]->getMean() << "\n";
	}
}

void MPAS_vis::measureTransectFlows()
{
	int i;

	float sh, dp;
	float x1, x2, y1, y2;
	for (i = 0; i < nTransects_S; ++i)layers[0]->computeTransectHR_2(i);
	for (i = 0; i < nTransects_D; ++i)layers[1]->computeTransectHR_2(i);

	float maxp, mlat, mlon;

	//layers[0]->computePeakTransectHR(11, maxp, mlon, mlat); to adjust the position



	// get the flow volumes  probably get rid of this
	layers[0]->getTransectFlow(0, sh, dp);	FV[1] = -sh; // Max Gulf Stream
	layers[0]->getTransectFlow(1, sh, dp);	FV[2] = sh;  // total Gulf Stream Return
	layers[0]->getTransectFlow(9, sh, dp);	FV[3] = sh;  // Gulf stream mid return
	layers[0]->getTransectFlow(2, sh, dp);	FV[4] = sh;
	layers[0]->getTransectFlow(6, sh, dp);	FV[5] = -sh;
	layers[0]->getTransectFlow(8, sh, dp);	FV[6] = -sh; // south return loop
	layers[0]->getTransectFlow(5, sh, dp);	FV[7] = -sh; // to UK

	layers[1]->getTransectFlow(0, sh, dp);	FV[8] = -sh; // deep return
	layers[1]->getTransectFlow(8, sh, dp);	FV[9] = sh; // deep return north loop
	layers[1]->getTransectFlow(7, sh, dp);	FV[10] = sh; // deep return
	layers[1]->getTransectFlow(9, sh, dp);	FV[11] = sh; // deep return

	/*
	FV[10] = FV[2] - FV[3]; // outer loop return
	FV[11] = FV[5] - FV[7]; // gulf stream south
	FV[12] = FV[1] - FV[2] - FV[7];  // to upper gyre
*/

	cerr << "FLOW VOLUMES \n";
	for (i = 1; i < 12; ++i)
			cerr << i << " " << FV[i] << "\n";

	// compute shallow flows
	/*
	// set the illustrating paths
	layers[0]->computePathline(5, 0, MID , true, FV[7]); // gulf stream to UK
	layers[0]->computePathLoop(1, 0, LAST, FV[10]); // Atlantic outer loop
	layers[0]->computePathLoop(1, 0, MID,FV[3]); // Atlantic middle loop
	layers[0]->computePathline(0, 6, MID, true, FV[11]); //Gulf steam south
	layers[0]->computePathline(2, 2, MID, true, FV[4]); // ATLANTIC UPPER GYRE
	layers[0]->computePath2part(8, 0, 6, MID, FV[6]); // LOWER ATLANTIC LOOP
	layers[0]->computePathline(3, 0, FIRST, true, FV[12]); //Gulf Stream to upper GYRE
//	cerr << "ATLANTIC PATHS STARTS \n";

	layers[1]->computePathline(7, 0, MID, true, 20.0f); // deep western current to south
	layers[1]->computePathline(4, 0, MID, false, 20.0f); // gulf stream to UK
	layers[1]->computePathline(4, 4, MID,true, 10.0f); // Upper atlantic deep loop

	*/

//	float pLon = -90.0f + layers[0]->startX[4] / 10.0f;
//	float pLat = layers[0]->startY[4] / 10.0f;

//	cerr << "XXXXXXXXXXXXXXXXXX " << flowDefPts_S[0].Lon << " " << pLon << "\n";
//	layers[0]->computeTransectHR(flowDefPts_S[0].Lon, flowDefPts_S[0].Lat, pLon, pLat, depthLayers, 9, MID_DEP,SHALLOW);

	layers[0]->getTransectAnchors(0, x1, y1, x2, y2);
	flowIndicator[0]->setAnchor(x1, y1);
	layers[0]->getTransectAnchors(2, x1, y1, x2, y2);
	flowIndicator[1]->setAnchor(x1, y1);
	layers[0]->getTransectAnchors(1, x1, y1, x2, y2);
	flowIndicator[2]->setAnchor(x2, y2);
	layers[0]->getTransectAnchors(4, x1, y1, x2, y2);
	flowIndicator[3]->setAnchor(x2, y2);
	layers[0]->getTransectAnchors(5, x1, y1, x2, y2);
	flowIndicator[4]->setAnchor(x1, y1);

	layers[1]->getTransectAnchors(0, x1, y1, x2, y2);
	flowIndicator[5]->setAnchor(x1, y1);
	layers[1]->getTransectAnchors(8, x1, y1, x2, y2);
	flowIndicator[6]->setAnchor(x1, y1);
	layers[1]->getTransectAnchors(7, x1, y1, x2, y2);
	flowIndicator[7]->setAnchor(x2, y2);
	layers[1]->getTransectAnchors(9, x1, y1, x2, y2);
	flowIndicator[8]->setAnchor(x2, y2);


}



void MPAS_vis::defineLand(int L)
{
	int i;
	double e, n, v;

	for (i = 0; i < YSIZE*XSIZE; ++i)
	{
		v = agF_S[0][i];
		if (v > 0)Land[i] = 1.0f;
		else Land[i] = 0.0f;
	}

	// wipe out small islands

	int r, c, ii, jj;
	for (r = 0; r <694; ++r) 
		for (c = 0; c < XSIZE - 6; ++c) //HA
		{
			ii = r*XSIZE + c;
			jj = ii + 6;
			if (Land[ii] > 0.5 && Land[jj] > 0.5)
			{
				Land[ii + 1] = 1.0;
				Land[ii + 2] = 1.0;
				Land[ii + 3] = 1.0;
			}
		}

	for (r = 0; r <694; ++r)
		for (c = 0; c < XSIZE - 6; ++c) //HA
		{
			ii = r*XSIZE + c;
			jj = (r+6)*XSIZE + c;
			if (Land[ii] > 0.5 && Land[jj] > 0.5)
			{
				Land[ii + XSIZE] = 1.0;
				Land[ii + 2*XSIZE] = 1.0;
				Land[ii + 3*XSIZE] = 1.0;
			}
		}





}


void MPAS_vis::makeTransportLayer(float **inVec, float *outVec, int top, int bottom)
{
	int i,L;
	double sum;
	for (i = 0; i < XSIZE*YSIZE; ++i) outVec[i] = 0.0f;

	for (i = 0; i < XSIZE*YSIZE; ++i) {
		sum = 0;
		for (L = top; L < bottom; ++L)
			sum += inVec[L][i] * depthLayers[L];
		tmp[i] = float(sum) / 1000.0;  // result in units of 1K
	}

	for (i = 0; i < XSIZE*YSIZE; ++i) outVec[i] = tmp[i];

}



void MPAS_vis::makeSmoothRadius()
{
	int i, L, ii;
	int r, c, rr,cc, rrr,ccc;
	int knwid, knhalf, kXhal;
	knwid = KSIZE;
	knhalf = knwid / 2;
	float minr2,rad2;

	cerr << "SMOOTH RADIUS START\n";

	smoothRadius = new float[700 * 1000];

	for (r = 0; r < 700 * 1000; ++r)smoothRadius[r] = 5.0f;

	for (r = 6; r <694; ++r) {
		for (c =6; c < XSIZE -6; ++c) //HACK should be knhalf
		{		
			minr2 = knhalf*knhalf;
			for (rr = -knhalf; rr <= knhalf; ++rr)
			{
		
				for (cc = -knhalf; cc <= knhalf; ++cc)
				{
					rad2 = rr*rr + cc*cc;

					rrr = r + rr; ccc = c + cc;
					ii = rrr*XSIZE + c + cc;
					if(rrr > 5 && rrr < 695)
						if(ccc > 5 && ccc < 695)
							if (Land[ii] < 0.1) {
								if (rad2 < minr2)minr2 = rad2;
							}
				}
			}
			ii = r*XSIZE + c;
			smoothRadius[ii] = sqrt(minr2);
		}
	}

	cerr << "SMOOTH RADIUS DONE\n";

}


// apply a cirular box filter using the land distance metric 

void MPAS_vis::makeSmoothVec_1(float *inVec, float *outVec, int top, int bottom)
{
	int i, L, ii, jj;
	int r, c, kr, kc;
	int knwid, knhalf, kXhalf;
	knwid = KSIZE;
	knhalf = knwid / 2;
	//double avgfac;
	double sum, v, wt;
	int cellcount;

	for (i = 0; i < XSIZE*YSIZE; ++i) outVec[i] = inVec[i];

	int rad2;
	int radSQ;

	//for (r = knhalf; r < (YSIZE - knhalf); ++r) {//HACK should be knhalf
	for (r = 10; r <680; ++r) {//HACK should be knhalf
		//for (c = kXhalf; c < XSIZE-kXhalf; ++c) //HACK should be knhalf
		for (c = 10; c < XSIZE - 10; ++c) //HACK should be knhalf
		{
			ii = r*XSIZE + c;
			knhalf = smoothRadius[ii];
			kXhalf = 2 * (int(knhalf / cosLat[r]) / 2) + 1; // make it odd	
			radSQ = knhalf*kXhalf;

			sum = 0;
			cellcount = 0;
			//for (kr = r - knhalf; kr <= r + knhalf; ++kr)  
			//	for (kc = c - kXhalf; kc <= c + kXhalf; ++kc)
			wt = 0.0001;
			for (kr = r - knhalf; kr <= r + knhalf; kr += 2) // fast smooth
				for (kc = c - kXhalf; kc <= c + kXhalf; kc += 2)
				{
					int rad2 = (kr - r)*(kr - r) + (kc - c)*(kc - c);
					if (rad2 <= radSQ)
					{
						ii = kr*XSIZE + kc;
						//sum = sum + inVec[ii]*Land[ii]; 
						if (ii > 0) {
							sum = sum + inVec[ii];// *Land[ii];
												  //wt += Land[ii];
							++cellcount;
						}
					}
				}
			jj = r*XSIZE + c;
			outVec[jj] = float(sum / cellcount);
			//outVec[jj] = float(sum / wt);
		}
	}

}



// smooth with land correction

void MPAS_vis::makeSmoothVec_2(float *inVec, float *outVec, int top, int bottom)
{
	int i, L, ii, jj;
	int r, c, kr, kc;
	int knwid, knhalf,kXhalf;
	knwid = KSIZE;
	knhalf = knwid / 2;
	//double avgfac;
	double sum, v,wt;
	int cellcount, landcount;

	for (i = 0; i < XSIZE*YSIZE; ++i) outVec[i] = inVec[i];

	int rad2;
	int radSQ;

	//for (r = knhalf; r < (YSIZE - knhalf); ++r) {//HACK should be knhalf
	for (r = knhalf; r <650; ++r) {//HACK should be knhalf
		kXhalf = 2 * (int(knhalf / cosLat[r]) / 2) + 1; // make it odd	

		radSQ = knhalf*kXhalf;

		for (c = kXhalf; c < XSIZE-kXhalf; ++c) //HACK should be knhalf
		//for (c = 100; c < XSIZE - 200; ++c) //HACK should be knhalf
		{
			sum = 0;
			cellcount = 0;
			landcount = 0;
			//for (kr = r - knhalf; kr <= r + knhalf; ++kr)  
			//	for (kc = c - kXhalf; kc <= c + kXhalf; ++kc)
			wt = 0.0001;
			for (kr = r - knhalf; kr <= r + knhalf; kr +=2) // fast smooth
				for (kc = c - kXhalf; kc <= c + kXhalf; kc+=2)
				{
					int rad2 = (kr-r)*(kr-r) + (kc-c)*(kc-c);
					if(rad2 <= radSQ)
					{
						ii = kr*XSIZE + kc;
						//sum = sum + inVec[ii]*Land[ii]; 
						sum = sum + inVec[ii];// *Land[ii];
						if(Land[ii] < 0.5)++landcount;
						++cellcount;
					}
				}
			jj = r*XSIZE + c;
			// correct for land -- add back values distributed over the land
			outVec[jj] = Land[jj] * (float(sum / cellcount) +inVec[jj] * float(landcount) / float(cellcount));
			//outVec[jj] = Land[jj] * (float(sum / cellcount) );
			//outVec[jj] = float(sum / wt);
		}
	}

}






void MPAS_vis::makeSpeed(float *Ev, float *Nv, float *Sv)
{
	double e, n;
	int i;

	for (i = 0; i < XSIZE*YSIZE; ++i)
	{
		e = Ev[i];
		n = Nv[i];
		Sv[i] = sqrt(e*e + n*n);
	}

}


void MPAS_vis::calcCellIO()
{
	int r, c, ir, ic, ii;
	int step;
	float min, max, v;
	double sum;
	step = 20;

	for (r = 0; r < 35; ++r)
		for (c = 0; c < 50; ++c)
			cellIOdelta[r][c] = 0.0f;

	for (r = 0; (r < 35-1); ++r)
		for (c = 0; c < 50-1; ++c)
		{
			sum = 0;
			// TOP AND BOTTOM
			for (ic = c * step; ic < (c * step + step); ++ic)
			{
				ii = r * step * NCOLS + ic;
				sum = sum + smoothF_N[0][ii]*cosLat[r*step]; // add inputs;
				ii = (r+1) * step * NCOLS + ic;
				sum =  sum - smoothF_N[0][ii]*cosLat[r * step + step]; // weight by lat cell wid
			}	

			// LEFT AND RIGHT
			for (ir = r * step; ir < (r * step + step); ++ir)
			{
				ii = ir *NCOLS + c*step;  // LEFT
				sum = sum + smoothF_E[0][ii];  // add inputs
				ii = ir *NCOLS + c * step + step;  //RIGHT
				sum = sum - smoothF_E[0][ii]; // subtract outputs
			}

			cellIOdelta[r][c] = sum;
		}

	FILE *fp;
	fp = fopen("Vert.csv", "w");
	min = 1000.0f; max = -1000.0f;

	for (r = 0; r < 34; ++r)
	{
		for (c = 0; c < 49; ++c)
		{
			v = cellIOdelta[r][c];
			if (v < -100000000.0f) v = 0.0f;
			if (v > 100000000.0f) v = 0.0f;
			if (v > max) max = v;
			if (v < min) min = v;

			cellIOdelta[r][c] = v;

			fprintf(fp, "%5.3f ,", v);
		}
		fprintf(fp, "\n");
	}
	fclose(fp);

	cerr << "MAX MIN 2 deg cells " << max << " " << min << "\n";
}




/*

void MPAS_vis::getDTM(float **(&dd),int &R, int &C)
{
	//dem->getDTM(dd, R,C);

	
	cerr << "TEST DEM GET 250 251 " << dd[250][251] << "\n";
}
*/
void MPAS_vis::printFlowGraph(float Lat)
{
	int ii, c, Lat10 , L;
	Lat10 = Lat*10.0f;

	float sum, nflow, x;
	double sum45, cum45;
	double sumDeep, cumDeep;

	char fname[30];

	sprintf(fname, "TransPlot_%d.csv", Lat10);

	sum = 0.0f;
	cum45 = 0.0;
	cumDeep = 0.0;
	cerr << fname << "\n";
	FILE *fp;

	fp = fopen(fname, "w");

	for (c = 0; c < XSIZE; c = c + 2)
	{
		x = float(c) / 10.0f -90.0f;
		ii = Lat10*XSIZE + c;

		nflow = Nvec[10][ii];
		if(nflow > -1.9f)
		{ 
			sum = sum + nflow;
		}
		else nflow = 0.0f;

		sum45 = 0.0;
		for (L = 0; L < 45; ++L)
		{
			nflow = Nvec[L][ii];
			if (nflow > -1.9f)
				sum45 = sum45 + nflow* depthLayers[L]; // aggregate flow down to Layer20;
				
		}
		cum45 = cum45 + sum45;

		sumDeep = 0.0;
		for (L = 45; L < 80; ++L)
		{
			nflow = Nvec[L][ii];
			if (nflow > -1.9f)
				sumDeep = sumDeep + nflow* depthLayers[L]; // aggregate flow down to Layer20;

		}
		cumDeep = cumDeep + sumDeep;

		fprintf(fp, "%4.2f, %4.3f, %4.3f, %4.3f, %4.3f \n", x, nflow, sum, cum45, cumDeep);
		//fprintf(fp, "%4.2f, %4.3f, %4.3f, %4.3f\n", x, nflow, sum, cum45);
	}
	fclose(fp);
}

void MPAS_vis::computeAMOCflows(float Lat)
{
	int ii, c, Lat10, L;

	Lat10 = Lat*10.0f;

	float sum, nflow, x;
	double sum45, cum45;
	double sumDeep, cumDeep,amoc;

	float cellwid;
	cellwid = cos(Lat*3.141592 / 180.0)*40000.0 / 360.0;
	
	cellwid = cellwid *100.0; // meters wid per grid cell
	//cerr << "Cell Width in meters " << cellwid << "\n";
	sum = 0.0f;
	cum45 = 0.0;
	cumDeep = 0.0;
	cerr << "-----------------------------------\n";
	FILE *fp;

	layers[0]->computePeakNorth(Lat10, depthLayers);
	//transects[0]->computePeakNorth_2(Lat10, agNf_S, c);
	int cpLat, cpLon;
	//transects[0]->findCriticalPoint(360, cpLat, cpLon);

	for (c = 0; c < XSIZE; c = ++c)
	{
		x = float(c) / 10.0f - 90.0f;
		ii = Lat10*XSIZE + c;
		nflow = Nvec[10][ii];
		if (nflow > -10.9f)
		{
			sum = sum + nflow;
		}
		else nflow = 0.0f;

		sum45 = 0.0;
		for (L = 0; L < MID_DEP; ++L)
		{
			nflow = Nvec[L][ii];
			if (nflow > -10.9f) // what is this??
				sum45 = sum45 + nflow* depthLayers[L]; // aggregate flow down to Layer20;

		}
		cum45 = cum45 + sum45;

		sumDeep = 0.0;
		for (L = MID_DEP; L < 80; ++L)
		{
			nflow = Nvec[L][ii];
			if (nflow > -10.9f)
				sumDeep = sumDeep + nflow* depthLayers[L]; // aggregate flow down to Layer20;
		}
		cumDeep = cumDeep + sumDeep;
	}
	double max = -10000000.00;
	int AMOC_Level ;
	amoc = 0.0; AMOC_Level = N_DEPTHS;
	for (L = N_DEPTHS-1; L > -1; --L)
	{

		for (c = 0; c < XSIZE; c = ++c)
		{
			ii = Lat10*XSIZE + c;
			amoc = amoc - Nvec[L][ii] * depthLayers[L];
		}
		if (amoc > max) {
			max = amoc;
			AMOC_Level = L;
		}
	}
	cout.precision(4);	
	cerr << "Lat " << Lat10 / 10.0f << "\n";
	cerr << "MAX AMOC " << AMOC_Level << " " << max*cellwid / 1000000.0 << "\n";

	cerr <<  "Surface Flow " << cum45*cellwid/1000000.0 << "\n";
	cerr  << "Deep Flow " << cumDeep*cellwid/ 1000000.0 << "\n";
}


void MPAS_vis::drawCurrentTransect()
{
	if(viewTransect == SHALLOW)
		layers[0]->drawOrthoFlowProfile(depthLayers,currentTransect_S);
	if (viewTransect == DEEP)
		layers[1]->drawOrthoFlowProfile(depthLayers, currentTransect_D);

}

void MPAS_vis::drawSlice(float Lat)
{
	int c,ii;
	int L;
	int Lat10;

	float x, y, v;
	float mfac;
	float red, green;

	glColor3f(1.0f, 1.0f, 1.0f);
	glEnable(GL_TEXTURE_2D);

	colormaps->setColorMap();
	//glBindTexture(GL_TEXTURE_2D, colormap->tex2dNames[1]);

	Lat10 = int(Lat*10.0f);

	for (L = 0; L < IN_LAYERS; ++L)
	{
		y = -float(L)*3.5f - 10.0f;

		glBegin(GL_QUAD_STRIP);

		for (c = 0; c < XSIZE; ++c)
		{
			mfac = depthLayers[L] *0.02f;
			x = float(c);
			ii = Lat10*XSIZE + c;
			v = Nvec[L][ii] * mfac + 0.5f - 0.015;

			glTexCoord2f( 0.1f, v);

			glVertex2f(x, y);
			glVertex2f(x, y-3.5f);
		}
		glEnd();
	}
	glDisable(GL_TEXTURE_2D);

}

void MPAS_vis::drawMask(unsigned char *buff)
{
	int r, c, ii;
	float xpos;
	float red;

	glColor3f(1.0f, 1.0f, 1.0f);
	for (r = 0; r < YSIZE; ++r)
	{

		glBegin(GL_POINTS);

		for (c = 0; c < XSIZE; ++c)
		{
			ii = (r*XSIZE + c)*3;
			if (buff[ii] > 0)
			{
				glColor4ub(buff[ii], buff[ii + 1], buff[ii + 2],127);

				//red = 3.5*smoothNv[ii]		
				xpos = 500.0f + float(c - 500)*cosLat[r];
				glVertex2f(xpos, float(r));
			}

		}
		glEnd();
	}
}

void MPAS_vis::drawBG()
{
	int r, c, ii;
	float xpos;
	float red;

	glColor3f(1.0f, 1.0f, 1.0f);
	for (r = 0; r < YSIZE; ++r)
	{

		glBegin(GL_POINTS);

		for (c = 0; c < XSIZE; ++c)
		{
			ii = r*XSIZE + c;
			red = 0.8f;

			if (Land[ii] > 0.0) glColor4f(red, red, red, 0.0f);
			else glColor4f(red, red, red, 1.0f);

			//red = 3.5*smoothNv[ii]		
			xpos = 500.0f + float(c - 500)*cosLat[r];

			glVertex2f(xpos, float(r));

		}
		glEnd();
	}
}

void MPAS_vis::drawGainLoss()
{
	float x1, x2, y;
	int r, c; 
	float red,v;
	for (r = 1; r < 33; ++r)
	{
		y = r* 20.0f;
		glBegin(GL_QUAD_STRIP);
		for (c = 1; c < 48; ++c)
		{

			v = cellIOdelta[r][c];
			red = v*5.0f + 0.5f;

			glColor3f(red, 0.1, 1.0 - red);

			x1 = 500.0f + float(c * 20 - 500)*cosLat[r * 20];

			glVertex2f(x1, y);

			v = cellIOdelta[r+1][c];
			red = v*5.0f + 0.5f;
			x2 = 500.0f + float(c * 20 - 500)*cosLat[r * 20 + 20];

			glColor3f(red, 0.1, 1.0 - red);
			glVertex2f(x2, y + 20.0f);
	
		}

		glEnd();
	}
}

void MPAS_vis::drawPathTraces_S(int rtype, float trans)
{
	int i;
	float p;


	glColor3f(1.0f, 0.0f, 0.0f);

//	cerr << "drawPaths \n";

	if(rtype == ANTS)glLineWidth(1.5f);
	if (rtype == MASK)glLineWidth(4.0f);
	//if (rtype == ANTS)trans = 0.3f;

	for (i = 0; i < 200; ++i) {
		p = float(i)*0.0045 + 0.01;
		layers[0]->computeDrawPathLoop(1, 0, p, 4, trans, rtype, colorScheme);// Atlantic gyre
	}

	//if (rtype == ANTS)trans = 0.25f;
	for (i = 0; i < 50; ++i) {
		p = float(i)*0.009 + 0.25;
		layers[0]->computeDrawPathLoop(8, 2, p, 5, trans, rtype, colorScheme); // North Atlantic reverse gyre
	}	

	//if (rtype == ANTS)trans = 0.5f;
	for (i = 0; i < 100; ++i) {
		p = float(i)*0.005 + 0.25;
		layers[0]->computeDrawPathABC(8, 2, 0, p, 5, trans, END, rtype, colorScheme); // Transfer to North Atlantic reverse gyre
	}

	//if (rtype == ANTS)trans = 0.25f;
	for (i = 0; i < 60; ++i) {
		p = float(i)*0.009 + 0.25;
		layers[0]->computeDrawPathLine(4, 0, p, true, 1, trans, END, rtype, colorScheme); // Gulf steam past UK
	}

	//if (rtype == ANTS)trans = 0.25f;

	glLineWidth(5.0f);
	for (i = 0; i < 400; ++i) {
		p = float(i)*0.0020 + 0.10;
		layers[0]->computeDrawPathLine(6, 0, p,false, 1, trans,END, rtype, colorScheme); // South out 
	}

}

void MPAS_vis::drawPathTraces_D(int rtype, float trans)
{
	int i;
	float p;

	glColor3f(1.0f, 0.0f, 0.0f);

	if (rtype == ANTS)glLineWidth(3.5f);
	if (rtype == MASK)glLineWidth(6.0f);


	for (i = 0; i < 90; ++i) {
		p = float(i)*0.009 + 0.05;
		if(rtype == MASK)
			layers[1]->computeDrawPathABC(8, 0, 4, p, 3, trans, START, rtype, colorScheme); // south Flow
		else
			layers[1]->computeDrawPathABC(8, 0, 4, p, 3, trans, START, rtype, colorScheme); // south Flow
	}


	for (i = 0; i < 90; ++i) {
		p = float(i)*0.009 + 0.25;
		if (rtype == MASK)
			layers[1]->computeDrawPathLine(7, 0, p, true, 3, trans, NONE, rtype, colorScheme); // From the North
		else
			layers[1]->computeDrawPathLine(7, 0, p, true, 3, trans, NONE, rtype, colorScheme); // From the North
	}


}



void MPAS_vis::drawVecField(int Layer)
{
	int r,c,ii,i;
	float x, y1,y2, xpos;
	float red;

	glColor3f(1.0f, 1.0f, 1.0f);

	//drawBG();
	
	if (drawTraces)
	{
		//flowVis_D->drawTracesStatic();
		//flowVis_S->drawTracesStatic();
		if (Layer == DEEP) {

			//glColor3f(0.2f, 0.2f, 1.0f);
			flowVis_D->drawPlets(4.5f);
			//flowVis_D->drawAnimatedTraces(4.5);
		}
		if (Layer == SHALLOW) {
			glColor3f(1.0f, 1.0f, 0.0f);
				flowVis_S->drawPlets(2.0f);
		}
			//flowVis_S->drawAnimatedTraces(2.0);
	}

	if(drawLongTrace) flowVis_S->drawLongTrace();

	//drawPathTraces_S();  // hack for now

}

void MPAS_vis::draw_Flow_Indicators_S(int ts)
{
	flowIndicator[0]->draw(FV[1],1, ts);
	flowIndicator[1]->draw(FV[4],1, ts);
	flowIndicator[2]->draw(FV[2],0, ts);
	flowIndicator[3]->draw(FV[3],0, ts);
	flowIndicator[4]->draw(-FV[5],1, ts);

	glColor3f(0.0, 0.0, 0.8);
	glLineWidth(2.0f);

	layers[0]->drawTransectLine(0);
	layers[0]->drawTransectLine(1);
	layers[0]->drawTransectLine(2);
	layers[0]->drawTransectLine(4);
	layers[0]->drawTransectLine(5);

}

void MPAS_vis::draw_Flow_Indicators_D(int ts)
{	//Deep
	flowIndicator[5]->draw(-FV[8], 1, ts);  // need to fix the volumen
	flowIndicator[6]->draw(FV[9], 1, ts);  // need to fix the volumen
	flowIndicator[7]->draw(FV[10], 1, ts);  // need to fix the volumen
	flowIndicator[8]->draw(FV[11], 1, ts);  // need to fix the volumen

	glEnable(GL_LINE_SMOOTH);
	layers[1]->drawTransectLine(0);
	layers[1]->drawTransectLine(8);
	layers[1]->drawTransectLine(7);
	layers[1]->drawTransectLine(9);
}

void MPAS_vis::drawTransects_Paths()  // draw the paths
{
	int i;

	//layers[0]->drawOrthoFlowProfile(depthLayers, currentTransect_S);

	//flowIndicator[0]->draw(FV[1]);
	//flowIndicator[1]->draw(FV[4]);
	//flowIndicator[2]->draw(FV[2]);

	if (drawTransct  && viewTransect == SHALLOW)
	{	
		//layers[0]->drawLongTrace();
		glColor3f(0.8, 0.8, 0.0);

		layers[0]->drawPaths();

		for (i = 0; i < nTransects_S; ++i)
		{ 	glColor3f(0.0, 0.0, 0.8);
			glLineWidth(2.0f);

			if (i == currentTransect_S)
			{
				glColor3f(0.0, 0.0, 0.5);
				glLineWidth(5.0f);
			}

			layers[0]->drawTransectLine(i);
		}
	}


	if (drawTransct  && viewTransect == DEEP)
	{
		//layer->drawLongTrace();
		glColor3f(0.0, 0.8, 0.8);
		layers[1]->drawPaths();

		for (i = 0; i < nTransects_D; ++i)
		{
			glColor3f(0.0, 0.0, 0.8);
			glLineWidth(2.0f);

			if (i == currentTransect_D)
			{
				glColor3f(0.0, 0.0, 0.5);
				glLineWidth(5.0f);
			}

			layers[1]->drawTransectLine(i);
		}

	}
}

void MPAS_vis::setSamplePts(unsigned char *buff, int L, bool first)
{
	if(L==0){
		flowVis_S->setPtColors(buff, first);
		flowVis_S->makeTraces();
	}
	if (L == 1) {
		flowVis_D->setPtColors(buff,first);
		flowVis_D->makeTraces();
	}

}





