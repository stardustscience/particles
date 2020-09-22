#include <GL/glut.h>
#include <math.h>
#include <iostream>
#include "Layer.h"
#include "flowPath.h"

#define IN_LAYERS 80
#define AMOC_DEP 45


using namespace std;

#define MAX_LEN 400  // the maximum length of a transect

Layer::Layer(int nr, int nc, int L)
{
	int i;
	double angl;

	float rad = 5.0f;
	nRows = nr;
	nCols = nc;

	for (i = 0; i < NCIRC; ++i)
	{
		angl = 2.0*3.141592*double(i)/double(NCIRC - 1);

		circX[i] = sin(angl)*rad;
		circY[i] = cos(angl)*rad;
	}

	cosCor = new float[nRows];
	for (i = 0; i < nRows; ++i)
		cosCor[i] = cos((double(i) / 10.0)*3.141592 / 180.0);
	for (i = 0; i < nRows; ++i)
		cosCor[i] = cos((double(i) / 10.0)*3.141592 / 180.0);

	flowPaths = new flowPath*[50];


	trLongX = new float[LONG_TRACE];
	trLongY = new float[LONG_TRACE];
	longTraceLen = 0;
	pathCount = 0;

	for (i = 0; i < 50; ++i) startX[i] = startY[i] = 100.0f;

	trace_1 = new traceLine;
	trace_2 = new traceLine;

	smX = new float[MAX_TRACE_LEN];
	smY = new float[MAX_TRACE_LEN];

	rn = 0.6f; gn = 0.0f; bn = 0.6f;

	if (L == 0){
		rn = 1.0f; gn = 0.0f; bn = 0.0f;
	}

	if (L == 1) {
		rn = 0.3f; gn = 0.3f; bn = 1.0f;
	}


	//transOrth = new float[MAX_LEN];


}

void Layer:: setFlows(float *fx, float *fy, float *spd, float **N, float **E)
{
	smFlowX = fx;
	smFlowY = fy;
	smSpeed = spd;
	fN = N;
	fE = E;
}

// change to define
float Layer::computeTransect(float x1, float y1, float x2, float y2, int n)
{
	// define the nth Transect
	float dx, dy, len;
	float px, py;
	float ddx, ddy;
	double cddx, cddy, lm; //  lengths in meters.
	int i, ix, iy, ii;
	double sum;
	FILE *fp;
	trSect tran;

	// set the endpoints
	transects_S[n].tx1 = (x1 + 90.0f)*10.0f; 
	transects_S[n].ty1 = y1*10.0f; // in model space
	transects_S[n].tx2 = (x2 + 90.0f)*10.0f;
	transects_S[n].ty2 = y2*10.0f;


	dx = transects_S[n].tx2 - transects_S[n].tx1; 
	dy = transects_S[n].ty2 - transects_S[n].ty1;
	len = sqrt(dx*dx + dy*dy); // in Lat Lon space
	VoX = dy / len; VoY = -dx / len;
	ddx = dx / float(MAX_LEN);
	ddy = dy / float(MAX_LEN);
	sum = 0.0;

//	if (print) fp = fopen("Trans.csv", "w");

	for (i = 0; i < MAX_LEN; ++i)
	{
		px = transects_S[n].tx1 + ddx*(float(i) + 0.5);
		py = transects_S[n].ty1 + ddy*(float(i) + 0.5);
		ix = px;
		iy = py;
		ii = ix + iy*nCols;

		cddx = ddx * cosCor[int(iy)] * 11132.0;  // meters
		cddy = ddy * 11132.0;
		lm = sqrt(cddx*cddx + cddy*cddy); // sample width

		transOrth[i] = vX_S[ii] * VoX + vY_S[ii] * VoY; // dot product

//		if (print) if (i % 5 == 3)fprintf(fp, "%5.3f, %5.3f, %5.3f \n", px, py, transOrth[i] * lm);
		sum += transOrth[i] * lm; 
	}	
//	cerr << "Len seg " << lm << "\n";
//	if (print)fclose(fp);

	cerr << "ORTHOGONAL FLOW Sv " << sum/1000.0 << "\n";
	return sum / 1000.0;
}

// do this using precomputed flows
float Layer::computePeakNorth(int L10, float *layers)
{
	int L,i,ii;
	double sum, max;
	float wid;
	int maxPos;

	wid = cosCor[L10] * 11132.0;

	maxPos = 0;
	max = sum = 0;
	for (i = 10; i < 900; ++i)
	{
		ii = i + L10*1000;
		for (L = 0; L < AMOC_DEP; ++L) 
		{
			sum = sum + fN[L][ii]*layers[L]*wid;
		}
		if (max < sum) {
			max = sum;
			maxPos = i;
		}

	}
	//cerr << "LAT " << L10/10 << "\n";
	cerr << "MAX N " << float(maxPos)/10.0f -90.0f << " " << max / 1000000.0 << "\n";
	return max / 1000000.0;
}


float Layer::computePeakNorth_2(int Lat10, int Lon10, float *agN, int &col, bool north) // remove second arg
{
	int L, i, ii;
	double sum, max, min;
	float wid;
	int maxPos, minPos;

	wid = cosCor[Lat10] * 11132.0;

	//cerr << "LON10 " << Lon10 << "\n";

	maxPos = minPos = col = 0;
	max = min = sum = 0;
	for (i = 10; i < Lon10+100; ++i)
	{
		ii = i + Lat10 * 1000;

		sum = sum + agN[ii]*wid;
		if (max < sum) {
			max = sum;
			maxPos = i;
		}
		if (min > sum) {
			min = sum;
			minPos = i;
		}

	}
	if (north) {
		col = maxPos;
		//cerr << "LAT " << Lat10 / 10 << " ";
		//cerr << "MAX NORTH " << float(maxPos) / 10.0f - 90.0f << " " << Lat10/10.00f << " " << max / 1000.0 << "\n";
	}
	else {
		col = minPos;
		//cerr << "LAT " << Lat10 / 10 << " ";
		//cerr << "MIN NORTH " << float(minPos) / 10.0f - 90.0f << " " << Lat10 / 10.00f << " " << min / 1000.0 << "\n";
	}
	if (north) return max / 1000.0;
	else return -min / 1000.0f;  // NOTE REVERSED SIGN
}

void Layer::findCriticalPoint(int sLat, int sLon, int &fLat, int &fLon, bool north)
{
	int r, c;
	int maxLatPos, maxLonPos;
	float v, max = 0;
	// find the maximum of N flow in the signal

	for (r = sLat - 20; r < sLat + 20; ++r) {
		v = computePeakNorth_2(r, 900 + sLon, vY_S, c, north);
		if (v > max) {
			max = v;
			maxLatPos = r;
			maxLonPos = c;
		}
		
	}
	

	fLat = maxLatPos; 

	fLon = maxLonPos;

	cerr << "Critical Pt R C " << fLat/10.0f << " " << fLon/10.0f - 90.0f << "\n";

}

// computes total, shallow and deep aggregate flows from hi-res raw data

void Layer::setTransect_S(float x1, float y1, float x2, float y2, int n)
{
	transects_S[n].tx1 = (x1 + 90.0f)*10.0f;// in model space
	transects_S[n].ty1 = y1*10.0f; 
	transects_S[n].tx2 = (x2 + 90.0f)*10.0f;
	transects_S[n].ty2 = y2*10.0f;
}


float Layer::computeTransectHR(float x1, float y1, float x2, float y2, float *layers, int n, int mid, int SD )
{
	float dx, dy, len;
	float px, py;
	float ddx, ddy;
	double cddx, cddy, lm; //  lengths in meters.
	int i,L, ix, iy, ii;
	double sum;
	float plus[80], minus[80];
	//double 
	FILE *fp;
	//float VoX, VoY; // unit vector orthogonal to transect
	double transOrt,v;
	float sumPlus, sumMinus;

	for (i = 0; i < 80; ++i) 
		plus[i] = minus[i] = netFlow[i]= 0.0f;	 
	if (SD == SHALLOW) {
		dx = transects_S[n].tx2 - transects_S[n].tx1; dy = transects_S[n].ty2 - transects_S[n].ty1;
	}
	/*
	if (SD == DEEP) {
		dx = transects_D[n].tx2 - transects_D[n].tx1; dy = transects_D[n].ty2 - transects_D[n].ty1;
	}*/

	len = sqrt(dx*dx + dy*dy); // in Lat Lon space
	VoX = dy / len; VoY = -dx / len;
	ddx = dx / float(MAX_LEN);
	ddy = dy / float(MAX_LEN);
	sum = 0.0;
	sumPlus = sumMinus = 0.0f;
	//cerr << "VOX VOY " << VoX << " " << VoY << "\n";
	//cerr << "DDY m " << ddy * 11132.0 << "\n";

	for (i = 0; i < MAX_LEN; ++i)
	{
		if (SD == SHALLOW) {
			px = transects_S[n].tx1 + ddx*(float(i) + 0.5); py = transects_S[n].ty1 + ddy*(float(i) + 0.5);
		}
		ix = px;
		iy = py;
		ii = ix + iy*nCols;

		cddx = ddx * cosCor[int(iy)] * 11132.0;  // meters
		cddy = ddy * 11132.0;
		lm = sqrt(cddx*cddx + cddy*cddy); // sample width in meters
		for (L = 0; L < 80; ++L) // e.g. // PUT in INNER LOOP
		{
			transOrt = fE[L][ii] * VoX + fN[L][ii] * VoY; // dot product
				//if (print) if (i % 5 == 3)fprintf(fp, "%5.3f, %5.3f, %5.3f \n", px, py, transOrth[i] * lm);
			v = (transOrt * lm*layers[L]);	
			if (v < 0.0) {
				sumMinus += v;
				minus[L] += v;
			}
			else {
				sumPlus += v;
				plus[L] += v;
			}
			netFlow[L] += v;
			sum += v;
		}
	}
	//	cerr << "Len seg " << lm << "\n";
	float sp, sm;
	int Lp, Lm;
	if (SD == SHALLOW)
		transects_S[n].shallowFlow = transects_S[n].deepFlow = 0.0;

	sp = sm = 0.0f;
	for (i = 0; i < 80; ++i) { 
		sp = sp + plus[i];
		sm = sm + minus[i];
		if (sp < sumPlus / 2) Lp = i;
		if (sm > sumMinus / 2) Lm = i;
		//if (SD == SHALLOW){

		// change this
		if(i<AMOC_DEP)transects_S[n].shallowFlow += netFlow[i] / 1000000.0;
			else transects_S[n].deepFlow += netFlow[i] / 1000000.0;

		netFlow[i] = netFlow[i]/1000000.0;
	}
	//cerr << "ORTHOGONAL FLOW Sv " << sum / 1000000.0 << "  " << sumPlus / 1000000.0 << " " << sumMinus / 1000000.0 << "\n";
	if (SD == SHALLOW)cerr << "Shallow Deep   " << transects_S[n].shallowFlow << " " << transects_S[n].deepFlow << "\n";

	return sum / 1000000.0;
}

float Layer::computeTransectHR_2(int n)
{
	float dx, dy, len;
	float px, py;
	float ddx, ddy;
	double cddx, cddy, lm; //  lengths in meters.
	int i, L, ix, iy, ii;
	double sumS, sumD;

	double transOrt, vS, vD;

	//if (SD == SHALLOW) {
	dx = transects_S[n].tx2 - transects_S[n].tx1; dy = transects_S[n].ty2 - transects_S[n].ty1;
		//cerr << "TR " << transects_S[n].tx1 << " " << transects_S[n].ty1 << " ";
		//cerr << "TR " << transects_S[n].tx2 << " " << transects_S[n].ty2 << "\n ";
//	}

	len = sqrt(dx*dx + dy*dy); // in Lat Lon space
	VoX = dy / len; VoY = -dx / len;
	ddx = dx / float(MAX_LEN);
	ddy = dy / float(MAX_LEN);
	sumS = sumD = 0.0;

	//cerr << "Ortho Vector " << VoX << " " << VoY << "\n";



	for (i = 0; i < MAX_LEN; ++i){
		px = transects_S[n].tx1 + ddx*(float(i) + 0.5); 
		py = transects_S[n].ty1 + ddy*(float(i) + 0.5);
		ix = px;
		iy = py;
		ii = ix + iy*nCols;

		cddx = ddx * cosCor[int(iy)] * 11132.0;  // meters
		cddy = ddy * 11132.0;
		lm = sqrt(cddx*cddx + cddy*cddy); // sample width in meters

		sumS += (vX_S[ii] * VoX + vY_S[ii]*VoY)*lm;// 
		//sumD += (vX_D[ii] * VoX + vY_D[ii]*VoY)*lm;
		
	}

	transects_S[n].shallowFlow = sumS/1000.0f; // layer flow.

	return sumS / 1000.0;
}

float Layer::computePeakTransectHR(int n, float &maxp, float &pxMax, float &pyMax)
{
	float dx, dy, len;
	float px, py;
	float ddx, ddy;
	double cddx, cddy, lm; //  lengths in meters.
	int i, L, ix, iy, ii;
	double sumS, sumD;
	double transOrt, vS, vD;
	float maxTrans;
	int maxIndex;

	dx = transects_S[n].tx2 - transects_S[n].tx1; dy = transects_S[n].ty2 - transects_S[n].ty1;

	len = sqrt(dx*dx + dy*dy); // in Lat Lon space
	VoX = dy / len; VoY = -dx / len;
	ddx = dx / float(MAX_LEN);
	ddy = dy / float(MAX_LEN);
	sumS = sumD = 0.0;
	maxTrans = 0.0;
	maxIndex = 0;
	//cerr << "Ortho Vector " << VoX << " " << VoY << "\n";

	i = 45;  //HACK

	//for (i = 0; i < MAX_LEN; ++i) {
		px = transects_S[n].tx1 + ddx*(float(i) + 0.5);
		py = transects_S[n].ty1 + ddy*(float(i) + 0.5);
		ix = px;
		iy = py;
		ii = ix + iy*nCols;

		cddx = ddx * cosCor[int(iy)] * 11132.0;  // meters
		cddy = ddy * 11132.0;
		lm = sqrt(cddx*cddx + cddy*cddy); // sample width in meters

		sumS -= (vX_S[ii] * VoX + vY_S[ii] * VoY)*lm;// 
		if (i % 40 == 39) cerr <<i << " " <<sumS / 1000.0 << "\n";

		if (sumS > maxTrans)
		{
			maxTrans = sumS;
			maxIndex = i;
			pxMax = px;
			pyMax = py;
		}

	//}

	transects_S[n].shallowFlow = sumS / 1000.0f; // layer flow.

	cerr << "MAX " << maxTrans/1000.0f << " " << pxMax << " " << pyMax <<"\n";

	transects_S[0].tx2 = pxMax;
	transects_S[0].ty2 = pyMax;
	transects_S[1].tx1 = pxMax;
	transects_S[1].ty1 = pyMax;

	return sumS / 1000.0;
}

void Layer::getTransectFlow(int n, float &s, float &d)
{
	
	//cerr << "\n Flow Volume " << transects_S[n].shallowFlow << "\n";
	//cerr << "DeepFlow    " << transects_S[n].deepFlow << "\n";
	s = transects_S[n].shallowFlow;
	d = transects_S[n].deepFlow;


}

void Layer::drawOrthoFlowProfile( float *layers, int n)
{
	float dx, dy, tLen;
	float px, py,x1,x2,px1,y, xp1,xp2;
	float ddx, ddy, tv;
	double cddx, cddy, lm; //  lengths in meters.
	int i, ix, iy, ii, L;
	int nCols = 1000;
	int steps;
	double sum;
	FILE *fp;
	float VoX, VoY; // unit vector orthogonal to transect
					// set the endpoints
	float xV, yV; // the east and North Vectors
	float OrthTrans;


	glColor3f(1.0f, 1.0f, 1.0f);
	glEnable(GL_TEXTURE_2D);
// step in model coorinates correct for the flow

	dx = transects_S[n].tx2 - transects_S[n].tx1; 
	dy = transects_S[n].ty2 - transects_S[n].ty1;
	tLen = sqrt(dx*dx + dy*dy); // len in lat lon space 1/10th of a deg
	steps = int(tLen);// *10.0f);

	VoX = dy / tLen; VoY = -dx / tLen;  // orthogonal vector
	ddx = dx / float(steps);
	ddy = dy / float(steps);
	for (i = 0; i < steps; ++i)	
	{			
		px1 = transects_S[n].tx1 + ddx*(float(i) + 0.5);
	//	px2 = tx1 + ddx*(float(i+1) + 0.5);
		//xp1 = 500.0f + float(rx1 - 500)*cosCor[int(ty1)];
		//xp2 = 500.0f + float(rx2 - 500)*cosCor[int(ty2)];
		//pdx = xp2-xp1; // the physical x distance 

		py = transects_S[n].ty1 + ddy*(float(i) + 0.5);
		ix = px1;
		iy = py;
		ii = ix + iy*nCols;

		x1 = float(i) + 100.0f;
		x2 = float(i+1) + 100.0f;

		glBegin(GL_QUAD_STRIP);

		for (L = 0; L < IN_LAYERS; ++L) 
		{	
			y = -float(L)*3.5f - 10.0f;

			yV = fN[L][ii];
			xV = fE[L][ii];
			//cddx = ddx * cosCor[int(iy)] * 11132.0;  // meters
			//cddy = ddy * 11132.0;
			//lm = sqrt(cddx*cddx + cddy*cddy); // sample width
			OrthTrans = xV * VoX + yV * VoY; // dot product
			tv = OrthTrans*layers[L] * 0.02f +0.5 - 0.015;
			glTexCoord2f(0.1f, tv);

			glVertex2f(x1, y);
			glVertex2f(x2, y);
		}
		glEnd();
	}

	glDisable(GL_TEXTURE_2D);
	glColor3f(0.0, 0.0, 0.0);
	glBegin(GL_LINE_STRIP);
	for (L = 0; L < IN_LAYERS; ++L)
	{
		y = -float(L)*3.5f - 10.0f;
		x1 = 600.0f + netFlow[L]*50.0f;
		glVertex2f(x1, y);
	}
	glEnd();
	glColor3f(0.5, 0.5, 0.2);
	glBegin(GL_LINES);
	for (L = 0; L < IN_LAYERS; ++L)
	{
		glVertex2f(600.0f, -10.0f);
		glVertex2f(600.0f, -350.0f);
	}
	glEnd();
	//cerr << "ORTHOGONAL FLOW Sv " << sum / 1000.0 << "\n";
}

void Layer::drawTransectLine(int n)
{
	int i;
	float xp1, xp2;


	xp1 = 500.0f + float(transects_S[n].tx1 - 500)*cosCor[int(transects_S[n].ty1)];
	xp2 = 500.0f + float(transects_S[n].tx2 - 500)*cosCor[int(transects_S[n].ty2)];

	//glColor3f(0.5f, 0.5f, 1.0f);
	glColor3f(rn,gn,bn);
	glBegin(GL_TRIANGLE_FAN);
	glVertex2f(xp1, transects_S[n].ty1);
	for (i = 0; i < NCIRC; ++i)
		glVertex2f(xp1+circX[i], transects_S[n].ty1+circY[i]);
	glEnd();

	glBegin(GL_TRIANGLE_FAN);
	glVertex2f(xp2, transects_S[n].ty2);
	for (i = 0; i < NCIRC; ++i)
		glVertex2f(xp2 + circX[i], transects_S[n].ty2 + circY[i]);
	glEnd();

	glBegin(GL_LINES);
	
	glVertex2f(xp1, transects_S[n].ty1);
	glColor3f(rn,gn,bn);
	glVertex2f(xp2, transects_S[n].ty2);
	glEnd();

	glColor3f(1.0, 1.0, 1.0);  // white circle
	glBegin(GL_LINE_STRIP);
	for (i = 0; i < NCIRC; ++i)
		glVertex2f(xp1 + circX[i], transects_S[n].ty1 + circY[i]);
	glEnd();
	glColor3f(1.0, 1.0, 1.0);
	glBegin(GL_LINE_STRIP);
	for (i = 0; i < NCIRC; ++i)
		glVertex2f(xp2 + circX[i], transects_S[n].ty2 + circY[i]);
	glEnd();


}

/*
void Layer::drawTransect_D(int n)  // consolidate
{
	int i;
	float xp1, xp2;


	xp1 = 500.0f + float(transects_D[n].tx1 - 500)*cosCor[int(transects_D[n].ty1)];
	xp2 = 500.0f + float(transects_D[n].tx2 - 500)*cosCor[int(transects_D[n].ty2)];

	glBegin(GL_TRIANGLE_FAN);
	glVertex2f(xp1, transects_D[n].ty1);
	for (i = 0; i < NCIRC; ++i)
		glVertex2f(xp1 + circX[i], transects_D[n].ty1 + circY[i]);
	glEnd();

	glColor3f(1.0, 0.0, 0.0);

	glBegin(GL_TRIANGLE_FAN);
	glVertex2f(xp2, transects_D[n].ty2);
	for (i = 0; i < NCIRC; ++i)
		glVertex2f(xp2 + circX[i], transects_D[n].ty2 + circY[i]);
	glEnd();

	glBegin(GL_LINES);
	glColor3f(0.0, 0.0, 0.0);
	glVertex2f(xp1, transects_D[n].ty1);
	glColor3f(1.0, 0.0, 0.0);
	glVertex2f(xp2, transects_D[n].ty2);
	glEnd();

	glBegin(GL_LINES);
	glColor3f(1.0, 0.0, 0.0);
	glVertex2f(xp1, transects_D[n].ty1);
	glVertex2f(xp1 + VoX*10.0f, transects_D[n].ty1 + VoY*10.0f);
	glEnd();

}
*/

void Layer::computeDrawPathLine(int t1, int t2, float p, bool dir, int colors, float tr, int startEndExt, int rtype, int scm)
{
	int i, is;
	float t1dx, t1dy;
	bool p1;
	float startx, starty;
	float avg1, c; // average flows
	float red, green, blue;

	t1dx = transects_S[t1].tx2 - transects_S[t1].tx1;
	t1dy = transects_S[t1].ty2 - transects_S[t1].ty1;
	setColors(colors, scm,red, green, blue);
	glColor4f(red, green, blue, tr);

	startx = transects_S[t1].tx1 + p*t1dx; // half way
	starty = transects_S[t1].ty1 + p*t1dy;
	//p1 = traceToTransect(startx, starty, trLongX, trLongY, t2, avg1, dir);
	p1 = traceToTransect_2(startx, starty, trace_1, t2, avg1, dir);

	if (p1)
	{
		if(startEndExt > 0) extendTrace_2(trace_1, startEndExt, 5, 20, dir);
		smoothTrace(trace_1);
		drawTrace(trace_1, rtype);
		//drawToTransect(startx, starty, trLongX, trLongY, t2, dir);
	}

}


bool Layer::computeDrawPathLoop(int t1, int t2, float p, int colors, float tr, int rtype, int scm)
{
	int i, is;
	float t1dx, t1dy;
	bool p1, p2, found;
	float startx, starty;
	float avg1, avg2,c; // average flows
	float red, green, blue;

	t1dx = transects_S[t1].tx2 - transects_S[t1].tx1;
	t1dy = transects_S[t1].ty2 - transects_S[t1].ty1;

	setColors(colors, scm,red, green, blue);
	glColor4f(red, green, blue, tr);


	found = false;

		startx = transects_S[t1].tx1 + p*t1dx; // half way
		starty = transects_S[t1].ty1 + p*t1dy;
		//p1 = traceToTransect(startx, starty, trLongX, trLongY, t2, avg1, false);
		p1 = traceToTransect_2(startx, starty, trace_1, t2, avg1, false);
		p2 = traceToTransect_2(startx, starty, trace_2, t2, avg1, true);
		//p2 = traceToTransect(startx, starty, trLongX, trLongY, t2, avg2, true);

		//c = 1.0 - 5.0*(avg1 + avg2);
		//if (c < 0.0f) c = 0.0f;
		//if (c > 1.0f) c = 1.0f;

		if (p1 && p2)
		{
			smoothTrace(trace_1);
			drawTrace(trace_1, rtype);
			smoothTrace(trace_2);
			drawTrace(trace_2, rtype);

			found = true;

			//drawTrace_2(trace_1);
			//drawTrace_2(trace_2);
			/*

			
			glColor4f(0.0f,0.0f,0.0f, tr);
			glLineWidth(2.5f);
			drawTrace(trace_1);
			glColor4f(red, green, blue, 1.0);
			glLineWidth(1.2f);
			drawTrace(trace_1);

glPushMatrix();
			glColor4f(0.0f, 0.0f, 0.0f, tr);
			glLineWidth(2.5f);
			drawTrace(trace_2);
			glColor4f(red, green, blue, 1.0);
			glLineWidth(1.2f);
			drawTrace(trace_2);
			*/
			//drawToTransect(startx, starty, trLongX, trLongY, t2, false);
			//drawToTransect(startx, starty, trLongX, trLongY, t2, true);
			
		}

		return found;

}

void Layer::computeDrawPathABC(int t1, int t2, int t3, float p, int colors, float tr, int startEndExt, int rtype,int scm)
{
	int i, is;
	float t1dx, t1dy;
	bool p1, p2, found;
	float startx, starty;
	float avg1, avg2, c; // average flows
	float red, green, blue;

	t1dx = transects_S[t1].tx2 - transects_S[t1].tx1;
	t1dy = transects_S[t1].ty2 - transects_S[t1].ty1;
	
	setColors(colors, scm, red, green, blue);
	glColor4f(red, green, blue,  tr);


	found = false;

	startx = transects_S[t1].tx1 + p*t1dx; // half way
	starty = transects_S[t1].ty1 + p*t1dy;
//	p1 = traceToTransect(startx, starty, trLongX, trLongY, t2, avg1, false);

	p1 = traceToTransect_2(startx, starty, trace_1, t2, avg1, false); //lower section
	p2 = traceToTransect_2(startx, starty, trace_2, t3, avg1, true);  // upper section
	//p2 = traceToTransect(startx, starty, trLongX, trLongY, t3, avg2, true);
	//p2 = true;

	//c = 1.0 - 5.0*(avg1 + avg2);
	//if (c < 0.0f) c = 0.0f;
	//if (c > 1.0f) c = 1.0f;

	if (p1 && p2)
	{
		//if (startEndExt == START) extendTrace_2(trace_1, START, 10, 80, true);  
		if (startEndExt == END) extendTrace_2(trace_1, END, 10, 80, false);  // works for the gulf stream n transferf
		smoothTrace(trace_1);
		drawTrace(trace_1,rtype);

		//if(startEndExt == END) extendTrace_2(trace_2, END, 10, 80, true);  
		if (startEndExt == START) extendTrace_2(trace_2, END, 10, 80, true);// works for the deep transfer
		smoothTrace(trace_2);
		drawTrace(trace_2,rtype);
		//drawToTransect(startx, starty, trLongX, trLongY, t2, false);
		//drawToTransect(startx, starty, trLongX, trLongY, t3, true);

	}

}


void Layer::setColors(int colors, int scheme, float &red,float  &green,float &blue)
{
	red = green = blue = 0.0;

	if(scheme == 1)
	{
	
		if (colors == 1) { red = 1.0f; }
		if (colors == 2) { green = 1.0f; }
		if (colors == 3) { blue = 1.0f, red = 0.2f, green = 0.4f; }
		if (colors == 4) { red = 1.0f, green = 0.6f; }
		if (colors == 5) { red = 0.9f, green = 0.2f; blue = 1.0f; }
		if (colors == 6) { red = 0.0f, green = 0.1f; blue = 0.5f; }
	 }
	if (scheme == 2) {
		if (colors == 1) { red = 1.0f; green = 0.4; blue = 0.4; }
		if (colors == 2) { green = 1.0f; blue = 0.05; }
		if (colors == 3) { blue = 0.3f, red = 0.05f, green = 0.0f; }
		if (colors == 4) { red = 1.0, green = 0.6; blue = 0.05; }
		if (colors == 5) { red = 1.0, blue = 0.8; }
		if (colors == 6) { red = 0.0f, green = 0.1f; blue = 0.5f; }
	}
}


void Layer::computePathLoop(int t1, int t2, int tpe, float flowVol)
{
	int i,is;
	float t1dx, t1dy;
	bool p1, p2, found;
	float startx, starty;

	float avg1, avg2;

	t1dx = transects_S[t1].tx2 - transects_S[t1].tx1;
	t1dy = transects_S[t1].ty2 - transects_S[t1].ty1;

	t1dx *= 0.05f; t1dy *= 0.05f;
	found = false;


	if (tpe == MID) is = 9;
	if (tpe == LAST) is = 19;
	for (i = is; i >0; --i) // start in the middle
	{
		startx = transects_S[t1].tx1 + t1dx*float(i); // half way
		starty = transects_S[t1].ty1 + t1dy*float(i);
		p1 = traceToTransect(startx, starty, trLongX, trLongY, t2, avg1, false);
		p2 = traceToTransect(startx, starty, trLongX, trLongY, t2, avg2,true);
		if (p1 && p2)
		{
			found = true;
			if (traceToTransect(startx, starty, trLongX, trLongY, t2, avg1, false));
			flowPaths[pathCount] = new flowPath(trLongX, trLongY, longTraceLen, flowVol, 0.0);
			flowPaths[pathCount]->makePath();
			flowPaths[pathCount]->smooth(7);
			++pathCount;
			if (traceToTransect(startx, starty, trLongX, trLongY, t2, avg2, true));
			flowPaths[pathCount] = new flowPath(trLongX, trLongY, longTraceLen, flowVol, 0.0);
			flowPaths[pathCount]->smooth(7);
			flowPaths[pathCount]->makePath();

			++pathCount;
		}
		if (found) break;
	}

	startX[pathCount - 1] = startx;
	startY[pathCount - 1] = starty;

}

void Layer::computePath2part(int t1, int t2, int t3, int tpe, float flow)
{
	int i, is;
	float t1dx, t1dy;
	bool p1, p2, found;
	float startx, starty;
	float avg1, avg2;

	t1dx = transects_S[t1].tx2 - transects_S[t1].tx1;
	t1dy = transects_S[t1].ty2 - transects_S[t1].ty1;

	t1dx *= 0.05f; t1dy *= 0.05f;
	found = false;

	if (tpe == MID) is = 9;
	if (tpe == LAST) is = 19;
	for (i = is; i >0; --i) // start in the middle
	{
		startx = transects_S[t1].tx1 + t1dx*float(i); // half way
		starty = transects_S[t1].ty1 + t1dy*float(i);
		p1 = traceToTransect(startx, starty, trLongX, trLongY, t2, avg1,false);
		p2 = traceToTransect(startx, starty, trLongX, trLongY, t3, avg2, true);
		if (p1 && p2)
		{
			found = true;
			if (traceToTransect(startx, starty, trLongX, trLongY, t2, avg1, false));
			flowPaths[pathCount] = new flowPath(trLongX, trLongY, longTraceLen, flow, 0.0);
			flowPaths[pathCount]->makePath();
			flowPaths[pathCount]->smooth(7);
			++pathCount;
			if (traceToTransect(startx, starty, trLongX, trLongY, t3, avg2, true));
			flowPaths[pathCount] = new flowPath(trLongX, trLongY, longTraceLen, flow, 0.0);
			flowPaths[pathCount]->smooth(7);
			flowPaths[pathCount]->makePath();

			++pathCount;
		}
		if (found) break;
	}
	startX[pathCount - 1] = startx;
	startY[pathCount - 1] = starty;
}

void Layer::computePathline(int t1, int t2, int tpe, bool reverse, float flow)
{
	int i,ip;
	float t1dx, t1dy;
	float startx, starty;
	bool isfirst;
	int firstPath, lastPath;
	float avg1;
	isfirst = false;

	lastPath = firstPath = 0;
	// find a middle path

	t1dx = transects_S[t1].tx2 - transects_S[t1].tx1;
	t1dy = transects_S[t1].ty2 - transects_S[t1].ty1;

	t1dx *= 0.05f; t1dy *= 0.05f;
	//reverse = true;

	ip = 1;
	if (tpe == MID) ip = 9;

	cerr << "COMPUTE PATH LINE \n";  // improve this
	for (i = ip; i <19; ++i)
	{
		startx = transects_S[t1].tx1 + t1dx*float(i); // half way
		starty = transects_S[t1].ty1 + t1dy*float(i);

		if (traceToTransect(startx, starty, trLongX, trLongY, t2, avg1, reverse))
		{
			if (!isfirst) { firstPath = i; isfirst = true; }
			lastPath = i;
		}
	}

	if(tpe == FIRST) {
		startx = transects_S[t1].tx1 + t1dx*float(firstPath); // half way
		starty = transects_S[t1].ty1 + t1dy*float(firstPath);
		if (traceToTransect(startx, starty, trLongX, trLongY, t2,avg1,  reverse));
		flowPaths[pathCount] = new flowPath(trLongX, trLongY, longTraceLen, flow, 0.0);
		flowPaths[pathCount]->makePath();
		++pathCount;
	}

	if (tpe == MID) {
		int mid = firstPath;
		startx = transects_S[t1].tx1 + t1dx*float(mid); // half way
		starty = transects_S[t1].ty1 + t1dy*float(mid);
		if (traceToTransect(startx, starty, trLongX, trLongY, t2, avg1, reverse)) {
			flowPaths[pathCount] = new flowPath(trLongX, trLongY, longTraceLen, flow, 0.0);
			flowPaths[pathCount]->smooth(7);
			flowPaths[pathCount]->makePath();

			++pathCount;
		}
	}

	if (tpe == LAST) {

		startx = transects_S[t1].tx1 + t1dx*float(lastPath); // half way
		starty = transects_S[t1].ty1 + t1dy*float(lastPath);
		if (traceToTransect(startx, starty, trLongX, trLongY, t2, avg1, reverse));
		flowPaths[pathCount] = new flowPath(trLongX, trLongY, longTraceLen, flow, 0.0);		
		flowPaths[pathCount]->smooth(7);
		flowPaths[pathCount]->makePath();

		++pathCount;
	}
	startX[pathCount - 1] = startx;
	startY[pathCount - 1] = starty;

	cerr << "FIRST LAST " << firstPath << " " << lastPath << "\n";
	
// logic for a loop 
	// compute a pathline from t1 to t2
}
/*
void Layer::drawLongTrace()
{
	int i;
	glLineWidth(5.0f);
	glColor3f(0.0f, 0.99f, 0.0f);
	glBegin(GL_LINE_STRIP);

	for (i = 0; i < longTraceLen; ++i)   // note this can be optimized
	{
		//xpos = 500.0f + float(trLongX[i] - 500)*cosCor[int(trLongY[i])];
		glVertex2f(trLongX[i], trLongY[i]);
	}
	glEnd();
	glLineWidth(1.3f);
}
*/
void Layer::drawPaths()
{
	int i;

	//cerr << "PATHS " << pathCount << "\n";
	for (i = 0; i < pathCount; ++i)
		flowPaths[i]->drawPath(false);
}


// compute line to stop at a particular transect.

bool Layer::traceToTransect(float xs, float ys, float *xt, float *yt, int t2, float &avgFlow, bool reverse)
{
	int i, ix, iy, ii;
	float xp, yp, xp2, yp2;
	float dx1, dy1, dx2, dy2;  // for predictor corrector
	float grad, constant; // the line through the transect
	bool success = false;
	float sign = 1.0f;
	if (reverse) sign = -1.0f;

	float Tscl = 5.0;
	float t2dx, t2dy, t2Len; 

	float sumFlow=0; // flow on trace;
	int count=0;

	float tSign, otSign, v;  // the sign of a point relative to the transect
	bool exit= false;
	tSign = otSign = 1.0;

	t2dx = transects_S[t2].tx2 - transects_S[t2].tx1;
	t2dy = transects_S[t2].ty2 - transects_S[t2].ty1;
	//t2len = sqrt(t2dx*t2dx + d2dy*d2dy);
	grad = t2dy / t2dx; // doesnt work for vertical transects_S
	constant = transects_S[t2].ty1 - grad*transects_S[t2].tx1;
	// eqn through transect pts y = const + grad*x
	float cc; // cosine correction
	//	TRACE FORWARD
	// need to interpolate at some point
	xp = xs;
	yp = ys;

	xt[0] = xp;
	yt[0] = yp;

	for (i = 0; i<LONG_TRACE; ++i)   // note this can be optimized
	{
		ix = int(xp);
		iy = int(yp);

		if (ix < 3)ix = 3;// exit = true;
		if (iy < 4) break;// exit = true;
		if (ix > nCols - 2)  break;// exit = true;
		if (iy > nRows - 2) break;// exit = true;
		
		// NOW do the Transect Test
		v = yp - grad*xp;
		if (v > constant) tSign = 1.0;
		else tSign = -1.0;
		if (otSign*tSign < 0.0) {
			//cerr << "Sign Change " << i << "\n";

			// note the i> 10 is to allow for full loops.
			if(xp > transects_S[t2].tx1 && xp < transects_S[t2].tx2 && i>10){
				success = exit = true;}
		}
		otSign = tSign;
		if (exit) break;

		cc = cosCor[iy];
		ii = iy*nCols + ix;

		dx1 = sign*smFlowX[ii] / (cc*(smSpeed[ii] + 0.01));
		dy1 = sign*smFlowY[ii] / (smSpeed[ii] + 0.01);

		xp2 = xp + dx1*Tscl;
		yp2 = yp + dy1*Tscl;

		ix = int(xp2);
		iy = int(yp2);  if (iy > 699) iy = 699;
		cc = cosCor[iy];
		ii = iy*nCols + ix;

		if (ii < 0) ii = 10;
		dx2 = sign*smFlowX[ii] / (cc*(smSpeed[ii] + 0.01));
		dy2 = sign*smFlowY[ii] / (smSpeed[ii] + 0.01);

		xp = xp + 0.5*(dx2 + dx1)*Tscl;
		yp = yp + 0.5*(dy1 + dy2)*Tscl;

		sumFlow += cc*smSpeed[ii];

		xt[i+1] = xp;
		yt[i+1] = yp;
	}
	longTraceLen = i+1 ;

	for (i = 0; i<longTraceLen; ++i)
		xt[i] = 500.0f + float(xt[i] - 500)*cosCor[int(yt[i])];

	avgFlow = sumFlow/float(longTraceLen);

	return success;
}


// new version includes line extension
bool Layer::traceToTransect_2(float xs, float ys, traceLine *tr, int t2, float &avgFlow, bool reverse)
{
	int i, ix, iy, ii;
	float xp, yp, xp2, yp2;
	float dx1, dy1, dx2, dy2;  // for predictor corrector
	float grad, constant; // the line through the transect
	bool success = false;
	float sign = 1.0f;
	if (reverse) sign = -1.0f;

	float Tscl = 5.0;
	float t2dx, t2dy, t2Len;

	float sumFlow = 0; // flow on trace;
	int count = 0;

	float tSign, otSign, v;  // the sign of a point relative to the transect
	bool exit = false;
	tSign = otSign = 1.0;

	t2dx = transects_S[t2].tx2 - transects_S[t2].tx1;
	t2dy = transects_S[t2].ty2 - transects_S[t2].ty1;
	//t2len = sqrt(t2dx*t2dx + d2dy*d2dy);
	grad = t2dy / t2dx; // doesnt work for vertical transects_S
	constant = transects_S[t2].ty1 - grad*transects_S[t2].tx1;
	// eqn through transect pts y = const + grad*x
	float cc; // cosine correction
			  //	TRACE FORWARD
			  // need to interpolate at some point
	xp = xs;
	yp = ys;
	tr->start = TRACE_OFFSET;
	tr->traceX[TRACE_OFFSET] = xp;
	tr->traceY[TRACE_OFFSET] = yp;
	//xt[0] = xp;
	//yt[0] = yp;

	for (i = 0; i<LONG_TRACE; ++i)   // note this can be optimized
	{
		ix = int(xp);
		iy = int(yp);

		if (ix < 3)ix = 3;// exit = true;
		if (iy < 4) break;// exit = true;
		if (ix > nCols - 2)  break;// exit = true;
		if (iy > nRows - 2) break;// exit = true;
								  // NOW do the Transect Test
		v = yp - grad*xp;
		if (v > constant) tSign = 1.0;
		else tSign = -1.0;
		if (otSign*tSign < 0.0) {

			// note the i> 10 is to allow for full loops.
			if (xp > transects_S[t2].tx1 && xp < transects_S[t2].tx2 && i>10) {
				success = exit = true;
			}
		}
		otSign = tSign;
		if (exit) break;

		cc = cosCor[iy];
		ii = iy*nCols + ix;

		dx1 = sign*smFlowX[ii] / (cc*(smSpeed[ii] + 0.01));
		dy1 = sign*smFlowY[ii] / (smSpeed[ii] + 0.01);

		xp2 = xp + dx1*Tscl;
		yp2 = yp + dy1*Tscl;

		ix = int(xp2);
		iy = int(yp2);  if (iy > 699) iy = 699;
		cc = cosCor[iy];
		ii = iy*nCols + ix;

		if (ii < 0) ii = 10;
		dx2 = sign*smFlowX[ii] / (cc*(smSpeed[ii] + 0.01));
		dy2 = sign*smFlowY[ii] / (smSpeed[ii] + 0.01);

		xp = xp + 0.5*(dx2 + dx1)*Tscl;
		yp = yp + 0.5*(dy1 + dy2)*Tscl;

		sumFlow += cc*smSpeed[ii]; 

		tr->traceX[TRACE_OFFSET +i + 1] = xp;

		tr->traceY[TRACE_OFFSET + i + 1] = yp;

	//	xt[i + 1] = xp;
	//	yt[i + 1] = yp;
	}

	tr->end  = TRACE_OFFSET + i + 1;

	avgFlow = sumFlow / float(longTraceLen);

	return success;
}

void Layer::extendTrace_2(traceLine *tr, int startEnd, int n, int randval, bool reverse)
{
	int i, ii, ti, ix, iy, isign;
	int nsteps;
	float xp, yp, xp2, yp2;
	float dx1, dy1, dx2, dy2;  // for predictor corrector
	float cc;
	float sign = 1.0f;
	float Tscl = 5.0;
	if (reverse) sign = -1.0f;

	if (n == 0)return;


	nsteps = n + rand() % randval;


	if (startEnd == START)
	{
		isign = -1;
		ti = tr->start;
		tr->start = tr->start - nsteps;
	}
	else {
		isign = 1;
		ti = tr->end - 1;
		tr->end = tr->end + nsteps - 2;
	}

	xp = tr->traceX[ti];
	yp = tr->traceY[ti];

	for (i = 0; i < nsteps; ++i)
	{
		ix = int(xp);
		iy = int(yp);
		cc = cosCor[iy];
		ii = iy*nCols + ix;

		dx1 = sign*smFlowX[ii] / (cc*(smSpeed[ii] + 0.01));
		dy1 = sign*smFlowY[ii] / (smSpeed[ii] + 0.01);

		xp2 = xp + dx1*Tscl;
		yp2 = yp + dy1*Tscl;

		ix = int(xp2);
		iy = int(yp2);  if (iy > 699) iy = 699;
		cc = cosCor[iy];
		ii = iy*nCols + ix;

		if (ii < 0) ii = 10;
		dx2 = sign*smFlowX[ii] / (cc*(smSpeed[ii] + 0.01));
		dy2 = sign*smFlowY[ii] / (smSpeed[ii] + 0.01);

		xp = xp + 0.5*(dx2 + dx1)*Tscl;
		yp = yp + 0.5*(dy1 + dy2)*Tscl;

		tr->traceX[ti] = xp;
		tr->traceY[ti] = yp;

		ti = ti + isign;
	}
	

}

void Layer::drawTrace_2(traceLine *tr)
{
	int i,j;
	float xx, yy;


	for (i = tr->start; i < tr->end; i=i+10)
	{		
		if (i < (tr->end - 6)) {
			glBegin(GL_LINE_STRIP);
			for (j = 0; j < 6; ++j) {
				yy = tr->traceY[i + j];
				xx = 500.0f + (tr->traceX[i + j] - 500.0f)*cosCor[int(yy)];
				glVertex2f(xx, yy);
			}
			glEnd();
		}
	}
}


void Layer::smoothTrace(traceLine *tr)
{
	int i,j, k2;
	float sumX, sumY;

	k2 = 11;


	for (i = tr->start; i < tr->end - k2; i = ++i)
	{
		sumX = sumY = 0.0;
		for (j = 0; j < k2; ++j) {
			sumX += tr->traceX[i + j];
			sumY += tr->traceY[i + j];
		}
		smX[i + 5] = sumX / 11.0f;
		smY[i + 5] = sumY / 11.0f;
	}


	for (i = tr->start + 5; i < tr->end - 6; i = ++i)
	{
		tr->traceX[i] = smX[i];
		tr->traceY[i] = smY[i];
	}

}


void Layer::drawTrace(traceLine *tr, int rtype)
{
	int i;
	float xx,yy;

	glBegin(GL_LINE_STRIP);

	if(rtype == ANTS)
		for ( i = tr->start; i < tr->end; ++i)
		{
			yy = tr->traceY[i];
			xx = 500.0f + (tr->traceX[i] - 500.0f)*cosCor[int(yy)];
			glVertex2f(xx, yy);
		}

	if (rtype == MASK)
		for (i = tr->start; i < tr->end; ++i)
		{
			yy = tr->traceY[i];
			xx = tr->traceX[i];
			glVertex2f(xx, yy);
		}
	glEnd();
}

bool Layer::drawToTransect(float xs, float ys, float *xt, float *yt, int t2, bool reverse)
{
	int i, ix, iy, ii;
	float xp, yp, xp2, yp2,xx;
	float dx1, dy1, dx2, dy2;  // for predictor corrector
	float grad, constant; // the line through the transect
	bool success = false;
	float sign = 1.0f;
	if (reverse) sign = -1.0f;

	float Tscl = 5.0;
	float t2dx, t2dy, t2Len; //

	float tSign, otSign, v;  // the sign of a point relative to the transect
	bool exit = false;
	tSign = otSign = 1.0;

	t2dx = transects_S[t2].tx2 - transects_S[t2].tx1;
	t2dy = transects_S[t2].ty2 - transects_S[t2].ty1;
	//t2len = sqrt(t2dx*t2dx + d2dy*d2dy);
	grad = t2dy / t2dx; // doesnt work for vertical transects_S
	constant = transects_S[t2].ty1 - grad*transects_S[t2].tx1;
	// eqn through transect pts y = const + grad*x

	float cc;

	xp = xs;
	yp = ys;
	ix = int(xp); iy = int(yp);
	if (ix < 0) ix = 0; if (iy < 0) iy = 0;
	if (ix >= nCols) ix = nCols - 1; if (iy >= nRows) iy = nRows - 1;

	glBegin(GL_LINE_STRIP);

	xx = 500.0f + float(xp - 500.0)*cosCor[int(yp)];

	glVertex2f(xx, yp);

	for (i = 0; i<LONG_TRACE; ++i)   // note this can be optimized
	{
		ix = int(xp);
		iy = int(yp);

		if (ix < 0) exit = true;
		if (iy < 0)exit = true;
		if (ix >= nCols) exit = true;
		if (iy >= nRows) exit = true;

		// NOW do the Transect Test
		v = yp - grad*xp;
		if (v > constant) tSign = 1.0;
		else tSign = -1.0;
		if (otSign*tSign < 0.0) {
			//cerr << "Sign Change " << i << "\n";
			// note the i> 10 is to allow for full loops.
			if (xp > transects_S[t2].tx1 && xp < transects_S[t2].tx2 && i>10) {
				success = exit = true;
			}
		}
		otSign = tSign;
		if (exit) break;

		cc = cosCor[iy];
		ii = iy*nCols + ix;

		dx1 = sign*smFlowX[ii] / (cc*(smSpeed[ii] + 0.01));
		dy1 = sign*smFlowY[ii] / (smSpeed[ii] + 0.01);

		xp2 = xp + dx1*Tscl;
		yp2 = yp + dy1*Tscl;

		ix = int(xp2);
		iy = int(yp2);  if (iy > 699) iy = 699;
		cc = cosCor[iy];
		ii = iy*nCols + ix;

		if (ii < 0) ii = 10;
		dx2 = sign*smFlowX[ii] / (cc*(smSpeed[ii] + 0.01));
		dy2 = sign*smFlowY[ii] / (smSpeed[ii] + 0.01);
		xp = xp + 0.5*(dx2 + dx1)*Tscl;
		yp = yp + 0.5*(dy1 + dy2)*Tscl;

		xx = 500.0f + float(xp - 500.0)*cosCor[int(yp)];

		glVertex2f(xx, yp);


		//xt[i + 1] = 500.0f + float(xt[i] - 500)*cosCor[int(yt[i])];
		//yt[i + 1] = yp;
	}

	glEnd();


	return success;
}

// find the point at which a path intersects a transect
// Used to calculate flow volumes
void Layer:: intersectTransect(int ip, int t2)
{
	int i, len;
	float t2dx, t2dy, t2Len; 
	float xp, yp, v, grad, constant;

	float tSign, otSign;  // the sign of a point relative to the transect
	bool exit = false;
	bool success;
	tSign = otSign = 1.0;
	t2dx = transects_S[t2].tx2 - transects_S[t2].tx1;
	t2dy = transects_S[t2].ty2 - transects_S[t2].ty1;
	grad = t2dy / t2dx; // doesnt work for vertical transects_S
	constant = transects_S[t2].ty1 - grad*transects_S[t2].tx1;

	for (i = 0; i < longTraceLen; ++i)   // note this can be optimized
	{
		xp = flowPaths[ip]->xPath[i];
		yp = flowPaths[ip]->yPath[i];

		v = yp - grad*xp;
		if (v > constant) tSign = 1.0;
		else tSign = -1.0;

		if (otSign*tSign < 0.0) {
			//cerr << "Sign Change " << i << "\n";
			// note the i> 10 is to allow for full loops.
			if (xp > transects_S[t2].tx1 && xp < transects_S[t2].tx2 && i>10) {
				success = exit = true;
			}
		}
		otSign = tSign;


		if (exit) break;

		// note the i> 10 is to allow for full loops.
	}
	otSign = tSign;

	cerr << " Transect Intersection X Y " << xp << " " << yp << " \n";
}

void Layer::getTransectAnchors(int i, float &x1, float &y1, float &x2, float &y2)
{
	float xx;
	xx = 500.0f + (transects_S[i].tx1 - 500.0f)*cosCor[int(transects_S[i].ty1)];
	x1 = xx;// transects_S[i].tx1;
	y1 = transects_S[i].ty1;
	xx = 500.0f + (transects_S[i].tx2 - 500.0f)*cosCor[int(transects_S[i].ty2)];
	x2 = xx; // transects_S[i].tx2;
	y2 = transects_S[i].ty2;

}




