#include "pathlets.h"
#include <iostream>
#include <GL/glut.h>
#include <math.h>
#include "flowPath.h"
//#include "GlyphArray.h"

using namespace std;

// SCHEME_1  SCHEME_2


//#define STRLEN	30  // was 75
//#define MID 15  // this is the backward trace part.

#define STRLEN	80  // was 75
#define MID 40  // this is the backward trace part.
#define PLET 35

#define SAMPLE_SPACE 7 //10

#define MAXAGE	80 // was 10

pathlets::pathlets(int nr, int nc, int level, int colScm)
{
	int i;
	nCols = nc;
	nRows = nr;

	colorScheme = colScm;
	cerr << "TRACE R C " << nRows << " " << nCols << "\n";

	thisLayer = level;

	cosCor = new float[nRows];
	for (i = 0; i < nRows; ++i)
		cosCor[i] = cos((double(i) / 10.0)*3.141592 / 180.0);

	genPts(SAMPLE_SPACE, level);
	traceScale = 18.0f;// 12.0;

	traceScale_2 = 0.5f;// 0.2f;
	traceLen = LONG_TRACE;

	//genPlets();

}



void pathlets::drawAnimatedTraces(float lineWid) // A COMPLET SET OF LINE TRACES
{
	int i,j,k;
	int a;
	int T;  // the tail
	float xpos;
	float *tx, *ty;
	float f;
	unsigned char r, g, b, pc;

	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);

	glEnable(GL_LINE_SMOOTH);

	glLineWidth(lineWid);

	glPushMatrix();	// put it in geographic coords 

	for(i=0;i<nStreakTraces;++i)
	{
		++streakletAge[i];
		++streakletColorAge[i];
		if (streakletAge[i] > MAXAGE)  // we restart if the particle
		{			
			//trace(Xpaths[i],Ypaths[i],Xstart[i], Ystart[i],timer);
			//traceSide(Xpaths[i],Ypaths[i],sideX[i],sideY[i],Xstart[i], Ystart[i],timer);
			Tails[i] = 0;
			streakletAge[i] = 0;
		}
		

		T = Tails[i]%(STRLEN-PLET);

		tx = Xpaths[i];
		ty = Ypaths[i];		

		glBegin(GL_LINE_STRIP);
		a = streakletAge[i];
		f = fade[a]*0.75f;

		r = traceRed[i]; g = traceGreen[i]; b = traceBlue[i];

		for (j = T;j< T+PLET;++j )   // note this can be optimized
		{				
			//pc = 0.1+ float(j-T)/PLET;  // what should this be?
			//glColor4f(1.0,0.3,0.3,pc*f);//0.4f);
			pc = 255 * (j-T) / PLET;  // what should this be?

			glColor4ub(r, g, b, pc);
			//glColor4f(1.0, 0.3, 0.3, pc);//0.4f);

			//k=j%STRLEN; // dont need this?
			k = j;

			xpos = 500.0f + float(tx[k] - 500)*cosCor[int(ty[k])];

			glVertex3f(xpos , ty[k] ,0.0);
		}	
		glEnd();	
		
		float px,py;
		px = tx[k]; py = ty[k]; // head pos

		
		++Tails[i];
		
	}
	glPopMatrix();

}



void pathlets::drawTracesStatic() // A COMPLET SET OF LINE TRACES
{
	int i,j;
	float *tx, *ty;
	float xpos;
	float strlen;
	unsigned char r, g, b,pc;

	strlen = STRLEN;

	glDisable(GL_TEXTURE_2D);

	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);

	glEnable(GL_LINE_SMOOTH);

	glLineWidth(1.3f);

	glPushMatrix();	// put it in geographic coords 

	for(i=0;i<nStreakTraces;++i)
	{
		tx = Xpaths[i];
		ty = Ypaths[i];	
		r = traceRed[i]; g = traceGreen[i];b = traceBlue[i];
		//r = g = 255; b = 0;
		pc = 255;

		glBegin(GL_LINE_STRIP);

		for (j = 0;j< STRLEN;++j )   // note this can be optimized
		{				
			pc =  255*j/ STRLEN;  // what should this be?

			//pcf = float(j) / STRLEN;
			xpos = 500.0f + float(tx[j] - 500)*cosCor[int(ty[j])];
			//glColor4f(pc*0.9f,0.9*pc,0.9*pc,pc);//0.4f)

			glColor4ub(r,g,b, pc);

			glVertex3f(xpos , ty[j] ,0.0f);
		}	
		glEnd();	
		
	}
	glPopMatrix();

}


void pathlets::drawLongTrace()
{
	int i;
	glLineWidth(3.0f);
	glColor3f(0.0f, 0.99f, 0.0f);
	glBegin(GL_LINE_STRIP);

	for (i = 0; i < traceLen; ++i)   // note this can be optimized
	{
		//xpos = 500.0f + float(trLongX[i] - 500)*cosCor[int(trLongY[i])];
		glVertex2f(trLongX[i], trLongY[i]);
	}
	glEnd(); 
	glLineWidth(1.3f);
}
/*
void pathlets::drawArrowHead(int si, float *xp, float *yp, float wid)
{
	float vx, vy;
	float x1, x2;
	float dx, dy;
	float ww;
	ww = 0.3 + wid*0.05;

	//x1 = 500.0f + float(xp[si] - 500)*cosCor[int(yp[si])];
	//x2 = 500.0f + float(xp[si+8] - 500)*cosCor[int(yp[si+8])];
	x1 = xp[si];
	x2 = xp[si+8];

	dx = x2-x1;
	dy = yp[si+8] - yp[si];
	glBegin(GL_TRIANGLES);
	vx = x1 + dy*ww;
	vy = yp[si] - dx*ww;
	glVertex2f(vx, vy);
	vx = x1 - dy*ww;
	vy = yp[si] + dx*ww;
	glVertex2f(vx, vy);
	vy = yp[si+8];
	glVertex2f(x2, vy);

	glEnd();


}
*/
/*
void pathlets::drawPaths(bool top, bool bottom)
{
	int i,ti, tlen;
	float xpos;

	//cerr << "draw currents\n";

	if (bottom)
	{
		//cerr << "bottom\n";
		glColor3f(0.2f, 0.2f, 0.8f);
		for (ti = nPaths_S; ti < nPaths_S + nPaths_D; ++ti)
		{
			glLineWidth(pathWids[ti]);
			tlen = pathLens[ti];

			glBegin(GL_LINE_STRIP);

			for (i = 0; i < tlen; i = i + 3)   // note this can be optimized
			{
				//xpos = 500.0f + float(pathsX[ti][i] - 500)*cosCor[int(pathsY[ti][i])];
				//glVertex2f(xpos, pathsY[ti][i]);
				glVertex2f(pathsX[ti][i], pathsY[ti][i]);

			}

			glEnd();

			for (i = 0; i < tlen; ++i)
				if (i % 100 == 49 && i < (tlen - 20)) drawArrowHead(i, pathsX[ti], pathsY[ti], pathWids[ti]);
		}
	}

	if (top)
	{
		//tcerr << "top\n";
		glColor3f(1.0f, 0.0f, 0.0f);

		for (ti = 0; ti < nPaths_S; ++ti)
		{ 
			if(ti > 1 && ti < 5)
				flowPaths[ti]->drawPath(true);
			else
				flowPaths[ti]->drawPath(false);
		}
		for (ti = 0; ti < nPaths_S; ++ti)
		{
			tlen = pathLens[ti];

			for (i = 0; i < tlen; ++i)
				if (i % 100 == 49 && i < (tlen - 20)) drawArrowHead(i, pathsX[ti], pathsY[ti], pathWids[ti]);
		}
		

	}

	glLineWidth(1.3f);

}

*/

void pathlets::makeTraces()  // add data
{
	int i;

	for(i=0;i<nStreakTraces;++i)
	{
		//traceSide(Xpaths[i],Ypaths[i], sideX[i],sideY[i],Xstart[i],Ystart[i],0);
		trace(Xpaths[i],Ypaths[i], Xstart[i],Ystart[i]);
		streakletAge[i] = rand()%MAXAGE;
		streakletColorAge[i] = rand() % (MAXAGE ); //*2
		Tails[i] = 0;
	}
}


void pathlets::longTrace(float xs, float ys, bool reverse) // SD shallow vs deep
{

	longTracer(xs, ys, trLongX, trLongY, reverse); // reverse
//	if (TB == 2)longTracer(xs, ys, trLongX, trLongY, vX_D, vY_D, speed_S, reverse);


}


void pathlets::longTracer(float xs, float ys, float *xt, float *yt, bool reverse)
{
	int i, ix, iy, ii;
	float xp, yp,xp2,yp2;
	float dx1, dy1, dx2, dy2;  // for predictor corrector
	float sign = 1.0f;
	if(reverse) sign = -1.0f;
	float Tscl = 5.0;
	float cc;

	//	TRACE FORWARD
	// need to interpolate at some point
	xp = xs;
	yp = ys;
	ix = int(xp); iy = int(yp);
	if (ix < 0) ix = 0; if (iy < 0) iy = 0;
	if (ix >= nCols) ix = nCols - 1; if (iy >= nRows) iy = nRows - 1;

	for (i = 0; i<LONG_TRACE; ++i)   // note this can be optimized
	{
		ix = int(xp);
		iy = int(yp);

		if (ix < 3) ix = 3;
		if (iy < 4) iy = 4;
		if (ix >= nCols) ix = nCols - 1;
		if (iy >= nRows) iy = nRows - 1;

		cc = cosCor[iy];
		ii = iy*nCols + ix;

		dx1 = sign*vX_S[ii] / (cc*(speed_S[ii] + 0.01));
		dy1 = sign*vY_S[ii]/(speed_S[ii]+0.01);


		xp2 = xp + dx1*Tscl;
		yp2 = yp + dy1*Tscl;

		ix = int(xp2);
		iy = int(yp2);  if (iy > 699) iy = 699;
		cc = cosCor[iy];
		ii = iy*nCols + ix;

		if (ii < 0) ii = 10;
		dx2 = sign*vX_S[ii] / (cc*(speed_S[ii] + 0.01));
		dy2 = sign*vY_S[ii]/(speed_S[ii] + 0.01);

		xp = xp + 0.5*(dx2 + dx1)*Tscl;
		yp = yp + 0.5*(dy1 + dy2)*Tscl;
		xt[i] = xp;
		yt[i] = yp;
	}

	// now fix the x values
	for (i = 0; i<LONG_TRACE; ++i)   
		xt[i]= 500.0f + float(xt[i] - 500)*cosCor[int(yt[i])];
}

void pathlets::trace(float *tx, float *ty,  float xs, float ys)
{  // we compute a new path for streak s at time t;
		int i, ix,iy, ii;
		float xp, yp;
		float dx1,dy1;	

		float cc;	 
		//	TRACE FORWARD
		// need to interpolate at some point
		xp = xs;
		yp = ys;
		ix = int(xp); iy = int(yp);
		if (ix < 0) ix = 0; if (iy < 10) iy =10;
		if (ix >= nCols) ix = nCols-1; if (iy >= nRows) iy = nRows-1;

		for (i = MID; i < STRLEN; ++i)   // note this can be optimized
		{
			ii = 0;
			ix = int(xp);
			if (yp < 10.0f) yp = 10.0f;
			iy = int(yp);

			if (ix < 0) ix = 0;
			if (ix >= nCols-1) ix = nCols - 1;
			if (iy >= nRows - 4) {  // stop tracing
				dx1 = dy1 = 0.0f;
			}
			else {
				cc = cosCor[iy];
				ii = iy*nCols + ix;
				dx1 = vX_S[ii]/cc;
				dy1 = vY_S[ii];
				//dx1 = vX_S[ii] /( cc * speed_S[ii] +0.0001);
				//dy1 = vY_S[ii]/speed_S[ii]+0.0001;
			}

			//xp = xp + dx1*traceScale;
			//yp = yp + dy1*traceScale; 

			xp = xp + dx1*traceScale + dx1*traceScale_2/ (speed_S[ii] + 0.0001);
			yp = yp + dy1*traceScale + dy1*traceScale_2 / (speed_S[ii] + 0.0001);

			
			tx[i] = xp;
			ty[i] = yp;
		}

		xp = xs; yp = ys;
		ix = int(xp); iy = int(yp);
		if (ix < 0) ix = 0;
		if (iy < 10) iy = 10;
		if (ix >= nCols-1) ix = nCols-1;
		if (iy >= nRows-4) iy = nRows-4;

		for (i = MID;i> -1;--i)   // note this can be optimized
		{	

			ix = int(xp);
			if (yp < 10.0f) yp = 10.0f;
			iy = int(yp);

			if (ix < 0) ix = 0;
			if (ix >= nCols) ix = nCols-1;
			if (iy >= nRows-5) iy = nRows-6;

			cc = cosCor[iy];

			ii = iy*nCols + ix;
		    dx1 = vX_S[ii]/cc;
			dy1 = vY_S[ii];
//
	//		dx1 = vX_S[ii] / (cc * speed_S[ii] + 0.0001);
	//		dy1 = vY_S[ii] / speed_S[ii] + 0.0001;

				///dx1 = fData->getUVelocity(time1,rParams->PressureLayer,iy,ix)*traceScale;
		        ///dy1 = fData->getVVelocity(time1,rParams->PressureLayer,iy,ix)*traceScale;
			//xp = xp - dx1*traceScale;
			//yp = yp - dy1*traceScale; //HACK


			xp = xp - dx1*traceScale - dx1*traceScale_2 / (speed_S[ii] + 0.0001);
			yp = yp - dy1*traceScale - dy1*traceScale_2 / (speed_S[ii] + 0.0001);

				
			tx[i] = xp;
			ty[i] = yp;
		}
	
}

// backup the start positions
void pathlets::backupStart(int i)  
{
	float xp, yp, cc, dx1, dy1;
	int  ix, iy, j,  ii;

	xp = Xstart[i];
	yp = Ystart[i];

	streakletAge[i] = 0;

	for (j = 0; j < PLET; ++j)
	{
		if (yp < 0.0f) yp = 0.0f;	
		if (yp > 699.0f) yp = 699.0f;

		ix = int(xp);
		iy = int(yp);
		cc = cosCor[iy];
		ii = iy*nCols + ix;
		dx1 = vX_S[ii] / cc;
		dy1 = vY_S[ii];

		xp = xp - dx1*traceScale - dx1*traceScale_2 / (speed_S[ii] + 0.0001);
		yp = yp - dy1*traceScale - dy1*traceScale_2 / (speed_S[ii] + 0.0001);

	}


	//float xpos = 500.0f + float(xp - 500)*cosCor[int(yp)];
	//	glColor3f(0.0, 1.0, 1.0);
	//glRectf(xpos - 5.0f, yp - 5.0f, xpos + 5.0f, yp + 5.0f);
	newXstart[i] = xp;
	newYstart[i] = yp;

}




// start with a simple loop

void pathlets::genPlets()
{
	int i;
	Xplets = new float*[nStreakTraces];
	Yplets = new float*[nStreakTraces];
	Heads = new int[nStreakTraces];
	//Tails = new int[nStreakTraces];
	for (i = 0; i < nStreakTraces; ++i)
	{
		Xplets[i] = new float[PLET];
		Yplets[i] = new float[PLET];
	}

	for (i = 0; i < nStreakTraces; ++i)backupStart(i);

	for (i = 0; i < nStreakTraces; ++i)
	{
		Heads[i] = rand()% PLET;
		Xplets[i][Heads[i]] = newXstart[i];
		Yplets[i][Heads[i]] = newYstart[i];
		isActive[i] = true;

	}
	// initialize them all
	for (i = 0; i < PLET; ++i)updateAllPlets();  

	cerr << "\n***********************************PLETS************************** \n";

	cerr << Xplets[0][Heads[0]] << " " << Yplets[0][Heads[0]] << "\n";
	int k;
	k = Heads[0];
/*
	for (i = 0; i < 10; ++i)
	{ 
		
		cerr << k << " " << Xplets[0][k] << " " << Yplets[0][k] << "\n";
		k = (k + 1) % PLET;
	}
	*/
}

void pathlets::drawPlets(float Lwid)
{
	int i, j, k, ii,a;
	float xpos,*tx, *ty;
	float f;
	unsigned char r,g,b,pc;
	int trans;

	if (Lwid < 3.0f)trans = 160;  // a hack to set the transparency 
	else trans = 255;

	glLineWidth(Lwid);
	for (i = 0; i < nStreakTraces; ++i)
	{
		if (isActive[i]) {
			k = Heads[i];
			tx = Xplets[i];
			ty = Yplets[i];

			a = streakletAge[i];
			f = fade[a] * 0.75f;  // not used

			r = traceRed[i]; g = traceGreen[i]; b = traceBlue[i];

			glBegin(GL_LINE_STRIP);
			for (j = 1; j <= PLET; ++j) {

				ii = (j + k) % PLET;
				xpos = 500.0f + float(tx[ii] - 500)*cosCor[int(ty[ii])];

				pc = trans * j / PLET;  // what should this be?

				glColor4ub(r, g, b, pc);

				glVertex2f(xpos, ty[ii]);

			}
			glEnd();

			++streakletAge[i];
			if (streakletAge[i] >= MAXAGE)  // we restart if the particle
			{
				restart(i);  //WORK_1
				//isActive[i] = false;//only do this if advancing in time
			}
		}

	}
	updateAllPlets();
}

void pathlets::restart(int i)  // trace a single plet to get it restarted
{
	float xp, yp, cc,dx1, dy1;
	int  ix, iy, j, k, ii,h;

	h = Heads[i];
	backupStart(i);

	Xplets[i][h] = newXstart[i];
	Yplets[i][h] = newYstart[i];
	xp = Xplets[i][h];
	yp = Yplets[i][h];	

	streakletAge[i] = 0;

	for (j = 0; j < PLET; ++j)
	{
		if (yp < 0.0f) yp = 0.0f;
		if (yp > 690.0f) yp = 690.0f;
		if (xp < 0.0f) xp = 0.0f;
		k = (Heads[i] + 1) % PLET;  // next position
		Heads[i] = k;

		ix = int(xp);
		iy = int(yp);		
		cc = cosCor[iy];
		ii = iy*nCols + ix;
		dx1 = vX_S[ii] / cc;
		dy1 = vY_S[ii];

		//xp = xp + dx1*traceScale;
		//yp = yp + dy1*traceScale;
		xp = xp + dx1*traceScale + dx1*traceScale_2 / (speed_S[ii] + 0.0001);
		yp = yp + dy1*traceScale + dy1*traceScale_2 / (speed_S[ii] + 0.0001);
		Xplets[i][k] = xp;
		Yplets[i][k] = yp;
	}

}

void pathlets::updateAllPlets()  // updates them all used in rendering
{
	int i,ix,iy, j,k,ii;
	float xp, yp,cc, dx1,dy1;
	for (i = 0; i < nStreakTraces; ++i)
	{
		if (isActive[i]) {
			xp = Xplets[i][Heads[i]];
			yp = Yplets[i][Heads[i]];

			k = (Heads[i] + 1) % PLET;  // next position
			Heads[i] = k;

			if (xp > float(nCols - 4)) xp = nCols - 4;
			if (yp < 10.0f) yp = 10.0f;

			if (yp > 690.0f) yp = 690.0f;
			ix = int(xp);
			iy = int(yp);

			cc = cosCor[iy];
			ii = iy*nCols + ix;
			dx1 = vX_S[ii] / cc;
			dy1 = vY_S[ii];
			//xp = xp + dx1*traceScale + dx1*traceScale_2 / (speed_S[ii] + 0.0001);

			xp = xp + dx1*traceScale + dx1*traceScale_2 / (speed_S[ii] + 0.0001);
			yp = yp + dy1*traceScale + dy1*traceScale_2 / (speed_S[ii] + 0.0001);
			Xplets[i][k] = xp;
			Yplets[i][k] = yp;

			//Xplets[i][k] = xp + dx1*traceScale;
			//Yplets[i][k] = yp + dy1*traceScale;
		}
	}

}

void pathlets::genPts(int sample, int L)  // sample density
{
	int i,j;
	int count;
	int border = 40; // border with no samples
	int right = nCols - 30;
	int top = nRows - 30;
	float jfac,hjfac;
	count = 0;
	int xsample;

	thisLayer = L;

	for (j= border;j< right; j +=sample)
		for (i=border;i<top; i = i+sample)
		{		
					++count;
		}
		cerr << "COUNT " << count << "\n";

	trLongX = new float[LONG_TRACE];
	trLongY = new float[LONG_TRACE];

	nStreakTraces = count;
	cerr << "NUMBER OF TRACES " << nStreakTraces << "\n";

	Xstart = new float[nStreakTraces];
	Ystart = new float[nStreakTraces];
	isActive = new bool[nStreakTraces];

	traceRed = new unsigned char[nStreakTraces];
	traceGreen = new unsigned char[nStreakTraces];
	traceBlue = new unsigned char[nStreakTraces];
	// default colors
	if(L==1)
		for (i = 0; i < nStreakTraces; ++i) {
			if (colorScheme == 1) { traceRed[i] = 0; traceGreen[i] = 00; traceBlue[i] = 20; }//scheme 1
			else { traceRed[i] = 80; traceGreen[i] = 80; traceBlue[i] = 100; }
		}
	if (L ==0)
		for (i = 0; i < nStreakTraces; ++i) {
			if (colorScheme == 1) {
				traceRed[i] = 250; traceGreen[i] = 50; traceBlue[i] = 0;}//scheme 1
			else { traceRed[i] = 140; traceGreen[i] = 80; traceBlue[i] = 80; }
		}

	newXstart = new float[nStreakTraces];  // what is this for
	newYstart = new float[nStreakTraces];

	count = 0;
	jfac = sample*0.75f;
	hjfac = jfac/2.0f;

	cerr << "FACS " << jfac << " " << hjfac << "\n";

	float maxX = -1000.0f;

	for (i=border;i<top; i = i+sample)
	{
		xsample = sample/cos(0.5*3.141592*double(i)/double(nRows));
		//if(i%10 == 0) cerr << "XSAMPLE"
		//xsample = sample;
		if (i % 100 == 0)cerr << "XSAMPLE " << xsample << "\n";
		//if (i > nRows / 2) xsample = sample *3;
		for (j= border;j< right; j +=xsample)
		{		
			if(i < 400 || j<850)  // should use land as a mask
			{ 
			Xstart[count] = j + jfac*float(rand()%1000)/999.0 - hjfac;
			Ystart[count] = i + jfac*float(rand()%1000)/999.0 - hjfac;

			newXstart[count] = Xstart[count];
			newYstart[count] = Ystart[count];				
			
			++ count;
			}
		}
 }
	//cerr << "MAX X " << maxX << "\n";
	// Hand Seed  to fill in some gaps
	Xstart[count] = 450; Ystart[count] = 40;++count;
	Xstart[count] = 455; Ystart[count] = 40; ++count;
	Xstart[count] = 460; Ystart[count] = 40; ++count;
	Xstart[count] = 465; Ystart[count] = 40; ++count;
	Xstart[count] = 450; Ystart[count] = 45; ++count;
	Xstart[count] = 455; Ystart[count] = 45; ++count;
	Xstart[count] = 460; Ystart[count] = 45; ++count;
	Xstart[count] = 465; Ystart[count] = 45; ++count;
	Xstart[count] = 450; Ystart[count] = 43; ++count;
	Xstart[count] = 455; Ystart[count] = 43; ++count;
	Xstart[count] = 460; Ystart[count] = 43; ++count;
	Xstart[count] = 465; Ystart[count] = 43; ++count;

	Xstart[count] = 128; Ystart[count] = 239; ++count;
	Xstart[count] = 130; Ystart[count] = 239; ++count;
	Xstart[count] = 132; Ystart[count] = 239; ++count;
	Xstart[count] = 134; Ystart[count] = 239; ++count;


	RNDcols = new float[nStreakTraces];
	Xpaths = new float*[nStreakTraces];
	Ypaths = new float*[nStreakTraces];

	Tails = new int[nStreakTraces];
	for(i=0;i<nStreakTraces;++i)
	{
		Xpaths[i] = new float[STRLEN];
		Ypaths[i] = new float[STRLEN];
		//sideX[i] = new float[STRLEN];
		//sideY[i] = new float[STRLEN];
	}	
	cerr << "GEN XStart " << Xstart[0] << "\n";

	streakletAge = new int[nStreakTraces];// was NPTS
	for(i=0;i<nStreakTraces;++i)  // initialize the ages randomly
	{
		streakletAge[i] = rand()%MAXAGE;
		RNDcols[i] = 0.5+0.5*float(rand()%100)/99.0f;
	}

	streakletColorAge = new int[nStreakTraces];// was NPTS

	fade = new float[MAXAGE];
	for(i=0;i<MAXAGE;++i) 
	{
		fade[i] = 1.0;
	}

	fade[MAXAGE-1] = 0.25f;
	fade[MAXAGE-2] = 0.5f;
	fade[MAXAGE-3] = 0.75f;  // fade it out
	fade[0] = 0.1f;  // fade in
	fade[1] = 0.3f;
	fade[2] = 0.5f;
	fade[3] = 0.7f;
	fade[4] = 0.9f;

}



void pathlets::setPtColors( unsigned char *mask,bool first)  // sample density vased on screen mask
{
	int i, j, ix,iy, indx, iv;
	int count;
	int border = 30; // border with no samples
	int right = nCols - 30;
	int top = nRows - 30;
	count = 0;

	cerr << "R C " << nRows << " " << nCols << "\n";
		cerr << "NUMBER OF TRACES " << nStreakTraces << "\n";

	count = 0;

	for (i = 0; i<nStreakTraces; ++i)  // 
	{
		
		// find the seed point in pixels
		iy = Ystart[i];
		ix = 500.0f + (Xstart[i] - 500.0f)*cosCor[iy];
		indx = (iy*nCols + ix) * 3;

		//iv = mask[indx+2];
		iv = mask[indx]+mask[indx+1];

		if (first&&iv > 0)  //initialize for frame 1
		{
			traceRed[i] = mask[indx];
			traceGreen[i] = mask[indx + 1];
			traceBlue[i] = mask[indx + 2];

		}


		if (streakletAge[i] == 1)
		{
			if (iv > 0) // color if not black
			{
				traceRed[i] = mask[indx];
				traceGreen[i] = mask[indx + 1];
				traceBlue[i] = mask[indx + 2];
				//isActive[i] = true;
			}
			else { // set the default colors

				if (streakletColorAge[i] < MAXAGE * 3 / 4) //*2
				{


					if (thisLayer == 1) {
						if (colorScheme == 1) { traceRed[i] = 0; traceGreen[i] = 00; traceBlue[i] = 20; }//scheme 1
						else { traceRed[i] = 80; traceGreen[i] = 80; traceBlue[i] = 100; }
					}
					if (thisLayer == 0) {
						if (colorScheme == 1) {
							traceRed[i] = 250; traceGreen[i] = 50; traceBlue[i] = 0;
						}//scheme 1
						else { traceRed[i] = 180; traceGreen[i] = 80; traceBlue[i] = 80; }
					}

					streakletColorAge[i] = 0;
				}
			
			}
		}

			/*
			if (iv > 0) // color if not black
			{
				traceRed[i] = mask[indx];
				traceGreen[i] = mask[indx + 1];
				traceBlue[i] = mask[indx + 2];
				//isActive[i] = true;
			}
			else {
				if (streakletColorAge[i] == (MAXAGE - 2)) //*2
				{

					if (thisLayer == 1) {
						if (colorScheme == 1) { traceRed[i] = 0; traceGreen[i] = 00; traceBlue[i] = 20; }//scheme 1
						else { traceRed[i] = 80; traceGreen[i] = 80; traceBlue[i] = 100; }
					}
					if (thisLayer == 0) {
						if (colorScheme == 1) {
							traceRed[i] = 250; traceGreen[i] = 50; traceBlue[i] = 0;
						}//scheme 1
						else { traceRed[i] = 180; traceGreen[i] = 80; traceBlue[i] = 80; }
					}
				}
			}
			*/
		}

}




