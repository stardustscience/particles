//class flowPath
//{
//public:
//	flowPath(float *xp, float *yp, int len, float flow);
//	makePath();
//	drawPath();
//private:
//	float flowRate;
//	float *xPath, *yPath;
//	float *leftX, *leftY, *rightX, rightY;
//	int Len;
//	float red, green, blue;
//};
#include <math.h> 
#include <GL/glut.h>
#include "flowPath.h"


flowPath::flowPath(float *xp, float *yp, int len, float vol, float scl)
{
	xPath = xp;
	yPath = yp;
	pLen = len;
	leftX = new float[pLen];
	rightX = new float[pLen];
	leftY = new float[pLen];
	rightY = new float[pLen];
	tmp = new float[pLen];

	pathScale = scl;
	wid = 0.3f*vol;

}

void flowPath::smooth(int sm)
{
	int i, k, ii, sm2;
	float avg;

	sm2 = sm / 2;
	for (i = sm2; i < pLen - sm2; ++i)
	{
		avg = 0;
		for (k = 0; k < sm; ++k)
		{
			ii = (i - sm2 + k) % pLen;
			avg = avg + xPath[ii];
		}
		tmp[i%pLen] = avg / float(sm);
	}
	for (i = sm2; i < pLen - sm2; ++i)
		xPath[i] = tmp[i];

	for (i = sm2; i < pLen + sm2; ++i)
	{
		avg = 0;
		for (k = 0; k < sm; ++k)
		{
			ii = (i - sm2 + k) % pLen;
			avg = avg + yPath[ii];
		}
		tmp[i%pLen] = avg / float(sm);
	}
		for (i = sm2; i < pLen - sm2; ++i)
		yPath[i] = tmp[i];

}

/*
void flowPath::smooth(int sm)
{
	int i, k,ii,sm2;
	float avg;

	sm2 = sm / 2;
	for (i = sm2; i < pLen + sm2; ++i)
	{
		avg = 0;
		for (k = 0; k < sm; ++k)
		{
			ii = (i - sm2 + k)%pLen;
			avg = avg + xPath[ii];
		}
		tmp[i%pLen] = avg / float(sm);
	}
	for (i = 0; i < pLen; ++i)
		xPath[i] = tmp[i];

	for (i = sm2; i < pLen + sm2; ++i)
	{
		avg = 0;
		for (k = 0; k < sm; ++k)
		{
			ii = (i - sm2 + k) % pLen;
			avg = avg + yPath[ii];
		}
		tmp[i%pLen] = avg / float(sm);
	}
	for (i = 0; i < pLen; ++i)
		yPath[i] = tmp[i];
	
}
*/

void flowPath::makePath()
{
	int i;
	float sLn;
	float dx, dy;
	float wSC = 1.0;
	for (i = 0; i < pLen-1; ++i)
	{
		wSC = 1.0 - pathScale*float(i)/ float(pLen);
		dx = xPath[i + 1] - xPath[i];
		dy = yPath[i + 1] - yPath[i];
		sLn = sqrt(dx*dx + dy*dy);

		leftX[i] = xPath[i] + wSC*wid*dy / sLn;
		leftY[i] = yPath[i] - wSC*wid*dx / sLn;

		rightX[i] = xPath[i] - wSC*wid*dy / sLn;
		rightY[i] = yPath[i] + wSC*wid*dx / sLn;

	}
	leftX[pLen -1] = xPath[pLen -1] + wSC*wid*dy / sLn;
	leftY[pLen - 1] = yPath[pLen - 1] - wSC*wid*dx / sLn;

	rightX[pLen - 1] = xPath[pLen - 1] - wSC*wid*dy / sLn;
	rightY[pLen - 1] = yPath[pLen - 1] + wSC*wid*dx / sLn;

}

void flowPath::drawPath(bool close)
{
	int i;
	glBegin(GL_TRIANGLE_STRIP);
	
	for (i = 0; i < pLen - 1; i=i+2)
	{
		glVertex2f(leftX[i], leftY[i]);
		glVertex2f(rightX[i], rightY[i]);

	}

	if (close)
	{
		glVertex2f(leftX[0], leftY[0]);
		glVertex2f(rightX[0], rightY[0]);
	}
	glEnd();

}
