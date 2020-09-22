#include <GL/glut.h>

#include "flowBar.h"

#define BAR_WID 20.0f
#define MAX_HT 105.0f  // equivalent to 100 Sv

flowBar:: flowBar(float xp, float yp, int ci, int between, float mn)
{
	int i;
	xPos = xp;
	yPos = yp;
	left = xp + 3.5f;
	right = xp + BAR_WID - 3.0f;

	if (ci == 1) { red = 1.0f, green = 0.3f, blue = 0.3f; }
	if (ci == 2) { red = 0.35f, green = 0.5f, blue = 0.95f; }

	anchorX = anchorY = 100.0f;

	for (i = 0; i < MAX_STEPS; ++i) history[i] = 0.0f;
	nsteps = 0;
	mean = mn;

	betweens = between;

}

void flowBar::draw(float flow, int RL,  int ts)
{
	int i;
	float offset, x,y,v;
	if (RL == 1) offset = BAR_WID;
	else offset = 0;

	nsteps = ts / 4;

	history[ts/4] = flow;

	glLineWidth(1.0f);

	glColor3f(0.5f, 0.5f, 0.5f);
	glRectf(xPos, yPos, xPos + BAR_WID, yPos + MAX_HT);

	glColor3f(0.85f, 0.85f, 0.85f);
	glBegin(GL_LINES);
	for (i = 0; i < 11; ++i)
	{
		glVertex2f(xPos, yPos + i*10);
		glVertex2f(xPos+BAR_WID, yPos + i*10);
	}
	glEnd();

	glColor3f(red, green, blue);
	glRectf(left, yPos, right, yPos + flow);

	glColor3f(0.0f, 0.55f, 0.5f);
	glBegin(GL_LINES);	
		glVertex2f(xPos + offset, yPos );
		glVertex2f(anchorX, anchorY);
	glEnd();

	glColor4f(red, green, blue, 0.55f);

	// draw the History Graph

	glBegin(GL_TRIANGLE_STRIP);
	for (i = 0; i <= nsteps; ++i) {
		x = xPos  - i;
		y = yPos + history[nsteps - i];
		glVertex2f(x, y);
		glVertex2f(x, yPos);
	}
	glEnd();

	// anomalies plot

	glColor4f(0.8, 0.8, 0.8, 1.0);  // background for anomalies
	glRectf(xPos, yPos - 32.0f, xPos - nsteps, yPos);

	glBegin(GL_TRIANGLE_STRIP);
	for (i = 0; i <= nsteps; ++i) {
		x = xPos - i;
		v = history[nsteps - i] - mean;
		if (v < 0.0f)glColor3f(0.3, 0.3, 1.0);
		else glColor3f(1.0, 0.2, 0.2);
		y = yPos - 16.0f + v*1.4f;
		glVertex2f(x, y);
		glVertex2f(x, yPos-16.0f);
	}
	glEnd();

	int yrs;
	float start,gap;
	start = xPos - nsteps;
	gap = float(12.0f*float(betweens)/4.0f);  // FIX this  12*betwweens/4

	glColor4f(red*0.4, green*0.4, blue*0.4, 0.8f);

	glBegin(GL_LINES);
	yrs = ts /(12*betweens); // = 12*BETWEENS
	for (i = 0; i <= yrs; ++i) {
		x = start + i*gap;
		glVertex2f(x, yPos-5.0f);
		glVertex2f(x, yPos+20.0f);

	}
	glEnd();

}

float flowBar::getMean()
{
	int i;

	float sum=0.0f;
	for (i = 0; i <= nsteps; ++i)
	{
		sum += history[i];
	}

	return sum / float(nsteps);
	
}