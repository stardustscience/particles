#include <GL/glut.h>
#include <iostream>

#include "arrow.h"
using namespace std;


Arrow::Arrow(float xp, float yp, float dx, float dy, float r, float g, float b, float p)
{
	float xOrth, yOrth;
	x1 = xp + dx;
	y1 = yp + dy;
	xOrth = dy *0.25;
	yOrth = -dx*0.25;

	x2 = xp + xOrth;
	y2 = yp + yOrth;
	x3 = xp - xOrth;
	y3 = yp - yOrth;

	red = r; green = g; blue = b;
	plus = p;

	cerr << "\n"<< x1 << " " << y1 << "\n";
	cerr << x2 << " " << y2 << "\n";
	cerr << x3 << " " << y3 << "\n";
}

void Arrow::draw()
{

	glLineWidth(1.5);
	glColor4f(red, green, blue, 0.5);
	glBegin(GL_TRIANGLES);
	glVertex2f(x1, y1);
	glVertex2f(x2, y2);
	glVertex2f(x3, y3);
	glEnd();

	glColor4f(red+plus, green + plus, blue+plus, 0.75 );

	glBegin(GL_LINE_STRIP);
	glVertex2f(x1, y1);
	glVertex2f(x2, y2);
	glVertex2f(x3, y3);
	glVertex2f(x1, y1);
	glEnd();

//		glRectf(100.0, 100.0, 200.0, 200.0);


}
