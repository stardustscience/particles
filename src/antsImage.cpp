#include <GL/glut.h>
#include "antsImage.h"

antsImage::antsImage()
{

	current = new unsigned char[1000 * 700 * 3];
	ants = new unsigned char[1000 * 700 * 4];

}

void antsImage::capture()
{
	glReadPixels(140, 320,1000, 700, GL_RGB, GL_UNSIGNED_BYTE, current);
}

void antsImage::captureAnts()
{
	glReadPixels(140, 320, 1000, 700, GL_RGBA, GL_UNSIGNED_BYTE, ants);
}

void antsImage::drawPrev()
{
	glRasterPos2i(140, 320);

	glDrawPixels(1000, 700, GL_RGB, GL_UNSIGNED_BYTE, current);
}

void antsImage::drawTransparentAnts()
{
	int i;
	for (i = 0; i < 1000 * 700;++i) // make transparent
		ants[i * 4 + 3] = 30;
	glRasterPos2i(140, 320);

	glDrawPixels(1000, 700, GL_RGBA, GL_UNSIGNED_BYTE, ants);
}



