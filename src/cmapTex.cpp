#include <GL/glut.h>
#include <iostream>
#include <math.h>
using namespace std;

#include "cmapTex.h"

#define TEX_MAP_ROWS 1024


cmapTex::cmapTex()
{
	basePos = 0;
	cmapCount = 0;
}	

void cmapTex::load(char *fname,bool reverse, int cmapno)
{
	float v;
	int count,rows;
	FILE *cf;

	int i,j;

	//cerr << "LOAD Color Map " << fname << "\n";
	count = 0;
	cf = fopen(fname,"r");
	while(fscanf(cf,"%f",&v)!=EOF)
		++count;

	fclose(cf);
	cf = fopen(fname,"r");
	rows = count/4;
	//cerr << "Count " << rows << "\n";

	if(reverse)
		for(i=0;i<rows;++i)
		{
			fscanf(cf,"%f",&v); vals[i] = v;
			for(j=0;j<3;++j)
				fscanf(cf,"%f",&original_rgb[rows-i-1][j]);// reverse
		}
	else
		for(i=0;i<rows;++i)
		{
			fscanf(cf,"%f",&v); vals[i] = v;
			for(j=0;j<3;++j)
				fscanf(cf,"%f",&original_rgb[i][j]);// reverse
		}
	fclose(cf);

	float v1,v2,dist,p;
	int i1,i2,id;
	float r1,r2,g1,g2,b1,b2;

	for(i=0;i<rows-1;++i)
	{
		v1 = vals[i]; v2 = vals[i+1];
		i1 = v1*1000.0 + 0.5;
		i2 = v2*1000.0 + 0.5;
		id = i2-i1;
		r1 = original_rgb[i][0];r2 = original_rgb[i+1][0];
		g1 = original_rgb[i][1];g2 = original_rgb[i+1][1];
		b1 = original_rgb[i][2];b2 = original_rgb[i+1][2];

		for(j=i1;j<i2;++j)
		{
			p = float(j-i1)/float(id);
			rgb[cmapno][j][0] = (1.0f - p)*r1 + p*r2;
			rgb[cmapno][j][1] = (1.0f - p)*g1 + p*g2;
			rgb[cmapno][j][2] = (1.0f - p)*b1 + p*b2;
		}

		dist= (v2-v1)*1000.0f;
		//cerr << i << " " << vals[i] << " " << dist << "\n";
	}
	for(i=1000;i<MAXCOLS;++i)
	{
		rgb[cmapno][i][0] = original_rgb[rows-1][0]; //HACK for DEEPS
		rgb[cmapno][i][1] = original_rgb[rows-1][1];
		rgb[cmapno][i][2] = original_rgb[rows-1][2];

	}

	++cmapCount;

	//mkTex2Dcmap();

}

void cmapTex::mkTex2Dcmap()// move this out later
{
	// a simple mult-color heat map
	unsigned char *texCmap;
	int i,j,index;
	int n,r,g,b,a; 
	n= TEX_MAP_ROWS*8*4;
	texCmap = new unsigned char[n];

	for (i = 0; i < n; ++i) texCmap[i] = 255;

	cerr << "GEN Colormap\n";

	for(i=0;i<TEX_MAP_ROWS;++i)
		for(j=0;j<2;++j)
		{
			index = i*8*4+ j*4;  // first two slots
			texCmap[index] = int(rgb[0][i][0]*255.5f);
			texCmap[index+1] = int(rgb[0][i][1]*255.5f);
			texCmap[index+2] = int(rgb[0][i][2]*255.5f);
			texCmap[index+3] = 255;
		}
	if(cmapCount > 1 ) 
		for(i=0;i<TEX_MAP_ROWS;++i)
		for(j=2;j<4;++j)
		{
			index = i*8*4+ j*4; // second two slots
			texCmap[index] = int(rgb[1][i][0]*255.5f);
			texCmap[index+1] = int(rgb[1][i][1]*255.5f);
			texCmap[index+2] = int(rgb[1][i][2]*255.5f);
			texCmap[index+3] = 255;
		}
	if (cmapCount > 2)
		for (i = 0; i<TEX_MAP_ROWS; ++i)
			for (j = 4; j<6; ++j)
			{
				index = i * 8 * 4 + j * 4; // second two slots
				texCmap[index] = int(rgb[2][i][0] * 255.5f);
				texCmap[index + 1] = int(rgb[2][i][1] * 255.5f);
				texCmap[index + 2] = int(rgb[2][i][2] * 255.5f);
				texCmap[index + 3] = 255;
			}
	if (cmapCount > 3)
		for (i = 0; i<TEX_MAP_ROWS; ++i)
			for (j = 6; j<8; ++j)
			{
				index = i * 8 * 4 + j * 4; // second two slots
				texCmap[index] = int(rgb[3][i][0] * 255.5f);
				texCmap[index + 1] = int(rgb[3][i][1] * 255.5f);
				texCmap[index + 2] = int(rgb[3][i][2] * 255.5f);
				texCmap[index + 3] = 255;
			}


	glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
	glGenTextures(5, tex2dNames);

	glEnable(GL_TEXTURE_2D);

	//glTexImage3D = (PFNGLTEXIMAGE3DPROC) wglGetProcAddress("glTexImage3D");
	//if(glTexImage3D == NULL) cerr << "GetTexImage3D failed \n";
	glBindTexture(GL_TEXTURE_2D, tex2dNames[1]);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S,  GL_CLAMP);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T,  GL_CLAMP);
	//glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_R, GL_REPEAT);

	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	glTexEnvi(GL_TEXTURE_ENV,GL_TEXTURE_ENV_MODE,GL_MODULATE);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, 8,1024,0, GL_RGBA, GL_UNSIGNED_BYTE, texCmap); 
//glTexImage3D(GL_TEXTURE_3D, 0, GL_RGB8, WIDTH, HEIGHT, DEPTH, 0, GL_RGB,
	glDisable(GL_TEXTURE_2D);

	delete [] texCmap;
}

void cmapTex::drawTexColorMap()
{

	glColor3f(1.0,1.0,1.0);

	glEnable(GL_TEXTURE_2D);
	glBindTexture(GL_TEXTURE_2D, tex2dNames[1]);
	glBegin(GL_TRIANGLE_STRIP);

	glTexCoord2f(0.0,0.0);
//glColor3f(1.0f,1.0f,0.0f);
	glVertex2f(-2.0,-2.0);

	glTexCoord2f(0.0,1.0);
	glVertex2f(-2.0,2.0);

//glColor3f(0.5f,0.5f,1.0f);
glTexCoord2f(1.0,0.0);
	glVertex2f(2.0,-2.0);

	glTexCoord2f(1.0,1.0);
	glVertex2f(2.0,2.0);
	glEnd();
}

void cmapTex::setColorMap()
{

	glEnable(GL_TEXTURE_2D);
	glBindTexture(GL_TEXTURE_2D, tex2dNames[1]);

}

void cmapTex::loadStepped(char *fname)
{
	float v;
	int count,rows;
	FILE *cf;

	int i,j;

	cerr << "LOAD Stepped Color Map " << fname << "\n";
	count = 0;
	cf = fopen(fname,"r");
	while(fscanf(cf,"%f",&v)!=EOF)
		++count;
	fclose(cf);
	cf = fopen(fname,"r");
	rows = count/4;
	cerr << "ColorMap " << rows << "\n";

	for(i=0;i<rows;++i)
	{
		fscanf(cf,"%f",&v); vals[i] = v;
		for(j=0;j<3;++j)
			fscanf(cf,"%f",&original_rgb[i][j]);
	}
	fclose(cf);

	//float r1,g1,b1,
	float bandwid;

	int ind;
	bandwid = float(1000)/rows;

	for(i=0;i<1000;++i)
	{
		ind = int(float(i)/bandwid);

		rgb[0][i][0] = original_rgb[ind][0];
		rgb[0][i][1] = original_rgb[ind][1];
		rgb[0][i][2] = original_rgb[ind][2];

	}
	for(i=1000;i<MAXCOLS;++i)
	{
		rgb[0][i][0] = original_rgb[rows-1][0];
		rgb[0][i][1] = original_rgb[rows-1][1];
		rgb[0][i][2] = original_rgb[rows-1][2];
	}
}

void cmapTex::stepMap(int nSteps)
{
	// turn a loaded colormap into a stepped map.

	int ind,i,ii;
	float bandwid;
	bandwid = float(1000)/nSteps;

	for(i=0;i<1000;++i)
	{
		ind = int(float(i)/bandwid);
		ii = int(ind*(bandwid) + bandwid*0.5f); 

		//rgb[i][0] = backupMap[ii][0];
		//rgb[i][1] = backupMap[ii][1];
		//rgb[i][2] = backupMap[ii][2];
	}

}

void cmapTex::setColor(float v)
{
	float r,g,b;

	r = rgb[0][int(v*1000.0f + 0.5)][0];
	g = rgb[0][int(v*1000.0f + 0.5)][1];
	b = rgb[0][int(v*1000.0f + 0.5)][2];
	glColor4f(r,g,b,1.0);

}

void cmapTex::getColor(float v,float &r, float &g, float &b)
{
	//if (v < 0) cerr << "Less than zero " << v << "\n";
	r = rgb[0][int(v*1000.0f + 0.5)][0];
	g = rgb[0][int(v*1000.0f + 0.5)][1];
	b = rgb[0][int(v*1000.0f + 0.5)][2];
}



	

