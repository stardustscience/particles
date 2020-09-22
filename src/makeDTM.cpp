
#include <math.h>
#include <iostream>
#include <GL/glut.h>
#include "makeDTM.h"

using namespace std;
/*

#define MKER 15

class makeDTM
{
public:
	makeDTM(int r, int c, float LLx, float LLy, float URx, float URy);
	void addPoint(float x, float y, float val);
	void finalStep();
private:
	float kdiam; // the kernel diameter
	float scale, trans; // the transformation to the grid.
	int nker;
	float **dem;
	float **weights;
	float kernel[MKER][MKER];// max kernel
	int rows, cols; //the rows and cols with padding
	int arows, acols;
};
*/
makeDTM::makeDTM(int c, int r, float LLx, float LLy, float URx, float URy)
{
	int i,j;
	arows = r;
	acols = c;

	//compress in x by cos y in this case 40 deg.0.76

	Xstep = (URx - LLx)/float(c);
	Ystep = (URy - LLy)/float(r);

	cerr << "X step deg " << Xstep << "\n";

	rows = r+MKER+2;
	cols = c+MKER+2;

	LLat = LLy; LLon = LLx;


	dem = new float*[rows];
	weights = new float*[rows];
	for(i=0;i<rows;++i)
	{
		dem[i] = new float[cols];
		weights[i] = new float[cols];
	}

	for(i=0;i<rows;++i)
		for(j=0;j<cols;++j)
		{
			dem[i][j] = 0.0;
			weights[i][j] = 0.000001;
		}

	transX = -LLx;
	transY = -LLy;
	scaleX = float(acols)/(URx - LLx);
	scaleY = float(arows)/(URy - LLy);

	count = 0;

}

void makeDTM::setupLighting()
{
	/*
	glEnable(GL_LIGHTING);
	float light_position[] = {-10.0,50.0,-50.0,0.0};
	glLightfv(GL_LIGHT0,GL_POSITION, light_position);
	float specReflection[] = { 1.0f, 1.0f, 1.0f, 1.0f };
    glMaterialfv(GL_FRONT, GL_SPECULAR, specReflection);
    glMateriali(GL_FRONT, GL_SHININESS, 30);

	glDisable(GL_LIGHTING);
*/
}

void makeDTM::getDTM(float **(&d),int &R, int &C)
{

	//cerr << "TEST GET 250 250 " << dem[250][250] << "\n";
	d = dem;
	R = arows; C = acols;
}

void makeDTM::addPoint(float x, float y, float val)
{
	int ix, iy;
	ix = int((x + transX)*scaleX);
	iy = int((y + transY)*scaleY);

	if((ix < 0)|| (iy < 0)) return;
	if((ix >= acols) || (iy >= arows)) return;

	addKernel(iy,ix,val);
	++count;
}

void makeDTM::divide()
{
	int i,j;

	cerr << "N POINTS " << count << "\n";
	for(i=0;i<rows;++i)
		for(j=0;j<cols;++j)
		{

			dem[i][j] = dem[i][j]/weights[i][j];
		}

}

void makeDTM::saveDEM()
{
	int i,j;
	float max, min, v, range;

	FILE *fout;

	max = -10000.0f;
	min = 10000.0f;
/*
	fscanf(infile, "%s %lf", DummyChars,&cellSize);// cell size
	fscanf(infile, "%s %d", DummyChars,&noData);
*/

	cerr << "SAVE " << count << "\n";
	for(i=0;i<rows;++i)
		for(j=0;j<cols;++j)
		{
			v = dem[i][j] ;
			if(v > max) max = v;
			if(v>0.01 && v < min) min = v;	
		}
	range = max-min;

	cerr << "min max " << min << " " << max << "\n";
	fout = fopen("SSH_Pacific.txt","w");	
	fprintf(fout,"cols %d\n",acols);
	fprintf(fout,"rows %d\n",arows);
	fprintf(fout,"LonCorner %f \n",LLon);
	fprintf(fout,"LatCorner %f \n",LLat);
	fprintf(fout,"xStep %f \n",Xstep);
	fprintf(fout,"yStep %f \n",Ystep);
	fprintf(fout,"NODAT -99999.0 \n" );

	for(i=0;i<arows;++i)
		for(j=0;j<acols;++j)
		{
			v = (dem[i+13][j+13]-min)/range ;
			if(v<0.0f)v = 0.0f;
			fprintf(fout,"%5.3f ",v);
			if(j%30 == 29) fprintf(fout,"\n");
		}
	fclose(fout);

}

void makeDTM::makeNormals()
{
	int r,c;
	float h1,h2,h3,dhx,dhy,dhz,len;

	float vscale = 80.0;

	normals = new Point3f *[rows];
	for(r=0;r<rows;++r)
		normals[r] = new Point3f[cols];

	for(r=0;r<rows;++r)
		for(c=0;c<cols;++c)
		{
			normals[r][c][0] = 0.0; 
			normals[r][c][1] = 1.0;
			normals[r][c][2] = 0.0;
		}

	for(r=0;r<rows-1;++r)
	{
		for(c=0;c<cols-1;++c)
		{
			h1 = dem[r][c]; h2 = dem[r][c+1];
			h3 = dem[r+1][c];
			dhx = vscale*(h1-h2); dhy = vscale*(h3-h1); // arbitrary v scale
			len = sqrt(dhx*dhx + dhy*dhy + 1.0);
			dhx /=len; dhy /= len; dhz = 1.0/len;

			normals[r][c][0] = dhx;
			normals[r][c][1] = dhy;
			normals[r][c][2] = dhz;
		}
		normals[r][cols-1][0] = dhx; //use the last value
		normals[r][cols-1][1] = dhy;
		normals[r][cols-1][2] = dhz;
	}
	/*
	for(c=0;c<cols;++c)
	{
		normals[rows-1][c][0] = 0.0;//normals[nRows-2][c][0]; //copy the last row
		normals[rows-1][c][1] = 1.0;//normals[nRows-2][c][1];
		normals[rows-1][c][2] = 0.0;//normals[rows-2][c][2];
	}
	*/
}

void makeDTM::draw()
{
	int i,j;
	float c;

	glEnable(GL_COLOR_MATERIAL);
	glEnable(GL_NORMALIZE);

	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);

	//glEnable(GL_DEPTH_TEST);
// the bal
	
	//glDisable(GL_LIGHTING);

	glPushMatrix();
	glTranslatef(-(float)cols/2.0,-(float)rows/2.0,-251.0);
	
	for(i=0;i<rows-1;++i)
	{
		glBegin(GL_TRIANGLE_STRIP);
		for(j=0;j<cols;++j)
		{
			c = dem[i][j]*0.9f;
			glColor3f(c,0.6*c,0.6f-c);

			//glColor3f(1.0,1.0,1.0);
			//glNormal3f(0.0f,0.0f,1.0f);
			glNormal3fv(normals[i][j]);
			glVertex3f(float(j),float(i),dem[i][j]*80.0);
			c = dem[i+1][j]*0.9f;
			glColor3f(c,0.6*c,0.6f-c);
			glNormal3fv(normals[i+1][j]);
			//glVertex2i(j,i+1);
			glVertex3f(float(j),float(i+1),dem[i+1][j]*80.0);
			
		}
		glEnd();
	}

	glPopMatrix();

	glPushMatrix();
	glTranslatef(0.0,20.0,25.0);
	//glutSolidSphere(50.0,15,15);
	
	glPopMatrix();



}



void makeDTM::mkKernel(float rad)
{
	float center,sq;
	float dx,dy,dist,sum,d2,h;
	float posVec[5]; // used in averaging
	int i,j,k,l,mid;
	kdiam = rad*2.0f;
	nker = int(kdiam+0.99f);
	if(nker%2 == 0)  nker = nker+1; // make sure it is odd

	cerr << "Kernel\n";

	posVec[0] = 0.1f; posVec[4] = 0.9f;
	posVec[1] = 0.3f; posVec[3] = 0.7f;
	posVec[2] = 0.5f;

	cerr << "Kernel " << nker << "\n"; ;

	//center = float(nker/2);// - 0.5f;
	center = float(nker/2) + 0.5;
	mid = nker/2;

	for(i=0;i<MKER;++i)
		for(j=0;j<MKER;++j)
			kernel[i][j] = 0.0f;

	for(i=0;i<nker;++i)
	{
		for(j=0;j<nker;++j)
		{
			sum = 0.0f;
			for(k=0;k<5;++k)// super sample
				for(l=0;l<5;++l)
				{
					dx = float(i)+posVec[k]-center;
					dy = float(j)+posVec[l]-center;
					d2 = dx*dx + dy*dy;
					dist = 0.0f;
					if(d2 > 0.0001)dist = sqrt(d2);
					h = 1.0-dist/rad;
					if(h>0.0) sum = sum+h;

				}
			kernel[i][j] = sum/25.0f;
		}
	}
	for(j=0;j<nker;++j)
		cerr <<kernel[6][j] << " ";
	cerr << "\n";
}

void makeDTM::addKernel(int iy, int ix,float val)
{
	int i,j;

	for(i=0;i<nker;++i)
		for(j=0;j<nker;++j)
		{
			dem[iy+i][ix+j] += kernel[i][j]*val;
			weights[iy+i][ix+j] += kernel[i][j];

		}
}