

#include <GL/glut.h>
#include <iostream>

#include "MPAS_vis.h" 
#include "savePNG.h"
#include "antsImage.h"
#include "arrow.h"
//#include "readpng.h"
//#include "Streak.h"

#define SCR_WID 1280
#define MOV_HT 720  // for 1280x720 format


using namespace std;

MPAS_vis *MP;

float winWid,winHeight;
float rx,ry;
int timer;
int startYr, currentYr;


unsigned char *readbuf;

int Layer;

int frameCount;
int maxFrames;

float Nslice;

savePNG *imageWriter;
antsImage *antsFrame;

//Streak *streams;

bool drawAnts;
bool resam;
bool drawTransects;

bool shallow;
bool deep;
Arrow *a1, *a2, *a3, *a4,*a5, *a6, *a7, *a8;



void drawTimeBar()
{
	float left, right,x, yrStep;
	int i,yrs;

	char date[200];

	currentYr = startYr + frameCount / (BETWEENS * 12);

	sprintf(date, "%d",currentYr);
	left = 120.0f; 
	right = left + 1020.0f*frameCount / maxFrames;
	glColor3f(0.2f, 0.2f, 0.2f);
	glRectf(0.0f, 300.0f, 1280.0f, 320.0f);
	glColor3f(0.2f, 0.7f, 0.7f);
	glRectf(left, 305.0f, right, 315.0f);

	yrs = 1+maxFrames / (BETWEENS*12);

	yrStep = 1020.0f / float(yrs - 1);

	glColor3f(1.0f, 1.0f, 1.0f);
	glBegin(GL_LINES);
	for (i = 0; i < yrs; ++i)
	{
		x = yrStep*i + left;
		glLineWidth(1.0f);

		if(i%10==0){
			glLineWidth(2.0f);
			glVertex2f(x, 295);
			glVertex2f(x, 325);
		}
		else {

			glVertex2f(x, 300);
			glVertex2f(x, 320);
		}
	}
	glEnd();

	glRasterPos2f(60.0f, 310.0f);
	for (i = 0; i < 4; ++i) glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18, date[i]);

}

void draw( void )
{
	glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
	glColor3f(1.0, 1.0, 0.0);

	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glPushMatrix();
	glTranslatef(140.0f, 320.0f, 0.0f);
	//glRotatef(-30.0f, 1.0f, 0.0f, 0.0f);
		//MP->drawSlice(Nslice);
		//MP->drawVecField(DEEP);  // the streamlets	

		if (drawAnts) {
			MP->drawPathTraces_D(ANTS, 0.15f);
			a1->draw();
			a2->draw();
			MP->drawPathTraces_S(ANTS, 0.15f);  // the smart traces illustration
			a3->draw();
			a4->draw();
			a5->draw();
			a6->draw();
			a7->draw();
			a8->draw();
			//MP->drawMask(readbuf);
		}

		// the fireflies
		MP->drawVecField(DEEP);  // the pathlets
		MP->drawVecField(SHALLOW);  // the pathlets

	

		glColor3f(1.0f, 1.0f, 0.0f);
	//	glRectf(0.0, 0.0, 200.0, 200.0);

		MP->drawBG();	
		if (drawTransects) {
			MP->draw_Flow_Indicators_D(frameCount);
			MP->draw_Flow_Indicators_S(frameCount);
		}

	//	MP->drawCurrentTransect();  // draws the ocean transect below need hi res data	
		
	glPopMatrix();

	drawTimeBar();


/*
	glColor3f(0.5,0.99,0.99);
	glBegin(GL_LINES);
	glVertex2f(0.0f, Nslice*10.0f + 300.0f);
	glVertex2f(1000.0f, Nslice*10.0f + 300.0f);
	glEnd();
*/
	glColor3f(0.2f, 0.8f, 0.2f);
	glBegin(GL_LINES);
	//glVertex2f(0.0f, 290.0f-Layer*3.5f) ;
	//glVertex2f(1000.0f, 290.0f-Layer*3.5f);
	glVertex2f(0.0f, 290.0f - 45.0f*3.5f);
	glVertex2f(1000.0f, 290.0f - 45.0f*3.5f);
	glEnd();

	glutSwapBuffers();
}

void ddd()
{
	int i, j, iv, indx;

	glColor3f(1.0, 1.0, 0.0);

	glBegin(GL_POINTS);
	for (j = 0; j < 700; ++j)
		for (i = 0; i < 1000; ++i)
		{	
			indx = (j*1000 + i) * 3;
			iv = readbuf[indx];// +mask[indx + 1] + mask[indx + 2];
			if (iv > 0)glVertex2i(i,j);
		}

	glEnd();

}

void makeAntMask(bool first)
{
	int ct=0;

	glClearColor(0.0f , 0.0f, 0.0f,1.0f);

	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	MP->drawPathTraces_S(MASK,1.0f);  // note: this is a lat, lon draw
	glReadBuffer(GL_BACK);
	glReadPixels(0, 0, 1000, 700, GL_RGB, GL_UNSIGNED_BYTE, readbuf);
	glutSwapBuffers();
	MP->setSamplePts(readbuf,0,first);	

	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	MP->drawPathTraces_D(MASK, 1.0f);
	glReadBuffer(GL_BACK);
	glReadPixels(0, 0, 1000, 700, GL_RGB, GL_UNSIGNED_BYTE, readbuf);
	glutSwapBuffers();
	MP->setSamplePts(readbuf, 1,first);
	if (MP->colorScheme == SCHEME_2)glClearColor(0.35f, 0.35f, 0.35f, 1.0);
	
	resam = true;
}

void redraw()
{
	//drawMask();
	//
	//if (resam)ddd();
	//else
	draw();

	//glutSwapBuffers();
}



void setLighting()
{
	glEnable(GL_LIGHTING);
	float light_position[] = {-20.0,10.0,40.0,0.0};
	glLightfv(GL_LIGHT0,GL_POSITION, light_position);
	float specReflection[] = { 0.8f, 0.8f, 0.8f, 1.0f };
    glMaterialfv(GL_FRONT, GL_SPECULAR, specReflection);
    glMateriali(GL_FRONT, GL_SHININESS, 150);

	glDisable(GL_LIGHTING);


}

void writePNG(char *fn)
{
	unsigned char *pixels;


	//pixels = new unsigned char[1000 * 680 * 3];
	pixels = new unsigned char[1280 * 720 * 3];

	glReadBuffer(GL_FRONT);
	glReadPixels(0, 300, 1280, 720, GL_RGB, GL_UNSIGNED_BYTE,pixels);

	cerr << "OUTPUT frame\n";
	imageWriter->writeImage(fn, 1280, 720, pixels, "nm");

	delete[] pixels;

}

void makeAntsMovie()
{
	int b, m;
	float p;
	int Months, Start;
	bool first = true;

	antsFrame = new antsImage();
	char frame[40];
	int mn; // model frame number
	startYr = 1930;
	Months =  599;  // 239 maxshould be n-1
	Start = (startYr - 1900) * 12 + 1; // 1945
	maxFrames = Months*BETWEENS;



	glClear(GL_COLOR_BUFFER_BIT);
	glutSwapBuffers();
	glClear(GL_COLOR_BUFFER_BIT );
	antsFrame->capture();


	for (m = 0; m < Months; ++m)
	{
		mn = m + Start;
		MP->loadLayers(mn);
		for (b = 0; b < BETWEENS; ++b) {
			p = float(b) / float(BETWEENS);
			MP->interplolateLayers(p, 0);
			MP->interplolateLayers(p, 1);

			glClear(GL_COLOR_BUFFER_BIT);	
			glReadBuffer(GL_BACK);

			glPushMatrix();
			glTranslatef(140.0f, 320.0f, 0.0f);
			if (deep) MP->drawPathTraces_D(ANTS, 1.0f);
				
			if (shallow)MP->drawPathTraces_S(ANTS, 1.0f);  // the smart traces illustration
			glPopMatrix();

			antsFrame->captureAnts();
			glClear(GL_COLOR_BUFFER_BIT);
			antsFrame->drawPrev();
			glPushMatrix();
				glTranslatef(140.0f, 320.0f, 0.0f);
				MP->drawBG();			
			glPopMatrix();
			antsFrame->drawTransparentAnts();
			antsFrame->capture();


			if (drawTransects) {	
				glPushMatrix();
				glTranslatef(140.0f, 320.0f, 0.0f);
					MP->measureTransectFlows();
					if (deep) {
						MP->draw_Flow_Indicators_D(frameCount);
						a1->draw();
						a2->draw();
					}
					if (shallow) {
						MP->draw_Flow_Indicators_S(frameCount);
						a3->draw();
						a4->draw();
						a5->draw();
						a6->draw();
						a7->draw();
						a8->draw();
					}

					
				glPopMatrix();
			}
			drawTimeBar();

			glutSwapBuffers();
			sprintf(frame, "frames_1/F_%04d.png", frameCount);
			cerr << frame << "\n";
			writePNG(frame);


			++frameCount;

		}
	}
}

void makeMovie() 
{
	int b,m;
	float p;
	int Months, Start;
	bool first=true;

	char frame[40];
	int mn; // model frame number
	startYr = 1930;
	//Months = 359;  // 239 maxshould be n-1
	Months = 599;  // 239 maxshould be n-1
	Start = (startYr-1900)*12+1; // 1945
	maxFrames = Months*BETWEENS;

	for (m = 0; m < Months; ++m)
	{
		mn = m + Start;
		MP->loadLayers(mn);

		for (b = 0; b < BETWEENS; ++b) {
			p = float(b) / float(BETWEENS);
			MP->interplolateLayers(p, 0);
			MP->interplolateLayers(p, 1);
		//	computePeakTransectHR(int n, float &maxp, float &pxMax, float &pyMax);

			MP->measureTransectFlows();
			makeAntMask(first);
			first = false;

			draw();// redraw();

			//if (frameCount < 10) sprintf(frame, "frames/F_000%d.png", frameCount);
			//if ((frameCount < 10) && (frameCount < 100)) sprintf(frame, "frames/F_00%d.png", frameCount);
			//if (frameCount > 99) sprintf(frame, "frames/F_0%d.png", frameCount);
			sprintf(frame, "frames_1/F_%04d.png", frameCount);
			cerr << frame << "\n";
			++frameCount;
			writePNG(frame);
		}
	}
	MP->getTransectMeans();
}

void motion(int x, int y) 
// called when a mouse is in motion with a button down
{
 	rx = x; ry = winHeight - y;
}


void mousebutton(int button, int state, int x, int y)
{
	float Lat, Lon;
	if (button == GLUT_LEFT_BUTTON && state == GLUT_DOWN)

	{
		rx = x; ry = winHeight - y;

		MP->setPoint(float(x-140), y);

		redraw();
	}
}

void keyboard(unsigned char key, int x, int y)
// x and y givethe mouse pos
{
	//cerr << "Key " << key << " " << int(key) << "\n";

	//cerr << "\n";

	if (key == '.') {
		Nslice = Nslice + 1.0f; MP->computeAMOCflows(Nslice);
		cerr << "Lat " << Nslice << "\n";
	}
	if (key == ',') {
		Nslice = Nslice - 1.0f; MP->computeAMOCflows(Nslice);
		cerr << "Lat " << Nslice << "\n";
		if (Nslice < 0.1f) Nslice = 0.1;
	}

	if (key == '[') {
		Layer = Layer + 1;
		if (Layer >= 80) Layer = 79;
		MP->defineLand(Layer);
	}
	if (key == ']')
	{
		Layer = Layer - 1;
		if (Layer < 0) Layer = 0;
		MP->defineLand(Layer);
	}

	if (key == 'd')
	{
		redraw();
	}


	if (key == 'M') makeAntMask(true);

	if (key == 'a')  drawAnts = !drawAnts;
	if (key == 'f') MP->drawTraces = !MP->drawTraces;
	if (key == 't') drawTransects = !drawTransects;// MP->drawTransct = !MP->drawTransct;
	if (key == 'c') MP->drawCurrents = !MP->drawCurrents;
	if (key == 'p') MP->drawLongTrace = !MP->drawLongTrace;
	if (Layer >= IN_LAYERS) Layer = IN_LAYERS - 1;

	//cerr << "Layer " << Layer << " " << MP->layerDepths[Layer]<< "\n";

	//if (key == 's') showUV=!showUV;
	//if (key == 'p') MP->printFlowGraph(Nslice);

	if (key == 's') MP->changeTraceLen(-5);
	if (key == 'l') MP->changeTraceLen(5);
	if (key == '1')
	{
		MP->topPaths = !MP->topPaths;
		if(MP->topPaths)
			MP->makeTraces(1);
		MP -> viewTransect = SHALLOW;
	}
	if (key == '2')
	{
		MP->bottomPaths = !MP->bottomPaths;
		if (MP->bottomPaths)
			MP->makeTraces(2);
		MP->viewTransect = DEEP;
	}

	if (key == '3')
		MP->loadLayers(300);
	if (key == '4')
		MP->loadLayers(400);
	if (key == '5')
		MP->loadLayers(500);
	if (key == '6')
		MP->loadLayers(600);
	if (key == '7')
		MP->loadLayers(700);
	
	if (key == 'h') MP->genSaveCompact(75);
		//MP->saveBinaryLayers(51*12+3);

	if (key == 'x') MP->incrementTransect();

	if(key == 'w')writePNG("test.png");

	if (key == 'm')
		makeAntsMovie();// 
		//makeMovie();
	redraw();

}

int main(int argc, char *argv[])
{
	cerr << "hello world\n";
	frameCount = 1;
	currentYr = 0;

	shallow = true;
	deep = true;

	readbuf = new unsigned char[1000 * 700 * 3];
	for (int i = 0; i < 1000 * 700 * 3; ++i)readbuf[i] = 100;

	Layer = 0;
	frameCount = 0;
	Nslice = 0.1f;
	maxFrames = 30 * 11;

	imageWriter = new savePNG(600, 600);	
	a1 = new Arrow(250.0f, 200.0f, 43.0f, -15.0f, 0.1, 0.1, 0.4,0.25);
	a2 = new Arrow(620.0f, 500.0f, 8.0f, -60.0f, 0.1, 0.1, 0.4, 0.25);

	a3 = new Arrow(400.0f, 250.0f, -70.0f, -15.0f, 0.8, 0.5, 0.3, 0.15);  // big gyre
	a7 = new Arrow(340.0f, 384.0f, 70.0f, 0.0f, 0.8, 0.5, 0.3, 0.15);  // big gyre
	a4 = new Arrow(427.0f, 550.0f, 18.0f, -55.0f, 0.8, 0.1, 0.8, 0.15);
	a5 = new Arrow(527.0f, 550.0f, 6.0f, 55.0f, 0.8, 0.1, 0.8, 0.15);
	a6 = new Arrow(220.0f, 153.0f, -60.0f, -5.0f, 0.8, 0.3, 0.3, 0.15);
	a8 = new Arrow(645.0f, 591.0f, 30.0f,45.0f, 0.8, 0.3, 0.3, 0.15);

	int r,c;

	drawAnts = true;
	drawTransects = true;
	resam = false;

	//timer = 0;

	MP = new MPAS_vis();
	MP->readMPAS(60*12+1);
	//MP->loadLayers(100);
	MP->loadTransects_2();

	winWid = SCR_WID;//720.0;
	winHeight = 1020.0f;//360.0;

	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE );
	glutCreateWindow("MPAS N Atlantic");
	glutPositionWindow(200,50);
	glutReshapeWindow(winWid,winHeight);
	
	if(MP->colorScheme == SCHEME_1)	glClearColor(0.0,0.0,0.0,1.0);


	if (MP->colorScheme == SCHEME_2)glClearColor(0.35f,0.35f,0.35f, 1.0);

	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glOrtho(0.0,1280.0,0.0,1020.0, -300.0, 500.0);

//	glFrustum(-300.0,300.0,-200.0,200.0,460.0,1500.0);
//	glTranslatef(0.0,0.0,-560.0);

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

	//MP->loadColorMaps();

	glutDisplayFunc(redraw);
	//glutIdleFunc(redraw);
	glutMotionFunc( motion);	
	glutMouseFunc( mousebutton);
	glutKeyboardFunc( keyboard );
	glutMainLoop();

	return 0;
}