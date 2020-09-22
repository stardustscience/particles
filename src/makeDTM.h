#define MKER 25

class makeDTM
{
public:
	makeDTM(int c, int r, float LLx, float LLy, float URx, float URy);

	void setupLighting();
	void addPoint(float x, float y, float val);
	void mkKernel(float rad);
	void addKernel(int ix, int iy,float val);
	void divide();
	void saveDEM();
	void draw();
	float **dem;
	void getDTM(float **(&d),int &R, int &C);
	void makeNormals();
	void finalStep();	
	int rows, cols; //the nominal rows and cols
private:
	
	float kdiam; // the kernal diameter
	float LLat, LLon;  // the lat and lon corners
	float Xstep,Ystep;
	float scaleX, scaleY, transX,transY; // the transformation to the grid.
	int nker;

	typedef float Point3f[3];
	Point3f **normals;  // for shaded bathy (CW)
	
	float **weights;
	float kernel[MKER][MKER];// max kernel

	int arows, acols;
	int count;
};