#define NCOLS 1000
#define MAXCOLS 1200
#define SMOOTH 1
#define STEPPED 2


class cmapTex
{
public:
	cmapTex();
	void load(char *fname,bool reverse, int cmapno);
	void mkTex2Dcmap();
	void loadStepped(char *fname);
	void stepMap(int nSteps);
	void moveMap(int step);
	void setColor(float v);

	void setColorMap();
	void getColor(float v,float &r, float &g, float &b);

	void drawTexColorMap();
	unsigned int tex2dNames[6];
private:
	int cmapCount;
	int basePos;  // the colormap start position
	float rgb[4][MAXCOLS][3]; // can handle 4 colormaps
	//ColorUCS *colDiffs;

	float original_rgb[400][3];
	float vals[400];
};
