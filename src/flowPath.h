
class flowPath
{
public:
	flowPath(float *xp, float *yp, int len, float flow, float lenscale);
	void makePath();
	void smooth(int sm);
	void drawPath(bool close);	
	float *xPath, *yPath;	
	int pLen;
private:
	float flowRate;
	float wid;
	float pathScale;

	float *leftX, *leftY, *rightX, *rightY;
	float *tmp;
	float red, green, blue;

};

	