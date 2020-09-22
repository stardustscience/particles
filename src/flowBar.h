#define MAX_STEPS 400

class flowBar
{
public:
	flowBar(float xp, float yp, int ci, int bet, float mean);
	void setAnchor(float ax, float ay) { anchorX = ax; anchorY = ay; };
	void draw(float flow, int RL,int ts);
	float getMean();
private:
	float xPos, yPos, left, right;
	float red, green, blue;
	float anchorX, anchorY;  // these denote the associated transect

	int betweens;

	float mean;

	float history[MAX_STEPS];
	int nsteps;

};


// every 2 frames = 11*30/3  == 110

