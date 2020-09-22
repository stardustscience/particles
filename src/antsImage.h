

class antsImage
{

public:
	antsImage();
	void capture();
	void captureAnts();
	void drawPrev();
	void drawTransparentAnts();

private:
	unsigned char *current;
	unsigned char *ants;

};
