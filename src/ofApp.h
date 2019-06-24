#pragma once

#include "ofMain.h"
#include "../Cell.h"
class ofApp : public ofBaseApp{

	public:
		void setup();
		void update();
		void draw();

		void keyPressed(int key);
		void keyReleased(int key);
		void mouseMoved(int x, int y );
		void mouseDragged(int x, int y, int button);
		void mousePressed(int x, int y, int button);
		void mouseReleased(int x, int y, int button);
		void mouseEntered(int x, int y);
		void mouseExited(int x, int y);
		void windowResized(int w, int h);
		void dragEvent(ofDragInfo dragInfo);
		void gotMessage(ofMessage msg);
		void applyCoulombsLaw(cell *a, cell *b);
		void applyHookesLaw(cell *a, cell *b);
		void attractToCenter(cell *p);
		bool wayToSort(const cell &i,const cell &j);
		float totalEnergy(vector<cell>* cells_);
		void grow();
		ofVboMesh plane;
		
		int width = 10;
		int height = 10;
		float dim = 20;
		
		vector<cell> cells;
		vector<cell> newCells;
		vector<int> adj;

		ofEasyCam camera;

		float repulsion = 400;
		float springConstant = 200;

		bool doRender = false;
		float foceTimer = false;
		float rotSpeed = 0;
		
		bool run = true;
};
