#pragma once

#include "ofMain.h"

//each cell corresponds to a vertex in the mesh. this cell processes the physics of each vertex

class cell {
public:

	ofMesh * base;

	vector<int> neighbors;

	int index;

	bool bud;

	float radius;

	ofVec3f velocity;
	ofVec3f acceleration;
	
	cell(ofMesh* base_ ,int index_,float radius, bool bud_);

	void run();

	void applyForce(ofVec3f force);

};