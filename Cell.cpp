#include "cell.h"

cell::cell(ofMesh*base_,int index_,float radius_, bool bud_) :base(base_),index(index_),radius(radius_), bud(bud_) {

}

void cell::run() {

	velocity += acceleration;
	
	velocity *= .5;
	
	if (velocity.length() > 1) {
		velocity.normalize();
		//velocity *= 5;
	}

	(*base).setVertex(index, (*base).getVertex(index) + velocity);
	acceleration.set(0, 0, 0);
}

void cell::applyForce(ofVec3f force) {

	acceleration = acceleration + (force / radius);

}