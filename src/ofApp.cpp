#include "ofApp.h"

//A program that simulates growth of a biological organism by extruding and deforming a base mesh

/*
The following functions (applyCoulombsLaw, applyHookesLaw, attractToCenter, TotalEnergy) are functions adapted
from a javascript library for creating force directed graphs (graphs of vertices and edges where edges are modeled as springs) called springy.js
*/

//a function that treats each input cell as a charged particles and attracts or repels them proportional to the squared distance between them
void ofApp::applyCoulombsLaw(cell *a, cell *b) {

	ofVec3f d = plane.getVertex(a->index) - plane.getVertex(b->index);
	float dist = d.length() + 0.1;
	ofVec3f direction = d.getNormalized();

	a->applyForce((direction*repulsion) / (dist*dist*0.5));
	b->applyForce((direction*repulsion) / (dist*dist*-0.5));

}

//a function that simulates hooke's law between two cells, simulating a spring between each node based on the radius of each cell (averaged to produce an ideal spring length) attracting or repelling based on that
void ofApp::applyHookesLaw(cell *a, cell *b) {
	ofVec3f d = plane.getVertex(b->index) - plane.getVertex(a->index);
	float disp = ((b->radius + a->radius) / 2) - d.length();
	ofVec3f dir = d.getNormalized();

	a->applyForce(dir*(springConstant*disp*-0.5));
	b->applyForce(dir*(springConstant*disp*0.5));

}

//attracts cell to the origin to counteract constant reduction in size
void ofApp::attractToCenter(cell *p) {

	ofVec3f dir = plane.getVertex(p->index)*-1;
	p->applyForce(dir*(repulsion / 50));


}

//calculates the total energy of the system. this is used to stop the jittering effect but i havent gotten it to work yet
float ofApp::totalEnergy(vector<cell>*cells_) {
	float energy = 0;


	for (int i = 0; i < cells_->size(); i++) {
		float speed = (*cells_)[i].velocity.length();
		energy += 0.5*(*cells_)[i].radius*speed*speed;
	}


	return energy;
}

//--------------------------------------------------------------
void ofApp::setup(){


	
	plane.enableIndices();
	
	//this following series of for loops generates a triangulated mesh

	for (int y = 0; y < width; y++) {
	
		for (int x = 0; x < height; x++) {
			
			//creating vertices and indices

			plane.addVertex(ofVec3f((x*dim+(dim*.5*(y%2)))-(dim*width/2), y*dim - (dim*height/2),ofRandom(-1,1)));

			if ((y % 2 == 0)&&(x<width-1)) {
				plane.addIndex(x + y*width);
				plane.addIndex(x + y*width+1);
				plane.addIndex(x + (y+1)*width);
				
				plane.addIndex(x + y*width + 1);
				plane.addIndex(x + (y + 1)*width);
				plane.addIndex(x + (y + 1)*width + 1);
			}
			else if((y%2==1)&&(x<width-1)&&(y<height-1)){

				plane.addIndex(x + y*width);
				plane.addIndex(x + (y + 1)*width);
				plane.addIndex(x + (y + 1)*width + 1);

				plane.addIndex(x + y*width);
				plane.addIndex(x + y*width + 1);
				plane.addIndex(x + (y + 1)*width + 1);

			}

			//generating new cells for each vertex in the mesh. some of these meshes are buds with a lower mass than the average
			
			bool b = false;
			bool r = dim;

			if (ofRandom(100) > 80) {
				b = true;
				r = dim /2;
			}

			cells.push_back(cell(&plane, x + y*width,r, b));

		}
	}
	
	cout << plane.getNumVertices()<<"\n";

	//all this commented code has been implemented in the grow function of ofApp. It has some differences that might be important so I'm keeping it here for now 
	
	/*
	for (int i = 0; i < cells.size(); i++) {

		if (cells[i].bud) {

			adj.clear();

			//removing unneeded indices

			for (int j = plane.getNumIndices() - 1; j >= 0; j -= 3) {


				if (plane.getIndex(j) == cells[i].index) {

					adj.push_back(plane.getIndex(j - 1));
					adj.push_back(plane.getIndex(j - 2));

					plane.removeIndex(j);
					plane.removeIndex(j - 1);
					plane.removeIndex(j - 2);

				}
				else if (plane.getIndex(j - 1) == cells[i].index) {

					adj.push_back(plane.getIndex(j));
					adj.push_back(plane.getIndex(j - 2));

					plane.removeIndex(j);
					plane.removeIndex(j - 1);
					plane.removeIndex(j - 2);

				}
				else if (plane.getIndex(j - 2) == cells[i].index) {

					adj.push_back(plane.getIndex(j));
					adj.push_back(plane.getIndex(j - 1));

					plane.removeIndex(j);
					plane.removeIndex(j - 1);
					plane.removeIndex(j - 2);

				}


			}
			//removing duplicates
			std::sort(adj.begin(), adj.end());
			adj.erase(unique(adj.begin(), adj.end()), adj.end());

			//generating center point of adjacent vertices
			ofVec3f avg = plane.getVertex(cells[i].index);

			//ordering adjacent vertices

			ofVec3f norm;

			norm.set((plane.getVertex(adj[0]) - avg).getCrossed(plane.getVertex(adj[1]) - avg).getNormalized());

			for (int j = 0; j < adj.size() - 2; j++) {
				ofVec3f edge0;
				edge0.set((plane.getVertex(adj[j]) - avg).getNormalized());

				int bestIndex = j + 1;
				float bestValue = -1.0;

				for (int k = j + 1; k < adj.size(); k++) {
					ofVec3f edge1;
					edge1.set((plane.getVertex(adj[k]) - avg).getNormalized());

					if (edge0.getCrossed(edge1).dot(norm) > 0) {
						float currentValue = edge0.dot(edge1);

						if (currentValue > bestValue) {
							bestIndex = k;
							bestValue = currentValue;
						}
					}
				}

				iter_swap(adj.begin() + j + 1, adj.begin() + bestIndex);
			}

			ofVec3f sum;
			sum.set(0, 0, 0);

			for (int j = 0; j < floor(adj.size() / 2) + 1; j++) {
				sum += plane.getVertex(adj[j]);
			}
			sum /= floor(adj.size() / 2) + 1;

			ofVec3f sum2;
			sum2.set(0, 0, 0);

			for (int j = floor(adj.size() / 2); j <= adj.size(); j++) {
				sum2 += plane.getVertex(adj[j%adj.size()]);
			}
			sum2 /= floor(adj.size() / 2) + 1;


			//plane.setVertex(cells[i], plane.getVertex(cells[i]) - dim*.5);
			plane.setVertex(cells[i].index, sum);

			plane.addVertex(sum2);
			newCells.push_back(cell(&plane, plane.getNumVertices() - 1, cells[i].radius*.5, true));

			//linking new vertices to ordered adjacent vertices
			for (int j = 0; j < floor(adj.size() / 2); j++) {
				//plane.addIndex(plane.getNumVertices() - 2);
				plane.addIndex(cells[i].index);
				plane.addIndex(adj[j]);
				plane.addIndex(adj[(j + 1)]);
			}

			plane.addIndex(cells[i].index);
			plane.addIndex(plane.getNumVertices() - 1);
			plane.addIndex(adj[0]);

			plane.addIndex(cells[i].index);
			plane.addIndex(plane.getNumVertices() - 1);
			plane.addIndex(adj[floor(adj.size() / 2)]);

			for (int j = floor(adj.size() / 2); j < adj.size(); j++) {
				plane.addIndex(plane.getNumVertices() - 1);
				plane.addIndex(adj[j]);
				plane.addIndex(adj[(j + 1) % adj.size()]);
			}
		}

	}
	cells.insert(cells.end(), newCells.begin(), newCells.end());
	

	cout << plane.getNumVertices();

	
	sort(cells.begin(), cells.end(), [](const cell &a, const cell&b) {
		return a.index < b.index;
	});
	*/
	

}


//--------------------------------------------------------------
void ofApp::update(){

	if (run) {

		//applies coulombs law between every vertex and every other vertex. this ensures the mesh generally tries to avoid self intersection
		for (int i = 0; i < plane.getNumVertices(); i++) {
			for (int j = 0; j < plane.getNumVertices(); j++) {
				if (i != j) {
					applyCoulombsLaw(&cells[i], &cells[j]);
				}
			}
		}


		//gets each triangle in the mesh and applies hooke's law between each edge
		for (int i = 0; i < plane.getNumIndices(); i += 3) {

			applyHookesLaw(&cells[plane.getIndex(i)], &cells[plane.getIndex(i + 1)]);
			applyHookesLaw(&cells[plane.getIndex(i + 1)], &cells[plane.getIndex(i + 2)]);
			applyHookesLaw(&cells[plane.getIndex(i + 2)], &cells[plane.getIndex(i)]);

		}


		//attracts every point to center
		for (int i = 0; i < cells.size(); i++) {
			attractToCenter(&cells[i]);
		}
	}

	//this would stop the simulation when the mesh is in its ideal shape to stop jittering but it doesnt work correctly yet
	if (totalEnergy(&cells) < 0.01) {
		//run = false;
	}

	//runnin each cell
	for (int i = 0; i < cells.size(); i++) {
		cells[i].run();
	}


}

//--------------------------------------------------------------
void ofApp::draw(){

	camera.begin();
	ofBackground(0);

	ofSetColor(100);
	plane.drawFaces();

	ofSetColor(255);
	plane.drawWireframe();
	
	//drawing an ellipse at every seed point
	for (int i = 0; i < cells.size(); i++) {
		if (cells[i].bud) {
			ofSetColor(0, 255, 0);
			ofDrawSphere(plane.getVertex(cells[i].index), .2);
		}
	}
	
	camera.end();
}

void ofApp::grow() {

	//clearing the temporary vector of cells

	newCells.clear();
	
	//goes throughe each cell and grows the mesh based on which cells are buds
	for (int i = 0; i < cells.size(); i++) {

		if (cells[i].bud) {

			adj.clear();

			//removing the links between the bud vertex and the mesh, and storing the adjacent vertices

			for (int j = plane.getNumIndices() - 1; j >= 0; j -= 3) {


				if (plane.getIndex(j) == cells[i].index) {

					adj.push_back(plane.getIndex(j - 1));
					adj.push_back(plane.getIndex(j - 2));

					plane.removeIndex(j);
					plane.removeIndex(j - 1);
					plane.removeIndex(j - 2);

				}
				else if (plane.getIndex(j - 1) == cells[i].index) {

					adj.push_back(plane.getIndex(j));
					adj.push_back(plane.getIndex(j - 2));

					plane.removeIndex(j);
					plane.removeIndex(j - 1);
					plane.removeIndex(j - 2);

				}
				else if (plane.getIndex(j - 2) == cells[i].index) {

					adj.push_back(plane.getIndex(j));
					adj.push_back(plane.getIndex(j - 1));

					plane.removeIndex(j);
					plane.removeIndex(j - 1);
					plane.removeIndex(j - 2);

				}


			}
			
			//removing duplicates in the adjacent vertex vector
			std::sort(adj.begin(), adj.end());
			adj.erase(unique(adj.begin(), adj.end()), adj.end());

			//generating center point of adjacent vertices
			ofVec3f avg = plane.getVertex(cells[i].index);

			//generating the normal to the plane of the vertices
			ofVec3f norm;
			norm.set((plane.getVertex(adj[0]) - avg).getCrossed(plane.getVertex(adj[1]) - avg).getNormalized());

			//ordering the adjacent vertices in the array based on their angles from a randomly chosen start point, and whether their dot product is in the same direction or opposite to the normal vector
			//this is based on an algorithm i found online to order vertices in a circle
			for (int j = 0; j < adj.size() - 2; j++) {
				ofVec3f edge0;
				edge0.set((plane.getVertex(adj[j]) - avg).getNormalized());

				int bestIndex = j + 1;
				float bestValue = -1.0;

				for (int k = j + 1; k < adj.size(); k++) {
					ofVec3f edge1;
					edge1.set((plane.getVertex(adj[k]) - avg).getNormalized());

					if (edge0.getCrossed(edge1).dot(norm) > 0) {
						float currentValue = edge0.dot(edge1);

						if (currentValue > bestValue) {
							bestIndex = k;
							bestValue = currentValue;
						}
					}
				}

				iter_swap(adj.begin() + j + 1, adj.begin() + bestIndex);
			}
			//choosing a random index to rotate the vector of indices by. this means the new vertex will be added in a random position around the original vertex
			int r = ofRandom(0, adj.size());

			std::rotate(adj.begin(), adj.begin()+r, adj.end());

			//creating two new vertices inside the circle of vertices
			ofVec3f sum;
			sum.set(0, 0, 0);

			for (int j = 0; j < floor(adj.size() / 2) + 1; j++) {
				sum += plane.getVertex(adj[j]);
			}
			sum /= floor(adj.size() / 2) + 1;

			ofVec3f sum2;
			sum2.set(0, 0, 0);

			for (int j = floor(adj.size() / 2); j <= adj.size(); j++) {
				sum2 += plane.getVertex(adj[j%adj.size()]);
			}
			sum2 /= floor(adj.size() / 2) + 1;


			//setting the position of the bud vertex based on one of the previously generated vertices
			plane.setVertex(cells[i].index, sum + norm);

			//adding a new vertex
			plane.addVertex(sum2 + norm);
			newCells.push_back(cell(&plane, plane.getNumVertices() - 1, cells[i].radius*.5, false));

			//linking new vertices to ordered adjacent vertices
			for (int j = 0; j < floor(adj.size() / 2); j++) {
				plane.addIndex(cells[i].index);
				plane.addIndex(adj[j]);
				plane.addIndex(adj[(j + 1)]);
			}

			plane.addIndex(cells[i].index);
			plane.addIndex(plane.getNumVertices() - 1);
			plane.addIndex(adj[0]);

			plane.addIndex(cells[i].index);
			plane.addIndex(plane.getNumVertices() - 1);
			plane.addIndex(adj[floor(adj.size() / 2)]);

			for (int j = floor(adj.size() / 2); j < adj.size(); j++) {
				plane.addIndex(plane.getNumVertices() - 1);
				plane.addIndex(adj[j]);
				plane.addIndex(adj[(j + 1) % adj.size()]);
			}
		}

	}

	//insering the newCells into the original cells vector
	cells.insert(cells.end(), newCells.begin(), newCells.end());
	//sorting the cells so the correspond to the order of the vertices in the mesh
	sort(cells.begin(), cells.end(), [](const cell &a, const cell&b) {
		return a.index < b.index;
	});
}

//--------------------------------------------------------------
void ofApp::keyPressed(int key){
	//grows the mesh every time space is pressed
	if (key == ' ') {
		grow();
	}
	

}

//--------------------------------------------------------------
void ofApp::keyReleased(int key){

}

//--------------------------------------------------------------
void ofApp::mouseMoved(int x, int y ){

}

//--------------------------------------------------------------
void ofApp::mouseDragged(int x, int y, int button){

}

//--------------------------------------------------------------
void ofApp::mousePressed(int x, int y, int button){

}

//--------------------------------------------------------------
void ofApp::mouseReleased(int x, int y, int button){

}

//--------------------------------------------------------------
void ofApp::mouseEntered(int x, int y){

}

//--------------------------------------------------------------
void ofApp::mouseExited(int x, int y){

}

//--------------------------------------------------------------
void ofApp::windowResized(int w, int h){

}

//--------------------------------------------------------------
void ofApp::gotMessage(ofMessage msg){

}

//--------------------------------------------------------------
void ofApp::dragEvent(ofDragInfo dragInfo){ 

}
