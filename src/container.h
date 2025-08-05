#ifndef CONTAINER_H
#define CONTAINER_H

#pragma once

#include "CGL/CGL.h"
using namespace std;
using namespace CGL;

class Container {
public:
	float length;
	float height;
	float width;

	std::vector<float> vertices;

	Vector3D p1;
	Vector3D p2;
	Vector3D p3;
	Vector3D p4;
	Vector3D p5;
	Vector3D p6;
	Vector3D p7;
	Vector3D p8;


	Container() {
		this->length = 1;
		this->height = 1;
		this->width = 1;
		this->p1 = Vector3D(-0.5f, 0.5f, 0.5f);
		this->p2 = Vector3D(0.5f, 0.5f, 0.5f);
		this->p3 = Vector3D(0.5f, -0.5f, 0.5f);
		this->p4 = Vector3D(-0.5f, -0.5f, 0.5f);

		this->p5 = Vector3D(-0.5f, 0.5f, -0.5f);
		this->p6 = Vector3D(0.5f, 0.5f, -0.5f);
		this->p7 = Vector3D(0.5f, -0.5f, -0.5f);
		this->p8 = Vector3D(-0.5f, -0.5f, -0.5f);

		vertices_setup();
	}

	void setLength(float l);
	void setWidth(float w);
	void setHeight(float h);
protected:
	void vertices_setup() {
		//front face edges
		this->vertices.push_back(this->p1.x);
		this->vertices.push_back(this->p1.y);
		this->vertices.push_back(this->p1.z);
		this->vertices.push_back(this->p2.x);
		this->vertices.push_back(this->p2.y);
		this->vertices.push_back(this->p2.z);
		this->vertices.push_back(this->p2.x);
		this->vertices.push_back(this->p2.y);
		this->vertices.push_back(this->p2.z);
		this->vertices.push_back(this->p3.x);
		this->vertices.push_back(this->p3.y);
		this->vertices.push_back(this->p3.z);
		this->vertices.push_back(this->p3.x);
		this->vertices.push_back(this->p3.y);
		this->vertices.push_back(this->p3.z);
		this->vertices.push_back(this->p4.x);
		this->vertices.push_back(this->p4.y);
		this->vertices.push_back(this->p4.z);
		this->vertices.push_back(this->p4.x);
		this->vertices.push_back(this->p4.y);
		this->vertices.push_back(this->p4.z);
		this->vertices.push_back(this->p1.x);
		this->vertices.push_back(this->p1.y);
		this->vertices.push_back(this->p1.z);

		//back face edges
		this->vertices.push_back(this->p5.x);
		this->vertices.push_back(this->p5.y);
		this->vertices.push_back(this->p5.z);
		this->vertices.push_back(this->p6.x);
		this->vertices.push_back(this->p6.y);
		this->vertices.push_back(this->p6.z);
		this->vertices.push_back(this->p6.x);
		this->vertices.push_back(this->p6.y);
		this->vertices.push_back(this->p6.z);
		this->vertices.push_back(this->p7.x);
		this->vertices.push_back(this->p7.y);
		this->vertices.push_back(this->p7.z);
		this->vertices.push_back(this->p7.x);
		this->vertices.push_back(this->p7.y);
		this->vertices.push_back(this->p7.z);
		this->vertices.push_back(this->p8.x);
		this->vertices.push_back(this->p8.y);
		this->vertices.push_back(this->p8.z);
		this->vertices.push_back(this->p8.x);
		this->vertices.push_back(this->p8.y);
		this->vertices.push_back(this->p8.z);
		this->vertices.push_back(this->p5.x);
		this->vertices.push_back(this->p5.y);
		this->vertices.push_back(this->p5.z);

		//other edges
		this->vertices.push_back(this->p1.x);
		this->vertices.push_back(this->p1.y);
		this->vertices.push_back(this->p1.z);
		this->vertices.push_back(this->p5.x);
		this->vertices.push_back(this->p5.y);
		this->vertices.push_back(this->p5.z);

		this->vertices.push_back(this->p2.x);
		this->vertices.push_back(this->p2.y);
		this->vertices.push_back(this->p2.z);
		this->vertices.push_back(this->p6.x);
		this->vertices.push_back(this->p6.y);
		this->vertices.push_back(this->p6.z);

		this->vertices.push_back(this->p3.x);
		this->vertices.push_back(this->p3.y);
		this->vertices.push_back(this->p3.z);
		this->vertices.push_back(this->p7.x);
		this->vertices.push_back(this->p7.y);
		this->vertices.push_back(this->p7.z);

		this->vertices.push_back(this->p4.x);
		this->vertices.push_back(this->p4.y);
		this->vertices.push_back(this->p4.z);
		this->vertices.push_back(this->p8.x);
		this->vertices.push_back(this->p8.y);
		this->vertices.push_back(this->p8.z);


	}
}; // class Container

#endif
