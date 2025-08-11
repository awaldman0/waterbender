#ifndef PARTICLE_H
#define	PARTICLE_H

#pragma once
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include "CGL/CGL.h"
#include "container.h"
#include <iostream>
using namespace std;
using namespace CGL;

const float GRAV_CONST = 6.674e-11;
const float GRAV_MULTIPLIER = 400;
const float INFLUENCE_MULTIPLIER = 1.2f;
const float BOUNCE_MULTIPLIER = -0.6;

glm::vec3 sphericalToCartesian(float r, float theta, float phi) {
	float x = r * sin(theta) * cos(phi);
	float y = r * cos(theta);
	float z = r * sin(theta) * sin(phi);
	return glm::vec3(x, y, z);
};


class Particle {
public:
	Vector3D center;
	Vector3D velocity = Vector3D(0.0f, 0.0f, 0.0f);
	float radius;
	float drawing_radius;

	Particle() {
		this->center = Vector3D(0.0f, 0.0f, 0.0f);
		this->radius = 0.05f;
		this->drawing_radius = 0.03f;
	}

	Particle(float x, float y, float z, float radius, float drawing_radius) {
		this->center = Vector3D(x, y, z);
		this->radius = radius;
		this->drawing_radius = drawing_radius;
	}

	void draw(int num_segments, std::vector<float>* circle_pts, glm::mat4* inv) {
		//for (int i = 0; i < num_segments; i++) {
		//	float theta = 2.0f * M_PI * float(i) / float(num_segments); //get the current angle
		//	float x = this->radius * cosf(theta); //calculate the x component
		//	float y = this->radius * sinf(theta); //calculate the y component
		//	glm::vec4 pt = glm::vec4(x + this->center.x, y + this->center.y, this->center.z, 1.0);

		//	circle_pts->push_back(pt.x);
		//	circle_pts->push_back(pt.y);
		//	circle_pts->push_back(pt.z);
		//}
		int stacks = 10;
		int sectors = 10;

		// generate circumference points using integer steps
		for (float i = 0.0f; i <= stacks; i++) {
			float theta1 = (i / stacks) * glm::pi<float>();
			float theta2 = (i + 1) / stacks * glm::pi<float>();
			for (float j = 0.0f; j < sectors; j++) {
				float phi1 = j / sectors * 2 * glm::pi<float>();
				float phi2 = (j + 1) / sectors * 2 * glm::pi<float>();
				glm::vec3 c = glm::vec3(this->center.x, this->center.y, this->center.z);

				glm::vec3 v1 = sphericalToCartesian(this->drawing_radius, theta1, phi1) + c;
				glm::vec3 v2 = sphericalToCartesian(this->drawing_radius, theta1, phi2) + c;
				glm::vec3 v3 = sphericalToCartesian(this->drawing_radius, theta2, phi1) + c;
				glm::vec3 v4 = sphericalToCartesian(this->drawing_radius, theta2, phi2) + c;

				// Triangle 1: v1-v2-v3
				circle_pts->insert(circle_pts->end(), { v1.x, v1.y, v1.z });
				circle_pts->insert(circle_pts->end(), { v2.x, v2.y, v2.z });
				circle_pts->insert(circle_pts->end(), { v3.x, v3.y, v3.z });

				// Triangle 2: v2-v4-v3
				circle_pts->insert(circle_pts->end(), { v2.x, v2.y, v2.z });
				circle_pts->insert(circle_pts->end(), { v4.x, v4.y, v4.z });
				circle_pts->insert(circle_pts->end(), { v3.x, v3.y, v3.z });

			}
		}
	}

	void updatePosition(glm::vec4* gravity, Container* c, glm::mat4* invRotation, vector<Particle*>* particles, float deltaTime) {
		//bounds checking
		int x = 0;
		int y = 0;
		int z = 0;
		if (boundsCollision(c, &x, &y, &z)) {
			//reverse particle's direction in the axis of each boundary it violated
			if (x == 1) {
				this->velocity.x *= BOUNCE_MULTIPLIER;
			}
			if (y == 1) {
				this->velocity.y *= BOUNCE_MULTIPLIER;
			}
			if (z == 1) {
				this->velocity.z *= BOUNCE_MULTIPLIER;
			}
			this->velocity.x += gravity->x;
			this->velocity.y += gravity->y;
			this->velocity.z += gravity->z;
			this->center += this->velocity;
		}
		else {
			this->center += this->velocity;
			this->velocity.x += gravity->x;
			this->velocity.y += gravity->y;
			this->velocity.z += gravity->z;
		}
		//check for particle-particle collisions
		particleCollision(particles);


	}

	bool boundsCollision(Container* c, int* x, int* y, int* z) {
		bool ret = false;
		if (this->center.x + this->radius > c->width / 2 - FLT_EPSILON) {
			this->center.x = c->width / 2 - FLT_EPSILON - this->radius;
			*x = 1;
			ret = true;
		}
		if (this->center.x - this->radius < -(c->width / 2) + FLT_EPSILON) {
			this->center.x = -(c->width / 2) + FLT_EPSILON + this->radius;
			*x = 1;
			ret = true;
		}
		if (this->center.y + this->radius > c->height / 2 - FLT_EPSILON) {
			this->center.y = c->height / 2 - FLT_EPSILON - this->radius;
			*y = 1;
			ret = true;
		}
		if (this->center.y - this->radius < -(c->height / 2) + FLT_EPSILON) {
			this->center.y = -(c->height / 2) + FLT_EPSILON + this->radius;
			*y = 1;
			ret = true;
		}
		if (this->center.z + this->radius > c->length / 2 - FLT_EPSILON) {
			this->center.z = c->length / 2 - FLT_EPSILON - this->radius;
			*z = 1;
			ret = true;
		}
		if (this->center.z - this->radius < -(c->length / 2) + FLT_EPSILON) {
			this->center.z = -(c->length / 2) + FLT_EPSILON + this->radius;
			*z = 1;
			ret = true;
		}
		return ret;
	}

	void particleCollision(vector<Particle*>* particles) {
		for (Particle* p : (*particles)) {
			float distsquared = (this->center.x - p->center.x) * (this->center.x - p->center.x) +
				(this->center.y - p->center.y) * (this->center.y - p->center.y) +
				(this->center.z - p->center.z) * (this->center.z - p->center.z);
			float rsqaured = (this->radius + p->radius) * (this->radius + p->radius);

			if (p != this && distsquared < rsqaured - FLT_EPSILON) {
				Vector3D v_delta = this->velocity - p->velocity;
				Vector3D n = this->center - p->center;
				n.normalize();
				float lambda = -dot(n, v_delta);
				Vector3D impulse = lambda * n;
				this->velocity += impulse;
				p->velocity -= impulse;
				this->center = p->center + n * (this->radius + p->radius + FLT_EPSILON);
			}
			if (p != this && distsquared < (rsqaured * INFLUENCE_MULTIPLIER) - FLT_EPSILON) {
				Vector3D n = this->center - p->center;
				n.normalize();
				Vector3D grav_force = n * (GRAV_CONST / rsqaured) * GRAV_MULTIPLIER;
				this->velocity += grav_force;
				p->velocity -= grav_force;
			}
		}
	}
};

#endif // !PARTICLE_H


#pragma once