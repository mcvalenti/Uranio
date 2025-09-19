/*
 * mission.h
 *
 *  Created on: 30 mar. 2024
 *      Author: ceci
 */

#ifndef MISSION_H_
#define MISSION_H_

#include "LDVector.h"


class SpaceVehicle
{
	public:
		LDVector init_sv;
		float mass;
		SpaceVehicle();
		SpaceVehicle(LDVector& init_sv, float mass);
		//~SpaceVehicle();
		friend ostream& operator<<(ostream& os, const SpaceVehicle& sat);
		
};		


#endif /* MISSION_H_ */
