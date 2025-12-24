/*
    SymIntegration is branching from SymbolicC++ 3.35
    SymbolicC++ : An object oriented computer algebra system written in C++

    Copyright (C) 2008 Yorick Hardy and Willi-Hans Steeb

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along
    with this program; if not, write to the Free Software Foundation, Inc.,
    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
*/
// THANKS SENTINEL!!! and Freya too

// g++ -o result main.cpp -lsymintegration

#include "symintegrationc++.h"
#include <bits/stdc++.h>
#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>

//#include <glm/glm.hpp>
//#include <glm/gtc/matrix_transform.hpp>
//#include <glm/gtc/type_ptr.hpp>

#define π 3.1415926535897f
#include <chrono>

using namespace std::chrono;
using namespace std;
using namespace SymbolicConstant;


int main() {
	// Get starting timepoint
	auto start = high_resolution_clock::now();

	/*glm::mat4 transrotate = glm::mat4(1.0f);
	for (int i = 0; i < 4; ++i) 
	{
		for (int j=0; j<4; ++j) 
		{
			cout << std::fixed << setw(23) << setprecision(6) << transrotate[i][j] << setw(23);
		}
		cout << endl;
	}*/
	dvec x = loadVectorFromFile("vector.txt");
	cout << "\nVector x:"<< endl;
	printVector(x);

	dvec reflection = reflection_yaxis(x);
	cout << "\nReflection about the y-axis:"<< endl;
	printVector(reflection);
	
	dvec reflection2 = reflection_xaxis(x);
	cout << "\nReflection about the x-axis:"<< endl;
	printVector(reflection2);
	
	dvec reflection3 = reflection_linex(x);
	cout << "\nReflection about the line y=x:"<< endl;
	printVector(reflection3);
	
	dvec op1 = orthogonalprojection_xaxis(x);
	cout << "\nOrthogonal projection on the x-axis:"<< endl;
	printVector(op1);

	dvec op2 = orthogonalprojection_yaxis(x);
	cout << "\nOrthogonal projection on the y-axis:"<< endl;
	printVector(op2);

	dvec rotation1 = rotation_ccw(x, π/6);
	cout << "\nRotation through an angle π/6 :"<< endl;
	printVector(rotation1);
	
	dvec x3d = loadVectorFromFile("vector3d.txt");
	dvec x3d1 = loadVectorFromFile("vector3d1.txt");
	cout << "\nVector x in 3d:"<< endl;
	printVector(x3d);
	
	dvec reflection4 = reflection_xyplane(x3d);
	cout << "\nReflection about the xy-plane:"<< endl;
	printVector(reflection4);
	
	dvec reflection5 = reflection_xzplane(x3d);
	cout << "\nReflection about the xz-plane:"<< endl;
	printVector(reflection5);
	
	dvec reflection6 = reflection_yzplane(x3d);
	cout << "\nReflection about the yz-plane:"<< endl;
	printVector(reflection6);

	dvec op3 = orthogonalprojection_xyplane(x3d);
	cout << "\nOrthogonal projection on the xy-plane:"<< endl;
	printVector(op3);

	dvec op4 = orthogonalprojection_xzplane(x3d);
	cout << "\nOrthogonal projection on the xz-plane:"<< endl;
	printVector(op4);

	dvec op5 = orthogonalprojection_yzplane(x3d);
	cout << "\nOrthogonal projection on the yz-plane:"<< endl;
	printVector(op5);

	dvec rotation2 = rotation3d_xaxis_ccw(x3d, π/6);
	cout << "\nRotation about the positive x-axis through an angle π/6 :"<< endl;
	printVector(rotation2);

	dvec rotation5 = rotation3d_xaxis_ccw(x3d1, π/6);
	cout << "\nRotation about the negative x-axis through an angle π/6 :"<< endl;
	printVector(rotation5);

	dvec rotation3 = rotation3d_yaxis_ccw(x3d, π/6);
	cout << "\nRotation about the positive y-axis through an angle π/6 :"<< endl;
	printVector(rotation3);

	dvec rotation4 = rotation3d_zaxis_ccw(x3d, π/6);
	cout << "\nRotation about the positive z-axis through an angle π/6 :"<< endl;
	printVector(rotation4);

	dvec compress1 = compressionexpansion3d_xdirection(x3d, 0.3);
	cout << "\nCompression in the x-direction with factor 0.3:"<< endl;
	printVector(compress1);

	dvec expansion1 = compressionexpansion3d_xdirection(x3d, 2.3);
	cout << "\nExpansion in the x-direction with factor 2.3:"<< endl;
	printVector(expansion1);

	dvec dilation1 = contractiondilation(x, 2.3);
	cout << "\nDilation in 2d with factor 2.3:"<< endl;
	printVector(dilation1);

	dvec dilation2 = contractiondilation(x3d, 2.3);
	cout << "\nDilation in 3d with factor 2.3:"<< endl;
	printVector(dilation2);

	dvec shear1 = shear_xdirection(x, 3.3);
	cout << "\nShear in 2d in the x-direction with factor 3.3:"<< endl;
	printVector(shear1);

	dvec shear2 = shear_ydirection(x, 3.3);
	cout << "\nShear in 2d in the y-direction with factor 3.3:"<< endl;
	printVector(shear2);
	
	// Get ending timepoint
	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<microseconds>(stop - start);

	cout << "\nTime taken by function: " << duration.count() << " microseconds" << endl;

	return 0;
}
