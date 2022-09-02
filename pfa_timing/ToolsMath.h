#ifndef TOOLSMATH_H
#define TOOLSMATH_H

#include "includes.h"
#include "Tools.C"

std::vector<double> CrossProduct(std::vector<double> const &a, std::vector<double> const &b);
/*
	return cross product of two vector 3d
*/
double Norm(double v1[3]);
/*
	return the norm
*/

double Norm(std::vector<double>  const &v1);
/*
	return the norm
*/

double MinimalDistanceLines(double vector1[3], double vector2[3], double point1[3], double point2[3]);
/*
	return the distance minimal between two lines 
*/

double CosAngleBetweenVectors(double vector1[3], double vector2[3]);
/*
	return cosinus
*/

double MinimalDistanceLinePoint(double vectorA[3], double vectorB[3], double A[3], double B[3], double distanceMin, double error);
/*
	return the distance minimal between a line and a point
*/


#endif