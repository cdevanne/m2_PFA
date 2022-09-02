#ifndef TOOLSMATH_C
#define TOOLSMATH_C

#include "ToolsMath.h"

std::vector<double> CrossProduct(std::vector<double> const &a, std::vector<double> const &b)
{
	std::vector<double> r (a.size());  
	r[0] = a[1]*b[2]-a[2]*b[1];
	r[1] = a[2]*b[0]-a[0]*b[2];
	r[2] = a[0]*b[1]-a[1]*b[0];
	return r;
}

double Norm(double v1[3])
{
	double norm = 0;

	for (int i = 0; i < 3; ++i)
	{
		norm += v1[i] * v1[i];
	}

	return sqrt(norm);
}

double Norm(std::vector<double>  const &v1)
{
	double norm = 0;

	for (int i = 0; i < 3; ++i)
	{
		norm += v1[i] * v1[i];
	}

	return sqrt(norm);
}

double MinimalDistanceLines(double vector1[3], double vector2[3], double point1[3], double point2[3])
{
	//distance minimum entre deux droites
	double distance;
	distance = 
		abs(
			(
				(point2[0] - point1[0]) * ((vector1[1] * vector2[2]) - (vector1[2] * vector2[1]))
				+(point2[1] - point1[1]) * ((vector1[2] * vector2[0]) - (vector1[0] * vector2[2]))
				+(point2[2] - point1[2]) * ((vector1[0] * vector2[1]) - (vector1[1] * vector2[0]))
			)
			/
			(
				sqrt(pow(((vector1[1] * vector2[2]) - (vector1[2] * vector2[1])), 2.)
					+pow (((vector1[2] * vector2[0]) - (vector1[0] * vector2[2])), 2.)
					+pow (((vector1[0] * vector2[1]) - (vector1[1] * vector2[0])), 2.))
			)
		);
		
	return distance;
}

double CosAngleBetweenVectors(double vector1[3], double vector2[3])
{
	//cos(angle)
	double product = 0;
	double norm1 = 0;
	double norm2 = 0;

	for (int i = 0; i < 3; i++)
	{
		product += vector1[i] * vector2[i];
		norm1 += vector1[i] * vector1[i];
		norm2 += vector2[i] * vector2[i];				
	}

	double cos = product / (sqrt(norm1) * sqrt(norm2));

	return cos;
}

double MinimalDistanceLinePoint(double vectorA[3], double vectorB[3], double A[3], double B[3], double distanceMin, double error)
{
	//trouve la coordonnée z ou la distance entres deux lignes est minimal
	double normV1 = Norm(vectorA) + 0.0000001;
	
	double originalB[3];
	double originalA[3];
	double distance = 10000;
	double oldDistance = 10000000.;
	double sign = 1.;
	double step = 8.;
	std::vector<double> mainVector (3);
	std::vector<double> AB (3);

	for (int i = 0; i < 3; ++i)
	{
		originalB[i] = B[i];
		originalA[i] = A[i];
		mainVector[i] = vectorA[i];
	}

	do
	{
		if (oldDistance < distance)
		{
			sign *= -1.;
			step *= 1./2.;
		}
		oldDistance = distance;


		for (int i = 0; i < 3; ++i)
		{
			B[i] = B[i] += vectorB[i]*0.5*sign;
		}

		for (int i = 0; i < 3; ++i)
		{
			AB[i] = B[i] - A[i];
		}

		std::vector<double> ABxU(3);

		double norm;
		norm = Norm(ABxU) + 0.0000001;
		distance = norm/(normV1+ 0.0000001);
	
		
		

	}while(distance > distanceMin + error && step >= .5);

	double d1 = B[2] - originalB[2];

	return d1;
}



#endif