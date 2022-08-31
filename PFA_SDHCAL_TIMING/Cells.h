#ifndef CELLS_H
#define CELLS_H


#include "includes.h"



class Cell
{
/*
cells of the detector, SDHCAL is 96*96*48
*/

public:

// Constructor

	Cell();


	Cell(int xId, int yId, int zId);
	/*
		Create cells not hited
	*/


	Cell(int xId, int yId, int zId, double x, double y, double z, int particleId, int thr, double timing);
	/*
		Create cells hited and save informations
	*/


	~Cell();


// Setter

	void SetHited(bool isHited);
	/*
		set attribute hited
	*/
	void SetClusterised(bool isClusterised);
	/*
		set attribute clusterised
	*/
	void SetArborId(int id);
	/*
		set attribute arborId
	*/


// Getter

	int* GetId();
	/*
		Return the xId, yId, zId of the cell
	*/

	double* GetPosition();
	/*
		Return the position in a int[3] of the cell
	*/

	int GetParticleId();
	/*
		Return real particle Id
	*/

	int GetThreshold();
	/*
		Return the threshold
	*/
	double GetTiming();
	/*
		Return the timing
	*/

	int GetArborId();
	/*
		Return the attribute arborId
	*/

	bool GetHited();
	/*
		Return the the attribute hited
	*/

	bool GetClusterised();
	/*
		Return the the attribute clusterised
	*/



// Method

	double Distance(Cell* cell,  double xCoef = 1., double yCoef = 1., double zCoef = 1.);
	/*
		Return the distance between two cells, default coef are cartesian coordinate system
	*/

private:

// Attribute

	int positionId[3];	//position on the cell the 48*48*98 cells disposition
	double position[3];	//real position
	int particle;
	int threshold;
	double timing;
	int arborId;

	bool hited;
	bool clusterised;
};



#endif