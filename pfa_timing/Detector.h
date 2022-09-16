#ifndef DETECTOR_H
#define DETECTOR_H

#include "includes.h"

#include "Cells.C"


class Detector
{
	
public:
	Detector();
	~Detector();

	void AddCell(Cell* newCell,int xId, int yId, int zId,  bool hited = false);
	/*
		add a cell in the detector
	*/

	int GetNumberOfHits();
	/*
		Return the number of hits
	*/
	int GetSize(int i);
	/*
		Return the size of the detector (96, 96, 48)
	*/
	Cell* GetCell(int xId, int yId, int zId);
	/*
		Return the cell in location 
	*/

	std::vector<Cell*> GetHits(); 
	/*
		Return the a vector of all cells hited
	*/

private:
	Cell* cells[96][96][48];		//map of all particles
	std::vector<Cell*> cellsHited;	//list of all hits
	int detectorSize[3];
};


#endif