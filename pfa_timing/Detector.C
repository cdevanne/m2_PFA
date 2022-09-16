#ifndef DETECTOR_C
#define DETECTOR_C

#include "Detector.h"


Detector::Detector()
{
	detectorSize[0] = 96;
	detectorSize[1] = 96;
	detectorSize[2] = 48;

	cellsHited.clear();
}

Detector::~Detector() 
{
	for (int i = 0; i < GetSize(0); ++i)
	{
		for (int j = 0; j < GetSize(1); ++j)
		{
			for (int k = 0; k < GetSize(2); ++k)
			{
				delete this->GetCell(i, j, k);
			}
		}	
	}
}


void Detector::AddCell(Cell* newCell, int xId, int yId, int zId, bool hited = false)
{
	cells[xId][yId][zId] = newCell;
	if (hited == true) cellsHited.push_back(newCell);
}

Cell* Detector::GetCell(int xId, int yId, int zId)
{
	return cells[xId][yId][zId];
}

std::vector<Cell*> Detector::GetHits()
{
	return cellsHited;
}

int Detector::GetSize(int i)
{
	return detectorSize[i];
}

int Detector::GetNumberOfHits()
{
	return cellsHited.size();
}


#endif