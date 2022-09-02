#ifndef HELPER_FOR_PROGRAMMING_C
#define HELPER_FOR_PROGRAMMING_C

#include "HelperForProgramming.h"





void optimalData(std::vector<Arbor*> arbors)
{

		int id1 = -1;
		int id2 = -1;

		int iteration = 0;

		while(id1 == -1 && iteration < arbors.size())
		{
			if (arbors[iteration]->GetMainParticleId() == 1)
			{
				id1 = iteration;
			}
			++iteration;
			
		}
		iteration = 0;

		while(id2 == -1 && iteration < arbors.size())
		{
			if (arbors[iteration]->GetMainParticleId() == 2)
			{
				id2 = iteration;
			}
			++iteration;
		}

		if (id2 == -1)
		{
			id2 = 1;
		}
		if (id1 == -1)
		{
			id2 = 0;
		}

		for (int i = 0; i < arbors.size(); ++i)
		{
			if (arbors[i]->GetMainParticleId() == 2)
			{
				arbors[i]->SetArborId(id2);
			}else
			{
				arbors[i]->SetArborId(id1);
			}
		}
}



void FeedSecondParticle(std::vector<Arbor*> arbors, double range, double timing) 
{
	if (arbors.size() == 1) return arbors;
	int id2 = arborsMaxSizeTop3(arbors)[1];
	int id1 = arborsMaxSizeTop3(arbors)[0];
	for (int i = 0; i < arbors.size(); ++i)
	{	
		
		double dt = abs(arbors[id2]->GetTiming() -arbors[i]->GetTiming());
		double dtMin = distanceSeed(arbors[i], arbors[id2]) * timing;
		if (i != id1 && arbors[id2]->Distance(arbors[i], 1., 1., 0.05) < range && dt > dtMin)
		{
			arbors[i]->SetArborId(id2);
		}
	}
}


#endif