#include <cstdlib>
#include "PseudoRandom.h"

namespace GenericLatticeNoiseAlgorithm
{
	void PseudoRandom::Initiate(int s)
	{
		useCSRandom = true;
		scale = 1.f;
		//initiate c# pseudo random
		srand(s);
		for (int i=0; i<251; i++)
		{
			for (int j=0; j<251; j++)
			{
				for (int k = 0; k < 5; k++)
				{
					//hashField[i][j][k] = ((float)rnd.Next(0, 1000000) - 500000f) / 1000000f; //prn ranging -1 to +1
					int rndValue = rand() - RAND_MAX/2;
					hashField[i][j][k] = ((float)rndValue * 2.0f) / RAND_MAX;
				}
			}
		}
		//initiate custom prn hash
		seed = s;
		seed = (40000 + s % 100000) * 4000;
		seedA = seed + 43876431;
		seedB = seed + 81256937;
		seedC = seed + 124532173;
		seedD = seed + 159683467;
	}

	float PseudoRandom::HashRandom(float x, float y, int index)
	{
		x *= scale; x += 25100;
		y *= scale; y += 25100;
		if (useCSRandom == true)
		{
			return hashField[(int)x % 251][(int)y % 251][index];
		}
		else
		{
			return HashF((int)x, (int)y, index);
		}
	}

	int PseudoRandom::Hash(int i, int j, int index) //Returns a integer 0 to 2147483647
	{
		int a = i;
		int b = j;
		for (int r = 0; r < 3; r++)
		{
			a = Rotate((a ^ seedA) + (b ^ seedC), 25 - index - index);
			b = Rotate((a ^ seedB) + (b ^ seedD), 3 + index + index);
		}
		return a ^ b;
	}
}
