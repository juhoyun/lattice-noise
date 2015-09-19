#include <stdio.h>
#include <math.h>
#include "PseudoRandom.h"

float heightmap[512][512];

using namespace GenericLatticeNoiseAlgorithm;

typedef float (*ProximityFunction)(PseudoRandom& pRnd, float x, float y, float lambda, float xHash, float yHash);
typedef float (*FadeFunction)(float x, float y);

float Gradient(PseudoRandom& pRnd, float x, float y, float lambda, float xHash, float yHash)
{
    return (x * pRnd.HashRandom(xHash, yHash, 0) + y * pRnd.HashRandom(xHash, yHash, 1)) / (lambda);
}

float Hermite(float x, float y)
{
    return (1 - (x * x * (3 - 2 * x))) * (1 - (y * y * (3 - 2 * y)));
}

int main()
{
	PseudoRandom pRnd;
	int size = 512;
	int numberOfOctaves = 3;
    int initialStep = size / 4;
    FadeFunction Fade = Hermite;
    ProximityFunction Proximity = Gradient;
	pRnd.SetScale(1024.f / size);

	for (int x = 0; x < size + 1; x += 1)
	{
		for (int y = 0; y < size + 1; y += 1)
		{
			heightmap[x][y] = 0.5f;
		}
	}

	for (int octave = 0; octave < numberOfOctaves; octave++)
	{
		int sizeGrid;
		float sizeGridF;
		int init = 0;
		sizeGrid = initialStep / (1 << (octave));
		init = 0;
		//if (rb_G01.IsChecked == true) { init = 0; }
		//if (rb_G02.IsChecked == true) { init = - sizeGrid / 2; }
		//if (rb_G03.IsChecked == true) { init = - sizeGrid / 4; }

		sizeGridF = (float)sizeGrid;
		float relativeSize = sizeGridF / size;
		//float persistence = (float)pow(relativeSize, 1.0 - Persistence_slider.Value);
		float persistence = relativeSize;

		int gridU, gridV;
		float hBase;                          

		if (sizeGrid >= 1)
		{
			for (int x = init; x < size + 1; x += sizeGrid)
			{
				if (x == size - sizeGrid)
				{
					gridU = sizeGrid + 1;
				} else 
				{ 
					gridU = sizeGrid; 
				}
				for (int y = init; y < size + 1; y += sizeGrid)
				{
					if (y == size - sizeGrid) { gridV = sizeGrid + 1; } else { gridV = sizeGrid; }

					for (float u = 0; u < gridU; u++)
					{
						for (float v = 0; v < gridV; v++)
						{
							if ((x + (int)u > size - 1 || x + (int)u < 0 || y + (int)v > size -1 || y + (int)v < 0) == false)
							{
								float us = u / sizeGridF;
								float vs = v / sizeGridF;

								hBase =
											Proximity(pRnd, u, v, sizeGridF, (float)x, (float)y)
											* Fade(u / sizeGridF, v / sizeGridF)

											+ Proximity(pRnd, u - sizeGridF, v, sizeGridF, (float)x + sizeGridF, (float)y)
											* Fade(1 - u / sizeGridF, v / sizeGridF)

											+ Proximity(pRnd, u, v - sizeGridF, sizeGridF, (float)x, (float)y + sizeGridF)
											* Fade(u / sizeGridF, 1 - v / sizeGridF)

											+ Proximity(pRnd, u - sizeGridF, v - sizeGridF, sizeGridF, (float)x + sizeGridF, (float)y + sizeGridF)
											* Fade(1 - u / sizeGridF, 1 - v / sizeGridF);

								heightmap[x + (int)u][y + (int)v] += hBase * persistence;
							}
						}
					}
				}
			}
		}
	}

	return 0;
}
