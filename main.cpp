#include <stdio.h>
#include <math.h>
#include <cstdlib>
#include "PseudoRandom.h"
#include "ProximityFunctions.h"

#define MAP_SIZE	512

float heightmap[MAP_SIZE+1][MAP_SIZE+1];

using namespace GenericLatticeNoiseAlgorithm;

typedef float (*ProximityFunction)(PseudoRandom& pRnd, float x, float y, float lambda, float xHash, float yHash);
typedef float (*FadeFunction)(float x, float y);

ProximityFunction ProximityFunctions[] =
{
	Gradient,						// 0 - Perlin
	Hill,							// 1
	HillAndSlope,
	Curvature,
	GradientRigged,
	GradientSemiRigged,				// 5
	GradientSemiRiggedIncreased,
	GradientBillowed,
	CurvatureRigged,
	CurvatureBillowed,
	CurvatureSemiRigged,			// 10
	CurvatureSemiRiggedDisplaced,
	DoublePlain,
	VerticalEdge,
	VerticalEdgeInverse,
	TriangularEdgeOne,				// 15
	TriangularEdgeTwo,
	MonkeySaddle,
	Hiperbolic,
	HiperbolicPlains,
	HiperbolicPlainsDisplaced,		// 20
	VerticalEdgeDisplaced,
	VerticalEdgeInverseDisplaced,
	Parabolic,
	ParabolicInverse,
	ParabolicComposed,				// 25
	ParabolicComposedII,
	DiplacedParabole
};

const char *ProximityFunctionNames[] =
{
	"Gradient",						// 0 - Perlin
	"Hill",							// 1
	"HillAndSlope",
	"Curvature",
	"GradientRigged",
	"GradientSemiRigged",			// 5
	"GradientSemiRiggedIncreased",
	"GradientBillowed",
	"CurvatureRigged",
	"CurvatureBillowed",
	"CurvatureSemiRigged",			// 10
	"CurvatureSemiRiggedDisplaced",
	"DoublePlain",
	"VerticalEdge",
	"VerticalEdgeInverse",
	"TriangularEdgeOne",			// 15
	"TriangularEdgeTwo",
	"MonkeySaddle",
	"Hiperbolic",
	"HiperbolicPlains",
	"HiperbolicPlainsDisplaced",	// 20
	"VerticalEdgeDisplaced",
	"VerticalEdgeInverseDisplaced",
	"Parabolic",
	"ParabolicInverse",
	"ParabolicComposed",			// 25
	"ParabolicComposedII",
	"DiplacedParabole"
};

float Hermite(float x, float y)
{
    return (1 - (x * x * (3 - 2 * x))) * (1 - (y * y * (3 - 2 * y)));
}

int main(int n, char *av[])
{
	PseudoRandom pRnd;
	int seed;
	int funcIndex;
	if (n == 1)
		seed = 100;
	else
	{
		seed = atoi(av[1]);
	}

	if (n > 2)
		funcIndex = atoi(av[2]);
	else
		funcIndex = 0;
	pRnd.Initiate(seed);
	int size = MAP_SIZE;
	int numberOfOctaves = 3;
    int initialStep = size / 4;
    FadeFunction Fade = Hermite;
    ProximityFunction Proximity = ProximityFunctions[funcIndex];
	pRnd.SetScale(1024.f / size);

	for (int x = 0; x < size + 1; x += 1)
	{
		for (int y = 0; y < size + 1; y += 1)
		{
			heightmap[x][y] = 0.5f;
		}
	}

	float minH = 10.0f;
	float maxH = -10.0f;

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
					if (y == size - sizeGrid) 
					{ 
						gridV = sizeGrid + 1; 
					} else 
					{ 
						gridV = sizeGrid; 
					}

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

								float h = heightmap[x + (int)u][y + (int)v] += hBase * persistence;

								if (h > maxH)
									maxH = h;
								if (h < minH)
									minH = h;
							}
						}
					}
				}
			}
		}
	}
	printf("Min: %f, Max: %f\n", minH, maxH);

	char filename[1000];
	sprintf(filename, "%s_%d.raw", ProximityFunctionNames[funcIndex], seed);
	FILE* fp = fopen(filename, "wb");
	if (!fp)
	{
		printf("File Open Error: %s\n", filename);
		return -1;
	}
	for (int x = 0; x < size + 1; ++x)
	{
		for (int y = 0; y < size + 1; ++y)
		{
			unsigned short color = (unsigned short)((heightmap[x][y] - 0.3f) * 65535.0f / 0.4f);
			fwrite(&color, 1, 2, fp);
		}
	}
	fclose(fp);
	printf("Saved %s\n", filename);

	return 0;
}
