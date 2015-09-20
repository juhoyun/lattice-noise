#ifndef _ProximityFunctions_H_
#define _ProximityFunctions_H_

#include "PseudoRandom.h"

namespace GenericLatticeNoiseAlgorithm
{
	float Gradient(PseudoRandom& pRnd, float x, float y, float lambda, float xHash, float yHash);
	float Hill(PseudoRandom& pRnd, float x, float y, float lambda, float xHash, float yHash);
	float HillAndSlope(PseudoRandom& pRnd, float x, float y, float lambda, float xHash, float yHash);
	float Curvature(PseudoRandom& pRnd, float x, float y, float lambda, float xHash, float yHash);
	float GradientRigged(PseudoRandom& pRnd, float x, float y, float lambda, float xHash, float yHash);
	float GradientSemiRigged(PseudoRandom& pRnd, float x, float y, float lambda, float xHash, float yHash);
	float GradientSemiRiggedIncreased(PseudoRandom& pRnd, float x, float y, float lambda, float xHash, float yHash);
	float GradientBillowed(PseudoRandom& pRnd, float x, float y, float lambda, float xHash, float yHash);
	float CurvatureRigged(PseudoRandom& pRnd, float x, float y, float lambda, float xHash, float yHash);
	float CurvatureBillowed(PseudoRandom& pRnd, float x, float y, float lambda, float xHash, float yHash);
	float CurvatureSemiRigged(PseudoRandom& pRnd, float x, float y, float lambda, float xHash, float yHash);
	float CurvatureSemiRiggedDisplaced(PseudoRandom& pRnd, float x, float y, float lambda, float xHash, float yHash);
	float DoublePlain(PseudoRandom& pRnd, float x, float y, float lambda, float xHash, float yHash);
	float VerticalEdge(PseudoRandom& pRnd, float x, float y, float lambda, float xHash, float yHash);
	float VerticalEdgeInverse(PseudoRandom& pRnd, float x, float y, float lambda, float xHash, float yHash);
	float TriangularEdgeOne(PseudoRandom& pRnd, float x, float y, float lambda, float xHash, float yHash);
	float TriangularEdgeTwo(PseudoRandom& pRnd, float x, float y, float lambda, float xHash, float yHash);
	float MonkeySaddle(PseudoRandom& pRnd, float x, float y, float lambda, float xHash, float yHash);
	float Hiperbolic(PseudoRandom& pRnd, float x, float y, float lambda, float xHash, float yHash);
	float HiperbolicPlains(PseudoRandom& pRnd, float x, float y, float lambda, float xHash, float yHash);
	float HiperbolicPlainsDisplaced(PseudoRandom& pRnd, float x, float y, float lambda, float xHash, float yHash);
	float VerticalEdgeDisplaced(PseudoRandom& pRnd, float x, float y, float lambda, float xHash, float yHash);
	float VerticalEdgeInverseDisplaced(PseudoRandom& pRnd, float x, float y, float lambda, float xHash, float yHash);
	float Parabolic(PseudoRandom& pRnd, float x, float y, float lambda, float xHash, float yHash);
	float ParabolicInverse(PseudoRandom& pRnd, float x, float y, float lambda, float xHash, float yHash);
	float ParabolicComposed(PseudoRandom& pRnd, float x, float y, float lambda, float xHash, float yHash);
	float ParabolicComposedII(PseudoRandom& pRnd, float x, float y, float lambda, float xHash, float yHash);
	float DiplacedParabole(PseudoRandom& pRnd, float x, float y, float lambda, float xHash, float yHash);
};

#endif /* _ProximityFunctions_H_ */