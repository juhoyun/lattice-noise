#include <math.h>
#include <algorithm>
#include "PseudoRandom.h"

namespace GenericLatticeNoiseAlgorithm
{
	float Gradient(PseudoRandom& pRnd, float x, float y, float lambda, float xHash, float yHash)
	{
		return (x * pRnd.HashRandom(xHash, yHash, 0) + y * pRnd.HashRandom(xHash, yHash, 1)) / (lambda);
	}

	float Hill(PseudoRandom& pRnd, float x, float y, float lambda, float xHash, float yHash)
	{
		return pRnd.HashRandom(xHash, yHash, 0) * 0.5f ;
	}

	float HillAndSlope(PseudoRandom& pRnd, float x, float y, float lambda, float xHash, float yHash)
	{
		return pRnd.HashRandom(xHash, yHash, 2) * 0.5f
			+ std::min(0.f, x * pRnd.HashRandom(xHash, yHash, 0) + y * pRnd.HashRandom(xHash, yHash, 1))
			/ lambda;
	}

	float Curvature(PseudoRandom& pRnd, float x, float y, float lambda, float xHash, float yHash)
	{
		return  
			(x * x * pRnd.HashRandom(xHash, yHash, 0) + y * y * pRnd.HashRandom(xHash, yHash, 1) + 2.f * x * y * pRnd.HashRandom(xHash, yHash, 2)) /
			(lambda * lambda);
	}

	float GradientRigged(PseudoRandom& pRnd, float x, float y, float lambda, float xHash, float yHash)
	{
		return  
			(-fabs(x * pRnd.HashRandom(xHash, yHash, 0) + y * pRnd.HashRandom(xHash, yHash, 1))) / (lambda);
	}

	float GradientSemiRigged(PseudoRandom& pRnd, float x, float y, float lambda, float xHash, float yHash)
	{
		return 
			(std::min(0.f, x * pRnd.HashRandom(xHash, yHash, 0) + y * pRnd.HashRandom(xHash, yHash, 1))) / (lambda);
	}

	float GradientSemiRiggedIncreased(PseudoRandom& pRnd, float x, float y, float lambda, float xHash, float yHash)
	{
		return
			(std::min(0.f, x * pRnd.HashRandom(xHash, yHash, 0) + y * pRnd.HashRandom(xHash, yHash, 1) - pRnd.HashRandom(xHash, yHash, 2))) / 
			(lambda);
	}

	float GradientBillowed(PseudoRandom& pRnd, float x, float y, float lambda, float xHash, float yHash)
	{
		return 
			(fabs(x * pRnd.HashRandom(xHash, yHash, 0) + y * pRnd.HashRandom(xHash, yHash, 1))) / (lambda);
	}

	float CurvatureRigged(PseudoRandom& pRnd, float x, float y, float lambda, float xHash, float yHash)
	{
		return 
			(-fabs(x * x * pRnd.HashRandom(xHash, yHash, 0) + y * y * pRnd.HashRandom(xHash, yHash, 1) + 2.f * x * y * pRnd.HashRandom(xHash, yHash, 2))) /
			(lambda * lambda);
	}

	float CurvatureBillowed(PseudoRandom& pRnd, float x, float y, float lambda, float xHash, float yHash)
	{
		return 
			(fabs(x * x * pRnd.HashRandom(xHash, yHash, 0) + y * y * pRnd.HashRandom(xHash, yHash, 1) + 2.f * x * y * pRnd.HashRandom(xHash, yHash, 2))) /
			(lambda * lambda);
	}

	float CurvatureSemiRigged(PseudoRandom& pRnd, float x, float y, float lambda, float xHash, float yHash)
	{
		return 
			(std::min(0.f, x * x * pRnd.HashRandom(xHash, yHash, 0) + y * y * pRnd.HashRandom(xHash, yHash, 1) + 2.f * x * y * pRnd.HashRandom(xHash, yHash, 2))) /
			(lambda * lambda);
	}

	float CurvatureSemiRiggedDisplaced(PseudoRandom& pRnd, float x, float y, float lambda, float xHash, float yHash)
	{
		return 
			(std::min(
			0.f,
			(x + pRnd.HashRandom(xHash, yHash, 1) * 0.75f * lambda) * (x + pRnd.HashRandom(xHash, yHash, 1) * 0.75f * lambda) * pRnd.HashRandom(xHash, yHash, 0) 
			+ (y + pRnd.HashRandom(xHash, yHash, 0) * 0.75f * lambda) * (y + pRnd.HashRandom(xHash, yHash, 0) * 0.75f * lambda) * pRnd.HashRandom(xHash, yHash, 1)
			+ 2.f * (x + pRnd.HashRandom(xHash, yHash, 1) * 0.75f * lambda) * (y + pRnd.HashRandom(xHash, yHash, 0) * 0.75f * lambda) * pRnd.HashRandom(xHash, yHash, 2))) /
			(lambda * lambda);
	}

	float DoublePlain(PseudoRandom& pRnd, float x, float y, float lambda, float xHash, float yHash)
	{
		return 
			(std::min(
			x * pRnd.HashRandom(xHash, yHash, 0) + y * pRnd.HashRandom(xHash, yHash, 1),
			x * pRnd.HashRandom(xHash, yHash, 1) + y * pRnd.HashRandom(xHash, yHash, 2))) 
			/ lambda;
	}

	float VerticalEdge(PseudoRandom& pRnd, float x, float y, float lambda, float xHash, float yHash)
	{
		return 
			(std::min(
			x * pRnd.HashRandom(xHash, yHash, 0) + y * pRnd.HashRandom(xHash, yHash, 1) + pRnd.HashRandom(xHash, yHash, 2),
			-x * pRnd.HashRandom(xHash, yHash, 0) - y * pRnd.HashRandom(xHash, yHash, 1) + pRnd.HashRandom(xHash, yHash, 2))) 
			/ lambda;
	}

	float VerticalEdgeInverse(PseudoRandom& pRnd, float x, float y, float lambda, float xHash, float yHash)
	{
		return 
			(std::max(
			x * pRnd.HashRandom(xHash, yHash, 0) + y * pRnd.HashRandom(xHash, yHash, 1) + pRnd.HashRandom(xHash, yHash, 2),
			-x * pRnd.HashRandom(xHash, yHash, 0) - y * pRnd.HashRandom(xHash, yHash, 1) + pRnd.HashRandom(xHash, yHash, 2))) 
			/ lambda;
	}

	float TriangularEdgeOne(PseudoRandom& pRnd, float x, float y, float lambda, float xHash, float yHash)
	{
		return 
			(std::min(
			std::max(-x * pRnd.HashRandom(xHash, yHash, 0) - y * pRnd.HashRandom(xHash, yHash, 1), 
			-x * pRnd.HashRandom(xHash, yHash, 1) + y * pRnd.HashRandom(xHash, yHash, 0)),
			std::max(x * pRnd.HashRandom(xHash, yHash, 0) + y * pRnd.HashRandom(xHash, yHash, 1), 
			x * pRnd.HashRandom(xHash, yHash, 1) - y * pRnd.HashRandom(xHash, yHash, 0)))) 
			/ lambda;
	}

	float TriangularEdgeTwo(PseudoRandom& pRnd, float x, float y, float lambda, float xHash, float yHash)
	{
		return 
			(std::min(
			std::min(std::max(-x * pRnd.HashRandom(xHash, yHash, 0) - y * pRnd.HashRandom(xHash, yHash, 1), -x * pRnd.HashRandom(xHash, yHash, 1) + y * pRnd.HashRandom(xHash, yHash, 0)), std::max(x * pRnd.HashRandom(xHash, yHash, 0) + y * pRnd.HashRandom(xHash, yHash, 1), x * pRnd.HashRandom(xHash, yHash, 1) - y * pRnd.HashRandom(xHash, yHash, 0))),
			std::min(std::max(-x * pRnd.HashRandom(xHash, yHash, 0) - y * pRnd.HashRandom(xHash, yHash, 1), -x * pRnd.HashRandom(xHash, yHash, 1) + y * pRnd.HashRandom(xHash, yHash, 0)), std::max(x * pRnd.HashRandom(xHash, yHash, 0) + y * pRnd.HashRandom(xHash, yHash, 1), -x * pRnd.HashRandom(xHash, yHash, 1) + y * pRnd.HashRandom(xHash, yHash, 0))))) 
			/ (lambda);
	}

	float MonkeySaddle(PseudoRandom& pRnd, float x, float y, float lambda, float xHash, float yHash)
	{
		return
			(pRnd.HashRandom(xHash, yHash, 0) * (x * x * x - 3 * x * y * y)) 
			/ (lambda * lambda * lambda);
	}

	float Hiperbolic(PseudoRandom& pRnd, float x, float y, float lambda, float xHash, float yHash)
	{
		return
			(1.f / (30.f + fabs(x * pRnd.HashRandom(xHash, yHash, 0) + y * pRnd.HashRandom(xHash, yHash, 1))))
			/ lambda;
	}

	float HiperbolicPlains(PseudoRandom& pRnd, float x, float y, float lambda, float xHash, float yHash)
	{
		return
			(std::min(
			std::max((-x * pRnd.HashRandom(xHash, yHash, 0) - y * pRnd.HashRandom(xHash, yHash, 1)) * (-x * pRnd.HashRandom(xHash, yHash, 0) - y * pRnd.HashRandom(xHash, yHash, 1)) * (-x * pRnd.HashRandom(xHash, yHash, 0) - y * pRnd.HashRandom(xHash, yHash, 1)),
			(-x * pRnd.HashRandom(xHash, yHash, 1) + y * pRnd.HashRandom(xHash, yHash, 0)) * (-x * pRnd.HashRandom(xHash, yHash, 1) + y * pRnd.HashRandom(xHash, yHash, 0)) * (-x * pRnd.HashRandom(xHash, yHash, 1) + y * pRnd.HashRandom(xHash, yHash, 0))),
			std::max((x * pRnd.HashRandom(xHash, yHash, 0) + y * pRnd.HashRandom(xHash, yHash, 1)) * (x * pRnd.HashRandom(xHash, yHash, 0) + y * pRnd.HashRandom(xHash, yHash, 1)) * (x * pRnd.HashRandom(xHash, yHash, 0) + y * pRnd.HashRandom(xHash, yHash, 1)),
			(x * pRnd.HashRandom(xHash, yHash, 1) - y * pRnd.HashRandom(xHash, yHash, 0)) * (x * pRnd.HashRandom(xHash, yHash, 1) - y * pRnd.HashRandom(xHash, yHash, 0)) * (x * pRnd.HashRandom(xHash, yHash, 1) - y * pRnd.HashRandom(xHash, yHash, 0)))))
			/ (lambda * lambda * lambda);
	}

	float HiperbolicPlainsDisplaced(PseudoRandom& pRnd, float x, float y, float lambda, float xHash, float yHash)
	{
		x += lambda * pRnd.HashRandom(xHash, yHash, 1) * 0.25f;
		y += lambda * pRnd.HashRandom(xHash, yHash, 2) * 0.25f;
		return 
			(std::min(
			std::max((-x * pRnd.HashRandom(xHash, yHash, 0) - y * pRnd.HashRandom(xHash, yHash, 1)) * (-x * pRnd.HashRandom(xHash, yHash, 0) - y * pRnd.HashRandom(xHash, yHash, 1)) * (-x * pRnd.HashRandom(xHash, yHash, 0) - y * pRnd.HashRandom(xHash, yHash, 1)),
			(-x * pRnd.HashRandom(xHash, yHash, 1) + y * pRnd.HashRandom(xHash, yHash, 0)) * (-x * pRnd.HashRandom(xHash, yHash, 1) + y * pRnd.HashRandom(xHash, yHash, 0)) * (-x * pRnd.HashRandom(xHash, yHash, 1) + y * pRnd.HashRandom(xHash, yHash, 0))),
			std::max((x * pRnd.HashRandom(xHash, yHash, 0) + y * pRnd.HashRandom(xHash, yHash, 1)) * (x * pRnd.HashRandom(xHash, yHash, 0) + y * pRnd.HashRandom(xHash, yHash, 1)) * (x * pRnd.HashRandom(xHash, yHash, 0) + y * pRnd.HashRandom(xHash, yHash, 1)),
			(x * pRnd.HashRandom(xHash, yHash, 1) - y * pRnd.HashRandom(xHash, yHash, 0)) * (x * pRnd.HashRandom(xHash, yHash, 1) - y * pRnd.HashRandom(xHash, yHash, 0)) * (x * pRnd.HashRandom(xHash, yHash, 1) - y * pRnd.HashRandom(xHash, yHash, 0)))))
			/ (lambda * lambda * lambda);
	}

	float VerticalEdgeDisplaced(PseudoRandom& pRnd, float x, float y, float lambda, float xHash, float yHash)
	{
		return 
			(std::min(
			((x + 0.5f * lambda) * pRnd.HashRandom(xHash, yHash, 1)) * pRnd.HashRandom(xHash, yHash, 0) + ((y + 0.5f * lambda) * pRnd.HashRandom(xHash, yHash, 0)) * pRnd.HashRandom(xHash, yHash,
			1) + pRnd.HashRandom(xHash, yHash, 2), -((x + 0.5f * lambda) * pRnd.HashRandom(xHash, yHash, 1)) * pRnd.HashRandom(xHash, yHash, 0) - ((y + 0.5f * lambda) * pRnd.HashRandom(xHash, yHash, 0)) * pRnd.HashRandom(xHash, yHash, 1) + pRnd.HashRandom(xHash, yHash, 2)))
			/ (lambda);
	}

	float VerticalEdgeInverseDisplaced(PseudoRandom& pRnd, float x, float y, float lambda, float xHash, float yHash)
	{
		return (std::max(
			(x + pRnd.HashRandom(xHash, yHash, 1)) * pRnd.HashRandom(xHash, yHash, 0) + (y + pRnd.HashRandom(xHash, yHash, 0)) * pRnd.HashRandom(xHash, yHash, 1) + pRnd.HashRandom(xHash, yHash, 2), 
			-(x + pRnd.HashRandom(xHash, yHash, 1)) * pRnd.HashRandom(xHash, yHash, 0) - (y + pRnd.HashRandom(xHash, yHash, 0)) * pRnd.HashRandom(xHash, yHash, 1) + pRnd.HashRandom(xHash, yHash, 2))) 
			/ (lambda);
	}

	float Parabolic(PseudoRandom& pRnd, float x, float y, float lambda, float xHash, float yHash)
	{
		return
			(x * x * pRnd.HashRandom(xHash, yHash, 0) + y * y * pRnd.HashRandom(xHash, yHash, 1) + 2 * x * y * pRnd.HashRandom(xHash, yHash, 2))
			/ (lambda * lambda);
	}

	float ParabolicInverse(PseudoRandom& pRnd, float x, float y, float lambda, float xHash, float yHash)
	{
		return
			(pRnd.HashRandom(xHash, yHash, 2) * lambda - x * x * pRnd.HashRandom(xHash, yHash, 0) - y * y * pRnd.HashRandom(xHash, yHash, 1) - 2 * x * y * pRnd.HashRandom(xHash, yHash, 2))
			/ (lambda * lambda);
	}

	float ParabolicComposed(PseudoRandom& pRnd, float x, float y, float lambda, float xHash, float yHash)
	{
		return
			std::max(
			pRnd.HashRandom(xHash, yHash, 2) * lambda - x * x * pRnd.HashRandom(xHash, yHash, 0) - y * y * pRnd.HashRandom(xHash, yHash, 1) - 2 * x * y * pRnd.HashRandom(xHash, yHash, 2),
			(x * x * pRnd.HashRandom(xHash, yHash, 0) + y * y * pRnd.HashRandom(xHash, yHash, 1) + 2 * x * y * pRnd.HashRandom(xHash, yHash, 2)))
			/ (lambda * lambda);
	}
	float ParabolicComposedII(PseudoRandom& pRnd, float x, float y, float lambda, float xHash, float yHash)
	{
		return
			std::min(
			0.f,
			(x * x * pRnd.HashRandom(xHash, yHash, 0) + y * y * pRnd.HashRandom(xHash, yHash, 1) + 2 * x * y * pRnd.HashRandom(xHash, yHash, 2)))
			/ (lambda * lambda);
	}

	float DiplacedParabole(PseudoRandom& pRnd, float x, float y, float lambda, float xHash, float yHash)
	{
		return
			((x + lambda * pRnd.HashRandom(xHash, yHash, 0)) * (x + lambda * pRnd.HashRandom(xHash, yHash, 1)) 
			+ (y + lambda * pRnd.HashRandom(xHash, yHash, 2)) * (y + lambda * pRnd.HashRandom(xHash, yHash, 3)))
			/ (lambda * lambda);
	}
}
