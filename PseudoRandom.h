#ifndef _PseudoRandom_H_
#define _PseudoRandom_H_

namespace GenericLatticeNoiseAlgorithm
{
	class PseudoRandom
	{
	private:
		float hashField[251][251][5];
		int seed;
		int seedA, seedB, seedC, seedD;

		bool useCSRandom;
		float scale;

	public:
		PseudoRandom() : useCSRandom(true)
		{
		}

		void SetCSharpRandomClasUsed(bool o)
		{
			useCSRandom = o;
		}
		void SetScale(float s)
		{
			scale = s;
		}

		float HashF(int i, int j, int index) //Returns a float -1 to 1. Not necessarily well weighted all over range
		{
			return (float)Hash(i, j, index) / 2147483648u;
		}

		int Rotate(int x, int b)
		{
			return (x << b) ^ (x >> (32 - b));
		}

		void Initiate(int s);
		float HashRandom(float x, float y, int index);
		int Hash(int i, int j, int index); //Returns a integer 0 to 2147483647
	};
}

#endif /* _PseudoRandom_H_ */