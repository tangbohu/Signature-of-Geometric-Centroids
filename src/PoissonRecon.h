#pragma once
#if defined DLL_EXPORT
#define DECLDIR __declspec(dllexport)
#else
#define DECLDIR __declspec(dllimport)
#endif

class DECLDIR PoissonRecon
{
public:
	 int MPoisson(char* In,char* Out,int Depth=8,int Degree=2,bool Verbose=0,int Refine=3,int IsoDivide=8,
		int SolverDivide=8,float Scale=1.25,float SamplesPerNode=1.0,
		bool Binary=0,bool Confidence=0,bool NoResetSamples=0,bool NoClipTree=0);

};