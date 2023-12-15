#pragma once

#include <complex>
#include <vector>

class DynamicRangeCompressor
{
private:
	// System parameters
	int SampleRate;
	int NumBands;
	int FftSize;
	// Crossover parameters
	std::vector<double> OctaveFreqs;
	int NumOctaves;
	int NumOctavesPerBand;
	double* CrossOverFreqs;
	int NyquistBin;
	double BinBandwidth;
	int* CrossOverBins;
	// Compression curve parameters
	double* Targets;
	double* Slopes;
	double MaxBoost;
	double* Thresholds;
	// Averaging coefficients
	double AttackTime;
	double AttackCoef;
	double ReleaseTime;
	double ReleaseCoef;
	const double NoiseCoef;
	const double NoiseLeak;
	// Estimators
	double* BandEnergy;
	double* AvgEnergy;
	double* NoiseEnergy;
	double* Gain;
	// Helper functions
	void ComputeOctaves(double startFrequency);
	void ComputeThresholds();
public:
	DynamicRangeCompressor(int numBands = 3, double attackTime = 0.01, double releaseTime = 0.3, double maxBoost = 30.0, int fftSize = 256, int fs = 48000);
	~DynamicRangeCompressor();
	int SetParameter(int paramId, double paramValue);
	void RunDrc(std::complex<double>* fftIn, double* gains);
};