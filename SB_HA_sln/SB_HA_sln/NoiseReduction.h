#pragma once

#include <complex>

class NoiseReduction
{
private:
	// System parameters
	int SampleRate;
	int FftSize;
	int NyquistBin;

	// Estimators
	double* BinEnergy;
	double* NoiseEnergy;
	double* NoiseReductionCoef;
	double* SnrBias;

	// Timer to init noise floor estimates
	int timer;
	static const int TIMEOUT = 10;

	// Estimator coefficients 
	double BinEnergyCoef;
	double NoiseEnergyCoef;
	const double NoiseLeak;

	// Gain array
	double* NrGains;
	// Snr estimator
	double* Snr;

	// Helper function
	void CalculateBiases();
public:
	NoiseReduction(int fs = 48000, int fftSize = 256);
	~NoiseReduction();
	void RunNr(std::complex<double>* fftIn, double* gains);
};