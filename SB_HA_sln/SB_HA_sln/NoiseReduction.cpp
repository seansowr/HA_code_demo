#include <cmath>
#include <algorithm>
#include "NoiseReduction.h"

/*	
*	Noise Reduction Constructor - Initializes parameters and estimator arrays
*	Inputs:
*			fs - Sample rate of the HA system
*			fftSize - number of fft points for HA system
*/
NoiseReduction::NoiseReduction(int fs, int fftSize) : NoiseLeak(1.0001)
{
	SampleRate = fs;
	FftSize = fftSize;
	NyquistBin = FftSize / 2 + 1;
	BinEnergy = new double[NyquistBin];
	memset(BinEnergy, 0, sizeof(double) * NyquistBin);
	NoiseEnergy = new double[NyquistBin];
	memset(NoiseEnergy, 0, sizeof(double) * NyquistBin);
	NoiseReductionCoef = new double[NyquistBin];
	NoiseReductionCoef[0] = 0.2;
	for (int k = 1; k < 30; k++)
	{
		NoiseReductionCoef[k] = 0.1;
	}
	for (int k = 30; k < NyquistBin; k++)
	{
		NoiseReductionCoef[k] = 0.5;
	}
	SnrBias = new double[NyquistBin];
	CalculateBiases();
	timer = 0;
	BinEnergyCoef = 1.0 - 1.0 / ((double(SampleRate) / double(FftSize)) * 0.01);
	NoiseEnergyCoef = 1.0 - 1.0 / ((double(SampleRate) / double(FftSize)) * 10.0);
	NrGains = new double[NyquistBin];
	Snr = new double[NyquistBin];
	for (int k = 0; k < NyquistBin; k++)
	{
		NrGains[k] = 1.0;
		Snr[k] = 1.0;
	}
}

/*
*	Noise Reduction Destructor - Deallocates memory from heap
*/
NoiseReduction::~NoiseReduction()
{
	delete[] BinEnergy;
	delete[] NoiseEnergy;
	delete[] NoiseReductionCoef;
	delete[] SnrBias;
	delete[] NrGains;
}

/*
*	Calculate Biases - calculates the SNR bias to be applied in each DFT bin up to the Nyquist rate
*	according to the formula bias(binIndex) = binIndex^(-1/3)
*/
void NoiseReduction::CalculateBiases()
{
	for (int k = 0; k < NyquistBin; k++)
	{
		SnrBias[k] = pow((k + 1), -1.0 / 3.0);
	}
}

/*
*	RunNr - Runs the noise reduction algorithm
*	Inputs:
*			fftIn - complex array of DFT values for the current frame
*	Outputs:
*			gains - array of noise reduction gains to be applied in each DFT bin
*/
void NoiseReduction::RunNr(std::complex<double>* fftIn, double* gains)
{
	if (timer < TIMEOUT)
	{
		// use simple average to initialize estimators over 100 frames
		double energySum = 0.0;
		for (int k = 0; k < NyquistBin; k++)
		{
			BinEnergy[k] += real(fftIn[k] * conj(fftIn[k])) / TIMEOUT;
			NoiseEnergy[k] = BinEnergy[k]; // std::max(BinEnergy[k], 1.0e-12);			// Set minimum noise floor to -120dB
			Snr[k] = BinEnergy[k] / (NoiseEnergy[k] + 1.0e-8);
			gains[k] = NrGains[k];
			energySum += NoiseEnergy[k];
		}
		timer++;
		if (energySum == 0.0) timer = 0.0;
	}
	else
	{
		// Run algorithm in each bin
		for (int k = 0; k < NyquistBin; k++)
		{
			// Estimate the energy in each bin
			BinEnergy[k] = BinEnergyCoef * BinEnergy[k] + (1 - BinEnergyCoef) * real(fftIn[k] * conj(fftIn[k]));
			// Estimate the noise in each bin
			if (BinEnergy[k] > NoiseEnergy[k])
			{
				// slowly increase if bin energy is greater than previous noise estimate
				NoiseEnergy[k] *= NoiseLeak;
			}
			else
			{
				// Otherwise use leaky integrator
				NoiseEnergy[k] = NoiseEnergyCoef * NoiseEnergy[k] + (1 - NoiseEnergyCoef) * BinEnergy[k];
			}
			NoiseEnergy[k] = std::max(NoiseEnergy[k], 1.0e-12);
			// Estimate the SNR
			Snr[k] = BinEnergy[k] / (NoiseEnergy[k] + 1.0e-8);
			// Apply bias and calculate the modified SNR^3
			double biasSnrCubed = pow(SnrBias[k] * Snr[k], 3.0);
			// Crossfade NR gains from previous frame to current frame
			NrGains[k] = 0.35 * biasSnrCubed / (biasSnrCubed + NoiseReductionCoef[k]) + 0.65 * NrGains[k];
			gains[k] = NrGains[k];
		}
	}
}