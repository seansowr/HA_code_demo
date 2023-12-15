#include <cmath>
#include "DynamicRangeCompressor.h"
#include "DynamicRangeCompressorParams.h"

/*
*	Dynamic Range Compressor Constructor - Initializes estimators and sets/calculates compression parameters
*	Inputs:
*			numBands - number of compression bands
*			attackTime - time constant for the attack average coefficient
*			releaseTime - time constant for the release average coefficient
*			maxBoost - the maximum amount of gain (in dB) that can be applied to each band
*			fftSize - the number of points in the DFT of the HA system
*			fs - sampling rate of the HA system
*/
DynamicRangeCompressor::DynamicRangeCompressor(int numBands, double attackTime, double releaseTime, double maxBoost, int fftSize, int fs) : NoiseCoef(0.99947), NoiseLeak(1.0001)
{
	SampleRate = fs;
	NumBands = numBands;
	FftSize = fftSize;
	ComputeOctaves(20.0);
	NumOctavesPerBand = NumOctaves / NumBands;
	CrossOverFreqs = new double[NumBands - 1];
	for (int k = 0; k < NumBands - 1; k++)
	{
		CrossOverFreqs[k] = OctaveFreqs[(k + 1) * NumOctavesPerBand - 1];
	}
	
	NyquistBin = FftSize / 2 + 1;
	BinBandwidth = (double)SampleRate / (double)FftSize;
	CrossOverBins = new int[NumBands];
	double f = 0.0;
	int binIdx = 0;
	for (int k = 0; k < NumBands - 1; k++)
	{
		while (f < CrossOverFreqs[k])
		{
			f += BinBandwidth;
			binIdx++;
		}
		CrossOverBins[k] = binIdx;
	}

	Targets = new double[NumBands];
	memset(Targets, 0, sizeof(double) * NumBands);
	Slopes = new double[NumBands];
	memset(Slopes, 0, sizeof(double) * NumBands);
	if (NumBands == 3)
	{
		Targets[0] = -20.0;
		Targets[1] = -8.0;
		Targets[2] = -12.0;
		Slopes[0] = 1.0 / 8.0;
		Slopes[1] = 1.0 / 8.0;
		Slopes[2] = 1.0 / 8.0;
	}
	MaxBoost = maxBoost;
	Thresholds = new double[NumBands];
	ComputeThresholds();
	AttackTime = attackTime;
	ReleaseTime = releaseTime;
	AttackCoef = 1.0 - 1.0 / (BinBandwidth * AttackTime);
	ReleaseCoef = 1.0 - 1.0 / (BinBandwidth * ReleaseTime);

	BandEnergy = new double[NumBands];
	memset(BandEnergy, 0, sizeof(double) * NumBands);
	AvgEnergy = new double[NumBands];
	memset(AvgEnergy, 0, sizeof(double) * NumBands);
	NoiseEnergy = new double[NumBands];
	Gain = new double[NumBands];
	for (int k = 0; k < NumBands; k++)
	{
		NoiseEnergy[k] = 1.0e-8;
		Gain[k] = 1.0;
	}
}

/*
*	Dynamic Range Compressor Constructor - deallocates memory on the heap
*/
DynamicRangeCompressor::~DynamicRangeCompressor()
{
	delete[] CrossOverFreqs;
	delete[] CrossOverBins;
	delete[] Targets;
	delete[] Slopes;
	delete[] Thresholds;
	delete[] BandEnergy;
	delete[] AvgEnergy;
	delete[] NoiseEnergy;
	delete[] Gain;
}

/*
*	Compute Octaves - Calculates the frequencies of the octaves of some start frequency up to the Nyquist rate
*	Inputs:
*			startFrequency - frequency from which octave frequencies are computed
*/
void DynamicRangeCompressor::ComputeOctaves(double startFrequency)
{
	OctaveFreqs.push_back(startFrequency);
	while (*(OctaveFreqs.end() - 1) < SampleRate / 2)
	{
		OctaveFreqs.push_back(*(OctaveFreqs.end() - 1) * 2);
	}
	NumOctaves = OctaveFreqs.size();
}

/*
*	Compute Thresholds - Computes the compression thresholds based user set parameters
*/
void DynamicRangeCompressor::ComputeThresholds()
{
	for (int k = 0; k < NumBands; k++)
	{
		Thresholds[k] = (Targets[k] - (Slopes[k] * Targets[k]) - MaxBoost) / (1 - Slopes[k]);
	}
}

/*
*	Set Parameter - Sets user parameters
*	Inputs:
*			paramId - the enumerated parameter ID
*			paramValue - the value of the parameter to be set
*	Outputs:
*			returns 0 on success
*			returns -1 on failure
*/
int DynamicRangeCompressor::SetParameter(int paramId, double paramValue)
{
	switch (paramId)
	{
	case OverallTargetLevel:
		if (paramValue <= 0.0 && paramValue >= -40.0)
		{
			Targets[0] = paramValue;
			Targets[1] = paramValue;
			Targets[2] = paramValue;
			ComputeThresholds();
			return 0;
		}
		else return -1;
	case BassTargetLevel:
		if (paramValue <= 0.0 && paramValue >= -40.0)
		{
			Targets[0] = paramValue;
			ComputeThresholds();
			return 0;
		}
		else return -1;
	case MidsTargetLevel:
		if(paramValue <= 0.0 && paramValue >= -40.0)
		{
			Targets[1] = paramValue;
			ComputeThresholds();
			return 0;
		}
		else return -1;
	case TrebleTargetLevel:
		if (paramValue <= 0.0 && paramValue >= -40.0)
		{
			Targets[2] = paramValue;
			ComputeThresholds();
			return 0;
		}
		else return -1;
	case MaxBoostParam:
		if (paramValue <= 40.0 && paramValue >= 0.0)
		{
			MaxBoost = paramValue;
			ComputeThresholds();
			return 0;
		}
		else return -1;
	case BassRatio:
		if (paramValue <= 20.0 && paramValue >= 1.0)
		{
			Slopes[0] = 1 / paramValue;
			ComputeThresholds();
			return 0;
		}
		else return -1;
	case MidsRatio:
		if (paramValue <= 20.0 && paramValue >= 1.0)
		{
			Slopes[1] = 1 / paramValue;
			ComputeThresholds();
			return 0;
		}
		else return -1;
	case TrebleRatio:
		if (paramValue <= 20.0 && paramValue >= 1.0)
		{
			Slopes[2] = 1 / paramValue;
			ComputeThresholds();
			return 0;
		}
		else return -1;
	default: return -1;		// Unsupported parameter
	}
}

/*
*	RunDrc - Runs the multi-band dynamic range compression algorithm
*	Inputs:
*			fftIn - complex array of DFT values for the current frame
*	Outputs:
*			gains - the array of gains to be applied in each DFT bin
*/
void DynamicRangeCompressor::RunDrc(std::complex<double>* fftIn, double* gains)
{
	// Clear band energy estimator
	memset(BandEnergy, 0, sizeof(double) * NumBands);
	int m = 0;
	// Calculate the power in each band occording to bin crossovers
	for (int k = 0; k < NyquistBin; k++)
	{
		BandEnergy[m] += 2 * std::real(fftIn[k] * conj(fftIn[k])) / FftSize;
		if (k == 0 || k == NyquistBin - 1)
		{
			BandEnergy[m] -= std::real(fftIn[k] * conj(fftIn[k])) / FftSize;
		}
		if (k == CrossOverBins[m])
		{
			m++;
		}
	}
	
	for (m = 0; m < NumBands; m++)
	{
		// Estimate the average energy in each bin according to attack and release behavior
		// Attack when the new estimate is > than the average, release when the new estimate is < than the average
		if (BandEnergy[m] > AvgEnergy[m])
		{
			AvgEnergy[m] = AttackCoef * AvgEnergy[m] + (1.0 - AttackCoef) * BandEnergy[m];
		}
		else
		{
			AvgEnergy[m] = ReleaseCoef * AvgEnergy[m] + (1.0 - ReleaseCoef) * BandEnergy[m];
		}
		// Estimate the noise energy in each band
		if (AvgEnergy[m] > NoiseEnergy[m])
		{
			NoiseEnergy[m] *= NoiseLeak;
		}
		else
		{
			NoiseEnergy[m] = NoiseCoef * NoiseEnergy[m] + (1.0 - NoiseCoef) * AvgEnergy[m];
		}
		// If the average energy is 5dB > than the noise floor, apply gain according to compression parameters
		if (AvgEnergy[m] > 3.1623 * NoiseEnergy[m])
		{
			Gain[m] = std::pow(10.0, MaxBoost / 20.0);
			double energy_dB = 10.0 * std::log10(AvgEnergy[m]);
			if (energy_dB > Thresholds[m])
			{
				double gain_dB = MaxBoost + (energy_dB - Thresholds[m]) * (Slopes[m] - 1.0);
				Gain[m] = std::pow(10.0, gain_dB / 20.0);
			}
		}
	}
	// Map the gains in the compression bands to the corresponding DFT bins
	m = 0;
	for (int k = 0; k < NyquistBin; k++)
	{
		gains[k] = Gain[m];
		if (k == CrossOverBins[m]) m++;
	}
}