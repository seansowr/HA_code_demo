#pragma once

#include "DynamicRangeCompressor.h"
#include "FeedbackCanceller.h"
#include "NoiseReduction.h"
extern "C" {
#include "kiss_fftr.h"
}

#define DEFAULT_SAMPLE_RATE 48000
#define DEFAULT_FFT_LENGTH 256
#define DEFAULT_FRAME_LENGTH 128
#define IFFT_SCALE_FACTOR 1/DEFAULT_FFT_LENGTH

class HearingAugmentation
{
private:
	// KISSFFT states
	kiss_fftr_cfg fftState;
	kiss_fftr_cfg ifftState;
	// some useful arrays
	double* inputBuffer;
	double* outputFrame;
	float* fftIn;
	float* ifftOut;
	kiss_fft_cpx* complexDft;
	std::complex<double>* dft;
	// HA filter arrays
	double* EqualizationGains;
	double* MinPhaseIr;
	double FirBuffer[DEFAULT_FRAME_LENGTH];
	// 3 main algorithms
	NoiseReduction NrObj;
	DynamicRangeCompressor DrcObj;
	FeedbackCanceller FbcObj;
	// Anti-rumble filter
	double DirectForm2States[2] = { 0 };
	static const double HpfCoefficientsA[];
	static const double HpfCoefficientsB[];
	// EQ constants
	static const double OctaveEdgeFreqs[];
	static const int NUM_EQ_OCTAVES = 8;
	// Enable/Disable algorithms
	bool EnableFbc;
	bool EnableNr;
	bool EnableDrc;
	// System sample rate
	double SampleRate;

	/*
	*	RunFft - Runs the FFT on the input buffer
	*/
	void RunFft();

	/*
	*	MinPhaseDesign - Designs the minimum phase FIR that applies DRC and NR
	*		Inputs:
	*				magResponseSqrd - the magnitude-squared response of the HA
	*								filter in the frequency domain
	*/
	void MinPhaseDesign(double* magResponseSqrd);

	/*
	*	RunFir - Runs the HA FIR filter
	*		Inputs:
	*				input - audio sample going into the HA filter
	*		Return:
	*				output sample after filtering
	*/
	double RunFir(double input);

	/*
	*	RunAntiRumbleFilter - Runs the anti-rumble IIR filter
	*		Inputs:
	*				input - audio sample going into the anti-rumble filter
	*		Return:
	*				output sample after filtering
	*/
	double RunAntiRumbleFilter(double input);
public:
	/*
	*	Constructor - sets up HA object
	*		Inputs:
	*				fs - sampling rate of the system (default to 48kHz)
	*/
	HearingAugmentation(double fs = DEFAULT_SAMPLE_RATE);

	/*
	*	Destructor - cleans up memory
	*/
	~HearingAugmentation();

	/*
	*	BufferInput - shifts old data to front of buffer, and puts new data at the back
	*				e.g. [x(0), x(1), ... , x(n), x(n+1), ... x(n + N - 1)] where N is
	*				the frame length and n indexes the first new sample
	*		Inputs:
	*				inputFrame - array holding new input data (must be of length 128)
	*/
	void BufferInput(double* inputFrame);

	/*
	*	RunHA - runs all the algorithms in the HA solution for the current frame
	*/
	void RunHA();

	/*
	*	Get output - grants access to the output frame
	*		Return:
	*				output frame after HA processing
	*/
	double* GetOutput();

	/*
	*	SetParameter - sets a parameter enumerated in HearingAugmentationParams.h
	*		Inputs:
	*				paramId - enumerated ID for the parameter to set
	*				paramValue - value to set for the chosen parameter
	*		Return:
	*				0 on success, -1 on failure
	*/
	int SetParameter(int paramId, double paramValue);
	int SetParameter(int paramId, bool paramValue);

	/*
	*	SetEqualizationGains - sets the gains for the 8 octaves as determined by calibration
	*		Inputs:
	*				octaveGains - array containing the linear (i.e. NOT in dB) gains for each octave
	*/
	void SetEqualizationGains(double* octaveGains);

	//debug only
	double* getfbcir()
	{
		return FbcObj.GetFgIr();
	}
};