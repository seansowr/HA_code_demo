#include <algorithm>
#include "HearingAugmentation.h"
#include "HearingAugmentationParams.h"

// Butterworth filter coefficients for a HPF with cutoff at 100Hz
const double HearingAugmentation::HpfCoefficientsA[] = { -1.981488509144574, 9.816582826171344e-01 };
const double HearingAugmentation::HpfCoefficientsB[] = { 9.907866979404267e-01, -1.981573395880853, 9.907866979404267e-01 };
// Frequencies for octave edges for calibrated EQ
const double HearingAugmentation::OctaveEdgeFreqs[] = { 176.777, 353.553, 707.107, 1414.214, 2828.427, 5656.854, 11313.708, 22627.417 };

/*
*	Constructor - sets up HA object
*		Inputs:
*				fs - sampling rate of the system (default to 48kHz)
*/
HearingAugmentation::HearingAugmentation(double fs)
{
	// Allocate some memory
	fftState = kiss_fftr_alloc(DEFAULT_FFT_LENGTH, 0, 0, 0);
	ifftState = kiss_fftr_alloc(DEFAULT_FFT_LENGTH, 1, 0, 0);
	fftIn = new float[DEFAULT_FFT_LENGTH];
	memset(fftIn, 0, sizeof(float) * DEFAULT_FFT_LENGTH);
	ifftOut = new float[DEFAULT_FFT_LENGTH];
	memset(ifftOut, 0, sizeof(float) * DEFAULT_FFT_LENGTH);
	inputBuffer = new double[DEFAULT_FFT_LENGTH];
	memset(inputBuffer, 0, sizeof(double) * DEFAULT_FFT_LENGTH);
	outputFrame = new double[DEFAULT_FRAME_LENGTH];
	memset(outputFrame, 0, sizeof(double) * DEFAULT_FRAME_LENGTH);
	complexDft = new kiss_fft_cpx[DEFAULT_FFT_LENGTH];
	memset(complexDft, 0, sizeof(kiss_fft_cpx) * DEFAULT_FFT_LENGTH);
	dft = new std::complex<double>[DEFAULT_FFT_LENGTH];
	memset(dft, 0, sizeof(std::complex<double>) * DEFAULT_FFT_LENGTH);
	EqualizationGains = new double[DEFAULT_FFT_LENGTH];
	MinPhaseIr = new double[DEFAULT_FRAME_LENGTH];
	// Set EQ gains to unity
	for (int k = 0; k < DEFAULT_FFT_LENGTH; k++)
	{
		EqualizationGains[k] = 1.0;
	}
	for (int n = 0; n < DEFAULT_FRAME_LENGTH; n++)
	{
		MinPhaseIr[n] = 0.0;
		FirBuffer[n] = 0.0;
	}
	MinPhaseIr[0] = 1.0;	// Initialize min phase FIR to passthrough
	// All algorithms on by default
	EnableFbc = true;
	EnableNr = true;
	EnableDrc = true;

	SampleRate = fs;
}

/*
*	Destructor - cleans up memory
*/
HearingAugmentation::~HearingAugmentation()
{
	// clean up allocated memory
	free(fftState);
	free(ifftState);
	delete[] fftIn;
	delete[] inputBuffer;
	delete[] complexDft;
	delete[] dft;
}

/*
*	BufferInput - shifts old data to front of buffer, and puts new data at the back
*				e.g. [x(0), x(1), ... , x(n), x(n+1), ... x(n + N - 1)] where N is
*				the frame length and n indexes the first new sample
*		Inputs:
*				inputFrame - array holding new input data (must be of length 128)
*/
void HearingAugmentation::BufferInput(double* inputFrame)
{
	for (int n = 0; n < DEFAULT_FRAME_LENGTH; n++)
	{
		// Shift old data to front of buffer
		inputBuffer[n] = inputBuffer[n + DEFAULT_FRAME_LENGTH];

		// Put new data in second half of buffer
		inputBuffer[n + DEFAULT_FRAME_LENGTH] = inputFrame[n];
	}
}

/*
*	Get output - grants access to the output frame
*		Return:
*				output frame after HA processing
*/
double* HearingAugmentation::GetOutput()
{
	return outputFrame;
}

/*
*	RunFft - Runs the FFT on the input buffer
*/
void HearingAugmentation::RunFft()
{
	// convert to float for kissfft
	for (int k = 0; k < DEFAULT_FFT_LENGTH; k++)
	{
		fftIn[k] = (float)inputBuffer[k];
	}

	// run FFT
	kiss_fftr(fftState, fftIn, complexDft);

	// put data into useable array
	for (int k = 0; k < DEFAULT_FFT_LENGTH; k++)
	{
		dft[k] = { (double)complexDft[k].r, (double)complexDft[k].i };
	}
}

/*
*	MinPhaseDesign - Designs the minimum phase FIR that applies DRC and NR
*		Inputs:
*				magResponseSqrd - the magnitude-squared response of the HA
*								filter in the frequency domain
*/
void HearingAugmentation::MinPhaseDesign(double* magResponseSqrd)
{
	// Take the log of the magnitude-squared response
	kiss_fft_cpx cpxDft[DEFAULT_FFT_LENGTH];
	for (int k = 0; k < DEFAULT_FFT_LENGTH; k++)
	{
		cpxDft[k].r = logf((float)magResponseSqrd[k] + 1.0e-10f);
		cpxDft[k].i = 0.0f;
	}
	// Calculate the cepstrum
	float cepstrum[DEFAULT_FFT_LENGTH];
	kiss_fftri(ifftState, cpxDft, cepstrum);
	// Fold the cepstrum onto itself
	float holomorph[DEFAULT_FFT_LENGTH];
	for (int n = 0; n < DEFAULT_FFT_LENGTH/2; n++)
	{
		holomorph[n] = cepstrum[n] * 2.0f * IFFT_SCALE_FACTOR;
		holomorph[n + DEFAULT_FFT_LENGTH / 2] = 0.0f;
	}
	holomorph[0] /= 2.0f;
	// Take the FFT of the folded cepstrum
	kiss_fftr(fftState, holomorph, cpxDft);
	// Take the exponential of the holomorphic filter
	for (int k = 0; k < DEFAULT_FFT_LENGTH; k++)
	{
		cpxDft[k].r = expf(cpxDft[k].r / 2.0f);
		cpxDft[k].i = 0.0f;
	}
	// Take the IFFT to generate the minimum phase FIR
	kiss_fftri(ifftState, cpxDft, holomorph);
	for (int n = 0; n < DEFAULT_FRAME_LENGTH; n++)
	{
		MinPhaseIr[n] = (double)holomorph[n] * IFFT_SCALE_FACTOR;
	}
}

/*
*	RunFir - Runs the HA FIR filter
*		Inputs:
*				input - audio sample going into the HA filter
*		Return:
*				output sample after filtering
*/
double HearingAugmentation::RunFir(double input)
{
	// Put the newest sample at the back of the buffer
	FirBuffer[0] = input;
	// Run the convolution and shift the buffer forward
	double output = FirBuffer[0]*MinPhaseIr[0];
	for (int n = DEFAULT_FRAME_LENGTH-1; n > 0; n--)
	{
		output += FirBuffer[n] * MinPhaseIr[n];
		FirBuffer[n] = FirBuffer[n - 1];
	}
	return output;
}

/*
*	RunAntiRumbleFilter - Runs the anti-rumble IIR filter
*		Inputs:
*				input - audio sample going into the anti-rumble filter
*		Return:
*				output sample after filtering
*/
double HearingAugmentation::RunAntiRumbleFilter(double input)
{
	// Run direct-form II IIR
	double state = input - HpfCoefficientsA[0] * DirectForm2States[0] - HpfCoefficientsA[1] * DirectForm2States[1];
	double output = HpfCoefficientsB[0] * state + HpfCoefficientsB[1] * DirectForm2States[0] + HpfCoefficientsB[2] * DirectForm2States[1];
	// Update the saved states
	DirectForm2States[1] = DirectForm2States[0];
	DirectForm2States[0] = state;
	return output;
}

/*
*	RunHA - runs all the algorithms in the HA solution for the current frame
*/
void HearingAugmentation::RunHA()
{
	// Run FBC and get output
	double fbcOutput[DEFAULT_FRAME_LENGTH];
	if (EnableFbc)
	{
		for (int n = 0; n < DEFAULT_FRAME_LENGTH; n++)
		{
			FbcObj.BufferReference(outputFrame[n]);
			fbcOutput[n] = FbcObj.RunFbc(inputBuffer[DEFAULT_FRAME_LENGTH + n]);
		}
	}
	// If FBC is enabled copy the output into the input buffer for further processing
	if (EnableFbc)
	{
		memcpy(inputBuffer + DEFAULT_FRAME_LENGTH, fbcOutput, sizeof(double) * DEFAULT_FRAME_LENGTH);
	}
	// Arrays to hold the gains from NR and DRC
	double nrgains[DEFAULT_FFT_LENGTH / 2 + 1];
	double drcgains[DEFAULT_FFT_LENGTH / 2 + 1];
	// FFT of input buffer
	RunFft();
	// Run NR algorithm
	NrObj.RunNr(dft, nrgains);
	// Run DRC algorithm
	DrcObj.RunDrc(dft, drcgains);
	// Array to hold the magnitude-squared of the HA filter
	double magSqr[DEFAULT_FFT_LENGTH];
	double mag;
	for (int k = 0; k < DEFAULT_FFT_LENGTH / 2 + 1; k++)
	{
		// Apply EQ
		mag = EqualizationGains[k];
		// If enabled, apply NR
		if (EnableNr) mag *= nrgains[k];
		// If enabled, apply DRC
		if (EnableDrc) mag *= drcgains[k];
		// Calculate the magnitude-squared
		magSqr[k] = mag * mag;
	}
	// Copy the first half of the magnitude-squared into the second to ensure symmetry
	int idx = DEFAULT_FFT_LENGTH / 2 - 1;
	for (int k = DEFAULT_FFT_LENGTH / 2 + 1; k < DEFAULT_FFT_LENGTH; k++)
	{
		magSqr[k] = magSqr[idx--];
	}
	// Design the filter
	MinPhaseDesign(magSqr);
	for (int n = 0; n < DEFAULT_FRAME_LENGTH; n++)
	{
		// Apply the HA filter
		outputFrame[n] = RunFir(inputBuffer[DEFAULT_FRAME_LENGTH + n]);
		// Apply the anti-rumble filter
		outputFrame[n] = RunAntiRumbleFilter(outputFrame[n]);
		
	}
}

/*
*	SetParameter - sets a parameter enumerated in HearingAugmentationParams.h
*		Inputs:
*				paramId - enumerated ID for the parameter to set
*				paramValue - value to set for the chosen parameter
*		Return:
*				0 on success, -1 on failure
*/
int HearingAugmentation::SetParameter(int ParamId, double ParamValue)
{
	switch (ParamId)
	{
	case DrcOverallTargetLevel:
		return DrcObj.SetParameter(ParamId, ParamValue);
	case DrcBassTargetLevel:
		return DrcObj.SetParameter(ParamId, ParamValue);
	case DrcMidsTargetLevel:
		return DrcObj.SetParameter(ParamId, ParamValue);
	case DrcTrebleTargetLevel:
		return DrcObj.SetParameter(ParamId, ParamValue);
	case DrcMaxBoost:
		return DrcObj.SetParameter(ParamId, ParamValue);
	case DrcBassRatio:
		return DrcObj.SetParameter(ParamId, ParamValue);
	case DrcMidsRatio:
		return DrcObj.SetParameter(ParamId, ParamValue);
	case DrcTrebleRatio:
		return DrcObj.SetParameter(ParamId, ParamValue);
	case FbcPreDelay:
		if (ParamValue >= 0.0 && ParamValue <= 0.3)
		{
			int samples = int(SampleRate * ParamValue);
			FbcObj.SetPreDelay(samples);
			return 0;
		}
		else
		{
			return -1;
		}
	default:
		return -1;	// unsupported parameter
	}
}

/*
*	SetParameter - sets a parameter enumerated in HearingAugmentationParams.h
*		Inputs:
*				paramId - enumerated ID for the parameter to set
*				paramValue - value to set for the chosen parameter
*		Return:
*				0 on success, -1 on failure
*/
int HearingAugmentation::SetParameter(int ParamId, bool ParamValue)
{
	switch (ParamId)
	{
	case DrcEnable:
		EnableDrc = ParamValue;
		return 0;
	case NrEnable:
		EnableNr = ParamValue;
		return 0;
	case FbcEnable:
		EnableFbc = ParamValue;
		return 0;
	default:
		return -1;	// unsupported parameter
	}
}

/*
*	SetEqualizationGains - sets the gains for the 8 octaves as determined by calibration
*		Inputs:
*				octaveGains - array containing the linear (i.e. NOT in dB) gains for each octave
*/
void HearingAugmentation::SetEqualizationGains(double* octavegains)
{
	int octIdx = 0;
	double binCenterFreq;
	// Map octave gains to DFT bins
	for (int k = 0; k < DEFAULT_FFT_LENGTH / 2 + 1; k++)
	{
		binCenterFreq = k * SampleRate / DEFAULT_FFT_LENGTH;
		if (binCenterFreq > OctaveEdgeFreqs[octIdx])
		{
			octIdx++;
			octIdx = std::min(octIdx, NUM_EQ_OCTAVES - 1);
		}
		EqualizationGains[k] = octavegains[octIdx];
	}
}