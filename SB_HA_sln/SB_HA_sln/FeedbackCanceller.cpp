#include <cstring>
#include <algorithm>
#include "FeedbackCanceller.h"

#define MAX_FBC_PRE_DELAY 14400
#define FB_PATH_IR_LEN 128
#define NLMS_STEP_SIZE_MIN 0.01

FeedbackCanceller::FeedbackCanceller() : AvgPwrCoef(0.995), NoisePwrCoef(0.999997), NoisePwrLeak(1.000002)
{
	InputPower = 0.0;
	ReferencePower = 0.0;
	CrossCov = 0.0;
	ErrorPowerBg = 0.0;
	ErrorPowerFg = 0.0;
	NoisePower = 1.0e-8;

	ReferenceDelayBuffer = new double[MAX_FBC_PRE_DELAY];
	PreDelay = 1;
	WriteIdx = PreDelay;
	ReadIdx = 0;

	BackgroundIr = new double[FB_PATH_IR_LEN];
	memset(BackgroundIr, 0, sizeof(double) * FB_PATH_IR_LEN);
	ForegroundIr = new double[FB_PATH_IR_LEN];
	memset(ForegroundIr, 0, sizeof(double) * FB_PATH_IR_LEN);
	FirBuffer = new double[FB_PATH_IR_LEN];
	memset(FirBuffer, 0, sizeof(double) * FB_PATH_IR_LEN);

	timer = 0;
}

FeedbackCanceller::~FeedbackCanceller()
{
	delete[] ReferenceDelayBuffer;
	delete[] BackgroundIr;
	delete[] ForegroundIr;
	delete[] FirBuffer;
}

double FeedbackCanceller::RunFbc(double input)
{
	
	double DelayedRef = ReferenceDelayBuffer[0];
	//ReadIdx %= PreDelay;

	InputPower = AvgPwrCoef * InputPower + (1.0 - AvgPwrCoef) * input * input;
	ReferencePower = AvgPwrCoef * ReferencePower + (1.0 - AvgPwrCoef) * DelayedRef * DelayedRef;
	CrossCov = AvgPwrCoef * CrossCov + (1.0 - AvgPwrCoef) * DelayedRef * input;
	double crosscorr = CrossCov * CrossCov / (InputPower * ReferencePower + 1.0e-10);
	if (timer < TIMEOUT)
	{
		NoisePower = std::max(InputPower, 1.0e-8);
		timer++;
	}
	else if (InputPower > 2.0*NoisePower)
	{
		NoisePower *= NoisePwrLeak;
	}
	else
	{
		NoisePower = NoisePwrCoef * NoisePower + (1.0 - NoisePwrCoef) * InputPower;
	}
	
	for (int n = 0; n < FB_PATH_IR_LEN - 1; n++)
	{
		FirBuffer[n] = FirBuffer[n + 1];
	}
	FirBuffer[FB_PATH_IR_LEN - 1] = DelayedRef;
	double BgFirOut = 0.0;
	double FgFirOut = 0.0;
	for (int n = 0; n < FB_PATH_IR_LEN; n++)
	{
		BgFirOut += BackgroundIr[n] * FirBuffer[n];
		FgFirOut += ForegroundIr[n] * FirBuffer[n];
	}

	double BgError = input - BgFirOut;
	double FgError = input - FgFirOut;
	ErrorPowerBg = AvgPwrCoef * ErrorPowerBg + (1.0 - AvgPwrCoef) * BgError * BgError;
	ErrorPowerFg = AvgPwrCoef * ErrorPowerFg + (1.0 - AvgPwrCoef) * FgError * FgError;
	double fbrleBg = InputPower / (ErrorPowerBg + 1.0e-10);
	double fbrleFg = InputPower / (ErrorPowerFg + 1.0e-10);

	if (ErrorPowerBg < ErrorPowerFg && fbrleBg > 1.0 && InputPower < 1.14 * NoisePower && fbrleBg > fbrleFg)
	{
		memcpy(ForegroundIr, BackgroundIr, sizeof(double) * FB_PATH_IR_LEN);
		ErrorPowerFg = ErrorPowerBg;
	}

	double stepSize = NLMS_STEP_SIZE_MIN * fabs(1.0 - sqrt(InputPower) / (sqrt(ErrorPowerBg) + 1.0e-10));
	stepSize = std::min(stepSize, NLMS_STEP_SIZE_MIN);
	
	if (InputPower < 3.16*NoisePower && crosscorr < 0.2 && timer >= TIMEOUT)// && InputPower > NoisePower)
	{
		double referenceNorm = 0.0;
		for (int n = 0; n < FB_PATH_IR_LEN; n++)
		{
			referenceNorm += FirBuffer[n] * FirBuffer[n];
		}

		for (int n = 0; n < FB_PATH_IR_LEN; n++)
		{
			BackgroundIr[n] += stepSize * BgError * FirBuffer[n] / (referenceNorm + 1.0e-10);
		}
	}

	/*for (int n = FB_PATH_IR_LEN - 1; n > 0; n--)
	{
		FirBuffer[n] = FirBuffer[n - 1];
	}*/

	return FgError;
}

void FeedbackCanceller::SetPreDelay(int preDelaySamples)
{
	memset(ReferenceDelayBuffer, 0, sizeof(double) * MAX_FBC_PRE_DELAY);
	PreDelay = preDelaySamples + 1;

}

void FeedbackCanceller::BufferReference(double reference)
{
	for (int n = 0; n < PreDelay - 1; n++)
	{
		ReferenceDelayBuffer[n] = ReferenceDelayBuffer[n + 1];
	}
	ReferenceDelayBuffer[PreDelay - 1] = reference;
}