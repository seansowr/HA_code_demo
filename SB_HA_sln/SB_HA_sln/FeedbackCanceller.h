#pragma once

class FeedbackCanceller
{
private:
	double InputPower;
	double ErrorPowerFg;
	double ErrorPowerBg;
	double NoisePower;
	double ReferencePower;
	double CrossCov;
	const double AvgPwrCoef;
	const double NoisePwrCoef;
	const double NoisePwrLeak;

	double* ReferenceDelayBuffer;
	int PreDelay;
	int WriteIdx;
	int ReadIdx;

	double* BackgroundIr;
	double* ForegroundIr;
	double* FirBuffer;
	
	int timer;
	static const int TIMEOUT = 100;

public:
	FeedbackCanceller();
	~FeedbackCanceller();
	double RunFbc(double input);
	void SetPreDelay(int preDelaySamples);
	void BufferReference(double reference);

	// debug only
	double* GetFgIr()
	{
		return ForegroundIr;
	}
};