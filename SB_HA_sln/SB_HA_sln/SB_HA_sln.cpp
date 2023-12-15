// SB_HA_sln.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <fstream>
#include "audioDataReader.h"
#include "HearingAugmentation.h"
#include "HearingAugmentationParams.h"

#define FS 48000
#define FBK_DELAY_SEC 0.2
//#define FBK

double RunFeedbackFir(double input, double* inputBuffer, double* impulseResponse, int bufferLength);

int main()
{
    // Create some user parameters for the Hearing Augmentation processing
    bool enableDrc = true;
    bool enableNr = true;
    bool enableFbc = false;
    // Target levels are in units of dBFS (FS = full scale) where 0dBFS is the maximum possible level in a digital system
    double overallTargetLeveldBFS = -20.0;  // This will set bass, mids, and treble target levels to -20dBFS. 
                                            // This is the parameter to control the overall loudness the DRC will aim for

    double midsTargetLeveldBFS = -8.0;      // This will set the mids target level 12dB higher than the overall target
                                            // i.e. -20 + 12 = -8

    double trebleTargetLeveldBFS = -12.0;   // This will set the treble target level 8dB higher than the overall target
                                            // i.e. -20 + 8 = -12;

    double maxBoostdB = 30.0;               // This controls the maximum amount of positive gain (boost) that can be
                                            // applied to quiet signals. This should be set relatively high value for
                                            // any amplification (e.g. between 20dB and 40dB. 30dB is a good default)

    double ratio = 10.0;                    // This will set the compression ratio to 10:1. Any signals that exceed the
                                            // DRC threshold (calculated internally) will be attenuated by this ratio
                                            // (e.g. if the close talker's level exceeds the threshold by 10dB, it will be
                                            // attenuated by 9dB). A higher ratio will make the level sound more even,
                                            // while a lower ratio will make the audio sound more natural.



    // Read in some audio data from a txt file
    const char* filepath1 = "C:\\Users\\Owner\\Documents\\OctaveCode\\example1_input_BEATS.txt";
    std::vector<double> data1;
    audioDataReader dataReader;
    dataReader.parseData(filepath1);
    dataReader.getData(&data1);
    std::vector<double> output;
    
#ifdef FBK
    const int FbkDlySamples = int(FS * FBK_DELAY_SEC) + 128;
    double* fbkFir = new double[FbkDlySamples];
    memset(fbkFir, 0, sizeof(double) * FbkDlySamples);
    for (int n = int(FS * FBK_DELAY_SEC); n < FbkDlySamples; n++)
    {
        fbkFir[n] = double(rand() % 100 - 50) * 0.001;
    }
    fbkFir[int(FS * FBK_DELAY_SEC)] = 0.1;
    double* fbkBuffer = new double[FbkDlySamples];
    memset(fbkBuffer, 0, sizeof(double) * FbkDlySamples);
    double fbkPathInput[DEFAULT_FRAME_LENGTH] = { 0 };
    double fbkPathOutput = 0.0;
#endif

    // Create hearing augmentation object
    HearingAugmentation ha;

    // Enable/disable internal algorithms
    ha.SetParameter(DrcEnable, enableDrc);
    ha.SetParameter(NrEnable, enableNr);
    ha.SetParameter(FbcEnable, enableFbc);

    // Set some parameters (described above)
    ha.SetParameter(DrcOverallTargetLevel, overallTargetLeveldBFS);
    //ha.SetParameter(DrcBassTargetLevel, overallTargetLeveldBFS);      <-- This is redundant if I want the bass target
                                                                        //  to be the same as the overall target
    ha.SetParameter(DrcMidsTargetLevel, midsTargetLeveldBFS);
    ha.SetParameter(DrcTrebleTargetLevel, trebleTargetLeveldBFS);
    ha.SetParameter(DrcMaxBoost, maxBoostdB);
    ha.SetParameter(DrcBassRatio, ratio);
    ha.SetParameter(DrcMidsRatio, ratio);
    ha.SetParameter(DrcTrebleRatio, ratio);     // I've set all bands to the same ratio but they could all be different

#ifdef FBK
    double* fbcIrEst;
#endif

    for (int n = 0; n < data1.size() - DEFAULT_FRAME_LENGTH; n += DEFAULT_FRAME_LENGTH)
    {
#ifdef FBK
        for (int nn = 0; nn < DEFAULT_FRAME_LENGTH; nn++)
        {
            fbkPathOutput = RunFeedbackFir(fbkPathInput[nn], fbkBuffer, fbkFir, FbkDlySamples);
            data1[n + nn] += fbkPathOutput;
        }
#endif
        ha.BufferInput(&data1[n]);
        ha.RunHA();
        double* outframe = ha.GetOutput();
        
        
        for (int nn = 0; nn < DEFAULT_FRAME_LENGTH; nn++)
        {
            output.push_back(outframe[nn]);
#ifdef FBK
            fbkPathInput[nn] = outframe[nn];
            fbcIrEst = ha.getfbcir();
#endif
        }
        
    }
    std::ofstream outfile1;
    outfile1.open("C:\\Users\\Owner\\Documents\\HAoutput1.txt");
    int outSize = output.size();
    for (int n = 0; n < outSize; n++)
    {
        outfile1 << output[n] << std::endl;
    }
    outfile1.close();

    return 0;

}

double RunFeedbackFir(double input, double* inputBuffer, double* impulseResponse, int bufferLength)
{
    inputBuffer[0] = input;
    double output = inputBuffer[0] * impulseResponse[0];
    for (int n = bufferLength - 1; n > 0; n--)
    {
        output += inputBuffer[n] * impulseResponse[n];
        inputBuffer[n] = inputBuffer[n - 1];
    }
    return output;
}

// Run program: Ctrl + F5 or Debug > Start Without Debugging menu
// Debug program: F5 or Debug > Start Debugging menu

// Tips for Getting Started: 
//   1. Use the Solution Explorer window to add/manage files
//   2. Use the Team Explorer window to connect to source control
//   3. Use the Output window to see build output and other messages
//   4. Use the Error List window to view errors
//   5. Go to Project > Add New Item to create new code files, or Project > Add Existing Item to add existing code files to the project
//   6. In the future, to open this project again, go to File > Open > Project and select the .sln file
