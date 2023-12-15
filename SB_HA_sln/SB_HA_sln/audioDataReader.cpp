#include <iostream>
#include <fstream>
#include <string>
#include "audioDataReader.h"

audioDataReader::audioDataReader()
{

}

audioDataReader::~audioDataReader()
{

}

void audioDataReader::parseData(const char* filename)
{
	std::ifstream datafile;
	datafile.open(filename);
	std::string data;
	double value;

	if (datafile)
	{
		while (datafile >> data)
		{
			value = std::stod(data);
			audioData.push_back(value);
		}
	}

	datafile.close();

}

void audioDataReader::getData(std::vector<double>* data)
{
	data->assign(audioData.begin(), audioData.end());
}