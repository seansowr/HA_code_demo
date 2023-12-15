#pragma once

#include <vector>

class audioDataReader
{
private:
	std::vector<double> audioData;
public:
	audioDataReader();
	~audioDataReader();
	void parseData(const char* filename);
	void getData(std::vector<double>* data);
};