#include <iostream>
#include <chrono>
#include "mcaFilter.h"
#include "filtering.h"
#include "resampler.h"
#pragma warning(disable:4996) 
#define NUM_DATA 19



class Manager
{
private:
    McaFilter mc;
    Resampler<float>        rsample[12];
	std::chrono::time_point<std::chrono::steady_clock> time_zero; // Start Time
	std::string filename;

public:
    Manager():time_zero(std::chrono::steady_clock::now()), filename("src/xpsp7.log"){}
	Manager(std::string& filename_) :time_zero(std::chrono::steady_clock::now()), filename(filename_) {}

    void run()
    {
		
		// Open file 
		FILE* pFile;
		float data_[NUM_DATA];

		float tstamp_prev = 0.0;
		bool first{ true };

		pFile = fopen(filename.c_str(), "r");
		if (pFile == NULL)
			perror("Error opening file");
		else {
			while (!feof(pFile))
			{
				// Read a line of data
				int res = fscanf(pFile, "%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f", &data_[0], &data_[1],
					&data_[2], &data_[3], &data_[4], &data_[5], &data_[6], &data_[7],
					&data_[8], &data_[9], &data_[10], &data_[11], &data_[12], &data_[13], &data_[14], &data_[15], &data_[16], &data_[17], &data_[18], &data_[19]);

				float dt = (float)std::chrono::duration_cast<std::chrono::microseconds>(
					std::chrono::steady_clock::now() - time_zero)
					.count();

				// PUNNING
				float dataCueInput[12];
				memcpy(dataCueInput, data_, sizeof(dataCueInput));

				for (int i = 0; i < 12; ++i)
					rsample[i].addSample(dataCueInput[i], dt);

				//Re sample for smooth data
				float data[12]{};
				for (int i = 0; i < 12; ++i)
					data[i] = rsample[i].filt(dt);

				//Filtering
				double filteredData[12];
				mc.filtering(data);
				mc.getData(filteredData);
				mc.logFilteredData(filteredData, data[0]);
			}
			fclose(pFile);
		}
    }

};

int main()
{
	Manager mg;
	mg.run();
}

