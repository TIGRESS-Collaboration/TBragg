/**
 * g++ -g -O0 -o TBraggSimulation TBraggSimulation.cpp      ------------    compiler line
 * To run program, use following command:
 *		cat fileName (multiple files) |./TBraggSimulation 
 * or:
 *		./TBraggSimulation fileName (multiple files)
 * Can handle SRIM collision data from multiple isotopes
 * Ability to change (multiple) parameters from the command line: 
 *		./TBraggSimulation --parameter value fileName
**/

/**
 * Program Name: TBragg Simulation, v7
 * File: TBraggSimulation.cpp
 * Last Modified By: Owen Paetkau
 * Date: 01/12/2016
 * Purpose: This program is meant to take ion collision data (as simulated in SRIM/TRIM) and apply  
 * 	short and long energy filters to approximate energy and shape of the Bragg peak. The idea is to 
 * 	be able to simulate particle identification using the TBragg detector with multiple isotopes of 
 * 	similar atomic mass number (example: Sr-94, Rb-94, and Mo-94). The advantage to running the 
 * 	simulation is to determine if the expected constituents of a cocktail beam will be 
 *	distinguishable, and thus will help determine what to expect for beam tuning and the 
 * 	performance of diagnostics.
 * Changes from Previous Version (v6) (Owen Paetkau):
 *	- Commented out the impulseResponse due to error in averaging.
 * 	- Changed window2 in energyFilters to create proper distance between windows.
 *	- Also widened the range of the windows so it reaches before and after the simulated data set.
 *		--> This is more accurate for simulation and more closely matched measured values.
 *	- It is important to ensure the proper drift velocity. Any data taken from before October 5th, 2016
 *	  should be simulated with a drift velocity of 4.8cm/us while anything afterwards uses 6.6cm/us.
 * Changes from Previous Version (v5) (Alex Kurkijan):
 * 	- Adjusted maxDepth to represent proper depth for the TBragg detector (10.00 cm)
 *	- Implemented a new effDepth which is meant to simulate the different drift velocities between the 
 *	  Frisch grid and anode, in comparison to between the Mylar window and the Frisch grid
 *		--> Should note that while this may more accurately simulate the TBragg, it has little
 *			noticeable effect in the simulation results
 *  	- Added an impulse response function more effectively simulate the more gradual curve that results 
 * 	  in the PID charts from the TBragg detector.
 *		--> Implements a simple moving average over the energy values. The number of energy values 
 * 			sampled in the moving average can be changed. 
 *	- Arrays storing energy values are now explicitly overwritten before opening and reading collision
 *	  data from a new file
 *	- Added the capability to output simulation data (short and long filters values) to file for 
 *	  potential further analysis
 *	- Tweaked some formatting, added/changed comments as necessary
**/

#include <iostream>
#include <istream>
#include <algorithm>
#include <string>
#include <fstream>
#include <vector>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <cstring>
#include <math.h>
#include <cmath>
#include <utility>
#include <string>

using namespace std;
float	adcSampleRate = 0.02;				// Sample rate in microseconds
float 	l_peakt = 0.5;					// Long peak time in microseconds (isobutane - 4, CH4 - 0.5, P-8 - 0.5)
float 	l_gapt = 4.5;					// Long gap time in microseconds (isobutane - 2, CH4 - 4, P-8 - 4.5)
float 	s_peakt = 0.14;					// Short peak tme in microseconds (isobutane - 1, CH4 - 0.14, P-8 - 0.14)
float 	s_gapt = 0.24;					// Short gap time in microseconds (isobutane - 0.2, CH4 - 0.10, P-8 - 0.24)
//float 	driftVel = 6.6;					// Drift velocity in cm/usec for most of the tube (isobutane - 2.5, CH4 - 9, P-8 - 4.8)
float 	driftVel = 6.77;				// NOTE: Above value is prior to Oct 5th, 2016, below is afterwards.
//float 	driftVelFschAn = 6.6; 				// Drift velocity in cm/usec between Frisch grid and Anode (isobutane - 2.5, CH4 - 9, P-8 - 4.8)
float 	driftVelFschAn = 6.77;			// NOTE: Above value is prior to Oct 5th, 2016, below is afterwards.
float 	minDepth = 15000; 				// Beginning of active volume in angstroms
float 	maxDepth = 10; 					// Chamber length in centimetres
float 	fschDist = 0.01; 				// Distance from Frisch Grid to Anode in centimeters 
float 	effDepth = maxDepth-fschDist;		 	// Effective Depth of the ionization chamber
int   	conver = 100000000;	 			// Converting centimeters to angstroms
int   	arraySize = 10000; 				// Size of all arrays in this program 
							// (potential segmentation fault if not large enough)
int 	respRange = 35;					// Impulse response sample range
int   	maxSamples = (int) (effDepth/(driftVel*adcSampleRate) 
			+ fschDist/(driftVelFschAn*adcSampleRate)); // Maximum samples allowed due to the adc



float average (float energyValues[], int pos, float peakt) { 
	float  average = 0;
	int    count = 0;

	//printf("%d\t%d\n", pos, pos + (int)(peakt/adcSampleRate));
	// Loop through the array indexes as specified by pos and peakt
	for (int i = pos; i < (pos + (int)(peakt / adcSampleRate)); i++) { 
		if (i < 0) {
			average += energyValues[0];
			count++;
			//printf("\tsum : %.2f\taverage : %.2fkeV\t[%d]: %.2f\t count: %d\n", average, average / count, i, energyValues[i], count);
		} else if (i > (maxSamples - 1)) {
			average += energyValues[maxSamples - 1];
			count++;
			//printf("\tsum : %.2f\taverage : %.2fkeV\t[%d]: %.2f\t count: %d\n", average, average / count, i, energyValues[i], count);
		} else {
			average += energyValues[i];
			count++;
			//printf("\tsum : %.2f\taverage : %.2fkeV\t[%d]: %.2f\t count: %d\n",average, average / count, i, energyValues[i], count);
		}	
	}
	//printf("average = %.2f / %d = %.2f\n", average, count, average / count);	
	average = average / count;
	return average;
}

float sampledEnergies (float ionDistance[], float ionEnergy[], float closeEnergy[], int collisionNum) { 
	float 	simPos[arraySize]; 	// Simulated position according to sample rate and chamber length
	float	samTime[arraySize]; 	// Simulated time according to sample rate and chamber length
	float	anodeDistance;		// Distance from the anode
	int 	fschCount, loopCount;

	for (int i = 0; i < maxSamples; i++) {
		fschCount = 0;
		loopCount = 0;
		samTime[i] = i*adcSampleRate; // Simulates elapsed time from sample rate and iteration

		if (effDepth <= (driftVel * samTime[i]) ) { // If the ion is past the Frisch grid
			// Simulate the position from the time and drift velocity between Frisch grid and Anode
			simPos[i] = (maxDepth - effDepth - (driftVelFschAn*(samTime[i]) - samTime[fschCount]))*conver; 
			//printf("%.2f\t%.2f\n",samTime[i],simPos[i]);
		} else { // Ion is still in the main body of the detector
			// Simulate position based off time and drift velocity
			simPos[i] = (maxDepth - (driftVel*samTime[i]))*conver; 
			fschCount = i;
		}
	
		
		// Ensures that the energy drawn from text file is ONLY if an ion gets past the Mylar
		if (simPos[i] > minDepth) { 
			// Loop through energy values for one near-matching the simulated position			
			for (loopCount; loopCount < collisionNum; loopCount++) { 				
				anodeDistance = (maxDepth*conver) - ionDistance[loopCount];				

				if (simPos[i] > anodeDistance) { // Once found, end loop and write energy to array
					closeEnergy[i] = ionEnergy[loopCount];
					//printf("%.2f\t%.2f\n",anodeDistance,closeEnergy[i]);					
					loopCount = collisionNum;
				}		
			}
		} else { // Removes points below minDepth by ending the loop
			// Fill remaining points in the array before terminating the loop
			for (int j = i; j < maxSamples; j++) { 
				closeEnergy[j] = closeEnergy[i - 1];
			}			

			i = maxSamples;
		}
		//printf("%.2f\t%.2f\t%f\n",samTime[i],simPos[i],closeEnergy[i]);	
	}
} 


/*float impulseResponse (float energyArray[], float responseArray[], int responseRange, int nsamples) { 
	float	sum = 0;	
	
		// Simple moving average			
		for (int k = 0; k < nsamples; k++) { // Loop for each possible index of the impulse response
			// The following performs a similar function to the above average() function
			for (int l = k - responseRange; l < k; l++) { // loop through and summ the energy values within the range of the impulse response
				if (l < 0) {
					sum = sum;
					//sum += energyArray[0];
				} else if (l > nsamples - 1) {
					//sum += energyArray[nsamples - 1];
					sum = sum;
				} else {
					sum += energyArray[l];
					//printf("i = %d\t j = %d\t sum = %.2f\n", i, j, sum);				
				}
			}			
		
			responseArray[k] = (sum / responseRange);
			//printf("i = %d\t Actual Energy = %.2f\t Impulse average = %.2f\n", k, energyArray[k], responseArray[k]);
			sum = 0;
		} 	
	

} */

float energyFilters (float energyArray[], float peakt, float gapt, int nsamples) { // Long or short filter is determined from peak and gap times
	float	difference = 0;
	float 	maxValue = 0;
	float	window1;
	float 	window2;	

	//printf("%.2f\t%.2f\n",peakt,gapt);
	for (int i = - 1000; i < nsamples + 500; i++) { // Subtracted a value to give more complete samples on either side of the curve
		window1 = average(energyArray, i, peakt); //creates the first window
		window2 = average(energyArray, i + (int)((gapt + peakt)/adcSampleRate), peakt); //creates the second window

		difference = window1 - window2;

		//printf("%f\t%f\n",energyArray[i],window1);
		//printf("%d window1 - window2 = %.2f - %.2f = %.2f\n", i, window1, window2, difference);

		if (difference > maxValue) {
			maxValue = difference;
			//printf("%.2f\n", maxValue);
		}
	}	
	return maxValue;
}

float processSRIMData (istream& srimCollisionDataSource, float shortFilter[], float longFilter[], string ionName[], float ionMass[], float initialEnergy[], int ionNum[], int *totalIonsP, int *isotopeCountP) { // IonMass, ionNum, initialEnergy and ionName written to arrays to account for multiple isotopes (ie., files)

	float 	X[arraySize]; // Array of posi(int)(s_gapt / adcSampleRate)tions from text file
	float 	E[arraySize]; // Array of energies from text file
	float	closeEnergy[arraySize]; // Energies drawn from the text file corresponding to the decided distances
	float 	measuredEnergy[arraySize]; // Energies "measured" after convultions from the Frisch grid or preamp
	float	simPos[arraySize];

	int 	colNum = 0; // Number of collisions
	int 	index, ionCount;
	int		maxNum = 1;
	bool 	readIon = false;
	string 	line;

	while (getline(srimCollisionDataSource, line)) {
		if ( (line.length() > 30) && (line.substr (6,8).compare("Ion Name") == 0)) { // Determining element
			//(*ionNameP) = line.substr (23,2);
			ionName[*isotopeCountP] = line.substr (23,2); // Array to allow for multiple isotopes (ie., using cat command)
		}

		if ( (line.length() > 30) && (line.substr (6,8).compare("Ion Mass") == 0)) { // Determining ion mass
			//(*ionMassP) = strtof((line.substr (22,7)).c_str(),NULL);
			ionMass[*isotopeCountP] = strtof((line.substr (22,7)).c_str(),NULL);
		}

		if ( (line.length() > 30) && (line.substr (6,10).compare("Ion Energy") == 0)) { // Determining initial energy
			//(*initialEnergyP) = strtof((line.substr (18,11)).c_str(),NULL);
			initialEnergy[*isotopeCountP] = strtof((line.substr (18,11)).c_str(),NULL);
			*totalIonsP = *totalIonsP + ionNum[*isotopeCountP]; // Keeps track of the total number of ions for storing purposes	
			*isotopeCountP += 1; // Keeps track of total number of isotopes in the system
		}
		
		if ( line.substr (1,1).compare("0") == 0 || line.substr (1,1).compare("1") == 0 ) { // Identifying an ion line.
			readIon = true;
			ionCount = strtof((line.substr(1,5)).c_str(),NULL);

			X[colNum]  = atof ( (line.substr (17,10)).c_str() );
			E[colNum]  = atof ( (line.substr (7,9)).c_str() ) ;
			//printf("%.2f\t%.2f\n",X[colNum],E[colNum]);

			colNum += 1;
		} else {
			if ( readIon == true ) { // When end of ion is reached
				sampledEnergies(X, E, closeEnergy, colNum); // Used to simulate closeEnergy and store appropriate values within closeEnergy				
				
				//impulseResponse(closeEnergy, measuredEnergy, respRange, maxSamples); // Applies the TBragg's mystery impulse response
				index = ionCount + *totalIonsP;	// Stores each ion in a unique location in the short/long filters					

				shortFilter[index] = energyFilters(closeEnergy, s_peakt, s_gapt, maxSamples); // applies the short filter  (without the impulse response)
				longFilter[index] = energyFilters(closeEnergy, l_peakt, l_gapt, maxSamples); // applies the long filter  (without the impulse response)
				//shortFilter[index] = energyFilters(measuredEnergy, s_peakt, s_gapt, maxSamples); // applies the short filter (with the impulse response)
				//longFilter[index] = energyFilters(measuredEnergy, l_peakt, l_gapt, maxSamples); // applies the long filter (with the impulse response)
							

				ionNum[*isotopeCountP] += 1; // keeping track of number of ions in the current file
		
				//printf("%d\t%f\t%f\t%d\n",index, shortFilter[index], longFilter[index], ionNum[*isotopeCountP]);
								
				for (int k = 0; k < arraySize; k++) { // Wipes clean any array that stores energy values to zero
					// This is absolutely necessary for multiple isotopes!
					// Higher-energy ions (which have more collisions and more data to store) will 
					// not be overwritten by following lower-energy ions
					E[k] = 0;
					closeEnergy[k] = 0;
					measuredEnergy[k] = 0;
				}

				readIon = false;
				colNum = 0;
			} else {
				readIon = false;
			}
		}	
	}
}

float printValues (string ionName[], float ionMass[], float initialEnergy[], float shortFilter[], float longFilter[], int ionNum[], int isotopeCount) { 

	float 	averageShort;
	float	averageLong;
	int 	index, total = 0; // Used to keep track of the total number of ions over all the files

	ofstream fileOutData ("output_filter_data.txt"); // Raw number data from each run
	ofstream fileOutInfo ("output_ion_info.txt"); // Ion Summary output information to be matched with the raw number data

	fileOutData << "###### TBragg Simulation v6 -- Energy Filters Output ######" << endl;
	fileOutData << "# To match the below index number with its respective isotope, please see output_ion_info.txt \n\n" << endl;
	fileOutData << "Index \t Long Filter \t Short Filter" << endl;

	fileOutInfo << "###### TBragg Simulation v6 -- Isotope information Summary ######\n\n" << endl;

	for (int k = 0; k < isotopeCount; k++) { // Loops over values to allow for each initial conditions to be printed.
		averageShort = 0;
		averageLong = 0;

		if (fileOutInfo.is_open()) {
			fileOutInfo << "# Index #" << k + 1 << endl; // This number will match w/output_filter_data
			fileOutInfo << "# Ion Name : " << ionName[k] << endl;
			fileOutInfo << "# Ion Mass : " << ionMass[k] << "u" << endl;
			fileOutInfo << "# Initial Energy : " << initialEnergy[k] << " keV" << endl;
			fileOutInfo << "# Energy per nucleon : " << initialEnergy[k] / ionMass[k] << " keV/u" << endl;
		}

		cout << "# Ion Name : " << ionName[k] << endl;
		cout << "# Ion Mass : " << ionMass[k] << "u" << endl;
		cout << "# Initial Energy : " << initialEnergy[k] << " keV" << endl;
		cout << "# Energy per nucleon : " << initialEnergy[k] / ionMass[k] << " keV/u" << endl;
		
		//printf("#Chamber Length : %f cm\n",maxLength);
	
		for (int j = 1; j < ionNum[k + 1] + 1; j++) { // Averages and prints the short and long values
			index =	j + total; // The unique index where each value is stored
			averageShort += shortFilter[index];
			averageLong  += longFilter[index];
			printf("%.2f, %.2f\n", longFilter[index], shortFilter[index]);

			if (fileOutData.is_open()) {
				fileOutData << k + 1 << "\t\t\t" << longFilter[index] << "\t\t\t" << shortFilter[index] << endl;
			}
		}

		averageShort = averageShort / ionNum[k + 1]; // Is ionNum-1, as ionNum is originally defined as 1 instead of 0.
		averageLong = averageLong / ionNum[k + 1];
		total = total + ionNum[k + 1]; // To keep track of the total number of ions over all files

		printf("#There are %d ions in this file.\n", ionNum[k + 1]); // Assumes that there is the same number of ions in each file.
		printf("#AVG SHORT FILTER : %.2fkeV\n", averageShort);	
		printf("#AVG LONG FILTER : %.2fkeV\n\n", averageLong);

		if (fileOutData.is_open()) {
			fileOutData << endl; // Add an extra space between ions (easier to read later)
		}

		if (fileOutInfo.is_open()) {
			fileOutInfo << "# There are " << ionNum[k + 1] << " ions in this file." << endl;
			fileOutInfo << "# AVG SHORT FILTER: \t" << averageShort << " keV" << endl;
			fileOutInfo << "# AVG LONG FILTER: \t" << averageLong << " keV \n" << endl;
		}
	}

	fileOutData.close();
	fileOutInfo.close();
}

int setConstants (char *itemName[], int *iP, float *adcSampleRateP, float *l_peaktP, float *l_gaptP, float *s_peaktP, float *s_gaptP, float *driftVelP, float *driftVelFschAnP, float *minDepthP, float *maxDepthP, float *fschDistP, int *respRangeP) { 

	*iP += 1; // Increases i to skip next loop as well as store value that comes afterwards	
	string parameter = itemName[*iP - 1]; // converts the argv into a string again

	if (parameter.substr(2,13) == "adcSampleRate")	{
		*adcSampleRateP = atof(itemName[*iP]);
		//printf("%f\n",*adcSampleRateP);
	}

	if (parameter.substr(2,7) == "l_peakt") {
		*l_peaktP = atof(itemName[*iP]);
		//printf("%f\n",*l_peaktP);
	}

	if (parameter.substr(2,6) == "l_gapt")	{
		*l_gaptP = atof(itemName[*iP]);
		//printf("%f\n",*l_gaptP);
	}

	if (parameter.substr(2,7) == "s_peakt") {
		*s_peaktP = atof(itemName[*iP]);
		//printf("%f\n",*s_peaktP);
	}


	if (parameter.substr(2,6) == "s_gapt")	{
		*s_gaptP = atof(itemName[*iP]);
		//printf("%f\n",*s_gaptP);
	}

	if (parameter.substr(2,8) == "driftVel") {
		*driftVelP = atof(itemName[*iP]);
		//printf("%f\n",*driftVelP);
	}

	if (parameter.substr(2,14) == "driftVelFschAn") {
		*driftVelFschAnP = atof(itemName[*iP]);
		//printf("%f\n",*driftVelFschAnP);
	}

	if (parameter.substr(2,8) == "minDepth") {
		*minDepthP = atof(itemName[*iP]);
		//printf("%f\n",*minDepthP);
	}

	if (parameter.substr(2,8) == "maxDepth") {
		*maxDepthP = atof(itemName[*iP]);
		//printf("%f\n",*maxDepthP);
	}

	if (parameter.substr(2,8) == "fschDist") {
		*fschDistP = atof(itemName[*iP]);
		//printf("%f\n",*fschDistP);
	}

	if (parameter.substr(2,8) == "respRange") {
		*respRangeP = atof(itemName[*iP]);
		//printf("%f\n",*respRangeP);
	}
}

int main(int argc, char* argv[]) {
	float	shortFilter[arraySize]; // Returned short values from filters for each ion
	float 	longFilter[arraySize]; // Returned long values from filters for each ion

	int 	ionNum[arraySize]; // Number of ions (total value over all files in the case of cat command)
	int	isotopeCount = 0; // Used to count number of isotopes when using the cat command
	int	totalIons = 0; // Used to keep track of the total ions over all of the files. Needs to be set here so it doesn't reset when using "./TBraggSimulation" input method (see arg loop)
	string 	fileName;
	float 	initialEnergy[arraySize]; // Initial energy for all of the ions
	float 	ionMass[arraySize]; // Mass of the ion
	string *ionName = new string[arraySize];

	if (argc > 1) { // Reads through the number of files after the executable call of format "./TBraggSimulation Si28.txt Mo95.txt ..."
		for (int i = 1; i < argc; i++) { // Looping over each argument in the above executable call
			string argv_ = argv[i]; // Converts the argument to a string so it is easier to manipulate 
		
			if (strcmp(argv_.substr(0,2).c_str(), "--") == 0) { // Determining if there is a command to set a parameter
			 	setConstants(argv, &i, &adcSampleRate, &l_peakt, &l_gapt, &s_peakt, &s_gapt, &driftVel, &driftVelFschAn, &minDepth, &maxDepth, &fschDist, &respRange); 
				// Accesses subroutine to define parameter called upon. More parameters can be easily added to this subroutine.		
			} else { 
				fileName = argv_;	
				ifstream myfile; 	
				myfile.open(fileName.c_str());
		
				if(myfile.is_open()) { // Reads that the file is open and completes the subroutine
					processSRIMData(myfile, shortFilter, longFilter, ionName, ionMass, initialEnergy, ionNum, &totalIons, &isotopeCount);
				}
			}	
		}

		printValues(ionName, ionMass, initialEnergy, shortFilter, longFilter, ionNum, isotopeCount); // To allow for loop to perform properly with stored values in the array

	} else { // Processes and prints values from the line "cat filename (multiple) | ./TBraggSimulation"
		processSRIMData(cin, shortFilter, longFilter, ionName, ionMass, initialEnergy, ionNum, &totalIons, &isotopeCount);
		printValues(ionName, ionMass, initialEnergy, shortFilter, longFilter, ionNum, isotopeCount);
	}
}


