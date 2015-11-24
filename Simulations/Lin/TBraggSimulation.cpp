//g++ -g -O0 -o TBraggSimulation TBraggSimulation_v5.cpp      ------------    compiler line
//To run program, use following command : cat fileName (multiple files) |./TBraggSimulation OR ./TBraggSimulation fileName (multiple files)
//It can do multiple isotopes at a time where said files are Collision files from SRIM simulations.
//Ability to change parameters from the command line. Use the syntax of "./TBraggSimulation --parameter value fileName". This can be done for multiple parameters


//------------------------------------------------------Owen Paetkau----------------------------------------------------------------------------------------------------------------------------
//----------------------TBragg Simulation---------------31/07/2015-------------------v5 Changes-------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//-----This program is able to take the Collision files found -----------------------05/08/2015 - Added a count to keep track of the total ions as well as an array that stores the number of---
//-----in the SRIM output and apply short and long energy filters--------------------ions in each file with a corresponding isotope (or file) number. This removes the need to assume the-------
//-----to determine the total energy and shape of the Bragg peak.--------------------number of ions in each file making the readout more complete.----------------------------------------------
//-----This allows you to do proper particle identification in-----------------------NOTE : Any changes to parameters must come before the file names so they affect the analysis.--------------
//-----combination with the SRIM software.------------------------------------------------------------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//-----This version has complete functionality and will allow-----------------------------------------------------------------------------------------------------------------------------------
//-----up to 10000 ions and 10000 collisions per ion. It will-----------------------------------------------------------------------------------------------------------------------------------
//-----output a list of short and long filters in keV for each----------------------------------------------------------------------------------------------------------------------------------
//-----ion as well as the initial conditions of the simulations.--------------------------------------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

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
float adcSampleRate = 0.02;//sample rate in microseconds
float l_peakt = 0.5; //long peak time in microseconds
float l_gapt = 4.5; //long gap time in microseconds
float s_peakt = 0.14; //short peak tme in microseconds
float s_gapt = 0.24; // short gap time in microseconds
float driftVel = 4.8;//drift velocity in cm/micron
float minDepth = 15000; //beginning of active volume in angstroms
float maxDepth = 12.28; //chamber length in centimetres
int   conver = 100000000; //converting centimeters to angstroms
int   arraySize = 10000; //size of all arrays in this program (potential location of segmentation fault if not large enough)
int   maxSamples = maxDepth / (driftVel * adcSampleRate); //maximum samples allowed due to the adc


float average (float energyValues[], int pos, float peakt) //determine the average value in some window
{

float  average = 0;
int    count = 0;

//printf("%d\t%d\n", pos, pos + (int)(peakt/adcSampleRate));
	for (int i = pos; i < (pos + (int)(peakt/ adcSampleRate)); i++)
	{
		if ( i < 0)
		{
			average = average + energyValues[0];
			count++;
			//printf("\tsum : %.2f\taverage : %.2fkeV\t[%d]: %.2f\t count: %d\n", average, average / count, i, energyValues[i], count);
		}
		
		else if (i > (maxSamples - 1))
		{
			average = average + energyValues[maxSamples - 1];
			count++;
			//printf("\tsum : %.2f\taverage : %.2fkeV\t[%d]: %.2f\t count: %d\n", average, average / count, i, energyValues[i], count);
		}
		
		else
		{
			average = average + energyValues[i];
			count++;
			//printf("\tsum : %.2f\taverage : %.2fkeV\t[%d]: %.2f\t count: %d\n",average, average / count, i, energyValues[i], count);
		}
		
	}
	//printf("average = %.2f / %d = %.2f\n", average, count, average / count);	
	average = average / count;
	return average;
}

float sampledEnergies ( float ionDistance[], float ionEnergy[], float closeEnergy[], int collisionNum) //simulates the position and time and then draws an energy value from the text file.
{												   

float 	simPos[arraySize];//simulated position according to sample rate and chamber length
float	samTime[arraySize];//simulated time according to sample rate and chamber length
float	anodeDistance;

	for (int i = 0; i < maxSamples; i++)
	{
		samTime[i] = i * adcSampleRate; //simulates the time due to the sample rate and length of the chamber
		simPos[i] = (maxDepth - (driftVel * samTime[i])) * conver;//simulates the position from the time and drift velocity.
		//printf("%.2f\t%.2f\n",samTime[i],simPos[i]);
	
		int loopcount = 0;

		if ( simPos[i] > minDepth)//to account for a minimum active depth (in most cases the mylar window
		{	
			for (loopcount; loopcount < collisionNum; loopcount++)//searches through the file for an energy corresponding to a value close to the simulated position
			{
				anodeDistance = (maxDepth * conver) - ionDistance[loopcount];
		
				if (simPos[i] > anodeDistance)//once value is found, terminates loop and writes energy to array
				{
					closeEnergy[i] = ionEnergy[loopcount];
					//printf("%.2f\t%.2f\n",anodeDistance,closeEnergy[i]);
					loopcount = collisionNum;
				}
			}
		}

		else//removes points below minDepth by ending the loop
		{
			for (int j = i; j < maxSamples; j++)//fills remaining points in the array before terminating the loop
			{
				closeEnergy[j] = closeEnergy[i - 1];
			}			

			i = maxSamples;
		}
		//printf("%.2f\t%.2f\t%f\n",samTime[i],simPos[i],closeEnergy[i]);	
	}
} 

float energyFilters (float energyArray[], float peakt, int gapt, int nsamples)//long or short filter is determined from peak and gap times
{

float	difference = 0;
float 	maxValue = 0;
float	window1;
float 	window2;	

	for (int i = - 5; i < nsamples; i++)//subtracted a value to give more complete samples on either side of the curve
	{	
		window1 = average(energyArray, i, peakt);
		window2 = average(energyArray, i + gapt, peakt);

		difference = window1 - window2;

		//printf("%f\t%f\n",energyArray[i],lwindow1);
		//printf("window1 - window2 = %.2f - %.2f = %.2f\n", window1, window2, difference);

		if (difference > maxValue)
		{
			maxValue = difference;
			//printf("%.2f\n", maxValue);
		}
	}	
	return maxValue;
}

float processSRIMData (istream& srimCollisionDataSource, float shortFilter[], float longFilter[], string ionName[], float ionMass[], float initialEnergy[], int ionNum[], int *totalIonsP, int *isotopeCountP)
{//ionMass, ionNum, initialEnergy and ionName written to arrays to account for multiple isotopes (ie., files)

float 	X[arraySize];//array of positions from text file
float 	E[arraySize];//array of energies from text file
float	closeEnergy[arraySize];//energies drawn from the text file corresponding to the decided distances

int 	colNum = 0; //number of collisions
int 	index, ionCount;
int	maxNum = 1;
bool 	readIon = false;
string 	line;

	while ( getline (srimCollisionDataSource, line) )
	{
		if ( (line.length() > 30) && (line.substr (6,8).compare("Ion Name") == 0)) //determining element
		{
			//(*ionNameP) = line.substr (23,2);
			ionName[*isotopeCountP] = line.substr (23,2);//array to allow for multiple isotopes (ie., using cat command)
		}

		if ( (line.length() > 30) && (line.substr (6,8).compare("Ion Mass") == 0)) //determining ion mass
		{
			//(*ionMassP) = strtof((line.substr (22,7)).c_str(),NULL);
			ionMass[*isotopeCountP] = strtof((line.substr (22,7)).c_str(),NULL);
		}

		if ( (line.length() > 30) && (line.substr (6,10).compare("Ion Energy") == 0)) //determining initial energy
		{
			//(*initialEnergyP) = strtof((line.substr (18,11)).c_str(),NULL);
			initialEnergy[*isotopeCountP] = strtof((line.substr (18,11)).c_str(),NULL);
			*totalIonsP = *totalIonsP + ionNum[*isotopeCountP]; //keeps track of the total number of ions for storing purposes	
			*isotopeCountP += 1;//keeps track of total number of isotopes in the system
		}
		
		if ( line.substr (1,1).compare("0") == 0 || line.substr (1,1).compare("1") == 0 ) //identifying an ion line.
		{
			readIon = true;
			ionCount = strtof((line.substr(1,5)).c_str(),NULL);
			//printf("%d\n",colNum + 1);

			X[colNum]  = atof ( (line.substr (17,10)).c_str() );
			E[colNum]  = atof ( (line.substr (7,9)).c_str() ) ;

			//printf("%.2f\t%.2f\n",X[colNum],E[colNum]);

			//ionCount = ionCompare;
			colNum += 1;
		}
		
		else
		{
			if ( readIon == true ) // when end of ion is reached
			{

				sampledEnergies(X, E, closeEnergy, colNum);//used to simulate closeEnergy and store appropriate values within closeEnergy

				
				index = ionCount + *totalIonsP;	//stores each ion in a unique location in the short/long filters			

				shortFilter[index] = energyFilters(closeEnergy, s_peakt, (int)(s_gapt / adcSampleRate), maxSamples);//applies the long filter 
				longFilter[index] = energyFilters(closeEnergy, l_peakt, (int)(l_gapt / adcSampleRate), maxSamples);//applies the short filter
				ionNum[*isotopeCountP] += 1;//keeping track of number of ions in the current file
		
				//printf("%d\t%f\t%f\t%d\n",index, shortFilter[index], longFilter[index], ionNum[*isotopeCountP]);

				readIon = false;
				//(*ionNumP) += 1;
				colNum = 0;
			}
			else 
			{
				readIon = false;
			}
		}	
	}

}

float printValues (string ionName[], float ionMass[], float initialEnergy[], float shortFilter[], float longFilter[], int ionNum[], int isotopeCount)
{//goes through all isotopes and prints both diagnostic values and filters

float 	averageShort;
float	averageLong;
int 	index, total = 0;//used to keep track of the total number of ions over all the files

for ( int k = 0; k < isotopeCount; k++) //loops over values to allow for each initial conditions to be printed.
{ 

averageShort = 0;
averageLong = 0;

cout <<"#Ion Name : " << ionName[k] << endl;
printf("#Ion Mass : %.3f u\n", ionMass[k]); 
printf("#Initial Energy : %.0f keV\n",initialEnergy[k]);
printf("#Energy per nucleon : %.3f keV/u\n", initialEnergy[k] / ionMass[k]);
//fprintf("#Chamber Length : %f cm\n",maxLength);
	
		for (int j = 1; j < ionNum[k + 1] + 1; j++) //averages and prints the short and long values
		{
			index =	j + total;//the unique index where each value is stored

			averageShort += shortFilter[index];
			averageLong  += longFilter[index];
			printf("%.2f, %.2f\n", longFilter[index], shortFilter[index]);
		}

averageShort = averageShort / ionNum[k + 1];// it is ionNum-1 since ionNum is originally defined as 1 instead of 0.
averageLong = averageLong / ionNum[k + 1];
total = total + ionNum[k + 1]; //used to keep track of the total number of ions over all the files

printf("#There are %d ions in this file.\n", ionNum[k + 1]); //This assumes that there is the same number of ions in each file.
printf("#AVG SHORT FILTER : %.2fkeV\n", averageShort);	
printf("#AVG LONG FILTER : %.2fkeV\n\n", averageLong);

}
}

int setConstants (char* itemName[], int *iP, float *adcSampleRateP, float *l_peaktP, float *l_gaptP, float *s_peaktP, float *s_gaptP, float *driftVelP, float *minDepthP, float *maxDepthP)
{//takes in constant values from command lines and sets them appropriately

*iP += 1;//increases i to skip next loop as well as store value that comes afterwards	
string parameter = itemName[*iP - 1];//converts the argv into a string again
	if ( parameter.substr(2,13) == "adcSampleRate")
	{
		*adcSampleRateP = atof(itemName[*iP]);
		//printf("%f\n",*adcSampleRateP);
	}

	if ( parameter.substr(2,7) == "l_peakt")
	{
		*l_peaktP = atof(itemName[*iP]);
		//printf("%f\n",*l_peaktP);
	}

	if ( parameter.substr(2,6) == "l_gapt")
	{
		*l_gaptP = atof(itemName[*iP]);
		//printf("%f\n",*l_gaptP);
	}

	if ( parameter.substr(2,7) == "s_peakt")
	{
		*s_peaktP = atof(itemName[*iP]);
		//printf("%f\n",*s_peaktP);
	}

	if ( parameter.substr(2,6) == "s_gapt")
	{
		*s_gaptP = atof(itemName[*iP]);
		//printf("%f\n",*s_gaptP);
	}

	if ( parameter.substr(2,8) == "driftVel")
	{
		*driftVelP = atof(itemName[*iP]);
		//printf("%f\n",*driftVelP);
	}

	if ( parameter.substr(2,8) == "minDepth")
	{
		*minDepthP = atof(itemName[*iP]);
		//printf("%f\n",*minDepthP);
	}

	if ( parameter.substr(2,8) == "maxDepth")
	{
		*maxDepthP = atof(itemName[*iP]);
		//printf("%f\n",*maxDepthP);
	}
}

int main(int argc, char* argv[])
{

float	shortFilter[arraySize];//returned short values from filters for each ion
float 	longFilter[arraySize];//returned long values from filters for each ion

int 	ionNum[arraySize]; //number of ions (total value over all files in the case of cat command)
int	isotopeCount = 0; //used to count number of isotopes when using the cat command
int	totalIons = 0;//used to keep track of the total ions over all of the files. Needs to be set here so it doesn't reset when using "./TBraggSimulation" input method (see arg loop)
string 	fileName;
float 	initialEnergy[arraySize]; //initial energy for all of the ions
float 	ionMass[arraySize]; //mass of the ion
string 	ionName[arraySize]; //name of the ion (element)

	if (argc > 1) //reads through the number of files after the executable call of format "./TBraggSimulation Si28.txt Mo95.txt ..."
	{
		for (int i = 1; i < argc; i++)//looping over each argument in the above executable call
		{			
			string argv_ = argv[i]; //converts the argument to a string so it is easier to manipulate 
		
			if (strcmp (argv_.substr(0,2).c_str(), "--") == 0) //determining if there is a command to set a parameter
			{
			 	setConstants(argv, &i, &adcSampleRate, &l_peakt, &l_gapt, &s_peakt, &s_gapt, &driftVel, &minDepth, &maxDepth); 
				//accesses subroutine to define parameter called upon. More parameters can be easily added to this subroutine.		
			}
	
			else
			{ 
				fileName = argv_;	
				ifstream myfile; 	
				myfile.open(fileName.c_str());
		
				if( myfile.is_open() ) //reads that the file is open and completes the subroutine
				{
					processSRIMData(myfile, shortFilter, longFilter, ionName, ionMass, initialEnergy, ionNum, &totalIons, &isotopeCount);
				}
			}	
		}
		printValues(ionName, ionMass, initialEnergy, shortFilter, longFilter, ionNum, isotopeCount);//to allow for loop to perform properly with stored values in the array
	}

	else //processes and prints values from the line "cat filename (multiple) | ./TBraggSimulation"
	{
		processSRIMData(cin, shortFilter, longFilter, ionName, ionMass, initialEnergy, ionNum, &totalIons, &isotopeCount);
		printValues(ionName, ionMass, initialEnergy, shortFilter, longFilter, ionNum, isotopeCount);
	}
}

















