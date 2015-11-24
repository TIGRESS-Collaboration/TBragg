//g++ -o TBraggSimulation_v1.exe TBraggSimulation_v1.cpp
#include <iostream>
#include <algorithm>
#include <string>
#include <fstream>
#include <vector>
#include <stdlib.h>
#include <cstring>
#include <math.h>
#include <cmath>

using namespace std;
const float adcSampleRate = 0.02;//sample rate in microseconds
const float l_peakt = 0.5; //long peak time in microseconds
const float l_gapt = 4.5; //long gap time in microseconds
const float s_peakt = 0.14; //short peak tme in microseconds
const float s_gapt = 0.24; // short gap time in microseconds
const float mylarSize = 1.5; //mylar size in micrometers
const float driftVel = 4.8;//drift velocity in cm/micron
const float chamberLength = 10; //chamber length in centimetres

int main ( )
{

float samTime[1000];
float simPos[1000];
int maxSamples = chamberLength / (driftVel * adcSampleRate);

for (int i = 0; i < maxSamples; i++)
{
	samTime[i] = i * adcSampleRate; //simulates the time due to the sample rate and length of the chamber
	simPos[i] = chamberLength - (driftVel * samTime[i]);//simulates the position from the time and drift velocity.
	printf("%.2f\t%.2f\n",samTime[i],simPos[i]);
}

/*for (int j = 0; j < maxSamples; j++)
{
	simPos[j] = chamberLength - (driftVel * samTime[j]);//simulates the position from the time and drift velocity.
	printf("%.2f\n",simPos[j]);
}*/

} 
