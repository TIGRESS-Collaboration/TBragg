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
int main()
{

const float adcSampleRate = 0.02;//sample rate in microseconds
const float l_peakt = 0.5; //long peak time in microseconds
const float l_gapt = 4.5; //long gap time in microseconds
const float s_peakt = 0.14; //short peak tme in microseconds
const float s_gapt = 0.24; // short gap time in microseconds
const float mylarSize = 1.5; //mylar size in micrometers
const float driftVel = 4.8;//drift velocity in cm/micron
const float chamberLength = 10; //chamber length in centimetres
const int   conver = 100000000; //converting angstroms to centimeters
const int   arraySize = 1000; //size of all arrays in this program
const int   maxSamples = chamberLength / (driftVel * adcSampleRate);
const int   maxSamples2 = chamberLength * ((int)(1 / (driftVel * adcSampleRate)));

printf("%d\n",maxSamples);
printf("%d\n",maxSamples2);

int gapt = (int)(s_gapt / adcSampleRate);
float gapt2 =(int)(s_gapt / adcSampleRate);

printf("%d\n", gapt);
printf("%f\n",gapt2);

}
