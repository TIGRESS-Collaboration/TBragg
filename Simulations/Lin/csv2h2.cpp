
/**
 * g++ csv2h2.cpp -o csv2h2 `root-config --cflags --libs`    ----------  compile line
 * Must be run with data piped in from TBragg Simulation:
 * ./TBraggSimulation <collision file> | ./csv2h2
*/


/**
 * File: csv2h2.cpp
 * Last Modified By: Alex Kurkjian
 * Date: 10/08/2016
 * Purpose: This program simply pulls the data from the output of the TBragg simulation and plots
 * 	on a histogram using ROOT.
**/

#include <Rtypes.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TApplication.h>
#include <TH2.h>
#include <TF2.h>
#include <TGraph.h>
#include <TMarker.h>
#include <TLegend.h>

#include <iostream>
#include <algorithm>
#include <string>
#include <fstream>
#include <vector>
#include <stdlib.h>
#include <cstring>
#include <math.h>
#include <string>

using namespace std;

float minx = 0;
float miny = 0;
float maxx = 200000;
float maxy = 100000; // Some default values which may be adjusted later
int xbins  = 5000;
int ybins  = 5000;

int main(int argc, const char* argv[]) { 
	TApplication *app = new TApplication("app",0,0);
	TFile* outFile = new TFile("graph.root", "RECREATE");
	TCanvas *c1 = new TCanvas("c1","TBRAGG Simulation",200,10,700,500);
	TH2F* histo;

	
	if (argc == 2)
	{
		maxx = atof(argv[1]);
		maxy = atof(argv[1]);
	}	

	if (argc == 3)
	{
		maxx = atof(argv[1]);
		maxy = atof(argv[2]);
	}	

	if (argc == 4)
	{
		maxx = atof(argv[1]);
		maxy = atof(argv[2]);
		xbins = atoi(argv[3]);
		ybins = atoi(argv[4]);
	}	

	if (argc == 5)
	{
		maxx = atof(argv[1]);
		maxy = atof(argv[2]);
		xbins = atoi(argv[3]);
		ybins = atoi(argv[4]);
	}

	histo = new TH2F("PID","Particle Identification Plot",xbins,minx,maxx,ybins,miny,maxy);
	//histo->SetXTitle("Long Filter"); 
	//histo->SetYTitle("Short Filter");
	histo->GetXaxis()->CenterTitle();
	histo->GetYaxis()->CenterTitle();
	
	{ // Put this in braces so I can break it out somewhere else eventually
		string line;
		while (getline(cin,line)) {
			float x;
			float y;
			char* cpline = new char [line.length()+1];
			strcpy(cpline,line.c_str());
			char* pcy;
			char* pcx; // After commas and/or white space
			pcx = strtok(cpline,",");
			if (pcx!=0) {
				pcy = strtok(NULL,",");
				//printf("%s\t%s\t%s\n",pcx,pcy, cpline);
				if (pcy!=0) { // If the data entry is legitimate with x and y enties
					x = atof(pcx);
					y = atof(pcy);
					histo->Fill(x,y);
				} // End of pcy!=0 condition
			} // End of pcx!=0 condition
		} // End of while getline condition
	}
	histo->DrawCopy("COLZ");
	c1->Draw();
	c1->SaveAs("PID_canvas.ps");
	outFile->Write();
	outFile->Close();

	app->Run(0);

	return 0;
}

	


	
