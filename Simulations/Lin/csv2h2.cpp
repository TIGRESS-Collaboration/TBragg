


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
float maxx = 500000;
float maxy = 125000; // Some defautl values to be changed later
int xbins  = 1024;
int ybins  = 256;

int main(int argc, const char* argv[]) { // Not sure if I'll ever use this
	TApplication *app = new TApplication("app",0,0);
	TFile* outFile = new TFile("graph.root", "RECREATE");
	TCanvas *c1 = new TCanvas("c1","TBRAGG Simulation",200,10,700,500);
	TH2I* histo;
	TGraph* graph; // May upgrade to multiple graphs later
	// gStyle->SetPalette(1);
	
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

	histo = new TH2I("PID","Particle Identification Plot",xbins,minx,maxx,ybins,miny,maxy);
	graph = new TGraph(); // hope this works
	graph->SetMarkerStyle(8);
	graph->SetMarkerSize(0.2);
	graph->SetMarkerColor(3);
	graph->GetXaxis()->SetTitle("Long Filter (MeV)");
	graph->GetYaxis()->SetTitle("Short Filter (MeV)");
	graph->GetXaxis()->CenterTitle();
	graph->GetYaxis()->CenterTitle();
	graph->SetTitle("Particle Identification Plot");

	
	{ // Put this in braces so I can break it out somewhere else eventually
		string line;
		while (getline(cin,line)) {
			float x;
			float y;
			char* cpline = new char [line.length()+1];
			strcpy(cpline,line.c_str());
			char* pcy;
			char* pcx; // after commas and/or white space
			pcx = strtok(cpline,",");
			if (pcx!=0) {
				pcy = strtok(NULL,",");
				//printf("%s\t%s\t%s\n",pcx,pcy, cpline);
				if (pcy!=0) { // If this is a legitmate thing with two entries
					x = atof(pcx);
					y = atof(pcy);
					// This needs to be thought through a bit better ....d
					histo->Fill(x,y);
					graph->SetPoint(graph->GetN(),x,y);
				} // end of pcy!=0 condition
			} // end of pcx!=0 condition
		} // end of while getline condition
	}
	histo->DrawCopy("COLZ");
	graph->SetMarkerStyle(8);
	graph->SetMarkerSize(1);
	//graph->SetMarkerColorAlpha(1,0.5);
	//graph->GetXaxis()->SetTitle("Long Filter (MeV)");
	//graph->GetYaxis()->SetTitle("Short Filter (MeV)");
	//graph->GetXaxis()->CenterTitle();
	//graph->GetYaxis()->CenterTitle();
//	graph->Draw("PSAME");
	c1->Draw();
	c1->SaveAs("blah");
	outFile->Write();
	outFile->Close();

	app->Run(0);

	return 0;
}

	


	
