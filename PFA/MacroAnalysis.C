#include <TH2.h>
#include <TStyle.h>
#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>

void MacroAnalysis() 
{
    fstream myfile;
    myfile.open("datas.txt",ios::in);

    std::vector<double> gapSize;
    std::vector<double> purity;
    std::vector<double> efficiency;

    if ( myfile.is_open() ) 
    {
        double data = 0.0;
        while ( myfile >> data ) 
        {
            gapSize.push_back((double)data);
            myfile >> data;
            purity.push_back((double)data);
            myfile >> data;
            efficiency.push_back((double)data);
        }
    }
    int const n = gapSize.size();

    TCanvas *c1 = new TCanvas("c1","Purity",200,10,500,300);
    double x[n], ypurity[n], yefficiency[n];

    for (int i = 0; i < n; ++i)
    {
        x[i] = gapSize[i];
        yefficiency[i] = efficiency[i];
        ypurity[i] = purity[i];
        
    }

    double ypurityR[6];
	ypurityR[0] = 37.5;
	ypurityR[1] = 73.;
	ypurityR[2] = 85.;
	ypurityR[3] = 88.;
	ypurityR[4] = 90.;
	ypurityR[5] = 91.;

	double yefficiencyR[6];
	yefficiencyR[0] = 38.;
	yefficiencyR[1] = 75.;
	yefficiencyR[2] = 88.;
	yefficiencyR[3] = 94.;
	yefficiencyR[4] = 96.;
	yefficiencyR[5] = 97.;

    TCanvas *c_pl = new TCanvas("c","Graph2D example",0,0,1200,600);

    c_pl->Divide(2,1);

    c_pl->cd(1);
    TGraph* gr = new TGraph(n,x,ypurity);
    TAxis *axis = gr->GetXaxis();
    TAxis *axisY = gr->GetYaxis();

    axis->SetLimits(0.,40.);   
    axisY->SetLimits(0.,100.);              // along X 
    gr->GetHistogram()->SetMaximum(101.);   // along Y        
    gr->GetHistogram()->SetMinimum(0.);    

    
    gr->SetLineColor(2);
    gr->SetLineWidth(1.5);
    gr->SetMarkerColor(4);
    gr->SetMarkerSize(1.);
    gr->SetMarkerStyle(21);


	TGraph* grR = new TGraph(n,x,ypurityR);
	grR->SetLineColor(1);
    grR->SetLineWidth(1.5);
    grR->SetMarkerColor(3);
    grR->SetMarkerSize(1.);
    grR->SetMarkerStyle(21);
grR->GetHistogram()->SetMaximum(101.);   // along Y        
    grR->GetHistogram()->SetMinimum(0.);   

    TMultiGraph *mg = new TMultiGraph();
    mg->Add(gr,"lp");
    mg->Add(grR,"lp");
    mg->SetTitle("efficiency for neutral particle");
    mg->Draw("a");

  

    c_pl->cd(2);
    TGraph* gr2 = new TGraph(n,x,yefficiency);
    TAxis *axis2 = gr2->GetXaxis();
    TAxis *axisY2 = gr2->GetXaxis();

    axis2->SetLimits(0.,40.);     
    axisY2->SetLimits(0.,100.);            // along X 
    gr2->GetHistogram()->SetMaximum(101.);   // along Y        
    gr2->GetHistogram()->SetMinimum(0.);    

    
    gr2->SetLineColor(2);
    gr2->SetLineWidth(1.5);
    gr2->SetMarkerColor(4);
    gr2->SetMarkerSize(1.);
    gr2->SetMarkerStyle(21);

    TGraph* gr2R = new TGraph(n,x,ypurityR);
	gr2R->SetLineColor(1);
    gr2R->SetLineWidth(1.5);
    gr2R->SetMarkerColor(3);
    gr2R->SetMarkerSize(1.);
    gr2R->SetMarkerStyle(21);
gr2R->GetHistogram()->SetMaximum(101.);   // along Y        
    gr2R->GetHistogram()->SetMinimum(0.);   

    TMultiGraph *mg2 = new TMultiGraph();
    mg2->SetTitle("purity for neutral particle");
    mg2->Add(gr2,"lp");
    mg2->Add(gr2R,"lp");
    mg2->Draw("a");
}
