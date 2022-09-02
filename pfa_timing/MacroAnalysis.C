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
            efficiency.push_back((double)data);
            myfile >> data;
            purity.push_back((double)data);
        }
    }
    int const n = gapSize.size();

    //TCanvas *c1 = new TCanvas("c1","Purity",200,10,500,300);
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
	yefficiencyR[0] = 36.;
	yefficiencyR[1] = 74.;
	yefficiencyR[2] = 88.;
	yefficiencyR[3] = 92.;
	yefficiencyR[4] = 95.;
	yefficiencyR[5] = 97.;

    double ypurityWithoutTiming[6];
    ypurityWithoutTiming[0] = 32.9126;
    ypurityWithoutTiming[1] = 78.0036;
    ypurityWithoutTiming[2] = 87.0555;
    ypurityWithoutTiming[3] = 91.0181;
    ypurityWithoutTiming[4] = 91.0181;
    ypurityWithoutTiming[5] = 96.1097;

    double yefficiencyWithoutTiming[6];
    yefficiencyWithoutTiming[0] =  34.4264;
    yefficiencyWithoutTiming[1] =  64.3777;
    yefficiencyWithoutTiming[2] =  78.213;
    yefficiencyWithoutTiming[3] =  82.4602;
    yefficiencyWithoutTiming[4] =  88.3716;
    yefficiencyWithoutTiming[5] =  88.3716;

    double V1P[6];
    V1P[0] = 50.;
    V1P[1] = 75.;
    V1P[2] = 85.;
    V1P[3] = 87.;
    V1P[4] = 88.;
    V1P[5] = 88.;

    double V1E[6];
    V1E[0] = 60.;
    V1E[1] = 75.;
    V1E[2] = 85.;
    V1E[3] = 90.;
    V1E[4] = 92.;
    V1E[5] = 92.;




    TCanvas *c_pl = new TCanvas("c","Graph2D example",0,0,1200,600);

    c_pl->Divide(2,1);

    c_pl->cd(1);
    TGraph* gr = new TGraph(n,x,ypurity);
    TAxis *axis = gr->GetXaxis();
    TAxis *axisY = gr->GetYaxis();



    
    gr->SetLineColor(2);
    gr->SetLineWidth(1.5);
    gr->SetMarkerColor(2);
    gr->SetMarkerSize(1.);
    gr->SetMarkerStyle(21);


	TGraph* grR = new TGraph(n,x,ypurityR);
	grR->SetLineColor(1);
    grR->SetLineWidth(1.5);
    grR->SetMarkerColor(1);
    grR->SetMarkerSize(1.);
    grR->SetMarkerStyle(21);
    grR->GetHistogram()->SetMaximum(101.);   // along Y        
    grR->GetHistogram()->SetMinimum(0.);

    TGraph* grWT = new TGraph(n,x,ypurityWithoutTiming);
    grWT->SetLineColor(4);
    grWT->SetLineWidth(1.5);
    grWT->SetMarkerColor(4);
    grWT->SetMarkerSize(1.);
    grWT->SetMarkerStyle(21);
    grWT->GetHistogram()->SetMaximum(101.);   // along Y        
    grWT->GetHistogram()->SetMinimum(0.);    

    // TGraph* aaa = new TGraph(n,x,V1P);
    // aaa->GetHistogram()->SetLineColor(3);
    // aaa->SetLineWidth(1.5);
    // aaa->SetMarkerColor(3);
    // aaa->SetMarkerSize(1.);
    // aaa->SetMarkerStyle(21);
    // aaa->GetHistogram()->SetMaximum(101.);   // along Y        
    // aaa->GetHistogram()->SetMinimum(0.);  

    TMultiGraph *mg = new TMultiGraph();
    mg->Add(gr,"lp");
    mg->Add(grR,"lp");
    mg->Add(grWT,"lp");
    // mg->Add(aaa,"lp");

    mg->SetTitle("purity for neutral particle");
    mg->Draw("AC*");

    axis->SetLimits(0.,40.);   
    axisY->SetLimits(0.,101.);              // along X 
    gr->GetHistogram()->SetMaximum(101.);   // along Y        
    gr->GetHistogram()->SetMinimum(0.);    

    auto legend = new TLegend();
    legend->AddEntry(gr,"with timing", "p");
    legend->AddEntry(grWT,"without timing", "p");
    // legend->AddEntry(aaa,"Goal", "p");
    legend->AddEntry(grR,"Remi's data", "p");

    legend->Draw();


  

    c_pl->cd(2);
    TGraph* gr2 = new TGraph(n,x,yefficiency);
    TAxis *axis2 = gr2->GetXaxis();
    TAxis *axisY2 = gr2->GetYaxis();

    
    gr2->SetLineColor(2);
    gr2->SetLineWidth(1.5);
    gr2->SetMarkerColor(2);
    gr2->SetMarkerSize(1.);
    gr2->SetMarkerStyle(21);

    TGraph* gr2R = new TGraph(n,x,yefficiencyR);
	gr2R->SetLineColor(1);
    gr2R->SetLineWidth(1.5);
    gr2R->SetMarkerColor(1);
    gr2R->SetMarkerSize(1.);
    gr2R->SetMarkerStyle(21);
    gr2R->GetHistogram()->SetMaximum(101.);   // along Y        
    gr2R->GetHistogram()->SetMinimum(0.);   

    TGraph* gr2WT = new TGraph(n,x,yefficiencyWithoutTiming);
    gr2WT->SetLineColor(4);
    gr2WT->SetLineWidth(1.5);
    gr2WT->SetMarkerColor(4);
    gr2WT->SetMarkerSize(1.);
    gr2WT->SetMarkerStyle(21);
    gr2WT->GetHistogram()->SetMaximum(101.);   // along Y        
    gr2WT->GetHistogram()->SetMinimum(0.);    

    // TGraph* bbb = new TGraph(n,x,V1E);
    // bbb->GetHistogram()->SetLineColor(3);
    // bbb->SetLineWidth(1.5);
    // bbb->SetMarkerColor(3);
    // bbb->SetMarkerSize(1.);
    // bbb->SetMarkerStyle(21);
    // bbb->GetHistogram()->SetMaximum(101.);   // along Y        
    // bbb->GetHistogram()->SetMinimum(0.);  

    TMultiGraph *mg2 = new TMultiGraph();
    mg2->SetTitle("efficiency for neutral particle");
    mg2->Add(gr2,"lp");
    mg2->Add(gr2R,"lp");
    mg2->Add(gr2WT,"lp");
    // mg2->Add(bbb,"lp");
    mg2->Draw("AC*");

    axis2->SetLimits(0.,40.);     
    axisY2->SetLimits(0.,100.);            // along X 
    gr2->GetHistogram()->SetMaximum(101.);   // along Y        
    gr2->GetHistogram()->SetMinimum(0.);    
}
