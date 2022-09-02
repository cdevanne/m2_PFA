#include <TH2.h>
#include <TStyle.h>
#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>

void MacroAnalysis2() 
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
    //perfect timing
    double ypurityR[6];
	ypurityR[0] = 51.0178;
	ypurityR[1] = 71.1151;
	ypurityR[2] = 82.4624;
	ypurityR[3] = 89.1579;
	ypurityR[4] = 95.752;
	ypurityR[5] = 95.752;

	double yefficiencyR[6];
	yefficiencyR[0] = 50.9506   ;
	yefficiencyR[1] = 71.481 ;
	yefficiencyR[2] = 81.5025;
	yefficiencyR[3] = 87.1525 ;
	yefficiencyR[4] = 91.1384;
	yefficiencyR[5] = 91.1384;

  
    // 100ps
    double ypurityWithoutTiming[6];
    ypurityWithoutTiming[0] = 49.5959;
    ypurityWithoutTiming[1] = 71.0698;
    ypurityWithoutTiming[2] = 83.2603;
    ypurityWithoutTiming[3] = 90.0437;
    ypurityWithoutTiming[4] = 94.9026;
    ypurityWithoutTiming[5] = 94.9026;

    double yefficiencyWithoutTiming[6];
    yefficiencyWithoutTiming[0] = 47.157;
    yefficiencyWithoutTiming[1] = 72.6219;
    yefficiencyWithoutTiming[2] = 82.2274;
    yefficiencyWithoutTiming[3] = 87.6613;
    yefficiencyWithoutTiming[4] = 90.0509;
    yefficiencyWithoutTiming[5] = 90.0509;


 


    // 1ns
    double V1P[6];
    V1P[0] = 27.1953;
    V1P[1] = 43.8638;
    V1P[2] = 44.298;
    V1P[3] = 54.2774;
    V1P[4] = 52.1929;
    V1P[5] = 52.1929;

    double V1E[6];
    V1E[0] = 29.1426 ;
    V1E[1] = 45.0306 ;
    V1E[2] = 46.5951;
    V1E[3] = 56.1832;
    V1E[4] = 52.8783;
    V1E[5] =  52.8783;




    TCanvas *c_pl = new TCanvas("c","Graph2D example",0,0,1200,600);

    c_pl->Divide(2,1);

    c_pl->cd(1);
    TGraph* gr = new TGraph(n,x,ypurity);
    TAxis *axis;
    TAxis *axisY; 



    
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

    TGraph* aaa = new TGraph(n,x,V1P);
    aaa->GetHistogram()->SetLineColor(3);
    aaa->SetLineWidth(1.5);
    aaa->SetMarkerColor(3);
    aaa->SetMarkerSize(1.);
    aaa->SetMarkerStyle(21);
    aaa->GetHistogram()->SetMaximum(101.);   // along Y        
    aaa->GetHistogram()->SetMinimum(0.);  

    TMultiGraph *mg = new TMultiGraph();
    mg->Add(gr,"lp");
    mg->Add(grR,"lp");
    mg->Add(grWT,"lp");
    mg->Add(aaa,"lp");

    mg->SetTitle("purity for neutral particle");
    mg->Draw("AC*");

    mg->GetXaxis()->SetLimits(3.,32.);  
    mg->SetMinimum(0.);
    mg->SetMaximum(101.);

    mg->GetXaxis()->SetTitle("distance between particles (mm)");
    mg->GetYaxis()->SetTitle("purity (%)");


    auto legend = new TLegend();
    legend->AddEntry(gr,"500ps", "p");
    legend->AddEntry(grWT,"100ps", "p");
    legend->AddEntry(aaa,"1ns", "p");
    legend->AddEntry(grR,"perfect timing", "p");

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

    TGraph* bbb = new TGraph(n,x,V1E);
    bbb->GetHistogram()->SetLineColor(3);
    bbb->SetLineWidth(1.5);
    bbb->SetMarkerColor(3);
    bbb->SetMarkerSize(1.);
    bbb->SetMarkerStyle(21);
    bbb->GetHistogram()->SetMaximum(101.);   // along Y        
    bbb->GetHistogram()->SetMinimum(0.);  

    TMultiGraph *mg2 = new TMultiGraph();
    mg2->Add(gr2,"lp");
    mg2->Add(gr2R,"lp");
    mg2->Add(gr2WT,"lp");
    mg2->Add(bbb,"lp");

    mg2->SetTitle("efficiency for neutral particle");
    

    mg2->GetXaxis()->SetLimits(3.,32.);  
    mg2->SetMinimum(0.);
    mg2->SetMaximum(101.);

    mg2->GetXaxis()->SetTitle("distance between particles (mm)");
    mg2->GetYaxis()->SetTitle("efficiency (%)");

    mg2->Draw("AC*");
}
