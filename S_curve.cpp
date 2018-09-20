#include <iostream>
#include <fstream>
#include <vector>
#include <cstring>
#include <stdio.h>
#include <iomanip>
#include <typeinfo>

#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include <TH1F.h>
#include <TH1D.h>
#include <TF1.h>
#include <TFitResultPtr.h>
#include <TFitResult.h>
#include <TMath.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TString.h>

using namespace std;

int num_histos(const char *fname) { //calculate number of histograms in the root file
	TFile *f = TFile::Open(fname, "READ");
	int total = 0; 
	TIter next(f->GetListOfKeys());  //iterator for keys
    while ((TKey *)next()) {
	 	total++;
	}
	return total;
};

double charge_transf(int val){ //charge in ke^-
	double coefficient = 15.67;
	double charge = coefficient * val;	
	return charge;
}

template <typename T> string NumberToString ( T Number ) //make string out of any number entry
  {
     ostringstream ss;
     ss << Number;
     return ss.str();
  }


//logistic sigmoid functions
Double_t logisticF(Double_t *x,Double_t *par){
	Double_t func = 1/(1+TMath::Exp(-par[0]*(x[0]-par[1])));
	return func;

}


//error functions
Double_t erfF(Double_t *x,Double_t *par){
	Double_t func = (1+TMath::Erf((x[0]-par[0])/par[1]))/2.0;
	return func;
}


//calculation of errors for Middle point calculation
pair<double, double> errorXMV(double x,  double x1, double x2  ){
	return make_pair(x-x1, x2-x); //down, up
}

//calculation of errors for fits if position of transiotion is different from 0.5
pair<double, double> errorXlogistic( TF1 *name, double x, double level){
	double par0 = name -> GetParameter(0);
	double par1 = name -> GetParameter(1);


	Double_t err0 = name -> GetParError(0);
	Double_t err1 = name -> GetParError(1);

	double errorY = TMath::Sqrt(TMath::Power(-1*(TMath::Exp(-par0*err0*(-par1+ x))*(par1 - x))/TMath::Power(1 + TMath::Exp(-par0*(-par1+ x)),2),2) + TMath::Power( -err1*(par0*TMath::Exp(-par0*(-par1 + x)))/TMath::Power(1 + TMath::Exp(-par0*(-par1+ x)),2),2)); 

	double valueX = name -> GetX(level);

	double valueX1 = name -> GetX(level + errorY);
	double valueX2 = name -> GetX(level - errorY);

	return make_pair(valueX1-valueX, valueX-valueX2);
}

pair<double, double> errorXerf( TF1 *name, double x, double level){
	double par0 = name -> GetParameter(0);
	double par1 = name -> GetParameter(1);


	Double_t err0 = name -> GetParError(0);
	Double_t err1 = name -> GetParError(1);
	//cout << err0 << " " << err1 << endl;
	double errorY = TMath::Sqrt(TMath::Power(err0*TMath::Exp((par0-x)*(-par0+x)/(par1*par1))/(par1*TMath::Sqrt(TMath::Pi())), 2)+ TMath::Power(err1*TMath::Exp((par0 - x)*(-par0 + x)/(par1*par1))*(x-par0)/(par1*par1*TMath::Sqrt(TMath::Pi())),2));
	//cout << errorY << endl;
	
	double valueX = name -> GetX(level);

	double valueX1, valueX2;
	if(level + errorY >1 && level - errorY >0 ) valueX1 = name -> GetX(0.99);
	if(level + errorY <1 && level - errorY <0 ) valueX2 = name -> GetX(0.01);
	if(level + errorY <1 && level - errorY >0) { cout << "1" << endl; valueX1 = name -> GetX(level + errorY); valueX2 = name -> GetX(level - errorY);}
	else {cout << "2" << endl; valueX1 = name -> GetX(0.99); valueX2 = name -> GetX(0.01);}
	
	return make_pair(valueX1-valueX, valueX-valueX2); //up and down
}


 
void S_curve(){

    int offset[15]={1,1,1,1,1,1,1,1,1,1,1,1,0,0,0};
    int attenuation[15]={3,3,3,2,2,2,1,1,1,0,0,0,0,0,0};
    int thresholds[15]={42,38,36,42,52,62,42,52,62,42,52,62,42,52,62};
    string names[15]={"histo_run271", "histo_run272", "histo_run273", "histo_run98","histo_run100", "histo_run102", "histo_run64", "histo_run93", "histo_run96", "histo_run63", "histo_run59", "histo_run61","histo_run114","histo_run142","histo_run146"};

    for(int p = 0; p<2; p++){ //looping over two boards
    	string dir;
	if(p==0) dir = "/afs/cern.ch/user/v/valukash/RICH-threshold-analysis/DAC/Board1_13/";
   	if(p==1) dir = "/afs/cern.ch/user/v/valukash/RICH-threshold-analysis/DAC/Board1_06/";

	for(int q = 0; q < 15; q++){ //looping over data files

    		string end = ".root";
		string name_histo = names[q];
		string name_buf = dir+name_histo+end;
		string savedir = dir + "S_curves/";

    		TFile *in = new TFile(name_buf.c_str()); 

    		TH1F *histo_logistic = new TH1F("histo_logistic", "histo_logistic", 500, 0, 500); 
    		TH1F *histo_erf = new TH1F("histo_erf", "histo_erf", 500, 0, 500);
    		TH1D *fit_status = new TH1D("fit_status", "fit_status", 5, 0, 5); 

    		double level_th = 0.5; //choose a level of transition point
    
		int att = attenuation[q];
    		int off = offset[q];
    		int th = thresholds[q];

    		int n_trig = 10000; // number of triggers

    		int ncounts;        // number of DAC counts per DAC scan

    		ofstream file;
    		file.open(savedir +name_histo+"/"+TString::Format("thresholds.txt"), ios::out);
   
    		file << "The threshold level:" << level_th  << endl;
    		file << "Attenuation: " << att <<  endl;
    		file << "Offset: " << off << endl;
    		file << "Threshold: " << th << endl;
    		file << "Channel ID" << setw(30) << "Middle value, DAC st" << setw(30) << "Error DOWN, DAC st" << setw(30) << "Error UP, DAC st" << setw(30) << "Middle value, ke-" << setw(30) << "logistic, DAC st" << setw(30) << "Error DOWN, DAC st" << setw(30) << "Error UP, DAC st" << setw(30) << "Chi2" << setw(30) << "Fit status" <<  setw(30) << "logistic, ke" << setw(30) << "erf, DAC st" << setw(30) << "Error DOWN, DAC st" << setw(30) << "Error UP, DAC st" << setw(30) << "chi2" << setw(30) << "Fit status" <<  setw(30) << "erf, ke" << endl;    

    		int nPMT = 4; 
    		int nEC = 4;
    		int ngraphs = num_histos(name_buf.c_str())/32; //number of graphs
	
                TGraphAsymmErrors* graphs[nEC][nPMT][ngraphs];
                TCanvas *c[nEC*nPMT];

                for(int i =0; i < nEC; i++){
    	    		for(int j = 0; j < nPMT-2; j++){
      				c[4*i+j] = new TCanvas(TString::Format("c%d", i), TString::Format("c%d", i), 16000, 16000);
      				c[4*i+j] -> Divide(8, 8);
				vector <double> valuesxone;
				vector <double> valuesxzero;
				vector <double> valuesxmiddle;
				vector <int> noisecode;

    				for(int l = 1; l < ngraphs+1; l++){ 
					c[4*i+j] -> cd(l);
					graphs[i][j][l] = new TGraphAsymmErrors();
    					graphs[i][j][l] = (TGraphAsymmErrors*) in ->FindObjectAny(TString::Format("S-Curve EC%d PMT%d Ch %d_g;1", i, j, l));
					if (graphs[i][j][l] == nullptr) cout << "Message [data file] Requested object does not exist or has an incompatible type\n";
     					if (graphs[i][j][l]) cout << TString::Format("Graph named \"%s\" has the title \"%s\"", graphs[i][j][l]-> GetName(),graphs[i][j][l]-> GetTitle()) << endl;

					gStyle -> SetLineWidth(1);
    					graphs[i][j][l]-> GetYaxis() -> SetTitle("Normilized counts");
    					graphs[i][j][l]-> GetXaxis() -> SetTitle("DAC counts");

					ncounts = graphs[i][j][l]->GetN();
					//cout << "Message [DAC scans] number of steps:" << ncounts << endl;

					Double_t pointsx[ncounts];
					Double_t pointsy[ncounts];
					for(int u = 10; u < ncounts; u++) graphs[i][j][l]-> GetPoint(u, pointsx[u], pointsy[u]); //HACK: int u = 10 for bad data files where first counts are on.
					if(pointsy[10]>0.1 || pointsy[ncounts-1] < 0.9) { //preselection of empty/ON channels HACK: pointsy[10]>0.1 for bad data files, where first counts are on
						cout << TString::Format("Message [DAC scans] S-Curve EC%d PMT%d Ch %d is empty or on", i, j, l) << endl;
						noisecode.push_back(l);
						valuesxzero.push_back(0);
						valuesxone.push_back(0);
						valuesxmiddle.push_back(0);
						continue;
					}
    					else {  //channel is okay 
						for(int u = 10; u < ncounts; u++){
							if(pointsy[u]!=0) { //find last point = 0
								double value0 = pointsx[u-1]-0.5;
								valuesxzero.push_back(value0);
								//cout << "Zero point:" << value0 << " " << valuesxzero[l-1] << endl;
								break;
							}
						}
						for(int u = 10; u < ncounts; u++){ // find first point = 1
							if(pointsy[u]>0.9) {
								double value1 = pointsx[u]-0.5;
								valuesxone.push_back(value1);
								//cout <<"One point:"<< valuesxone[l-1] << endl;
								double valuem = (valuesxone[l-1]+valuesxzero[l-1])/2;
								valuesxmiddle.push_back(valuem);
								//cout << "Message [DAC scans] middle point:" << valuesxmiddle[l-1] << endl;
								break;
							}
						}

					        //cout << "Message [transition point] transition point from middle point:" << valuesxmiddle[l-1]+off*off_value << "(" << charge_transf(valuesxmiddle[l-1]) << " in ke)" << endl;
						
						auto logistic = new TF1("logistic", logisticF, 0, 256, 2); //256 - max number of DAC counts 
						auto erf = new TF1("erf", erfF, 0 , 256, 2); 
						
						logistic -> SetParameters(1., valuesxmiddle[l-1]);
						erf -> SetParameters(valuesxmiddle[l-1], 1.0);

						logistic -> SetLineColor(kBlue);
						erf -> SetLineColor(kRed);
				
						//cout << "Message [Fit] fitting logistic ..." << endl;
						TFitResultPtr r_log = graphs[i][j][l]-> Fit(logistic, "s");
						Double_t chi2_log = r_log -> Chi2();
						//graphs[i][j][l]-> DrawClone("AC*"); //if one wants to draw results from logistic sigmoid

						//cout << "Message [Fit] logistic fit status: " << r_log << endl;

						//cout << "Message [transition point] transition point from logistic:" << logistic -> GetX(level_th) << "(" << charge_transf(logistic -> GetX(level_th)) << " in ke)" << endl;
					
						cout << "Message [Fit] fitting erf ..." << endl;
						TFitResultPtr r_erf = graphs[i][j][l]-> Fit(erf, "s");
						Double_t chi2_erf = r_erf -> Chi2();
						graphs[i][j][l]-> Fit(erf);
						graphs[i][j][l]-> Draw("AC*");

						fit_status -> Fill(r_erf);

						cout << "Message [Fit] Erf fit status: " << r_erf << endl;
						cout << "Message [transition point] transition point from Erf:" << erf -> GetX(level_th)  << "(" << charge_transf(erf -> GetX(level_th)) << " in ke)" << endl;
				
						double Merrx1, Merrx2, Lerrx1, Lerrx2, Eerrx1, Eerrx2;
				
						tie(Merrx1, Merrx2) = errorXMV(valuesxmiddle[l-1], valuesxzero[l-1], valuesxone[l-1]);
						//tie(Lerrx1, Lerrx2) = errorXlogistic(logistic, logistic -> GetX(level_th), level_th); //uncomment if transition poit is not 0.5
						//tie(Eerrx1, Eerrx2) = errorXerf(erf, erf -> GetX(level_th), level_th); //uncomment if transition poit is not 0.5

						Lerrx1 = logistic -> GetParError(1); //comment if transition poit is not 0.5
						Eerrx1 = erf -> GetParError(0);  //comment if transition poit is not 0.5

						//cout << "---------------------[INFO] Errors-------------------- " << endl;
						//cout << "Message [error] middle point: " << Merrx1 << " " << Merrx2 << endl;
						//cout << "Message [error] logistic: " << Lerrx1 << " " << Lerrx2 << endl;
						//cout << "Message [error] Erf: " << Eerrx1 << " " << Eerrx2 << endl;

						//cout << "Message [Fit] logistic fit Chi2:" <<  logistic -> GetChisquare() << endl;
				                //cout << "Message [Fit] Erf fit Chi2: " << erf -> GetChisquare() << endl;

						histo_logistic -> Fill(chi2_log);
						histo_erf -> Fill(chi2_erf);

						file << "EC" << i <<" PMT" << j <<" Ch"<< l << setw(30) << valuesxmiddle[l-1] << setw(30) << Merrx1 << setw(30) << Merrx2 << setw(30) << charge_transf(valuesxmiddle[l-1]) << setw(30) << logistic -> GetX(level_th)<< setw(30) << Lerrx1 << setw(30) << Lerrx1 << setw(30) << chi2_log << setw(30) << r_log << setw(30) << charge_transf(logistic -> GetX(level_th)) << setw(30)<< erf -> GetX(level_th) << setw(30) << Eerrx1 << setw(30) << Eerrx1 << setw(30)<< chi2_erf << setw(30)  << r_erf << setw(30) << charge_transf(erf -> GetX(level_th)) << endl;

		//for (auto i = valuesxone.begin(); i != valuesxone.end(); ++i) cout << *i << endl;	
		//for (auto i = valuesxzero.begin(); i != valuesxzero.end(); ++i) cout << *i << endl;	
		//for (auto i = valuesxmiddle.begin(); i != valuesxmiddle.end(); ++i) cout << *i << endl;
		//cout << "Sizes of vectors:" << valuesxone.size() << " " << "=" << " " << valuesxzero.size() << " " << "=" << " " << valuesxmiddle.size() << ".";
				}
			}

		valuesxone.clear();
		valuesxzero.clear();
		valuesxmiddle.clear();
		noisecode.clear();

		c[4*i+j] -> SaveAs(savedir +name_histo+"/"+TString::Format("S-Curve_EC%d_PMT%d.pdf", i, j));
		c[4*i+j] -> Close();

		}


	}
		
	
		file.close();

		//Draw chi2 and fit_status

		/*c1 -> Divide(1,2);

		c1 -> cd(1);
		histo_logistic -> GetXaxis() -> SetTitle("#chi^{2}");
		histo_logistic -> GetYaxis() -> SetTitle("counts");
		histo_logistic -> Draw();
		TLatex t1(0.8, 400,"#frac{1}{1+exp(-a(x-x_{m}))}");
		t1.Draw();*/

		TCanvas *c1 = new TCanvas();
		c1 -> cd();
		histo_erf -> GetXaxis() -> SetTitle("#chi^{2}");
		histo_erf -> GetYaxis() -> SetTitle("counts");
		histo_erf -> Draw();
		TLatex t2(0.8, 200,"0.5(1+Erf[a(x-x_{m})])");
		t2.Draw();

		c1 -> SaveAs(savedir +name_histo+TString("/fit_properties.pdf"));
		c1 -> Close();

		TCanvas *c2 = new TCanvas();
		c2 -> cd();
		fit_status -> Draw();
		c2 -> SaveAs(savedir +name_histo+TString("/fit_status.pdf"));
		c2 -> Close();	

		}
	
	}	



};
