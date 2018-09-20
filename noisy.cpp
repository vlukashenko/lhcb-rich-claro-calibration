#include <iostream>
#include <fstream>
#include <vector>
#include <cstring>
#include <stdio.h>
#include <iomanip>
#include <regex>

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

//readout values of transition points
vector <string> values(string name_file, int histo_num, int opt, int header = 6){//opt==1 middle point calculation; opt==2 logistic sigmoid; opt==3 Erf fit;
//Note: one can choose average/max values of chi2 as threshold, by default coefficient is 2.5 chi^2_{max or av}. 
	string board;
	if(opt == 1) board = "Board1_13";
	if(opt == 0) board = "Board1_06";
	int i = 0;		
	int k = 0;

	ifstream file;
	file.open(name_file, ios::in);

	vector <string> names;
	vector <double> chi2_log;
	vector <double> chi2_erf;

	TH1D *histo_chi2 = new TH1D("histo_chi2", "histo_chi2", 1000, 0, 1000);

	double aver_chi2_log;
	double aver_chi2_erf;

	if (file.is_open()){

		while(i<header){
			file.ignore(256, '\n');
			i++;
			}
		string line;

		while (getline(file, line)) {
    			istringstream is(line);    // construct temporary istringstream
			//cout << line << endl;
    			double col2, col3, col4, col5, col6, col7, col8, col9, col10, col11, col12, col13, col14, col15, col16, col17;
    			string cols0, cols1, cols2;
			
			//col9 -chi2 logistic, col15 - chi2 erf
			if (is >> cols0  >> cols1 >> cols2 >> col2 >> col3 >> col4 >> col5 >> col6 >> col7 >> col8 >> col9 >> col10 >> col11 >> col12 >> col13>> col14 >> col15 >> col16 >> col17) {
				
   			chi2_log.push_back(col9);
			chi2_erf.push_back(col15);
	
			}
    			else cout << "This line didn't meet the expected format." << endl;
		}
		
		for(auto i = chi2_log.begin(); i != chi2_log.end(); ++i) aver_chi2_log += *i/chi2_log.size();	
		for(auto i = chi2_erf.begin(); i != chi2_erf.end(); ++i) {aver_chi2_erf += *i/chi2_erf.size(); histo_chi2 -> Fill(*i);}


		TCanvas *c = new TCanvas();
		c-> cd();
		histo_chi2 -> Draw();
  		c->Update();
 		double chimax = histo_chi2 -> GetMaximumBin();
  		double maxY = 10;
		//aver_chi2_erf = chimax;
		
		TLine l1(chimax, 1e-8, chimax, maxY); 
   		TLatex latex;  
		latex.SetTextSize(0.05);
		//latex.DrawLatex(chimax+1, 8,"#chi^{2}_{max}");
   		//latex.DrawLatex(2.5*chimax+1, 8,"2.5*#chi^{2}_{max}");


		TLine l2(2.5*chimax,1e-8, 2.5*chimax, maxY);

		l1.SetLineStyle(9);
		l1.SetLineColor(2); 
		l1.SetLineWidth(3);
		l2.SetLineStyle(9);
		l2.SetLineColor(2); 
		l2.SetLineWidth(3);
		l1.Draw("same");
		l2.Draw("same");

		latex.DrawLatex(chimax+5, 10,"#chi^{2}_{max}");
   		latex.DrawLatex(2.5*chimax+5, 10,"2.5*#chi^{2}_{max}");

		c -> SaveAs("/afs/cern.ch/user/v/valukash/RICH-threshold-analysis/DAC/"+TString(board)+"/chi2_run"+TString(to_string(histo_num))+".pdf");

	}
 	
	file.close();

	file.open(name_file, ios::in);
	i = 0;
	if (file.is_open()){
		int sum = 0;
		while(i<header){
			file.ignore(256, '\n');
			i++;
			}
		string line;

		while (getline(file, line)) {
    			istringstream is(line);    // construct temporary istringstream
    			double col2, col3, col4, col5, col6, col7, col8, col9, col10, col11, col12, col13, col14, col15, col16, col17;
    			string cols0, cols1, cols2;

			//col9 -chi2 logistic, col15 - chi2 erf
			if (is >> cols0  >> cols1 >> cols2 >> col2 >> col3 >> col4 >> col5 >> col6 >> col7 >> col8 >> col9 >> col10 >> col11 >> col12 >> col13>> col14 >> col15 >> col16 >> col17) {
   			if(col9>2.5*aver_chi2_log || col15 > 2.5*aver_chi2_erf) {names.push_back(cols0+" "+cols1+" "+cols2); sum+=1;}

			//cout << cols0 << " " << cols1 << " " << cols2 << " " << col2 << " " << col3 << " " << col4 << " " << col5 << " " <<  col6 << " " << col7 << " " << col8 <<  " " << col9 <<  " " << col10  << " " << col11 << col12 << " " << col13 << " " << col14 << " " << col15 << " " << col16 << " " << col17 << endl;
  				
			}
    			else cout << "Message [file txt]: this line didn't meet the expected format" << endl;
		}
		cout << "Message [noisy channels] total number: " << sum << endl;

	}
 
	file.close();

	return names;	

}

//readout file header
void header(string name_file, int header = 4) //read values from file header; int nheadlines - number of header lines
{
	ifstream file;
	file.open(name_file, ios::in);
	string line;
	int i = 1;
	if (file.is_open()){
		file.ignore(256, '\n');
		while(i<header){
			getline(file, line);
			cout << line << endl;
			i++;
			}

 	}
	file.close();


}

void noisy(){
	//DATA
	int nruns = 15;
	string dir1 = "/afs/cern.ch/user/v/valukash/RICH-threshold-analysis/DAC/Board1_13/S_curves/";
	string dir2 = "/afs/cern.ch/user/v/valukash/RICH-threshold-analysis/DAC/Board1_06/S_curves/";
 	int num[] = {114, 142, 146, 63, 59, 61, 64, 93, 96, 98, 100, 102, 273, 272, 271};

	vector <string> histo_names;

	cout << "======================================= BOARD 13 =======================================" << endl;
	for(int k =0; k < nruns; k++) {
		vector <string> noised;
		vector <int> val;
		histo_names.push_back("histo_run"+to_string(num[k]));
		cout << "____________________________________________" << endl;
		header(dir1+histo_names[k]+"/thresholds.txt");
		noised = values(dir1+histo_names[k]+"/thresholds.txt", num[k], 1);
		for(auto i = noised.begin(); i != noised.end(); ++i) cout << *i << endl;

	}
	cout << "======================================= BOARD 06 =======================================" << endl;

	for(int k=0 ; k < nruns; k++) {
		vector <string> noised;
		vector <int> val;
		histo_names.push_back("histo_run"+to_string(num[k]));
		cout << "____________________________________________" << endl;
		header(dir1+histo_names[k]+"/thresholds.txt");
		noised = values(dir2+histo_names[k]+"/thresholds.txt", num[k], 0);
		for(auto i = noised.begin(); i != noised.end(); ++i) cout << *i << endl;

	}

};