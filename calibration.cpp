#include <iostream>
#include <fstream>
#include <vector>
#include <cstring>
#include <stdio.h>
#include <algorithm>
#include <functional>
#include <cstdlib>
#include <cmath>


#include <TCanvas.h>
#include <TMultiGraph.h>
#include <TMath.h>
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TStyle.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TPad.h>
#include <TString.h>
#include <TF1.h>
#include <TGraphAsymmErrors.h>
#include <TH1F.h>
#include <TLatex.h>

using namespace std;

struct value{ //structure for transition point + errors.

	vector <string> name;

	vector <double> arrDAC; 
	vector <double> arrCh;

	vector <double> err1DAC;
	vector <double> err2DAC;

    	void clear(){
       		name.clear();    //name of channel
		arrDAC.clear();  //transition point in DAC counts
		arrCh.clear();   //transition point in charge
		err1DAC.clear(); //transition point error down in DAC counts
		err2DAC.clear(); //transition point error up in DAC counts
   	 }
};


vector <string> header(string name_file, int nheadlines = 4) //read values from file header; int nheadlines - number of header lines
{
	ifstream file;
	file.open(name_file, ios::in);
	string line;
 	vector <string> result;
	int i = 1;
	if (file.is_open()){
		file.ignore(256, '\n');
		while(i<nheadlines){
			getline(file, line);
			result.push_back(line);
			i++;
			}
 	}
	file.close();
	return result;
}

value values(string name_file, int opt = 3, int nheadlines = 6){ //get transition points; opt: 1 - middle point calculation; 2 - logistic signoid; 3 - Erf fit

	int i = 0;		
	int k = 0;

	ifstream file;
	file.open(name_file, ios::in);

	value results;

	if (file.is_open()){

		while(i<nheadlines){
			file.ignore(256, '\n');
			i ++;
			}
		string line;

		while (getline(file, line)) {
    			istringstream is(line);    // construct temporary istringstream
    			double col2, col3, col4, col5, col6, col7, col8, col9, col10, col11, col12, col13, col14, col15, col16, col17;
    			string cols0, cols1, cols2;
   
			if (is >> cols0  >> cols1 >> cols2 >> col2 >> col3 >> col4 >> col5 >> col6 >> col7 >> col8 >> col9 >> col10 >> col11 >> col12 >> col13>> col14 >> col15 >> col16 >> col17) {
				if(opt==1){                                       // MV
					results.name.push_back(cols0+"_"+cols1+"_"+cols2);
					results.arrDAC.push_back(col2*1.0);
					results.arrCh.push_back(col5*1.0);
					results.err1DAC.push_back(col3*1.0);
					results.err2DAC.push_back(col4*1.0);
					
				}
				if(opt==2){ 					  // logistic
					results.name.push_back(cols0 +"_"+cols1+"_"+cols2);
					results.arrDAC.push_back(col6);
					results.arrCh.push_back(col11);
					results.err1DAC.push_back(col7);
					results.err2DAC.push_back(col8);
				
				}

				if(opt==3){                                       // erf
					results.name.push_back(cols0+"_"+cols1+"_"+cols2);
					results.arrDAC.push_back(col12);
					results.arrCh.push_back(col17);
					results.err1DAC.push_back(col13);
					results.err2DAC.push_back(col14);
				
				}
				++k;
			}
    			else cout << "Message [txt file] this line didn't meet the expected format" << endl;
		}
 	

	}
	file.close();
	return results;	
}


vector <double> change_one(string name_file, string name_PMT, int nheadlines = 6){ //change one value if middle point calculation should be used

	int i = 0;		
	int k = 0;

	ifstream file;
	file.open(name_file, ios::in);

	vector <double> results;

	if (file.is_open()){

		while(i<nheadlines){
			file.ignore(256, '\n');
			i ++;
			}
		string line;

		while (getline(file, line)) {
    			istringstream is(line);    // construct temporary istringstream
    			double col2, col3, col4, col5, col6, col7, col8, col9, col10, col11, col12, col13, col14, col15, col16, col17;
    			string cols0, cols1, cols2;
   
			if (is >> cols0  >> cols1 >> cols2 >> col2 >> col3 >> col4 >> col5 >> col6 >> col7 >> col8 >> col9 >> col10 >> col11 >> col12 >> col13>> col14 >> col15 >> col16 >> col17) {
				                           
				if((cols0 +"_"+cols1 +"_"+cols2) == name_PMT) {results.push_back(col5); results.push_back(col3); results.push_back(col4);}
				++k;
				//cout << cols0 << " " << cols1 << " " << cols2 << " " << col2 << " " << col3 << " " << col4 << " " << col5 << " " <<  col6 << " " << col7 << " " << col8 <<  " " << col9 <<  " " << col10  << " " << col11 << col12 << " " << col13 << " " << col14 << " " << col15 << " " << col16 << " " << col17 <<endl;
			}
    			else cout << "Message [txt file] this line didn't meet the expected format." << endl;
    
		}
 	

	}
	file.close();
	return results;	
}


double av_chi2(string name_file, int opt = 3, int nheadlines = 6){ //average chi2 value; opt similar as before

	int i = 0;		
	int k = 0;

	ifstream file;
	file.open(name_file, ios::in);

	vector <string> names;
	vector <double> chi2_log;
	vector <double> chi2_erf;

	double aver_chi2_log;
	double aver_chi2_erf;

	if (file.is_open()){

		while(i<nheadlines){
			file.ignore(256, '\n');
			i++;
			}
		string line;

		while (getline(file, line)) {
    			istringstream is(line);    // construct temporary istringstream
    			double col2, col3, col4, col5, col6, col7, col8, col9, col10, col11, col12, col13, col14, col15, col16, col17;
    			string cols0, cols1, cols2;
			if (is >> cols0  >> cols1 >> cols2 >> col2 >> col3 >> col4 >> col5 >> col6 >> col7 >> col8 >> col9 >> col10 >> col11 >> col12 >> col13>> col14 >> col15 >> col16 >> col17) {
				
   				chi2_log.push_back(col9);
				chi2_erf.push_back(col15);
				//cout << cols0 << " " << cols1 << " " << cols2 << " " << col2 << " " << col3 << " " << col4 << " " << col5 << " " <<  col6 << " " << col7 << " " << col8 <<  " " << col9 <<  " " << col10  << " " << col11 << col12 << " " << col13 << " " << col14 << " " << col15 << " " << col16 << " " << col17 << endl;
				}
    			else cout << "Message [txt file] this line didn't meet the expected format." << endl;
		}
		
		for(auto i = chi2_log.begin(); i != chi2_log.end(); ++i) aver_chi2_log += *i/chi2_log.size();	
		for(auto i = chi2_erf.begin(); i != chi2_erf.end(); ++i) aver_chi2_erf += *i/chi2_erf.size();	

	}
 	
	file.close();
        if(opt == 1) {cout << "Message [c++ function] wrong opt in av_chi2; middle point calculation is not a fit" << endl; return 1;}
	if(opt == 2) return aver_chi2_log;
	if(opt == 3) return aver_chi2_erf;
};

//functions to shift vectors to the left/right
vector <int> shift_vector_int_R(vector <int> vec){
	vector<int> new_(vec);
	for (size_t i = 1; i < vec.size(); i++) vec.at(i) = new_.at(i - 1);
	vec.push_back(new_.at(vec.size() - 1));
	return vec;
	}

vector <int> shift_vector_int_L(vector <int> vec){
	vector<int> new_(vec);
	for (size_t i = 1; i < vec.size()-1; i++) vec.at(i) = new_.at(i + 1);
	vec.push_back(new_.at(0));
	return vec;
	}

vector <double> shift_vector_R(vector <double> vec){
			
	vector<double> new_(vec);
	for (size_t i = 1; i < vec.size(); i++) vec.at(i) = new_.at(i - 1);
	vec.push_back(new_.at(vec.size() - 1));
	return vec;
	}

vector <double> shift_vector_L(vector <double> vec){
			
	vector<double> new_(vec);
	for (size_t i = 1; i < vec.size()-1; i++) vec.at(i) = new_.at(i + 1);
	vec.push_back(new_.at(0));
	return vec;
	}

vector <string> shift_vector_name_R(vector <string> vec){
			
	vector<string> new_(vec);
	for (size_t i = 1; i < vec.size(); i++) vec.at(i) = new_.at(i - 1);
	vec.push_back(new_.at(vec.size() - 1));
	return vec;
	}

vector <string> shift_vector_name_L(vector <string> vec){
			
	vector<string> new_(vec);
	for (size_t i = 1; i < vec.size()-1; i++) vec.at(i) = new_.at(i + 1);
	vec.push_back(new_.at(0));
	return vec;
	}

//linear fit
Double_t fitF(Double_t *x, Double_t *par) {return par[0]+par[1]*x[0];} 	

bool decide(string name_file, string name_PMT, double chi2, int opt = 3, int nheadlines = 6){ //decide is noisy  opt similar as before

	ifstream file;
	file.open(name_file, ios::in); 

	int i = 0;
	bool decision;
	if (file.is_open()){

		while(i<nheadlines){
			file.ignore(256, '\n');
			++i;
			}
		string line;

		while (getline(file, line)) {
    			istringstream is(line);    // construct temporary istringstream
    			double col2, col3, col4, col5, col6, col7, col8, col9, col10, col11, col12, col13, col14, col15, col16, col17;
    			string cols0, cols1, cols2;

			if (is >> cols0  >> cols1 >> cols2 >> col2 >> col3 >> col4 >> col5 >> col6 >> col7 >> col8 >> col9 >> col10 >> col11 >> col12 >> col13>> col14 >> col15 >> col16 >> col17) {
				
   			if((cols0 +"_"+cols1 +"_"+cols2) == name_PMT){
				if(opt == 2 && col9>2.5*chi2)  decision = true;
				if(opt == 3 && col15 > 2.5*chi2)  decision = true;
				else decision = false;
				}
			//cout << cols0 << " " << cols1 << " " << cols2 << " " << col2 << " " << col3 << " " << col4 << " " << col5 << " " <<  col6 << " " << col7 << " " << col8 <<  " " << col9 <<  " " << col10  << " " << col11 << col12 << " " << col13 << " " << col14 << " " << col15 << " " << col16 << " " << col17 << endl;
  				
			}
    			else cout << "Message [txt file] this line didn't meet the expected format." << endl;
		}

	}
 
	file.close();
	return decision;

}



int calibration(){

	vector <int> EC1s, EC2s, EC3s;	
	vector <int> PMT1s, PMT2s, PMT3s;
	vector <int> ch1s, ch2s, ch3s;
	vector <string> bad_calibr_names;
	string dir;

	string name_histo1, name_histo2, name_histo3;

	double slope_13[5][4][2][64];
	double offset_13[5][4][2][64];
	double slope_06[5][4][2][64];
	double offset_06[5][4][2][64];

	int count_MV = 0;
	int count_HE_sl = 0;
	int count_HE_off = 0;
	double max_error = 0.1;
	double coefficient = 15.67;

	int nentries_off;
	int nentries_sl;

	double av_slope_all = 0;
	double av_at_zero_all = 0;
	double av_slope_18 = 0;
	double av_at_zero_18 = 0;
	double bad_calibr = 0;

	int option = 3;        // 1 - MV, 2 - logistic, 3 - erf;
	double scale = 1.;
	int up = 10000;        //att3: 10, att2: 5, att1: 3, att0: 1 if scale = 1000 (Me-)
	int up_slope = 450;    //0.45 if scale = 1000 (Me-)

	TString name_graph;
 	string name_canvas;
	TGraphAsymmErrors *calibration = new TGraphAsymmErrors();

	TFile calibration_file("calibration_new.root", "RECREATE");

  	TTree *calibration06 = new TTree("calibration06", "calibration06");
  	TTree *calibration13 = new TTree("calibration13", "calibration13");

	TH1D * histo_count_MV = new TH1D("histo_count_MV", "histo_count", 100, 0, 100);                 //count of noisy channels
	TH1D * histo_count_HE_off = new TH1D("histo_count_HE_off", "histo_count_HE_off", 100, 0, 100);  //how many more than 10% error
	TH1D * histo_count_HE_sl = new TH1D("histo_count_HE_sl", "histo_count_HE_sl", 100, 0, 100);     //how many more than 10% error

   	Int_t EC06_v, PMT06_v, CH06_v;
   	Double_t thresholdstep_06_att0_v, thresholdstep_06_att1_v, thresholdstep_06_att2_v, thresholdstep_06_att3_v, thresholdstep_06_off0_v;
	Double_t offset_06_att0_v, offset_06_att1_v, offset_06_att2_v, offset_06_att3_v, offset_06_off0_v;

   	Int_t EC13_v, PMT13_v, CH13_v;
   	Double_t thresholdstep_13_att0_v, thresholdstep_13_att1_v, thresholdstep_13_att2_v, thresholdstep_13_att3_v, thresholdstep_13_off0_v;
	Double_t offset_13_att0_v, offset_13_att1_v, offset_13_att2_v, offset_13_att3_v, offset_13_off0_v;

	int ECb, PMTb, CHb;

	Int_t sizeb = 64000;
	Int_t split = 0;

   	calibration06->Branch("EC06", &EC06_v, sizeb, split);
   	calibration06->Branch("PMT06", &PMT06_v, sizeb, split);
   	calibration06->Branch("CH06", &CH06_v, sizeb, split);

  	calibration06->Branch("thresholdstep_06_att0", &thresholdstep_06_att0_v, sizeb, split);
    	calibration06->Branch("thresholdstep_06_att1", &thresholdstep_06_att1_v, sizeb, split);
  	calibration06->Branch("thresholdstep_06_att2", &thresholdstep_06_att2_v, sizeb, split);
  	calibration06->Branch("thresholdstep_06_att3", &thresholdstep_06_att3_v, sizeb, split); 	
  	calibration06->Branch("thresholdstep_06_off0", &thresholdstep_06_off0_v, sizeb, split);
	
  	calibration06->Branch("offset_06_att0", &offset_06_att0_v, sizeb, split);
  	calibration06->Branch("offset_06_att1", &offset_06_att1_v, sizeb, split);
  	calibration06->Branch("offset_06_att2", &offset_06_att2_v, sizeb, split);
  	calibration06->Branch("offset_06_att3", &offset_06_att3_v, sizeb, split);
  	calibration06->Branch("offset_06_off0", &offset_06_off0_v, sizeb, split);

   	calibration13->Branch("EC13", &EC13_v, sizeb, split);
   	calibration13->Branch("PMT13", &PMT13_v, sizeb, split);
   	calibration13->Branch("CH13", &CH13_v, sizeb, split);

  	calibration13->Branch("thresholdstep_13_att0", &thresholdstep_13_att0_v, sizeb, split);
    	calibration13->Branch("thresholdstep_13_att1", &thresholdstep_13_att1_v, sizeb, split);
  	calibration13->Branch("thresholdstep_13_att2", &thresholdstep_13_att2_v, sizeb, split);
  	calibration13->Branch("thresholdstep_13_att3", &thresholdstep_13_att3_v, sizeb, split);
  	calibration13->Branch("thresholdstep_13_off0", &thresholdstep_13_off0_v, sizeb, split); 	
	
  	calibration13->Branch("offset_13_att0", &offset_13_att0_v, sizeb, split);
 	calibration13->Branch("offset_13_att1", &offset_13_att1_v, sizeb, split);
  	calibration13->Branch("offset_13_att2", &offset_13_att2_v, sizeb, split);
  	calibration13->Branch("offset_13_att3", &offset_13_att3_v, sizeb, split);
  	calibration13->Branch("offset_13_off0", &offset_13_off0_v, sizeb, split);

for(int g=0; g<2; g++){

	double offset;
	double slope;

	for(int q=0;q<5;q++){

		gStyle->SetOptStat(0);
		TH1D *histo_slope_all_06 = new TH1D("histo_slope_all_06", "histo_slope_all", 225, 0, up_slope);  
		TH1D *histo_at_zero_all_06 = new TH1D("histo_at_zero_all_06", "histo_at_zero_all", 500, 0, up); 
		histo_slope_all_06 -> GetYaxis() -> SetTitle("counts");
		histo_at_zero_all_06 -> GetYaxis() -> SetTitle("counts");
		histo_slope_all_06 -> GetXaxis() -> SetTitle("#Delta Th, ke^{-}");
		histo_at_zero_all_06 -> GetXaxis() -> SetTitle("|Th(0)|, ke^{-}");

		TH1D *histo_slope_error_06 = new TH1D("histo_slope_error_06", "histo_slope_error", 3600, 0, up_slope);
		TH1D *histo_at_zero_error_06 = new TH1D("histo_at_zero_error_06", "histo_at_zero_error", 10000, 0, up);
		histo_slope_error_06 -> GetYaxis() -> SetTitle("counts");
		histo_at_zero_error_06 -> GetYaxis() -> SetTitle("counts");
		histo_slope_error_06 -> GetXaxis() -> SetTitle("#Delta(#Delta Th), ke^{-}");
		histo_at_zero_error_06 -> GetXaxis() -> SetTitle("#Delta|Th(0)|, ke^{-}");

		TH1D *histo_slope_all_13 = new TH1D("histo_slope_all_13", "histo_slope_all", 225, 0, up_slope); 
		TH1D *histo_at_zero_all_13 = new TH1D("histo_at_zero_all_13", "histo_at_zero_all", 500, 0, up);
		histo_slope_all_13 -> GetYaxis() -> SetTitle("counts");
		histo_at_zero_all_13 -> GetYaxis() -> SetTitle("counts");
		histo_slope_all_13 -> GetXaxis() -> SetTitle("#Delta Th, ke^{-}");
		histo_at_zero_all_13 -> GetXaxis() -> SetTitle("|Th(0)|, ke^{-}");

		TH1D *histo_slope_error_13 = new TH1D("histo_slope_error_13", "histo_slope_error", 3600, 0, up_slope); 
		TH1D *histo_at_zero_error_13 = new TH1D("histo_at_zero_error_13", "histo_at_zero_error", 10000, 0, up);
		histo_slope_error_13 -> GetYaxis() -> SetTitle("counts");
		histo_at_zero_error_13 -> GetYaxis() -> SetTitle("counts");
		histo_slope_error_13 -> GetXaxis() -> SetTitle("#Delta(#Delta Th), ke^{-}");
		histo_at_zero_error_13 -> GetXaxis() -> SetTitle("#Delta|Th(0)|, ke^{-}");
	
		//==============DIRECTOR=========================
		if(g==0){
			cout << "=======================================BOARD 13================================================" << endl;
			dir = "/afs/cern.ch/user/v/valukash/RICH-threshold-analysis/DAC/Board1_13/S_curves/";
		}
	
		if(g==1){
			cout << "=======================================BOARD 06================================================" << endl;
			dir = "/afs/cern.ch/user/v/valukash/RICH-threshold-analysis/DAC/Board1_06/S_curves/";
		}

		//=================DATA==========================
		//Attenuation: 3; Offset: 1
		if(q==0){
			cout<< "==================================Attenuation: 3; Offset: 1====================================" << endl;
			name_histo1 ="histo_run271";   //42 (271)
 			name_histo2 ="histo_run272";   //38
 			name_histo3 ="histo_run273";   //36
		}

		//Attenuation: 2; Offset: 1
		if(q==1){
			cout<< "==================================Attenuation: 2; Offset: 1====================================" << endl;
 			name_histo1 ="histo_run102";   //62 
			name_histo2 ="histo_run100";   //52
 			name_histo3 ="histo_run98";    //42
		}
 	
		//Attenuation: 1; Offset: 1
		if(q==2){
			cout<< "==================================Attenuation: 1; Offset: 1====================================" << endl;
			name_histo1 ="histo_run96";    //62
			name_histo2 ="histo_run93";    //52
 			name_histo3 ="histo_run64";    //42
		}

		//Attenuation: 0; Offset: 1
		if(q==3){
			cout<< "==================================Attenuation: 0; Offset: 1====================================" << endl;
			name_histo1 ="histo_run61";    //62
			name_histo2 ="histo_run59";    //52
 			name_histo3 ="histo_run63";    //42
		}

		//Attenuation: 0; Offset: 0	
		if(q==4){
			cout<< "==================================Attenuation: 0; Offset: 0====================================" << endl;
 			name_histo1 ="histo_run146";    //62
			name_histo2 ="histo_run142";    //52
 			name_histo3 ="histo_run114";    //42
		}

		if ( calibration_file.IsOpen() ) cout << "Message [ROOT file state]: file is open" << endl;

  		ofstream file_par;
    		file_par.open(dir +name_histo1+"/"+TString::Format("file_par.txt"), ios::out);
   		file_par << "Name" << setw(30) << "Attenuation" << setw(30) << "Offset" << setw(30) << "slope, ke-" << setw(30) << "|at zero|, ke-" << endl;   
		
		//read out transition point values for 3 threshold values;
		value thr1 = values(dir+name_histo1+"/thresholds.txt", option); 
		value thr2 = values(dir+name_histo2+"/thresholds.txt", option);
		value thr3 = values(dir+name_histo3+"/thresholds.txt", option);

		//get header of a "thresholds.txt" file
		vector <string> head = header(dir+name_histo1+"/thresholds.txt");
		size_t pos1 = head[0].find(":");
		size_t pos2 = head[1].find(":");
		int att = atoi(TString(head[0].substr(pos1+1, pos1+2)));
		int off = atoi(TString(head[1].substr(pos2+1, pos2+2)));

		//calculate averaged chi2
		double average_chi2_1 = av_chi2(dir+name_histo1+"/thresholds.txt");
		double average_chi2_2 = av_chi2(dir+name_histo2+"/thresholds.txt");
		double average_chi2_3 = av_chi2(dir+name_histo3+"/thresholds.txt");

		//make channel name lists same size
		auto max0 = max(thr1.name.size(), thr2.name.size());
		auto max_value  = max(max0, thr3.name.size());
		if(thr1.name.size()<max_value) thr1.name.insert(thr1.name.end(), max_value-thr1.name.size(), "EC0_PMT0_ch0"); 
		if(thr2.name.size()<max_value) thr2.name.insert(thr2.name.end(), max_value-thr2.name.size(), "EC0_PMT0_ch0"); 
		if(thr3.name.size()<max_value) thr3.name.insert(thr3.name.end(), max_value-thr3.name.size(), "EC0_PMT0_ch0");  
		
		cout << "========================[INFO]========================" << endl;

		for(int j = 0; j < max_value; j++){
			
			//decide if channel is noisy
			bool decide1 = decide(dir+name_histo1+"/thresholds.txt", thr1.name[j], average_chi2_1);
			bool decide2 = decide(dir+name_histo2+"/thresholds.txt", thr2.name[j], average_chi2_2);
			bool decide3 = decide(dir+name_histo3+"/thresholds.txt", thr3.name[j], average_chi2_3);
			
			//make arrays with EC, PMT, CH numbers
			size_t posEC1 = thr1.name[j].find("EC");
			size_t posEC2 = thr2.name[j].find("EC");
			size_t posEC3 = thr3.name[j].find("EC");

			EC1s.push_back(atoi(TString(thr1.name[j].substr(posEC1+2, posEC1+3))));
			EC2s.push_back(atoi(TString(thr2.name[j].substr(posEC2+2, posEC2+3))));
			EC3s.push_back(atoi(TString(thr3.name[j].substr(posEC3+2, posEC3+3))));

			size_t posPMT1 = thr1.name[j].find("PMT");
			size_t posPMT2 = thr2.name[j].find("PMT");
			size_t posPMT3 = thr3.name[j].find("PMT");

			PMT1s.push_back(atoi(TString(thr1.name[j].substr(posPMT1+3, posPMT1+5))));
			PMT2s.push_back(atoi(TString(thr2.name[j].substr(posPMT2+3, posPMT2+5))));
			PMT3s.push_back(atoi(TString(thr3.name[j].substr(posPMT3+3, posPMT3+5))));

			size_t pos1 = thr1.name[j].find("Ch");
			size_t pos2 = thr2.name[j].find("Ch");
			size_t pos3 = thr3.name[j].find("Ch");

			ch1s.push_back(atoi(TString(thr1.name[j].substr(pos1+2))));
			ch2s.push_back(atoi(TString(thr2.name[j].substr(pos2+2))));
			ch3s.push_back(atoi(TString(thr3.name[j].substr(pos3+2))));

			int jpos1 = find(ch1s.begin(), ch1s.end(), atoi(TString(thr1.name[j].substr(pos1+2)))) - ch1s.begin();
			int jpos2 = find(ch2s.begin(), ch2s.end(), atoi(TString(thr2.name[j].substr(pos2+2)))) - ch2s.begin();
			int jpos3 = find(ch3s.begin(), ch3s.end(), atoi(TString(thr3.name[j].substr(pos3+2)))) - ch3s.begin();

			//if(decide1==1) {EC1s.erase(EC1s.begin()+jpos1); PMT1s.erase(PMT1s.begin()+jpos1); ch1s.erase(ch1s.begin()+jpos1); thr1.name.push_back("EC0_PMT0_ch0");}
			//if(decide2==1) {EC2s.erase(EC2s.begin()+jpos2); PMT2s.erase(PMT2s.begin()+jpos2); ch2s.erase(ch2s.begin()+jpos2); thr2.name.push_back("EC0_PMT0_ch0");}
			//if(decide3==1) {EC3s.erase(EC3s.begin()+jpos3); PMT3s.erase(PMT3s.begin()+jpos3); ch3s.erase(ch3s.begin()+jpos3); thr3.name.push_back("EC0_PMT0_ch0");}

			vector <double> buf;
		
			if(g!=1 && q!=2){ //hack for run96 : not zero initial counts
				if(decide1==1) {count_MV++; buf = change_one(dir+name_histo1+"/thresholds.txt",  thr1.name[j]); thr1.arrCh[j] = buf[0]; thr1.err1DAC[j] = buf[1]; thr1.err2DAC[j] = buf[2];  cout << "Message [noise level]: noisy channel " << thr1.name[j] << " in " << name_histo1 << endl; buf.erase(buf.begin(), buf.end());}
				if(decide2==1) {count_MV++; buf = change_one(dir+name_histo2+"/thresholds.txt",  thr2.name[j]); thr2.arrCh[j] = buf[0]; thr2.err1DAC[j] = buf[1]; thr2.err2DAC[j] = buf[2];  cout << "Message [noise level]: noisy channel " << thr2.name[j] << " in " << name_histo2 << endl; buf.erase(buf.begin(), buf.end());}		
				if(decide3==1) {count_MV++; buf = change_one(dir+name_histo3+"/thresholds.txt",  thr3.name[j]); thr3.arrCh[j] = buf[0]; thr3.err1DAC[j] = buf[1]; thr3.err2DAC[j] = buf[2];  cout << "Message [noise level]: noisy channel " << thr3.name[j] << " in " << name_histo3 << endl; buf.erase(buf.begin(), buf.end());}		
				}
			buf.clear();
		}

		histo_count_MV -> Fill(count_MV); 

		for(int j = 0; j < max_value; j++){

			//cout << "In the general loop, names of graphs: " << thr1.name[j] << " " << thr2.name[j] << " " << thr3.name[j] << endl;
			//cout << "In the general loop: " << "EC" << EC1s[j] << " " << "PMT" << PMT1s[j] << " " << "Ch" << ch1s[j] << endl;
			//cout << "In the general loop: " << "EC" << EC2s[j] << " " << "PMT" << PMT2s[j] << " " << "Ch" << ch2s[j] << endl;
			//cout << "In the general loop: " << "EC" << EC3s[j] << " " << "PMT" << PMT3s[j] << " " << "Ch" << ch3s[j] << endl;

			//hack for att_2 B.13
		if(g==0 && q == 1 && EC3s[j]==3 && PMT3s[j]==0 && ch3s[j]==62) {

			vector <string> buffer_name;
			vector <int> buffer_EC, buffer_PMT, buffer_CH;
			vector <double> buffer1, buffer2, buffer3, buffer4;
		
			buffer_name = shift_vector_name_R(thr2.name); 
			thr2.name.swap(buffer_name); 

			buffer_EC = shift_vector_int_R(EC2s); 
			EC2s.swap(buffer_EC); 

			buffer_PMT = shift_vector_int_R(PMT2s); 
			PMT2s.swap(buffer_PMT); 

			buffer_CH = shift_vector_int_R(ch2s); 
			ch2s.swap(buffer_CH); 

  			buffer1 = shift_vector_R(thr2.arrCh);
			thr2.arrCh.swap(buffer1);

  			buffer2 = shift_vector_R(thr2.arrDAC);
			thr2.arrDAC.swap(buffer2);

  			buffer3 = shift_vector_R(thr2.err1DAC);
			thr2.err1DAC.swap(buffer3);

  			buffer4 = shift_vector_R(thr2.err2DAC);
			thr2.err2DAC.swap(buffer4);

			buffer_name.erase(buffer_name.begin(), buffer_name.end());
			buffer_EC.erase(buffer_EC.begin(), buffer_EC.end());
			buffer_PMT.erase(buffer_PMT.begin(), buffer_PMT.end());
			buffer_CH.erase(buffer_CH.begin(), buffer_CH.end());
			buffer1.erase(buffer1.begin(), buffer1.end());
			buffer2.erase(buffer1.begin(), buffer1.end());
			buffer3.erase(buffer1.begin(), buffer1.end());
			buffer4.erase(buffer1.begin(), buffer1.end());
	
			buffer_name = shift_vector_name_R(thr1.name); 
			thr1.name.swap(buffer_name); 

			buffer_EC = shift_vector_int_R(EC1s); 
			EC1s.swap(buffer_EC); 

			buffer_PMT = shift_vector_int_R(PMT1s); 
			PMT1s.swap(buffer_PMT); 

			buffer_CH = shift_vector_int_R(ch1s); 
			ch1s.swap(buffer_CH); 

  			buffer1 = shift_vector_R(thr1.arrCh);
			thr1.arrCh.swap(buffer1);

  			buffer2 = shift_vector_R(thr1.arrDAC);
			thr1.arrDAC.swap(buffer2);

  			buffer3 = shift_vector_R(thr1.err1DAC);
			thr1.err1DAC.swap(buffer3);

  			buffer4 = shift_vector_R(thr1.err2DAC);
			thr1.err2DAC.swap(buffer4);
			continue;} 
			
			
		//two points only
		if((EC1s[j]>EC2s[j] && EC2s[j]==EC3s[j] &&  PMT1s[j]<PMT2s[j] && PMT2s[j]==PMT3s[j] && ch2s[j]==ch3s[j] && ch1s[j]!=0)||( EC1s[j]==EC2s[j] && EC2s[j]==EC3s[j] && PMT1s[j]>PMT2s[j] && PMT2s[j]==PMT3s[j] && (PMT1s[j-1]!=PMT2s[j]  || ch3s[j-1]==64 || ch2s[j-1]==63) && ch2s[j]==ch3s[j] && ch1s[j]!=0)||(EC1s[j]==EC2s[j] && EC2s[j]==EC3s[j] && PMT1s[j]==PMT2s[j] && PMT2s[j]==PMT3s[j] && ch2s[j]==ch3s[j] && ch1s[j]!=0 && ch1s[j]>ch2s[j])){

			ECb = EC2s[j];
			PMTb = PMT2s[j];
			CHb = ch2s[j];

			cout << "Message [Number of point]: 2 points only! " << endl;

			const Int_t n = 2;
			Double_t x[n];
			if(q==0) {x[0] = 36.; x[1]=38.;}
			else   {x[0] = 42.; x[1]=52.;}

   	   		Double_t y[n]   = {thr3.arrCh[j]/scale,thr2.arrCh[j]/scale};
   	   		Double_t exl[n] = {0.5,0.5};
   	   		Double_t eyl[n] = {thr3.err1DAC[j]*coefficient/scale,thr2.err1DAC[j]*coefficient/scale};
   	   		Double_t exh[n] = {0.5,0.5};
   	   		Double_t eyh[n] = {thr3.err2DAC[j]*coefficient/scale, thr2.err2DAC[j]*coefficient/scale};
	
			name_graph = TString(thr2.name[j]);

			calibration =  new TGraphAsymmErrors(n,x,y,exl,exh,eyl,eyh);

			name_canvas.clear();	 
			name_canvas = thr2.name[j]; 
		
			vector <string> buffer_name = shift_vector_name_R(thr1.name); 
			thr1.name.swap(buffer_name); 

			vector <int> buffer_EC = shift_vector_int_R(EC1s); 
			EC1s.swap(buffer_EC); 

			vector <int> buffer_PMT = shift_vector_int_R(PMT1s); 
			PMT1s.swap(buffer_PMT); 

			vector <int> buffer_CH = shift_vector_int_R(ch1s); 
			ch1s.swap(buffer_CH); 

  			vector <double> buffer1 = shift_vector_R(thr1.arrCh);
			thr1.arrCh.swap(buffer1);

  			vector <double> buffer2 = shift_vector_R(thr1.arrDAC);
			thr1.arrDAC.swap(buffer2);

  			vector <double> buffer3 = shift_vector_R(thr1.err1DAC);
			thr1.err1DAC.swap(buffer3);

  			vector <double> buffer4 = shift_vector_R(thr1.err2DAC);
			thr1.err2DAC.swap(buffer4);
			
			goto plot;
		}

		if((EC1s[j]<EC2s[j] && EC2s[j]==EC3s[j] &&  PMT1s[j]>PMT2s[j] && PMT2s[j]==PMT3s[j] && ch2s[j]==ch3s[j] && ch1s[j]!=0)||(EC1s[j]==EC2s[j] && EC2s[j]==EC3s[j] && PMT1s[j]<PMT2s[j] && PMT2s[j]==PMT3s[j] && PMT1s[j+1]!=PMT2s[j] && ch2s[j]==ch3s[j] && ch1s[j]!=0)||(EC1s[j]==EC2s[j] && EC2s[j]==EC3s[j] && PMT1s[j]==PMT2s[j] && PMT2s[j]==PMT3s[j] && ch2s[j]==ch3s[j] && ch1s[j]!=0 && ch1s[j]<ch2s[j])){

			ECb = EC2s[j];
			PMTb = PMT2s[j];
			CHb = ch2s[j];

			cout << "Message [Number of point]: 2 points only! " << endl;

			const Int_t n = 2;

			Double_t x[n];
			if(q==0) {x[0] = 36.; x[1]=38.;}
			else  {x[0] = 42.; x[1]=52.;}

   	   		Double_t y[n]   = {thr3.arrCh[j]/scale,thr2.arrCh[j]/scale};
   	   		Double_t exl[n] = {0.5,0.5};
   	   		Double_t eyl[n] = {thr3.err1DAC[j]*coefficient/scale,thr2.err1DAC[j]*coefficient/scale};
   	   		Double_t exh[n] = {0.5,0.5};
   	   		Double_t eyh[n] = {thr3.err2DAC[j]*coefficient/scale, thr2.err2DAC[j]*coefficient/scale};
	
			name_graph = TString(thr2.name[j]);

			calibration =  new TGraphAsymmErrors(n,x,y,exl,exh,eyl,eyh);

			name_canvas.clear();	 
			name_canvas = thr2.name[j]; 

			vector <string> buffer_name;
			vector <int> buffer_EC, buffer_PMT, buffer_CH;
			vector <double> buffer1, buffer2, buffer3, buffer4;
		
			buffer_name = shift_vector_name_L(thr1.name); 
			thr1.name.swap(buffer_name); 

			buffer_EC = shift_vector_int_L(EC1s); 
			EC1s.swap(buffer_EC); 

			buffer_PMT = shift_vector_int_L(PMT1s); 
			PMT1s.swap(buffer_PMT); 

			buffer_CH = shift_vector_int_L(ch1s); 
			ch1s.swap(buffer_CH); 

  			buffer1 = shift_vector_L(thr1.arrCh);
			thr1.arrCh.swap(buffer1);

  			buffer2 = shift_vector_L(thr1.arrDAC);
			thr1.arrDAC.swap(buffer2);

  			buffer3 = shift_vector_L(thr1.err1DAC);
			thr1.err1DAC.swap(buffer3);

  			buffer4 = shift_vector_L(thr1.err2DAC);
			thr1.err2DAC.swap(buffer4);

			buffer_name.erase(buffer_name.begin(), buffer_name.end());
			buffer_EC.erase(buffer_EC.begin(), buffer_EC.end());
			buffer_PMT.erase(buffer_PMT.begin(), buffer_PMT.end());
			buffer_CH.erase(buffer_CH.begin(), buffer_CH.end());
			buffer1.erase(buffer1.begin(), buffer1.end());
			buffer2.erase(buffer1.begin(), buffer1.end());
			buffer3.erase(buffer1.begin(), buffer1.end());
			buffer4.erase(buffer1.begin(), buffer1.end());
		
			buffer_name = shift_vector_name_R(thr2.name); 
			thr2.name.swap(buffer_name); 

			buffer_EC = shift_vector_int_R(EC2s); 
			EC2s.swap(buffer_EC); 

			buffer_PMT = shift_vector_int_R(PMT2s); 
			PMT2s.swap(buffer_PMT); 

			buffer_CH = shift_vector_int_R(ch2s); 
			ch2s.swap(buffer_CH); 

  			buffer1 = shift_vector_R(thr2.arrCh);
			thr2.arrCh.swap(buffer1);

  			buffer2 = shift_vector_R(thr2.arrDAC);
			thr2.arrDAC.swap(buffer2);

  			buffer3 = shift_vector_R(thr2.err1DAC);
			thr2.err1DAC.swap(buffer3);

  			buffer4 = shift_vector_R(thr2.err2DAC);
			thr2.err2DAC.swap(buffer4);

			buffer_name.erase(buffer_name.begin(), buffer_name.end());
			buffer_EC.erase(buffer_EC.begin(), buffer_EC.end());
			buffer_PMT.erase(buffer_PMT.begin(), buffer_PMT.end());
			buffer_CH.erase(buffer_CH.begin(), buffer_CH.end());
			buffer1.erase(buffer1.begin(), buffer1.end());
			buffer2.erase(buffer1.begin(), buffer1.end());
			buffer3.erase(buffer1.begin(), buffer1.end());
			buffer4.erase(buffer1.begin(), buffer1.end());
		
			buffer_name = shift_vector_name_R(thr3.name); 
			thr3.name.swap(buffer_name); 

			buffer_EC = shift_vector_int_R(EC3s); 
			EC3s.swap(buffer_EC); 

			buffer_PMT = shift_vector_int_R(PMT3s); 
			PMT3s.swap(buffer_PMT); 

			buffer_CH = shift_vector_int_R(ch3s); 
			ch3s.swap(buffer_CH); 

  			buffer1 = shift_vector_R(thr3.arrCh);
			thr3.arrCh.swap(buffer1);

  			buffer2 = shift_vector_R(thr3.arrDAC);
			thr3.arrDAC.swap(buffer2);

  			buffer3 = shift_vector_R(thr3.err1DAC);
			thr3.err1DAC.swap(buffer3);

  			buffer4 = shift_vector_R(thr3.err2DAC);
			thr3.err2DAC.swap(buffer4);
		
			goto plot;

		}


		if((EC2s[j]>EC1s[j] && EC1s[j]==EC3s[j] &&  PMT2s[j]<PMT1s[j] && PMT1s[j]==PMT3s[j] && ch1s[j]==ch3s[j] && ch1s[j]!=0)||(EC1s[j]==EC2s[j] && EC2s[j]==EC3s[j] && PMT2s[j]>PMT1s[j] && PMT1s[j]==PMT3s[j] && (PMT2s[j-1]!=PMT1s[j] || ch2s[j-1]==63 || ch2s[j-1]==64) && ch1s[j]==ch3s[j] && ch1s[j]!=0)||(EC1s[j]==EC2s[j] && EC2s[j]==EC3s[j] && PMT2s[j]==PMT1s[j] && PMT1s[j]==PMT3s[j] && ch2s[j]>ch1s[j] && ch1s[j]==ch3s[j] && ch1s[j]!=0)){

			ECb = EC1s[j];
			PMTb = PMT1s[j];
			CHb = ch1s[j];

			cout << "Message [Number of point]: 2 points only! " << endl;

			const Int_t n = 2;
			Double_t x[n];
			if(q==0) {x[0] = 36.; x[1]=42.;}
			else {x[0] = 42.; x[1]=62.;}

   	   		Double_t y[n]   = {thr3.arrCh[j]/scale,thr1.arrCh[j]/scale};
   	   		Double_t exl[n] = {0.5,0.5};
   	   		Double_t eyl[n] = {thr3.err1DAC[j]*coefficient/scale, thr1.err1DAC[j]*coefficient/scale};
   	   		Double_t exh[n] = {0.5,0.5};
   	   		Double_t eyh[n] = {thr3.err2DAC[j]*coefficient/scale, thr1.err2DAC[j]*coefficient/scale};

			name_graph = TString(thr1.name[j]);

			calibration =  new TGraphAsymmErrors(n,x,y,exl,exh,eyl,eyh);

			name_canvas.clear();	 
			name_canvas = thr1.name[j]; 


			vector <string> buffer_name;
			vector <int> buffer_EC, buffer_PMT, buffer_CH;
			vector <double> buffer1, buffer2, buffer3, buffer4;
		
			buffer_name = shift_vector_name_R(thr2.name); 
			thr2.name.swap(buffer_name); 

			buffer_EC = shift_vector_int_R(EC2s); 
			EC2s.swap(buffer_EC); 

			buffer_PMT = shift_vector_int_R(PMT2s); 
			PMT2s.swap(buffer_PMT); 

			buffer_CH = shift_vector_int_R(ch2s); 
			ch2s.swap(buffer_CH); 

  			buffer1 = shift_vector_R(thr2.arrCh);
			thr2.arrCh.swap(buffer1);

  			buffer2 = shift_vector_R(thr2.arrDAC);
			thr2.arrDAC.swap(buffer2);

  			buffer3 = shift_vector_R(thr2.err1DAC);
			thr2.err1DAC.swap(buffer3);

  			buffer4 = shift_vector_R(thr2.err2DAC);
			thr2.err2DAC.swap(buffer4);
			
			goto plot;
		}

		if((EC2s[j]<EC1s[j] && EC1s[j]==EC3s[j] &&  PMT2s[j]>PMT1s[j] && PMT1s[j]==PMT3s[j] && ch1s[j]==ch3s[j] && ch1s[j]!=0)||(EC1s[j]==EC2s[j] && EC2s[j]==EC3s[j] && PMT2s[j]<PMT1s[j] && PMT1s[j]==PMT3s[j] && PMT2s[j+1]!=PMT1s[j] && ch1s[j]==ch3s[j] && ch1s[j]!=0)||(EC1s[j]==EC2s[j] && EC2s[j]==EC3s[j] && PMT2s[j]==PMT1s[j] && PMT1s[j]==PMT3s[j] && ch2s[j]<ch1s[j] && ch1s[j]==ch3s[j] && ch1s[j]!=0)){

			ECb = EC1s[j];
			PMTb = PMT1s[j];
			CHb = ch1s[j];

			cout << "Message [Number of point]: 2 points only! " << endl;

			const Int_t n = 2;
			Double_t x[n];
			if(q==0) {x[0] = 36.; x[1]=42.;}
			else {x[0] = 42.; x[1]=62.;}

   	   		Double_t y[n]   = {thr3.arrCh[j]/scale,thr1.arrCh[j]/scale};
   	   		Double_t exl[n] = {0.5,0.5};
   	   		Double_t eyl[n] = {thr3.err1DAC[j]*coefficient/scale, thr1.err1DAC[j]*coefficient/scale};
   	   		Double_t exh[n] = {0.5,0.5};
   	   		Double_t eyh[n] = {thr3.err2DAC[j]*coefficient/scale, thr1.err2DAC[j]*coefficient/scale};

			name_graph = TString(thr1.name[j]);

			calibration =  new TGraphAsymmErrors(n,x,y,exl,exh,eyl,eyh);

			name_canvas.clear();	 
			name_canvas = thr1.name[j]; 

			vector <string> buffer_name;
			vector <int> buffer_EC, buffer_PMT, buffer_CH;
			vector <double> buffer1, buffer2, buffer3, buffer4;

			buffer_name = shift_vector_name_L(thr2.name); 
			thr2.name.swap(buffer_name); 

			buffer_EC = shift_vector_int_L(EC2s); 
			EC2s.swap(buffer_EC); 

			buffer_PMT = shift_vector_int_L(PMT2s); 
			PMT2s.swap(buffer_PMT); 

			buffer_CH = shift_vector_int_L(ch2s); 
			ch2s.swap(buffer_CH); 

  			buffer1 = shift_vector_L(thr2.arrCh);
			thr2.arrCh.swap(buffer1);

  			buffer2 = shift_vector_L(thr2.arrDAC);
			thr2.arrDAC.swap(buffer2);

  			buffer3 = shift_vector_L(thr2.err1DAC);
			thr2.err1DAC.swap(buffer3);

  			buffer4 = shift_vector_L(thr2.err2DAC);
			thr2.err2DAC.swap(buffer4);

			buffer_name.erase(buffer_name.begin(), buffer_name.end());
			buffer_EC.erase(buffer_EC.begin(), buffer_EC.end());
			buffer_PMT.erase(buffer_PMT.begin(), buffer_PMT.end());
			buffer_CH.erase(buffer_CH.begin(), buffer_CH.end());
			buffer1.erase(buffer1.begin(), buffer1.end());
			buffer2.erase(buffer1.begin(), buffer1.end());
			buffer3.erase(buffer1.begin(), buffer1.end());
			buffer4.erase(buffer1.begin(), buffer1.end());
		
			buffer_name = shift_vector_name_R(thr1.name); 
			thr1.name.swap(buffer_name); 

			buffer_EC = shift_vector_int_R(EC1s); 
			EC1s.swap(buffer_EC); 

			buffer_PMT = shift_vector_int_R(PMT1s); 
			PMT1s.swap(buffer_PMT); 

			buffer_CH = shift_vector_int_R(ch1s); 
			ch1s.swap(buffer_CH); 

  			buffer1 = shift_vector_R(thr1.arrCh);
			thr1.arrCh.swap(buffer1);

  			buffer2 = shift_vector_R(thr1.arrDAC);
			thr1.arrDAC.swap(buffer2);

  			buffer3 = shift_vector_R(thr1.err1DAC);
			thr1.err1DAC.swap(buffer3);

  			buffer4 = shift_vector_R(thr1.err2DAC);
			thr1.err2DAC.swap(buffer4);

			buffer_name.erase(buffer_name.begin(), buffer_name.end());
			buffer_EC.erase(buffer_EC.begin(), buffer_EC.end());
			buffer_PMT.erase(buffer_PMT.begin(), buffer_PMT.end());
			buffer_CH.erase(buffer_CH.begin(), buffer_CH.end());
			buffer1.erase(buffer1.begin(), buffer1.end());
			buffer2.erase(buffer1.begin(), buffer1.end());
			buffer3.erase(buffer1.begin(), buffer1.end());
			buffer4.erase(buffer1.begin(), buffer1.end());
		
			buffer_name = shift_vector_name_R(thr3.name); 
			thr3.name.swap(buffer_name); 

			buffer_EC = shift_vector_int_R(EC3s); 
			EC3s.swap(buffer_EC); 

			buffer_PMT = shift_vector_int_R(PMT3s); 
			PMT3s.swap(buffer_PMT); 

			buffer_CH = shift_vector_int_R(ch3s); 
			ch3s.swap(buffer_CH); 

  			buffer1 = shift_vector_R(thr3.arrCh);
			thr3.arrCh.swap(buffer1);

  			buffer2 = shift_vector_R(thr3.arrDAC);
			thr3.arrDAC.swap(buffer2);

  			buffer3 = shift_vector_R(thr3.err1DAC);
			thr3.err1DAC.swap(buffer3);

  			buffer4 = shift_vector_R(thr3.err2DAC);
			thr3.err2DAC.swap(buffer4);
		
			goto plot;

		}

		if((EC3s[j]>EC2s[j] && EC1s[j]==EC3s[j] &&  PMT3s[j]<PMT2s[j] && PMT1s[j]==PMT3s[j] && ch1s[j]==ch3s[j] && ch1s[j]!=0)||(EC1s[j]==EC2s[j] && EC2s[j]==EC3s[j] && PMT3s[j]>PMT2s[j] && PMT1s[j]==PMT2s[j] && (PMT3s[j-1]!=PMT1s[j] || ch3s[j-1]==64 || ch3s[j-1]==63) && ch1s[j]==ch2s[j] && ch1s[j]!=0)||(EC1s[j]==EC2s[j] && EC2s[j]==EC3s[j] && PMT3s[j]==PMT1s[j] && PMT1s[j]==PMT2s[j] && ch3s[j]>ch2s[j] && ch1s[j]==ch2s[j] && ch1s[j]!=0)){

			ECb = EC1s[j];
			PMTb = PMT1s[j];
			CHb = ch1s[j];

			cout << "Message [Number of point]: 2 points only! " << endl;

			const Int_t n = 2;
			Double_t x[n];
			if(q==0) {x[0] = 38.; x[1]=42.;}
			else   {x[0] = 52.; x[1]=62.;}

   	   		Double_t y[n]   = {thr2.arrCh[j]/scale,thr1.arrCh[j]/scale};
   	   		Double_t exl[n] = {0.5,0.5};
   	   		Double_t eyl[n] = {thr2.err1DAC[j]*coefficient/scale,thr1.err1DAC[j]*coefficient/scale};
   	   		Double_t exh[n] = {0.5,0.5};
   	   		Double_t eyh[n] = { thr2.err2DAC[j]*coefficient/scale, thr1.err2DAC[j]*coefficient/scale};

			name_graph = TString(thr1.name[j]);

			calibration =  new TGraphAsymmErrors(n,x,y,exl,exh,eyl,eyh);

			name_canvas.clear();	 
			name_canvas = thr1.name[j]; 
		
			vector <string> buffer_name = shift_vector_name_R(thr3.name); 
			thr3.name.swap(buffer_name); 
		
			vector <int> buffer_EC = shift_vector_int_R(EC3s); 
			EC3s.swap(buffer_EC); 

			vector <int> buffer_PMT = shift_vector_int_R(PMT3s); 
			PMT3s.swap(buffer_PMT); 

			vector <int> buffer_CH = shift_vector_int_R(ch3s); 
			ch3s.swap(buffer_CH); 

  			vector <double> buffer1 = shift_vector_R(thr3.arrCh);
			thr3.arrCh.swap(buffer1);

  			vector <double> buffer2 = shift_vector_R(thr3.arrDAC);
			thr3.arrDAC.swap(buffer2);

  			vector <double> buffer3 = shift_vector_R(thr3.err1DAC);
			thr3.err1DAC.swap(buffer3);

  			vector <double> buffer4 = shift_vector_R(thr3.err2DAC);
			thr3.err2DAC.swap(buffer4);
			
			goto plot;
		}

		if((EC3s[j]<EC1s[j] && EC1s[j]==EC2s[j] &&  PMT1s[j]<PMT3s[j] && PMT2s[j]==PMT1s[j] && ch2s[j]==ch1s[j] && ch1s[j]!=0)||(EC1s[j]==EC2s[j] && EC2s[j]==EC3s[j] && PMT3s[j]<PMT1s[j] && PMT1s[j]==PMT2s[j] && PMT3s[j+1]!=PMT1s[j] && ch1s[j]==ch2s[j] && ch1s[j]!=0)||(EC1s[j]==EC2s[j] && EC2s[j]==EC3s[j] && PMT3s[j]==PMT1s[j] && PMT1s[j]==PMT2s[j] && ch3s[j]<ch2s[j] && ch1s[j]==ch2s[j])){

			ECb = EC1s[j];
			PMTb = PMT1s[j];
			CHb = ch1s[j];

			cout << "Message [Number of point]: 2 points only! " << endl;

			const Int_t n = 2;
			Double_t x[n];
			if(q==0) {x[0] = 38.; x[1]=42.;}
			else  {x[0] = 20.; x[1]=30.;}


   	   		Double_t y[n]   = {thr2.arrCh[j]/scale,thr1.arrCh[j]/scale};
   	   		Double_t exl[n] = {0.5,0.5};
   	   		Double_t eyl[n] = {thr2.err1DAC[j]*coefficient/scale,thr1.err1DAC[j]*coefficient/scale};
   	   		Double_t exh[n] = {0.5,0.5};
   	   		Double_t eyh[n] = { thr2.err2DAC[j]*coefficient/scale, thr1.err2DAC[j]*coefficient/scale};

			name_graph = TString(thr1.name[j]);

			calibration =  new TGraphAsymmErrors(n,x,y,exl,exh,eyl,eyh);

			name_canvas.clear();	 
			name_canvas = thr1.name[j]; 

			vector <string> buffer_name;
			vector <int> buffer_EC, buffer_PMT, buffer_CH;
			vector <double> buffer1, buffer2, buffer3, buffer4;
		
			buffer_name = shift_vector_name_L(thr3.name); 
			thr3.name.swap(buffer_name); 

			buffer_EC = shift_vector_int_L(EC3s); 
			EC3s.swap(buffer_EC); 

			buffer_PMT = shift_vector_int_L(PMT3s); 
			PMT3s.swap(buffer_PMT); 

			buffer_CH = shift_vector_int_L(ch3s); 
			ch3s.swap(buffer_CH); 

  			buffer1 = shift_vector_L(thr3.arrCh);
			thr3.arrCh.swap(buffer1);

  			buffer2 = shift_vector_L(thr3.arrDAC);
			thr3.arrDAC.swap(buffer2);

  			buffer3 = shift_vector_L(thr3.err1DAC);
			thr3.err1DAC.swap(buffer3);

  			buffer4 = shift_vector_L(thr3.err2DAC);
			thr3.err2DAC.swap(buffer4);

			buffer_name.erase(buffer_name.begin(), buffer_name.end());
			buffer_EC.erase(buffer_EC.begin(), buffer_EC.end());
			buffer_PMT.erase(buffer_PMT.begin(), buffer_PMT.end());
			buffer_CH.erase(buffer_CH.begin(), buffer_CH.end());
			buffer1.erase(buffer1.begin(), buffer1.end());
			buffer2.erase(buffer1.begin(), buffer1.end());
			buffer3.erase(buffer1.begin(), buffer1.end());
			buffer4.erase(buffer1.begin(), buffer1.end());
		
			buffer_name = shift_vector_name_R(thr1.name); 
			thr1.name.swap(buffer_name); 

			buffer_EC = shift_vector_int_R(EC1s); 
			EC1s.swap(buffer_EC); 

			buffer_PMT = shift_vector_int_R(PMT1s); 
			PMT1s.swap(buffer_PMT); 
	
			buffer_CH = shift_vector_int_R(ch1s); 
			ch1s.swap(buffer_CH); 

  			buffer1 = shift_vector_R(thr1.arrCh);
			thr1.arrCh.swap(buffer1);

  			buffer2 = shift_vector_R(thr1.arrDAC);
			thr1.arrDAC.swap(buffer2);

  			buffer3 = shift_vector_R(thr1.err1DAC);
			thr1.err1DAC.swap(buffer3);

  			buffer4 = shift_vector_R(thr1.err2DAC);
			thr1.err2DAC.swap(buffer4);

			buffer_name.erase(buffer_name.begin(), buffer_name.end());
			buffer_EC.erase(buffer_EC.begin(), buffer_EC.end());
			buffer_PMT.erase(buffer_PMT.begin(), buffer_PMT.end());
			buffer_CH.erase(buffer_CH.begin(), buffer_CH.end());
			buffer1.erase(buffer1.begin(), buffer1.end());
			buffer2.erase(buffer1.begin(), buffer1.end());
			buffer3.erase(buffer1.begin(), buffer1.end());
			buffer4.erase(buffer1.begin(), buffer1.end());
		
			buffer_name = shift_vector_name_R(thr2.name); 
			thr2.name.swap(buffer_name); 

			buffer_EC = shift_vector_int_R(EC2s); 
			EC2s.swap(buffer_EC); 

			buffer_PMT = shift_vector_int_R(PMT2s); 
			PMT2s.swap(buffer_PMT); 

			buffer_CH = shift_vector_int_R(ch2s); 
			ch2s.swap(buffer_CH); 

  			buffer1 = shift_vector_R(thr2.arrCh);
			thr2.arrCh.swap(buffer1);

  			buffer2 = shift_vector_R(thr2.arrDAC);
			thr2.arrDAC.swap(buffer2);

  			buffer3 = shift_vector_R(thr2.err1DAC);
			thr2.err1DAC.swap(buffer3);

  			buffer4 = shift_vector_R(thr2.err2DAC);
			thr2.err2DAC.swap(buffer4);
		
			goto plot;

		}

		//for 3 points!	
		if(EC1s[j]==EC2s[j] && EC2s[j]==EC3s[j] && PMT3s[j]==PMT1s[j] && PMT1s[j]==PMT2s[j] && ch1s[j]==ch2s[j] && ch1s[j]!=0 && ch3s[j]==ch2s[j]){

			ECb = EC1s[j];
			PMTb = PMT1s[j];
			CHb = ch1s[j];

			cout << "Message [Number of point]: all 3 points! " << endl;

			name_graph = TString(thr1.name[j]);
			name_canvas.clear();	 
			name_canvas = thr1.name[j]; 

			const Int_t n = 3;
			Double_t x[n];
			if(q==0) {x[0] = 36.; x[1]=38.; x[2]=42;}
			else   {x[0]=42; x[1] = 52.; x[2]=62.;}

			Double_t y[n]   = {thr3.arrCh[j]/scale,thr2.arrCh[j]/scale,thr1.arrCh[j]/scale};
   	   		Double_t exl[n] = {0.5,0.5,0.5};
   	   		Double_t eyl[n] = {thr3.err1DAC[j]*coefficient/scale,thr2.err1DAC[j]*coefficient/scale,thr1.err1DAC[j]*coefficient/scale};
   	   		Double_t exh[n] = {0.5,0.5,0.5};
   	   		Double_t eyh[n] = {thr3.err2DAC[j]*coefficient/scale, thr2.err2DAC[j]*coefficient/scale, thr1.err2DAC[j]*coefficient/scale};

			calibration = new TGraphAsymmErrors(n,x,y,exl,exh,eyl,eyh);

			goto plot;
		}



		plot: 
                
		cout << "Message [calibration]: " << TString(name_canvas) << endl;

		double x_p, y_p;
		calibration -> GetPoint(0, x_p, y_p);
		calibration -> GetPoint(1, x_p, y_p);
		calibration -> GetPoint(2, x_p, y_p);

	   	TCanvas *c = new TCanvas("c", TString("Calibration"+name_canvas), 600, 600);
		c -> cd();
		
		calibration -> SetTitle(name_graph);
		calibration -> GetXaxis() -> SetTitle("Threshold code"); 
		calibration -> GetYaxis() -> SetTitle("Charge, ke^{-}");
		calibration -> GetYaxis() -> SetMaxDigits(3);
		calibration -> Draw("AP");

		TF1 *f1 = new TF1("f1",fitF,0,63,2);
		f1 -> SetParameters(1, 0.2);
		calibration -> Fit("f1", "qs");

		TFitResultPtr r = calibration -> Fit("f1", "s");
		cout << "Message [calibration fit] fit status: " << r << endl;

		TF1 *myfit = (TF1*)calibration->GetFunction("f1");
		myfit -> SetLineWidth(3);
	
		offset = myfit  -> GetParameter(0);
		slope = myfit  -> GetParameter(1);

		if(slope<100000 && TMath::Abs(offset)<100000){  //100 for scale = 1000 (Me-) //Hack: normal adequate fit results
				double offset_error = myfit -> GetParError(0);
				double slope_error = myfit -> GetParError(1);
				
				if(offset_error/offset > max_error) count_HE_off++; //how many more than 10% error
				if(slope_error/slope > max_error) count_HE_sl++;    //how many more than 10% error
				
				av_slope_all+=slope;	
				av_at_zero_all+=TMath::Abs(offset);

				//Board_13
				if(g==0){
					histo_slope_all_13 -> Fill(slope);
					histo_at_zero_all_13 -> Fill(TMath::Abs(offset));
					histo_slope_error_13 -> Fill(slope_error);
					histo_at_zero_error_13 -> Fill(offset_error);
					slope_13[q][ECb][PMTb][CHb-1]= slope;
					offset_13[q][ECb][PMTb][CHb-1]= offset;
				}

				//Board_06
				if(g==1){
					histo_slope_all_06 -> Fill(slope);
					histo_at_zero_all_06 -> Fill(TMath::Abs(offset));
					histo_slope_error_06 -> Fill(slope_error);
					histo_at_zero_error_06 -> Fill(offset_error);
					slope_06[q][ECb][PMTb][CHb-1]= slope;
					offset_06[q][ECb][PMTb][CHb-1]= offset;
				}
			}

		else{
			TFitResultPtr r = calibration -> Fit("pol1", "s");
			cout << "Message [calibration fit] fit status: " << r << endl;
			TF1 *myfit2 = (TF1*)calibration->GetFunction("pol1");
			myfit2 -> SetLineWidth(3);
	
			offset = myfit2  -> GetParameter(0);
			slope = myfit2  -> GetParameter(1);
		
			if(slope<100000 && TMath::Abs(offset)<100000){  //100 for scale = 1000 (Me-)

				double offset_error = myfit2 -> GetParError(0);
				double slope_error = myfit2 -> GetParError(1);
			
				if(offset_error/offset > max_error) count_HE_off++; //how many more than 10% error
				if(slope_error/slope > max_error) count_HE_sl++;    //how many more than 10% error

				av_slope_all+=slope;	
				av_at_zero_all+=TMath::Abs(offset);

				//Board 13

				if(g==0){
					histo_slope_all_13 -> Fill(slope);
					histo_at_zero_all_13 -> Fill(TMath::Abs(offset));
					histo_slope_error_13 -> Fill(slope_error);
					histo_at_zero_error_13 -> Fill(offset_error);
					slope_13[q][ECb][PMTb][CHb-1]= slope;
					offset_13[q][ECb][PMTb][CHb-1]= offset;
				}


				//Board_06
				if(g==1){
					histo_slope_all_06 -> Fill(slope);
					histo_at_zero_all_06 -> Fill(TMath::Abs(offset));
					histo_slope_error_06 -> Fill(slope_error);
					histo_at_zero_error_06 -> Fill(offset_error);
					slope_06[q][ECb][PMTb][CHb-1]= slope;
					offset_06[q][ECb][PMTb][CHb-1]= offset;
				}	
			}

		 else { 
			cout << "Message [calibration fit]: calibration failed " << TString(name_canvas) << endl; 
			bad_calibr_names.push_back(thr1.name[j]); 
			++bad_calibr;
			if(g==0){slope_13[q][ECb][PMTb][CHb-1]= -1; offset_13[q][ECb][PMTb][CHb-1]= 0;}
			if(g==1){slope_06[q][ECb][PMTb][CHb-1]= -1; offset_06[q][ECb][PMTb][CHb-1]= 0;}
		}
	}
		cout << "Message [calibration fit] slope: " << slope << ", ke- " << endl;
		cout << "Message [calibration fit] offset: " << offset << ", ke-" << endl;

		auto legend = new TLegend(0.1,0.7,0.48,0.9);
   		legend->SetHeader(TString(head[0] + " " + head[1]),"C"); // option "C" allows to center the header
   		legend->AddEntry((TObject*)0,TString::Format("#Delta Th = %.1f ke^{-}", slope),"");
   		legend->AddEntry((TObject*)0,TString::Format("|Th(0)| = %.1f ke^{-}", TMath::Abs(offset)),"");
      		legend->Draw();
		
		c -> SaveAs(dir+name_histo1+ "/calibration/" + TString(name_canvas)+".pdf");
		c -> Close();

		if(file_par.is_open()) { file_par << name_canvas << setw(30) << att << setw(30) << off << setw(30) << slope << setw(30) << TMath::Abs(offset) << endl;}	
	}

	file_par.close();

	histo_count_HE_off -> Fill(count_HE_off);
	histo_count_HE_sl -> Fill(count_HE_sl);

	if ( calibration_file.IsOpen() ) {
					cout << "Message [ROOT file state]: writing down" << endl;
					string board_num;
					if(g==0) {nentries_sl =  histo_slope_all_13-> GetEntries(); nentries_off =  histo_at_zero_all_13-> GetEntries(); board_num = "13"; histo_slope_all_13->Write(TString("Board_" + board_num + head[0] + "_" + head[1] +"_"+"slope_all")); histo_at_zero_all_13->Write(TString("Board_" + board_num + head[0] + "_" + head[1] +"_"+"offset_all")); histo_slope_error_13->Write(TString("Board_" + board_num + head[0] + "_" + head[1] +"_"+"slope_error")); histo_at_zero_error_13->Write(TString("Board_" + board_num + head[0] + "_" + head[1] +"_"+"offset_error"));}
					if(g==1) {nentries_sl =  histo_slope_all_06-> GetEntries(); nentries_off =  histo_at_zero_all_06-> GetEntries(); board_num = "06"; histo_slope_all_06->Write(TString("Board_" + board_num + head[0] + "_" + head[1] +"_"+"slope_all")); histo_at_zero_all_06->Write(TString("Board_" + board_num + head[0] + "_" + head[1] +"_"+"offset_all")); histo_slope_error_06->Write(TString("Board_" + board_num + head[0] + "_" + head[1] +"_"+"slope_error")); histo_at_zero_error_06->Write(TString("Board_" + board_num + head[0] + "_" + head[1] +"_"+"offset_error"));}	
					}

	cout << "Message [calibration fit] number of bad calibrations: " << bad_calibr << endl;
	cout << "Message [calibration fit] list of bad calibrations: " << endl;
	for(auto j = bad_calibr_names.begin(); j!=bad_calibr_names.end(); j++) cout << *j << endl;

  	ofstream file;

    	file.open(dir +name_histo1+"/"+TString::Format("avaraged_results.txt"), ios::out);
   	file << "Name" << setw(30) << "Attenuation" << setw(30) << "Offset" << setw(30) << "all, ke-" << endl; 
	file << "slop" << setw(30) << att << setw(30) << off << setw(30) << av_slope_all/(nentries_sl-bad_calibr) << endl;    
	file << "offs" << setw(30) << att << setw(30) << off << setw(30) << av_at_zero_all/(nentries_off-bad_calibr) << endl;    
	file.close();

	EC1s.clear();	
	EC2s.clear();
	EC3s.clear();
	PMT1s.clear();	
	PMT2s.clear();
	PMT3s.clear();
	ch1s.clear();	
	ch2s.clear();
	ch3s.clear();
	bad_calibr_names.clear();
	head.clear();

	thr1.clear();
	thr2.clear();
	thr3.clear();

	nentries_sl = 0;
	nentries_off = 0;
	bad_calibr = 0;
	av_slope_all = 0;
	av_at_zero_all = 0;
	count_MV = 0;
	count_HE_off = 0;
	count_HE_sl = 0;

	delete histo_slope_all_13;
	delete histo_slope_all_06;
	delete histo_at_zero_all_13;
	delete histo_at_zero_all_06;
	delete histo_slope_error_06;
	delete histo_slope_error_13;
	delete histo_at_zero_error_06;
	delete histo_at_zero_error_13;
		}
		if ( calibration_file.IsOpen() ){
		cout << "Message [ROOT file state]: writing down" << endl;
		if(g==0){for(int t = 0; t<4; t++){for(int k = 0; k<2; k++){for(int l = 0; l<=64; l++){
				EC13_v=t; PMT13_v = k; CH13_v = l;
				for(int m = 0; m<5 ; m++){
					if(slope_13[m][t][k][l]<0.001 || slope_13[m][t][k][l] > 100000) slope_13[m][t][k][l] = -1; 
					if(m==0) {thresholdstep_13_att3_v = slope_13[m][t][k][l]; offset_13_att3_v = offset_13[m][t][k][l];}
					if(m==1) {thresholdstep_13_att2_v = slope_13[m][t][k][l]; offset_13_att2_v = offset_13[m][t][k][l];}					
					if(m==2) {thresholdstep_13_att1_v = slope_13[m][t][k][l]; offset_13_att1_v = offset_13[m][t][k][l];}
					if(m==3) {thresholdstep_13_att0_v = slope_13[m][t][k][l]; offset_13_att0_v = offset_13[m][t][k][l];}
					if(m==4) {thresholdstep_13_off0_v = slope_13[m][t][k][l]; offset_13_off0_v = offset_13[m][t][k][l];}
					}
				calibration13 -> Fill();
				}}}}
		if(g==1){for(int t = 0; t<4; t++){for(int k = 0; k<2; k++){for(int l = 0; l<64; l++){
				EC06_v=t; PMT06_v = k; CH06_v = l;
				for(int m = 0; m<5 ; m++){
					if(slope_06[m][t][k][l]<0.001 || slope_06[m][t][k][l] > 100000) slope_06[m][t][k][l] = -1; 
					if(m==0) {thresholdstep_06_att3_v = slope_06[m][t][k][l]; offset_06_att3_v = offset_06[m][t][k][l];}
					if(m==1) {thresholdstep_06_att2_v = slope_06[m][t][k][l]; offset_06_att2_v = offset_06[m][t][k][l];}					
					if(m==2) {thresholdstep_06_att1_v = slope_06[m][t][k][l]; offset_06_att1_v = offset_06[m][t][k][l];}
					if(m==3) {thresholdstep_06_att0_v = slope_06[m][t][k][l]; offset_06_att0_v = offset_06[m][t][k][l];}
					if(m==4) {thresholdstep_06_off0_v = slope_06[m][t][k][l]; offset_06_off0_v = offset_06[m][t][k][l];}
					}
				calibration06 -> Fill();
				}}}}
		}
		
 	}
        
	histo_count_MV -> Write();
	histo_count_HE_off -> Write();
	histo_count_HE_sl -> Write();
	calibration_file.Write();
	return 0;
};
