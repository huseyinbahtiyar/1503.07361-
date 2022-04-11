
#include <TSystem.h>
#include <TROOT.h>
#include "TRint.h"
#include <TFile.h>
#include <TTree.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TMultiGraph.h>
#include <TRandom3.h>
#include <TEventList.h>
#include <TChain.h>
#include <TH1D.h>
#include <TF1.h>
#include "TVector2.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TCanvas.h"
#include "TMatrixD.h"
#include "TMath.h"
#include "TStyle.h"
#include "time.h"
#include <map>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <math.h>
#include <string>
#include <vector>
#include <stddef.h>
#include <cstring>
#include <stdexcept>
#include <memory>
#include <new>
#include "RVersion.h"
#include "TRandom.h"
#include "TDirectory.h"
#include "TProcessID.h"
#include "TLegend.h"
#include "trQCD.h"
#include "TMatrixDEigen.h"
#include <TVector.h>
#include "TMatrixDSym.h"
#include "TMatrixDSymEigen.h"
#include "TMatrixDLazy.h"
#include "TVectorD.h"
#include "TDecompLU.h"
#include "TDecompSVD.h"
#include "TVectorDfwd.h"
#include <iostream>
#include <fstream>
//////////////////////////////////////////////////////////////////////////
//                                                                      //
// trQCD Class File                                                     //
//    Written by Huseyin Bahtiyar                                       //
// 			3.01.15                                                     //
//                                                                      //
//////////////////////////////////////////////////////////////////////////


	
trQCD::trQCD()
{

}
trQCD::~trQCD()
{

}

/// 	Find Probabilties of given q square value
/** Find Probabilties of given q square value, For 2pt functions */
int
trQCD::findProbs(int psqr, string data_name, vector < string > &output)
{
	int		divider = 0;
	//cout << "P square selected :" << psqr << endl;
	//cout << "Now looking for possibilities..." << endl;
	for (int px = -3; px < 4; px++) {
		for (int py = -3; py < 4; py++) {
			for (int pz = -3; pz < 4; pz++) {
				if (((px * px) + (py * py) + (pz * pz)) == psqr) {
					std::ostringstream oss;
					oss << data_name << "_px" << px << "_py" << py << "_pz" << pz;
					std::string tempr = oss.str();
					output.push_back(tempr);
					divider++;
				}
			}
		}
	}
	//cout << divider << " possibilities found." << endl;
	return (divider);
}


/// 	Find Probabilties of given q square value
/** Find Probabilties of given q square value, For 3pt functions */
int
trQCD::findProbs3pt(int psqr, string data_name,int pol,int snk ,int current,int src ,vector < string > &output)
{
	//char		tempr   [100];
	int		divider = 0;
	//cout << "P square selected :" << psqr << endl;
	//cout << "Now looking for possibilities..." << endl;
	for (int qx = -3; qx < 4; qx++) {
		for (int qy = -3; qy < 4; qy++) {
			for (int qz = -3; qz < 4; qz++) {
				if (((qx * qx) + (qy * qy) + (qz * qz)) == psqr) {
					std::ostringstream oss;
					oss << data_name <<"_p" << pol <<"_snk"<< snk<<"_g"<<current<<"_src"<<src<<"_qx" << qx << "_qy" << qy << "_qz" << qz;
					std::string tempr = oss.str();
					output.push_back(tempr);
					divider++;
				}
			}
		}
	}
	//cout << divider << " possibilities found." << endl;
	return (divider);
}


/// 	int to string conversion
/** int to string conversion for file namings */
string
trQCD::itos(int i)
{
	stringstream s;
	s << i;
	return s.str();
}

/// 	Read the "not averaged" data from vector
/** Read the "not averaged" data from vector  */
TMatrixD
trQCD::reader(vector < vector<int> > data, string possib_0, vector <string> path)
{
	int nconf=0;
	for(int vcount=0; vcount<data.size(); vcount++){
		nconf += data[vcount].size();
	}
	TMatrixD datamatrix2pt_0(Nt,nconf);
	ifstream	opentheFile2;
	//cout << "Path is: " << path << endl;
	//cout << "Now looking for files..." << endl;
	int		numberfiles = 0;
	int 	column = 0;
	for(int vcount=0; vcount<data.size(); vcount++){
		for (int folder = 0; folder < data[vcount].size(); folder++) {
			string		folder2 = itos(data[vcount][folder]);
			string		dosya   (path[vcount] + folder2 + possib_0);
			//cout<<dosya<<endl;
			opentheFile2.open(dosya.c_str());
			if (opentheFile2.is_open()) {
				for (int j = 0; j < Nt; j++) {
					Double_t	skip;
					opentheFile2 >> skip;
					opentheFile2 >> datamatrix2pt_0(j, column);
					opentheFile2 >> skip;
				}//for j
				column++;
			} //if
			else {
				cout << "File Error " << dosya << " could not open!" << endl;
				exit(1);
			}//else
			opentheFile2.close();
		}//for folder
	}//for vcount
	return datamatrix2pt_0;
}




/// 	Read the "not averaged" data from vector
/** Read the "not averaged" data from vector  */
TMatrixD
trQCD::reader(vector<int> data, string possib_0, string path)
{
	int nconf=data.size();
	TMatrixD datamatrix2pt_0(Nt,nconf);
	ifstream	opentheFile2;
	//cout << "Path is: " << path << endl;
	//cout << "Now looking for files..." << endl;
	int		numberfiles = 0;
	for (int folder = 0; folder < data.size(); folder++) {
		string		folder2 = itos(data[folder]);
		string		dosya   (path + folder2 + possib_0);
		opentheFile2.open(dosya.c_str());
		if (opentheFile2.is_open()) {
			for (int j = 0; j < Nt; j++) {
				Double_t	skip;
				opentheFile2 >> skip;
				opentheFile2 >> datamatrix2pt_0(j, folder);
				opentheFile2 >> skip;
			}
		} else {
			cout << "File Error " << dosya << " could not open!" << endl;
			exit(1);
		}
		opentheFile2.close();

	}
	return datamatrix2pt_0;
}


/// 	Read the "not averaged" data
/** Read the "not averaged" data  */
TMatrixD
trQCD::reader(int beginf, int endf, int increment, string possib_0, string path)
{
	int nconf=(endf-beginf)/increment;
	nconf++;
	TMatrixD datamatrix2pt_0(Nt,nconf);
	ifstream	opentheFile2;
	//cout << "Path is: " << path << endl;
	//cout << "Now looking for files..." << endl;
	int		numberfiles = 0;
	for (int folder = beginf; folder <= endf; folder = folder + increment) {
		int		i = (folder - beginf) / increment;
		string		folder2 = itos(folder);
		string		dosya   (path + folder2 + possib_0);
		opentheFile2.open(dosya.c_str());
		if (opentheFile2.is_open()) {
			for (int j = 0; j < Nt; j++) {
				double		skip;
				opentheFile2 >> skip;
				opentheFile2 >> datamatrix2pt_0(j, i);
				opentheFile2 >> skip;
			}
		} else {
			cout << "File Error " << dosya << " could not open!" << endl;
			exit(1);
		}
		opentheFile2.close();

	}
	return datamatrix2pt_0;
}


/// Read the data from vector
/** Read the data from vector continue from  int cont */
TMatrixD
trQCD::readerC(vector<int> data, string  dataname, string path, int cont, TMatrixD & datamatrix2pt_0)
{
	ifstream	opentheFile2;
	cout << "Path: " << path <<" File : "<<dataname <<endl;
	int		numberfiles = 0;
	cont=cont+1;	//fix
	for (int folder = 0; folder < data.size(); folder++) {
		string		folder2 = itos(data[folder]);
        string		dosya   (path + folder2 + dataname);
		int say = folder + cont;
		opentheFile2.open(dosya.c_str());
        if (opentheFile2.is_open())
        {
            for (int j = 0; j < Nt; j++) {
                double		skip;
                opentheFile2 >> skip;
                opentheFile2 >> datamatrix2pt_0(j, say);
                opentheFile2 >> skip;
            }
        }
        else {
            cout<<"File Error "<<dosya<<" could not open!"<<endl;
            exit(1);
        }
        opentheFile2.close();
	}
	return datamatrix2pt_0;
}


/// Read the data
/** Read the data and continue from  int cont */
TMatrixD
trQCD::readerC(int beginf, int endf, int increment, int cont, string  dataname, string path, TMatrixD & datamatrix2pt_0)
{
	int nconf=(endf-beginf)/increment;
	ifstream	opentheFile2;
	cout << "Path: " << path <<" File : "<<dataname <<endl;
	int		numberfiles = 0;
	cont=cont+1;//fix
	for (int folder = beginf; folder <= endf; folder = folder + increment) {
		int         i = cont+(folder - beginf) / increment;
        string		folder2 = itos(folder);
        string		dosya   (path + folder2 + dataname);
		opentheFile2.open(dosya.c_str());
        if (opentheFile2.is_open())
        {
            for (int j = 0; j < Nt; j++) {
                double		skip;
                opentheFile2 >> skip;
                opentheFile2 >> datamatrix2pt_0(j, i);
                opentheFile2 >> skip;
            }
        }
        else {
            cout<<"File Error "<<dosya<<" could not open!"<<endl;
            exit(1);
        }
        opentheFile2.close();
	}
	return datamatrix2pt_0;
}


/// Reader functions
/** Read the data and calcualte average from vector */
TMatrixD
trQCD::readerAveraged(vector <vector <int> > data, int psqr, vector < string > &possib2pt, int divider, vector <string> path2pt)
{
	int nconf=0;
	for(int vcount=0; vcount<data.size(); vcount++){
		nconf += data[vcount].size();
	}
	TMatrixD datamatrix2pt(Nt,nconf);
	ifstream	opentheFile3;
	int		nbcol = datamatrix2pt.GetNcols();
	//cout << "Path is: " << path2pt << endl;
	//cout << "Now looking for files then calculate the average..." << endl;
	int		numberfiles = 0;
	int		column = 0;
	
	for(int vcount=0; vcount<data.size(); vcount++){
		for (int folder = 0; folder < data[vcount].size(); folder++) {
			for (int it = 0; it < possib2pt.size(); it++) {
				string		folder2 = itos(data[vcount][folder]);
				string		dosya3  (path2pt[vcount] + folder2 + possib2pt[it]);
				opentheFile3.open(dosya3.c_str());
				if (opentheFile3.is_open()) {
					//cout << "Opened " << dosya3 << endl;
					for (int j = 0; j < Nt; j++) {
						double		skip;
						Double_t	val;
						opentheFile3 >> skip;
						opentheFile3 >> val;
						opentheFile3 >> skip;
						datamatrix2pt(j, column) += val;
					}//j
				}//if 
				else {
					cout << "File Error " << dosya3 << " could not open!" << endl;
					exit(1);
				}//else
				opentheFile3.close();
			}//possib
			column++;
		}//folder
	}//vcount
	double		scale = 1. / divider;
	datamatrix2pt = datamatrix2pt * scale;
	return datamatrix2pt;
}


/// Reader functions
/** Read the data and calcualte average from vector */
TMatrixD
trQCD::readerAveraged(vector <int> data, int psqr, vector < string > &possib2pt, int divider, string path2pt)
{
	int nconf = data.size();
	TMatrixD datamatrix2pt(Nt,nconf);
	ifstream	opentheFile3;
	int		nbcol = datamatrix2pt.GetNcols();
	//cout << "Path is: " << path2pt << endl;
	//cout << "Now looking for files then calculate the average..." << endl;
	int		numberfiles = 0;
	for (int folder = 0; folder < data.size(); folder++) {
		for (int it = 0; it < possib2pt.size(); it++) {
			string		folder2 = itos(data[folder]);
			string		dosya3  (path2pt + folder2 + possib2pt[it]);
			opentheFile3.open(dosya3.c_str());
			if (opentheFile3.is_open()) {
				//cout << "Opened " << dosya3 << endl;
				for (int j = 0; j < Nt; j++) {
					double		skip;
					Double_t	val;
					opentheFile3 >> skip;
					opentheFile3 >> val;
					opentheFile3 >> skip;
					datamatrix2pt(j, folder) += val;
				}
			} else {
				cout << "File Error " << dosya3 << " could not open!" << endl;
				exit(1);
			}

			opentheFile3.close();
		}
	}
	double		scale = 1. / divider;
	datamatrix2pt = datamatrix2pt * scale;
	return datamatrix2pt;
}

/// Reader functions
/** Read the data and calcualte average */
TMatrixD
trQCD::readerAveraged(int beginf, int endf, int increment, int psqr, vector < string > &possib2pt, int divider, string path2pt)
{
	int nconf=(endf-beginf)/increment;
	nconf++;
	TMatrixD datamatrix2pt(Nt,nconf);
	ifstream	opentheFile3;
	int		nbcol = datamatrix2pt.GetNcols();
	//cout << "Path is: " << path2pt << endl;
	//cout << "Now looking for files then calculate the average..." << endl;
	int		numberfiles = 0;
	for (int folder = beginf; folder <= endf; folder = folder + increment) {
		int		i = (folder - beginf) / increment;
		for (int it = 0; it < possib2pt.size(); it++) {
			string		folder2 = itos(folder);
			string		dosya3  (path2pt + folder2 + possib2pt[it]);
			opentheFile3.open(dosya3.c_str());
			if (opentheFile3.is_open()) {
				//cout << "Opened " << dosya3 << endl;
				for (int j = 0; j < Nt; j++) {
					double		skip;
					Double_t	val;
					opentheFile3 >> skip;
					opentheFile3 >> val;
					opentheFile3 >> skip;
					datamatrix2pt(j, i) += val;
				}
			} else {
				cout << "File Error " << dosya3 << " could not open!" << endl;
				exit(1);
			}

			opentheFile3.close();
		}
	}
	double		scale = 1. / divider;
	datamatrix2pt = datamatrix2pt * scale;
	return datamatrix2pt;
}


/// Reader functions
/** Read the data and calcualte average from vector */
TMatrixD
trQCD::readerAveragedC(vector <int> data, int psqr, vector < string > &possib2pt, int divider, string path2pt,TMatrixD & datamatrix2pt,int cont)
{
	ifstream	opentheFile3;
	int		nbcol = datamatrix2pt.GetNcols();
	double		scale = 1. / divider;
	//cout << "Path is: " << path2pt << endl;
	cout << "Now looking for files then calculate the average..." << endl;
	int		numberfiles = 0;
	cont=cont+1;//fix
	for (int folder = 0; folder < data.size(); folder++) {
		int         i = cont+folder;
		for (int it = 0; it < possib2pt.size(); it++) {
			string		folder2 = itos(folder);
			string		dosya3  (path2pt + folder2 + possib2pt[it]);
			opentheFile3.open(dosya3.c_str());
			if (opentheFile3.is_open()) {
				//cout << "Opened " << dosya3 << endl;
				for (int j = 0; j < Nt; j++) {
					double		skip;
					Double_t	val;
					opentheFile3 >> skip;
					opentheFile3 >> val;
					opentheFile3 >> skip;
					datamatrix2pt(j, i) += val*scale;
				}
			} else {
				cout << "File Error " << dosya3 << " could not open!" << endl;
				exit(1);
			}

			opentheFile3.close();
		}
	}
	return datamatrix2pt;
}


/// Reader functions
/** Read the data and calcualte average  */
TMatrixD
trQCD::readerAveragedC(int beginf, int endf, int increment, int psqr, vector < string > &possib2pt, int divider, string path2pt,TMatrixD & datamatrix2pt,int cont)
{
	int nconf=(endf-beginf)/increment;
	nconf++;
	ifstream	opentheFile3;
	int		nbcol = datamatrix2pt.GetNcols();
	double		scale = 1. / divider;
	
	//cout << "Path is: " << path2pt << endl;
	cout << "Now looking for files then calculate the average..." << endl;
	int		numberfiles = 0;
	cont=cont+1;//fix
	for (int folder = beginf; folder <= endf; folder = folder + increment) {
		int         i = cont+(folder - beginf) / increment;
		for (int it = 0; it < possib2pt.size(); it++) {
			string		folder2 = itos(folder);
			string		dosya3  (path2pt + folder2 + possib2pt[it]);
			opentheFile3.open(dosya3.c_str());
			if (opentheFile3.is_open()) {
				//cout << "Opened " << dosya3 << endl;
				for (int j = 0; j < Nt; j++) {
					double		skip;
					Double_t	val;
					opentheFile3 >> skip;
					opentheFile3 >> val;
					opentheFile3 >> skip;
					datamatrix2pt(j, i) += val*scale;
				}
			} else {
				cout << "File Error " << dosya3 << " could not open!" << endl;
				exit(1);
			}

			opentheFile3.close();
		}
	}
	return datamatrix2pt;
}


/// Magnetic Form Fac Calculation
/** Magnetic Form Fac Calculation using vector 1310.5915 eq:13  */
TMatrixD
trQCD::MagneticFormfac(vector <int> data, int snk, int src, int sink , int q2, vector <Double_t> constant, string name3pt, string path3pt, TMatrixD datamat2pt, TMatrixD datamat2pt_0, TMatrixD datamat2pt_shwl,bool diag)
{
	int nbfiles = data.size();
    int		g[4] = {8, 1, 2, 4};
    int		divider = 0;
    double		ort_deltaE_Jack[Nt];
    TMatrixD	Tot_Ratio(Nt, nbfiles);
    for (int j = 1; j < 4; j++) {
    	for (int i = 1; i < 4; i++) {
            for (int qz = 3; qz > -4; qz=qz-1) {
		        for (int qy = 3; qy > -4; qy=qy-1) {
					for (int qx = 3; qx > -4; qx=qx-1) {
			            if (((qx * qx) + (qy * qy) + (qz * qz)) == q2) {
							if ((qx!=0) && (levi_civita(i,j,1)!=0)) {
                                divider+=1;
                                TMatrixD	Ratio(Nt, nbfiles);
                                int		levi = levi_civita(i,j,1);
                                std::ostringstream oss;
                                oss << name3pt << "_p" << j << "_snk" << snk << "_g" << g[i] << "_src" << src << "_qx" << qx << "_qy" << qy << "_qz" << qz;
                                std::string data_3pt = oss.str();
                                Ratio = reader(data, data_3pt.c_str(), path3pt);
                                Ratio = (1.0/(levi*qx)) * Ratio;
                                Tot_Ratio = Tot_Ratio + Ratio;
                                Ratio.Clear();
                            } //if qx
                            if ((qy != 0) && (levi_civita(i, j, 2) != 0) ) {
                                divider+=1;
                                TMatrixD	Ratio(Nt, nbfiles);
                                int		levi = levi_civita(i,j,2);
                                std::ostringstream oss;
                                oss << name3pt << "_p" << j << "_snk" << snk << "_g" << g[i] << "_src" << src << "_qx" << qx << "_qy" << qy << "_qz" << qz;
                                std::string data_3pt = oss.str();
                                Ratio = reader(data, data_3pt.c_str(), path3pt);
                                Ratio = (1.0/(levi*qy)) * Ratio;
                                Tot_Ratio = Tot_Ratio + Ratio;
                                Ratio.Clear();
                            } //if qy
                            if ((qz != 0) && (levi_civita(i,j,3) != 0)) {
                                divider+=1;
                                TMatrixD	Ratio(Nt, nbfiles);
                                int		levi = levi_civita(i,j,3);
                                std::ostringstream oss;
                                oss << name3pt << "_p" << j << "_snk" << snk << "_g" << g[i] << "_src" << src << "_qx" << qx << "_qy" << qy << "_qz" << qz;
                                std::string data_3pt = oss.str();
                                Ratio = reader(data, data_3pt.c_str(), path3pt);
                                Ratio = (1.0/(levi*qz)) * Ratio;
                                Tot_Ratio = Tot_Ratio + Ratio;
                                Ratio.Clear();
                            } //if qz
                        } //j
                    } //i
                } //if q2
            } //qx
        } //qy
    } //qz
    double		scale = 1. / divider;
    Tot_Ratio = scale * Tot_Ratio;
    TMatrixD	Tot_Ratio_scale(Nt, nbfiles);
    TMatrixD	Tot_Ratio_jack(Nt, nbfiles);
    for (int col = 0; col < nbfiles; col++) {
        for(int row = 0; row < Nt; row++) {
            Tot_Ratio_scale (row,col) =  Tot_Ratio(row,col) * constant[col];
        }
    }
    Tot_Ratio_jack = JackKnife(Tot_Ratio_scale);
    return Ratiohalf(datamat2pt_0, datamat2pt, datamat2pt_shwl, Tot_Ratio_jack, sink, diag);
}

/// Pauli Form Factor Calculation
/** Pauli form factor calculation using arxiv:1603.04762 eq(12) */
vector <Double_t> 
trQCD::PauliFormFac(vector < vector <Double_t> > G_E, vector < vector <Double_t> > G_M, vector <Double_t> Qsqr, vector <Double_t>  par1_GE_err,vector <Double_t>  par1_GM_err, vector <Double_t> M1,  vector <Double_t> M2,int select) 
{
	vector <Double_t> F2;
	vector <Double_t> f2fit;
	vector <Double_t> f2err;
	Qsqr.erase(Qsqr.begin());
		for(int i=0; i<G_E[1].size(); i++){ //conf
			for(int k=0; k<G_M.size(); k++){//qsqr 1-8
			  Double_t Msqr = pow((inva*(M1[k]+M2[k])),2.0);
		  	  f2fit.push_back(((Msqr/(Msqr-Qsqr[k]))*(G_M[k][i]-G_E[k+1][i])));
		  	  f2err.push_back(TMath::Sqrt(pow(par1_GE_err[k+1],2.0)+pow(par1_GM_err[k], 2.0)));
	    	}// vectors are defined
			TF1 *f1;
			switch(select)
			{
				case 0 : //Monopole form
					f1 = new TF1("f1","[0]/(1+x/[1]^2)",0.0, Qsqr[6]);
				break;
				case 1 : //Dipole form
			       	 f1 = new TF1("f1","[0]/(1+x/[1]^2)^2",0.0, Qsqr[6]);
		    	break;
		    	case 2 : //Exponential form
			       	 f1 = new TF1("f1","[0]*TMath::Exp(-x/[1]^2)",0.0, Qsqr[6]);
		    	break;
			}
		    TRandom* random1 = new TRandom;
		    double par1 = random1->Rndm();
		    double par2 = random1->Rndm();
			TGraphErrors *grav = new TGraphErrors(Qsqr.size(),&(Qsqr[0]),&(f2fit[0]),0,&(f2err[0]));
		    f1->SetParameter(0,par1);
		    f1->SetParameter(1,par2);
			grav->Fit(f1, "RQN+","S");
			grav->Fit(f1, "RQN+","S");
			grav->Fit(f1, "RQN+","S");
			grav->Fit(f1, "RQN+","S");
			grav->Fit(f1, "RQN+","S");
			F2.push_back(f1->GetParameter(0));
			f1->SetLineColor(2);
			grav->SetMarkerColor(4);
			grav->SetMarkerStyle(21);
			grav->GetXaxis()->SetTitle("Q^{2} [GeV^{2}]");
			grav->GetYaxis()->SetTitle("F_{2}");
			//TCanvas * Plat = new TCanvas("f2", "f2");
		    grav->Draw("ALP");
			f1->Draw("SAME");
			//Plat->Write();
			f2fit.clear();
			f2err.clear();
			grav->Clear();
			//Plat->Clear();
		}//conf ends;
	return F2;
}

/// Pauli Form Factor Calculation
/** Pauli form factor calculation using arxiv:1603.04762 eq(12) */
vector <Double_t> 
trQCD::PauliFormFac(vector < vector <Double_t> > G_E, vector < vector <Double_t> > G_M, vector <Double_t> Qsqr, vector <Double_t>  par1_GE_err,vector <Double_t>  par1_GM_err,Double_t M1, Double_t M2,int select) 
{
	vector <Double_t> F2;
	Double_t Msqr = (M1 + M2) * (M1+ M2);
	vector <Double_t> f2fit;
	vector <Double_t> f2err;
	Qsqr.erase(Qsqr.begin());
		for(int i=0; i<G_E[1].size(); i++){ //conf
			for(int k=0; k<G_M.size(); k++){//qsqr 1-8
		  	  f2fit.push_back(((Msqr/(Msqr-Qsqr[k]))*(G_M[k][i]-G_E[k+1][i])));
		  	  f2err.push_back(TMath::Sqrt(pow(par1_GE_err[k+1],2.0)+pow(par1_GM_err[k], 2.0)));
	    	}// vectors are defined
			TF1 *f1;
			switch(select)
			{
				case 0 : //Monopole form
					f1 = new TF1("f1","[0]/(1+x/[1]^2)",0.0, Qsqr[6]);
				break;
				case 1 : //Dipole form
			       	 f1 = new TF1("f1","[0]/(1+x/[1]^2)^2",0.0, Qsqr[6]);
		    	break;
		    	case 2 : //Exponential form
			       	 f1 = new TF1("f1","[0]*TMath::Exp(-x/[1]^2)",0.0, Qsqr[6]);
		    	break;
			}
		    TRandom* random1 = new TRandom;
		    double par1 = random1->Rndm();
		    double par2 = random1->Rndm();
			TGraphErrors *grav = new TGraphErrors(Qsqr.size(),&(Qsqr[0]),&(f2fit[0]),0,&(f2err[0]));
		    f1->SetParameter(0,par1);
		    f1->SetParameter(1,par2);
			grav->Fit(f1, "RQN+","S");
			grav->Fit(f1, "RQN+","S");
			grav->Fit(f1, "RQN+","S");
			grav->Fit(f1, "RQN+","S");
			grav->Fit(f1, "RQN+","S");
			F2.push_back(f1->GetParameter(0));
			f1->SetLineColor(2);
			grav->SetMarkerColor(4);
			grav->SetMarkerStyle(21);
			grav->GetXaxis()->SetTitle("Q^{2} [GeV^{2}]");
			grav->GetYaxis()->SetTitle("F_{2}");
			//TCanvas * Plat = new TCanvas("f2", "f2");
		    grav->Draw("ALP");
			f1->Draw("SAME");
			//Plat->Write();
			f2fit.clear();
			f2err.clear();
			grav->Clear();
			//Plat->Clear();
		}//conf ends;
	return F2;
}

/// Calculation of magnetic moment 
/** Calculation of magnetic moment in nuclear magneton*/
vector <Double_t> 
	trQCD::mu_N(vector <Double_t> GM0, vector <Double_t> mass1,vector <Double_t> mass2){
		vector <Double_t> muN;
		for(int i=0; i<GM0.size(); i++){
			muN.push_back(GM0[i]*(0.938272/(((mass1[i]+mass2[i])/2)*inva)));
		}
		return muN;
	}

/// Decayrate
/** Decayrate calculation using B1 -> B2 gamma KeV*/
vector <Double_t> 
	trQCD::DecayRate(vector <Double_t> F2, Double_t M1, Double_t M2){
		
		Double_t q = TMath::Abs((TMath::Power(M1,2)-TMath::Power(M2,2))/(2*M1));
		Double_t Mm= (M1+M2);
		Double_t constant=(4*TMath::Power(q,3))/(137*TMath::Power(Mm,2))*TMath::Power(10,6);
		vector <Double_t> Decay;
		for(int i=0; i<F2.size(); i++){
			Decay.push_back(constant*F2[i]*F2[i]);
		}
		return Decay;
}
/// Lifetime
/** Lifetime calculation hbar/Decay */
vector <Double_t> 
	trQCD::Lifetime(vector <Double_t> Decay){
		vector <Double_t> Life;
		for(int i=0; i<Decay.size(); i++){
			Life.push_back(6.582119514*TMath::Power(10,-19)/Decay[i]);
		}
		return Life;
}
	


/// Overloaded Magnetic Form Fac Calculation supper
/** Magnetic Form Fac Calculation returns the raw matrix, user have to jackknife it 1310.5915 eq:13  */
TMatrixD
trQCD::MagneticFormfac(vector <vector <int> > data, int snk, int src, int sink , int q2, string name3pt, vector <string> path3pt)
{
	int nbfiles=0;
	for(int vcount=0; vcount<data.size(); vcount++){
		nbfiles += data[vcount].size();
	}
    int		g[4] = {8, 1, 2, 4};
    int		divider = 0;
    double		ort_deltaE_Jack[Nt];
    TMatrixD	Tot_Ratio(Nt, nbfiles);
//cout<<"{";
    for (int j = 1; j < 4; j++) {
    	for (int i = 1; i < 4; i++) {
            for (int qz = 3; qz > -4; qz--) {
		        for (int qy = 3; qy > -4; qy--) {
					for (int qx = 3; qx > -4; qx--) {
			            if (((qx * qx) + (qy * qy) + (qz * qz)) == q2) {
							if ((qx!=0) && (levi_civita(i,j,1)!=0)) {
                                divider+=1;
                                TMatrixD	Ratio(Nt, nbfiles);
                                int		levi = levi_civita(i,j,1);
                                std::ostringstream oss;
                                oss << name3pt << "_p" << j << "_snk" << snk << "_g" << g[i] << "_src" << src << "_qx" << qx << "_qy" << qy << "_qz" << qz;

                                std::string data_3pt = oss.str();
//cout<<"{"<<j<<","<<g[i]<<","<<qx<<","<<qy<<","<<qz<<"},";
				//cout<<data_3pt<<endl;
                                Ratio = reader(data, data_3pt.c_str(), path3pt);
                                Ratio = (1.0/(levi*qx)) * Ratio;
                                Tot_Ratio = Tot_Ratio + Ratio;
                                Ratio.Clear();
                            } //if qx
                            if ((qy != 0) && (levi_civita(i, j, 2) != 0) ) {
                                divider+=1;
                                TMatrixD	Ratio(Nt, nbfiles);
                                int		levi = levi_civita(i,j,2);
                                std::ostringstream oss;
                                oss << name3pt << "_p" << j << "_snk" << snk << "_g" << g[i] << "_src" << src << "_qx" << qx << "_qy" << qy << "_qz" << qz;
//cout<<"{"<<j<<","<<g[i]<<","<<qx<<","<<qy<<","<<qz<<"},";
                                std::string data_3pt = oss.str();
				//cout<<data_3pt<<endl;
                                Ratio = reader(data, data_3pt.c_str(), path3pt);
                                Ratio = (1.0/(levi*qy)) * Ratio;
                                Tot_Ratio = Tot_Ratio + Ratio;
                                Ratio.Clear();
                            } //if qy
                            if ((qz != 0) && (levi_civita(i,j,3) != 0)) {
                                divider+=1;
                                TMatrixD	Ratio(Nt, nbfiles);
                                int		levi = levi_civita(i,j,3);
                                std::ostringstream oss;
                                oss << name3pt << "_p" << j << "_snk" << snk << "_g" << g[i] << "_src" << src << "_qx" << qx << "_qy" << qy << "_qz" << qz;
                                std::string data_3pt = oss.str();
				//cout<<data_3pt<<endl;
//cout<<"{"<<j<<","<<g[i]<<","<<qx<<","<<qy<<","<<qz<<"},";
                                Ratio = reader(data, data_3pt.c_str(), path3pt);
                                Ratio = (1.0/(levi*qz)) * Ratio;
                                Tot_Ratio = Tot_Ratio + Ratio;
                                Ratio.Clear();
                            } //if qz
                        } //j
                    } //i
                } //if q2
            } //qx
        } //qy
    } //qz
//cout<<"}"<<endl;
    double		scale = 1.0 / divider;
    Tot_Ratio = scale * Tot_Ratio;
    return Tot_Ratio;
}


/// Overloaded Magnetic Form Fac Calculation
/** Magnetic Form Fac Calculation returns the raw matrix, user have to jackknife it 1310.5915 eq:13  */
TMatrixD
trQCD::MagneticFormfac(vector <int> data, int snk, int src, int sink , int q2, string name3pt, string path3pt)
{
	int nbfiles = data.size();
    int		g[4] = {8, 1, 2, 4};
    int		divider = 0;
    double		ort_deltaE_Jack[Nt];
    TMatrixD	Tot_Ratio(Nt, nbfiles);

    for (int j = 1; j < 4; j++) {
    	for (int i = 1; i < 4; i++) {
            for (int qz = 3; qz > -4; qz=qz-1) {
		        for (int qy = 3; qy > -4; qy=qy-1) {
					for (int qx = 3; qx > -4; qx=qx-1) {
			            if (((qx * qx) + (qy * qy) + (qz * qz)) == q2) {
							if ((qx!=0) && (levi_civita(i,j,1)!=0)) {
                                divider+=1;
                                TMatrixD	Ratio(Nt, nbfiles);
                                int		levi = levi_civita(i,j,1);
                                std::ostringstream oss;
                                oss << name3pt << "_p" << j << "_snk" << snk << "_g" << g[i] << "_src" << src << "_qx" << qx << "_qy" << qy << "_qz" << qz;

                                std::string data_3pt = oss.str();
				//cout<<data_3pt<<endl;
                                Ratio = reader(data, data_3pt.c_str(), path3pt);
                                Ratio = (1.0/(levi*qx)) * Ratio;
                                Tot_Ratio = Tot_Ratio + Ratio;
                                Ratio.Clear();
                            } //if qx
                            if ((qy != 0) && (levi_civita(i, j, 2) != 0) ) {
                                divider+=1;
                                TMatrixD	Ratio(Nt, nbfiles);
                                int		levi = levi_civita(i,j,2);
                                std::ostringstream oss;
                                oss << name3pt << "_p" << j << "_snk" << snk << "_g" << g[i] << "_src" << src << "_qx" << qx << "_qy" << qy << "_qz" << qz;

                                std::string data_3pt = oss.str();
				//cout<<data_3pt<<endl;
                                Ratio = reader(data, data_3pt.c_str(), path3pt);
                                Ratio = (1.0/(levi*qy)) * Ratio;
                                Tot_Ratio = Tot_Ratio + Ratio;
                                Ratio.Clear();
                            } //if qy
                            if ((qz != 0) && (levi_civita(i,j,3) != 0)) {
                                divider+=1;
                                TMatrixD	Ratio(Nt, nbfiles);
                                int		levi = levi_civita(i,j,3);
                                std::ostringstream oss;
                                oss << name3pt << "_p" << j << "_snk" << snk << "_g" << g[i] << "_src" << src << "_qx" << qx << "_qy" << qy << "_qz" << qz;
                                std::string data_3pt = oss.str();
				//cout<<data_3pt<<endl;
                                Ratio = reader(data, data_3pt.c_str(), path3pt);
                                Ratio = (1.0/(levi*qz)) * Ratio;
                                Tot_Ratio = Tot_Ratio + Ratio;
                                Ratio.Clear();
                            } //if qz
                        } //j
                    } //i
                } //if q2
            } //qx
        } //qy
    } //qz
    double		scale = 1.0 / divider;
    Tot_Ratio = scale * Tot_Ratio;
    return Tot_Ratio;
}

Double_t
trQCD::ReNormVec(Double_t kappaqrk_Real, Double_t kappacrit_Real)
{
	Double_t c1=-0.331;
	Double_t Beta=1.9;
	Double_t Nf=3.0;
	Double_t a=0.0907;
	Double_t c0=1.0-8.0*c1;
	Double_t g=6/Beta;
	Double_t P=1.0-0.1402*g;
	Double_t R=1.0-0.2689*g;
	Double_t Mu=1/a;
	Double_t gMSinv=(((c0*P)+(8*c1*R))/g)-0.1006+(0.03149*Nf)+((11-(2*Nf)/3)*TMath::Log(Mu*a))/(8*(TMath::Pi())*(TMath::Pi()));
	Double_t gMS= 1.0/gMSinv;
	Double_t u0=TMath::Power((1-(0.8412/Beta)),0.25);
	Double_t ZV=1.0-0.0277*gMS;
	Double_t bV=1.0+0.0382*gMS;
	Double_t m=0.5*(1.0/kappaqrk_Real-1/kappacrit_Real);
	Double_t rAxC= u0*ZV*(1.0+(bV*m)/u0);
	return rAxC;
}




/// Overloaded Magnetic Form Fac Calculation
/** Magnetic Form Fac Calculation 1310.5915 eq:13  */
TMatrixD
trQCD::MagneticFormfac(int beginfolder, int endfolder, int increment, int snk, int src, int sink , int q2, vector <Double_t> constant, string name3pt, string path3pt, TMatrixD datamat2pt, TMatrixD datamat2pt_0, TMatrixD datamat2pt_shwl,bool diag)
{
    int		nbfiles = (endfolder - beginfolder) / increment + 1;
    int		g[4] = {8, 1, 2, 4};
    int		divider = 0;
    double		ort_deltaE_Jack[Nt];
    TMatrixD	Tot_Ratio(Nt, nbfiles);
    for (int qx = -3; qx < 4; qx++) {
        for (int qy = -3; qy < 4; qy++) {
            for (int qz = -3; qz < 4; qz++) {
                if (((qx * qx) + (qy * qy) + (qz * qz)) == q2) {
                    for (int i = 1; i < 4; i++) {
                        for (int j = 1; j < 4; j++) {
                            if ((qx!=0) && (levi_civita(i,j,1)!=0)) {
                                divider+=1;
                                TMatrixD	Ratio(Nt, nbfiles);
                                int		levi = levi_civita(i, j, 1);
                                std::ostringstream oss;
                                oss << name3pt << "_p" << j << "_snk" << snk << "_g" << g[i] << "_src" << src << "_qx" << qx << "_qy" << qy << "_qz" << qz;
                                std::string data_3pt = oss.str();
                                Ratio = reader(beginfolder, endfolder, increment, data_3pt.c_str(), path3pt);
                                Ratio = (1.0/(levi*qx)) * Ratio;
                                Tot_Ratio = Tot_Ratio + Ratio;
								//Ratio.Print();
                                Ratio.Clear();
                            } //if qx
                            if ((qy!=0) && (levi_civita(i,j,2)!=0)) {
                                divider+=1;
                                TMatrixD	Ratio(Nt, nbfiles);
                                int		levi = levi_civita(i, j, 2);
                                std::ostringstream oss;
                                oss << name3pt << "_p" << j << "_snk" << snk << "_g" << g[i] << "_src" << src << "_qx" << qx << "_qy" << qy << "_qz" << qz;
                                std::string data_3pt = oss.str();
                                Ratio = reader(beginfolder, endfolder, increment, data_3pt.c_str(), path3pt);
                                Ratio = (1.0/(levi*qy)) * Ratio;
                                Tot_Ratio = Tot_Ratio + Ratio;
                                Ratio.Clear();
                            } //if qy
                            if ((qz!=0) && (levi_civita(i,j,3)!=0)) {
                                divider+=1;
                                TMatrixD	Ratio(Nt, nbfiles);
                                int		levi = levi_civita(i, j, 3);
                                std::ostringstream oss;
                                oss << name3pt << "_p" << j << "_snk" << snk << "_g" << g[i] << "_src" << src << "_qx" << qx << "_qy" << qy << "_qz" << qz;
                                std::string data_3pt = oss.str();
                                Ratio = reader(beginfolder, endfolder, increment, data_3pt.c_str(), path3pt);
                                Ratio = (1.0/(levi*qz)) * Ratio;
                                Tot_Ratio = Tot_Ratio + Ratio;
                                Ratio.Clear();
                            } //if qz
                        } //j
                    } //i
                } //if q2
            } //qx
        } //qy
    } //qz
    double		scale = 1. / divider;
    Tot_Ratio = scale * Tot_Ratio;
    TMatrixD	Tot_Ratio_scale(Nt, nbfiles);
    TMatrixD	Tot_Ratio_jack(Nt, nbfiles);
    for (int col = 0; col < nbfiles; col++) {
        for(int row = 0; row < Nt; row++) {
            Tot_Ratio_scale (row,col) =  Tot_Ratio(row,col) * constant[col];
        }
    }
    Tot_Ratio_jack = JackKnife(Tot_Ratio_scale);
    return Ratiohalf(datamat2pt_0, datamat2pt, datamat2pt_shwl, Tot_Ratio_jack, sink, diag);
}


/// Combine two matrices.
/** Combine two matricex x1 and x2 and return as thrid matrix which size is x1.size() + x2.size() . */
TMatrixD
trQCD::CombineMatrices(TMatrixD x1, TMatrixD x2) {
	const Int_t nrow1 = x1.GetNrows();
	const Int_t nrow2 = x2.GetNrows();
	const Int_t ncol1 = x1.GetNcols();
	const Int_t ncol2 = x2.GetNcols();

	const Int_t nrow3 = nrow1;
	const Int_t ncol3 = ncol1+ncol2;

	TMatrixD x3(nrow3,ncol3);
	TMatrixDSub(x3,0,nrow1-1,0,ncol1-1) = x1;
	TMatrixDSub(x3,0,nrow1-1,ncol1,ncol3-1) = x2;
	return x3;
}


/// Combine three matrices.
/** Combine three matricex x1 and x2 and return as thrid matrix which size is x1.size() + x2.size() . */
TMatrixD
trQCD::CombineMatrices(TMatrixD x1, TMatrixD x2, TMatrixD x3) {
	const Int_t nrow1 = x1.GetNrows();
	const Int_t nrow2 = x2.GetNrows();
	const Int_t nrow3 = x3.GetNrows();
	
	const Int_t ncol1 = x1.GetNcols();
	const Int_t ncol2 = x2.GetNcols();
	const Int_t ncol3 = x3.GetNcols();
	

	const Int_t nrow4 = nrow1;
	const Int_t tempcol_1 = ncol1+ncol2;
	const Int_t new_col = ncol1+ncol2+ncol3;
	TMatrixD x4(nrow4,new_col);

	TMatrixDSub(x4,0,nrow1-1,0,ncol1-1) = x1;
	TMatrixDSub(x4,0,nrow1-1,ncol1,tempcol_1-1) = x2;
	TMatrixDSub(x4,0,nrow1-1,tempcol_1,new_col-1) = x3;


	return x4;
}

/// Combine four matrices.
/** Combine four matricex x1 and x2 and return as thrid matrix which size is x1.size() + x2.size() . */
TMatrixD
trQCD::CombineMatrices(TMatrixD x1, TMatrixD x2, TMatrixD x3, TMatrixD x4) {
	const Int_t nrow1 = x1.GetNrows();
	const Int_t nrow2 = x2.GetNrows();
	const Int_t nrow3 = x3.GetNrows();
	const Int_t nrow4 = x4.GetNrows();
	
	const Int_t ncol1 = x1.GetNcols();
	const Int_t ncol2 = x2.GetNcols();
	const Int_t ncol3 = x3.GetNcols();
	const Int_t ncol4 = x4.GetNcols();
	
	const Int_t tempcol_1 = ncol1+ncol2;
	const Int_t tempcol_2 = ncol1+ncol2+ncol3;
	const Int_t new_col = ncol1+ncol2+ncol3+ncol4;
	TMatrixD x5(nrow4,new_col);
	TMatrixDSub(x5,0,nrow1-1,0,ncol1-1) = x1;
	TMatrixDSub(x5,0,nrow1-1,ncol1,tempcol_1-1) = x2;
	TMatrixDSub(x5,0,nrow1-1,tempcol_1,tempcol_2-1) = x3;
	TMatrixDSub(x5,0,nrow1-1,tempcol_2,new_col-1) = x4;

	return x5;
}

/// 	form factor calculation for spin 1/2 particles
/** 1503.07361 eq:15. 
There is no scale here 																						*/
TMatrixD
trQCD::Ratio1503(TMatrixD datamat2pt_0, TMatrixD datamat2pt, TMatrixD datamat2pt_shwl, TMatrixD datamat3pt, int sink,bool diagonal)
{
    int			nbcol = datamat2pt.GetNcols();
    
    TMatrixD	datamatrix_Sum(Nt, nbcol);
    
    TMatrixD	datamatrix_division(sink+1, nbcol);
    
    datamatrix_Sum = (datamat3pt) * pow(Nl, 3.0);
    for (int row = 0; row <= sink; row++) {
	for (int col = 0; col < nbcol; col++) {
		int place=2*row;
		datamatrix_division(row,col) = (datamatrix_Sum(row, col) / datamat2pt_shwl(sink, col)) * TMath::Sqrt(TMath::Abs(datamat2pt_0(place, col)/datamat2pt(place, col)));
				}
			}
}
/// 	form factor calculation for spin 1/2 particles
/** The diagonal equation from 1310.5915 eq:12 (transition) equation from 1503.07361 eq:15. 
There is no scale here
*/
TMatrixD
trQCD::Ratiohalf (TMatrixD datamat2pt_0, TMatrixD datamat2pt, TMatrixD datamat2pt_shwl, TMatrixD datamat3pt, int sink,bool diagonal)
{
    int			nbcol = datamat2pt.GetNcols();
    
    TMatrixD	datamatrix_Sum(Nt, nbcol);
    
    TMatrixD	datamatrix_division(sink, nbcol);
    
    datamatrix_Sum = (datamat3pt) * pow(Nl, 3.0);
		for(int row = 0; row < sink; row++) {
		int place= sink-row;
	        for (int col = 0; col < nbcol; col++) {
			Double_t temp1 =(datamat2pt_0(row, col)/datamat2pt(row, col));
			Double_t temp2 =(datamat2pt(place, col)/datamat2pt_0(place, col));
			Double_t temp3 =(datamat2pt_0(sink, col)/datamat2pt(sink,col));
			datamatrix_division(row,col) = (datamatrix_Sum(row, col) / datamat2pt_shwl(sink, col))*(TMath::Sqrt(TMath::Abs(temp1*temp2*temp3 )));
	        }
	    }	

    return datamatrix_division;
}

/// 	form factor calculation for spin 1/2 particles
/** equation from 1503.07361 eq:15. 
There is scale vector for electric form factor calculation */
TMatrixD
trQCD::Ratio1503 (TMatrixD datamat2pt_0, TMatrixD datamat2pt, TMatrixD datamat2pt_shwl, TMatrixD datamat3pt, vector <Double_t> scale ,int sink, bool diagonal)
{
    int			nbcol = datamat2pt.GetNcols();
    
    TMatrixD	datamatrix_scale(Nt, nbcol);
    
    TMatrixD	datamatrix_Sum(Nt, nbcol);
    
    TMatrixD	datamatrix_division(sink+1, nbcol);
        for (int col = 0; col < nbcol; col++) {
        	for(int row = 0; row < Nt; row++) {
            		datamatrix_scale (row,col) =  datamat3pt(row,col) * scale[col];
        	}
    }
    datamatrix_Sum = datamatrix_scale * pow(Nl, 3.0);

		for (int row = 0; row <= sink; row++) {
				for (int col = 0; col < nbcol; col++) {
					int place=2*row;
		                        Double_t temp1 =TMath::Sqrt((datamat2pt_0(place, col)/datamat2pt(place, col)));
	                                Double_t temp4 =(datamatrix_Sum(row, col) / datamat2pt_shwl(sink, col));
					datamatrix_division(row,col) = temp4*temp1;
				}
			}
}

/// 	form factor calculation for spin 1/2 particles
/** The diagonal equation from 1310.5915 eq:12 the non-diagonal (transition) equation from 1503.07361 eq:15. 
There is scale vector for electric form factor calculation */
TMatrixD
trQCD::Ratiohalf (TMatrixD datamat2pt_0, TMatrixD datamat2pt, TMatrixD datamat2pt_shwl, TMatrixD datamat3pt, vector <Double_t> mass_1, int qsqr ,int sink, bool electric)
{
    int			nbcol = datamat2pt.GetNcols();
    
    TMatrixD	datamatrix_scale(sink, nbcol);
    
    TMatrixD	datamatrix_Sum(sink, nbcol);
    
    TMatrixD	datamatrix_division(sink, nbcol);

		for(int row = 0; row < sink; row++) {
		int place=sink-row;
	        	for (int col = 0; col < nbcol; col++) {
				Double_t temp1 =TMath::Sqrt((datamat2pt_0(row, col)/datamat2pt(row, col)));
				Double_t temp2 =TMath::Sqrt(datamat2pt(place, col)/datamat2pt_0(place, col));
				Double_t temp3 =TMath::Sqrt(datamat2pt_0(sink, col)/datamat2pt(sink,col));
				Double_t temp4 = datamat3pt(row, col) / datamat2pt_shwl(sink, col);
				datamatrix_division(row,col) = temp4*temp1*temp2*temp3;
	    		}	
		}
    Double_t Constant=0.0;
    for (int col = 0; col < nbcol; col++) {
                if(electric) {
                        //electric
                        Constant= ConstEFF(qsqr, mass_1[col]);  
                }    
                else {
                        //magnetic
                        Constant= ConstMFF(qsqr, mass_1[col]);
                }    
                for(int row = 0; row < sink; row++) {
                        datamatrix_scale (row,col) =  datamatrix_division(row,col) * Constant;
                }    
	Constant=0.0;
    } 
    datamatrix_Sum = datamatrix_scale * pow(Nl, 3.0);

    return datamatrix_Sum;
}


/// 	TGraph function
/** Because I am lazy */
void
trQCD::draw_TGraph(vector < Double_t > &data, vector < Double_t > &xarray, vector < Double_t > &stdDev, string grname, string x_name, string y_name)
{
    std::vector<Double_t> zero (data.size(),0);
    TGraphErrors *gr = new TGraphErrors(data.size(),&(xarray[0]),&(data[0]),&(zero[0]),&(stdDev[0]));
    gr->SetTitle(grname.c_str());
    gr->SetMarkerColor(4);
    gr->SetMarkerStyle(21);
    gr->Draw("ALP");
    gr->GetXaxis()->SetTitle(x_name.c_str());
    gr->GetYaxis()->SetTitle(y_name.c_str());
    gr->Write(grname.c_str());
}


/// 	Fit function
/** It fits the constant [0] */
double
trQCD::FitYap(vector < Double_t > &data, vector < Double_t > &xarray, vector < Double_t > &stdDev,int *aralik, string grname, bool write)
{
    int		fitbas = aralik[0];
    int		fitbit = aralik[1];
    //TCanvas * Plat = new TCanvas(grname.c_str(), grname.c_str());
    TRandom* random1 = new TRandom;
    TGraphErrors *gr = new TGraphErrors(data.size(),&(xarray[0]),&(data[0]),0,&(stdDev[0]));
    //TF1 * f1 = new TF1("f1", "[0]", fitbas, fitbit);
    double par1 = random1->Rndm();
    double par2 = random1->Rndm();
    TF1 * f1 = new TF1("f1", "([0]/(2*[1]))*TMath::Exp(-[1]*x)", fitbas, fitbit);
    f1->SetParameter(0,par1);
    f1->SetParameter(1,par2);
    gr->Fit(f1, "RQN+");
    gr->Fit(f1, "RQN+");
    gr->Fit(f1, "RQN+");
    gr->Fit(f1, "RQN+");
    gr->Fit(f1, "RQN+");
    gr->Fit(f1, "RQN+");
    gr->Fit(f1, "RQ");
    char		parameter [100];
    double		parz;
    parz = f1->GetParameter(1);
    sprintf(parameter, "p0=%e", parz);
    TLegend        *legend = new TLegend(0.70, 0.90, 0.99, 0.99);
    legend->SetTextFont(50);
    legend->SetTextSize(0.04);
    legend->AddEntry(f1, parameter, "l");
    legend->Draw();
    gr->SetMarkerColor(4);
    gr->SetMarkerStyle(21);
    gr->Draw("ALP");
    f1->SetLineColor(2);
    f1->Draw("SAME");
    gr->SetTitle(grname.c_str());
    gr->GetXaxis()->SetTitle("Timestep");
    if(write){
    	gr->Write(grname.c_str());
    	//Plat->Update();
    	//Plat->Write(grname.c_str());
	}
    
    return parz;
}

/// 	Fitter for mass calculation
/** It fits the constant [0] */
Double_t
trQCD::FitterLog(vector < Double_t > &data, vector < Double_t > &xarray, vector < Double_t > &stdDev,int *aralik)
{
    int		fitbas = aralik[0];
    int		fitbit = aralik[1];
    std::vector<Double_t> zero (data.size(),0);
    TGraphErrors *gr = new TGraphErrors(data.size(),&(xarray[0]),&(data[0]),&(zero[0]),&(stdDev[0]));
    TF1 * f1 = new TF1("f1", "[0]", fitbas, fitbit);
    gr->Fit(f1, "RQN+");
    gr->Fit(f1, "RQN+");
    gr->Fit(f1, "RQN+");
    gr->Fit(f1, "RQN+");
    gr->Fit(f1, "RQN+");
    gr->Fit(f1, "RQN+");
    char		parameter [100];
    double		parz;
    f1->Draw();
    parz = f1->GetParameter(0);
    sprintf(parameter, "p0=%e", parz);
    TLegend        *legend = new TLegend(0.70, 0.90, 0.99, 0.99);
    legend->SetTextFont(50);
    legend->SetTextSize(0.04);
    legend->AddEntry(f1, parameter, "l");
    legend->Draw();
    gr->SetMarkerColor(4);
    gr->SetMarkerStyle(21);
    gr->Draw("ALP");
    gr->GetXaxis()->SetTitle("Timestep");
    gr->GetYaxis()->SetTitle("Average of  #frac{C(n+1)}{C(n)}");
    //gr->Write("Log Mass-Fit");
    
    return parz;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
/// 	Fitter for mass calculation
/** It fits the constant [0] */
Double_t
trQCD::Fitter(vector < Double_t > &data, vector < Double_t > &xarray, vector < Double_t > &stdDev,int *aralik)
{
    int		fitbas = aralik[0];
    int		fitbit = aralik[1];
    std::vector<Double_t> zero (data.size(),0);
    TGraphErrors *gr = new TGraphErrors(data.size(),&(xarray[0]),&(data[0]),0,&(stdDev[0]));
    TRandom* random1 = new TRandom;
    double par1 = random1->Rndm();
    double par2 = random1->Rndm();
    TF1 * f1 = new TF1("f1", "([0]/(2*[1]))*TMath::Exp(-[1]*x)", fitbas, fitbit);
    f1->SetParameter(0,par1);
    f1->SetParameter(1,2);
    gr->Fit(f1, "RQN+");
    gr->Fit(f1, "RQN+");
    gr->Fit(f1, "RQN+");
    gr->Fit(f1, "RQN+");
    gr->Fit(f1, "RQN+");
    gr->Fit(f1, "RQ");
    char		parameter [100];
    double		parz;
    f1->Draw();
    parz = f1->GetParameter(1);
    sprintf(parameter, "p0=%e", parz);
    TLegend        *legend = new TLegend(0.70, 0.90, 0.99, 0.99);
    legend->SetTextFont(50);
    legend->SetTextSize(0.04);
    legend->AddEntry(f1, parameter, "l");
    legend->Draw();
    gr->SetMarkerColor(4);
    gr->SetMarkerStyle(21);
    gr->Draw("ALP");
    gr->GetXaxis()->SetTitle("Timestep");
    gr->GetYaxis()->SetTitle("Average of  #frac{C(n+1)}{C(n)}");
    //gr->Write("Mass-Fit");
    
    return parz;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////
/****************		Another Ratio						*************/
/***************  	It is not written well so I disabled it	*************/
//////////////////////////////////////////////////////////////////////////
void
trQCD::RijRatio(int snkp, int srcp,int current,int polarization,string Path3pt,string Path2pt,string quark,int beginfolder,int endfolder,int beginfolder2,int endfolder2,int sink,int qx,int qy,int qz,TMatrixD & datamatrix2_jk,TMatrixD & datamatrix2_shwl_0_jk,TMatrixD & Ratio)
{
/*	int			nbnum = (endfolder - beginfolder) / 20;
	int			nbnum2 = (endfolder2 - beginfolder2) / 20;
	int 		nbfiles = nbnum+nbnum2+2;
	TMatrixD	Ratio1(Nt, nbfiles);
	TMatrixD	Ratio2(Nt, nbfiles);
	/****************************************/
/*	TMatrixD	datamatrix2_nuc(Nt, nbfiles);
	TMatrixD	datamatrix2_nucm(Nt, nbfiles);
	/****************************************/
/*	TMatrixD	datamatrix3_shwl(Nt,nbfiles);
	TMatrixD	datamatrix3_shwlm(Nt,nbfiles);
	/****************************************/
/*	int			nbcol = datamatrix2_jk.GetNcols();
	char		names2pt_nuc[100];
	char		names2pt_nucm[100];
	char		names3pt_shwl_cst[100];
	char		names3pt_shwl_mst[100];

	sprintf(names2pt_nuc, "/DATA_Sigma_c_shsh_px%d_py%d_pz%d", qx,qy,qz);
	sprintf(names2pt_nucm, "/DATA_Sigma_c_shsh_px%d_py%d_pz%d", -qx,-qy,-qz);
	sprintf(names3pt_shwl_cst, "/nonlocal_cur3ptfn_wallDeltaPFormFacH_%s_p%d_snk%d_g%d_src%d_qx%d_qy%d_qz%d",quark.c_str(),polarization,snkp,current,srcp,qx,qy,qz);
	sprintf(names3pt_shwl_mst, "/nonlocal_cur3ptfn_wallDeltaPFormFacH_%s_p%d_snk%d_g%d_src%d_qx%d_qy%d_qz%d",quark.c_str(),polarization,snkp,current,srcp,-qx,-qy,-qz);
	/****************************************/
/*	reader(beginfolder, endfolder, names2pt_nuc, Path2pt, datamatrix2_nuc);
	reader(beginfolder, endfolder, names2pt_nucm, Path2pt, datamatrix2_nucm);

   	reader(beginfolder,endfolder, names3pt_shwl_cst, Path3pt, datamatrix3_shwl);
   	reader(beginfolder,endfolder, names3pt_shwl_mst, Path3pt, datamatrix3_shwlm);

	readerC(beginfolder2, endfolder2, nbnum, names2pt_nuc, Path2pt, datamatrix2_nuc);
	readerC(beginfolder2, endfolder2, nbnum, names2pt_nucm, Path2pt, datamatrix2_nucm);

   	readerC(beginfolder2, endfolder2, nbnum, names3pt_shwl_cst, Path3pt, datamatrix3_shwl);
   	readerC(beginfolder2, endfolder2, nbnum, names3pt_shwl_mst, Path3pt, datamatrix3_shwlm);

	/*************************************************/
/*	TMatrixD	datamatrix2_nuc_jk(Nt, nbfiles);
	TMatrixD	datamatrix2_nucm_jk(Nt, nbfiles);
	TMatrixD	datamatrix3_shwl_jk(Nt, nbfiles);
	TMatrixD	datamatrix3_shwlm_jk(Nt, nbfiles);
	/****************************************/
	/****************************************/
/*	JackKnife(datamatrix2_nuc, datamatrix2_nuc_jk);
	JackKnife(datamatrix3_shwl, datamatrix3_shwl_jk);
	JackKnife(datamatrix3_shwlm, datamatrix3_shwlm_jk);
	JackKnife(datamatrix2_nucm, datamatrix2_nucm_jk);
	/****************************************/
/*	for (int row = 0; row <= sink; row++) {
		for (int col = 0; col < nbcol; col++) {
			int place=2*row;
			Double_t Temp=datamatrix2_nuc_jk(place,col)/datamatrix2_jk(place,col);
			Double_t Temp2=datamatrix2_nucm_jk(place,col)/datamatrix2_jk(place,col);
			if(Temp<0){
				Temp=-Temp;
			}
			if(Temp2<0){
				Temp2=-Temp2;
			}
			Ratio1(row,col)=datamatrix3_shwl_jk(row,col)/(datamatrix2_shwl_0_jk(sink,col) * TMath::Sqrt(Temp));
			Ratio2(row,col)=datamatrix3_shwlm_jk(row,col)/(datamatrix2_shwl_0_jk(sink,col) * TMath::Sqrt(Temp2));

		}
	}
	Ratio1 = Ratio1 *TMath::Power(32,3);
	Ratio2 = Ratio2 *TMath::Power(32,3);
	Ratio = (Ratio1-Ratio2);
*/
}



/// 	Good old jackknife fella
/** Jack knife function */
TMatrixD
trQCD::JackKnife(TMatrixD datamat) {
    int		N_col = datamat.GetNcols();
    int		N_row = datamat.GetNrows();
	if(N_col<2) {
		cout<<"Error! Jackknife cannot be made"<<endl;
		return datamat;
	}
	else {
		TMatrixD outputMat(N_row,N_col);
	    Double_t	Sum = 0.0;
	    for             (int rows = 0; rows < N_row; rows++) {
	        Sum = GetSumRows (datamat, rows);
	        for (int colz = 0; colz < N_col; colz++) {
	            outputMat(rows, colz) = (Sum - datamat(rows, colz)) / (N_col - 1);
	        }
	    }
	    return outputMat;		
	}

}
/// Calculate q^2 of specific particle
/** Calculating -q^2=2*m *(m-E) */
TMatrixD
trQCD::calculateQsqr(vector <Double_t> mass_par1, vector <Double_t> mass_par2){
	double q = TMath::TwoPi()/Nl;
	TMatrixD calqsqr(9,mass_par1.size());
	for(int qsqr=0; qsqr<9;qsqr++) {
		if(qsqr==7)
			qsqr++;
		for (unsigned i=0; i<mass_par1.size(); i++){
			Double_t E_par1= Energy(qsqr,mass_par1[i]);
			Double_t temp_qs=(mass_par1.at(i)*mass_par1.at(i))+(mass_par2.at(i)*mass_par2.at(i)) - (2*mass_par2.at(i)*E_par1);	
			calqsqr(qsqr,i)=-inva*inva*temp_qs;
		}
	}
	return calqsqr;	
}

/// Calculate Energy of specific particle
/** Calculating sqrt(m^2 + q^2 )*/
Double_t
trQCD::Energy(int qsqr,Double_t mass_par1)
{
	Double_t E_par;
	Double_t q = TMath::TwoPi()/Nl;
	Double_t qq=q*q;
	E_par=TMath::Sqrt((mass_par1*mass_par1)+(qsqr*(qq)));
	return E_par;
}

/// Calculate q^2 of specific particle
/** Calculating -q^2=2*m *(m-E) */
Double_t
trQCD::qq2(int qsqr, Double_t mass)
{
	Double_t qqsqr=2*mass*(mass-Energy(qsqr,mass));
	return qqsqr;
}

/// Calculate Q^2 of specific particle
/** Calculating Q^2=-a^2 q^2 */
Double_t
trQCD::Q2(int qsqr, Double_t mass)
{
	Double_t QQ=-inva*inva*qq2(qsqr,mass);
}

/// Electric FormFactor Constant of Ratio
/** sqrt(2E/(E+m)) */
Double_t
trQCD::ConstEFF(int qsqr, Double_t mass)
{
	Double_t ConstE= TMath::Sqrt((2*Energy(qsqr,mass))/(Energy(qsqr,mass)+mass));
	return ConstE;
}

/// Electric FormFactor Constant of Ratio
/** 2*sqrt(E2*E1/(m2+E2)*(m1+E1)) */
Double_t
trQCD::ConstEFF(int qsqr, Double_t mass1, Double_t mass2)
{
	Double_t ConstE= 2*TMath::Sqrt((Energy(qsqr,mass2)*Energy(qsqr,mass1))/((mass2+Energy(qsqr,mass2))*(mass1+Energy(qsqr,mass1))));
	return ConstE;
}


/// Magnetic FormFactor Constant of Ratio
/** sqrt(2E(E+m))/q */
Double_t
trQCD::ConstMFF(int qsqr, Double_t mass)
{
	Double_t q = TMath::TwoPi()/Nl;
	Double_t ConstM = TMath::Sqrt((2*Energy(qsqr,mass)*(Energy(qsqr,mass)+mass)))/q;
	return ConstM;
}

/// Magnetic FormFactor Constant of Ratio
/** sqrt(2E(E+m))/q */
Double_t
trQCD::ConstMFF(int qsqr, Double_t mass1,Double_t mass2)
{
	Double_t q = TMath::TwoPi()/Nl;
	Double_t ConstM =(2/q)*TMath::Sqrt((Energy(qsqr,mass1)*Energy(qsqr,mass2)*(Energy(qsqr,mass1)+mass1))/(Energy(qsqr,mass2)+mass2));
	return ConstM;
}

/// 	returns rows of matrix
/** Returns rows of matrix */
std::vector<Double_t>
trQCD::GetRows(TMatrixD Mat, int rows)
{
	int		N_colz = Mat.GetNcols();
	std::vector <Double_t> rowz;
	for (int i = 0; i < N_colz; i++) {
		rowz.push_back(Mat(rows,i));
	}
	return rowz;
}


/// 	returns cols of matrix
/** Returns cols of matrix */
std::vector<Double_t>
trQCD::GetCols (TMatrixD Mat, int column){
    int		N_row = Mat.GetNrows();
    std::vector <Double_t> cols;
    for (int i = 0; i < N_row; i++) {
        cols.push_back(Mat(i, column));
    }
    return cols;
}

///     returns cols of matrix till sink point
/** Returns cols of matrix */
std::vector<Double_t>
trQCD::GetColsSink (TMatrixD Mat, int column,int sink){
    std::vector <Double_t> cols;
    for (int i = 0; i < sink; i++) {
        cols.push_back(Mat(i, column));
    }
    return cols;
}

/// 	returns sum of the rows of matrix
/** it is useful for jackknife */
Double_t
trQCD::GetSumRows(TMatrixD Mat, int rows)
{
	/**************************************************/
	int		N_colz = Mat.GetNcols();
	Double_t	SumofRow = 0.0;
	for (int i = 0; i < N_colz; i++) {
		SumofRow += Mat(rows,i);
	}
	return SumofRow;
	/**************************************************/
    
}


/// 	Calculate average to the sink point
/** Calculate average of every time point seperately over all configurations to the given time point */
std::vector<Double_t>
trQCD::average_sink  (TMatrixD M, int sink){
	/**************************************************/ 
    int		N_max = sink;
    int		N_cf_max = M.GetNcols();
    std::vector <Double_t> avg_S;
    for (int i = 0; i < N_max; i++) {
        double		averagetot = 0.0;
        for (int j = 0; j < N_cf_max; j++) {
            averagetot += M(i, j);
        }
        avg_S.push_back(averagetot/N_cf_max);
    }
    return avg_S;
    /**************************************************/

}

///     Standard deviation
/**  Calculate the standard deviation, returns the average and error */
Double_t
trQCD::Average (std::vector<Double_t> the_Mat)
{
    /**************************************************/
    int          N_col = the_Mat.size();
    Double_t    tot = 0.0;
    Double_t    average=0.0;
    Double_t fix = N_col-1;
        fix = fix/N_col;
    for (unsigned i=0; i<the_Mat.size(); i++){
        average += the_Mat.at(i);
        }
        average = average/N_col;
    return average;
    /**************************************************/
}
///     Standard deviation
/**  Calculate the standard deviation, returns the average and error */
Double_t
trQCD::StdDev (std::vector<Double_t> the_Mat)
{
    /**************************************************/
    int          N_col = the_Mat.size();
    Double_t    tot = 0.0;
    Double_t    average=0.0;
    Double_t fix = N_col-1;
        fix = fix/N_col;
    for (unsigned i=0; i<the_Mat.size(); i++){
        average += the_Mat.at(i);
        }
        average = average/N_col;
    for (unsigned j=0; j<the_Mat.size(); j++) {
        tot += TMath::Power((the_Mat.at(j) - average), 2);
    }
        Double_t error = TMath::Sqrt(tot*fix);
    return error;
    /**************************************************/
}

/// 	Standard deviation
/**  Calculate the standard deviation, returns the average and error */
Double_t
trQCD::StdDev (std::vector<Double_t> the_Mat, Double_t & error)
{
    /**************************************************/
    int		 N_col = the_Mat.size();
    Double_t	tot = 0.0;
    Double_t 	average=0.0;
    Double_t fix = N_col-1;
	fix = fix/N_col;
    for (unsigned i=0; i<the_Mat.size(); i++){
        average += the_Mat.at(i);
	}
	average = average/N_col;
    for (unsigned j=0; j<the_Mat.size(); j++) {
        tot += TMath::Power((the_Mat.at(j) - average), 2);
    }
	error = TMath::Sqrt(tot*fix);
    return average;
    /**************************************************/
}


/// 	Standard deviation
/**  Calculate the standard deviation, returns the average and error */
std::vector<Double_t>
trQCD::StdDev (TMatrixD the_Mat, vector <Double_t> & the_Mat_ort, int sink){
    /**************************************************/
    int		N_col = the_Mat.GetNcols();
    int		N_row = the_Mat.GetNrows();
    std::vector<Double_t> stdSapma;
    Double_t	tot = 0.0;
    Double_t	fix = 0.0;
    fix =N_col - 1;
    fix =(fix) / N_col;
    the_Mat_ort	= average_sink(the_Mat, sink);
    for (int i = 0; i < sink; i++) {
        tot = 0.0;
        for (int j = 0; j < N_col; j++) {
            tot += (TMath::Power((the_Mat(i,j) - the_Mat_ort.at(i)),2));
        }
        stdSapma.push_back(TMath::Sqrt(fix*tot));
    }
    return stdSapma;
    /**************************************************/
}

/// 	Calculate average of Qsquare
/**  Calculate gm vs qsquare */
void
trQCD::central_fit (vector <Double_t>  average, vector <Double_t>  aver_err, vector <Double_t> qsquare, int select, string graphname)
{
	int boyut= qsquare.size()-1;
	
	TRandom* random1 = new TRandom;
	TF1 *f1, *f2;
	qsquare.erase(qsquare.begin());
	switch(select)
	{
		case 0 : //Monopole form
			
			f1 = new TF1("f1","[0]/(1+x/[1]^2)",qsquare[0], qsquare[boyut]);
	    		f2= new TF1("f2","[0]/(1+x/[1]^2)",-0.1,qsquare[boyut]+0.1);
		break;
		case 1 : //Dipole form
       	 		f1 = new TF1("f1","[0]/(1+x/[1]^2)^2",qsquare[0], qsquare[boyut]);
       	 		f2 = new TF1("f2","[0]/(1+x/[1]^2)^2",-0.1, qsquare[boyut]+0.1);
		break;
		case 2 : //Exponential form
       	 		f1 = new TF1("f1","[0]*TMath::Exp(-x/[1]^2)",qsquare[0], qsquare[boyut]);
       	 		f2 = new TF1("f2","[0]*TMath::Exp(-x/[1]^2)",-0.1, qsquare[boyut]+0.1);
		break;
	}
	double par1 = random1->Rndm();
	double par2 = random1->Uniform(0,5);
	f1->SetParameters(par1, par2);
	TGraphErrors *grav = new TGraphErrors(average.size(),&(qsquare[0]),&(average[0]),0,&(aver_err[0]));
	grav->Fit(f1, "RQN+","S");
	grav->Fit(f1, "RQN+","S");
	grav->Fit(f1, "RQN+","S");
	grav->Fit(f1, "RQN+","S");
	grav->Fit(f1, "RQN+","S");
	grav->SetMarkerColor(4);
	grav->SetMarkerStyle(21);
	grav->GetXaxis()->SetTitle("Q^2[GeV^-1]");
	grav->GetYaxis()->SetTitle("Average");
	grav->Draw("ALP");
	//TCanvas * Plat = new TCanvas(graphname.c_str(), graphname.c_str());
	f2->SetNpx(300);
	f2->SetParameters(f1->GetParameter(0), f1->GetParameter(1));
	f2->SetTitle(graphname.c_str());
	f2->SetLineColor(2);
	f2->SetLineStyle(4);
	//f2->Draw("SAME");
	//TCanvas * Plat2 = new TCanvas("dummy", "dummy");
        f2->Draw();
	grav->Draw("SAME P");
	grav->Write(graphname.c_str());
	//Plat->Write(graphname.c_str());
	//Plat->Clear();
	//Plat2->Write("dummmy");
}


/// 	Calculate average of Qsquare
/**  Calculate ge vs qsquare */
void
trQCD::central_fit (vector <Double_t>  average, vector <Double_t>  aver_err, vector <Double_t> qsquare, double charge ,int select, string graphname)
{
	int boyut= qsquare.size()-1;
	
    /**************************************************/
	TRandom* random1 = new TRandom;
	TF1 *f1;
	switch(select)
	{
		case 0 : //Monopole form
		f1 = new TF1("f1","[0]/(1+x/[1]^2)",0.0, qsquare[boyut]);
		break;
		case 1 : //Dipole form
       	 	f1 = new TF1("f1","[0]/(1+x/[1]^2)^2",0.0, qsquare[boyut]);
    		break;
    		case 2 : //Exponential form
       	 	f1 = new TF1("f1","[0]*TMath::Exp(-x/[1]^2)",0.0, qsquare[boyut]);
    		break;
	}
   	double par1 = charge;
	double par2= random1->Uniform(0,5);
	f1->SetParameters(par1, par2);
	TGraphErrors *grav = new TGraphErrors(average.size(),&(qsquare[0]),&(average[0]),0,&(aver_err[0]));
	f1->FixParameter(0,par1);
	grav->Fit(f1, "RQN+","S");
	grav->SetMarkerColor(4);
	grav->SetMarkerStyle(21);
	grav->SetTitle(graphname.c_str());
	grav->GetXaxis()->SetTitle("Q^2[GeV^-1]");
	grav->GetYaxis()->SetTitle("Average");
	//TCanvas * Plat = new TCanvas(graphname.c_str(), graphname.c_str());
	f1->SetTitle(graphname.c_str());
	f1->SetLineColor(2);
	f1->SetLineStyle(4);
	//f2->Draw("SAME");
	//TCanvas * Plat2 = new TCanvas("dummy", "dummy");
	f1->Draw();
	grav->Draw("SAME P");
	grav->Write(graphname.c_str());
	//Plat->Write(graphname.c_str());
	//Plat->Clear();
    /**************************************************/
}
/// 	For cosmetics
/**  Drawing q^2 graph with some cosmetics */
void
trQCD::Draw_and_Fit(vector <Double_t> data, vector <Double_t> Gm, vector <Double_t> lambda, double charge , vector <Double_t> qsquare, vector <Double_t> stdsapma,Double_t fitconst, int select ,string graphname)
{
	int boyut= qsquare.size()-1;
	
	TF1 *f1;
	switch(select)
	{
		case 0 : //Monopole form
			f1 = new TF1("f1","[0]/(1+x/[1]^2)",qsquare[0], qsquare[boyut]);
		break;
		case 1 : //Dipole form
       	 	f1 = new TF1("f1","[0]/(1+x/[1]^2)^2",qsquare[0], qsquare[boyut]);
    	break;
    	case 2 : //Exponential form
       	 	f1 = new TF1("f1","[0]*exp(-x/([1])^2)",qsquare[0], qsquare[boyut]);
    	break;
	}
	TRandom* random1 = new TRandom;;
	double par2 = random1->Uniform(0,5);
	f1->SetParameters(fitconst, par2);
	TGraphErrors *grav = new TGraphErrors(data.size(),&(qsquare[0]),&(data[0]),0,&(stdsapma[0]));
	grav->Fit(f1, "RQN+","S");
	grav->Fit(f1, "RQN+","S");
	grav->Fit(f1, "RQN+","S");
	grav->Fit(f1, "RQN+","S");
	grav->Fit(f1, "RN+","S");
	f1->SetLineColor(2);
	cout<<graphname<<" Chi square: "<<f1->GetChisquare()<<" P Value: "<<f1->GetProb()<<endl;
	grav->SetMarkerColor(4);
	grav->SetMarkerStyle(21);
	TGraph *errorband= ErrorBand(Gm,lambda, qsquare,0.00001,select);
	TCanvas * Plat = new TCanvas(graphname.c_str(), graphname.c_str());
	errorband->SetFillStyle(3003);
	errorband->SetTitle(" ");
	errorband->Draw("AF SAME");
    errorband->GetXaxis()->SetTitle("Q^{ 2} [GeV^{ 2}]");
    errorband->GetYaxis()->SetTitle("G_{E}");
	f1->Draw("SAME");
	grav->Draw("P SAME");
	Plat->Write(graphname.c_str());
	Plat->Clear();
}
/// 	For cosmetics
/**  Drawing q^2 graph with some cosmetics */
void
trQCD::Draw_and_Fit(vector <Double_t> data, vector <Double_t> Gm, vector <Double_t> lambda, vector <Double_t> qsquare, vector <Double_t> stdsapma,Double_t fitconst, int select ,string graphname)
{
	int boyut= qsquare.size()-1;
	
	TF1 *f1, *f2;
	switch(select)
	{
		case 0 : //Monopole form
			f1 = new TF1("f1","[0]/(1+x/[1]^2)",qsquare[1], qsquare[boyut]);
	    	f2= new TF1("f2","[0]/(1+x/[1]^2)",-0.1,qsquare[boyut]+0.1);
		break;
		case 1 : //Dipole form
       	 	f1 = new TF1("f1","[0]/(1+x/[1]^2)^2",qsquare[1], qsquare[boyut]);
       	 	f2 = new TF1("f2","[0]/(1+x/[1]^2)^2",-0.1, qsquare[boyut]+0.1);
    	break;
    	case 2 : //Exponential form
       	 	f1 = new TF1("f1","[0]*exp(-x/([1])^2)",qsquare[1], qsquare[boyut]);
       	 	f2 = new TF1("f2","[0]*exp(-x/([1])^2)",-0.1, qsquare[boyut]+0.1);
    	break;
	}
	TRandom* random1 = new TRandom;;
	double par2 = random1->Uniform(0,5);
	f1->SetParameters(fitconst, par2);
//	f1->FixParameter(0,fitconst);
//	f1->FixParameter(1,avv_lambda);
	TGraphErrors *grav = new TGraphErrors(data.size(),&(qsquare[0]),&(data[0]),0,&(stdsapma[0]));
	grav->Fit(f1, "RQN+","S");
	grav->Fit(f1, "RQN+","S");
	grav->Fit(f1, "RQN+","S");
	grav->Fit(f1, "RQN+","S");
	grav->Fit(f1, "RN+","S");
	f2->SetParameters(f1->GetParameter(0), f1->GetParameter(1));
	f1->SetLineColor(2);
   	f2->SetLineColor(2);
	f2->SetLineStyle(4);
	cout<<graphname<<" Chi square: "<<f1->GetChisquare()<<" P Value: "<<f1->GetProb()<<endl;
	grav->SetMarkerColor(4);
	grav->SetMarkerStyle(21);
  //      TGraphErrors *grav2 = new TGraphErrors(data2.size(),&(qsqr2[0]),&(data2[0]),0,&(stdsapma2[0]));
	TGraph *errorband= ErrorBand(Gm,lambda, qsquare,0.0001,select);
	TCanvas * Plat = new TCanvas(graphname.c_str(), graphname.c_str());
	//f2->SetTitle(" ");
	errorband->SetFillStyle(3003);
	errorband->SetTitle(" ");
	errorband->Draw("AF SAME");
    errorband->GetXaxis()->SetTitle("Q^{ 2} [GeV^{ 2}]");
    errorband->GetYaxis()->SetTitle("G_{M}");
	f2->Draw("SAME");
	f1->Draw("SAME");
	grav->RemovePoint(0.0);
	grav->Draw("P SAME");
	Plat->Write(graphname.c_str());
	Plat->Clear();
}



/// Error band function
/**   */
TGraph* 
trQCD::ErrorBand(vector <Double_t> data,vector <Double_t> lambda, vector <Double_t> qsquare, Double_t stepsize, int select)
{
	int boyut= qsquare.size()-1;
	int nbdata = (qsquare[boyut]+0.1)/stepsize;
	Double_t qsqr[2*nbdata];
	Double_t errors[2*nbdata];
	int counter=0;
	for(Double_t j=0.0; j<qsquare[boyut]+0.1; j=j+stepsize){
        vector<Double_t> yvector;
        Double_t error = 0.0;
		for(int i = 0; i<data.size(); i++){
			switch(select)
        		{
        		case 0 : //Monopole form
					yvector.push_back(data[i]/((1+j/(lambda[i]*lambda[i]))));
        		break;
        		case 1 : //Dipole form
					yvector.push_back(data[i]/((1+j/(lambda[i]*lambda[i]))*(1+j/(lambda[i]*lambda[i]))));
		        break;
        		case 2 : //Exponential form
		        	yvector.push_back(data[i]*exp(-j/(lambda[i]*lambda[i])));
		        break;
        		}
		}
		Double_t ort=StdDev (yvector, error);
		yvector.clear();
		qsqr[2*nbdata-counter-1]=j;
		qsqr[counter]=j;
		errors[2*nbdata-counter-1]=ort+error;
		errors[counter]=ort-error;
		counter++;
	}
	TGraph *gr = new TGraph(2*nbdata,qsqr,errors);
	return gr;
}

/// 	Calculate average of Qsquare
/**  Calculate ge qsquare */
void
trQCD::central_fit (TMatrixD datamatrix_jkknf, vector <Double_t> qsquare, double charge , int select, string graphname)
{
	int boyut= qsquare.size()-1;
    int		nbcol = datamatrix_jkknf.GetNcols();
    int		nbrow = datamatrix_jkknf.GetNrows();
    double parz=0;
	std::vector<Double_t> aver;
	std::vector<Double_t> aver_err;
	std::vector<Double_t> qsqr2;
	for(int i=0; i<qsquare.size();i++){
		qsqr2.push_back(qsquare[i]);
	}
	for (int k = 0; k < nbrow; k++) {
		if(k==7)
			k++;
		double hata= 0.0;
		double ort= 0.0;
		std::vector<Double_t> row;
       	row = GetRows(datamatrix_jkknf, k);
		ort=StdDev (row, hata);
		aver.push_back(ort);
		aver_err.push_back(hata);
	}
	TRandom* random1 = new TRandom;
	TF1 *f1;
	switch(select)
	{
		case 0 : //Monopole form
			f1 = new TF1("f1","[0]/(1+x/[1]^2)",0.0, qsqr2[boyut]);
		break;
		case 1 : //Dipole form
       	 	f1 = new TF1("f1","[0]/(1+x/[1]^2)^2",0.0, qsqr2[boyut]);
	    	break;
    		case 2 : //Exponential form
       	 	f1 = new TF1("f1","[0]*TMath::Exp(-x/[1]^2)",0.0, qsqr2[boyut]);
	    	break;
	}
	double par1 = charge;
	double par2= random1->Uniform(0,5);
	f1->SetParameters(par1, par2);
	TGraphErrors *grav = new TGraphErrors(aver.size(),&(qsqr2[0]),&(aver[0]),0,&(aver_err[0]));
	f1->FixParameter(0,par1);
	grav->Fit(f1, "RQN+","S");
	grav->Fit(f1, "RQN+","S");
	grav->Fit(f1, "RQN+","S");
	grav->Fit(f1, "RQN+","S");
	grav->Fit(f1, "RQN+","S");
	grav->SetMarkerColor(4);
	grav->SetMarkerStyle(21);
	grav->SetTitle(graphname.c_str());
	grav->GetXaxis()->SetTitle("Q^2[GeV^-1]");
	grav->GetYaxis()->SetTitle("Average");
	grav->Draw("ALP");
	TCanvas * Plat = new TCanvas(graphname.c_str(), graphname.c_str());
	f1->SetTitle(graphname.c_str());
	f1->SetLineColor(2);
	f1->SetLineStyle(4);
	f1->Draw();
	grav->Draw("P");
	Plat->Write(graphname.c_str());
	Plat->Clear();
}

/// 	Calculate average of Qsquare
/**  Calculate gm qsquare */
void
trQCD::central_fit (TMatrixD datamatrix_jkknf, vector <Double_t> qsquare, int select, string graphname)
{
	int boyut= qsquare.size()-1;
        /**************************************************/
    int		nbcol = datamatrix_jkknf.GetNcols();
    int		nbrow = datamatrix_jkknf.GetNrows();
    double parz=0;
	std::vector<Double_t> aver;
	std::vector<Double_t> aver_err;
	std::vector<Double_t> qsqr2;
	for(int i=0; i<qsquare.size();i++){
		qsqr2.push_back(qsquare[i]);
	}
	for (int k = 0; k < nbrow; k++) {
		if(k==7)
			k++;
		double hata= 0.0;
		double ort= 0.0;
		std::vector<Double_t> row;
       	row = GetRows(datamatrix_jkknf, k);
		ort=StdDev (row, hata);
		aver.push_back(ort);
		aver_err.push_back(hata);
	}
	TRandom* random1 = new TRandom;
	TF1 *f1, *f2;
	switch(select)
	{
		case 0 : //Monopole form
			f1 = new TF1("f1","[0]/(1+x/[1]^2)",0.0, qsqr2[boyut]);
	    	f2= new TF1("f2","[0]/(1+x/[1]^2)",-0.1,qsqr2[boyut]+0.1);
		break;
		case 1 : //Dipole form
       	 	f1 = new TF1("f1","[0]/(1+x/[1]^2)^2",0.0, qsqr2[boyut]);
       	 	f2 = new TF1("f2","[0]/(1+x/[1]^2)^2",-0.1, qsqr2[boyut]+0.1);
    	break;
    	case 2 : //Exponential form
       	 	f1 = new TF1("f1","[0]*TMath::Exp(-x/[1]^2)",0.0, qsqr2[boyut]);
       	 	f2 = new TF1("f2","[0]*TMath::Exp(-x/[1]^2)",-0.1, qsqr2[boyut]+0.1);
    	break;
	}
	qsqr2.erase(qsqr2.begin());
	aver.erase(aver.begin());
	aver_err.erase(aver_err.begin());
	double par1 = random1->Rndm();
	double par2 = random1->Uniform(0,5);
	f1->SetParameters(par1, par2);
	TGraphErrors *grav = new TGraphErrors(aver.size(),&(qsqr2[0]),&(aver[0]),0,&(aver_err[0]));
	grav->Fit(f1, "RQN+","S");
	grav->Fit(f1, "RQN+","S");
	grav->Fit(f1, "RQN+","S");
	grav->Fit(f1, "RQN+","S");
	grav->Fit(f1, "RQN+","S");
	grav->SetMarkerColor(4);
	grav->SetMarkerStyle(21);
	grav->GetXaxis()->SetTitle("Q^2[GeV^-1]");
	grav->GetYaxis()->SetTitle("Average");
	TCanvas * Plat = new TCanvas(graphname.c_str(), graphname.c_str());
	f2->SetParameters(f1->GetParameter(0), f1->GetParameter(1));
	f2->SetLineColor(2);
	f2->SetLineStyle(4);
	f1->SetTitle(graphname.c_str());
    	f2->Draw();
	grav->Draw("P");
	Plat->Write(graphname.c_str());
	Plat->Clear();			
	aver.clear();
    /**************************************************/
}



/// 	Calculate average of Qsquare
/**  Calculate ge vs qsquare */
vector <Double_t>
trQCD::stat_fit (vector <vector <Double_t> > datamatrix_jkknf, vector<Double_t> stddev, TMatrixD qsquare, vector<Double_t>& lambda ,double charge ,int select, string graphname)
{
	int boyut= stddev.size()-1;
    /**************************************************/
    int		nbcol = qsquare.GetNcols();
    int		nbrow = qsquare.GetNrows();
    std::vector<Double_t> x;
    std::vector<Double_t> G_QSQR;
    double parz=0;

	Double_t temp_fit_1=0.0;
	Double_t temp_fit_2=0.0;
    for (int i = 0; i < nbcol; i++) {
        std::vector<Double_t> colon;
        std::vector<Double_t> colon_x;
		for(int j=0; j<datamatrix_jkknf.size(); j++)
        colon.push_back(datamatrix_jkknf[j][i]);
		colon_x = GetCols(qsquare, i);
        string          isim   (graphname + itos(i));
		TRandom* random1 = new TRandom;
		Double_t par1=0.0;
		Double_t par2=0.0;	
		if(i==0) {
			par1 = random1->Uniform(0,5);
			par2 = random1->Rndm();
		}
		else {
			par1=temp_fit_1;
			par2=temp_fit_2;
		}
		colon_x.erase(colon_x.begin()+7);
		TF1 *f1;
		switch(select)
		{
			case 0 : //Monopole form
					f1 = new TF1("f1","[0]/(1+x/[1]^2)",0.0, colon_x[boyut]);
			break;
			case 1 : //Dipole form
		       	 	f1 = new TF1("f1","[0]/(1+x/[1]^2)^2",0.0, colon_x[boyut]);
	    		break;
	    		case 2 : //Exponential form
		       	 	f1 = new TF1("f1","[0]*TMath::Exp(-x/[1]^2)",0.0, colon_x[boyut]);
	    		break;
		}
   		f1->SetParameters(par1, charge);
		//f1->FixParameter(0,charge);
		TCanvas * Plat = new TCanvas(isim.c_str(), isim.c_str());
		TGraphErrors *gr = new TGraphErrors(colon_x.size(),&(colon_x[0]),&(colon[0]),0,&(stddev[0]));
		gr->Fit(f1, "REQN+","S");
		gr->Fit(f1, "REQN+","S");
		gr->Fit(f1, "REQN+","S");
		gr->Fit(f1, "REQN+","S");
		gr->Fit(f1, "REQN+","S");
		gr->Fit(f1, "REQN+","S");
		gr->SetMarkerColor(4);
		gr->SetMarkerStyle(21);
		f1->SetLineColor(2);
		f1->SetTitle(graphname.c_str());
		f1->Draw();
		gr->Draw("P");
		Plat->Write(isim.c_str());
		G_QSQR.push_back(f1->GetParameter(0));
		lambda.push_back(f1->GetParameter(1));
		temp_fit_1=f1->GetParameter(0);
		temp_fit_2=f1->GetParameter(1);
		colon_x.clear();
		colon.clear();
		gr->Clear();
		Plat->Clear();
	}
	return G_QSQR;
}

/// 	Calculate average of Qsquare
/**  Calculate gm vs qsquare */
vector <Double_t>
trQCD::stat_fit (vector <vector <Double_t> > datamatrix_jkknf, vector<Double_t> stddev, TMatrixD qsquare, vector<Double_t>& lambda ,int select, string graphname)
{
	int boyut= stddev.size()-1;
    /**************************************************/
	
    int		nbcol = qsquare.GetNcols();
    int		nbrow = qsquare.GetNrows();
    std::vector<Double_t> x;
    std::vector<Double_t> G_QSQR;
    double parz=0;
	Double_t temp_fit_1=0.0;
	Double_t temp_fit_2=0.0;
    for (int i = 0; i < nbcol; i++) {
        std::vector<Double_t> colon;
        std::vector<Double_t> colon_x;
		for(int j=0; j<datamatrix_jkknf.size(); j++)
       		colon.push_back(datamatrix_jkknf[j][i]);
		colon_x = GetCols(qsquare, i);
        string          isim   (graphname + itos(i));
		TRandom* random1 = new TRandom;
		Double_t par1=0.0;
		Double_t par2=0.0;	
		if(i==0) {
			par1 = random1->Uniform(0,5);
			par2 = random1->Rndm();
		}
		else {
			par1=temp_fit_1;
			par2=temp_fit_2;
		}
		colon_x.erase(colon_x.begin()+7);
		colon_x.erase(colon_x.begin());
		TF1 *f1;
		switch(select)
		{
			case 0 : //Monopole form
				f1 = new TF1("f1","[0]/(1+x/[1]^2)",0.0, colon_x[boyut]);
				break;
			case 1 : //Dipole form
		       	 f1 = new TF1("f1","[0]/(1+x/[1]^2)^2",0.0, colon_x[boyut]);
	    			break;
	    	case 2 : //Exponential form
		       	 f1 = new TF1("f1","[0]*TMath::Exp(-x/[1]^2)",0.0, colon_x[boyut]);
	    	break;
		}
	   	f1->SetParameters(par1, par2);
		TCanvas * Plat = new TCanvas(isim.c_str(), isim.c_str());
		TGraphErrors *gr = new TGraphErrors(colon_x.size(),&(colon_x[0]),&(colon[0]),0,&(stddev[0]));
		gr->Fit(f1, "REQN+","S");
		gr->Fit(f1, "REQN+","S");
		gr->Fit(f1, "REQN+","S");
		gr->Fit(f1, "REQN+","S");
		gr->Fit(f1, "REQN+","S");
		gr->Fit(f1, "REQN+","S");
		gr->SetMarkerColor(4);
		gr->SetMarkerStyle(21);
		f1->SetLineColor(2);
		f1->SetTitle(graphname.c_str());
		f1->Draw();
        	gr->Draw("P");
		Plat->Write(isim.c_str());
		G_QSQR.push_back(f1->GetParameter(0));
		lambda.push_back(f1->GetParameter(1));
		temp_fit_1=f1->GetParameter(0);
		temp_fit_2=f1->GetParameter(1);
		colon_x.clear();
		x.clear();
		colon.clear();
		Plat->Clear();
		gr->Clear();
		
    }
	
return G_QSQR;
    /**************************************************/
}


/// 	Calculate average of Qsquare fit function
/**  Calculate gm and ge vs qsquare */
vector <Double_t>
trQCD::GetQsqrAver (TMatrixD datamatrix_jkknf,int sink ,int qsqr ,int vfirst, int vlast,string graphname)
{
    /**************************************************/
    //cout<<"........................Calculating Q^2 Average plot "<< graphname <<" "<<qsqr<<"........................"<<endl;
    int		nbcol = datamatrix_jkknf.GetNcols();
    int		nbrow = datamatrix_jkknf.GetNrows();
    std::vector<Double_t> x;
    int		fitarashwl [] = {vfirst, vlast};
    std::vector<Double_t> aver;
    std::vector<Double_t> stddev;
	
    for (int i = 0; i < nbrow; i++) {
        x.push_back(i);
    }
	stddev=StdDev(datamatrix_jkknf,aver,sink);
    	string          isim2   (graphname+"_"+itos(qsqr)+"_average");
	TGraphErrors *grav = new TGraphErrors(aver.size(),&(x[0]),&(aver[0]),0,&(stddev[0]));
	TF1 * f1 = new TF1("f1", "[0]", vfirst, vlast);
	grav->Fit(f1, "RQN+");
	grav->Fit(f1, "RQN+");
	grav->Fit(f1, "RQN+");
	grav->Fit(f1, "RQN+");
	grav->Fit(f1, "RQN+");
	grav->Fit(f1, "RQN+");
	grav->Fit(f1, "RQN");
	char		parameter [100];
	sprintf(parameter, "p0=%e", f1->GetParameter(0));
	TLegend        *legend = new TLegend(0.70, 0.90, 0.99, 0.99);
	legend->SetTextFont(50);
	legend->SetTextSize(0.04);
	legend->AddEntry(f1, parameter, "l");
	legend->Draw();
    	grav->SetTitle(isim2.c_str());
	grav->SetMarkerColor(4);
	grav->SetMarkerStyle(21);
	grav->Draw("ALP");
	f1->SetLineColor(2);
	f1->Draw("SAME");
	grav->Write(isim2.c_str());
	aver.clear();
	grav->Clear();
	return stddev;
    /**************************************************/
}


/// 	Calculate Qsquare fit function
/**  Calculate gm and ge vs qsquare */
vector <Double_t>
trQCD::GetQsqr (TMatrixD datamatrix_jkknf, int qsqr ,int vfirst, int vlast, string graphname)
{
    /**************************************************/
    //cout<<"........................Calculating Q^2 plots "<<graphname<<" "<<qsqr<<"........................"<<endl;
    int		nbcol = datamatrix_jkknf.GetNcols();
    int		nbrow = datamatrix_jkknf.GetNrows();
    std::vector<Double_t> x;
    std::vector<Double_t> resultz;
    std::vector<Double_t> aver;
    std::vector<Double_t> stddev;
	stddev=StdDev(datamatrix_jkknf,aver,nbrow);
    int		fitarashwl [] = {vfirst, vlast};
    string          isim   (graphname + itos(qsqr));
    //datamatrix_jkknf.Write(isim.c_str());
    for (int i = 0; i < nbrow; i++) {
        x.push_back(i);
    }
	string          isim2   (graphname+"_"+itos(qsqr)+"_average");
	TGraphErrors *grav = new TGraphErrors(aver.size(),&(x[0]),&(aver[0]),0,&(stddev[0]));
	TF1 * f1 = new TF1("f1", "[0]", vfirst, vlast);
	grav->Fit(f1, "RQN+");
	grav->Fit(f1, "RQN+");
	grav->Fit(f1, "RQN+");
	grav->Fit(f1, "RQN+");
	grav->Fit(f1, "RQN+");
	grav->Fit(f1, "RQN+");
	grav->Fit(f1, "RQN");
	grav->SetMarkerColor(4);
	grav->SetMarkerStyle(21);
	grav->Draw("ALP");
	f1->SetLineColor(2);
	f1->Draw("SAME");
	grav->Write(isim2.c_str());
	Double_t fitosko = f1->GetParameter(0);
    for (int i = 0; i < nbcol; i++) {
    	string          isim3   (graphname+"_"+itos(qsqr)+"_"+itos(i));
    	std::vector<Double_t> colon;
		colon = GetCols(datamatrix_jkknf, i);
		TGraphErrors *gr = new TGraphErrors(colon.size(),&(x[0]),&(colon[0]),0,&(stddev[0]));
		TF1 * f2 = new TF1("f2", "[1]", vfirst, vlast);
		f2->SetParameter(1,fitosko);
		gr->Fit(f2, "REQN+");
		gr->Fit(f2, "REQN+");
		gr->Fit(f2, "REQN+");
		gr->Fit(f2, "REQN+");
		gr->Fit(f2, "REQN+");
		gr->Fit(f2, "REQN+");
		gr->Fit(f2, "RQ");
		char		parameter [100];
		resultz.push_back(f2->GetParameter(1));
		gr->SetMarkerColor(4);
		gr->SetMarkerStyle(21);
		gr->Draw("ALP");
		f2->SetLineColor(2);
		f2->Draw("SAME");
		//gr->Write(isim2.c_str());
		colon.clear();
		gr->Clear();
	}
	
	return resultz;
    /**************************************************/
}


/// 	Calculate Qsquare fit function
/**  Calculate gm and ge vs qsquare */
vector <Double_t>
trQCD::GetQsqr (TMatrixD datamatrix_jkknf, vector<Double_t> stddev ,int qsqr ,int vfirst, int vlast,string graphname)
{
    /**************************************************/
    //cout<<"........................Calculating Q^2 plots "<<graphname<<" "<<qsqr<<"........................"<<endl;
    int		nbcol = datamatrix_jkknf.GetNcols();
    int		nbrow = datamatrix_jkknf.GetNrows();
    std::vector<Double_t> x;
    std::vector<Double_t> resultz;
    int		fitarashwl [] = {vfirst, vlast};
    string          isim   (graphname + itos(qsqr));
    //datamatrix_jkknf.Write(isim.c_str());
    for (int i = 0; i < nbrow; i++) {
        x.push_back(i);
    }
    for (int i = 0; i < nbcol; i++) {
    	string          isim2   (graphname+"_"+itos(qsqr)+"_"+itos(i));
    	std::vector<Double_t> colon;
		colon = GetCols(datamatrix_jkknf, i);
		TGraphErrors *gr = new TGraphErrors(colon.size(),&(x[0]),&(colon[0]),0,&(stddev[0]));
		TF1 * f1 = new TF1("f1", "[0]", vfirst, vlast);
		gr->Fit(f1, "REQN+");
		gr->Fit(f1, "REQN+");
		gr->Fit(f1, "REQN+");
		gr->Fit(f1, "REQN+");
		gr->Fit(f1, "REQN+");
		gr->Fit(f1, "REQN+");
		gr->Fit(f1, "RQ");
		char		parameter [100];
		resultz.push_back(f1->GetParameter(0));
		gr->SetMarkerColor(4);
		gr->SetMarkerStyle(21);
		gr->Draw("ALP");
		f1->SetLineColor(2);
		f1->Draw("SAME");
		//gr->Write(isim2.c_str());
		colon.clear();
		gr->Clear();

	}
	return resultz;
    /**************************************************/
}



/// 	TvectorD to TMatrixD  conversion
/**  TVectorD to TMatrixD  conversion 1 row dimension colon matrix */
//TMatrixD
//trQCD::VectortoMatrix (TVectorD v1, int size){
//    TMatrixD son(1,size);
//    for(int t=0;t<size;t++)
//        son(0,t) = v1(t);
//    return son;
//}
///// 	TvectorD to TMatrixD  conversion
///**  TVectorD to TMatrixD  conversion dim row 1 dimension colon matrix */
//TMatrixD
//trQCD::VectortoMatrixT (TVectorD v1, int size){
//    TMatrixD son(size,1);
//    for(int t=0;t<size;t++)
//        son(t,0) = v1(t);
//    return son;
//}

/// 	Write Routine to text
/**  Write the matrices to files */
void
trQCD::Write_File (vector <TMatrixD> Mat_t, vector <TMatrixD> Mat_t0, vector <TMatrixD> Mat_dt,string t_file,string t0_file,string dt_file)
{
	int dim = Mat_t0.at(0).GetNcols();
	
	
	for(unsigned i=0; i<Mat_t0.size(); i++){
	ofstream	Matrixt,Matrixt0,Matrixdt;
 	ostringstream os_1;
  	os_1 <<"Corrmat/t/"<< t_file <<"_c"<<i;
  	string t_file_c= os_1.str();
 	ostringstream os_2;
  	os_2 <<"Corrmat/dt/"<< dt_file <<"_c"<<i;
  	string dt_file_c= os_2.str();
 	ostringstream os_3;
  	os_3 <<"Corrmat/t0/"<< t0_file <<"_c"<<i;
  	string t0_file_c = os_3.str();

	Matrixt.open(t_file_c.c_str());
	Matrixt0.open(t0_file_c.c_str());
	Matrixdt.open(dt_file_c.c_str()); 
		
		for(int fff=0; fff<dim; fff++){
	    	for(int eee=0; eee<dim; eee++){
		    	Double_t tempnum1 =Mat_t.at(i)(fff,eee);
			Double_t tempnum2 =Mat_t0.at(i)(fff,eee);
			Double_t tempnum3= Mat_dt.at(i)(fff,eee);
			if(eee<dim-1){
  		    	Matrixt << tempnum1 <<"\t";
		    	Matrixt0 << tempnum2 <<"\t";
			Matrixdt <<tempnum3<<"\t";
			}
			else{
  		    	Matrixt << tempnum1 ;
		    	Matrixt0 << tempnum2;
			Matrixdt <<tempnum3;
			}
			}
		    if(fff<dim-1){	
  		    Matrixt << endl;
		    Matrixt0 << endl;
		    Matrixdt <<endl;
		    }
		}	
	
Matrixt   <<endl;
			Matrixt0   <<endl;
			Matrixdt  <<endl;
	}

}



/// 	Calculates Eigenvalues and eigenvectors of correlation matrices
/**  This function can be used for excited states calculations it calculates 0905.3616 eq(21) G_{ij} */
vector < vector  <Double_t> >
trQCD::Eigen_Mass (vector <TMatrixD> Mat_t, vector <TMatrixD> Mat_t0, vector <TMatrixD> Mat_dt)
{
	vector <vector  <Double_t> > Sonuc;
	int dim = Mat_t0.at(0).GetNcols();
	for(unsigned i=0; i<Mat_t0.size(); i++){
		TMatrixD Mat_invt0(dim,dim);
		Double_t minim = Mat_t0.at(i)(0,0);
		for(int rowz=0; rowz<dim;rowz++){
			for(int colz=0;colz<dim;colz++){
				if(Mat_t0.at(i)(rowz,colz)<minim){
					minim=Mat_t0.at(i)(rowz,colz);
				}
				if(Mat_t0.at(i)(rowz,colz)==0)
					cout<<"WARNING MATRIX HAS VALUE 0 PLEASE CHECK MATRIX, Inversion might be wrong!!!!!"<<endl;
			}
		}
		//TMatrixD Matrix_t0 = Mat_t0.at(i);
		
		TMatrixD H_square = Mat_t0.at(i);
		TDecompLU lu(H_square);
		lu.SetTol(1.e-30);
		TMatrixD H4 = lu.Invert();
						
		TMatrixD M1(dim,dim);
		TMatrixD M2(dim,dim);
		M1 = H4 * Mat_dt.at(i);
		M2 = Mat_dt.at(i) * H4;
		TMatrixD M3 (dim,dim);
		M3.Transpose(M2);
		
	/*	ofstream	Matrix1,Matrix2,Matrix_t;
		Matrix1.open("M1.txt");
		Matrix2.open("M2.txt");
		Matrix_t.open("Mat_t.txt"); 
		
		
		for(int fff=0; fff<dim; fff++){
		    for(int eee=0; eee<dim; eee++){
			    Double_t tempnum1 =M1(fff,eee);
				Double_t tempnum2 =M3(fff,eee);
				//Double_t tempnum3= Matrix_t(fff,eee);
	  		    Matrix1 << tempnum1 <<"\t";
			    Matrix2 << tempnum2 <<"\t";
				//Matrix_t <<tempnum3<<"\t";
			}	
			Matrix1   <<endl;
			Matrix2   <<endl;
			Matrix_t  <<endl;
		}*/
		cout<<"looop"<<endl;
		M1.Print();

        TMatrixDEigen eigen(M1);
       	TMatrixD eigenVal = eigen.GetEigenValues();
       // TMatrixD eigenVecUT = eigen.GetEigenVectors();
		TMatrixD eigenVecU = eigen.GetEigenVectors();
        TMatrixDEigen eigen2(M3);
       	// TMatrixD eigenVal2 = eigen2.GetEigenValues();
       // TMatrixD eigenVecVT = eigen2.GetEigenVectors();
		TMatrixD eigenVecV = eigen2.GetEigenVectors();
       // TMatrixD eigenVecU(dim,dim);
		//TMatrixD eigenVecV(dim,dim);
		
		//eigenVecU.Transpose(eigenVecUT);
		//eigenVecV.Transpose(eigenVecVT);
		
		vector < TMatrixD > eigenU;
		vector < TMatrixD > eigenV;
		//eigenVecU.Print();
		//eigenVecV.Print();
		
		for (Int_t ir = 0; ir < eigenVecU.GetNrows(); ir++)
		{
		  const TVectorD v = TMatrixDRow_const(eigenVecU,ir);
		  const Double_t norm = TMath::Sqrt(v.Norm2Sqr());
		  TMatrixDRow(eigenVecU,ir) *= 1/norm;
		}
		
	
		
		eigenVecU.Print();
		//eigenVal.Print();
		cout<<"looop ends"<<endl;
		
		
		
		
		
		
				
		for(int alpha=0; alpha<dim; alpha++){
			vector <Double_t> U = GetCols(eigenVecU,alpha);
			vector <Double_t> V = GetRows(eigenVecV,alpha);
			TMatrixD tem_eigU(dim,1);
			TMatrixD tem_eigV(1,dim);
			for(int ss=0; ss<dim; ss++){	
				tem_eigU(ss,0) = U.at(ss);
				tem_eigV(0,ss) = V.at(ss);				
			}
			eigenU.push_back(tem_eigU);
			eigenV.push_back(tem_eigV);
        }
		
		vector <Double_t> soonuc;
	
		for(int alpha=0; alpha<dim; alpha++){
			//TMatrixD son = eigenV.at(alpha)*Mat_t.at(i)*eigenU.at(alpha); boyle olmasi gerekmez mi?
			TMatrixD son = eigenV.at(alpha)*Mat_t.at(i)*eigenU.at(alpha);
        	Double_t sonD=son(0,0);
			soonuc.push_back(sonD);
		}
		Sonuc.push_back(soonuc);
	}		
	return Sonuc;
  /*myfile.close();
  myfile2.close();
  myfile3.close();*/
  
}

//vector <vector < > > to TMatrixD conversion.
TMatrixD
trQCD::VectortoMatrix (vector <vector < vector < Double_t > > >vec, int alpha){
	int colz = vec.size();
	int rowz = vec[0].size();
    TMatrixD son(colz,rowz);
    for(int kol=0;kol<colz;kol++)
		for(int rowzz=0;rowzz<rowz;rowzz++)
        	son(kol,rowzz) = vec[kol][rowzz][alpha];
    return son;
}



/// 	Creates 2x2 correlation matrix
/**  This function can be used for excited states calculations it calculates 0905.3616 eq(21) G_{ij} */
vector <TMatrixD>
trQCD::Correlation_Matrix (TMatrixD Mat1_1, TMatrixD Mat1_g5,TMatrixD Matg5_1,TMatrixD Matg5_g5, TFile *filename, int begin,int end, string matrix_name)
{
	//TFile	*cfile = new TFile("corrmats.root", "RECREATE");	
	vector <TMatrixD> Corr_Mats;
	int		nbconf = Mat1_1.GetNcols();
	for (int t=begin; t<end; t++)
	{
		for (int j = 0; j < nbconf; j++) {
		std::string tempr (matrix_name+"_t"+itos(t)+"_c"+itos(j)); 
		TMatrixD TempMatrix(2,2);
			TempMatrix(0,0) = Mat1_1(t,j);
			TempMatrix(0,1) = Mat1_g5(t,j);
			//TempMatrix(0,2) = Mat1_g5g4(t,j);
			TempMatrix(1,0) = Matg5_1(t,j);
			TempMatrix(1,1) = Matg5_g5(t,j);
			/*TempMatrix(1,2) = Matg5_g5g4(t,j);
			TempMatrix(2,0) = Matg5g4_1(t,j);
			TempMatrix(2,1) = Matg5g4_g5(t,j);
			TempMatrix(2,2) = Matg5g4_g5g4(t,j);*/
			Corr_Mats.push_back(TempMatrix);
			//TempMatrix.Print();
			TempMatrix.Write(tempr.c_str());
			filename->Write();
		}
	}
	return Corr_Mats;
	filename->Write();
	//cfile->Close();
}





/// 	Creates 3x3 correlation matrix
/**  This function can be used for excited states calculations it calculates 0905.3616 eq(21) G_{ij} */
vector <TMatrixD>
trQCD::Correlation_Matrix (TMatrixD Mat1_1, TMatrixD Mat1_g5,TMatrixD Mat1_g5g4,TMatrixD Matg5_1,TMatrixD Matg5_g5,TMatrixD Matg5_g5g4,TMatrixD Matg5g4_1,TMatrixD Matg5g4_g5,TMatrixD Matg5g4_g5g4, TFile *filename, int begin,int end, string matrix_name)
{
	//TFile	*cfile = new TFile("corrmats.root", "RECREATE");	
	vector <TMatrixD> Corr_Mats;
	int		nbconf = Mat1_1.GetNcols();
	for (int t=begin; t<end; t++)
	{
		for (int j = 0; j < nbconf; j++) {
		std::string tempr (matrix_name+"_t"+itos(t)+"_c"+itos(j)); 
		TMatrixD TempMatrix(3,3);
			TempMatrix(0,0) = Mat1_1(t,j);
			TempMatrix(0,1) = Mat1_g5(t,j);
			TempMatrix(0,2) = Mat1_g5g4(t,j);
			TempMatrix(1,0) = Matg5_1(t,j);
			TempMatrix(1,1) = Matg5_g5(t,j);
			TempMatrix(1,2) = Matg5_g5g4(t,j);
			TempMatrix(2,0) = Matg5g4_1(t,j);
			TempMatrix(2,1) = Matg5g4_g5(t,j);
			TempMatrix(2,2) = Matg5g4_g5g4(t,j);
			Corr_Mats.push_back(TempMatrix);
			//TempMatrix.Print();
			TempMatrix.Write(tempr.c_str());
			filename->Write();
		}
	}
	return Corr_Mats;
	filename->Write();
	//cfile->Close();
}


/// 	GetMass function, calculate mass function
/**  Calculate the effective mass using logarithm of two point functions */
Double_t
trQCD::GetMassLog (TMatrixD datamatrix_jkknf, int vfirst, int vlast,string name)
{
    /**************************************************/
    int		nbcol = datamatrix_jkknf.GetNcols();
	int		nbrow = datamatrix_jkknf.GetNrows();
	
    TMatrixD	DeltaE_Jack2(nbrow, nbcol);
    std::vector<Double_t> x;
    std::vector<Double_t> colon;
    Double_t mass;
	std::vector<Double_t> weight;
    int		fitarashwl [] = {vfirst, vlast};
    
    for (int i = 0; i < nbrow-1; i++) {
        x.push_back(i+1);
        for (int j = 0; j < nbcol; j++) {
            DeltaE_Jack2(i, j) = (TMath::Log(TMath::Abs(datamatrix_jkknf(i, j) / datamatrix_jkknf(i+1, j))));
		}
    }
	for (int i = 0; i < nbrow-1; i++) {
		std::vector<Double_t> Row;
		Row = GetRows(DeltaE_Jack2,i);
		Double_t error = StdDev(Row);
		weight.push_back(error);
		Row.clear();
	}
	
	colon= average_sink(DeltaE_Jack2,nbrow-1);
    mass= FitterLog(colon, x, weight , fitarashwl, name);
    
    return mass;
    /**************************************************/
}

/// 	Fitter for mass calculation
/** It fits the constant [0] */
Double_t
trQCD::FitterLog(vector < Double_t > &data, vector < Double_t > &xarray, vector < Double_t > &stdDev,int *aralik,string name)
{
    int		fitbas = aralik[0];
    int		fitbit = aralik[1];
    std::vector<Double_t> zero (data.size(),0);
    TGraphErrors *gr = new TGraphErrors(data.size(),&(xarray[0]),&(data[0]),&(zero[0]),&(stdDev[0]));
    TF1 * f1 = new TF1("f1", "[0]", fitbas, fitbit);
    gr->Fit(f1, "RQN+");
    gr->Fit(f1, "RQN+");
    gr->Fit(f1, "RQN+");
    gr->Fit(f1, "RQN+");
    gr->Fit(f1, "RQN+");
    gr->Fit(f1, "RQN+");
    char		parameter [100];
    double		parz;
    f1->Draw();
    parz = f1->GetParameter(0);
    sprintf(parameter, "p0=%e", parz);
    TLegend        *legend = new TLegend(0.70, 0.90, 0.99, 0.99);
    legend->SetTextFont(50);
    legend->SetTextSize(0.04);
    legend->AddEntry(f1, parameter, "l");
    legend->Draw();
    gr->SetMarkerColor(4);
    gr->SetMarkerStyle(21);
    gr->Draw("ALP");
	gr->SetTitle(name.c_str());
    gr->GetXaxis()->SetTitle("Timestep");
    gr->GetYaxis()->SetTitle("Average of  #frac{C(n+1)}{C(n)}");
    gr->Write(name.c_str());
    
    return parz;
}


/// 	GetMass function, calculate mass function
/**  Calculate the effective mass using logarithm of two point functions */
std::vector<Double_t>
trQCD::GetMassLog (TMatrixD datamatrix_jkknf, vector <Double_t> weight, int vfirst, int vlast)
{
    /**************************************************/
    int		nbcol = datamatrix_jkknf.GetNcols();
	int		nbrow = datamatrix_jkknf.GetNrows();
	
    TMatrixD	DeltaE_Jack2(nbrow, nbcol);
    std::vector<Double_t> x;
    std::vector<Double_t> colon;
    std::vector<Double_t> mass;
    int		fitarashwl [] = {vfirst, vlast};
    
    for (int i = 0; i < nbrow; i++) {
        x.push_back(i);
        for (int j = 0; j < nbcol; j++) {
            DeltaE_Jack2(i, j) = (TMath::Log(TMath::Abs(datamatrix_jkknf(i, j) / datamatrix_jkknf((i + 1) % nbrow, j))));
		}
    }
    for (int i = 0; i < nbcol; i++) {
        colon = GetCols(DeltaE_Jack2, i);
        mass.push_back(FitterLog(colon, x, weight , fitarashwl));
    }
    return mass;
    /**************************************************/
}

/// 	GetMass function, calculate mass function
/**  Calculate the effective mass using two point functions */
std::vector<Double_t>
trQCD::GetMass (TMatrixD datamatrix_jkknf, vector <Double_t> weight, int vfirst, int vlast)
{
    /**************************************************/
    int		nbcol = datamatrix_jkknf.GetNcols();
    TMatrixD	DeltaE_Jack2(Nt, nbcol);
    std::vector<Double_t> x;
    std::vector<Double_t> colon;
    std::vector<Double_t> mass;
    int		fitarashwl [] = {vfirst, vlast};
    
    for (int i = 0; i < Nt; i++) {
        x.push_back(i);
        //for (int j = 0; j < nbcol; j++) {
            //DeltaE_Jack2(i, j) = (TMath::Log(TMath::Abs(datamatrix_jkknf(i, j) / datamatrix_jkknf((i + 1) % Nt, j))));
        	//DeltaE_Jack2(i, j) = datamatrix_jkknf(i, j) / datamatrix_jkknf((i + 1) % Nt, j);
	//}
    }
    for (int i = 0; i < nbcol; i++) {
        colon = GetCols(datamatrix_jkknf, i);
        mass.push_back(Fitter(colon, x, weight , fitarashwl));
    }
    return mass;
    /**************************************************/
}




/// 	GetMass function, calculate average of mass function
/**  Calculate the effective mass using two point functions */
Double_t 
trQCD::GetMassAver(TMatrixD datamatrix_jkknf, vector<Double_t> & stddev, int vfirst, int vlast) {
    /**************************************************/
	int			nbcol = datamatrix_jkknf.GetNcols();
	int			nbrow = datamatrix_jkknf.GetNrows();
	TMatrixD	DeltaE_Jack2(nbrow, nbcol);
        std::vector<Double_t> x;
        std::vector<Double_t> colon;
        std::vector<Double_t> mass;
        std::vector<Double_t> colon2;
        std::vector<Double_t> stddev2;

	int 		fitarashwl[]={vfirst,vlast};
	
	for (int i = 0; i < nbrow; i++) {
		x.push_back(i);
		//for (int j = 0; j < nbcol; j++) {
			//DeltaE_Jack2(i, j) = (TMath::Log(TMath::Abs(datamatrix_jkknf(i, j) / datamatrix_jkknf((i + 1) % Nt, j))));
			//DeltaE_Jack2(i, j) = datamatrix_jkknf(i, j) / datamatrix_jkknf((i + 1) % Nt, j);
		//}
	}
	stddev = StdDev(datamatrix_jkknf,colon, nbrow);	
	return FitYap(colon, x, stddev , fitarashwl,"Average Mass",true);
    /**************************************************/
}

///  Vector Sorting
/**  For fitting and finding median */
bool 
trQCD::sortVector(vector<Double_t>& realM, vector<Double_t>& M){
    /**************************************************/
    M=realM;
    //bubble sort algorithm
    for(int i=0; i<M.size(); i++){
        for(int j=i+1; j<M.size(); j++){
            if(M[i]>M[j]){
                Double_t tempr=M[i];
                M[i]=M[j];
                M[j]=tempr;
            }        
        }
    }
    return true;
    /**************************************************/
}

///     Clear bad fits
/**  this codes looks upper and lover percent band of the median and clear remaining from analysis vector */
vector <Double_t> 
trQCD::clear_bad_fits(vector<Double_t>& M, vector<Int_t>& bad_fits, Double_t percent){
    /**************************************************/
    Double_t median=Median(M);
    vector <Double_t> cleaned;
    Double_t upper;
    Double_t lower;
    cout<<"Median is : "<< median<<endl;
	Double_t fix=0.0;
    fix =M.size() - 1;
    fix =(fix) / M.size();
    if(median>0){
        upper =  fix*median+(percent*median/100); 
        lower =  fix*median-(percent*median/100); 
    }else {
        upper=  fix*median-(percent*median/100); 
        lower=  fix*median+(percent*median/100); 
    }
    for(int i=0; i<M.size(); i++){
        if((M[i]<upper) and (M[i]>lower))
		{
            cleaned.push_back(M[i]);
		}else {
			bad_fits.push_back(i);
		}
    }
	if(cleaned.size()<bad_fits.size()){
		cerr<<"Too many bad fits the analysis might be wrong!"<<endl;
		exit (EXIT_FAILURE);
	}
    return cleaned;
    /**************************************************/
}

///     Find median
/**  Find median */
Double_t 
trQCD::Median(vector<Double_t>& M){
    /**************************************************/
    std::vector<Double_t> TempM;
    if(sortVector(M,TempM)){
        if((TempM.size()%2)==0){
            Int_t temp= (TempM.size())/2;
            return TempM[temp];
        }else {
            Int_t temp= (TempM.size())/2;
            return (TempM[temp]+TempM[temp+1])/2;
        }    
    }
    /**************************************************/
}

/// Remove the selected column of matrix.
/** Remove the selected column of matrix using the TMatrixD   */
TMatrixD
trQCD::RemoveColumn(TMatrixD x1, Int_t column) {
	const Int_t nrow1 = x1.GetNrows();
	const Int_t ncol1 = x1.GetNcols();

	TMatrixD x2(nrow1,ncol1-1);
	TMatrixDSub(x2,0,nrow1-1,0,column-1) = TMatrixDSub(x1,0,nrow1-1,0,column-1);
	TMatrixDSub(x2,0,nrow1-1,column,ncol1-2) = TMatrixDSub(x1,0,nrow1-1,column+1,ncol1-1);
	return x2;
}

/// Remove the selected columns of matrix.
/** Remove the selected columns of matrix in given vector   */
TMatrixD
trQCD::RemoveColumns(TMatrixD x1, vector<Int_t> to_be_cleaned) {
	std::vector<TMatrixD> tempmatrix; 

	Int_t sayi = to_be_cleaned.size()-1;
	for(int i=0; i<to_be_cleaned.size(); i++) {
		if(i==0){			
			tempmatrix.push_back(RemoveColumn(x1,to_be_cleaned[sayi-i]));			
		}else{
			tempmatrix.push_back(RemoveColumn(tempmatrix[i-1],to_be_cleaned[sayi-i]));
		}
	}
	return tempmatrix[sayi];
}

/// Remove the selected data.
/** Remove the selected columns of vector in given vector   */
vector<Double_t>
trQCD::RemoveColumns(vector<Double_t> data, vector<Int_t> to_be_cleaned) {
	std::vector<Double_t> temp; 
	temp=data;
	int sayi = to_be_cleaned.size()-1;
        for(int i=0; i<to_be_cleaned.size(); i++) {
		int say=to_be_cleaned[sayi-i];
		temp.erase(temp.begin()+say);
	}
	return temp;
}


/// 	Good old fella bootstrap
/**  bootstrap the good old fella */
void
trQCD::bootstrap(TMatrixD G, TMatrixD & G_bootstrap)
{
    /**************************************************/
	TRandom        *r2 = new TRandom();
	int		N_cf = G.GetNcols();
	for (int i = 0; i < N_cf; i++) {
		int		alpha = int (r2->Uniform(0, N_cf));
		for (int j = 0; j < Nt; j++) {
			G_bootstrap(j, i) += G(j, alpha);
		}
	}
    /**************************************************/
}



///  Levi civita function for Magnetic form factor
/**  levi civita function  */
int
trQCD::levi_civita (int i, int j, int k){
    /**************************************************/
    if((i==j)||(i==k)||(j==k)) return 0;
    else return ( (i-j)*(j-k)*(k-i)/2 );
    /**************************************************/
}




