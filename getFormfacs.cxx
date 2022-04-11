#include <stdlib.h>
#include "Riostream.h"
#include "TFile.h"
#include "TObject.h"
#include "TNetFile.h"
#include "TRandom.h"
#include "TTree.h"
#include "TBranch.h"
#include "TClonesArray.h"
#include "TStopwatch.h"
#include "TCanvas.h"
#include "TMatrix.h"
#include "TMatrixDEigen.h"
#include <TVector.h>
#include "TMatrixDSym.h"
#include "TMatrixDSymEigen.h"
#include "TVectorDfwd.h"
#include "trQCD.h"
#include <map>
#include <sstream>
using namespace std;
int main(int argc, char **argv)
	{
	trQCD *analyze= new trQCD();
	TStopwatch watch;
	watch.Start();
    int it =0;
	int u=0;
	int qsqrlimit=9;
	int		fitara_nuc [2];
	double q = TMath::TwoPi()/Nl;
	Double_t charge = 0.0;
	int		source=0;
	int		sink = 12;
	/************************** vector namings *******************************/
	Double_t Ks=0.13640;
	Double_t Kc = 0.1246;
	Double_t Kcrit=0.137868;

	vector<string> isim;
	isim.push_back("omega_exct");
	
	std::vector <bool> diagonal;
	diagonal.push_back(false);
	
	std::vector <int> fit_bas;
	fit_bas.push_back(3);

	
	std::vector <int> fit_bit;
	fit_bit.push_back(5);


    std::vector <int> error_percent;
	error_percent.push_back(90);
	
	vector<string> FitForm;
	FitForm.push_back("Monopole"); 
	FitForm.push_back("Dipole"); 	
	FitForm.push_back("Exponential"); 
	

    if((argc>1)) {
        it = std::atoi (argv[1]);
		u = std::atoi (argv[2]);
	if(it>isim.size()){
		cout<<"usage: ./getFormfacs int baryon(";
		for(int bars=0; bars<isim.size(); bars++)
	     cout<<bars<<":"<<isim[bars];
		cout<<") int lightquark (0:dsc 1:usc)"<<endl;
	     exit(1);
	}
     }else {
 		cout<<"usage: ./getFormfacs int baryon(";
 		for(int bars=0; bars<isim.size(); bars++)
 	     cout<<bars<<":"<<isim[bars];
 		cout<<") int lightquark (0:dsc 1:usc)"<<endl;
 	     exit(1);
     }
	if(u==1){
		charge = 1.0;
	}
	/************************** vector namings *******************************/
	vector<TMatrixD> GE_par1_qsqr;
	vector<TMatrixD> GM_par1_qsqr;
	std::vector <Double_t> xvec;
	std::vector <Double_t> avcalqsqr;
	
	std::vector<Double_t> mass_par1;
	std::vector<Double_t> mass_par2;
	std::vector <Double_t> GM_par1,GE_par1;
	std::ostringstream dosyaadi2;
	if(u==1){
		dosyaadi2<<isim[it]<<"_usc_sonuc.root";
	}else {
		dosyaadi2<<isim[it]<<"_dsc_sonuc.root";
	}
	std::string dosyaadi = dosyaadi2.str();
	string okuncak (isim[it]+"_matrix.root");
	TFile   *f = new TFile(okuncak.c_str());
	TFile	*hfile = new TFile(dosyaadi.c_str(), "RECREATE");				
	/**************	Mass Calculations	**************/

   	string strpar1_g1_1_g2_1 (isim[it]+"par1_g1_1_g2_1"); 
   	string strpar1_g1_1_g2_g5(isim[it]+"par1_g1_1_g2_g5"); 
   	string strpar1_g1_1_g2_g4g5(isim[it]+"par1_g1_1_g2_g4g5"); 
   	string strpar1_g1_g5_g2_1(isim[it]+"par1_g1_g5_g2_1"); 
   	string strpar1_g1_g5_g2_g5 (isim[it]+"par1_g1_g5_g2_g5"); 
   	string strpar1_g1_g5_g2_g4g5 (isim[it]+"par1_g1_g5_g2_g4g5"); 
   	string strpar1_g1_g4g5_g2_1 (isim[it]+"par1_g1_g4g5_g2_1");  
   	string strpar1_g1_g4g5_g2_g5 (isim[it]+"par1_g1_g4g5_g2_g5"); 
   	string strpar1_g1_g4g5_g2_g4g5(isim[it]+"par1_g1_g4g5_g2_g4g5"); 
		
	TMatrixD* par1_g1_1_g2_1= (TMatrixD*) f->Get(strpar1_g1_1_g2_1.c_str());
	TMatrixD* par1_g1_1_g2_g5= (TMatrixD*) f->Get(strpar1_g1_1_g2_g5.c_str());
	TMatrixD* par1_g1_1_g2_g4g5= (TMatrixD*) f->Get(strpar1_g1_1_g2_g4g5.c_str());
	TMatrixD* par1_g1_g5_g2_1= (TMatrixD*) f->Get(strpar1_g1_g5_g2_1.c_str());
	TMatrixD* par1_g1_g5_g2_g5= (TMatrixD*) f->Get(strpar1_g1_g5_g2_g5.c_str());
	TMatrixD* par1_g1_g5_g2_g4g5= (TMatrixD*) f->Get(strpar1_g1_g5_g2_g4g5.c_str());
	TMatrixD* par1_g1_g4g5_g2_1= (TMatrixD*) f->Get(strpar1_g1_g4g5_g2_1.c_str());
	TMatrixD* par1_g1_g4g5_g2_g5= (TMatrixD*) f->Get(strpar1_g1_g4g5_g2_g5.c_str()); 
	TMatrixD* par1_g1_g4g5_g2_g4g5= (TMatrixD*) f->Get(strpar1_g1_g4g5_g2_g4g5.c_str());
	
	
	int		nbcol = par1_g1_1_g2_1->GetNcols();
	int 		nbfiles= par1_g1_1_g2_1->GetNcols();	
		

	TMatrixD par1_g1_1_g2_1_jack(Nt,nbfiles);
	TMatrixD par1_g1_1_g2_g5_jack(Nt,nbfiles);
	TMatrixD par1_g1_1_g2_g4g5_jack(Nt,nbfiles);
	TMatrixD par1_g1_g5_g2_1_jack(Nt,nbfiles);
	TMatrixD par1_g1_g5_g2_g5_jack(Nt,nbfiles);
	TMatrixD par1_g1_g5_g2_g4g5_jack(Nt,nbfiles);
	TMatrixD par1_g1_g4g5_g2_1_jack(Nt,nbfiles);
	TMatrixD par1_g1_g4g5_g2_g5_jack(Nt,nbfiles);
	TMatrixD par1_g1_g4g5_g2_g4g5_jack(Nt,nbfiles);
	
	par1_g1_1_g2_1_jack = analyze->JackKnife(*par1_g1_1_g2_1);
	par1_g1_1_g2_g5_jack = analyze->JackKnife(*par1_g1_1_g2_g5);
	par1_g1_1_g2_g4g5_jack = analyze->JackKnife(*par1_g1_1_g2_g4g5);
	par1_g1_g5_g2_1_jack = analyze->JackKnife(*par1_g1_g5_g2_1);
	par1_g1_g5_g2_g5_jack = analyze->JackKnife(*par1_g1_g5_g2_g5);
	par1_g1_g5_g2_g4g5_jack = analyze->JackKnife(*par1_g1_g5_g2_g4g5);
	par1_g1_g4g5_g2_1_jack = analyze->JackKnife(*par1_g1_g4g5_g2_1);
	par1_g1_g4g5_g2_g5_jack = analyze->JackKnife(*par1_g1_g4g5_g2_g5);
	par1_g1_g4g5_g2_g4g5_jack = analyze->JackKnife(*par1_g1_g4g5_g2_g4g5);
	
	vector <TMatrixD> Corrmat_t;
	vector <TMatrixD> Corrmat_dt;
	Corrmat_t = analyze->Correlation_Matrix(par1_g1_1_g2_1_jack,par1_g1_1_g2_g5_jack,par1_g1_1_g2_g4g5_jack,par1_g1_g5_g2_1_jack,par1_g1_g5_g2_g5_jack,par1_g1_g5_g2_g4g5_jack,par1_g1_g4g5_g2_1_jack,par1_g1_g4g5_g2_g5_jack,par1_g1_g4g5_g2_g4g5_jack,10);
	Corrmat_dt = analyze->Correlation_Matrix(par1_g1_1_g2_1_jack,par1_g1_1_g2_g5_jack,par1_g1_1_g2_g4g5_jack,par1_g1_g5_g2_1_jack,par1_g1_g5_g2_g5_jack,par1_g1_g5_g2_g4g5_jack,par1_g1_g4g5_g2_1_jack,par1_g1_g4g5_g2_g5_jack,par1_g1_g4g5_g2_g4g5_jack,12);

       // analyze->Eigen_Mass(Corrmat_t,Corrmat_dt);

	hfile->Write();
	hfile->Close();
	watch.Stop();
			cout <<	"Cpu Time: " << watch.CpuTime() << " Real Time: " << watch.RealTime() << endl;

			return 0;		
}

