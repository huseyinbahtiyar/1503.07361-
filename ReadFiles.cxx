#include "trQCD.h"
/// 	The Ntuple generator
/** Reads the files, fills the raw data in TMatrixDs, generates root file  */
using namespace std;
int main(int argc, char **argv)
	{
	trQCD *analyze= new trQCD();
	TStopwatch watch;
	watch.Start();
    	int 		it =0;
	int 		qsqrlimit=1;
	int		source=0;
	int		sink = 12;
    if(argc==2) {
        it = std::atoi (argv[1]);
     }//if
	 vector <string> Paths2pt;
	 vector < vector <int> > confs;


	
	/************************** vector namings *******************************/
	
	vector<string> isim;
	isim.push_back("omega");

    // Initialize String Array
    string gamma_1[3] = {"1_g5", "g4g5_1",	"g5_1"};
	string gamma_2[3] = {"1_g5", "g4g5_1",	"g5_1"};
	string smearing[3] = {"50","70","100"};
	string smearing2[3] = {"50","70","100"};


	string dosyaadi (isim[it]+"_matrix.root");
	TFile          *hfile = new TFile(dosyaadi.c_str(), "RECREATE");


	for(int qsqr=0; qsqr<qsqrlimit; qsqr++){
		for(int g1=0; g1<sizeof(gamma_1)/sizeof(gamma_1[0]); g1++){
			for(int g2=0; g2<sizeof(gamma_2)/sizeof(gamma_2[0]); g2++){
				for(int s1=0; s1<sizeof(smearing)/sizeof(smearing[0]); s1++){
					string path = ("13700_spec_all/spec/shpt/"+ gamma_1[g1]+"/"+ gamma_2[g2]+ "/"+ smearing[s1]+"/stripped/"); 
					string strpar1_g1_1_g2_1 (isim[it]+"shpt_q"+analyze->itos(qsqr)+ "_s"+smearing[s1]+"_g1" + gamma_1[g1]+"_g2"+ gamma_2[g2]); 
					if(qsqr==0) {
					    	TMatrixD par1_g1_1_g2_1(Nt, 79);
						TMatrixD par1_g1_1_g2_1_2(Nt, 59);
						TMatrixD par1_g1_1_g2_1_3(Nt, 138);
						par1_g1_1_g2_1 = analyze->reader(3510,4290,10, "/DATA_General_ProtonLike_shpt_px0_py0_pz0", path);
						par1_g1_1_g2_1_2 = analyze->reader(4310,4890,10, "/DATA_General_ProtonLike_shpt_px0_py0_pz0", path);
						par1_g1_1_g2_1_3=analyze->CombineMatrices(par1_g1_1_g2_1,par1_g1_1_g2_1_2);
						par1_g1_1_g2_1_3.Write(strpar1_g1_1_g2_1.c_str());
						hfile->Write();
						
					}
					else{
					    TMatrixD par1_g1_1_g2_1(Nt, 63);
						
						std::vector <string >	shpt_2pt;
						int scale = analyze->findProbs(qsqr,"/DATA_General_ProtonLike_shpt", shpt_2pt);
						par1_g1_1_g2_1 = analyze->readerAveraged(3510,4130,10, qsqr, shpt_2pt, scale, path);
						par1_g1_1_g2_1.Write(strpar1_g1_1_g2_1.c_str());
						hfile->Write();	
								
					}					
				}//s1
			}//gammas2
		}//gammas
		
		/*for(int g1=0; g1<sizeof(gamma_1)/sizeof(gamma_1[0]); g1++){
			for(int g2=0; g2<sizeof(gamma_2)/sizeof(gamma_2[0]); g2++){
				for(int s1=0; s1<sizeof(smearing)/sizeof(smearing[0]); s1++){
						for(int s2=0; s2<sizeof(smearing2)/sizeof(smearing2[0]); s2++){
							cout<<g1<<" "<<g2<<endl;
							string path = ("shsh/"+ gamma_1[g1]+"/"+ gamma_2[g2]+ "/"+ smearing[s1]+"/"+ smearing2[s2]+"/stripped/"); 
							string strpar1_g1_1_g2_1 (isim[it]+"shsh_q"+analyze->itos(qsqr)+ "_g1" + gamma_1[g1]+"_g2"+ gamma_2[g2]); 
							if(qsqr==0) {
  							  	
								TMatrixD par1_shsh_all(Nt, 43);						
								
								par1_shsh_all = analyze->reader(1200,2040,20, "/DATA_General_ProtonLike_shsh_px0_py0_pz0", path);
								par1_shsh_all.Write(strpar1_g1_1_g2_1.c_str());
								hfile->Write();
							
							}
							else{
								TMatrixD par1_shsh_all(Nt, 43);						
							
								std::vector <string >	shpt_2pt;
								int scale = analyze->findProbs(qsqr,"/DATA_General_ProtonLike_shsh", shpt_2pt);
								par1_shsh_all = analyze->readerAveraged(1200,2040,20, qsqr, shpt_2pt, scale, path);
								par1_shsh_all.Write(strpar1_g1_1_g2_1.c_str());
								hfile->Write();	
									
							}					
						}//s2
					}//s1
				
				}//gammas2
			}//gammas
		*/
		
		
		
		hfile->Write();
		
	}//qsqr
	
	
	hfile->Close();

}//main

