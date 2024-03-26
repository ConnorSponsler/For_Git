#include <vector>
#include <iostream>
#include <fstream>
#include <nuSQuIDS/nuSQuIDS.h>
#include "Constants.c"

using namespace nusquids;


int main()
{
  squids::Const units;
  unsigned int numneu = 3;

  double Emin=0.1*units.GeV;
  double Emax=1e6*units.GeV;
  int Ebin=1000;
  double czmin=0;
  double czmax=0.1;
  int czbin=2;


  nuSQUIDSAtm<nuSQUIDS,EmittingEarthAtm> nus(linspace(czmin,czmax,czbin), logspace(Emin,Emax,Ebin), numneu, both,false);
  nus.Set_AtmHeight(60);
  nus.Set_NeutrinoSources(true);
  
  //IMPORTANT
  nus.Set_IncludeOscillations(false);
  //nus.Set_AtmProducedFlavors(1,2);
    
  //Initialize all mixing angles to 0
  for(int i=0;i<numneu;i++){
    for(int j=0;j<numneu;j++){
      if(i<j){
        nus.Set_MixingAngle(i,j,0);
      }
    }    
  }
 
  
  nus.Set_MixingAngle( numneu-3 , numneu-2 ,Th12);
  nus.Set_MixingAngle( numneu-3 , numneu-1 ,Th13);
  nus.Set_MixingAngle( numneu-2 , numneu-1 ,Th23);
  if(numneu == 3){
    nus.Set_SquareMassDifference(1, M2);
    nus.Set_SquareMassDifference(2, M3);  
  }


  if(numneu == 4){
    nus.Set_SquareMassDifference(1, M1-QM2);
    nus.Set_SquareMassDifference(2, M2-QM2);
    nus.Set_SquareMassDifference(3, M3-QM2);  
    nus.Set_MixingAngle(0,2,QTh2);
  }
  
  if(numneu == 5){
    if(quasi1==1){ 
        nus.Set_SquareMassDifference(1, QM2-QM1); 
        nus.Set_SquareMassDifference(2, M1-QM1);  
        nus.Set_SquareMassDifference(3, M2-QM1);
        nus.Set_SquareMassDifference(4, M3-QM1);    
        nus.Set_MixingAngle(0,2,QTh1);
        nus.Set_MixingAngle(1,3,QTh2);
        
      }else if(quasi3==1){        
        nus.Set_SquareMassDifference(1, QM3-QM2);  
        nus.Set_SquareMassDifference(2, M1-QM2);  
        nus.Set_SquareMassDifference(3, M2-QM2);
        nus.Set_SquareMassDifference(4, M3-QM2); 
        nus.Set_MixingAngle(0,3,QTh2);
        nus.Set_MixingAngle(1,4,QTh3);        
      }
  
  }
  
  if(numneu==6){      
    nus.Set_SquareMassDifference(1, QM2-QM1);
    nus.Set_SquareMassDifference(2, QM3-QM1); 
    nus.Set_SquareMassDifference(3, M1-QM1);  
    nus.Set_SquareMassDifference(4, M2-QM1);
    nus.Set_SquareMassDifference(5, M3-QM1);  
    nus.Set_MixingAngle(1,4,QTh2);
    nus.Set_MixingAngle(0,3,QTh1);
    nus.Set_MixingAngle(2,5,QTh3); 
  }
  
  //We set the GSL step function
  nus.Set_GSL_step(gsl_odeiv2_step_rk4);

  //Setting the numerical precision of gsl integrator.
  nus.Set_rel_error(1.0e-6);
  nus.Set_abs_error(1.0e-6);

  //Set true the progress bar during the evolution.
  nus.Set_ProgressBar(true);

  marray<double,1> E_range = nus.GetERange();
  //Array that contains the initial state of the system, fist component is energy and second every one of the flavors
  marray<double,4> inistate{nus.GetNumCos(),nus.GetNumE(),2,numneu};
  std::fill(inistate.begin(),inistate.end(),0);
  nus.Set_initial_state(inistate,flavor);

  nus.EvolveState();

  std::ofstream file("fluxes_from_60.txt");


  int Nen=1000;
  int Ncz=100;
  double lEmin=log10(Emin);
  double lEmax=log10(Emax);;


  //Writing to the file!  
  for(double cz=czmin;cz<czmax;cz+=(czmax-czmin)/(double)Ncz){
    for(double lE=lEmin; lE<lEmax; lE+=(lEmax-lEmin)/(double)Nen){
      double row_sum = 0;
      double row_values[3] = {0,0,0};
      double E=pow(10.0,lE);
      file << cz << " " << E;
      for(int fl=0; fl<numneu; fl++){
        double curr_flav = nus.EvalFlavor(fl,cz, E, 0);
        row_sum += curr_flav;
        row_values[fl] = curr_flav;
      }
      //remove "/row_sum for non-normalized-flux"
      file << " " <<  row_values[0]/row_sum << " " << row_values[1]/row_sum << " " << row_values[2]/row_sum;
      file << std::endl;
    }
  }

  return 0;
}
