#include <vector>
#include <iostream>
#include <fstream>
#include <nuSQuIDS/nuSQuIDS.h>
#include "NSI.h"

using namespace nusquids;


int main()
{
  unsigned int numneu;
  squids::Const units;
  numneu = 6;



  double Emin=6.5e2*units.GeV;
  double Emax=8e2*units.GeV;
  int Ebin=1000;
  double czmin=0;
  double czmax=1;
  int czbin=2;

  //nuSQUIDSAtm<nuSQUIDSNSI, EmittingEarthAtm>  nus(linspace(czmin,czmax,czbin), logspace(Emin,Emax,Ebin), numneu, neutrino,false);
  
  nuSQUIDSNSI nus(logspace(Emin,Emax,Ebin), numneu, neutrino,false);
  //nus.Set_AtmHeight(60);
  //nus.Set_AtmEmission();
  //nus.Set_NeutrinoSources(true);
  
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
  
  
  std::shared_ptr<EarthAtm> earth_atm = std::make_shared<EarthAtm>();
  earth_atm->SetAtmosphereHeight(60);
  std::shared_ptr<EarthAtm::Track> earth_atm_track = std::make_shared<EarthAtm::Track>(earth_atm->MakeTrack(acos(0)));
  
  
  const double layer_1 = 500*units.km; //Baseline corresponds to 20km height in the atmosphere at cos(zen)=0 to the earths surface 
  
  std::shared_ptr<Vacuum> vacuum = std::make_shared<Vacuum>();
  std::shared_ptr<Vacuum::Track> track_vac = std::make_shared<Vacuum::Track>(layer_1);
  
  std::shared_ptr<ConstantDensity> cdens = std::make_shared<ConstantDensity>(0.001225, 0.5);
  std::shared_ptr<ConstantDensity::Track> track_cdens = std::make_shared<ConstantDensity::Track>(layer_1);

  // Then we set the body and trajectory in the nuSQUIDS object.
  nus.Set_Body(earth_atm);
  nus.Set_Track(earth_atm_track);
  //nus.Set_Body(cdens);
  //nus.Set_Track(track_cdens);
  //nus.Set_Body(vacuum);
  //nus.Set_Track(track_vac);

  //We set the GSL step function
  nus.Set_GSL_step(gsl_odeiv2_step_rk4);

  //Setting the numerical precision of gsl integrator.
  nus.Set_rel_error(1.0e-6);
  nus.Set_abs_error(1.0e-6);

  //Set true the progress bar during the evolution.
  nus.Set_ProgressBar(true);

  marray<double,1> E_range = nus.GetERange();
  //Array that contains the initial state of the system, fist component is energy and second every one of the flavors
  //marray<double,3> inistate{nus.GetNumE(),numneu};

  //Array that contains the initial state of the system, fist component is energy and second every one of the flavors
  marray<double,2> inistate{nus.GetNumE(),numneu};

  std::fill(inistate.begin(),inistate.end(),0);
  for ( int i = 0 ; i < inistate.extent(0); i++){
        for ( int j = 0 ; j < inistate.extent(1); j++){
          inistate[i][j] = (j == numneu-2) ? 1 : 0;
        }
      }
  //Set the initial state in nuSQuIDS object
  nus.Set_initial_state(inistate,flavor);
  //Propagate the neutrinos in the earth for the path defined in path

  nus.EvolveState();
  //std::cout << "end evolve state";
  std::ofstream file("fluxes_flavor.txt");
    
int Nen =1000;
  int Ncz=5;
  double lEmin=log10(Emin);
  double lEmax=log10(Emax);
    //Writing to the file!  
  for(double lE=lEmin; lE<lEmax; lE+=(lEmax-lEmin)/(double)Nen){
    double E=pow(10.0,lE);
    file << E;
    for(int fl=0; fl<numneu; fl++){
	    file << " " <<  nus.EvalFlavor(fl,E);
    }
    file << std::endl;
  }


  return 0;
}
