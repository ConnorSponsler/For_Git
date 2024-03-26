#include <vector>
#include <iostream>
#include <fstream>
#include <nuSQuIDS/nuSQuIDS.h>
#include "Constants.c"

#include <iostream>
#include <fstream>
#include <vector>
#include <stdexcept> // For std::runtime_error


// Helper function to read data from a file and populate vectors
void ReadDataFromFile(const std::string& filename,
                      std::vector<double>& x_values,
                      std::vector<double>& y1_values,
                      std::vector<double>& y2_values,
                      std::vector<double>& y3_values) {
    std::ifstream file(filename);
    if (!file) {
        throw std::runtime_error("Cannot open file");
    }

    double x, y1, y2, y3;
    while (file >> x >> y1 >> y2 >> y3) {
        x_values.push_back(x);
        y1_values.push_back(y1);
        y2_values.push_back(y2);
        y3_values.push_back(y3);
    }
}


using namespace nusquids;

double flux_function(double enu, double cz){
  return 0;
}

int main()
{
  unsigned int numneu;
  squids::Const units;
  numneu = 3;


  double Emin=0.1*units.GeV;
  double Emax=1e6*units.GeV;
  int Ebin=1000;
  double czmin=0;
  double czmax=1;
  int czbin=2;
  auto energy_list = logspace(Emin,Emax,Ebin);
  //nuSQUIDSAtm<nuSQUIDSNSI, EmittingEarthAtm>  nus(linspace(czmin,czmax,czbin), logspace(Emin,Emax,Ebin), numneu, neutrino,false);
  
  nuSQUIDS nus(energy_list, numneu, neutrino,false);
  
  //IMPORTANT
  nus.Set_IncludeOscillations(false);
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
  nus.Set_SquareMassDifference(1, M2);
  nus.Set_SquareMassDifference(2, M3);  

  std::shared_ptr<EarthAtm> emit_earth_atm = std::make_shared<EarthAtm>();
  emit_earth_atm->SetAtmosphereHeight(20);
  std::shared_ptr<EarthAtm::Track> earth_atm_track = std::make_shared<EarthAtm::Track>(emit_earth_atm->MakeTrack(acos(0)));
  
  
  const double layer_1 = 0*units.km; //Baseline corresponds to 20km height in the atmosphere at cos(zen)=0 to the earths surface 
  
  std::shared_ptr<Vacuum> vacuum = std::make_shared<Vacuum>();
  std::shared_ptr<Vacuum::Track> track_vac = std::make_shared<Vacuum::Track>(layer_1);
  
  std::shared_ptr<ConstantDensity> cdens = std::make_shared<ConstantDensity>(0.001225, 0.5);
  std::shared_ptr<ConstantDensity::Track> track_cdens = std::make_shared<ConstantDensity::Track>(layer_1);

  // Then we set the body and trajectory in the nuSQUIDS object.
  nus.Set_Body(emit_earth_atm);
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

  std::string filename = "surf_flux.txt"; // Correctly setting the filename
  std::cout << "Starting setting initial flux\n";
  
  
  
  std::vector<double> x, y2, y3, y4;
  ReadDataFromFile(filename, x, y2, y3, y4);
  for (int i=0; i<121; i++){
    std::cout << x[i] << "\n";
  }

  // Directly use vectors with constructor
  AkimaSpline spline2(x, y2);
  AkimaSpline spline3(x, y3);
  AkimaSpline spline4(x, y4);

  for (size_t i = 0; i < inistate.extent(0); i++) {
      auto curr_E = energy_list[i]*1e-9;
      // Adjust the function call according to your definition
      inistate[i][0] = spline2(curr_E); // For 2nd column
      inistate[i][1] = spline3(curr_E); // For 3rd column
      inistate[i][2] = spline4(curr_E); // For 4th column
  }
  std::cout << "Finished setting initial flux\n";
  nus.Set_initial_state(inistate,flavor);
  //Propagate the neutrinos in the earth for the path defined in path

  nus.EvolveState();
  //std::cout << "end  std::cout << "end setting init flux" << std::endl; evolve state";
  std::ofstream file("fluxes_from_20.txt");
    
  int Nen =1000;
  int Ncz=5;
  double lEmin=log10(Emin);
  double lEmax=log10(Emax);
    //Writing to the file!  
  for(double lE=lEmin; lE<lEmax; lE+=(lEmax-lEmin)/(double)Nen){
    double row_sum = 0;
    double row_values[3] = {0,0,0};
    double E=pow(10.0,lE);
    file << E;
    for(int fl=0; fl<numneu; fl++){
            double curr_flav = nus.EvalFlavor(fl,E);
            row_sum += curr_flav;
            row_values[fl] = curr_flav;
    }
    //remove "/row_sum for non-normalized-flux"
    file << " " <<  row_values[0]/row_sum << " " << row_values[1]/row_sum <<" " << row_values[2]/row_sum ;
    file << std::endl;
  }


  return 0;
}
