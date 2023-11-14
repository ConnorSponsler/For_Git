#include "Constants.c"


using namespace nusquids;


class nuSQUIDSNSI: public nuSQUIDS {
  public:
    squids::SU_vector NSI;
    squids::SU_vector SI_NC;
    squids::SU_vector SI_CC;
    squids::SU_vector custom_HI;
    std::vector<squids::SU_vector> NSI_evol;
    std::vector<squids::SU_vector> SI_NC_evol;
    std::vector<squids::SU_vector> SI_CC_evol;
    double HI_prefactor;
    double DMP_prefactor;
    void AddToPreDerive(double x){
      for(int ei = 0; ei < ne; ei++){
        //asumming same hamiltonian for neutrinos/antineutrinos
        //SU_vector h0 = H0(E_range[ei],0);
        //NSI_evol[ei] = NSI.Evolve(h0,(x-Get_t_initial()));
        NSI_evol[ei] = NSI.Evolve(H0_array[ei],(x-Get_t_initial()));
        SI_NC_evol[ei] = SI_NC.Evolve(H0_array[ei],(x-Get_t_initial()));
        SI_CC_evol[ei] = SI_CC.Evolve(H0_array[ei],(x-Get_t_initial()));
      }
    }
    

        
    squids::SU_vector HI(unsigned int ei,unsigned int index_rho) const{ 
      double CC = HI_prefactor*current_density*current_ye;
      double NC= -0.5*HI_prefactor*current_density*(1.0-current_ye);
      double BSMC = DMP_prefactor*current_density;    
      
      return BSMC*NSI_evol[ei]+CC*SI_CC_evol[ei]+NC*SI_NC_evol[ei];
    }
    
    double GetCurrentDensity() const {
        return HI_prefactor*current_density*current_ye;
    }
    
squids::SU_vector H0(double Enu, unsigned int irho) const{
  return nuSQUIDS::H0(Enu,irho);
}
  public:
   
  nuSQUIDSNSI(marray<double,1> Erange, unsigned int numneu, NeutrinoType NT,
              bool iinteraction) : nuSQUIDS(Erange,numneu,NT,iinteraction) 
  {    
  
  //Initialize all mixing angles to 0 to bypass nusquids' initialization of the first 3 mass states to active
  for(int i=0;i<numneu;i++){
    for(int j=0;j<numneu;j++){
      if(i<j){
        Set_MixingAngle(i,j,0);
      }
    }    
  }


    DMP_prefactor = BSMPotBase/DensityBase;
    // defining a complex matrix M which will contain our BSM Potential
    gsl_matrix_complex * M = gsl_matrix_complex_calloc(numneu,numneu);
    gsl_complex c1 {{ quasi1*1 , 0 }};
    gsl_complex c2 {{ quasi2*1 , 0 }};
    gsl_complex c3 {{ quasi3*1 , 0 }};

 
    HI_prefactor = params.sqrt2*params.GF*params.Na*pow(params.cm,-3); 
    
    gsl_matrix_complex * NCM = gsl_matrix_complex_calloc(numneu,numneu);
    gsl_matrix_complex * CCM = gsl_matrix_complex_calloc(numneu,numneu);
    gsl_complex d {{ 1 , 0.0 }};
    
    
    gsl_matrix_complex_set(NCM, numneu-3 , numneu-3 ,d);
    gsl_matrix_complex_set(NCM, numneu-2 , numneu-2 ,d);
    gsl_matrix_complex_set(NCM, numneu-1 , numneu-1 ,d);
    gsl_matrix_complex_set(CCM, numneu-3 , numneu-3 ,d);
    
    Set_MixingAngle( numneu-3 , numneu-2 ,Th12);
    Set_MixingAngle( numneu-3 , numneu-1 ,Th13);
    Set_MixingAngle( numneu-2 , numneu-1 ,Th23);
    
    if(numneu==4){
      gsl_matrix_complex_set(M,0,0,c2);    
      Set_MixingAngle(0,2,QTh2);
    }
    
    if(numneu==5){    
      if(quasi1==1){
        gsl_matrix_complex_set(M,0,0,c1);
        gsl_matrix_complex_set(M,1,1,c2);    
        Set_MixingAngle(0,2,QTh1);
        Set_MixingAngle(1,3,QTh2);
        
      }else if(quasi3==1){        
        gsl_matrix_complex_set(M,0,0,c2);
        gsl_matrix_complex_set(M,1,1,c3);
        Set_MixingAngle(0,3,QTh2);
        Set_MixingAngle(1,4,QTh3);        
      }
    }
    
    if(numneu==6){
      Set_MixingAngle(1,4,QTh2);
      gsl_matrix_complex_set(M,1,1,c2);        
      if(quasi!=0){
        gsl_matrix_complex_set(M,0,0,c1);
        gsl_matrix_complex_set(M,2,2,c3);
        Set_MixingAngle(0,3,QTh1);
        Set_MixingAngle(2,5,QTh3); 
      }
    }
    
        
    NSI = squids::SU_vector(M);
    SI_NC = squids::SU_vector(NCM);
    SI_CC = squids::SU_vector(CCM);
    
    // rotate to mass reprentation
    NSI.RotateToB1(params);
    SI_NC.RotateToB1(params);
    SI_CC.RotateToB1(params);
    
    NSI_evol.resize(ne);
    SI_NC_evol.resize(ne);
    SI_CC_evol.resize(ne);
    
    for(int ei = 0; ei < ne; ei++){
      NSI_evol[ei] = squids::SU_vector(nsun);
      SI_NC_evol[ei] = squids::SU_vector(nsun);
      SI_CC_evol[ei] = squids::SU_vector(nsun);

    }
    gsl_matrix_complex_free(M);
    gsl_matrix_complex_free(NCM);
    gsl_matrix_complex_free(CCM);

  }

  
  
};
