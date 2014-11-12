// Macro to reconstruct List Mode Data from GATE simulations in ROOT files
//
// 29 April 2014. V0.0 J.L. Herraiz and A. Sitek
// 30 Oct 2014,   V0.1
// --------------------------------------------------------------------
// This code was developed with the goal to be a useful tool for GATE users
//  and researchers working on PET image reconstruction
// Please contact the authors if you want to help in its development
// herraiz@mit.edu 

#include <stdio.h>
#include <math.h>
#include <iostream>
using namespace std;
#include <cstdlib>
#include <fstream>
#include <iomanip>

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TDirectory.h>

//---------- VARIABLES -------------------
Int_t    niter=99;
Float_t  Scanner = 800.0;                    // Default Diameter of the scanner (mm) (it is adapted automatically based on the root file)
Float_t  FOV_Z=130.0;                        // Default Axial FOV (Obtained later based on the root file)
Float_t  FOV_XY=250.0;                       // Size of Reconstructed Image (mm)
Int_t    RES=75, NZS=61;                     // Image Dimensions
Int_t    nvoxels=RES*RES*NZS;
Int_t*   IMG = (Int_t*)malloc(nvoxels*sizeof(Int_t));
Float_t* IMG_N = (Float_t*)malloc(nvoxels*sizeof(Float_t));
Float_t* SENS = (Float_t*)malloc(nvoxels*sizeof(Float_t));
Float_t  dxy = float(RES)/FOV_XY;
Float_t  dz = float(NZS)/FOV_Z;             // Default value
Float_t  RESM = 0.5*RES, NZM=0.5*NZS;
Int_t    ncoinc;
Float_t  source_X, source_Y, source_Z;
Float_t  det1_X, det1_Y, det1_Z;
Float_t  det2_X, det2_Y, det2_Z;
Int_t    b1, b2, c1, c2;
Int_t    ix,iy,iz,iV;
Int_t    indx;

//--------- FUNCTIONS ------------------
void define_branches(TTree *nt);
Int_t Arek(Int_t posic, Float_t det1_X,Float_t det1_Y,Float_t det1_Z,Float_t det2_X,Float_t det2_Y,Float_t det2_Z, Int_t* IMG, Float_t* SENS);
Float_t* MedianIMG(Float_t* IMG);
Float_t* normalization(TString namefile2, int nt);
Float_t* FiltNorm(Float_t* IMG, int nt);

//--------------------------------------
void FBI3D(TString namefile="", TString namefile2=""){

  Int_t    MaxNCoinc = 9e7;
  Int_t*   POS = (Int_t*)malloc(MaxNCoinc*sizeof(Int_t));
  FILE* fichero;
  FILE* fichero2;

  cout << namefile<<" "<<namefile2<<endl;

// === STEP 1 === NORMALIZATION ================

if (namefile2=="normf.raw"){   // Allows reusing previous normalization image
  cout << "Using precalculated Normalization file" << endl;
  fichero2 = fopen(namefile2,"rb");
  fread(SENS,sizeof(float),nvoxels,fichero2);
  fclose(fichero2); 
}else{          // Create normalization from ROOT file (filtered to reduce noise)  
  int nt = 10;  // Number of filter iterations
  SENS = normalization(namefile2,nt);
  // Storing the normalization image so that it can be used for other reconstructions of this scanner
  fichero2 = fopen("normf.raw","wb");    
  fwrite(SENS,sizeof(float),nvoxels,fichero2);
  fclose(fichero2);
}

// === STEP 2 === READ ROOT FILE  ======

 int posic; 
 TFile *f1 = new TFile(namefile,"READ");
 TTree *Coincidences = (TTree*)f1->Get("Coincidences"); 
 ncoinc = Coincidences->GetEntries();
 cout << "Simulated Coincidences = " << ncoinc << endl;
 define_branches(Coincidences);
 Scanner = (Coincidences->GetMaximum("globalPosX1")) - (Coincidences->GetMinimum("globalPosX1")); // Diameter of the scanner from the detected counts
 FOV_Z = (Coincidences->GetMaximum("globalPosZ1")) - (Coincidences->GetMinimum("globalPosZ1"));
 dz = float(NZS)/FOV_Z;  

// === STEP 3 === REFERENCE IMAGE (OPTIONAL) ====

for(Int_t iV=0; iV<nvoxels ; iV++){IMG_N[iV]=0.;}  // Initial Image
for(Int_t ic = 0; ic <= ncoinc; ic++){     // 
  Coincidences->GetEntry(ic);
  ix = int(source_X*dxy + RESM);
  iy = int(source_Y*dxy + RESM);
  iz = int(source_Z*dz + NZM);      
  if (ix>=0 && ix<RES && iy>=0 && iy<RES && iz>=0 && iz<NZS){   // Check voxel valid
   iV = iz*RES*RES+iy*RES+ix;
   IMG_N[iV]++;
  }
 }

fichero = fopen("image_ref.raw","wb");
fwrite(IMG_N,sizeof(float),nvoxels,fichero);
fclose(fichero);

Float_t sens_corr = 0.;
for(Int_t iV=0; iV<nvoxels ; iV++){
 sens_corr = SENS[iV];
 if (sens_corr>0){sens_corr = 1.0/sens_corr;}
 IMG_N[iV]*=sens_corr;
}  // Reference image corrected by sensitivity
fichero = fopen("image_ref_senscor.raw","wb");
fwrite(IMG_N,sizeof(float),nvoxels,fichero);
fclose(fichero);

// === STEP 4 === INITIAL IMAGE ================
for(Int_t iV=0; iV<nvoxels ; iV++){IMG[iV]=0;}  // Initial Image
for(Int_t iV=0; iV<ncoinc; iV++){POS[iV]=-1;}  // Initial Position

// === STEP 5 === ITERATIONS ==================
for(Int_t iter=1; iter<=niter ; iter++){        // Iterations (1 iteration consists of a whole loop over all data)
for(Int_t ic = 0;  ic<ncoinc;  ic++){     
 Coincidences->GetEntry(ic); 
 if ( ic % 100000 == 0 ) cout << "Iter= " << iter << " Processing " << ic << " coincidences " << endl;  
  posic = POS[ic]; 
  POS[ic] = Arek(posic,det1_X,det1_Y,det1_Z,det2_X,det2_Y,det2_Z,IMG,SENS);
 } // COINCIDENCES
} // ITERATIONS

f1->Close();

fichero = fopen("image.raw","wb");
if(fichero){
 fwrite(IMG,sizeof(int),nvoxels,fichero);
 fclose(fichero);
}

// === STEP 6 === FINAL SENSITIVY NORMALIZATION ====
for(Int_t iV=0; iV<nvoxels ; iV++){
  if (SENS[iV]>0){
    IMG_N[iV] = float(IMG[iV])/SENS[iV];
  }else{
    IMG_N[iV] = 0.;
  }
}  

fichero = fopen("image_norm.raw","wb");
if(fichero){
 fwrite(IMG_N,sizeof(float),nvoxels,fichero);
 fclose(fichero);
}

IMG_N = MedianIMG(IMG_N);
fichero = fopen("image_norm_med.raw","wb");
if(fichero){
 fwrite(IMG_N,sizeof(float),nvoxels,fichero);
 fclose(fichero);
}

}  // END

//-----------------------------------------------------
void define_branches(TTree *nt){
  nt->SetBranchAddress("sourcePosX1",&source_X);
  nt->SetBranchAddress("sourcePosY1",&source_Y);
  nt->SetBranchAddress("sourcePosZ1",&source_Z);
  nt->SetBranchAddress("globalPosX1",&det1_X);
  nt->SetBranchAddress("globalPosY1",&det1_Y);
  nt->SetBranchAddress("globalPosZ1",&det1_Z);
  nt->SetBranchAddress("globalPosX2",&det2_X);
  nt->SetBranchAddress("globalPosY2",&det2_Y);
  nt->SetBranchAddress("globalPosZ2",&det2_Z);
  nt->SetBranchAddress("blockID1",&b1);
  nt->SetBranchAddress("blockID2",&b2);
  nt->SetBranchAddress("crystalID1",&c1);
  nt->SetBranchAddress("crystalID2",&c2);
}

// --------------------------------------------------------------------
// --------------------------------------------------------------------

Int_t Arek(int ipos1, Float_t det1_X,Float_t det1_Y,Float_t det1_Z,Float_t det2_X,Float_t det2_Y,Float_t det2_Z, Int_t* IMG,Float_t* SENS){            

       Float_t deltax, deltay,deltaz,norm;
       Float_t x0,y0,z0;
       Float_t lambda;     
       Int_t ix,iy,iz;
       Int_t points;       
       Int_t iV;
       Float_t prob,ran;
       Int_t ipos2;
       Float_t val1,val2;    
       Float_t factor = 0.75*FOV_XY/Scanner;  // FOV/Scanner Diameter          
       Float_t sens1,sens2,sens_ratio;
       
       deltax = det2_X-det1_X;
       deltay = det2_Y-det1_Y;
       deltaz = det2_Z-det1_Z;

       lambda = 0.5 + factor*float(rand() % RES - RESM)/RES;   // Centered at 0.5 point between 0..1
    
       x0 = det1_X + lambda*deltax;    // Position of event in mm
       y0 = det1_Y + lambda*deltay;
       z0 = det1_Z + lambda*deltaz;
 
       ix = int(x0*dxy + RESM) + int((rand()%5 -2)/2);
       iy = int(y0*dxy + RESM) + int((rand()%5 -2)/2);
       iz = int(z0*dz + NZM)   + int((rand()%5 -2)/2);       

       if (ix>=0 && ix<RES && iy>=0 && iy<RES && iz>=0 && iz<NZS){   // Check voxel valid	 
        ipos2 = iz*RES*RES+iy*RES+ix;
        if (ipos1<0) { IMG[ipos2]++; ipos1=ipos2; return ipos1;}   // ONLY FOR FIRST PASS
       }else{
	 //	 cout<<"return "<<ix <<"  "<< iy <<"  " << iz << endl; 
	 return ipos1;
       }  	 
                    
       sens1 = SENS[ipos1];
       sens2 = SENS[ipos2];
       if (sens2>0.) { sens_ratio = sens1/sens2;} else { return ipos1;}

       val1 = IMG[ipos1];
       val2 = IMG[ipos2];

       prob =sens_ratio*float(val2+1.0)/float(val1);    
       ran = float(rand() % 1000)/1000.0;
       if (ran < prob) {
	IMG[ipos1]--;
        IMG[ipos2]++;
	ipos1 = ipos2;
       }      

       return ipos1;
}

// ------------------ NORMALIZATION ----------------------------

Float_t* normalization(TString namefile2, int nt){
 
 TFile *f2 = new TFile(namefile2,"READ");
 
 TTree *Coincidences2 = (TTree*)f2->Get("Coincidences"); 
 ncoinc = Coincidences2->GetEntries();
 cout << "Normalization file with " << ncoinc << " coincidences " <<endl;
 define_branches(Coincidences2);
 for(Int_t iV=0; iV<nvoxels ; iV++){IMG_N[iV]=0.;}  // Initial image
 
 for(Int_t ic = 0;  ic < ncoinc  ;  ic++){       // 
  Coincidences2->GetEntry(ic); 
  ix = int(source_X*dxy + RESM);
  iy = int(source_Y*dxy + RESM);
  iz = int(source_Z*dz + NZM);      
  if (ix>=0 && ix<RES && iy>=0 && iy<RES && iz>=0 && iz<NZS){    // Check voxel valid
   iV = iz*RES*RES+iy*RES+ix;
   IMG_N[iV]++;
  }
 }
 
 f2->Close();
 SENS = FiltNorm(IMG_N,nt);  // Optional, to avoid noise from normalization
 return SENS;
}

// -------------------- FILTER IMAGES ----------------------------
Float_t* FiltNorm(Float_t* IMG, int nt){

  Float_t* IMG_NF = (Float_t*)malloc(nvoxels*sizeof(Float_t));
  Int_t II,JJ,KK,iV,iVV;
  int num = 0;
  float sum = 0.;

  for (int it=0;it<nt;it++){

  cout << " Filtering. Iter "<< it<< endl;

  for (int K=0;K<NZS;K++){
  for (int J=0;J<RES;J++){
  for (int I=0;I<RES;I++){

   iV = K*RES*RES+J*RES+I;
   sum = 0.;
   num = 0 ;
   for(int dK = -1; dK <= 1; dK++) {		  
   for(int dJ = -1; dJ <= 1; dJ++) {		  
   for(int dI = -1; dI <= 1; dI++) {
    KK = K + dK;
    JJ = J + dJ;
    II = I + dI;		
    if (KK >=0 && KK<NZS && JJ>=0 && JJ<RES && II>=0 && II<RES) {   
     iVV = KK*RES*RES+JJ*RES+II;
     sum+=IMG[iVV];
     if (IMG[iVV]>0) { num++;}
    }          

   }
   }
   } // FILTER KERNEL'S LOOPS

   if (num>0) {IMG_NF[iV]=sum/num;}else{IMG_NF[iV]=0.;}

 }
 }
 } // VOXELS' LOOPS
 
  for (int iV=0;iV<nvoxels;iV++){IMG[iV]=IMG_NF[iV];}
 } // ITERATIONS' LOOP
  
return IMG;
}

// -------------------- MEDIAN FILTER IMAGES ----------------------------
Float_t* MedianIMG(Float_t* IMG){

  Float_t* IMG_NF = (Float_t*)malloc(nvoxels*sizeof(Float_t));
  Int_t II,JJ,KK,iV,iVV,iP,position;
  int num = 0;
  float sum = 0.;
  float swap =0.;
  int N = 27;
  Float_t* Array = (Float_t*)malloc(N*sizeof(Float_t));

  for (int K=0;K<NZS;K++){
  for (int J=0;J<RES;J++){
  for (int I=0;I<RES;I++){

   iV = K*RES*RES+J*RES+I;
   sum = 0.;
   num = 0 ;
   iP = 0;
   for(int dK = -1; dK <= 1; dK++) {		  
   for(int dJ = -1; dJ <= 1; dJ++) {		  
   for(int dI = -1; dI <= 1; dI++) {
    KK = K + dK;
    JJ = J + dJ;
    II = I + dI;	
    iP++;
    if (KK >=0 && KK<NZS && JJ>=0 && JJ<RES && II>=0 && II<RES) {   
     iVV = KK*RES*RES+JJ*RES+II;	
     Array[iP] = IMG[iVV];
    }else{
     Array[iP] = IMG[iV];
    }
   }
   }
   }

   for (int c = 0 ; c < N-1; c++ ){
    position = c;
    for (int d = c + 1 ; d < N ; d++ ){
     if ( Array[position] > Array[d] ) position = d;
    }
    if ( position != c ) {
     swap = Array[c];
     Array[c] = Array[position];
     Array[position] = swap;
    }
   }

   IMG_NF[iV]=Array[N/2];

 }
 }
 } // VOXELS' LOOPS
 
 for (int iV=0;iV<nvoxels;iV++){IMG[iV]=0.5*IMG_NF[iV] + 0.5*IMG[iV];}
  
return IMG;
}

// ---- MAIN --- ENABLES g++ COMPILATION ---
// http://en.wikibooks.org/wiki/ROOT/Getting_Started/Many_Ways_to_Use_ROOT

# ifndef __CINT__
int main(int argc, char * argv[]){
  TString namefile="";
  TString namefile2="";
  namefile = argv[1];
  namefile2 = argv[2];
  FBI3D(namefile,namefile2);
 return 0;
}
# endif
