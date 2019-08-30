//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Aug 12 18:32:20 2019 by ROOT version 6.14/02
// from TChain particles/
//////////////////////////////////////////////////////////

#ifndef particles_h
#define particles_h

#include "Qvector.h"
#include <TMath.h>
#include <TROOT.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <TMath.h>
#include <Rtypes.h>
#include <TString.h>
#include <TH1.h>
#include <TH2.h>
#include <TFile.h>
#include <TProfile.h>
#include <TChain.h>
#include <TString.h>
#include <string>

Float_t Get_bessel_resolution(Float_t a);

const Int_t   qbins = 100;
const Int_t  impmax = 18;
const Float_t minpt = 0.1;

enum usePart
{
  kPROTON = 0,
  kPROTONBAR,
  kNEUTRON,
  kNEUTRONBAR,
  kPIPLUS,
  kPIMINUS,
  kKPLUS,
  kKMINUS,
  kParticles
};

const struct useParticle
{
  Int_t id;
  Int_t pdg;
  Double_t mass;
  Int_t charge;
  TString name;
  TString displayName;
} gParticles[kParticles] =
{
  {kPROTON,       2212,  0.938,  1,  "_p",          "p"},
  {kPROTONBAR,   -2212,  0.938, -1,  "_pbar",       "#bar{p}"},
  {kNEUTRON,      2112,  0.938,  0,  "_n",          "n"},
  {kNEUTRONBAR,  -2112,  0.938,  0,  "_nbar",       "bar{n}"},
  {kPIPLUS,       211,   0.138,  1,  "_piplus",     "#pi^{+}"},
  {kPIMINUS,     -211,   0.138, -1,  "_piminus",    "#pi^{-}"},
  {kKPLUS,        321,   0.494,  1,  "_kplus",      "#Kappa^{+}"},
  {kKMINUS,      -321,   0.494, -1,  "_kminus",     "#Kappa^{-}"},
};




// Header file for the classes stored in the TTree if any.

class particles {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           npart;
   Double_t        impact_b;
   Int_t           ev;
   Int_t           tcounter;
   Int_t           pdgcode[2500];//[npart]
   Double_t        p0[2500];   //[npart]
   Double_t        px[2500];   //[npart]
   Double_t        py[2500];   //[npart]
   Double_t        pz[2500];   //[npart]
   Double_t        t[2500];   //[npart]
   Double_t        x[2500];   //[npart]
   Double_t        y[2500];   //[npart]
   Double_t        z[2500];   //[npart]

   // List of branches
   TBranch        *b_npart;   //!
   TBranch        *b_impact_b;   //!
   TBranch        *b_ev;   //!
   TBranch        *b_tcounter;   //!
   TBranch        *b_pdgcode;   //!
   TBranch        *b_p0;   //!
   TBranch        *b_px;   //!
   TBranch        *b_py;   //!
   TBranch        *b_pz;   //!
   TBranch        *b_t;   //!
   TBranch        *b_x;   //!
   TBranch        *b_y;   //!
   TBranch        *b_z;   //!

   //List of hists
   // Flow

   TProfile *hQxflow[3];			         //3 harm , 2 +-rap;
   TProfile *hQyflow[3];			         //3 harm , 2 +-rap;
   TProfile *hQx[3][2];			         //3 harm , 2 +-rap;
   TProfile *hQy[3][2];			         //3 harm , 2 +-rap;
   TProfile *hcos[3][2];			         //3 from v1 and vn
   TProfile *hrealre[3][2]; 		         //2 from v1 and vn
   TH1F *hpsi[3][2];
   TH1F *hpsi2[3][2];
   TH1F *histQxtest;
   TH1F *histQytest;
   TH1F *histQxrtest;
   TH1F *histQyrtest;
   
   TProfile *pVn_pT[3][2][kParticles];	//3 v1/v2/v3 ; 2 r from  v1 and vn; all particles
   TProfile *pVn_Y[3][2][kParticles];
   TProfile *hpVn_pT[3][2];         //all particles, 3 har, 2 res m==n and m!=n		
   TProfile *hpVn_Y[3][2];



   particles(TTree *tree=0);
   virtual ~particles();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
   virtual void     InitHists();
   virtual void     FillHists();
   virtual void     SaveHists();
   virtual void     HistsDell();
   virtual void     findQ(Float_t _minrap=1, Float_t _maxrap=2);      //make it better
   virtual void     findRes(Float_t _minrap=1, Float_t _maxrap=2);
   virtual void     findrealRes();
   virtual void     findvn(Float_t _minrap=1, Float_t _maxrap=2);
};

#endif

#ifdef particles_cxx
particles::particles(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {

#ifdef SINGLE_TREE
      // The following code should be used if you want this class to access
      // a single tree instead of a chain
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("Memory Directory");
      if (!f || !f->IsOpen()) {
         f = new TFile("Memory Directory");
      }
      f->GetObject("particles",tree);

#else // SINGLE_TREE

      // The following code should be used if you want this class to access a chain
      // of trees.
   TChain * chain = new TChain("particles","");
   std::ifstream in("filelist", ios::in);
	if (!in) {
		return;
	}
	std::string s;
   //int number=0;
	while (!in.eof()) {
		std::getline(in, s);
		TString mystr = s;
		chain->Add(mystr);
     // if(number++>10) break;
	}
	in.close();

   tree = chain;
#endif // SINGLE_TREE

   }
   Init(tree);
}

particles::~particles()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
   HistsDell();
}

Int_t particles::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t particles::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void particles::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("npart", &npart, &b_npart);
   fChain->SetBranchAddress("impact_b", &impact_b, &b_impact_b);
   fChain->SetBranchAddress("ev", &ev, &b_ev);
   fChain->SetBranchAddress("tcounter", &tcounter, &b_tcounter);
   fChain->SetBranchAddress("pdgcode", pdgcode, &b_pdgcode);
   fChain->SetBranchAddress("p0", p0, &b_p0);
   fChain->SetBranchAddress("px", px, &b_px);
   fChain->SetBranchAddress("py", py, &b_py);
   fChain->SetBranchAddress("pz", pz, &b_pz);
   fChain->SetBranchAddress("t", t, &b_t);
   fChain->SetBranchAddress("x", x, &b_x);
   fChain->SetBranchAddress("y", y, &b_y);
   fChain->SetBranchAddress("z", z, &b_z);
   Notify();
}

Bool_t particles::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void particles::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t particles::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
void particles::InitHists(){
   histQxtest = new TH1F("Qxhist", "Qx", 1000, -1, 1);
   histQytest = new TH1F("Qyhist", "Qy", 1000, -1, 1);
   histQxrtest = new TH1F("Qxhistr", "Qxr", 1000, -1, 1);
   histQyrtest = new TH1F("Qyhistr", "Qyr", 1000, -1, 1);
   for (int i=0;i<3;i++){
      hQxflow[i] = new TProfile(Form("hQx_h_%i", i+1), Form("hQx_h_%i;impact, [mm]; Qx",i+1),qbins,0,impmax);
      hQyflow[i] = new TProfile(Form("hQy_h_%i", i+1), Form("hQy_h_%i;impact, [mm]; Qy",i+1),qbins,0,impmax);
      for (int j=0; j<2;j++){
         hQx[i][j]     = new TProfile(Form("hQx_h_%i_s_%i",i+1,j),Form("hQx_harm_%i_s_%i;impact, [mm];hQx", i+1, j), qbins, 0, impmax);
         hQy[i][j]     = new TProfile(Form("hQy_h_%i_s_%i",i+1,j),Form("hQy_harm_%i_s_%i;impact, [mm];hQy", i+1, j), qbins, 0, impmax);
         hcos[i][j]    = new TProfile(Form("cos_h_%i_s_%i",i+1,j),Form("cos_harm_%i_s_%i;impact, [mm];cos", i+1, j), impmax, 0, impmax);
         hrealre[i][j] = new TProfile(Form("res_h_%i_s_%i",i+1,j),Form("res_harm_%i_s_%i;impact, [mm];res", i+1, j), impmax, 0, impmax);
         hpsi[i][j]    = new TH1F(Form("psi_%i_s_%i",i+1,j),Form("#Psi_%i_s_%i;#Psi, [rad]; nEntries", i+1, j),qbins,-3.2,3.2);
         hpsi2[i][j]    = new TH1F(Form("psi2_%i_s_%i",i+1,j),Form("#Psi_%i_s_%i;#Psi before r, [rad]; nEntries", i+1, j),qbins,-3.2,3.2);
         //add hists later.
         hpVn_pT[i][j] = new TProfile(Form("v_{%i}pt_%i", (i+1), j), Form("v_{%i} vs pt_%i; p_{t}, [GeV/c]; v_{%i}", (i+1), j,(i+1) ), 25, 0, 2.5);
         hpVn_Y[i][j]  = new TProfile(Form("v_{%i}Y_%i", (i+1), j), Form("v_{%i} vs Y_%i; Y; v_{%i}", i+1, j,i+1), 20, -1, 1);
         for (int k=0; k<kParticles; k++){
            pVn_pT[i][j][k] = new TProfile(Form("v_%ipt_%i_"+gParticles[k].name,i,j),Form("v_{%i} vs pt_%i_"+gParticles[k].displayName+"; p_{t}, [GeV/c]; v_{%i}",(i+1),j,(i+1)), 25, 0, 2.5 );
            pVn_Y[i][j][k]  = new TProfile(Form("v_%iY_%i_"+gParticles[k].name,i,j),Form("v_{%i} vs Y_%i_"+gParticles[k].displayName+"; rapidity; v_{%i}",(i+1),j,(i+1)), 20, -1, 1 );
         }
      }
   }
}
void particles::FillHists(){
   //empty
}
void particles::SaveHists(){
   //root save
   TFile* w = new TFile("Flows.root", "recreate");
   TDirectory *outputDir;
   w->cd();
   outputDir = w->mkdir("Qvectors");
   outputDir->cd ();
   histQxtest->Write();
   histQytest->Write();
   histQxrtest->Write();
   histQyrtest->Write();
   for (int i=0;i<3;i++){
      hQxflow[i]->Write();
      hQyflow[i]->Write();
      for (int j=0; j<2;j++){
         hQx[i][j]->Write();
         hQy[i][j]->Write();
      }
   }
   outputDir = w->mkdir("Psy");
   outputDir->cd ();
   for (int i=0;i<3;i++){
      for (int j=0;j<2;j++){
         hpsi[i][j]->Write();
         hpsi2[i][j]->Write();
         hcos[i][j]->Write();
         hrealre[i][j]->Write();
      }
   }
   outputDir = w->mkdir("flows");
   outputDir->cd ();
   for (int i=0;i<3;i++){
      for (int j=0;j<2;j++){
         hpVn_pT[i][j]->Write();
         hpVn_Y[i][j]->Write();
      }
   }
   TDirectory *outputDir2;
   for (int k=0; k<kParticles; k++){
      outputDir2 = outputDir->mkdir(gParticles[k].name);
      outputDir2->cd ();
      for (int i=0;i<3;i++){
         for (int j=0;j<2;j++){
            pVn_pT[i][j][k]->Write(); 
            pVn_Y[i][j][k]->Write(); 
         }
      } 
   }
   w->Close();
}

void particles::HistsDell(){
   // empty
}

void  particles::findQ(Float_t _minrap, Float_t _maxrap){   //revrite later
   for (int i=0;i<3;i++){  //i is EP harm -1
      Double_t Q[2]; Q[0]=0; Q[1]=0;
      Float_t Qxsubev[2], Qysubev[2]; Qxsubev[0]=0; Qxsubev[1]=0; Qysubev[0]=0; Qysubev[1]=0;
      Double_t rap;
      Int_t counter[2]={0,0};
      Double_t pt;
      Double_t psi;
      Int_t harm = i+1;
      for (Long64_t j = 0; j < npart; j++){
         Int_t mod;
         rap = 1.0/2.0 * log( (p0[j] + pz[j]) / (p0[j] - pz[j]) );
         if (fabs(rap) < _minrap || fabs(rap) > _maxrap) continue;
         pt = sqrt( px[j]*px[j] + py[j]*py[j] );
         if (pt < minpt) continue;
         psi = TMath::ATan2(py[j],px[j]);
         rap>0?mod=0:mod=1;
         Qxsubev[mod] += pt*cos(harm * psi);
         Qysubev[mod] += pt*sin(harm * psi);
         Q[mod] += pt;
         counter[mod]++;
      }
      if (!counter[0] || !counter[1] ) continue;
      hQxflow[i]->Fill(impact_b, (Qxsubev[0] + Qxsubev[1]) / (Q[0] + Q[1]) );
      hQyflow[i]->Fill(impact_b, (Qysubev[0] + Qysubev[1]) / (Q[0] + Q[1]) );
      for (int j=0; j<2; j++){
         Qxsubev[j] /= Q[j];
         Qysubev[j] /= Q[j];
         hQx[i][j]->Fill(impact_b, Qxsubev[j]);
         hQy[i][j]->Fill(impact_b, Qysubev[j]);
      }
      histQxtest->Fill(Qxsubev[0]);
      histQytest->Fill(Qysubev[0]);
   }
}
void  particles::findRes(Float_t _minrap, Float_t _maxrap){
   for (int i=0;i<3;i++){
      Double_t Q[2]; Q[0]=0; Q[1]=0;
      Float_t Qxsubev[2], Qysubev[2]; Qxsubev[0]=0; Qxsubev[1]=0; Qysubev[0]=0; Qysubev[1]=0;
      Double_t rap;
      Int_t counter[2]={0,0};
      Double_t pt;
      Int_t harm = i+1;
      for (Long64_t j = 0; j < npart; j++){
         Double_t psi;
         Int_t mod;
         rap = 1.0/2.0 * log( (p0[j] + pz[j]) / (p0[j] - pz[j]) );
         if ( fabs(rap) < _minrap || fabs(rap) > _maxrap) continue;
         pt = sqrt( px[j]*px[j] + py[j]*py[j] );
         if (pt < minpt) continue;
         psi = TMath::ATan2(py[j],px[j]);
         rap>0?mod=0:mod=1;
         Qxsubev[mod] += pt*cos(harm * psi);
         Qysubev[mod] += pt*sin(harm * psi);
         Q[mod] += pt;
         counter[mod]++;
      }
      if (!counter[0] || !counter[1]) continue;
      for (int j=0; j<2; j++){
         Qxsubev[j] /= Q[j];
         Qysubev[j] /= Q[j];
      }
      Float_t psi[2], curcos;
      psi[0] = 1.0 / harm * TMath::ATan2(Qysubev[0], Qxsubev[0]);
      psi[1] = 1.0 / harm * TMath::ATan2(Qysubev[1], Qxsubev[1]);
      hpsi2[i][0]->Fill( psi[0] );
      hpsi2[i][1]->Fill( psi[1] );
      Qxsubev[0] -= hQx[i][0]->GetBinContent((int) (impact_b * qbins) / impmax + 1);
      Qysubev[0] -= hQy[i][0]->GetBinContent((int) (impact_b * qbins) / impmax + 1);
      Qxsubev[1] -= hQx[i][1]->GetBinContent((int) (impact_b * qbins) / impmax + 1);
      Qysubev[1] -= hQy[i][1]->GetBinContent((int) (impact_b * qbins) / impmax + 1);
      psi[0] = 1.0 / harm * TMath::ATan2(Qysubev[0], Qxsubev[0]);
      psi[1] = 1.0 / harm * TMath::ATan2(Qysubev[1], Qxsubev[1]);
      curcos = cos( harm * (psi[1] - psi[0]) );
      hpsi[i][0]->Fill( psi[0] );
      hpsi[i][1]->Fill( psi[1] );
      hcos[i][0]->Fill( impact_b, curcos);
      if (i==0){
         curcos = cos( 2 * (psi[1] - psi[0]) );
         hcos[1][1]->Fill( impact_b, curcos );
         curcos = cos( 3 * (psi[1] - psi[0]) );
         hcos[2][1]->Fill( impact_b, curcos );
      }
      histQxrtest->Fill(Qxsubev[0]);
      histQyrtest->Fill(Qysubev[0]);
   } 
}
void  particles::findrealRes(){
   Int_t NbinsX = hcos[0][0]->GetNbinsX();
   Float_t widthofbin = hrealre[0][0]->GetXaxis()->GetBinWidth(1);
   Float_t res[3][2], err[3][2];
   for (int j=0; j<3; j++) {     // EP harm
      for (int k=0; k<2; k++) {  // flow harm
         if (j==0 && k==1) continue;
         for (int iBin=0; iBin<NbinsX; iBin++){
            res[j][k] = hcos[j][k]->GetBinContent(iBin + 1);
            err[j][k] = hcos[j][k]->GetBinError(iBin + 1);
            Int_t usingharm=1;
            if (k==1){
               usingharm++;
               if (j==2){
                  usingharm++;
               }
            }
           // res[j][k] = sqrt(Get_Res(Get_Chi(abs(res[j][k]),usingharm, 100)));  
            res[j][k] = Get_bessel_resolution_ap(abs(res[j][k]), usingharm);
            err[j][k] = Get_error_res2sub(abs(res[j][k]),err[j][k]);
            hrealre[j][k]->Fill( widthofbin * (iBin+0.5), res[j][k]+err[j][k] );
            hrealre[j][k]->Fill( widthofbin * (iBin+0.5), res[j][k]-err[j][k] );
         }
      }
   }
}
void  particles::findvn(Float_t _minrap, Float_t _maxrap) {
   for (int iHarm=0; iHarm<3; iHarm++){
      Float_t res;
      Double_t rap;
      Double_t pt;
      Int_t counter=0;
      res = hrealre[iHarm][0]->GetBinContent((int)impact_b + 1);
      Float_t PsiEP, _Qx=0, _Qy=0, _Q=0;
      // /*
      for (Long64_t i = 0; i < npart; i++){
         Float_t psi;
         rap = 1.0/2.0 * log( (p0[i] + pz[i]) / (p0[i] - pz[i]) );
         if (fabs(rap) < _minrap || fabs(rap) > _maxrap) continue;
         pt = sqrt( px[i]*px[i] + py[i]*py[i] );
         if (pt < minpt) continue;
         psi = TMath::ATan2(py[i],px[i]);
         _Qx += pt*cos( (iHarm+1) * psi);
         _Qy += pt*sin( (iHarm+1) * psi);
         _Q += pt;
         counter++;
      }
      if (counter < 2) continue;
      _Qx = _Qx/_Q;
      _Qy = _Qy/_Q;
      //_Qx -= hQxflow[iHarm]->GetBinContent(  (int) (impact_b * qbins) / impmax + 1  );
      //_Qy -= hQyflow[iHarm]->GetBinContent(  (int) (impact_b * qbins) / impmax + 1  );
      PsiEP = 1 / (iHarm + 1) * TMath::ATan2(_Qy,_Qx);   //EP from harm iHarm
      //*/
      PsiEP = 0;
      Float_t psy;
      Float_t vn;
      for (Long64_t i = 0; i < npart; i++){
         rap = 1.0/2.0 * log( (p0[i] + pz[i]) / (p0[i] - pz[i]) );
         if (fabs(rap)>1.0) continue;
         Int_t ppid = -1;
         for (int j=0;j<kParticles;j++){
            if (gParticles[j].pdg==pdgcode[i]) {
               ppid = j;
               break;
            }
            ppid = -1;
         }
         pt = sqrt( px[i]*px[i] + py[i]*py[i] );
         psy = TMath::ATan2(py[i],px[i]);
         vn = cos( (iHarm+1) * (PsiEP - psy) );
         hpVn_pT[iHarm][0]->Fill(pt, vn/res);         //harm for itself
         hpVn_Y[iHarm][0]->Fill(rap, vn/res);
         if (ppid!=-1) {
            pVn_pT[iHarm][0][ppid]->Fill(pt, vn/res);
            pVn_Y[iHarm][0][ppid]->Fill(rap, vn/res);
         }
         //v1 and v2 from EP1
         if (!iHarm) {
            for (int k=1;k<3;k++){     //harm of flow v_{k+1}
               vn = cos( (k+1) * (PsiEP-psy) );
               res=hrealre[k][1]->GetBinContent( (int)impact_b + 1 );
               hpVn_pT[k][1]->Fill(pt, vn/res);         //harm for itself
               hpVn_Y[k][1]->Fill(rap, vn/res);
               if (ppid!=-1) {
                  pVn_pT[k][1][ppid]->Fill(pt, vn/res);
                  pVn_Y[k][1][ppid]->Fill(rap, vn/res);
               }
            }
         }
      }
   }
}

#endif // #ifdef particles_cxx
