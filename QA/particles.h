//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Aug 12 18:32:20 2019 by ROOT version 6.14/02
// from TChain particles/
//////////////////////////////////////////////////////////

#ifndef particles_h
#define particles_h

#include <TROOT.h>
#include <iostream>
#include <fstream>
#include <math.h>
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
   TH1D *hB;
   TH2D *hppt;
   TH2D *hpEn;
   TH1D *hp[kParticles];
   TH1D *hpx[kParticles];
   TH1D *hpy[kParticles];
   TH1D *hpz[kParticles];
   TH1D *hpt[kParticles];
   TH1D *hrap[kParticles];
   TH1D *hpcr[kParticles];
   TH2D *hptrap[kParticles];
   TH2D *hptpsr[kParticles];
   TH2D *hrappsr[kParticles];
   // Flow
   TProfile *pVn_pT[2][kParticles];
   TProfile *pVn_Y[2][kParticles];
   TProfile *hpVn_pT[2];
   TProfile *hpVn_Y[2];


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
	   TChain * chain = new TChain("particles", "");
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
  // HistsDell();
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
   hB = new TH1D("hB","Impact parameter; B (fm);nEvents", 200, 0, 20);
   hpVn_Y[0] = new TProfile("v1_from_Y", "v1 from Y; Y; v1", 100, -1, 1);
   hpVn_Y[1] = new TProfile("v2_from_Y", "v2 from Y; Y; v1", 100, -1, 1);
   hpVn_pT[0] = new TProfile("v1_from_p_{t}", "v1 from p_{t}; p_{t}; v1", 100, 0, 5);
   hpVn_pT[1] = new TProfile("v2_from_p_{t}", "v1 from p_{t}; p_{t}; v1", 100, 0, 5);
   for (int i=0; i<kParticles; i++){
      hp[i] = new TH1D("hp"+gParticles[i].name, "P "+gParticles[i].displayName+"; p, [GeV/c]; nEntries",100,0,5);
      hpx[i] = new TH1D("hpx"+gParticles[i].name, "Px "+gParticles[i].displayName+"; px, [GeV/c]; nEntries",100, -1.5, 1.5);
      hpy[i] = new TH1D("hpy"+gParticles[i].name, "Py "+gParticles[i].displayName+"; py, [GeV/c]; nEntries",100, -1.5, 1.5);
      hpz[i] = new TH1D("hpz"+gParticles[i].name, "Pz "+gParticles[i].displayName+"; pz, [GeV/c]; nEntries",100, -5, 5);
      hpt[i] = new TH1D("hpt"+gParticles[i].name, "pt "+gParticles[i].displayName+"; pt, [GeV/c]; nEntries",100, 0, 2.5);
      hrap[i] = new TH1D("hrap"+gParticles[i].name, "rapidity "+gParticles[i].displayName+"; rapidity; nEntries",100, -2.5, 2.5);
      hpcr[i] = new TH1D("hpcr"+gParticles[i].name, "Pseudorapidity "+gParticles[i].displayName+"; pseudorapidity; nEntries",100, -8, 8);
      hptrap[i] = new TH2D("hptrap"+gParticles[i].name, "pt vs rapidity "+gParticles[i].displayName+"; pt, [GeV/c]; rapidity",100, 0, 2.5, 100, -2.5, 2.5);
      hptpsr[i] = new TH2D("hptpsr"+gParticles[i].name, "pt vs Pseudorapidity "+gParticles[i].displayName+"; pt, [GeV/c]; pseudorapidity",100, 0, 2.5, 100, -8, 8);
      hrappsr[i] = new TH2D("hrappsr"+gParticles[i].name, "rapidity vs Pseudorapidity "+gParticles[i].displayName+"; rapidity; pseudorapidity",100, -2.5, 2.5, 100, -8, 8);

         pVn_Y[0][i]  = new TProfile("v1_from_Y_"+gParticles[i].name, "v1 from Y "+gParticles[i].displayName+"; Y; v1", 100, -1, 1);
         pVn_Y[1][i]  = new TProfile("v2_from_Y_"+gParticles[i].name, "v2 from Y "+gParticles[i].displayName+"; Y; v2", 100, -1, 1);
         pVn_pT[0][i] = new TProfile("v1_from_p_{t}_"+gParticles[i].name, "v1 from p_{t} "+gParticles[i].displayName+"; p_{t}, [GeV/c]; v1", 100, 0, 2);
         pVn_pT[1][i] = new TProfile("v2_from_p_{t}_"+gParticles[i].name, "v2 from p_{t} "+gParticles[i].displayName+"; p_{t}, [GeV/c]; v2", 100, 0, 2);
      
   }
   hppt = new TH2D("hppt", "p_{t} vs p; p_{t}, [GeV/c]; p, [GeV/c]", 100, 0, 2.5, 100, 0, 5);
   hpEn = new TH2D("hpEn", "p_{t} vs Energy; p_{t}, [GeV/c]; E, [GeV]", 100, 0, 2.5, 100, 0, 5);
}
void particles::FillHists(){
   hB->Fill(impact_b);
   for (Long64_t i = 0; i < npart; i++){
      Int_t ppid = 10;
      for (int j=0;j<kParticles;j++){
         if (gParticles[j].pdg==pdgcode[i]) {
            ppid = j;
            break;
         }
         ppid = 10;
      }
      if (ppid==10) continue;
      Double_t rap = 1.0/2.0 * log( (p0[i] + pz[i]) / (p0[i] - pz[i]) );
      Double_t momentum = sqrt( px[i]*px[i] + py[i]*py[i] + pz[i]*pz[i] );
      Double_t psrap = 1.0/2.0 * log( (momentum + pz[i]) / (momentum - pz[i]) );
      Double_t pt = sqrt( px[i]*px[i] + py[i]*py[i] );
      if (ppid==kPROTON || ppid==kNEUTRON){      //kill the main parn of spectators.   to kill all |rap|<1, but...
         if (fabs(pz[i]) > 1.8 && fabs(pz[i]) < 3){  
            if (pt<0.5){
              // if (momentum>1.6 && momentum < 3.2){
                  continue;
             //  }
            }
         }
      }
      //std::cout<<pt<<std::endl;
      hpx[ppid]->Fill(px[i]);
      hpy[ppid]->Fill(py[i]);
      hpz[ppid]->Fill(pz[i]);
      hp[ppid]->Fill(momentum);
      hpt[ppid]->Fill(pt);
      hrap[ppid]->Fill(rap);
      hpcr[ppid]->Fill(psrap);
      hptrap[ppid]->Fill(pt,rap);
      hptpsr[ppid]->Fill(pt, psrap);
      hrappsr[ppid]->Fill(rap, psrap);
      hppt->Fill(pt, momentum);
      hpEn->Fill(pt, p0[i]);
      Double_t flow;
      Int_t flowSign=0;
      Double_t phi, psiRP;
      phi=TMath::ATan2(py[i],px[i]);
      psiRP=0;


      for (Int_t iHarm=0; iHarm<2; ++iHarm) {
         if (iHarm % 2 == 0 && rap < 0.) flowSign = -1;
         else flowSign = 1;
         flow = cos ( (iHarm + 1) * ( phi - psiRP) );
         hpVn_pT[iHarm]->Fill(pt, flowSign * flow );
         hpVn_Y[iHarm]->Fill(rap, flow );
    }
      for (Int_t iHarm=0; iHarm<2; iHarm++){   //1 and 2
         if (iHarm % 2 == 0 && rap < 0.) flowSign = -1;
         else flowSign = 1;
         flow = cos ( (iHarm + 1) * ( phi - psiRP) );
         pVn_pT[iHarm][ppid]->Fill(pt, flowSign * flow  );
          pVn_Y[iHarm][ppid]->Fill(rap, flow );
      }
   }
}
void particles::SaveHists(){
   //root save
   TFile* w = new TFile("QA.root", "recreate");
   TDirectory *outputDir;
   w->cd();
   hB->Write();
   hppt->Write();
   hpEn->Write();
   hpVn_Y[0]->Write();
   hpVn_Y[1]->Write();
   hpVn_pT[0]->Write();
   hpVn_pT[1]->Write();
   for (int i=0;i<kParticles;i++) {
    outputDir = w->mkdir(gParticles[i].name);
    outputDir->cd ();
    hpx[i]->Write();
    hpy[i]->Write();
    hpz[i]->Write();
    hp[i]->Write();
    hpt[i]->Write();
    hrap[i]->Write();
    hpcr[i]->Write();
    hptrap[i]->Write();
    hptpsr[i]->Write();
    hrappsr[i]->Write();
       pVn_Y[0][i]->Write();
       pVn_Y[1][i]->Write();
      pVn_pT[0][i]->Write();
      pVn_pT[1][i]->Write();
   }
   w->Close();
   delete w;
//  HistsDell();
}

void particles::HistsDell(){
   delete hB;
   delete hppt;
   delete hpEn;
   delete hpVn_Y[0];
   delete hpVn_Y[1];
   delete hpVn_pT[0];
   delete hpVn_pT[1];
   for (int i=0;i<kParticles;i++){
      delete hp[i];
      delete hpx[i];
      delete hpy[i];
      delete hpz[i];
      delete hpt[i];
      delete hrap[i];
      delete hpcr[i];
      delete hptrap[i];
      delete hptpsr[i];
      delete hrappsr[i];
      delete pVn_Y[0][i];
      delete pVn_Y[1][i];
      delete pVn_pT[0][i];
      delete pVn_pT[1][i];
   }
}

#endif // #ifdef particles_cxx
