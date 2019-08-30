#define particles_cxx
#include "particles.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "Qvector.h"

void particles::Loop()
{

   if (fChain == 0) return;

   const Float_t minrap = 1.0;
   const Float_t maxrap = 2.0;
   Long64_t nentries = fChain->GetEntries();
   std::cout<<"I will run " << nentries <<" events"<<std::endl;
   InitHists();

   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      if (!(jentry % 100))
         std::cout << "\r I COMPLETED " << jentry << " events out of " << nentries << std::flush;
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      Long64_t aaa = GetEntry(jentry);
      if (impact_b < 0) continue;    //only final steh has b
      findQ(minrap, maxrap);
   }
   std::cout<<"\n Mean Q found\n";
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      if (!(jentry % 100))
         std::cout << "\r I COMPLETED " << jentry << " events out of " << nentries << std::flush;
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      Long64_t aaa = GetEntry(jentry);
      if (impact_b < 0) continue;    //only final steh has b
      findRes(minrap, maxrap);
   }
   findrealRes();
   std::cout<<"\n Resolution found\n";
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      if (!(jentry % 100))
         std::cout << "\r I COMPLETED " << jentry << " events out of " << nentries << std::flush;
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      Long64_t aaa = GetEntry(jentry);
      if (impact_b < 0) continue;    //only final steh has b
      findvn(minrap, maxrap);
   }
   std::cout<<"\n flows are ready\n";
   SaveHists();   
}




//for correct work do:
/* root -l
   TChain ch("particles")
   ch.Add("/mnt/pool/rhic/4/vsnikolaev/Soft/Models/Basov/smashtr/part5/*.root")
   ch.MakeClass()
   .q
*/ //and cp /mnt/pool/rhic/4/vsnikolaev/Soft/Models/Basov/smashtr/part5/*.C ../

//some helpfull commands
/*
   root -l
 TFile *_file0 = TFile::Open("collisions_88386001_1.root")
_file0->ls()
collisions->Print()
particles->GetListOfLeaves()->ls()
particles->Scan()
particles->Draw("pdgcode","TMath::Abs(pdgcode)<5000")
particles->Draw("pz","pdgcode==2212")
particles->Draw("pz","pdgcode==211")
particles->Draw("pz","pdgcode==2212")
particles->Draw("pz","pdgcode==211")
particles->Draw("pz","pdgcode==321","pz<1.0")
 */


//   In a ROOT session, you can do:
//      root> .L particles.C
//      root> particles t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch