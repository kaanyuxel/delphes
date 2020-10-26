#ifdef __CLING__
R__LOAD_LIBRARY(libDelphes)
#include "classes/DelphesClasses.h"
#include "external/ExRootAnalysis/ExRootTreeReader.h"
#endif
///////////////////////////////////////////////////////////////////////
void analysis()
{
  gSystem->Load("libDelphes");
  TFile *outputFile;
  outputFile = new TFile("pileup_100.root","RECREATE");
  // Create chain of root trees
  TChain chain("Delphes");

  chain.Add("/eos/experiment/fcc/hh/generation/DelphesStandalone/Kaan/pp2HH_output/pileUp_hh_100_1.root");
  chain.Add("/eos/experiment/fcc/hh/generation/DelphesStandalone/Kaan/pp2HH_output/pileUp_hh_100_2.root");
  chain.Add("/eos/experiment/fcc/hh/generation/DelphesStandalone/Kaan/pp2HH_output/pileUp_hh_100_3.root");
  chain.Add("/eos/experiment/fcc/hh/generation/DelphesStandalone/Kaan/pp2HH_output/pileUp_hh_100_4.root");
  chain.Add("/eos/experiment/fcc/hh/generation/DelphesStandalone/Kaan/pp2HH_output/pileUp_hh_100_5.root");
  chain.Add("/eos/experiment/fcc/hh/generation/DelphesStandalone/Kaan/pp2HH_output/pileUp_hh_100_6.root");
  chain.Add("/eos/experiment/fcc/hh/generation/DelphesStandalone/Kaan/pp2HH_output/pileUp_hh_100_7.root");
  chain.Add("/eos/experiment/fcc/hh/generation/DelphesStandalone/Kaan/pp2HH_output/pileUp_hh_100_8.root");  
  chain.Add("/eos/experiment/fcc/hh/generation/DelphesStandalone/Kaan/pp2HH_output/pileUp_hh_100_9.root");
  chain.Add("/eos/experiment/fcc/hh/generation/DelphesStandalone/Kaan/pp2HH_output/pileUp_hh_100_10.root");        
  
  // Create object of class ExRootTreeReader
  ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);
  Long64_t numberOfEntries = treeReader->GetEntries();
  
  // Get pointers to branches used in this analysis
  TClonesArray *branchJet= treeReader->UseBranch("Jet");
  TClonesArray *branchPhoton = treeReader->UseBranch("Photon");

  // Decalare TLorentzVector for Z
  TLorentzVector firstPhotonFour(0,0,0,0);
  TLorentzVector secondPhotonFour(0,0,0,0); 
  TLorentzVector jetFour(0,0,0,0);
  TLorentzVector firstHiggs(0,0,0,0);
  TLorentzVector secondHiggs(0,0,0,0);  
  TLorentzVector noTagHiggs(0,0,0,0); 


  TH1 *higgsJets    = new TH1F("higgsJets","Reco Higgs From Jets", 60, 0, 300);    
  TH1 *higgsbJets    = new TH1F("higgsbJets","Reco Higgs From bJets", 60, 0, 300);  
  TH1 *higgsPhotons = new TH1F("higgsPhotons","Reco Higgs From Photons", 50, 60, 160);     
  Jet *jet;
  Jet *bJet1, bJet2;
  Photon *p1, *p2;
  vector<Jet*> allbJets;  
  vector<Jet*> allJets; 

  double dRJandP1 = 0.0;
  double dRJandP2 = 0.0;

  // Loop over all events
  for(Int_t entry = 0; entry < numberOfEntries; ++entry)
  {
    // print event number in the file
    if(entry%1000==0) cout << "Reading Event " << entry << endl;
    
    // Load selected branches with data from specified event
    treeReader->ReadEntry(entry);
    allJets.clear();
  
    // If event contains at least 2 electrons

    if(branchJet->GetEntries() >= 2 && branchPhoton->GetEntries() >= 2)
    {
      p1 = (Photon*) branchPhoton->At(0);
      p2 = (Photon*) branchPhoton->At(1);

      if(p1->PT > 20 && p2->PT > 20)
      {
        // isolation will be written here!!!!

        for(int i = 0; i < branchJet->GetEntries(); ++i)
        {
          jet = (Jet*) branchJet->At(i);
          jetFour = jet->P4();
          dRJandP1 = (p1->P4()).DeltaR(jetFour);
          dRJandP2 = (p2->P4()).DeltaR(jetFour);

          if(dRJandP1 > 0.3 && dRJandP2 > 0.3)
            allJets.push_back((Jet*)branchJet->At(i));
        }

        if(allJets.size() >= 2)
        {  
          noTagHiggs = allJets.at(0)->P4();
          noTagHiggs += allJets.at(1)->P4();
          higgsJets->Fill(noTagHiggs.M());

          if((allJets.at(0)->BTag & (1 << 0)) && (allJets.at(1)->BTag & (1 << 0)))
          {
            if(allJets.at(0)->PT > 30 && allJets.at(1)->PT > 30 && fabs(allJets.at(0)->Eta) < 2.5 && fabs(allJets.at(1)->Eta) < 2.5)
            {
              firstHiggs = p1->P4();
              firstHiggs += p2->P4();
              higgsPhotons->Fill(firstHiggs.M());
              /////////////////////////////////////
              secondHiggs = allJets.at(0)->P4();
              secondHiggs += allJets.at(1)->P4();
              higgsbJets->Fill(secondHiggs.M());          
            }  
          }
        }  
      }
    }    
  }

  outputFile->Write();
  outputFile->Close();
 // c1->Print("pileupsub_5.png");

}
