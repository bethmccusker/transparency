#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include<TMath.h>
#include<TTree.h>
#include<TFile.h>
#include<TROOT.h>
#include<iostream>
#include"event_display_test.C"
#include<sstream>

void chit_eventdisplaytest()
{
  gROOT->Reset();
  const TString SaveDir="/exp/sbnd/data/users/bethanym/wire_transparency/Selection_CRT";
  //  gStyle->SetOptStat(0); //Removing the stats box
  TColor::InvertPalette(); 
  TFile  myFile("/exp/sbnd/data/users/bethanym/wire_transparency/hd_variable_test/hists_decode_data_evb01_EventBuilder1_art1_run16740_10_20240912T082517-30bef869-9d29-42a0-aa77-091ad9c1620d_Reco1Comm-20241029T161501.root");
  TTree* myTree = (TTree*) myFile.Get("hitdumper/hitdumpertree");
  vector<double>*  ct_x1=0;
  vector<double>*  ct_x2=0;
  vector<double>*  ct_y1=0;
  vector<double>*  ct_y2=0;
  vector<double>*  ct_z1=0;
  vector<double>*  ct_z2=0;
  int track_colour;
  vector<int> track_type;
  myTree->SetBranchAddress("crt_track_x1", &ct_x1);
  myTree->SetBranchAddress("crt_track_x2", &ct_x2);
  myTree->SetBranchAddress("crt_track_y1", &ct_y1);
  myTree->SetBranchAddress("crt_track_y2", &ct_y2);
  myTree->SetBranchAddress("crt_track_z1", &ct_z1);
  myTree->SetBranchAddress("crt_track_z2", &ct_z2); 
  gROOT->cd(0);
  int nEntries = myTree->GetEntries(); // Get the number of entries in this tree
  vector<TVector3*> track_starts;
  vector<TVector3*> track_ends;
  for (int iEnt = 0; iEnt < nEntries; iEnt++) { //Big loop of entries
    myTree->GetEntry(iEnt); // Gets the next entry (filling the linked variables)
   
    int n = ct_x1->size();
    for (size_t c=0; c < n; ++c) {
      double dx;
      double dy;
      double dz;
      double CRT_theta_xz;
      double CRT_theta_yz;
      
      dx =ct_x2->at(c) - ct_x1->at(c);
      dy =ct_y2->at(c) - ct_y1->at(c);
      dz =ct_z2->at(c) - ct_z1->at(c);

      TVector3* track_start= new TVector3(ct_x1->at(c),ct_y1->at(c),ct_z1->at(c));
      TVector3* track_end= new TVector3(ct_x2->at(c),ct_y2->at(c),ct_z2->at(c));    
       
      track_starts.push_back(track_start);
      track_ends.push_back(track_end); 

      CRT_theta_xz= atan2(dx,dz)*(180/TMath::Pi());
      CRT_theta_yz= atan2(dy,dz)*(180/TMath::Pi());
        
      double modified_theta_xz_CRT = CRT_theta_xz;
      if ( CRT_theta_xz < -90 ) {
	modified_theta_xz_CRT = CRT_theta_xz + 180.;
      }
      if ( CRT_theta_xz > 90 ) {
	modified_theta_xz_CRT = CRT_theta_xz - 180.;
      }
      double modified_theta_yz_CRT = CRT_theta_yz;
      if ( CRT_theta_yz < -90 ) {
	modified_theta_yz_CRT = CRT_theta_yz + 180.;
      }
      if ( CRT_theta_yz > 90 ) {
	modified_theta_yz_CRT = CRT_theta_yz - 180.;
      }
     
      std::string index = Form("event%i_track%zu_angle%f",iEnt,c,modified_theta_yz_CRT);    
      event_display(track_start,track_end,index);

      if(modified_theta_yz_CRT>0){
	track_type.push_back(kBlue-6);
      }
      else if(modified_theta_yz_CRT<0){
	track_type.push_back(kTeal-5);
      }
    }  
  } //Big loop of entries
  event_display_multiple(track_starts,track_ends,track_type);
}
