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

//*******************************Function Definitions*******************************************//
bool Cut_NS_Function(double x1, double y1,double z1,double x2, double y2,double z2){
  bool z =(z1>-200 && z1<-150) || (z2>750 && z2<800) || (z1>750 && z1<800) || (z2>-200 && z2<-150);
  bool y = (y1>-360 && y1<360) || (y2>-360 && y2<360);
  bool x=((x1>-200 && x1<-150) && (x2>-200 && x2<-150)) || ((x1>150 && x1<200) && (x2>150 && x2<200));
  return x && y && z;
}
bool Cut_Flat_Function(double x1, double y1,double z1,double x2, double y2,double z2){
  bool z =(z1>-200 && z1<791) || (z2>-200 && z2<791);
  bool y = (y1>-380 && y1<-376) || (y2>-380 && y2<-376);
  bool x=((x1>-200 && x1<-150) && (x2>-200 && x2<-150)) || ((x1>150 && x1<200) && (x2>150 && x2<200));
  return x && y && z;
}
bool Cut_Induction_Angle_yz(double theta_yz){
  bool yz= (theta_yz>55 && theta_yz<65) || (theta_yz>-65 && theta_yz<-55);
  return yz;
}
bool Cut_Track_Length(double length){
  bool tracklength= (length>0 );
  return tracklength;
}
int Closest_Track(vector<double>* TPC_Angle, double CRT_Angle){
  int Closest_Index=-1;
  double Min_Difference=180;
  for(size_t i=0; i < TPC_Angle->size(); ++i){
    const double Difference=fabs(TPC_Angle->at(i)-CRT_Angle);
    if(Difference<Min_Difference && Difference<10){
      Min_Difference=Difference;
      Closest_Index=i;
    }
  }
  return Closest_Index;
}

//*****************************Main*********************************************************//
void Selection_induction()
{
  //******************************Accessing Tree and defining histograms*************************//
  gROOT->Reset();
  const TString SaveDir="/exp/sbnd/data/users/bethanym/wire_transparency/Selection_CRT";
  gStyle->SetOptStat(0); //Removing the stats box
  gStyle->SetPalette(kCandy);
  TColor::InvertPalette(); 
  TFile  myFile("/sbnd/data/users/arellano/tpc_comm/production_231127_brazilCM/sbnd_crt_lifetime/trees/muon_hitdumper_NS_elifetime15ms_sce.root");
  TTree* myTree = (TTree*) myFile.Get("hitdumper/hitdumpertree");
  vector<double>*  ct_x1=0;
  vector<double>*  ct_x2=0;
  vector<double>*  ct_y1=0;
  vector<double>*  ct_y2=0;
  vector<double>*  ct_z1=0;
  vector<double>*  ct_z2=0;
  vector<double>*  crt_space_point_x=0;
  vector<double>*  crt_space_point_y=0;
  vector<double>*  crt_space_point_z=0;
  vector<double>* muontrk_type=0;
  vector<double>*muontrk_theta_xz=0;
  vector<double>*muontrk_theta_yz=0;
  vector<double>*muontrk_tpc=0;
  myTree->SetBranchAddress("ct_x1", &ct_x1);
  myTree->SetBranchAddress("ct_x2", &ct_x2);
  myTree->SetBranchAddress("ct_y1", &ct_y1);
  myTree->SetBranchAddress("ct_y2", &ct_y2);
  myTree->SetBranchAddress("ct_z1", &ct_z1);
  myTree->SetBranchAddress("ct_z2", &ct_z2); 
  myTree->SetBranchAddress("crt_space_point_x", &crt_space_point_x);
  myTree->SetBranchAddress("crt_space_point_y", &crt_space_point_y);
  myTree->SetBranchAddress("crt_space_point_z", &crt_space_point_z); 
  myTree->SetBranchAddress("muontrk_type", &muontrk_type);
  myTree->SetBranchAddress("muontrk_theta_xz", &muontrk_theta_xz);
  myTree->SetBranchAddress("muontrk_theta_yz", &muontrk_theta_yz);
  myTree->SetBranchAddress("muontrk_tpc", &muontrk_tpc);
  gROOT->cd(0);
  TH2F *x_y_Hit_1 = new TH2F("xy1","xy1",70, -360, 360, 70, -360, 360);
  TH2F *x_y_Hit_2 = new TH2F("xy3","xy3",70, -360, 360, 70, -360, 360);
  TH1F *TPC_Theta_xz =new TH1F("theta_xz_1", "theta_xz_1", 50, -90, 90);
  TH1F *TPC_Theta_yz =new TH1F("theta_yz_1", "theta_yz_1", 50, -90,90);
  TH1F *CRT_Theta_xz =new TH1F("crttheta_xz_1", "crttheta_xz_1", 50, -90, 90);
  TH1F *CRT_Theta_yz =new TH1F("crttheta_yz_1", "crttheta_yz_1", 50, -90,90);
  TH2F *x_z_Hit_1 = new TH2F("xz1","xz1",70, -420, 420, 70, -200, 800);
  TH2F *x_z_Hit_2 = new TH2F("xz2","xz2",70, -420, 420, 70, -200, 800);

 
  vector<TVector3*> track_starts;
  vector<TVector3*> track_ends;
  vector<int> track_type;


  //******************************Entries loop**********************************//
  int nEntries = myTree->GetEntries(); // Get the number of entries in this tree  
  for (int iEnt = 0; iEnt < nEntries; iEnt++) { //Big loop of entries
    myTree->GetEntry(iEnt); // Gets the next entry (filling the linked variables)

    vector<int> Matched_Tracks;
    //*********************Modifying TPC Angles ****************************//
    vector<double>* Modified_TPC_theta_yz=muontrk_theta_yz;
    vector<double>* Modified_TPC_theta_xz=muontrk_theta_xz;
    for(size_t t=0; t<muontrk_type->size(); ++t){
      if ( muontrk_theta_xz->at(t) > 90 && muontrk_tpc->at(t) == 0 ) {
	Modified_TPC_theta_xz->push_back(muontrk_theta_xz->at(t) - 180.0);
      }
      else if ( muontrk_theta_xz->at(t) > 90 && muontrk_tpc->at(t) == 1 ) {
	Modified_TPC_theta_xz->push_back(180.0 - muontrk_theta_xz->at(t));
      }
      else if ( muontrk_theta_xz->at(t) < 90 && muontrk_tpc->at(t) == 1 ) {
	Modified_TPC_theta_xz->push_back(-muontrk_theta_xz->at(t));
      }
      
      if ( muontrk_theta_yz->at(t) < -90 ) {
	Modified_TPC_theta_yz->push_back( muontrk_theta_yz->at(t) + 180.);
      }
      if ( muontrk_theta_yz->at(t) > 90 ) {
	Modified_TPC_theta_yz->push_back( muontrk_theta_yz->at(t) - 180.);
      }
    }

    //***************CRT Track Loop ***************//                                                                                                                                                               
    for (size_t c=0; c < ct_x1->size(); ++c) {

      //Defining variables
      bool CutResult=false;
      double dx;
      double dy;
      double dz;
      double CRT_theta_xz;
      double CRT_theta_yz;
      double modified_theta_xz_CRT;
      double modified_theta_yz_CRT;
      double Track_Length; 
      int track_colour;
      //Calculating Variables
      dx =ct_x2->at(c) - ct_x1->at(c);
      dy =ct_y2->at(c) - ct_y1->at(c);
      dz =ct_z2->at(c) - ct_z1->at(c);
      Track_Length=sqrt((dx*dx)+(dy*dy)+(dz*dz));
      CRT_theta_xz= atan2(dx,dz)*(180/TMath::Pi());
      CRT_theta_yz= atan2(dy,dz)*(180/TMath::Pi());
      modified_theta_xz_CRT = CRT_theta_xz;
      if ( CRT_theta_xz < -90 ) {
	modified_theta_xz_CRT = CRT_theta_xz + 180.;
      }
      if ( CRT_theta_xz > 90 ) {
	modified_theta_xz_CRT = CRT_theta_xz - 180.;
      }
      modified_theta_yz_CRT = CRT_theta_yz;
      if ( CRT_theta_yz < -90 ) {
	modified_theta_yz_CRT = CRT_theta_yz + 180.;
      }
      if ( CRT_theta_yz > 90 ) {
	modified_theta_yz_CRT = CRT_theta_yz - 180.;
      }
      //Applying Cuts
      if(Cut_NS_Function(ct_x1->at(c),ct_y1->at(c),ct_z1->at(c),ct_x2->at(c),ct_y2->at(c),ct_z2->at(c)) && Cut_Flat_Function(ct_x1->at(c),ct_y1->at(c),ct_z1->at(c),ct_x2->at(c),ct_y2->at(c),ct_z2->at(c))&& Cut_Induction_Angle_yz(modified_theta_yz_CRT) && Cut_Track_Length(Track_Length)){
	CutResult=true;
  
	x_y_Hit_1->Fill(ct_x1->at(c),ct_y1->at(c));
	x_y_Hit_2->Fill(ct_x2->at(c),ct_y2->at(c));
	x_z_Hit_1->Fill(ct_x1->at(c),ct_z1->at(c));
	x_z_Hit_2->Fill(ct_x2->at(c),ct_z2->at(c));

	TVector3* track_start= new TVector3(ct_x1->at(c),ct_y1->at(c),ct_z1->at(c));
	TVector3* track_end= new TVector3(ct_x2->at(c),ct_y2->at(c),ct_z2->at(c));

	track_starts.push_back(track_start);
	track_ends.push_back(track_end);

	CRT_Theta_xz->Fill(modified_theta_xz_CRT);
	CRT_Theta_yz->Fill(modified_theta_yz_CRT);


	std::string index = Form("event%i_track%zu_angle%f",iEnt,c,modified_theta_yz_CRT);
	event_display(track_start,track_end,index);

	if(modified_theta_yz_CRT>0){                                                                                                                                                                       
	  track_type.push_back(kBlue-6);
	}                                                                                                                                                                                                    
	else if(modified_theta_yz_CRT<0){                                                                                                                                                                     
	  track_type.push_back(kTeal-5);                                                                                                                                                                   
	}                                                                                                                                                                                                     
     
      
	//******************TPC Track Matching *********************************//

	int Closest_Index= Closest_Track(Modified_TPC_theta_yz, modified_theta_yz_CRT);
	if(std::find(Matched_Tracks.begin(),Matched_Tracks.end(),Closest_Index)== Matched_Tracks.end()&& !(Closest_Index==-1)){
	  Matched_Tracks.push_back(Closest_Index);
	  Double_t theta_xz = Modified_TPC_theta_xz->at(Closest_Index);
	  Double_t theta_yz = Modified_TPC_theta_yz->at(Closest_Index);
	  TPC_Theta_xz->Fill(theta_xz);
	  TPC_Theta_yz->Fill(theta_yz);     

	}
      }// cut if loop
    }// end of CRT track loop                                                                                                                                                                                        
    
  } //Big loop of entries
  event_display_multiple(track_starts,track_ends,track_type);





  //*************** Making pretty histograms ********************//
  
  TCanvas*  x_y_1 = new TCanvas ("xy1", "xy1", 900, 700);
 
  x_y_Hit_1->Draw("COLZ");
  x_y_Hit_1->SetTitle("CRT Muon Hit 1");
  x_y_Hit_1->GetXaxis()->SetTitle("x Position (cm)");
  x_y_Hit_1->GetYaxis()->SetTitle("y Position (cm)");

  x_y_1->SaveAs(SaveDir + "/x_y_Hit_1.pdf");  


  TCanvas*  x_y_3 = new TCanvas ("xy3", "xy3", 900, 700);

  x_y_Hit_2->Draw("COLZ");
  x_y_Hit_2->SetTitle("CRT Muon Hit 2");
  x_y_Hit_2->GetXaxis()->SetTitle("x Position (cm)");
  x_y_Hit_2->GetYaxis()->SetTitle("y Position (cm)");

  x_y_3->SaveAs(SaveDir + "/x_y_Hit_2.pdf");
 
  TCanvas*  x_z_1 = new TCanvas ("xz1", "xz1", 900, 700);

  x_z_Hit_1->Draw("COLZ");
  x_z_Hit_1->SetTitle("CRT Muon Hit 1");
  x_z_Hit_1->GetXaxis()->SetTitle("x Position (cm)");
  x_z_Hit_1->GetYaxis()->SetTitle("z Position (cm)");

  x_z_1->SaveAs(SaveDir + "/x_z_Hit_1.pdf");


  TCanvas*  x_z_2 = new TCanvas ("xz3", "xz3", 900, 700);

  x_z_Hit_2->Draw("COLZ");
  x_z_Hit_2->SetTitle("CRT Muon Hit 2");
  x_z_Hit_2->GetXaxis()->SetTitle("x Position (cm)");
  x_z_Hit_2->GetYaxis()->SetTitle("z Position (cm)");

  x_z_2->SaveAs(SaveDir + "/x_z_Hit_2.pdf");

  
  
  TCanvas*  TPC_theta_xz = new TCanvas ("theta_xz1", "theta_xz1", 900, 700);

  TPC_Theta_xz->Draw("");
  TPC_Theta_xz->SetTitle("#theta_{xz} TPC");
  TPC_Theta_xz->GetXaxis()->SetTitle("Angle (degrees)");
  TPC_Theta_xz->GetYaxis()->SetTitle("");
  TPC_Theta_xz->SetLineWidth(4);
  TPC_Theta_xz->SetLineColor(kRed-6);

  TPC_theta_xz->SaveAs(SaveDir + "/TPC_Theta_xz.pdf");
  
  TCanvas*  TPC_theta_yz = new TCanvas ("theta_xy1", "theta_xy1", 900, 700);

  TPC_Theta_yz->Draw("");
  TPC_Theta_yz->SetTitle("#theta_{yz} TPC");
  TPC_Theta_yz->GetXaxis()->SetTitle("Angle (degrees)");
  TPC_Theta_yz->GetYaxis()->SetTitle("");
  TPC_Theta_yz->SetLineWidth(4);
  TPC_Theta_yz->SetLineColor(kRed-6);

  TPC_theta_yz->SaveAs(SaveDir + "/TPC_Theta_yz.pdf");
  
  
  TCanvas*  CRT_theta_xz = new TCanvas ("crttheta_xz1", "crttheta_xz1", 900, 700);

  CRT_Theta_xz->Draw("");
  CRT_Theta_xz->SetTitle("#theta_{xz} CRT");
  CRT_Theta_xz->GetXaxis()->SetTitle("Angle (degrees)");
  CRT_Theta_xz->GetYaxis()->SetTitle("");
  CRT_Theta_xz->SetLineWidth(4);
  CRT_Theta_xz->SetLineColor(kRed-6);

  CRT_theta_xz->SaveAs(SaveDir + "/CRT_Theta_xz.pdf");


  TCanvas*  CRT_theta_yz = new TCanvas ("crttheta_xy1", "crttheta_xy1", 900, 700);

  CRT_Theta_yz->Draw("");
  CRT_Theta_yz->SetTitle("#theta_{yz} CRT");
  CRT_Theta_yz->GetXaxis()->SetTitle("Angle (degrees)");
  CRT_Theta_yz->GetYaxis()->SetTitle("");
  CRT_Theta_yz->SetLineWidth(4);
  CRT_Theta_yz->SetLineColor(kRed-6);

  CRT_theta_yz->SaveAs(SaveDir + "/CRT_Theta_yz.pdf");

}
