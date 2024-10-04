#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

{
  gROOT->Reset();
  const TString SaveDir="/exp/sbnd/data/users/bethanym/wire_transparency/Selection_CRT";

  gStyle->SetOptStat(0); //Removing the stats box
  gStyle->SetPalette(kCandy);
  TColor::InvertPalette(); 
  TFile  myFile("/sbnd/data/users/arellano/tpc_comm/production_231127_brazilCM/sbnd_crt_lifetime/trees/muon_hitdumper_NS_elifetime15ms_sce.root");
  TTree* myTree = (TTree*) myFile.Get("hitdumper/hitdumpertree");
  vector<double>*  ct_x1;
  vector<double>*  ct_x2;
  vector<double>*  ct_y1;
  vector<double>*  ct_y2;
  vector<double>*  ct_z1;
  vector<double>*  ct_z2;
  vector<double>*  chit_x;
  vector<double>*  chit_y;
  vector<double>*  chit_z;
  myTree->SetBranchAddress("ct_x1", &ct_x1);
  myTree->SetBranchAddress("ct_x2", &ct_x2);
  myTree->SetBranchAddress("ct_y1", &ct_y1);
  myTree->SetBranchAddress("ct_y2", &ct_y2);
  myTree->SetBranchAddress("ct_z1", &ct_z1);
  myTree->SetBranchAddress("ct_z2", &ct_z2); 
  myTree->SetBranchAddress("chit_x", &chit_x);
  myTree->SetBranchAddress("chit_y", &chit_y);
  myTree->SetBranchAddress("chit_z", &chit_z); 
  gROOT->cd(0);
  TH2F *x_y_Hit_1_South = new TH2F("xy1","xy1",70, -360, 360, 70, -360, 360);
  TH2F *x_y_Hit_1_North = new TH2F("xy2","xy2",70, -360, 360, 70, -360, 360);
  TH2F *x_y_Hit_2_South = new TH2F("xy3","xy3",70, -360, 360, 70, -360, 360);
  TH2F *x_y_Hit_2_North = new TH2F("xy4","xy4",70, -360, 360, 70, -360, 360);
  TH2F *x_z_Hit_1_Bottom = new TH2F("xz1","xz1",70, -420, 420, 70, -200, 800);
  TH2F *x_z_Hit_2_Bottom = new TH2F("xz2","xz2",70, -420, 420, 70, -200, 800);
  TH2F *x_y_AllHits_South = new TH2F("xy5","xy5",70, -360, 360, 70, -360, 360);
  TH2F *x_y_AllHits_North = new TH2F("xy6","xy6",70, -360, 360, 70, -360, 360);
  TH2F *x_z_AllHits_Bottom = new TH2F("xz3","xz3",70, -200, 800, 70, -420, 420);
  int nEntries = myTree->GetEntries(); // Get the number of entries in this tree
  for (int iEnt = 0; iEnt < nEntries; iEnt++) {
    myTree->GetEntry(iEnt); // Gets the next entry (filling the linked variables)

    for (size_t i=0; i < ct_x1->size(); ++i) {
      if(ct_z1->at(i)>-200 && ct_z1->at(i)<-150){
	x_y_Hit_1_South->Fill(ct_x1->at(i),ct_y1->at(i));
      }
      if(ct_z2->at(i)>-200 && ct_z2->at(i)<-150){   
	x_y_Hit_2_South->Fill(ct_x2->at(i),ct_y2->at(i));
      }
      if(ct_z1->at(i)>750 && ct_z1->at(i)<800){
	x_y_Hit_1_North->Fill(ct_x1->at(i),ct_y1->at(i));
      }
      if(ct_z2->at(i)>750 && ct_z2->at(i)<800){
	x_y_Hit_2_North->Fill(ct_x2->at(i),ct_y2->at(i));
      }   
      if(ct_y1->at(i)>-400 && ct_y1->at(i)<-350){
	x_z_Hit_1_Bottom->Fill(ct_x1->at(i),ct_z1->at(i));
      }
      if(ct_y2->at(i)>-400 && ct_y2->at(i)<-350){
	x_z_Hit_2_Bottom->Fill(ct_x2->at(i),ct_z2->at(i));
      }
    }
    for(size_t i=0; i < chit_x->size(); ++i){
      if(chit_z->at(i)>-200 && chit_z->at(i)<-150){
	x_y_AllHits_South->Fill(chit_x->at(i),chit_y->at(i));
      }
      if(chit_z->at(i)>750 && chit_z->at(i)<800){
        x_y_AllHits_North->Fill(chit_x->at(i),chit_y->at(i));
      }
      if(chit_y->at(i)>-400 && chit_y->at(i)<-350){
        x_z_AllHits_Bottom->Fill(chit_z->at(i),chit_x->at(i));
      }

    }

  }
  /* TCanvas*  x_y_1 = new TCanvas ("xy1", "xy1", 900, 700);

  x_y_Hit_1_South->Draw("COLZ");
  x_y_Hit_1_South->SetTitle("South CRT Muon Hit 1");
  x_y_Hit_1_South->GetXaxis()->SetTitle("x Position (cm)");
  x_y_Hit_1_South->GetYaxis()->SetTitle("y Position (cm)");

  x_y_1->SaveAs(SaveDir + "/x_y_Hit_1_South.pdf");  

  TCanvas*  x_y_2 = new TCanvas ("xy2", "xy2", 900, 700);
  x_y_Hit_1_North->Draw("COLZ");
  x_y_Hit_1_North->SetTitle("North CRT Muon Hit 1");
  x_y_Hit_1_North->GetXaxis()->SetTitle("x Position (cm)");
  x_y_Hit_1_North->GetYaxis()->SetTitle("y Position (cm)");

  x_y_2->SaveAs(SaveDir + "/x_y_Hit_1_North.pdf");

  TCanvas*  x_y_3 = new TCanvas ("xy3", "xy3", 900, 700);

  x_y_Hit_2_South->Draw("COLZ");
  x_y_Hit_2_South->SetTitle("South CRT Muon Hit 2");
  x_y_Hit_2_South->GetXaxis()->SetTitle("x Position (cm)");
  x_y_Hit_2_South->GetYaxis()->SetTitle("y Position (cm)");

  x_y_3->SaveAs(SaveDir + "/x_y_Hit_2_South.pdf");

  TCanvas*  x_y_4 = new TCanvas ("xy4", "xy4", 900, 700);
  x_y_Hit_2_North->Draw("COLZ");
  x_y_Hit_2_North->SetTitle("North CRT Muon Hit 2");
  x_y_Hit_2_North->GetXaxis()->SetTitle("x Position (cm)");
  x_y_Hit_2_North->GetYaxis()->SetTitle("y Position (cm)");

  x_y_4->SaveAs(SaveDir + "/x_y_Hit_2_North.pdf");


  TCanvas*  x_z_1 = new TCanvas ("xz1", "xz1", 900, 700);
  x_z_Hit_1_Bottom->Draw("COLZ");
  x_z_Hit_1_Bottom->SetTitle("Bottom CRT Muon Hit 1");
  x_z_Hit_1_Bottom->GetXaxis()->SetTitle("x Position (cm)");
  x_z_Hit_1_Bottom->GetYaxis()->SetTitle("z Position (cm)");

  x_z_1->SaveAs(SaveDir + "/x_z_Hit_1_Bottom.pdf");

  TCanvas*  x_z_2 = new TCanvas ("xz2", "xz2", 900, 700);

  x_z_Hit_2_Bottom->Draw("COLZ");
  x_z_Hit_2_Bottom->SetTitle("Bottom CRT Muon Hit 2");
  x_z_Hit_2_Bottom->GetXaxis()->SetTitle("x Position (cm)");
  x_z_Hit_2_Bottom->GetYaxis()->SetTitle("z Position (cm)");

  x_z_2->SaveAs(SaveDir + "/x_z_Hit_2_Bottom.pdf");

  TCanvas*  x_y_5 = new TCanvas ("xy5", "xy5", 900, 700);

  x_y_AllHits_South->Draw("COLZ");
  x_y_AllHits_South->SetTitle("South CRT Muon Hits");
  x_y_AllHits_South->GetXaxis()->SetTitle("x Position (cm)");
  x_y_AllHits_South->GetYaxis()->SetTitle("y Position (cm)");

  x_y_5->SaveAs(SaveDir + "/x_y_AllHits_South.pdf");

  TCanvas*  x_y_6 = new TCanvas ("xy6", "xy6", 900, 700);

  x_y_AllHits_North->Draw("COLZ");
  x_y_AllHits_North->SetTitle("North CRT Muon Hits");
  x_y_AllHits_North->GetXaxis()->SetTitle("x Position (cm)");
  x_y_AllHits_North->GetYaxis()->SetTitle("y Position (cm)");

  x_y_6->SaveAs(SaveDir + "/x_y_AllHits_North.pdf");
  */
  TCanvas*  x_z_3 = new TCanvas ("xz3", "xz3", 900, 700);

  x_z_AllHits_Bottom->Draw("COLZ");
  x_z_AllHits_Bottom->SetTitle("Bottom CRT Muon Hits");
  x_z_AllHits_Bottom->GetXaxis()->SetTitle("z Position (cm)");
  x_z_AllHits_Bottom->GetYaxis()->SetTitle("x Position (cm)");

  x_z_3->SaveAs(SaveDir + "/x_z_AllHits_Bottom.pdf");




}
