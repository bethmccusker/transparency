#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include<TMath.h>
#include<TTree.h>
#include<TFile.h>
#include<TROOT.h>
#include<iostream>


void DrawCube(TCanvas *c1, double *rmin, double *rmax, int colour)
{
  c1->cd();
  TList *outline = new TList;
  TPolyLine3D *p1 = new TPolyLine3D(4);
  TPolyLine3D *p2 = new TPolyLine3D(4);
  TPolyLine3D *p3 = new TPolyLine3D(4);
  TPolyLine3D *p4 = new TPolyLine3D(4);
  p1->SetLineColor(colour);
  p1->SetLineWidth(2);
  p1->Copy(*p2);
  p1->Copy(*p3);
  p1->Copy(*p4);
  outline->Add(p1);
  outline->Add(p2);
  outline->Add(p3);
  outline->Add(p4); 
  TPolyLine3D::DrawOutlineCube(outline, rmin, rmax);
  p1->Draw();
  p2->Draw();
  p3->Draw();
  p4->Draw();

}

void DrawTrack(TCanvas *c1, TVector3* point1,TVector3* point2, int colour)
{
  TPolyLine3D *track = new TPolyLine3D(2);
  track->SetLineWidth(2);
  track->SetLineColor(colour);
  track->SetPoint(0,point1->X(), point1->Y(), point1->Z());
  track->SetPoint(1,point2->X(),point2->Y(),point2->Z());
  track->Draw();
}
void event_display(TVector3 *track_start,TVector3 *track_end,std::string index){

  TCanvas *c1 = new TCanvas("c1","",1200,1200);

  double  rmin[3] {-200,-200,0};
  double  rmax[3] {-200,200,500};
  double CRT_south_min[3] {-200,360,-179.05};
  double CRT_south_max[3] {-200,-360,-179.05};
  double CRT_north_min[3]{-200,360,770.95};
  double CRT_north_max[3]{-200,-360,770.95};
  double CRT_flat_min[3]{-200,-378,-199.05};
  double CRT_flat_max[3]{-200,-378,790.95};
  double CRT_lower_top_min[3]{-200,640.12,-154.05};
  double CRT_lower_top_max[3]{-200,640.12,745.95};
  double CRT_higher_top_min[3]{-200,928.6,-154.05};
  double CRT_higher_top_max[3]{-200,928.6,745.95};

  DrawCube(c1, rmin, rmax, kRed-6 );
  DrawCube(c1, CRT_north_min, CRT_north_max, kOrange+8 );
  DrawCube(c1, CRT_south_min, CRT_south_max, kOrange+8 );
  DrawCube(c1, CRT_flat_min, CRT_flat_max, kOrange+8 );
  DrawCube(c1, CRT_lower_top_min, CRT_lower_top_max, kOrange+8 );
  DrawCube(c1, CRT_higher_top_min, CRT_higher_top_max, kOrange+8 );

  DrawTrack(c1,track_start,track_end,kBlue-6);
  

  double arrow_start_x = - 0.7;
  double arrow_start_y = -0.8;
  double arrow_length = 0.2;
  TLatex coord_y;
  TLatex coord_z;
  coord_y.DrawLatex(arrow_start_x-0.015,arrow_start_y+0.2+0.03, "y");
  coord_z.DrawLatex(arrow_start_x+0.13 +0.02, arrow_start_y-0.02, "z");
  TArrow *arrow_y = new TArrow(arrow_start_x, arrow_start_y, arrow_start_x, arrow_start_y+0.2, 0.01, "|>");
  TArrow *arrow_z = new TArrow(arrow_start_x, arrow_start_y, arrow_start_x+0.13, arrow_start_y, 0.01, "|>");
  arrow_y->SetLineWidth(2);
  arrow_z->SetLineWidth(2);
  arrow_y->Draw();
  arrow_z->Draw();


  TView3D *view = (TView3D*) TView::CreateView(1);

  double c[3] = { 0, 0, 250 };
  double s[3] = { 1000, -1000, 1000 };
  view->SetRange(-600, -600, -600, 600, 600, 600);
  view->DefineViewDirection(s, c,
                            1, 0,
                            0, 1,
                            0, 1,
                            view->GetTnorm(),
                            view->GetTback());

  const std::string SaveDir="/exp/sbnd/data/users/bethanym/wire_transparency/Selection_CRT";
  std::string event = SaveDir + ("/event_display_" + index + ".pdf");
  c1->SaveAs(event.c_str());
  //   c1->SaveAs(SaveDir + "/event_display.pdf");                                                                                                                                                    

  //  delete c1;                                                                                                                                                                                    
}
void DrawMultipleTracks(TCanvas *c1, vector<TVector3*>& point_1, vector<TVector3*>& point_2 ,vector<int> colour) {
  size_t num_tracks=point_1.size();
  for (size_t i = 0; i < num_tracks; ++i) {
    DrawTrack(c1,point_1.at(i),point_2.at(i), colour.at(i));
  }
}


void event_display_multiple(vector<TVector3*>& track_start,vector<TVector3*>& track_end, vector<int> track_colour) {

  TCanvas *c1 = new TCanvas("c1","");
  double  rmin[3] {-200,-200,0};
  double  rmax[3] {-200,200,500};
  double CRT_south_min[3] {-200,360,-179.05};
  double CRT_south_max[3] {-200,-360,-179.05};
  double CRT_north_min[3]{-200,360,770.95};
  double CRT_north_max[3]{-200,-360,770.95};
  double CRT_flat_min[3]{-200,-378,-199.05};
  double CRT_flat_max[3]{-200,-378,790.95};
  double CRT_lower_top_min[3]{-200,640.12,-154.05};
  double CRT_lower_top_max[3]{-200,640.12,745.95};
  double CRT_higher_top_min[3]{-200,928.6,-154.05};
  double CRT_higher_top_max[3]{-200,928.6,745.95};
  DrawCube(c1, rmin, rmax, kRed-6 );
  DrawCube(c1, CRT_north_min, CRT_north_max, kOrange+6 );
  DrawCube(c1, CRT_south_min, CRT_south_max, kOrange+6 );
  DrawCube(c1, CRT_flat_min, CRT_flat_max, kOrange+6 );
  DrawCube(c1, CRT_lower_top_min, CRT_lower_top_max, kOrange+6 );
  DrawCube(c1, CRT_higher_top_min, CRT_higher_top_max, kOrange+6 );
 
  DrawMultipleTracks(c1, track_start, track_end, track_colour);


  double arrow_start_x = - 0.7;
  double arrow_start_y = -0.8;
  double arrow_length = 0.2;
  TLatex coord_y;
  TLatex coord_z;
  coord_y.DrawLatex(arrow_start_x-0.015,arrow_start_y+0.2+0.03, "y");
  coord_z.DrawLatex(arrow_start_x+0.13 +0.02, arrow_start_y-0.02, "z");
  TArrow *arrow_y = new TArrow(arrow_start_x, arrow_start_y, arrow_start_x, arrow_start_y+0.2, 0.01, "|>");
  TArrow *arrow_z = new TArrow(arrow_start_x, arrow_start_y, arrow_start_x+0.13, arrow_start_y, 0.01, "|>");
  arrow_y->SetLineWidth(2);
  arrow_z->SetLineWidth(2);
  arrow_y->Draw();
  arrow_z->Draw();

  TView3D *view = (TView3D*) TView::CreateView(1);
  double c[3] = { 0, 0, 250 };
  double s[3] = { 800, -800, 800 };
  view->SetRange(-600, -600, -300, 600, 600, 900);
  view->DefineViewDirection(s, c,
                            1, 0,
                            0, 1,
                            0, 1,
                            view->GetTnorm(),
                            view->GetTback());

  const TString SaveDir="/exp/sbnd/data/users/bethanym/wire_transparency/Selection_CRT";
  c1->SaveAs(SaveDir + "/event_display_multiple.pdf");                                                                     
  //  delete c1;                                                                                                     
}




void event_display_test(){

}
