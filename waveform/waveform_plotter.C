#include "TFile.h"
#include "TTree.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TH1F.h"
#include <iostream>
#include <map>
#include <vector>


std::pair<double, double> CalculateHalfWidthHeightAndAmplitudeCollection(const std::vector<double>& adc_vals, const std::vector<double>& time_vals) {
  if (adc_vals.empty() || time_vals.empty()) return {-1.0, -1.0};
    
  double max_adc = *std::max_element(adc_vals.begin(), adc_vals.end());
  double half_max = max_adc / 2.0;
    
  double t_left = -1, t_right = -1;
  for (size_t i = 1; i < adc_vals.size(); ++i) {
    if (adc_vals[i-1] < half_max && adc_vals[i] >= half_max) {
      t_left = time_vals[i-1] + (time_vals[i] - time_vals[i-1]) * 
	(half_max - adc_vals[i-1]) / (adc_vals[i] - adc_vals[i-1]);
    }
    if (adc_vals[i-1] >= half_max && adc_vals[i] < half_max) {
      t_right = time_vals[i-1] + (time_vals[i] - time_vals[i-1]) * 
	(half_max - adc_vals[i-1]) / (adc_vals[i] - adc_vals[i-1]);
      break;
    }
  }
  double half_width = (t_left >= 0 && t_right >= 0) ? (t_right - t_left) : -1.0;
  return {half_width, max_adc};
}


std::pair<double, double> CalculateHalfWidthHeightAndAmplitudeInduction(const std::vector<double>& adc_vals, const std::vector<double>& time_vals) {
  if (adc_vals.empty() || time_vals.empty()) return {-1.0, -1.0};
    
  double min_adc = *std::min_element(adc_vals.begin(), adc_vals.end());
  double half_min = min_adc / 2.0;
    
  double t_left = -1, t_right = -1;
  for (size_t i = 1; i < adc_vals.size(); ++i) {
    if (adc_vals[i-1] > half_min && adc_vals[i] <= half_min) {
      t_left = time_vals[i-1] + (time_vals[i] - time_vals[i-1]) * 
	(half_min - adc_vals[i-1]) / (adc_vals[i] - adc_vals[i-1]);
    }
    if (adc_vals[i-1] <= half_min && adc_vals[i] > half_min) {
      t_right = time_vals[i-1] + (time_vals[i] - time_vals[i-1]) * 
	(half_min - adc_vals[i-1]) / (adc_vals[i] - adc_vals[i-1]);
      break;
    }
  }
  double half_width = (t_left >= 0 && t_right >= 0) ? (t_right - t_left) : -1.0;
  return {half_width, min_adc};
}

void waveform_plotter() {
  gStyle->SetOptStat(0); //Removing the stats box
  gStyle->SetPalette(kCandy);
  TColor::InvertPalette();
  // Open the ROOT file
  TFile *file = TFile::Open("/exp/sbnd/data/users/bethanym/wire_transparency/filter_test/hists_decode_data_evb01_EventBuilder1_art1_run16740_10_20240912T082517-30bef869-9d29-42a0-aa77-091ad9c1620d_Reco1Comm-20250214T115006.root");

  // Access the tree containing your waveform data
  TTree *tree = (TTree*)file->Get("hitdumper/hitdumpertree");

  // Variables to store the data
  std::vector<int> *waveform_number = nullptr;
  std::vector<float> *adc_on_wire = nullptr;
  std::vector<float> *time_for_waveform = nullptr;
  std::vector<float> *wire_number = nullptr;

  // Set branch addresses to connect the variables to the tree
  tree->SetBranchAddress("waveform_number", &waveform_number);
  tree->SetBranchAddress("adc_on_wire", &adc_on_wire);
  tree->SetBranchAddress("time_for_waveform", &time_for_waveform);
  tree->SetBranchAddress("wire_number", &wire_number);

  const char *output_dir = "/exp/sbnd/data/users/bethanym/wire_transparency/templateFitting/";

  TH2D* heat_map = new TH2D("h_wire_vs_peak_time", "",  1700, 0, 1700, 3500,0, 3500);
  TH2D* width_map = new TH2D("h_wire_vs_width", "",  1700, 0, 1700, 500,0, 10);
  TH2D* peak_map = new TH2D("h_wire_vs_peak_amplitude", "",  1700, 0, 1700, 500,0, 300);
  TH2D* peak_vs_width = new TH2D("h_peak_vs_width", "",  500, 0,300, 500,0, 10);

  // Map to store waveform data indexed by waveform number
  std::map<int, std::vector<double>> waveform_adc_map;
  std::map<int, std::vector<double>> waveform_time_map;
  std::map<int, std::vector<int>> waveform_wire_map;

  // Loop over tree entries (events)
  Long64_t nEntries = tree->GetEntries();
  for (Long64_t i = 0; i < nEntries; i++) {
    tree->GetEntry(i); // Load the ith entry

    // Loop through each data point in the current entry and organize by waveform number
    for (size_t j = 0; j < waveform_number->size(); j++) {
      // Check for valid ADC and time values (filter out -9999 default)
      if (adc_on_wire->at(j) != -9999 && time_for_waveform->at(j) != -9999) {
	// Insert ADC and time values into the corresponding waveform number
	int wave_num = waveform_number->at(j);
	waveform_adc_map[wave_num].push_back(adc_on_wire->at(j));
	waveform_time_map[wave_num].push_back(time_for_waveform->at(j));
	waveform_wire_map[wave_num].push_back(wire_number->at(j));
      }
    }
  }

  // Now plot each waveform based on the unique waveform numbers
  for (const auto &entry : waveform_adc_map) {
    int wave_num = entry.first; // The current waveform number
    const std::vector<double> &adc_vals = entry.second;
    const std::vector<double> &time_vals = waveform_time_map[wave_num];
    const std::vector<int> &wire_num = waveform_wire_map[wave_num];
    // If no data, skip this waveform
    if (adc_vals.empty() || time_vals.empty()) {
      continue;
    }

    auto [half_width, amplitude] = CalculateHalfWidthHeightAndAmplitudeCollection(adc_vals, time_vals);
    //    std::cout << "Waveform " << wave_num << " Half-Width Height: " << half_width << " Amplitude: " << amplitude << "\n";
    // std::cout << " Wire " <<wire_num[0]<< endl; 
    // std::cout << " Time " <<time_vals[0]<< endl;
    // std::cout << " Wave " <<wave_num<< endl;

    //    auto [trough_half_width, trough_amplitude] = CalculateHalfWidthHeightAndAmplitudeInduction(adc_vals, time_vals);
    //    std::cout << " | Trough Half-Width: " << trough_half_width << " Trough Amplitude: " << trough_amplitude << "\n";

    heat_map->Fill(wire_num[0],time_vals[0]);
    width_map->Fill(wire_num[0],half_width);
    peak_map->Fill(wire_num[0],amplitude);
    peak_vs_width->Fill(half_width,amplitude);

    // Create a TGraph for the current waveform
    TGraph *graph = new TGraph(adc_vals.size(), &time_vals[0], &adc_vals[0]);
    graph->SetTitle(Form("Waveform Number %d;Time (ticks);ADC counts", wave_num));
    graph->SetLineColor(kRed-6);
    graph->SetLineWidth(3);

    // Create a canvas to draw the graph
    TCanvas *c1 = new TCanvas(Form("c1_waveform_%d", wave_num), "Waveform", 800, 600);
    graph->Draw("AL");

    // Save the plot as a PNG file
    //    c1->SaveAs(Form("%swaveform_%d.png",output_dir, wave_num));

    // Optionally, delete canvas and graph to free memory
    delete c1;
    delete graph;
  }

  TCanvas* heat = new TCanvas("heat", "2D Histogram Canvas", 800, 600);
  heat_map->GetXaxis()->SetTitle("Wire Number");
  heat_map->GetYaxis()->SetTitle("Time (ticks)");
  heat_map->Draw("COLZ");
  heat->SaveAs(Form("%sheat_map.pdf", output_dir));

  TCanvas* width = new TCanvas("width", "2D Histogram Canvas", 800, 600);
  width_map->GetXaxis()->SetTitle("Wire Number");
  width_map->GetYaxis()->SetTitle("Width at Half Height (ticks)");
  width_map->Draw("COLZ");
  width->SaveAs(Form("%swidth_map.pdf", output_dir));

  TCanvas* peak = new TCanvas("peak", "2D Histogram Canvas", 800, 600);
  peak_map->GetXaxis()->SetTitle("Wire Number");
  peak_map->GetYaxis()->SetTitle("Peak Amplitude (ADC)");
  peak_map->Draw("COLZ");
  peak->SaveAs(Form("%speak_map.pdf", output_dir));

  TCanvas* peakwidth = new TCanvas("peakwidth", "2D Histogram Canvas", 800, 600);
  peak_vs_width->GetXaxis()->SetTitle("Width at Half Height (ticks)");
  peak_vs_width->GetYaxis()->SetTitle("Peak Amplitude (ADC)");
  peak_vs_width->Draw("COLZ");
  peakwidth->SaveAs(Form("%sheat_map.pdf", output_dir));

  // Clean up
  delete file;
}
