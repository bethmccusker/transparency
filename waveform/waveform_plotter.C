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

  // Map to store waveform data indexed by waveform number
  std::map<int, std::vector<double>> waveform_adc_map;
  std::map<int, std::vector<double>> waveform_time_map;
  std::map<int, std::vector<double>> waveform_wire_map;

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

    // If no data, skip this waveform
    if (adc_vals.empty() || time_vals.empty()) {
      continue;
    }

    auto [half_width, amplitude] = CalculateHalfWidthHeightAndAmplitudeCollection(adc_vals, time_vals);
    std::cout << "Waveform " << wave_num << " Half-Width Height: " << half_width << " Amplitude: " << amplitude << "\n";


    auto [trough_half_width, trough_amplitude] = CalculateHalfWidthHeightAndAmplitudeInduction(adc_vals, time_vals);
    std::cout << " | Trough Half-Width: " << trough_half_width << " Trough Amplitude: " << trough_amplitude << "\n";

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

  // Clean up
  delete file;
}
