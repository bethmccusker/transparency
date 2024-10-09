#include "TFile.h"
#include "TTree.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TH1F.h"
#include <iostream>
#include <map>
#include <vector>

void waveform_plotter() {
  // Open the ROOT file
  TFile *file = TFile::Open("/exp/sbnd/data/users/bethanym/wire_transparency/testing_HD_pr/hists_decode_data_evb01_EventBuilder1_art1_run16740_10_20240912T082517-30bef869-9d29-42a0-aa77-091ad9c1620d_Reco1Comm-20241002T102140.root");

  // Access the tree containing your waveform data
  TTree *tree = (TTree*)file->Get("hitdumper/hitdumpertree");

  // Variables to store the data
  std::vector<int> *waveform_number = nullptr;
  std::vector<float> *adc_on_wire = nullptr;
  std::vector<float> *time_for_waveform = nullptr;

  // Set branch addresses to connect the variables to the tree
  tree->SetBranchAddress("waveform_number", &waveform_number);
  tree->SetBranchAddress("adc_on_wire", &adc_on_wire);
  tree->SetBranchAddress("time_for_waveform", &time_for_waveform);

  const char *output_dir = "/exp/sbnd/data/users/bethanym/wire_transparency/templateFitting/";

  // Map to store waveform data indexed by waveform number
  std::map<int, std::vector<double>> waveform_adc_map;
  std::map<int, std::vector<double>> waveform_time_map;

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

    // Create a TGraph for the current waveform
    TGraph *graph = new TGraph(adc_vals.size(), &time_vals[0], &adc_vals[0]);
    graph->SetTitle(Form("Waveform Number %d;Time (ticks);ADC counts", wave_num));
    graph->SetLineColor(kRed-6);
    graph->SetLineWidth(3);

    // Create a canvas to draw the graph
    TCanvas *c1 = new TCanvas(Form("c1_waveform_%d", wave_num), "Waveform", 800, 600);
    graph->Draw("AL");

    // Save the plot as a PNG file
    c1->SaveAs(Form("%swaveform_%d.png",output_dir, wave_num));

    // Optionally, delete canvas and graph to free memory
    delete c1;
    delete graph;
  }

  // Clean up
  delete file;
}
