#include "TFile.h"
#include "TTree.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TH1F.h"
#include <iostream>

// Function to calculate the average histogram from a vector of histograms
TH1D* CalculateAverageHistogram(const std::vector<TH1D*>& histograms) {
  if (histograms.empty()) {
    return nullptr;  // Return null pointer if vector is empty
  }

  int nBins = histograms[0]->GetNbinsX();
  double xMin = histograms[0]->GetXaxis()->GetXmin();
  double xMax = histograms[0]->GetXaxis()->GetXmax();

  TH1D* averageHistogram = new TH1D("AverageHistogram", "Average Waveform;Time (ticks);ADC counts", nBins, xMin, xMax);
  int numHistograms = histograms.size();

  for (int iBin = 1; iBin <= nBins; ++iBin) {
    double sumBinContents = 0.0;
    for (const auto& histo : histograms) {
      sumBinContents += histo->GetBinContent(iBin);
    }
    double averageBinContent = sumBinContents / numHistograms;
    averageHistogram->SetBinContent(iBin, averageBinContent);
  }
  std::cout<<numHistograms<<"    Number of histograms in template"<<std::endl;
  return averageHistogram;
}


void waveform_plotter() {
  gStyle->SetOptStat(0); //Removing the stats box
  TFile *file = TFile::Open("/exp/sbnd/data/users/bethanym/wire_transparency/testing_HD_pr/hists_decode_data_evb01_EventBuilder1_art1_run16740_10_20240912T082517-30bef869-9d29-42a0-aa77-091ad9c1620d_Reco1Comm-20\
241002T102140.root");
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

  // Vector to store histograms of each waveform
  std::vector<TH1D*> histograms;

  // Loop over tree entries (events)
  Long64_t nEntries = tree->GetEntries();
  for (Long64_t i = 0; i < nEntries; i++) {
    tree->GetEntry(i); // Load the ith entry

    // Map to store waveform data indexed by waveform number                                                                                                                                                         
    std::map<int, std::vector<double>> waveform_adc_map;
    std::map<int, std::vector<double>> waveform_time_map;

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
  
    // For each unique waveform, create a histogram
    for (const auto &entry : waveform_adc_map) {
      int wave_num = entry.first;
      const std::vector<double>& adc_vals = entry.second;
      const std::vector<double>& time_vals = waveform_time_map[wave_num];

      // Create a histogram for the waveform
      TH1D* hist = new TH1D(Form("waveform_%d", wave_num), Form("Waveform %d;Time (ticks);ADC counts", wave_num), time_vals.size(), time_vals.front(), time_vals.back());

      // Fill the histogram with ADC values at corresponding time bins
      for (size_t k = 0; k < adc_vals.size(); k++) {
	hist->SetBinContent(k+1, adc_vals[k]);
      }

      // Add the histogram to the vector
      histograms.push_back(hist);

      // Create a canvas to draw the histogram
      TCanvas *c1 = new TCanvas(Form("c1_waveform_%d", wave_num), "Waveform", 800, 600);
      hist->SetLineColor(kRed-6);
      hist->SetLineWidth(3);
      hist->Draw();

      // Save the plot as a PNG file with the output directory
      //      c1->SaveAs(Form("%swaveform_%d.jpg", output_dir, wave_num));

      // Clean up
      delete c1;
    }
  }
  // Calculate the average histogram (template)
  TH1D* averageHistogram = CalculateAverageHistogram(histograms);

  // Plot and save the averaged histogram
  TCanvas *c_avg = new TCanvas("c_avg", "Average Waveform", 800, 600);
  averageHistogram->SetLineColor(kRed-6);
  averageHistogram->SetLineWidth(3);
  averageHistogram->Draw();
  c_avg->SaveAs(Form("%stemplate_waveform.jpg", output_dir));

  // Clean up
  for (auto hist : histograms) {
    delete hist;
  }
  delete averageHistogram;
  delete file;
}
