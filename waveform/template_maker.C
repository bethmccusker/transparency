#include "TFile.h"
#include "TTree.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TH1F.h"
#include <iostream>

// Function to calculate the average histogram from a vector of histograms
TH1D* CalculateAverageHistogram(const std::vector<TH1D*>& histograms, int nBins) {
  TH1D* averageHistogram = new TH1D("AverageHistogram", "Average Waveform;Time (ticks);ADC counts", nBins, 1, 41);
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

// Function to calculate the scaling factor based on the integral of the histogram
double CalculateIntegral(TH1D* individualWaveform, TH1D* averageTemplate) {
  double integral_individual = individualWaveform->Integral();
  double integral_template = averageTemplate->Integral();
  // std::cout<<"Individual Integral :" <<integral_individual<<std::endl;
  // std::cout<<"Template Integral:"<<integral_template<<std::endl;

  if (integral_template == 0) {
    std::cerr << "Error: Template histogram has zero integral!" << std::endl;
    return 1.0;  // No scaling if template has zero integral
  }

  // Scaling factor is the ratio of integrals
  return integral_individual / integral_template;
}

// Function to calculate Chi-squared
double CalculateChiSquared(TH1D* individualWaveform, TH1D* scaledTemplate) {
  double chi2 = 0.0;
  int nBins = individualWaveform->GetNbinsX();

  for (int i = 1; i <= nBins; i++) {
    double O = individualWaveform->GetBinContent(i);  // Observed
    double E = scaledTemplate->GetBinContent(i);       // Expected

    if (E > 0) { // Avoid division by zero
      chi2 += (O - E) * (O - E) / E;
    }
  }

  return chi2;
}

// Function to calculate the scaling factor based on peak amplitude
double CalculateScalingFactor(TH1D* individualWaveform, TH1D* averageTemplate) {
  int centralBin = 21;

  double amplitude_individual = individualWaveform->GetBinContent(centralBin);
  double amplitude_template = averageTemplate->GetBinContent(centralBin);
  // std::cout<<"Individual amplitude"<<amplitude_individual<<std::endl;
  // std::cout<<"Template amplitude"<<amplitude_template<<std::endl;

  if (amplitude_template == 0) {
    std::cerr << "Error: Template histogram has zero amplitude at the central bin!" << std::endl;
    return 1.0;  // No scaling if template has zero amplitude
  }

  // Scaling factor is the ratio of amplitudes
  return amplitude_individual / amplitude_template;
}

// Function to compare an individual waveform to the template
void CompareWaveformToTemplate(TH1D* individualWaveform, TH1D* averageTemplate, const char* output_dir, int wave_num,int wire_num, TGraph* scalingGraph) {
  // Calculate scaling factor
  double scalingFactor = CalculateScalingFactor(individualWaveform, averageTemplate);
  std::cout << "Ratio of aplitudes for waveform " << wave_num << ": " << scalingFactor << std::endl;
  //Calculate integral scale factor
  double integralScale = CalculateIntegral(individualWaveform, averageTemplate);
  std::cout << "Ratio of integrals for waveform " << wave_num << ": " << integralScale << std::endl;
  scalingGraph->SetPoint(scalingGraph->GetN(), wire_num, scalingFactor);
  // Scale the template
  TH1D* Template = (TH1D*)averageTemplate->Clone(Form("scaled_template_%d", wave_num));
  Template->Scale(scalingFactor);  // Apply scaling factor

  // Calculate Chi-squared
  double chi2 = CalculateChiSquared(individualWaveform, Template);
  std::cout << "Chi-squared for waveform " << wave_num << ": " << chi2 << std::endl;

  // Create a canvas for comparison plot
  TCanvas *c_compare = new TCanvas(Form("c_compare_%d", wave_num), "Comparison with Template", 800, 600);

  // Draw the scaled template on the same canvas                                                                                                                                                                     
  Template->SetLineColor(kBlue-6);
  Template->SetLineWidth(2);
  Template->Draw("HIST");

  // Draw the individual waveform
  individualWaveform->SetLineColor(kRed-6);
  individualWaveform->SetLineWidth(2);
  individualWaveform->Draw("HIST SAME");

   // Add a legend to differentiate the plots
  TLegend *legend = new TLegend(0.7, 0.7, 0.9, 0.9);
  legend->AddEntry(individualWaveform, "Individual Waveform", "l");
  legend->AddEntry(Template, "Scaled Template", "l");
  legend->Draw();

  // Save the comparison plot
  c_compare->SaveAs(Form("%scomparison_waveform_%d.png", output_dir, wave_num));

  // Clean up
  delete c_compare;
  delete Template;
}


void template_maker() {
   gStyle->SetOptStat(0); //Removing the stats box
  TFile *file = TFile::Open("/exp/sbnd/data/users/bethanym/wire_transparency/testing_HD_pr/hists_decode_data_evb01_EventBuilder1_art1_run16740_10_20240912T082517-30bef869-9d29-42a0-aa77-091ad9c1620d_Reco1Comm-20\
241002T102140.root");
  TTree *tree = (TTree*)file->Get("hitdumper/hitdumpertree");

  // Variables to store the data
  std::vector<int> *waveform_number = nullptr;
  std::vector<float> *adc_on_wire = nullptr;
  std::vector<float> *time_for_waveform = nullptr;
  std::vector<int> *wire_number = nullptr;

  // Set branch addresses to connect the variables to the tree
  tree->SetBranchAddress("waveform_number", &waveform_number);
  tree->SetBranchAddress("adc_on_wire", &adc_on_wire);
  tree->SetBranchAddress("time_for_waveform", &time_for_waveform);
  tree->SetBranchAddress("wire_number", &wire_number);

  const char *output_dir = "/exp/sbnd/data/users/bethanym/wire_transparency/templateFitting/";

  // Vector to store histograms of each waveform
  int nBins = 41;
  std::vector<TH1D*> histograms;
  TGraph* scalingGraph = new TGraph();

  // Loop over tree entries (events)
  Long64_t nEntries = tree->GetEntries();
  for (Long64_t i = 0; i < nEntries; i++) {
    tree->GetEntry(i); // Load the ith entry

    // Map to store waveform data indexed by waveform number                                                                                                                                                         
    std::map<int, std::vector<double>> waveform_adc_map;
    // std::map<int, std::vector<double>> waveform_time_map;
    std::map<int, int> waveform_wire_map;

    // Loop through each data point in the current entry and organize by waveform number
    for (size_t j = 0; j < waveform_number->size(); j++) {
      // Check for valid ADC and time values (filter out -9999 default)
      if (adc_on_wire->at(j) != -9999 && time_for_waveform->at(j) != -9999) {
	// Insert ADC and time values into the corresponding waveform number
	int wave_num = waveform_number->at(j);
	waveform_adc_map[wave_num].push_back(adc_on_wire->at(j));
	//	waveform_time_map[wave_num].push_back(time_for_waveform->at(j));
	waveform_wire_map[wave_num] = wire_number->at(j);
      }
    }
  
    // For each unique waveform, create a histogram
    for (const auto &entry : waveform_adc_map) {
      int wave_num = entry.first;
      const std::vector<double>& adc_vals = entry.second;
      // const std::vector<double>& time_vals = waveform_time_map[wave_num];
      int wire_num = waveform_wire_map[wave_num];

      // Create a histogram for the waveform
      TH1D* hist = new TH1D(Form("waveform_%d", wave_num), Form("Waveform %d;Time (ticks);ADC counts", wave_num), nBins, 1, 41);

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
      //      c1->SaveAs(Form("%swaveform_%d.png", output_dir, wave_num));

      // Clean up
      delete c1;
    }
  }

  // Calculate the average histogram (template)
  TH1D* averageHistogram = CalculateAverageHistogram(histograms, nBins);

  // Plot and save the averaged histogram
  TCanvas *c_avg = new TCanvas("c_avg", "Average Waveform", 800, 600);
  averageHistogram->SetLineColor(kRed-6);
  averageHistogram->SetLineWidth(3);
  averageHistogram->Draw();
  c_avg->SaveAs(Form("%stemplate_waveform.png", output_dir));

  // Save the average histogram to a ROOT file
  TFile *outputFile = new TFile("/exp/sbnd/data/users/bethanym/wire_transparency/templateFitting/average_waveform.root", "RECREATE");
    averageHistogram->Write();
    outputFile->Close();
  
    // Compare individual waveforms to the average template and plot
    int wave_num=1;
    for (auto hist : histograms) {
      CompareWaveformToTemplate(hist, averageHistogram, output_dir, wave_num,wire_num,scalingGraph);
      wave_num++;
    }
    
    // Plot and save the scaling factor vs wire number
    TCanvas *c_scaling = new TCanvas("c_scaling", "Scaling Factor vs Wire Number", 800, 600);
    scalingGraph->SetTitle("Scaling Factor vs Wire Number;Wire Number;Scaling Factor");
    scalingGraph->SetMarkerStyle(20);
    scalingGraph->SetMarkerColor(kBlue);
    scalingGraph->Draw("AP");
    c_scaling->SaveAs(Form("%sscaling_factor_vs_wire_number.png", output_dir));

  // Clean up
  for (auto hist : histograms) {
    delete hist;
  }
  delete c_avg;
  delete outputFile;
  delete averageHistogram;
  delete file;
}
