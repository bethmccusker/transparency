#include "TFile.h"
#include "TTree.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TMath.h"
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
std::pair<double, double> CalculateIntegralPeakTrough(TH1D* individualWaveform, TH1D* averageTemplate) {
  int nBins = individualWaveform->GetNbinsX();
  int centerBin =21;

  // Calculate integral for peak region (1 to centerBin)
  double peak_individual = individualWaveform->Integral(1, centerBin);
  double peak_template = averageTemplate->Integral(1, centerBin);

  // Calculate integral for trough region (centerBin to nBins)
  double trough_individual = individualWaveform->Integral(centerBin, nBins);
  double trough_template = averageTemplate->Integral(centerBin, nBins);

  // Handle edge case: if any template region has zero integral, return 1.0 as the scaling factor
  if (peak_template == 0 || trough_template == 0) {
    std::cerr << "Error: Template histogram has zero integral in peak or trough region!" << std::endl;
    return std::make_pair(1.0, 1.0);   // No scaling if either region has zero integral
  }

  // Scaling factor is a weighted average of the ratios of peak and trough integrals
  double peak_scaling = peak_individual / peak_template;
  double trough_scaling = trough_individual / trough_template;

  return std::make_pair(peak_scaling, trough_scaling);
}


// Function to calculate Chi-squared
double CalculateChiSquared(TH1D* individualWaveform, TH1D* scaledTemplate) {
  double chi2 = 0.0;
  int nBins = individualWaveform->GetNbinsX();

  for (int i = 1; i <= nBins; i++) {
    double O = individualWaveform->GetBinContent(i);  // Observed
    double E = scaledTemplate->GetBinContent(i);       // Expected

    if (!(E == 0.0)) { // Avoid division by zero
      chi2 += (O - E) * (O - E) / TMath::Abs(E);
    }
  }

  return chi2;
}

// Function to calculate the scaling factor based on peak amplitude
std::pair<double, double> CalculatePeakTroughScalingFactor(TH1D* individualWaveform, TH1D* averageTemplate) {
  int nBins = individualWaveform->GetNbinsX();
  int centerBin = 21;  // Central bin for dividing peak/trough

  // Identify peak and trough bins in both individual and template histograms
  double peak_amplitude_individual = individualWaveform->GetMaximum();
  double peak_amplitude_template = averageTemplate->GetMaximum();

  double trough_amplitude_individual = individualWaveform->GetMinimum();
  double trough_amplitude_template = averageTemplate->GetMinimum();

  if (peak_amplitude_template == 0) {
    std::cerr << "Error: Template histogram has zero peak amplitude!" << std::endl;
    return std::make_pair(1.0, 1.0);  // Return default scaling factors in case of error
  }

  if (trough_amplitude_template == 0) {
    std::cerr << "Error: Template histogram has zero trough amplitude!" << std::endl;
    return std::make_pair(1.0, 1.0);  // Return default scaling factors in case of error
  }

  // Calculate scaling factor as the ratio of amplitudes for peak and trough
  double peak_scaling = peak_amplitude_individual / peak_amplitude_template;
  double trough_scaling = trough_amplitude_individual / trough_amplitude_template;

  // Return both scaling factors as a pair
  return std::make_pair(peak_scaling, trough_scaling);
}


// Function to compare an individual waveform to the template
void CompareWaveformToTemplate(TH1D* individualWaveform, TH1D* averageTemplate, const char* output_dir, int wave_num,int wire_num, TGraph* scalingGraph) {
  // Calculate separate scaling factors for peak (before central bin) and trough (after central bin)
  std::pair<double, double> scalingFactors = CalculatePeakTroughScalingFactor(individualWaveform, averageTemplate);

  // Extract the scaling factors
  double peakScalingFactor = scalingFactors.first;
  double troughScalingFactor = scalingFactors.second;

  // Log the scaling factors
  std::cout << "Scaling factor for peak region of waveform " << wave_num << ": " << peakScalingFactor << std::endl;
  std::cout << "Scaling factor for trough region of waveform " << wave_num << ": " << troughScalingFactor << std::endl;

  // Optionally update the scaling graph with the peak or trough scaling factor (you can decide which one to use here)
  scalingGraph->SetPoint(scalingGraph->GetN(), wire_num, peakScalingFactor);

  // Clone the template to apply different scaling for peak and trough
  TH1D* Template = (TH1D*)averageTemplate->Clone(Form("scaled_template_%d", wave_num));

  // Define the central bin
  int nBins = individualWaveform->GetNbinsX();
  int centerBin = 21;

  // Scale the peak region (1 to centerBin)
  for (int iBin = 1; iBin <= centerBin; ++iBin) {
    double originalContent = Template->GetBinContent(iBin);
    Template->SetBinContent(iBin, originalContent * peakScalingFactor);
  }

  // Scale the trough region (centerBin to nBins)
  for (int iBin = centerBin + 1; iBin <= nBins; ++iBin) {
    double originalContent = Template->GetBinContent(iBin);
    Template->SetBinContent(iBin, originalContent * troughScalingFactor);
  } 

  // Calculate Chi-squared
  double chi2 = CalculateChiSquared(individualWaveform, Template);
  std::cout << "Chi-squared for waveform " << wave_num << ": " << chi2 << std::endl;

  // Create a canvas for comparison plot
  TCanvas *c_compare = new TCanvas(Form("c_compare_%d", wave_num), "Comparison with Template", 800, 600);


  // Draw the individual waveform                                                                                                                                                                                    
  individualWaveform->SetLineColor(kRed-6);
  individualWaveform->SetLineWidth(2);
  individualWaveform->Draw("HIST");

  // Draw the scaled template on the same canvas                                                                                                                                                                     
  Template->SetLineColor(kBlue-6);
  Template->SetLineWidth(2);
  Template->Draw("HIST SAME");

   // Add a legend to differentiate the plots
  TLegend *legend = new TLegend(0.7, 0.7, 0.9, 0.9);
  legend->AddEntry(individualWaveform, "Individual Waveform", "l");
  legend->AddEntry(Template, "Scaled Template", "l");
  legend->Draw();

  // Save the comparison plot
      c_compare->SaveAs(Form("%scomparison_waveform_%d.pdf", output_dir, wave_num));

  // Clean up
  delete c_compare;
  delete Template;
}


void induction_template_maker() {
   gStyle->SetOptStat(0); //Removing the stats box
  TFile *file = TFile::Open("/exp/sbnd/data/users/bethanym/wire_transparency/hd_variable_test/hists_decode_data_evb01_EventBuilder1_art1_run16740_10_20240912T082517-30bef869-9d29-42a0-aa77-091ad9c1620d_Reco1Comm-20241017T220428.root");
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
  std::map<int, int> waveform_wire_map;
  // Loop over tree entries (events)
  Long64_t nEntries = tree->GetEntries();
  for (Long64_t i = 0; i < nEntries; i++) {
    tree->GetEntry(i); // Load the ith entry

    // Map to store waveform data indexed by waveform number                                                                                                                                                         
    std::map<int, std::vector<double>> waveform_adc_map;
    // std::map<int, std::vector<double>> waveform_time_map;
   
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
      //      c1->SaveAs(Form("%swaveform_%d.pdf", output_dir, wave_num));

      // Clean up
      delete c1;
    }
  }

  // Calculate the average histogram (template)
  TH1D* averageHistogram = CalculateAverageHistogram(histograms, nBins);

  // Plot and save the averaged histogram
  TCanvas *c_avg = new TCanvas("c_avg", "Average Waveform", 800, 600);
  averageHistogram->SetLineColor(kBlue-6);
  averageHistogram->SetLineWidth(3);
  averageHistogram->Draw();
  c_avg->SaveAs(Form("%stemplate_waveform.pdf", output_dir));

  
  // Compare individual waveforms to the average template and plot
  int wave_num=1;
  for (auto hist : histograms) {
    int wire_num = waveform_wire_map[wave_num];
    CompareWaveformToTemplate(hist, averageHistogram, output_dir, wave_num,wire_num,scalingGraph);
    wave_num++;
  }
    
  // Plot and save the scaling factor vs wire number
  TCanvas *c_scaling = new TCanvas("c_scaling", "Scaling Factor vs Wire Number", 800, 600);
  scalingGraph->SetTitle("Scaling Factor vs Wire Number;Wire Number;Scaling Factor");
  scalingGraph->SetMarkerStyle(8);
  scalingGraph->SetMarkerColor(kRed-6);
  scalingGraph->SetMarkerSize(0.75);
  scalingGraph->GetXaxis()->SetLimits(0, 1700);
  scalingGraph->Draw("AP");
  c_scaling->SaveAs(Form("%sscaling_factor_vs_wire_number.pdf", output_dir));
  
  // Clean up
  for (auto hist : histograms) {
    delete hist;
  }
  delete c_avg;
  delete averageHistogram;
  delete file;
}
