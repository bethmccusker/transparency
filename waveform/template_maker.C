#include "TFile.h"
#include "TTree.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TH1F.h"
#include <iostream>

//*********************************** Functions Start ******************************************//
TH1D* CalculateAverageHistogram(const std::vector<std::vector<TH1D*>>& histograms, int nBins) {
  // Create the final histogram to store averages
  auto averageHistogram = new TH1D("AverageHistogram", "Average Waveform;Time (ticks);ADC counts", nBins, 1, 41);

  int numEvents = histograms.size();
  if (numEvents == 0) {
    std::cerr << "No histograms provided!" << std::endl;
    return averageHistogram;
  }

  // Accumulator for bin-wise averages
  std::vector<double> averageOfAverages(nBins, 0.0);
  int totalHistograms = 0;

  for (const auto& eventHistograms : histograms) {
    int numHistograms = eventHistograms.size();
    if (numHistograms == 0) continue; // Skip empty events

    totalHistograms += numHistograms;

    // Temporary storage for this event's averages
    std::vector<double> eventAverages(nBins, 0.0);
    for (int iBin = 1; iBin <= nBins; ++iBin) {
      double sumBinContents = 0.0;
      for (const auto& histo : eventHistograms) {
        if (!histo) continue; // Handle null pointers gracefully
        sumBinContents += histo->GetBinContent(iBin);
      }
      eventAverages[iBin - 1] = sumBinContents / numHistograms;
    }
    // Accumulate into the global average
    for (int iBin = 0; iBin < nBins; ++iBin) {
      averageOfAverages[iBin] += eventAverages[iBin];
    }
  }

  // Normalize global averages and set the bin contents
  for (int iBin = 0; iBin < nBins; ++iBin) {
    averageHistogram->SetBinContent(iBin + 1, averageOfAverages[iBin] / numEvents);
  }

  std::cout << totalHistograms << "   Number of histograms in the full template" << std::endl;
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

    if (E>0) { // Avoid division by zero
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
void CompareWaveformToTemplate(TH1D* individualWaveform, TH1D* averageTemplate, const char* output_dir, int wave_num,int wire_num, TGraph* scalingGraph, double hit_tim, TH2D* heat_map) {
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

  heat_map->Fill(wire_num, hit_tim);

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


// Function to find the intersection of two lines in the form Z = mY + c
std::pair<double, double> findIntersection(double m1, double c1, double m2, double c2) {
  // Check if the lines are parallel
  if (m1 == m2) {
    std::cerr << "The lines are parallel and do not intersect." << std::endl;
    return {std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN()};
  }

  // Calculate the intersection point
  double y_intersect = (c2 - c1) / (m1 - m2); // Solve for Y
  double z_intersect = m1 * y_intersect + c1; // Solve for Z using either line equation

  return {y_intersect, z_intersect};
}


//*********************************************Functions End ***************************************************//


//******************************************** Main Body Start ***********************************************//
void template_maker() {
  gStyle->SetOptStat(0); //Removing the stats box
   gStyle->SetPalette(kCandy);
   TColor::InvertPalette(); 
  TFile *file = TFile::Open("/exp/sbnd/data/users/bethanym/wire_transparency/filter_test/hists_decode_data_evb04_process2_EventBuilder4_p2_art1_run16740_3_20240912T081019-2106a30f-6e4b-4c39-b710-5286b5565346_Reco1Comm-20241204T105223.root");
  TTree *tree = (TTree*)file->Get("hitdumper/hitdumpertree");

  // Variables to store the data
  Int_t *nhits = nullptr;
  std::vector<int> *waveform_number = nullptr;
  std::vector<float> *adc_on_wire = nullptr;
  std::vector<float> *time_for_waveform = nullptr;
  std::vector<int> *wire_number = nullptr;
  std::vector<double> *hit_time = nullptr;
 
  // Set branch addresses to connect the variables to the tree
  tree->SetBranchAddress("nhits", &nhits); 
  tree->SetBranchAddress("waveform_number", &waveform_number);
  tree->SetBranchAddress("adc_on_wire", &adc_on_wire);
  tree->SetBranchAddress("time_for_waveform", &time_for_waveform);
  tree->SetBranchAddress("wire_number", &wire_number);
  tree->SetBranchAddress("hit_time", &hit_time);

  const char *output_dir = "/exp/sbnd/data/users/bethanym/wire_transparency/templateFitting/";

  // Vector to store histograms of each waveform
  int nBins = 41;
  std::vector<std::vector<TH1D*>> histograms;
  std::vector<std::map<int, int>> waveform_wire_map;
  std::vector<std::map<int, double>> waveform_peak_time_map;

  // Loop over tree entries (events)
  Long64_t nEntries = tree->GetEntries();
  for (Long64_t i = 0; i < nEntries; i++) {
    tree->GetEntry(i); // Load the ith entry

    // Map to store waveform data indexed by waveform number                                       
    std::map<int, std::vector<double>> waveform_adc_map;
 
    std::map<int, int> temp_waveform_wire_map;
    std::map<int, double> temp_waveform_peak_time_map;  
    // Loop through each data point in the current entry and organize by waveform number
    for (size_t j = 0; j < waveform_number->size(); j++) {
      // Check for valid ADC and time values (filter out -9999 default)
      if (adc_on_wire->at(j) != -9999 && time_for_waveform->at(j) != -9999) {
	// Insert ADC and time values into the corresponding waveform number
	int wave_num = waveform_number->at(j);
	waveform_adc_map[wave_num].push_back(adc_on_wire->at(j));
	temp_waveform_wire_map[wave_num] = wire_number->at(j);
	temp_waveform_peak_time_map[wave_num] = hit_time->at(j);     
      }
    }

    waveform_wire_map.push_back(temp_waveform_wire_map);
    waveform_peak_time_map.push_back(temp_waveform_peak_time_map);
  
    std::vector<TH1D*> temp_histograms;
    // For each unique waveform, create a histogram
    for (const auto &entry : waveform_adc_map) {
      int wave_num = entry.first;
      const std::vector<double>& adc_vals = entry.second;
     
      // Create a histogram for the waveform
      TH1D* hist = new TH1D(Form("waveform_%d_entry_%lld", wave_num, i), Form("Waveform %d;Time (ticks);ADC counts", wave_num), nBins, 1, 41);

      // Fill the histogram with ADC values at corresponding time bins
      for (size_t k = 0; k < adc_vals.size(); k++) {
	hist->SetBinContent(k+1, adc_vals[k]);
      }

      // Add the histogram to the vector
      temp_histograms.push_back(hist);

      // Create a canvas to draw the histogram
      TCanvas *c1 = new TCanvas(Form("c1_waveform_%d", wave_num), "Waveform", 800, 600);
      hist->SetLineColor(kRed-6);
      hist->SetLineWidth(3);
      hist->Draw();

      // Save the plot as a PNG file with the output directory
      //          c1->SaveAs(Form("%swaveform_%d.pdf", output_dir, wave_num));
      std::cout<<" Waveform Number "<<wave_num<<std::endl;
      // Clean up
      delete c1;
    }
    histograms.push_back(temp_histograms);
  }

 
  TH1D* averageHistogram = CalculateAverageHistogram(histograms, nBins);

 
  // Plot and save the averaged histogram
  TCanvas *c_avg = new TCanvas("c_avg", "Average Waveform", 800, 600);
  averageHistogram->SetLineColor(kBlue-6);
  averageHistogram->SetLineWidth(3);
  averageHistogram->Draw();
  c_avg->SaveAs(Form("%stemplate_waveform.pdf", output_dir));


  std::vector<double> eventScalingFactors;  
  // Compare individual waveforms to the average template and plot                                                                                                                                          
for (int i = 0; i < nEntries; i++) {
  tree->GetEntry(i);
  if (nhits == 0) {
    std::cout << "no waveforms or hits for entry " << i << std::endl;
    continue;
  } 

  TGraph* scalingGraph = new TGraph();
  TH2D* heat_map = new TH2D("h_wire_vs_peak_time", "Wire Number vs Peak Time",  1700, 0, 1700, 3500,0, 3500);
          
  int wave_num=1;
  double sumScalingFactors = 0.0;
  int countScalingFactors = 0;
  for (auto hist : histograms.at(i)) {
    double hit_tim=waveform_peak_time_map.at(i)[wave_num];
    int wire_num = waveform_wire_map.at(i)[wave_num];
    CompareWaveformToTemplate(hist, averageHistogram, output_dir, wave_num,wire_num,scalingGraph,hit_tim, heat_map);
   
    double scalingFactor = CalculateScalingFactor(hist, averageHistogram);
    sumScalingFactors += scalingFactor;
    countScalingFactors++;

    wave_num++;
  }
    
  if (countScalingFactors > 0) {
    double averageScalingFactor = sumScalingFactors / countScalingFactors;
    eventScalingFactors.push_back(averageScalingFactor);
    std::cout << "Average scaling factor for entry " << i << ": " << averageScalingFactor << std::endl;
  }

  // Plot and save the scaling factor vs wire number
  TCanvas *c_scaling = new TCanvas("c_scaling", "Scaling Factor vs Wire Number", 800, 600);
  scalingGraph->SetTitle("Scaling Factor vs Wire Number;Wire Number;Scaling Factor");
  scalingGraph->SetMarkerStyle(8);
  scalingGraph->SetMarkerColor(kRed-6);
  scalingGraph->SetMarkerSize(0.75);
  scalingGraph->GetXaxis()->SetLimits(0, 1700);
  scalingGraph->Draw("AP");
  //  c_scaling->SaveAs(Form("%sscaling_factor_vs_wire_number_entry_%d.pdf", output_dir, i));

  TCanvas* c_wire_vs_time = new TCanvas("c_wire_vs_time", "Wire vs Peak Time", 900, 600);
  heat_map->GetXaxis()->SetTitle("Wire Number");
  heat_map->GetYaxis()->SetTitle("Waveform Time (Ticks)");
  // heat_map->GetYaxis()->SetRangeUser(500, 800);
  heat_map->Draw("COLZ");
  //  c_wire_vs_time->SaveAs(Form("%sheat_map_entry_%d.pdf", output_dir, i));

  // Clean up
  delete scalingGraph;
  delete heat_map;
 
 }
}


