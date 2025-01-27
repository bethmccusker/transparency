#include "TFile.h"
#include "TTree.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TMath.h"
#include <iostream>

//***********************************Functions Start ******************************************//
TH1D* CalculateAverageHistogram(const std::vector<std::vector<TH1D*>>& histograms, int nBins) {
  // Create the final histogram to store averages
  auto averageHistogram = new TH1D("AverageHistogram", "Average Waveform;Time (ticks);ADC counts", nBins, 1, 81);

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

TH1D* AlignToDownPeak(TH1D* waveform, int targetBin) {
  if (!waveform) return nullptr;

  // Find the bin with the minimum (down-peak) value
  int downPeakBin = waveform->GetMinimumBin();
  if (downPeakBin == 0) return waveform; // No valid minimum found, return original

  // Calculate the shift needed to align the down-peak with the target bin
  int shift = targetBin - downPeakBin;

  // Create a new histogram to store the shifted waveform
  TH1D* alignedWaveform = (TH1D*)waveform->Clone(Form("%s_aligned", waveform->GetName()));
  alignedWaveform->Reset(); // Clear contents

  // Shift bin contents
  int nBins = waveform->GetNbinsX();
  for (int iBin = 1; iBin <= nBins; ++iBin) {
    int shiftedBin = iBin + shift;
    if (shiftedBin >= 1 && shiftedBin <= nBins) {
      alignedWaveform->SetBinContent(shiftedBin, waveform->GetBinContent(iBin));
    }
  }

  return alignedWaveform;
}


// Function to calculate the scaling factor based on the integral of the histogram
std::pair<double, double> CalculateIntegralPeakTrough(TH1D* individualWaveform, TH1D* averageTemplate) {
  int nBins = individualWaveform->GetNbinsX();
  int centerBin =41;

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
  int centerBin = 41;  // Central bin for dividing peak/trough

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
void CompareWaveformToTemplate(TH1D* individualWaveform, TH1D* averageTemplate, const char* output_dir, int wave_num,int wire_num, TGraph* scalingGraph, double hit_tim, TH2D* heat_map, TH2D* scaling_factors ) {
  // Calculate separate scaling factors for peak (before central bin) and trough (after central bin)
  std::pair<double, double> scalingFactors = CalculatePeakTroughScalingFactor(individualWaveform, averageTemplate);

  // Extract the scaling factors
  double peakScalingFactor = scalingFactors.first;
  double troughScalingFactor = scalingFactors.second;

  // Log the scaling factors
   std::cout << "Scaling factor for peak region of waveform " << wave_num << ": " << peakScalingFactor << std::endl;
   std::cout << "Scaling factor for trough region of waveform " << wave_num << ": " << troughScalingFactor << std::endl;

  scaling_factors->Fill(peakScalingFactor, troughScalingFactor);

  // Optionally update the scaling graph with the peak or trough scaling factor (you can decide which one to use here)
  scalingGraph->SetPoint(scalingGraph->GetN(), wire_num, peakScalingFactor);

  // Clone the template to apply different scaling for peak and trough
  TH1D* Template = (TH1D*)averageTemplate->Clone(Form("scaled_template_%d", wave_num));

  // Define the central bin
  int nBins = individualWaveform->GetNbinsX();
  int centerBin = 41;

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
  /*
  TLegend *legend = new TLegend(0.7, 0.7, 0.9, 0.9);
  legend->AddEntry(individualWaveform, "Individual Waveform", "l");
  legend->AddEntry(Template, "Scaled Template", "l");
  legend->Draw();
  */
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

double vector_to_double(const std::vector<double>& vec) {
  if (vec.size() != 1) {
    // Handle the case where the vector doesn't have exactly one element
    // You can throw an exception, return a default value, or handle it differently
    throw std::runtime_error("Vector does not contain exactly one element."); 
  }
  return vec[0]; 
}

//**********************************Functions End********************************************//
//**********************************Main Body Start ********************************************//
void induction_1_template_maker() {
  gStyle->SetOptStat(0); //Removing the stats box
   gStyle->SetPalette(kCandy);
   TColor::InvertPalette();
  TFile *file = TFile::Open("/exp/sbnd/data/users/bethanym/wire_transparency/filter_test/hists_decode_data_evb04_process2_EventBuilder4_p2_art1_run16740_3_20240912T081019-2106a30f-6e4b-4c39-b710-5286b5565346_Reco1Comm-20241204T112223.root");
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

  TFile *file2=TFile::Open("/exp/sbnd/app/users/bethanym/wire_plane_transparency/code_repo/waveform/wires/wire_info.root", "READ");
  TTree *tree2 =(TTree*)file2->Get("data");
  std::vector<int> *chan=nullptr;
  std::vector<double> *wire_gradient=nullptr;
  std::vector<double> *wire_intercept=nullptr;

  tree2->SetBranchAddress("chan", &chan);
  tree2->SetBranchAddress("gradient", &wire_gradient);
  tree2->SetBranchAddress("intercept", &wire_intercept);

  const char *output_dir = "/exp/sbnd/data/users/bethanym/wire_transparency/templateFitting/";

  // Vector to store histograms of each waveform
  int nBins = 81;
  std::vector<std::vector<TH1D*>> histograms;
  std::vector<std::map<int, int>> waveform_wire_map;
  std::vector<std::map<int, int>> waveform_channel_map;
  std::vector<std::map<int, double>> waveform_peak_time_map;

  TH2D* h_sum = new TH2D("sum entriws", "2D Histogram;Y coordinate;Z coordinate;Scaling Probability",175, 0., 500.,100, -200., 200.);
  TH2D* h_count = new TH2D("count entries", "Histogram;Y coordinate;Z coordinate;Scaling Probability",175, 0., 500.,100, -200., 200.);

  // Loop over tree entries (events)
  Long64_t nEntries = tree->GetEntries();
  for (Long64_t i = 0; i < nEntries; i++) {
    tree->GetEntry(i); // Load the ith entry
  
    // Map to store waveform data indexed by waveform number                                                                                                                                                         
    std::map<int, std::vector<double>> waveform_adc_map;
    std::map<int, int> temp_waveform_wire_map;
    std::map<int, int> temp_waveform_channel_map;
    std::map<int, double> temp_waveform_peak_time_map;
    // Loop through each data point in the current entry and organize by waveform number
    for (size_t j = 0; j < waveform_number->size(); j++) {
      // Check for valid ADC and time values (filter out -9999 default)
      if (adc_on_wire->at(j) != -9999 && time_for_waveform->at(j) != -9999) {
	// Insert ADC and time values into the corresponding waveform number
	int wave_num = waveform_number->at(j);
	waveform_adc_map[wave_num].push_back(adc_on_wire->at(j));
	//	waveform_time_map[wave_num].push_back(time_for_waveform->at(j));
	temp_waveform_wire_map[wave_num] = wire_number->at(j);
	temp_waveform_channel_map[wave_num] = channel_number->at(j);
	temp_waveform_peak_time_map[wave_num] = hit_time->at(j);
      }
    }

    waveform_wire_map.push_back(temp_waveform_wire_map);
    waveform_channel_map.push_back(temp_waveform_channel_map);
    waveform_peak_time_map.push_back(temp_waveform_peak_time_map);
  
    std::vector<TH1D*> temp_histograms;
    // For each unique waveform, create a histogram
    for (const auto &entry : waveform_adc_map) {
      int wave_num = entry.first;
      const std::vector<double>& adc_vals = entry.second;
    

      // Create a histogram for the waveform
      TH1D* hist = new TH1D(Form("waveform_%d_entry_%lld", wave_num, i), Form("Waveform %d;Time (ticks);ADC counts", wave_num), nBins, 1, 81);

      // Fill the histogram with ADC values at corresponding time bins
      for (size_t k = 0; k < adc_vals.size(); k++) {
	hist->SetBinContent(k+1, adc_vals[k]);
      }
      // Align the waveform to the down-peak
      TH1D* alignedHist = AlignToDownPeak(hist, /* targetBin */ 47);

      // Add the aligned histogram to the list for averaging
      temp_histograms.push_back(alignedHist);

      // Create a canvas to draw the histogram
      TCanvas *c1 = new TCanvas(Form("c1_waveform_%d", wave_num), "Waveform", 800, 600);
      hist->SetLineColor(kRed-6);
      hist->SetLineWidth(3);
      hist->Draw();

      // Save the plot as a PNG file with the output directory
      //      c1->SaveAs(Form("%swaveform_%d.pdf", output_dir, wave_num));
      std::cout<<" Waveform Number "<<wave_num<<std::endl;
      // Clean up
      delete c1;
    }
    histograms.push_back(temp_histograms);
  }

  // Calculate the average histogram (template)
  TH1D* averageHistogram = CalculateAverageHistogram(histograms, nBins);
 
  // Plot and save the averaged histogram
  TCanvas *c_avg = new TCanvas("c_avg", "Average Waveform", 800, 600);
  averageHistogram->SetLineColor(kBlue-6);
  averageHistogram->SetLineWidth(3);
  averageHistogram->Draw();
  c_avg->SaveAs(Form("%stemplate_waveform.pdf", output_dir));

  std::vector<double> eventScalingFactorsPeak;
  std::vector<double> eventScalingFactorsTrough;  
  // Compare individual waveforms to the average template and plot
  for (int i = 0; i < nEntries; i++) {
    tree->GetEntry(i);
    if (nhits == 0) {
      std::cout << "no waveforms or hits for entry " << i << std::endl;
      continue;
    }

    TGraph* scalingGraph = new TGraph();
    TH2D* heat_map = new TH2D("h_wire_vs_peak_time", "Wire Number vs Peak Time",  500, 0, 2000, 500, 0, 3500);
    TH2D* scaling_factors = new TH2D("peak", "trough",  100, 0, 18, 100, 0, 6);
    double crt_grad=vector_to_double(*crt_gradient);
    double crt_int=vector_to_double(*crt_intercept);

    int wave_num=1;
    double sumScalingFactorsPeak=0.0;
    double sumScalingFactorsTrough = 0.0;
    int countScalingFactors = 0; 
   for (auto hist : histograms.at(i)) {
      double hit_tim=waveform_peak_time_map.at(i)[wave_num];
      int channel_num = waveform_channel_map.at(i)[wave_num]; 
      int wire_num = waveform_wire_map.at(i)[wave_num];
      //  std::cout << "wire number for  waveform " << wave_num << ": " << wire_num << std::endl;
      CompareWaveformToTemplate(hist, averageHistogram, output_dir, wave_num,wire_num,scalingGraph,hit_tim, heat_map, scaling_factors);
      std::pair<double, double> scalingFactors = CalculatePeakTroughScalingFactor(hist, averageHistogram);
      double peakScalingFactor = scalingFactors.first;
      double troughScalingFactor = scalingFactors.second;
      sumScalingFactorsPeak += peakScalingFactor;
      sumScalingFactorsTrough += troughScalingFactor;
      countScalingFactors++;
      int nentries=tree2->GetEntries();
      for (int entry = 0; entry < nentries; ++entry) {
	tree2->GetEntry(entry);
	for (size_t i = 0; i < chan->size(); ++i) {
	  if ((*chan)[i] == channel_num) {
	    std::cout<<"chan" <<chan->at(i)<<std::endl;
	    double gradient = wire_gradient->at(i);
	    double intercept = wire_intercept->at(i);
	    std::pair<double, double> yz_position=findIntersection(crt_grad,crt_int ,gradient,intercept);
	    double y_position = yz_position.first;
	    double z_position= yz_position.second;
	    std::cout<< "y=  "<< y_position<<std::endl;
	    std::cout<< "z=  "<< z_position<<std::endl;
	    h_sum->Fill(z_position, y_position, scalingFactor);
	    h_count->Fill(z_position,y_position,1);
	  }
	}
      }
      wave_num++;
    }
  
   if (countScalingFactors > 0) {
     double averageScalingFactorPeak = sumScalingFactorsPeak / countScalingFactors;
     double averageScalingFactorTrough = sumScalingFactorsTrough / countScalingFactors;
     eventScalingFactorsPeak.push_back(averageScalingFactorPeak);
     eventScalingFactorsTrough.push_back(averageScalingFactorTrough);
     std::cout << "Average peak scaling factor for entry " << i << ": " << averageScalingFactorPeak << std::endl;
     std::cout << "Average trough scaling factor for entry " << i << ": " << averageScalingFactorTrough << std::endl;
   }
    

   TH2D* h_avg = (TH2D*)h_sum->Clone("h_avg");
   h_avg->SetTitle("Averaged Values");
   h_avg->Divide(h_count);
   TCanvas* average = new TCanvas("averagw", "2D Histogram Canvas", 800, 600);
   h_avg->GetXaxis()->SetTitle("Z");
   h_avg->GetYaxis()->SetTitle("Y");
   h_avg->Draw("COLZ");
   average->SaveAs(Form("%saverage_scaling_across_dectector_%d.pdf", output_dir, i));
 
   TCanvas* sum = new TCanvas("sum", "2D Histogram Canvas", 800, 600);
   h_sum->GetXaxis()->SetTitle("Z");
   h_sum->GetYaxis()->SetTitle("Y");
   h_sum->Draw("COLZ");
   sum->SaveAs(Form("%ssum_scaling_across_dectector.pdf", output_dir));

    // Plot and save the scaling factor vs wire number
    TCanvas *c_scaling = new TCanvas("c_scaling", "Scaling Factor vs Wire Number", 800, 600);
    scalingGraph->SetTitle("Scaling Factor vs Wire Number;Wire Number;Scaling Factor");
    scalingGraph->SetMarkerStyle(8);
    scalingGraph->SetMarkerColor(kRed-6);
    scalingGraph->SetMarkerSize(0.75);
    scalingGraph->GetXaxis()->SetLimits(0, 1700);
    scalingGraph->Draw("AP");
    // c_scaling->SaveAs(Form("%sscaling_factor_vs_wire_number_entry_%d.pdf", output_dir, i));
    
    TCanvas* c_wire_vs_time = new TCanvas("c_wire_vs_time", "Wire vs Peak Time", 900, 600);
    heat_map->Draw("COLZ");
    heat_map->GetXaxis()->SetTitle("Wire Number");
    heat_map->GetYaxis()->SetTitle("Waveform Time (Ticks)");
    // c_wire_vs_time->SaveAs(Form("%sheat_map_entry_%d.pdf", output_dir, i));
  

    TCanvas* c_scaling_factors = new TCanvas("", "", 900, 600);
    scaling_factors->Draw("COLZ");
    scaling_factors->GetXaxis()->SetTitle("Peak Scaling");
    scaling_factors->GetYaxis()->SetTitle("Trough Scaling");
    // c_scaling_factors->SaveAs(Form("%sscaling_factors_entry_%d.pdf", output_dir, i));

    delete scalingGraph;
    delete heat_map;
    delete scaling_factors;

    //heat_map->Reset();
    //scaling_factors->Reset();
    //scalingGraph->Set(0);
  }

  //   Clean up
  // for (auto hist : histograms) {
    // delete hist;
    //  }

  //delete c_avg;
  //delete averageHistogram;
  //delete file;
}
