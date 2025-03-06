#include "TFile.h"
#include "TTree.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TF1.h"
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


double CalculateAreaUnderCurve(const std::vector<double>& adc_vals, const std::vector<double>& time_vals) {
  if (adc_vals.size() < 2 || time_vals.size() < 2 || adc_vals.size() != time_vals.size()) return -1.0;

  double area = 0.0;
  for (size_t i = 1; i < adc_vals.size(); ++i) {
    double delta_time = time_vals[i] - time_vals[i - 1];
    double average_adc = (adc_vals[i] + adc_vals[i - 1]) / 2.0;
    area += average_adc * delta_time;
  }

  return area;
}

bool HasDoublePeakFeature(const std::vector<double>& adc_vals, double threshold) {
  int peak_count = 0;
  for (size_t i = 1; i < adc_vals.size() - 1; ++i) {
    if (adc_vals[i] > threshold && adc_vals[i] > adc_vals[i - 1] && adc_vals[i] > adc_vals[i + 1]) {
      peak_count++;
      if (peak_count > 1) return true;
    }
  }
  return false;
}


void waveform_plotter() {
    gStyle->SetOptStat(0); //Removing the stats box
  gStyle->SetPalette(kCandy);
  TColor::InvertPalette();
  // Open the ROOT file
  TFile *file = TFile::Open("/exp/sbnd/data/users/bethanym/wire_transparency/filter_test/hists_decode_data_evb01_EventBuilder1_art1_run16740_10_20240912T082517-30bef869-9d29-42a0-aa77-091ad9c1620d_Reco1Comm-20250220T101021.root");

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
  TH2D* width_map = new TH2D("h_wire_vs_width", "",  250, 0, 1700, 100,0, 15);
  TH2D* peak_map = new TH2D("h_wire_vs_peak_amplitude", "",  250, 0, 1700, 100,0, 300);
  TH2D* peak_vs_width = new TH2D("h_peak_vs_width", "",  100, 0,300, 100,0, 15);
  TH2D* area_under_curve_map = new TH2D("h_wire_vs_area", "",  250, 0, 1700, 100,0, 2500);
  TH2D* area_under_curve_vs_width = new TH2D("h_wire_vs_area_vs_width", "",  100, 0, 15, 100,0, 2500);
  TH2D* area_under_curve_vs_amplitude = new TH2D("h_wire_vs_area_vs_amplitude", "",  100, 0, 300, 100,0, 2500);
  TH1F *amplitude_hist = new TH1F("amplitude_hist", "Amplitude Distribution;Amplitude (ADC counts);Entries", 100, 0,300);


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

  int pass_count=0;
  int fail_count=0;

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

    if (HasDoublePeakFeature(adc_vals, 10.0)){
      fail_count++;
      continue;
    }

    double max_adc = *std::max_element(adc_vals.begin(), adc_vals.end());
    if (max_adc > 140) {
      fail_count++;
      continue;
    }

  

    auto [half_width, amplitude] = CalculateHalfWidthHeightAndAmplitudeCollection(adc_vals, time_vals);
    if (half_width < 5.5 || half_width > 8.5) {
      fail_count++;
      continue;
    }

    pass_count++;
    double area_under_curve = CalculateAreaUnderCurve(adc_vals, time_vals);

    //    std::cout << "Waveform " << wave_num << " Half-Width Height: " << half_width << " Amplitude: " << amplitude << "\n";
    // std::cout << "area under curve" <<area_under_curve<< endl;
    //std::cout << " Wire " <<wire_num[0]<< endl; 
    // std::cout << " Time " <<time_vals[0]<< endl;
    // std::cout << " Wave " <<wave_num<< endl;

	//	auto [trough_half_width, trough_amplitude] = CalculateHalfWidthHeightAndAmplitudeInduction(adc_vals, time_vals);
	//  std::cout << " | Trough Half-Width: " << trough_half_width << " Trough Amplitude: " << trough_amplitude << "\n";

    heat_map->Fill(wire_num[0],time_vals[0]);
    width_map->Fill(wire_num[0],half_width);
    peak_map->Fill(wire_num[0],amplitude);
    peak_vs_width->Fill(amplitude, half_width);
    area_under_curve_map->Fill(wire_num[0],area_under_curve);
    area_under_curve_vs_width->Fill(half_width,area_under_curve);
    area_under_curve_vs_amplitude->Fill(amplitude,area_under_curve);
    amplitude_hist->Fill(amplitude);

    // Create a TGraph for the current waveform
    TGraph *graph = new TGraph(adc_vals.size(), &time_vals[0], &adc_vals[0]);
    graph->SetTitle(Form("Waveform Number %d;Time (ticks);ADC counts", wave_num));
    graph->SetLineColor(kRed-6);
    graph->SetLineWidth(3);

    // Create a canvas to draw the graph
    TCanvas *c1 = new TCanvas(Form("c1_waveform_%d", wave_num), "Waveform", 800, 600);
    graph->Draw("AL");


    if((half_width>8.5 || half_width<5.5)/* || (amplitude>140) || ((half_width>8.5 || half_width<5.5)&& amplitude>140)*/){
    c1->SaveAs(Form("%swaveform_%d.pdf",output_dir, wave_num));
    std::cout << "Waveform " << wave_num << " Half-Width Height: " << half_width << " Amplitude: " << amplitude << "\n";
    std::cout << "area under curve" <<area_under_curve<< endl;
    std::cout << " Wire " <<wire_num[0]<< endl;
    std::cout << " Time " <<time_vals[0]<< endl;
    std::cout << " Wave " <<wave_num<< endl;
    }
    //    std::cout << "Number of waveforms that passed: " << pass_count << std::endl;
    std::cout << "Number of waveforms that failed: " << fail_count << std::endl;
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
  peak_vs_width->GetYaxis()->SetTitle("Width at Half Height (ticks)");
  peak_vs_width->GetXaxis()->SetTitle("Peak Amplitude (ADC)");
  peak_vs_width->Draw("COLZ");
  peakwidth->SaveAs(Form("%speakwidth_map.pdf", output_dir));

  TCanvas* area = new TCanvas("area", "2D Histogram Canvas", 800, 600);
  area_under_curve_map->GetXaxis()->SetTitle("Wire Number");
  area_under_curve_map->GetYaxis()->SetTitle("Area Under Curve (ADC*ticks)");
  area_under_curve_map->Draw("COLZ");
  area->SaveAs(Form("%sarea_under_area_map.pdf", output_dir));

  TCanvas* areawidth = new TCanvas("areawidth", "2D Histogram Canvas", 800, 600);
  area_under_curve_vs_width->GetXaxis()->SetTitle("Width at Half Height (ticks)");
  area_under_curve_vs_width->GetYaxis()->SetTitle("Area Under Curve (ADC*ticks)");
  area_under_curve_vs_width->Draw("COLZ");
  areawidth->SaveAs(Form("%sarea_under_area_vs_width.pdf", output_dir));

  TCanvas* areapeak = new TCanvas("areapeak", "2D Histogram Canvas", 800, 600);
  area_under_curve_vs_amplitude->GetXaxis()->SetTitle("Peak Amplitude (ADC)");
  area_under_curve_vs_amplitude->GetYaxis()->SetTitle("Area Under Curve (ADC*ticks)");
  area_under_curve_vs_amplitude->Draw("COLZ");
  areapeak->SaveAs(Form("%sarea_under_area_vs_amplitude.pdf", output_dir));

  TCanvas *c2 = new TCanvas("c2", "Amplitude Histogram", 800, 600);
  amplitude_hist->Draw();
  double fit_min = amplitude_hist->GetXaxis()->GetXmax();
  double fit_max = amplitude_hist->GetXaxis()->GetXmin();
  for (int bin = 1; bin <= amplitude_hist->GetNbinsX(); ++bin) {
    if (amplitude_hist->GetBinContent(bin) > 65) {
      double bin_center = amplitude_hist->GetBinCenter(bin);
      if (bin_center < fit_min) fit_min = bin_center;
      if (bin_center > fit_max) fit_max = bin_center;
    }
  }

  if (fit_max > fit_min) {
    TF1 *gaus_fit = new TF1("gaus_fit", "gaus", fit_min, fit_max);
    amplitude_hist->Fit(gaus_fit, "R");
  }

  c2->SaveAs(Form("%samplitude_histogram.pdf", output_dir));


  // Clean up
  delete file;
}
