// norm plots with mean value and binning control

#include <TFile.h>
#include <TH1.h>
#include <THStack.h>
#include <TKey.h>
#include <TDirectory.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TEfficiency.h>
#include <TGraphAsymmErrors.h>
#include <iostream>
#include <vector>
#include <string>
#include <unordered_map>
#include <map>
#include <algorithm>
#include <TRatioPlot.h>
#include <TGaxis.h>
#include <sys/stat.h> // For mkdir

using namespace std;

struct ProcessInfo {
    string path;
    bool isSignal;
    int color; // stacked histogram color
    string legendName;
};

struct BinSettings {
    bool useFixedBinCount;  // true = use fixed number of bins, false = use rebinFactor
    int fixedBinCount;      // Desired number of bins (when useFixedBinCount = true)
    int rebinFactor;        // Traditional rebin factor (when useFixedBinCount = false)
};

struct HistogramSetting {
    string typeOfHisto;
    string title;
    string xAxisTitle;
    string saveName;
    pair<double, double> xRange; // x-axis range (min, max)
    BinSettings binSettings;     // Binning configuration
};

const string plotExtension = ".png"; // save file extension
const string savePath = "Norm_plots";

// Function to create the directory if it does not exist
void CreateDirectoryIfNotExists(const string& path) {
    struct stat info;
    if (stat(path.c_str(), &info)) {
        // The directory does not exist, so create it
        if (mkdir(path.c_str(), 0755)) {
            cerr << "Error: Could not create directory " << path << endl;
        } else {
            cout << "Directory created: " << path << endl;
        }
    } else if (!(info.st_mode & S_IFDIR)) {
        // The path exists, but it is not a directory
        cerr << "Error: " << path << " already exists, but it is not a directory." << endl;
    }
}

void setLatexSetting(TLatex& histoLatex, const string& histo_title) {
    histoLatex.SetNDC();
    histoLatex.SetTextSize(0.05);
    histoLatex.SetTextAlign(11);
    histoLatex.DrawLatex(0.1, 0.96, "CMS #scale[0.85]{#font[52]{Private Work}}");

    histoLatex.SetTextSize(0.04);
    histoLatex.SetTextAlign(31); // Align to the right
    histoLatex.DrawLatex(0.9, 0.945, "2024 year, 108.96 fb^{-1} [13.6 TeV]");

    histoLatex.SetTextSize(0.04);
    histoLatex.SetTextAlign(11);
    string title = "[ " + histo_title + " ]";
    histoLatex.DrawLatex(0.13, 0.92, title.c_str());
}

TGraphAsymmErrors* CreateRatioPlot(TH1* signal, THStack* background, const HistogramSetting& setting, double maxRatioLimit = 1.5) {
    if (!signal || !background) {
        cerr << "Error: Signal or background is null while creating the ratio plot." << endl;
        return nullptr;
    }

    auto First_Stacked_histo = (TH1*)background->GetHists()->First();
    if (!First_Stacked_histo) {
        cerr << "Error: No histogram found in the THStack." << endl;
        return nullptr;
    }

    // Use the user-defined x-range if specified, otherwise use the histogram's range
    double xAxisRange_low = setting.xRange.first != 0 || setting.xRange.second != 0 ? 
                           setting.xRange.first : First_Stacked_histo->GetXaxis()->GetXmin();
    double xAxisRange_high = setting.xRange.first != 0 || setting.xRange.second != 0 ? 
                            setting.xRange.second : First_Stacked_histo->GetXaxis()->GetXmax();

    double maxVal = -std::numeric_limits<double>::max();
    double minVal = std::numeric_limits<double>::min();

    TGraphAsymmErrors* grRatio = new TGraphAsymmErrors();

    // Calculate sum of background THStack
    TH1D* sumHist = static_cast<TH1D*>(signal->Clone("sumHist"));
    sumHist->Reset();
    TList* histList = background->GetHists();
    TIter next(histList);
    TH1D* hist;
    while ((hist = static_cast<TH1D*>(next()))) {
        sumHist->Add(hist);
    }

    // ratio histogram and its error
    for (int i = 1; i <= signal->GetNbinsX(); ++i) {
        double x = signal->GetBinCenter(i);
        
        // Skip points outside the specified X-axis range
        if ((setting.xRange.first != 0 || setting.xRange.second != 0) && 
            (x < setting.xRange.first || x > setting.xRange.second)) {
            continue;
        }

        double S = signal->GetBinContent(i);
        double B = sumHist->GetBinContent(i);
        double sigma_S = signal->GetBinError(i);
        double sigma_B = sumHist->GetBinError(i);

        if (B > 0 && S > 0) {
            double R = S / B;
            // Truncate R if it exceeds the maximum limit
            if (R > maxRatioLimit) {
                R = maxRatioLimit;
            }

            double errorLow = R * sqrt(pow(sigma_S / S, 2) + pow(sigma_B / B, 2));
            double errorHigh = errorLow;

            int iPoint = grRatio->GetN();
            grRatio->SetPoint(iPoint, x, R);
            grRatio->SetPointError(iPoint, 0.0, 0.0, errorLow, errorHigh);
            maxVal = max(maxVal, R + errorHigh);
            minVal = min(minVal, R - errorLow);
        }
    }

    // Configure the ratio plot
    grRatio->GetXaxis()->SetTitle(setting.xAxisTitle.c_str());
    grRatio->GetYaxis()->SetTitle("#bf{Ratio of [ #it{Sig / Bkg} ]}");
    grRatio->GetYaxis()->SetTitleOffset(0.7);
    grRatio->GetXaxis()->SetTitleSize(0.085);
    grRatio->GetYaxis()->SetLabelSize(0.08);
    grRatio->GetYaxis()->SetTitleSize(0.07);
    grRatio->GetYaxis()->SetNdivisions(505);

    // X-axis settings to match the main histogram
    grRatio->GetXaxis()->SetLimits(xAxisRange_low, xAxisRange_high); // Same x-axis range
    grRatio->GetXaxis()->SetNdivisions(505); // Same number of divisions
    grRatio->GetXaxis()->SetTickLength(0.03); // Same tick length
    grRatio->GetXaxis()->SetLabelSize(0.04); // Same label size

    // Adjust the y-axis scale to avoid very high values
    double yMin = 0.0; // Fixed lower limit
    double yMax = maxRatioLimit; // Adjustable upper limit

    grRatio->GetYaxis()->SetRangeUser(yMin, yMax); // Set y-axis limits

    // Configure the Y-axis to use scientific notation
    grRatio->GetYaxis()->SetMoreLogLabels(); // Enable more labels on the Y-axis
    grRatio->GetYaxis()->SetNoExponent(false); // Force scientific notation
    grRatio->GetYaxis()->SetMaxDigits(3); // Set the maximum number of digits in the exponent

    grRatio->SetLineWidth(4);
    grRatio->SetLineColor(kBlue + 2);
    grRatio->SetMarkerStyle(20);
    grRatio->SetMarkerSize(1);
    grRatio->SetMarkerColor(kBlue);
    grRatio->SetLineWidth(2);

    return grRatio;
}

void DrawStackedHistograms(const vector<pair<TH1*, ProcessInfo>>& histograms, const HistogramSetting& setting, const string& histName) {
    if (histograms.empty()) {
        cerr << "Error: No histograms provided to draw." << endl;
        return;
    }

    THStack* stack = new THStack("stack", "");
    if (!stack) {
        cerr << "Error: Failed to create THStack." << endl;
        return;
    }

    TH1* signalHist = nullptr;
    double maxY = 0;

    // Add background histograms to the stack
    for (const auto& histPair : histograms) {
        TH1* hist = histPair.first;
        if (!hist) {
            cerr << "Error: Null histogram found for process " << histPair.second.legendName << endl;
            continue;
        }

        TH1* processedHist = hist;
        
        // Apply rebinning according to settings
        if (setting.binSettings.useFixedBinCount) {
            // Fixed number of bins method
            int originalBins = hist->GetNbinsX();
            int rebinFactor = originalBins / setting.binSettings.fixedBinCount;
            if (rebinFactor < 1) rebinFactor = 1;
            
            processedHist = (TH1*)hist->Rebin(rebinFactor, Form("%s_rebinned", hist->GetName()));
        } else {
            // Traditional rebinning method
            if (setting.binSettings.rebinFactor > 1) {
                processedHist = (TH1*)hist->Rebin(setting.binSettings.rebinFactor, 
                                                Form("%s_rebinned", hist->GetName()));
            }
        }

        // Apply X-axis range if specified
        if (setting.xRange.first != 0 || setting.xRange.second != 0) {
            processedHist->GetXaxis()->SetRangeUser(setting.xRange.first, setting.xRange.second);
        }
        
        // Normalize to unit area
        if (processedHist->Integral() != 0) {
            processedHist->Scale(1.0 / processedHist->Integral());
        }

        if (histPair.second.isSignal) {
            signalHist = processedHist;
        } else {
            processedHist->SetFillColor(0); // No fill
            processedHist->SetLineColor(histPair.second.color);
            processedHist->SetLineWidth(2);
            stack->Add(processedHist);
            cout << "Background histogram added to stack: " << histPair.second.legendName << endl;
        }

        maxY = max(maxY, processedHist->GetMaximum());
    }

    if (!signalHist) {
        cerr << "Error: No signal histogram found." << endl;
        delete stack;
        return;
    }

    // Create canvas and pads
    TCanvas* canvas = new TCanvas("canvas", setting.title.c_str(), 1800, 1400);
    TPad* upperPad = new TPad("upperPad", "Upper Pad", 0.0, 0.30, 1.0, 1.0);
    TPad* lowerPad = new TPad("lowerPad", "Lower Pad", 0.0, 0.0, 1.0, 0.30);

    // Increase the right margin to create more space
    upperPad->SetRightMargin(0.20); // Increased right margin to 20%
    lowerPad->SetRightMargin(0.20); // Increased right margin to 20%

    upperPad->Draw();
    lowerPad->Draw();

    // Draw upper pad
    upperPad->cd();
    if (!stack->GetHists() || stack->GetHists()->GetEntries() == 0) {
        cerr << "Error: No histograms were added to the stack." << endl;
        delete stack;
        delete canvas;
        return;
    }

    stack->Draw("HIST");
    stack->GetXaxis()->SetTitle(setting.xAxisTitle.c_str());
    stack->GetYaxis()->SetTitle("#bf{Normalized Events / bin}");
    stack->SetMaximum(maxY * 10.0);
    stack->SetMinimum(1e-3);
    stack->GetXaxis()->SetLabelSize(0);

    // Use the user-defined x-range if specified, otherwise use the histogram's range
    double xAxisRange_low = setting.xRange.first != 0 || setting.xRange.second != 0 ? 
                           setting.xRange.first : stack->GetXaxis()->GetXmin();
    double xAxisRange_high = setting.xRange.first != 0 || setting.xRange.second != 0 ? 
                            setting.xRange.second : stack->GetXaxis()->GetXmax();

    // Set the x-axis range
    stack->GetXaxis()->SetLimits(xAxisRange_low, xAxisRange_high);
    stack->GetXaxis()->SetNdivisions(505);
    stack->GetXaxis()->SetTickLength(0.03);
    stack->GetXaxis()->SetLabelSize(0.04);

    // Draw the signal as a line (no fill)
    signalHist->SetLineColor(kBlack);
    signalHist->SetLineWidth(5);
    signalHist->SetFillStyle(0); // No fill
    signalHist->Draw("HIST same");

    upperPad->SetLogy(1);

    TLatex latex;
    setLatexSetting(latex, setting.title);

    // Add colored squares with the color and name of each process in the top-right corner
    double xText = 0.85; // Horizontal position of the text
    double yText = 0.85; // Initial vertical position of the text (higher up)
    double yStep = 0.05; // Spacing between text lines (increased to create more space)
    double xStep = 0.15; // Horizontal spacing between processes in the same line

    for (const auto& histPair : histograms) {
        TH1* hist = histPair.first;
        if (!hist) continue;

        double mean = hist->GetMean();
        double stdDev = hist->GetStdDev();

        // Format values to two decimal places
        stringstream meanStream, stdDevStream;
        meanStream << std::fixed << std::setprecision(2) << mean;
        stdDevStream << std::fixed << std::setprecision(2) << stdDev;

        // Create a TPaveText for the colored square
        double boxSize = 0.02; // Size of the square (height and width)
        TPaveText* box = new TPaveText(xText - 0.04, yText - boxSize / 2, xText - 0.04+ boxSize, yText + boxSize / 2, "NDC");
        box->SetFillColor(histPair.second.color); // Use the color defined in ProcessInfo
        box->SetLineColor(kBlack); // Set the border color to black
        box->SetBorderSize(1); // Border thickness
        box->Draw();

        // First line: process name and mean value
        string info1 = histPair.second.legendName + ": #mu = " + meanStream.str();
        latex.SetTextColor(kBlack); // Text color is black for better readability
        latex.SetTextSize(0.025); // Reduced font size to 0.025
        latex.DrawLatex(xText, yText, info1.c_str());
        yText -= yStep; // Move to the next line

        // Second line: standard deviation
        string info2 = "#sigma = " + stdDevStream.str();
        latex.DrawLatex(xText, yText, info2.c_str());
        yText -= yStep; // Move to the next line
    }

    // Draw lower pad
    lowerPad->cd();
    lowerPad->SetGrid(1, 1);
    lowerPad->SetTickx(1);

    auto* grRatio = CreateRatioPlot(signalHist, stack, setting);
    if (grRatio) {
        // X-axis settings to align with the main plot
        grRatio->GetXaxis()->SetLimits(xAxisRange_low, xAxisRange_high); // Same x-axis range
        grRatio->GetXaxis()->SetNdivisions(505); // Same number of divisions
        grRatio->GetXaxis()->SetTickLength(0.03); // Same tick length
        grRatio->GetXaxis()->SetLabelSize(0.04); // Same label size
        grRatio->Draw("AP");
    }

    // Save the plot
    string saveName = savePath + "/" + histName + "_" + setting.saveName + plotExtension;
    CreateDirectoryIfNotExists(savePath); // Create the directory if it does not exist
    canvas->SaveAs(saveName.c_str());

    // Clean up
    delete stack;
    delete canvas;
}

void Iteration_Directories_And_Histograms(const unordered_map<string, ProcessInfo>& processes, const unordered_map<string, vector<HistogramSetting>>& histogramSettings) {
    for (const auto& histSettingPair : histogramSettings) {
        const string& channelName = histSettingPair.first;
        const vector<HistogramSetting>& settings = histSettingPair.second;

        cout << "Processing folder: " << channelName << endl;

        for (const auto& setting : settings) {
            vector<pair<TH1*, ProcessInfo>> histogramsForStacking;

            cout << "Looking for histogram: " << setting.typeOfHisto << " in folder " << channelName << endl;

            for (const auto& processPair : processes) {
                const string& processName = processPair.first;
                const ProcessInfo& processInfo = processPair.second;

                cout << "Opening file: " << processInfo.path << " for process: " << processName << endl;

                TFile* file = TFile::Open(processInfo.path.c_str(), "READ");
                if (!file || file->IsZombie()) {
                    cerr << "Error: Could not open file " << processInfo.path << endl;
                    continue;
                }

                cout << "File opened successfully." << endl;

                // Navigate to the folder (Lepton or jet)
                TDirectory* dir = nullptr;
                file->GetObject(channelName.c_str(), dir);
                if (!dir) {
                    cerr << "Error: Folder not found: " << channelName << " in " << processInfo.path << endl;
                    file->Close();
                    continue;
                }

                cout << "Folder found: " << channelName << endl;

                // Look for the histogram inside the folder
                TH1* hist = nullptr;
                dir->GetObject(setting.typeOfHisto.c_str(), hist);
                if (!hist) {
                    cerr << "Error: Histogram not found: " << setting.typeOfHisto << " in " << channelName << endl;
                    file->Close();
                    continue;
                }

                cout << "Histogram found: " << setting.typeOfHisto << endl;

                hist->SetDirectory(0); // Unlink the histogram from the file
                histogramsForStacking.push_back(make_pair(hist, processInfo));
                file->Close();
            }

            if (histogramsForStacking.empty()) {
                cerr << "Error: No histograms were loaded for type: " << setting.typeOfHisto << endl;
                continue;
            }

            // Draw the histograms
            DrawStackedHistograms(histogramsForStacking, setting, channelName);

            // Clean up
            for (auto& histPair : histogramsForStacking) {
                delete histPair.first;
            }
        }
    }
}

int Ploter_Norm() {
    unordered_map<string, ProcessInfo> processes = {
{"TTSL", {"TTSL.root", false, kOrange+1, "TTSL"}},
//{"TTZ", {"TTZ.root", false, kAzure+1, "TTZ"}},
{"TTH", {"TTH.root", false, kPink+1, "TTH"}},
{"TT4b", {"TT4b.root", false, kSpring+3, "TT4b"}},
{"TTbb", {"TTbb.root", false, kViolet+1, "TTbb"}},
{"TTZH", {"TTZH.root", false, kTeal+1, "TTZH"}},
{"TTZZ", {"TTZZ.root", false, kRed+1, "TTZZ"}},
{"TTHH", {"TTHH.root", true, kBlack, "TTHH"}},
        // Add other processes as needed
    };

    // Define histogram settings with binning options
    unordered_map<string, vector<HistogramSetting>> histogramSettings = {
      {"Tree", 
            {
                // Format: {histName, title, xTitle, saveName, {xMin,xMax}, {useFixedBins, fixedBinCount, rebinFactor}}
                {"cutflow", "Cutflow", "Cuts", "cutflow", {0, 0}, {false, 10, 1}},
                {"cutflow_w", "Cutflow Weighted", "Cuts", "cutflow_w", {0, 0}, {false, 10, 1}}          
            }
        },
        {"Lepton", 
            {
                // Format: {histName, title, xTitle, saveName, {xMin,xMax}, {useFixedBins, fixedBinCount, rebinFactor}}
                {"lepCharge1", "Lepton Charge 1", "Charge", "lepCharge1", {0, 0}, {false, 10, 1}},
                {"lepCharge2", "Lepton Charge 2", "Charge", "lepCharge2", {0, 0}, {false, 10, 1}},
                {"LepNumber", "Lepton Number", "Number of Leptons", "LepNumber", {0, 0}, {false, 5, 1}},
                {"ElecNumber", "Electron Number", "Number of Electrons", "ElecNumber", {0, 0}, {false, 10, 1}},
                {"MuonNumber", "Muon Number", "Number of Muons", "MuonNumber", {0, 0}, {false, 5, 1}},
                {"elePT1", "Electron PT 1", "pT [GeV]", "elePT1", {0, 800}, {true, 10, 1}},
                {"elePT2", "Electron PT 2", "pT [GeV]", "elePT2", {0, 300}, {true, 10, 1}},
                {"muonPT1", "Muon PT 1", "pT [GeV]", "muonPT1", {0, 800}, {true, 10, 1}},
                {"muonPT2", "Muon PT 2", "pT [GeV]", "muonPT2", {0, 300}, {true, 10, 1}},
                {"leptonHT", "Lepton HT", "HT [GeV]", "leptonHT", {0, 3500}, {false, 10, 1}}
            }
        },
              {"jet", 
            {
                {"jetPT1", "Jet PT 1", "pT [GeV]", "jetPT1", {0, 1500}, {true, 10, 1}},
                {"jetPT2", "Jet PT 2", "pT [GeV]", "jetPT2", {0, 1500}, {true, 10, 1}},
                {"jetPT3", "Jet PT 3", "pT [GeV]", "jetPT3", {0, 1000}, {true, 10, 1}},
                {"jetPT4", "Jet PT 4", "pT [GeV]", "jetPT4", {0, 500}, {true, 10, 1}},
                {"jetPT5", "Jet PT 5", "pT [GeV]", "jetPT5", {0, 500}, {true, 10, 1}},
                {"jetPT6", "Jet PT 6", "pT [GeV]", "jetPT6", {0, 400}, {true, 10, 1}},
                {"bjetPT1", "B-Jet PT 1", "pT [GeV]", "bjetPT1", {0, 1500}, {false, 10, 1}},
                {"bjetPT2", "B-Jet PT 2", "pT [GeV]", "bjetPT2", {0, 1500}, {false, 10, 1}},
                {"bjetPT3", "B-Jet PT 3", "pT [GeV]", "bjetPT3", {0, 1000}, {false, 10, 1}},
                {"bjetPT4", "B-Jet PT 4", "pT [GeV]", "bjetPT4", {0, 400}, {false, 10, 1}},
                {"bjetPT5", "B-Jet PT 5", "pT [GeV]", "bjetPT5", {0, 300}, {false, 10, 1}},
                {"bjetPT6", "B-Jet PT 6", "pT [GeV]", "bjetPT6", {0, 200}, {false, 10, 1}},
                {"jetHT", "Jet HT", "HT [GeV]", "jetHT", {0, 2000}, {false, 10, 1}},
                {"jetBHT", "B-Jet HT", "HT [GeV]", "jetBHT", {0, 2000}, {false, 10, 1}},
                {"met", "Missing ET", "ET [GeV]", "met", {0, 1000}, {false, 10, 1}},
                {"jetNumber", "Jet Number", "Number of Jets", "jetNumber", {0, 0}, {false, 10, 1}},
                {"jetBNumber", "B-Jet Number", "Number of B-Jets", "jetBNumber", {0, 0}, {false, 10, 1}},
                {"invMass_HH1Matched", "Invariant Mass HH1 Matched", "Mass [GeV]", "invMass_HH1Matched", {0, 500}, {true, 10, 1}},
                {"invMass_HH2Matched", "Invariant Mass HH2 Matched", "Mass [GeV]", "invMass_HH2Matched", {0, 500}, {true, 10, 1}},
                {"invMass_HH1NotMatched", "Invariant Mass HH1 Not Matched", "Mass [GeV]", "invMass_HH1NotMatched", {0, 500}, {true, 10, 1}},
                {"invMass_HH2NotMatched", "Invariant Mass HH2 Not Matched", "Mass [GeV]", "invMass_HH2NotMatched", {0, 500}, {true, 10, 1}}
            }
        }
    };

    Iteration_Directories_And_Histograms(processes, histogramSettings);
    return 0;
}
//to run use: cmssw-el7  cmsenv  root -l -b -q Ploter_Norm.cpp
// to run use root -l -b -q Ploter_Norm.cpp
