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
const string savePath = "Plots_all_norm";

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
    stack->SetMinimum(1e-4);
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

void Iteration_Directories_And_Histograms(const vector<pair<string, ProcessInfo>>& processes, const unordered_map<string, vector<HistogramSetting>>& histogramSettings) {    for (const auto& histSettingPair : histogramSettings) {
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

int Plotter_all_norm() {
	// Substitua a sua definição de 'processes' por esta:
vector<pair<string, ProcessInfo>> processes = {
    // Backgrounds na ordem de empilhamento desejada (de baixo para cima)
    {"TTZZ", {"/eos/user/g/gvian/root_files/TTZZ.root", false, kRed+1, "TTZZ"}},
    {"TTZH", {"/eos/user/g/gvian/root_files/TTZH.root", false, kTeal+1, "TTZH"}},
    {"TTH", {"/eos/user/g/gvian/root_files/TTH.root", false, kPink+1, "TTH"}},
    {"TTbb", {"/eos/user/g/gvian/root_files/TTbb.root", false, kViolet+1, "TTbb"}},
    {"TT4b", {"/eos/user/g/gvian/root_files/TT4b.root", false, kSpring+3, "TT4b"}},
    {"TTSL", {"/eos/user/g/gvian/root_files/TTSL.root", false, kOrange+1, "TTSL"}},
    
    // Sinal (geralmente processado por último e não adicionado ao stack)
    {"TTHH", {"/eos/user/g/gvian/root_files/TTHH.root", true, kBlack, "TTHH"}}
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
        {"ST", "Scalar Transverse Momentum", "ST [GeV]", "ST", {0, 2000}, {false, 10, 1}},                  // ADICIONADO
        {"leptonPT1", "Lepton PT 1", "pT [GeV]", "leptonPT1", {0, 800}, {false, 10, 1}},
        {"leptonEta1", "Lepton Eta 1", "#eta", "leptonEta1", {-3, 3}, {false, 10, 1}},                       // ADICIONADO
        {"elePT1", "Electron PT 1", "pT [GeV]", "elePT1", {0, 800}, {true, 10, 1}},
        {"eleEta1", "Electron Eta 1", "#eta", "eleEta1", {-3, 3}, {false, 10, 1}},                           // ADICIONADO
        {"elePT2", "Electron PT 2", "pT [GeV]", "elePT2", {0, 300}, {true, 10, 1}},
        {"muonPT1", "Muon PT 1", "pT [GeV]", "muonPT1", {0, 800}, {true, 10, 1}},
        {"muonEta1", "Muon Eta 1", "#eta", "muonEta1", {-3, 3}, {false, 10, 1}},                           // ADICIONADO
        {"muonPT2", "Muon PT 2", "pT [GeV]", "muonPT2", {0, 300}, {true, 10, 1}},
        {"leptonHT", "Lepton HT", "HT [GeV]", "leptonHT", {0, 3500}, {false, 10, 1}}
    }        },
	{"jet", {
    // --- Jatos Gerais (1-8) ---
    {"jetPT1", "Jet PT 1", "pT [GeV]", "jetPT1", {0, 1500}, {true, 10, 1}},
    {"jetEta1", "Jet Eta 1", "#eta", "jetEta1", {-3, 3}, {false, 10, 1}},         // ADICIONADO
    {"jetBTagDisc1", "Jet 1 btagDisc", "btagDisc", "jetBTagDisc1", {0, 1}, {false, 10, 1}}, // ADICIONADO
    {"jetPT2", "Jet PT 2", "pT [GeV]", "jetPT2", {0, 1500}, {true, 10, 1}},
    {"jetEta2", "Jet Eta 2", "#eta", "jetEta2", {-3, 3}, {false, 10, 1}},         // ADICIONADO
    {"jetBTagDisc2", "Jet 2 btagDisc", "btagDisc", "jetBTagDisc2", {0, 1}, {false, 10, 1}}, // ADICIONADO
    {"jetPT3", "Jet PT 3", "pT [GeV]", "jetPT3", {0, 1000}, {true, 10, 1}},
    {"jetEta3", "Jet Eta 3", "#eta", "jetEta3", {-3, 3}, {false, 10, 1}},         // ADICIONADO
    {"jetBTagDisc3", "Jet 3 btagDisc", "btagDisc", "jetBTagDisc3", {0, 1}, {false, 10, 1}}, // ADICIONADO
    {"jetPT4", "Jet PT 4", "pT [GeV]", "jetPT4", {0, 500}, {true, 10, 1}},
    {"jetEta4", "Jet Eta 4", "#eta", "jetEta4", {-3, 3}, {false, 10, 1}},         // ADICIONADO
    {"jetBTagDisc4", "Jet 4 btagDisc", "btagDisc", "jetBTagDisc4", {0, 1}, {false, 10, 1}}, // ADICIONADO
    {"jetPT5", "Jet PT 5", "pT [GeV]", "jetPT5", {0, 500}, {true, 10, 1}},
    {"jetEta5", "Jet Eta 5", "#eta", "jetEta5", {-3, 3}, {false, 10, 1}},         // ADICIONADO
    {"jetBTagDisc5", "Jet 5 btagDisc", "btagDisc", "jetBTagDisc5", {0, 1}, {false, 10, 1}}, // ADICIONADO
    {"jetPT6", "Jet PT 6", "pT [GeV]", "jetPT6", {0, 400}, {true, 10, 1}},
    {"jetEta6", "Jet Eta 6", "#eta", "jetEta6", {-3, 3}, {false, 10, 1}},         // ADICIONADO
    {"jetBTagDisc6", "Jet 6 btagDisc", "btagDisc", "jetBTagDisc6", {0, 1}, {false, 10, 1}}, // ADICIONADO
    {"jetPT7", "Jet PT 7", "pT [GeV]", "jetPT7", {0, 400}, {true, 10, 1}},         // ADICIONADO
    {"jetEta7", "Jet Eta 7", "#eta", "jetEta7", {-3, 3}, {false, 10, 1}},         // ADICIONADO
    {"jetBTagDisc7", "Jet 7 btagDisc", "btagDisc", "jetBTagDisc7", {0, 1}, {false, 10, 1}}, // ADICIONADO
    {"jetPT8", "Jet PT 8", "pT [GeV]", "jetPT8", {0, 400}, {true, 10, 1}},         // ADICIONADO
    {"jetEta8", "Jet Eta 8", "#eta", "jetEta8", {-3, 3}, {false, 10, 1}},         // ADICIONADO
    {"jetBTagDisc8", "Jet 8 btagDisc", "btagDisc", "jetBTagDisc8", {0, 1}, {false, 10, 1}}, // ADICIONADO

    // --- Jatos Leves (1-6) ---
    {"lightJetPT1", "Light Jet PT 1", "pT [GeV]", "lightJetPT1", {0, 400}, {true, 10, 1}},
    {"lightJetEta1", "Light Jet Eta 1", "#eta", "lightJetEta1", {-3, 3}, {false, 10, 1}}, // ADICIONADO
    {"lightJetBTagDisc1", "Light Jet 1 btagDisc", "btagDisc", "lightJetBTagDisc1", {0, 1}, {false, 10, 1}}, // CORRIGIDO
    {"lightJetPT2", "Light Jet PT 2", "pT [GeV]", "lightJetPT2", {0, 400}, {true, 10, 1}},
    {"lightJetEta2", "Light Jet Eta 2", "#eta", "lightJetEta2", {-3, 3}, {false, 10, 1}}, // ADICIONADO
    {"lightJetBTagDisc2", "Light Jet 2 btagDisc", "btagDisc", "lightJetBTagDisc2", {0, 1}, {false, 10, 1}}, // ADICIONADO
    {"lightJetPT3", "Light Jet PT 3", "pT [GeV]", "lightJetPT3", {0, 400}, {true, 10, 1}},
    {"lightJetEta3", "Light Jet Eta 3", "#eta", "lightJetEta3", {-3, 3}, {false, 10, 1}}, // ADICIONADO
    {"lightJetBTagDisc3", "Light Jet 3 btagDisc", "btagDisc", "lightJetBTagDisc3", {0, 1}, {false, 10, 1}}, // ADICIONADO
    {"lightJetPT4", "Light Jet PT 4", "pT [GeV]", "lightJetPT4", {0, 400}, {true, 10, 1}},
    {"lightJetEta4", "Light Jet Eta 4", "#eta", "lightJetEta4", {-3, 3}, {false, 10, 1}}, // ADICIONADO
    {"lightJetBTagDisc4", "Light Jet 4 btagDisc", "btagDisc", "lightJetBTagDisc4", {0, 1}, {false, 10, 1}}, // ADICIONADO
    {"lightJetPT5", "Light Jet PT 5", "pT [GeV]", "lightJetPT5", {0, 400}, {true, 10, 1}},
    {"lightJetEta5", "Light Jet Eta 5", "#eta", "lightJetEta5", {-3, 3}, {false, 10, 1}}, // ADICIONADO
    {"lightJetBTagDisc5", "Light Jet 5 btagDisc", "btagDisc", "lightJetBTagDisc5", {0, 1}, {false, 10, 1}}, // ADICIONADO
    {"lightJetPT6", "Light Jet PT 6", "pT [GeV]", "lightJetPT6", {0, 400}, {true, 10, 1}},
    {"lightJetEta6", "Light Jet Eta 6", "#eta", "lightJetEta6", {-3, 3}, {false, 10, 1}}, // ADICIONADO
    {"lightJetBTagDisc6", "Light Jet 6 btagDisc", "btagDisc", "lightJetBTagDisc6", {0, 1}, {false, 10, 1}}, // ADICIONADO

    // --- B-Jets (1-6) ---
    {"bjetPT1", "B-Jet PT 1", "pT [GeV]", "bjetPT1", {0, 1500}, {true, 10, 1}},
    {"bjetEta1", "B-Jet Eta 1", "#eta", "bjetEta1", {-3, 3}, {false, 10, 1}},         // ADICIONADO
    {"bjetBTagDisc1", "B-Jet 1 btagDisc", "btagDisc", "bjetBTagDisc1", {0, 1}, {false, 10, 1}}, // ADICIONADO
    {"bjetPT2", "B-Jet PT 2", "pT [GeV]", "bjetPT2", {0, 1500}, {true, 10, 1}},
    {"bjetEta2", "B-Jet Eta 2", "#eta", "bjetEta2", {-3, 3}, {false, 10, 1}},         // ADICIONADO
    {"bjetBTagDisc2", "B-Jet 2 btagDisc", "btagDisc", "bjetBTagDisc2", {0, 1}, {false, 10, 1}}, // ADICIONADO
    {"bjetPT3", "B-Jet PT 3", "pT [GeV]", "bjetPT3", {0, 1000}, {true, 10, 1}},
    {"bjetEta3", "B-Jet Eta 3", "#eta", "bjetEta3", {-3, 3}, {false, 10, 1}},         // ADICIONADO
    {"bjetBTagDisc3", "B-Jet 3 btagDisc", "btagDisc", "bjetBTagDisc3", {0, 1}, {false, 10, 1}}, // ADICIONADO
    {"bjetPT4", "B-Jet PT 4", "pT [GeV]", "bjetPT4", {0, 400}, {true, 10, 1}},
    {"bjetEta4", "B-Jet Eta 4", "#eta", "bjetEta4", {-3, 3}, {false, 10, 1}},         // ADICIONADO
    {"bjetBTagDisc4", "B-Jet 4 btagDisc", "btagDisc", "bjetBTagDisc4", {0, 1}, {false, 10, 1}}, // ADICIONADO
    {"bjetPT5", "B-Jet PT 5", "pT [GeV]", "bjetPT5", {0, 300}, {true, 10, 1}},
    {"bjetEta5", "B-Jet Eta 5", "#eta", "bjetEta5", {-3, 3}, {false, 10, 1}},         // ADICIONADO
    {"bjetBTagDisc5", "B-Jet 5 btagDisc", "btagDisc", "bjetBTagDisc5", {0, 1}, {false, 10, 1}}, // ADICIONADO
    {"bjetPT6", "B-Jet PT 6", "pT [GeV]", "bjetPT6", {0, 200}, {true, 10, 1}},
    {"bjetEta6", "B-Jet Eta 6", "#eta", "bjetEta6", {-3, 3}, {false, 10, 1}},         // ADICIONADO
    {"bjetBTagDisc6", "B-Jet 6 btagDisc", "btagDisc", "bjetBTagDisc6", {0, 1}, {false, 10, 1}}, // ADICIONADO

    // --- Variáveis Globais do Evento ---
    {"met", "Missing E_T", "MET [GeV]", "met", {0, 1000}, {false, 10, 1}}, // ADICIONADO
    {"jetHT", "Jet HT", "HT [GeV]", "jetHT", {0, 3500}, {false, 10, 1}},
    {"jetBHT", "B-Jet HT", "HT [GeV]", "jetBHT", {0, 3500}, {false, 10, 1}},
    {"jetLightHT", "Light Jet HT", "HT [GeV]", "jetLightHT", {0, 1200}, {false, 10, 1}},
    {"jetHadronicHiggsHT", "Hadronic Higgs HT", "HT [GeV]", "jetHadronicHiggsHT", {0, 4000}, {false, 10, 1}},

    // --- Contagens de Objetos ---
    {"jetNumber", "Jet Number", "Number of Jets", "jetNumber", {0, 20}, {false, 10, 1}},
    {"jetBNumber", "B-Jet Number", "Number of B-Jets", "jetBNumber", {0, 15}, {false, 10, 1}},
    {"jetLightNumber", "Light Jet Number", "Number of Light Jets", "jetLightNumber", {0, 15}, {false, 10, 1}},
    {"jetHadronicHiggsNumber", "Hadronic Higgs Number", "Number of Had. Higgs", "jetHadronicHiggsNumber", {0, 10}, {false, 10, 1}},

    // --- Massas e Reconstrução ---
    {"invMass_HH1Matched", "Invariant Mass HH1 Matched", "Mass [GeV]", "invMass_HH1Matched", {0, 500}, {true, 10, 1}},
    {"invMass_HH2Matched", "Invariant Mass HH2 Matched", "Mass [GeV]", "invMass_HH2Matched", {0, 500}, {true, 10, 1}},
    {"invMass_HH1NotMatched", "Invariant Mass HH1 Not Matched", "Mass [GeV]", "invMass_HH1NotMatched", {0, 500}, {true, 10, 1}},
    {"invMass_HH2NotMatched", "Invariant Mass HH2 Not Matched", "Mass [GeV]", "invMass_HH2NotMatched", {0, 500}, {true, 10, 1}},
    {"invMass_Higgs1_mChi", "Invariant Mass H1 min(chi2)", "Mass [GeV]", "invMass_Higgs1_mChi", {0, 500}, {true, 10, 1}},
    {"invMass_Higgs2_mChi", "Invariant Mass H2 min(chi2)", "Mass [GeV]", "invMass_Higgs2_mChi", {0, 500}, {true, 10, 1}},
    {"invMass_hadW", "Invariant Mass Hadronic W", "Mass [GeV]", "invMass_hadW", {0, 400}, {true, 10, 1}},
    {"invMass_Z1", "Invariant Mass Z1", "Mass [GeV]", "invMass_Z1", {0, 400}, {true, 10, 1}},
    {"invMass_Z2", "Invariant Mass Z2", "Mass [GeV]", "invMass_Z2", {0, 400}, {true, 10, 1}},
    {"invMass_Higgs1", "Invariant Mass H1", "Mass [GeV]", "invMass_Higgs1", {0, 500}, {true, 10, 1}},
    {"invMass_Higgs2", "Invariant Mass H2", "Mass [GeV]", "invMass_Higgs2", {0, 500}, {true, 10, 1}},
    {"invMass_HiggsZ1", "Invariant Mass HZ1", "Mass [GeV]", "invMass_HiggsZ1", {0, 500}, {true, 10, 1}},
    {"invMass_HiggsZ2", "Invariant Mass HZ2", "Mass [GeV]", "invMass_HiggsZ2", {0, 500}, {true, 10, 1}},
    {"invMass_HiggsMatched", "Invariant Mass H Matched", "Mass [GeV]", "invMass_HiggsMatched", {0, 500}, {true, 10, 1}},
    {"invMass_HiggsNotMatched", "Invariant Mass H Not Matched", "Mass [GeV]", "invMass_HiggsNotMatched", {0, 500}, {true, 10, 1}},
    {"pT_Higgs1", "pT H1", "pT [GeV]", "pT_Higgs1", {0, 500}, {true, 10, 1}},
    {"pT_Higgs2", "pT H2", "pT [GeV]", "pT_Higgs2", {0, 500}, {true, 10, 1}},
    {"chi2Higgs", "chi2 HH", "#chi^{2}_{HH}", "chi2Higgs", {0, 5000}, {true, 10, 1}},
    {"chi2Z", "chi2 ZZ", "#chi^{2}_{ZZ}", "chi2Z", {0, 5000}, {true, 10, 1}},
    {"chi2HiggsZ", "chi2 HZ", "#chi^{2}_{HZ}", "chi2HiggsZ", {0, 5000}, {true, 10, 1}},
    {"chi2HiggsNotMatched", "chi2 H Not Matched", "#chi^{2}_{H,unmatched}", "chi2HiggsNotMatched", {0, 10}, {true, 10, 1}},
    {"chi2HiggsMatched", "chi2 H Matched", "#chi^{2}_{H,matched}", "chi2HiggsMatched", {0, 10}, {true, 10, 1}},
    {"chi2HHNotMatched", "chi2 HH Not Matched", "#chi^{2}_{HH,unmatched}", "chi2HHNotMatched", {0, 10}, {true, 10, 1}},
    {"chi2HHMatched", "chi2 HH Matched", "#chi^{2}_{HH,matched}", "chi2HHMatched", {0, 10}, {true, 10, 1}},
    {"jetAvgMass", "Average Jet Mass", "Mass [GeV]", "jetAvgMass", {0, 60}, {true, 10, 1}},
    {"jetBAvgMass", "Average B-Jet Mass", "Mass [GeV]", "jetBAvgMass", {0, 250}, {true, 10, 1}},
    {"higgsHadAvgMass", "Average Hadronic Higgs Mass", "Mass [GeV]", "higgsHadAvgMass", {0, 60}, {true, 10, 1}},
    {"jetLightAvgMass", "Average Light Jet Mass", "Mass [GeV]", "jetLightAvgMass", {0, 60}, {true, 10, 1}},
    {"jetBAvgMassSqr", "Average B-Jet Mass Sqr", "Mass^2 [GeV^2]", "jetBAvgMassSqr", {0, 2500}, {true, 10, 1}},
    {"higgsHadSoftDropMass1", "Hadronic Higgs SoftDrop Mass 1", "Mass [GeV]", "higgsHadSoftDropMass1", {0, 400}, {true, 10, 1}},
    {"higgsHadSoftDropMass2", "Hadronic Higgs SoftDrop Mass 2", "Mass [GeV]", "higgsHadSoftDropMass2", {0, 300}, {true, 10, 1}},
    {"maxPTmassjbb", "Mass jbb (max pT)", "Mass [GeV]", "maxPTmassjbb", {0, 100}, {true, 10, 1}},
    {"maxPTmassjjj", "Mass jjj (max pT)", "Mass [GeV]", "maxPTmassjjj", {0, 100}, {true, 10, 1}},

    // --- Variáveis Angulares ---
    {"deltaRavgjj", "Avg Delta R jj", "#DeltaR_{jj}^{avg}", "deltaRavgjj", {0, 5}, {true, 10, 1}},
    {"deltaRavgbb", "Avg Delta R bb", "#DeltaR_{bb}^{avg}", "deltaRavgbb", {0, 5}, {true, 10, 1}},
    {"deltaRavgbj", "Avg Delta R bj", "#DeltaR_{bj}^{avg}", "deltaRavgbj", {0, 5}, {true, 10, 1}},
    {"deltaEtaavgjj", "Avg Delta Eta jj", "#Delta#eta_{jj}^{avg}", "deltaEtaavgjj", {0, 3}, {true, 10, 1}},
    {"deltaEtaavgbb", "Avg Delta Eta bb", "#Delta#eta_{bb}^{avg}", "deltaEtaavgbb", {0, 3}, {true, 10, 1}},
    {"deltaEtaavgbj", "Avg Delta Eta bj", "#Delta#eta_{bj}^{avg}", "deltaEtaavgbj", {0, 3}, {true, 10, 1}},
    {"deltaRminjj", "Min Delta R jj", "#DeltaR_{jj}^{min}", "deltaRminjj", {0, 3}, {true, 10, 1}},
    {"deltaRminbb", "Min Delta R bb", "#DeltaR_{bb}^{min}", "deltaRminbb", {0, 3}, {true, 10, 1}},
    {"deltaRminbj", "Min Delta R bj", "#DeltaR_{bj}^{min}", "deltaRminbj", {0, 3}, {true, 10, 1}},
    {"pTdeltaRminjj", "pT Min Delta R jj", "pT [GeV]", "pTdeltaRminjj", {0, 600}, {true, 10, 1}},
    {"pTdeltaRminbb", "pT Min Delta R bb", "pT [GeV]", "pTdeltaRminbb", {0, 600}, {true, 10, 1}},
    {"pTdeltaRminbj", "pT Min Delta R bj", "pT [GeV]", "pTdeltaRminbj", {0, 600}, {true, 10, 1}},
    {"massDeltaRminjj", "Mass Min Delta R jj", "Mass [GeV]", "massDeltaRminjj", {0, 100}, {true, 10, 1}},
    {"massDeltaRminbb", "Mass Min Delta R bb", "Mass [GeV]", "massDeltaRminbb", {0, 100}, {true, 10, 1}},
    {"massDeltaRminbj", "Mass Min Delta R bj", "Mass [GeV]", "massDeltaRminbj", {0, 100}, {true, 10, 1}},
    {"deltaEtamaxbb", "Max Delta Eta bb", "#Delta#eta_{bb}^{max}", "deltaEtamaxbb", {0, 5}, {true, 10, 1}},
    {"deltaEtamaxjj", "Max Delta Eta jj", "#Delta#eta_{jj}^{max}", "deltaEtamaxjj", {0, 5}, {true, 10, 1}},
    {"deltaEtamaxbj", "Max Delta Eta bj", "#Delta#eta_{bj}^{max}", "deltaEtamaxbj", {0, 5}, {true, 10, 1}},

    // --- Formato do Evento (Event Shape) ---
    {"aplanarity", "Aplanarity", "A", "aplanarity", {0, 0.5}, {true, 10, 1}},
    {"sphericity", "Sphericity", "S", "sphericity", {0, 1}, {true, 10, 1}},
    {"transSphericity", "Transverse Sphericity", "S_{#perp}", "transSphericity", {0, 1}, {true, 10, 1}},
    {"C", "C value", "C value", "C", {0, 1}, {true, 10, 1}},
    {"D", "D value", "D value", "D", {0, 1}, {true, 10, 1}},
    {"centralityjb", "Centrality jb", "centrality_{jb}", "centralityjb", {0, 1}, {true, 10, 1}},
    {"centralityjl", "Centrality jl", "centrality_{jl}", "centralityjl", {0, 1}, {true, 10, 1}},
    {"H0", "H0", "H_{0}", "H0", {0.2, 0.45}, {true, 10, 1}},
    {"H1", "H1", "H_{1}", "H1", {-0.2, 0.45}, {true, 10, 1}},
    {"H2", "H2", "H_{2}", "H2", {-0.2, 0.3}, {true, 10, 1}},
    {"H3", "H3", "H_{3}", "H3", {-0.2, 0.3}, {true, 10, 1}},
    {"H4", "H4", "H_{4}", "H4", {-0.2, 0.3}, {true, 10, 1}},
    {"R1", "R1", "R_{1}", "R1", {0, 1}, {true, 10, 1}},
    {"R2", "R2", "R_{2}", "R2", {0, 1}, {true, 10, 1}},
    {"R3", "R3", "R_{3}", "R3", {0, 1}, {true, 10, 1}},
    {"R4", "R4", "R_{4}", "R4", {0, 1}, {true, 10, 1}},
    {"H0_bjet", "H0 bjet", "H_{0,bjet}", "H0_bjet", {-0.2, 0.45}, {true, 10, 1}},
    {"H1_bjet", "H1 bjet", "H_{1,bjet}", "H1_bjet", {-0.2, 0.45}, {true, 10, 1}},
    {"H2_bjet", "H2 bjet", "H_{2,bjet}", "H2_bjet", {-0.2, 0.3}, {true, 10, 1}},
    {"H3_bjet", "H3 bjet", "H_{3,bjet}", "H3_bjet", {-0.2, 0.3}, {true, 10, 1}},
    {"H4_bjet", "H4 bjet", "H_{4,bjet}", "H4_bjet", {-0.2, 0.3}, {true, 10, 1}},
    {"R1_bjet", "R1 bjet", "R_{1,bjet}", "R1_bjet", {0, 1}, {true, 10, 1}},
    {"R2_bjet", "R2 bjet", "R_{2,bjet}", "R2_bjet", {0, 1}, {true, 10, 1}},
    {"R3_bjet", "R3 bjet", "R_{3,bjet}", "R3_bjet", {0, 1}, {true, 10, 1}},
    {"R4_bjet", "R4 bjet", "R_{4,bjet}", "R4_bjet", {0, 1}, {true, 10, 1}},
    {"aplanarity_bjet", "Aplanarity bjet", "A_{bjet}", "aplanarity_bjet", {0, 0.5}, {true, 10, 1}},
    {"sphericity_bjet", "Sphericity bjet", "S_{bjet}", "sphericity_bjet", {0, 1}, {true, 10, 1}},
    {"transSphericity_bjet", "Transverse Sphericity bjet", "S_{#perp, bjet}", "transSphericity_bjet", {0, 1}, {true, 10, 1}},
    {"C_bjet", "C value bjet", "C value_{bjet}", "C_bjet", {0, 1}, {true, 10, 1}},
    {"D_bjet", "D value bjet", "D value_{bjet}", "D_bjet", {0, 1}, {true, 10, 1}}
}}
    };
    

    Iteration_Directories_And_Histograms(processes, histogramSettings);
    return 0;
}
//to run use: cmssw-el7  cmsenv  root -l -b -q Plotter_all_norm.cpp
// to run use root -l -b -q Plotter_all_norm.cpp
