# Implementation of the plotting step of the analysis
#
# The plotting combines the histograms to plots which allow us to study the
# inital dataset based on observables motivated through physics.


import ROOT
ROOT.gROOT.SetBatch(True)


# Declare a human-readable label for each variable on the plot axis
variable_labels = {
        "Higgs_mass": "Mass 4 leptons / GeV",
        "Z1_mass": "Mass Z_{1} / GeV",
        "Z2_mass": "Mass Z_{2} / GeV",
        }

# Retrieve a histogram from the input file based on the sample, the fina state
# and the variable name
def getHistogram(tfile, sample, final_state, variable):
    name = "{}_{}_{}".format(sample, final_state, variable)
    h = tfile.Get(name)
    if not h:
        raise Exception("Failed to load histogram {}.".format(name))
    return h


# Main function of the plotting step

# The plotting takes for each variable the histograms for each final state and sample
# and combines all signal and background processes, respectively. Then, the histgrams
# are plotted in a stacked manner overlain by the data points.
# This procedure is repeated with all final states combined to make the Higgs peak visible.
def main(variable):
    tfile = ROOT.TFile("histograms.root", "READ")

    # Styles
    ROOT.gStyle.SetOptStat(0)

    ROOT.gStyle.SetCanvasBorderMode(0)
    ROOT.gStyle.SetCanvasColor(ROOT.kWhite)
    ROOT.gStyle.SetCanvasDefH(600)
    ROOT.gStyle.SetCanvasDefW(600)
    ROOT.gStyle.SetCanvasDefX(0)
    ROOT.gStyle.SetCanvasDefY(0)

    ROOT.gStyle.SetPadTopMargin(0.08)
    ROOT.gStyle.SetPadBottomMargin(0.13)
    ROOT.gStyle.SetPadLeftMargin(0.16)
    ROOT.gStyle.SetPadRightMargin(0.05)

    ROOT.gStyle.SetHistLineColor(1)
    ROOT.gStyle.SetHistLineStyle(0)
    ROOT.gStyle.SetHistLineWidth(1)
    ROOT.gStyle.SetEndErrorSize(2)
    ROOT.gStyle.SetMarkerStyle(20)

    ROOT.gStyle.SetOptTitle(0)
    ROOT.gStyle.SetTitleFont(42)
    ROOT.gStyle.SetTitleColor(1)
    ROOT.gStyle.SetTitleTextColor(1)
    ROOT.gStyle.SetTitleFillColor(10)
    ROOT.gStyle.SetTitleFontSize(0.05)

    ROOT.gStyle.SetTitleColor(1, "XYZ")
    ROOT.gStyle.SetTitleFont(42, "XYZ")
    ROOT.gStyle.SetTitleSize(0.05, "XYZ")
    ROOT.gStyle.SetTitleXOffset(1.00)
    ROOT.gStyle.SetTitleYOffset(1.60)

    ROOT.gStyle.SetLabelColor(1, "XYZ")
    ROOT.gStyle.SetLabelFont(42, "XYZ")
    ROOT.gStyle.SetLabelOffset(0.007, "XYZ")
    ROOT.gStyle.SetLabelSize(0.04, "XYZ")

    ROOT.gStyle.SetAxisColor(1, "XYZ")
    ROOT.gStyle.SetStripDecimals(True)
    ROOT.gStyle.SetTickLength(0.03, "XYZ")
    ROOT.gStyle.SetNdivisions(510, "XYZ")
    ROOT.gStyle.SetPadTickX(1)
    ROOT.gStyle.SetPadTickY(1)

    ROOT.gStyle.SetPaperSize(20., 20.)
    ROOT.gStyle.SetHatchesLineWidth(5)
    ROOT.gStyle.SetHatchesSpacing(0.05)

    ROOT.TGaxis.SetExponentOffset(-0.08, 0.01, "Y")

    # Simulation
    signals = {}
    for final_state in ["FourMuons", "FourElectrons", "TwoMuonsTwoElectrons"]:
        signals[final_state] = getHistogram(tfile, "SMHiggsToZZTo4L", final_state, variable)

    def combineFinalStates(d):
        d["combined"] = d["FourMuons"].Clone()
        d["combined"].Add(d["FourElectrons"])
        d["combined"].Add(d["TwoMuonsTwoElectrons"])

    combineFinalStates(signals)

    backgrounds = {}
    backgrounds["FourMuons"] = getHistogram(tfile, "ZZTo4mu", "FourMuons", variable)
    backgrounds["FourElectrons"] = getHistogram(tfile, "ZZTo4e", "FourElectrons", variable)
    backgrounds["TwoMuonsTwoElectrons"] = getHistogram(tfile, "ZZTo2e2mu", "TwoMuonsTwoElectrons", variable)

    combineFinalStates(backgrounds)

    # Data
    data = {}
    for final_state, samples in [
                ["FourMuons", ["dataRunBMu", "dataRunCMu"]],
                ["FourElectrons", ["dataRunBElec", "dataRunCElec"]],
                ["TwoMuonsTwoElectrons", ["dataRunBMu", "dataRunCMu",
                                          "dataRunBElec", "dataRunCElec"]]
            ]:
        for sample in samples:
            h = getHistogram(tfile, sample, final_state, variable)
            if not final_state in data:
                data[final_state] = h
            else:
                data[final_state].Add(h)

    combineFinalStates(data)

    # Draw histograms

    for final_state in ["FourMuons", "FourElectrons", "TwoMuonsTwoElectrons", "combined"]:
        data_ = data[final_state]
        data_.SetMarkerStyle(20)
        data_.SetLineColor(ROOT.kBlack)

        background = backgrounds[final_state]
        background.SetLineWidth(3)
        background.SetFillColor(ROOT.TColor.GetColor(100, 192, 232))
        background.SetLineColor(ROOT.TColor.GetColor(100, 192, 232))

        signal = signals[final_state]
        signal.Add(background)
        signal.SetLineColor(ROOT.kRed)
        signal.SetLineWidth(3)
        signal.SetTitleSize(0.04, "XYZ")
        signal.SetTitleOffset(1.2, "XYZ")

        c = ROOT.TCanvas("", "", 600, 600)
        name = data_.GetTitle()
        if name in variable_labels:
            title = variable_labels[name]
        else:
            title = name
        signal.GetXaxis().SetTitle(title)
        signal.GetYaxis().SetTitle("N_{Events}")
        signal.SetMaximum(max(background.GetMaximum(), data_.GetMaximum()) * 1.4)

        signal.Draw("HIST")
        background.Draw("HIST SAME")
        data_.Draw("E1P SAME")

        # Add legend
        legend = ROOT.TLegend(0.6, 0.66, 0.90, 0.86)
        legend.AddEntry(background, "Z#gamma*, ZZ", "f")
        legend.AddEntry(signal, "m_{H} = 125 GeV", "l")
        legend.AddEntry(data_, "Data", "lep")
        legend.SetBorderSize(0)
        legend.Draw()

        # Add title
        latex = ROOT.TLatex()
        latex.SetNDC()
        latex.SetTextSize(0.04)
        latex.SetTextFont(42)
        latex.DrawLatex(0.6, 0.935, "11.6 fb^{-1} (2012, 8 TeV)")
        latex.DrawLatex(0.16, 0.935, "#bf{CMS Open Data}")

        # Save
        c.SaveAs("{}_{}.pdf".format(final_state, variable))
        c.SaveAs("{}_{}.png".format(final_state, variable))


# Loop over all variable names and make a plot for each
if __name__ == "__main__":
    for variable in variable_labels.keys():
        main(variable)

