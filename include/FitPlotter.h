/**
 * FitPlotter.h
 *
 * Publication-quality fit plots for the phi SDME closure test and real-data fit.
 *
 * Produces three sets of plots in <outputDir>/plots/:
 *
 *  1.  xsec_Q2_<i>.pdf   — dσ/dt vs t'  per Q² bin
 *                           Upper pad: data (black points) + MC fit (red line)
 *                           Lower pad: pull = (data − MC) / σ
 *
 *  2.  bsa_Q2_<i>.pdf    — A_LU vs t'   per Q² bin   (same layout)
 *
 *  3.  param_summary.pdf — fitted value ± σ vs truth for all 8 parameters
 *                           Colour-coded by |pull|: green < 1σ, yellow 1–2σ, red > 2σ
 *
 * Usage:
 *   FitPlotter plotter(bins, outputDir);
 *   plotter.plotAll(data, mc_best, result, truthParams);   // closure test
 *   plotter.plotAll(data, mc_best, result);                // real data (no truth)
 */

#pragma once

// ROOT
#include "TCanvas.h"
#include "TPad.h"
#include "TGraphErrors.h"
#include "TGraph.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TLine.h"
#include "TBox.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TAxis.h"
#include "TColor.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TError.h"    // gErrorIgnoreLevel, kWarning

// project
#include "BinDef.h"
#include "DataSet.h"
#include "AnalysisModule.h"
#include "Chi2Fitter.h"
#include "PhysicsParams.h"

#include <vector>
#include <string>
#include <cstdio>
#include <cmath>
#include <algorithm>
#include <filesystem>

namespace fs = std::filesystem;

// ── Style constants ──────────────────────────────────────────────────────────
namespace PlotStyle {
    static const int    kData     = kBlack;
    static const int    kMC       = kRed+1;
    static const int    kBand1s   = kGreen+2;    // ±1σ pull band
    static const int    kBand2s   = kYellow+1;   // ±2σ pull band

    static const double kMarkerSize = 0.9;
    static const double kLineWidth  = 2.0;
    static const double kPadSplit   = 0.30;     // lower pad fraction

    inline void apply() {
        gStyle->SetOptStat(0);
        gStyle->SetOptTitle(0);
        gStyle->SetPadTickX(1);
        gStyle->SetPadTickY(1);
        gStyle->SetPadLeftMargin(0.14);
        gStyle->SetPadRightMargin(0.05);
        gStyle->SetPadTopMargin(0.08);
        gStyle->SetPadBottomMargin(0.12);
        gStyle->SetLegendBorderSize(0);
        gStyle->SetLegendFillColor(0);
        gStyle->SetFrameLineWidth(1);
    }
}

// ── FitPlotter ────────────────────────────────────────────────────────────────
class FitPlotter {
public:

    FitPlotter(const BinDef& bins, const std::string& outputDir)
        : bins_(bins), outDir_(outputDir + "/plots")
    {
        fs::create_directories(outDir_);
        PlotStyle::apply();
        // Suppress ROOT info messages
        gROOT->SetBatch(true);
        gErrorIgnoreLevel = kWarning;
    }

    // ── Main entry point ─────────────────────────────────────────────────────
    // hasTruth = true for closure test (draws truth comparison in param plot)
    void plotAll(const DataSet&               data,
                 const std::vector<BinResult>& mc,
                 const FitResult&              result,
                 const PhysicsParams*          truthParams = nullptr)
    {
        printf("[FitPlotter] Writing plots to %s/\n", outDir_.c_str());

        for (int iQ2 = 0; iQ2 < bins_.nQ2(); ++iQ2) {
            plotXsec(iQ2, data, mc);
            plotBSA (iQ2, data, mc);
        }
        plotParamSummary(result, truthParams);

        printf("[FitPlotter] Done — %d canvases written.\n",
               2 * bins_.nQ2() + 1);
    }

private:

    // ═══════════════════════════════════════════════════════════════════════
    // 1. dσ/dt  vs  t'  per Q² bin  (upper: data+MC,  lower: pull)
    // ═══════════════════════════════════════════════════════════════════════
    void plotXsec(int iQ2,
                  const DataSet&               data,
                  const std::vector<BinResult>& mc) const
    {
        // Collect matching points
        std::vector<double> t, sigma, sig_err, mc_sigma, pull;
        for (const auto& dp : data.xsec) {
            if (dp.iQ2 != iQ2) continue;
            int idx = bins_.flatIndex(dp.iQ2, dp.it, dp.ixB);
            if (idx < 0 || idx >= (int)mc.size()) continue;
            const auto& r = mc[idx];
            if (r.dsigma_dt <= 0 || dp.error <= 0) continue;
            t       .push_back(bins_.tbins.center(dp.it));
            sigma   .push_back(dp.value);
            sig_err .push_back(dp.error);
            mc_sigma.push_back(r.dsigma_dt);
            pull    .push_back((dp.value - r.dsigma_dt) / dp.error);
        }
        if (t.empty()) return;

        double tMin = bins_.tbins.lo(0),  tMax = bins_.tbins.hi(bins_.nT()-1);
        double q2lo = bins_.Q2bins.lo(iQ2), q2hi = bins_.Q2bins.hi(iQ2);

        char name[128], title[256];
        snprintf(name,  sizeof(name),  "c_xsec_Q2_%d", iQ2);
        snprintf(title, sizeof(title), "d#sigma/dt  |  Q^{2} #in [%.2f, %.2f] GeV^{2}", q2lo, q2hi);

        TCanvas* c = makeCanvas(name, 600, 650);
        auto [padUp, padDn] = makeSplitPads(c);

        // ── Upper pad ─────────────────────────────────────────────────────
        padUp->cd();
        padUp->SetLogy();

        int n = (int)t.size();
        std::vector<double> tErr(n, 0.0);
        TGraphErrors* gData = new TGraphErrors(n,
            t.data(), sigma.data(), tErr.data(), sig_err.data());
        TGraph*       gMC   = new TGraph(n, t.data(), mc_sigma.data());

        styleData(gData);
        styleMC(gMC);

        // Frame
        double yMin = *std::min_element(sigma.begin(), sigma.end()) * 0.3;
        double yMax = *std::max_element(sigma.begin(), sigma.end()) * 4.0;
        TH1F* frame = padUp->DrawFrame(tMin, yMin, tMax, yMax);
        styleFrame(frame, "", "d#sigma/dt [nb/GeV^{2}]", true);

        gMC->Draw("L SAME");
        gData->Draw("P SAME");

        TLegend* leg = new TLegend(0.55, 0.72, 0.90, 0.88);
        leg->AddEntry(gData, "Pseudo-data", "pe");
        leg->AddEntry(gMC,   "Best fit",    "l");
        leg->SetTextSize(0.045);
        leg->Draw();

        drawLabel(title, 0.18, 0.88);
        drawCLAS12Label(0.18, 0.80);

        // ── Lower pad: pull ────────────────────────────────────────────────
        padDn->cd();
        drawPullPanel(padDn, n, t.data(), pull.data(), tErr.data(),
                      tMin, tMax, "t' [GeV^{2}]");

        std::string fname = outDir_ + "/xsec_Q2_" + std::to_string(iQ2) + ".pdf";
        c->SaveAs(fname.c_str());
        printf("  → %s\n", fname.c_str());
        delete c;
    }

    // ═══════════════════════════════════════════════════════════════════════
    // 2. A_LU  vs  t'  per Q² bin
    // ═══════════════════════════════════════════════════════════════════════
    void plotBSA(int iQ2,
                 const DataSet&               data,
                 const std::vector<BinResult>& mc) const
    {
        std::vector<double> t, bsa, bsa_err, mc_bsa, pull;
        for (const auto& dp : data.bsa) {
            if (dp.iQ2 != iQ2) continue;
            int idx = bins_.flatIndex(dp.iQ2, dp.it, dp.ixB);
            if (idx < 0 || idx >= (int)mc.size()) continue;
            const auto& r = mc[idx];
            if (dp.error <= 0) continue;
            t      .push_back(bins_.tbins.center(dp.it));
            bsa    .push_back(dp.value);
            bsa_err.push_back(dp.error);
            mc_bsa .push_back(r.A_LU);
            pull   .push_back((dp.value - r.A_LU) / dp.error);
        }
        if (t.empty()) return;

        double tMin = bins_.tbins.lo(0),  tMax = bins_.tbins.hi(bins_.nT()-1);
        double q2lo = bins_.Q2bins.lo(iQ2), q2hi = bins_.Q2bins.hi(iQ2);

        char name[128], title[256];
        snprintf(name,  sizeof(name),  "c_bsa_Q2_%d", iQ2);
        snprintf(title, sizeof(title), "A_{LU}  |  Q^{2} #in [%.2f, %.2f] GeV^{2}", q2lo, q2hi);

        TCanvas* c = makeCanvas(name, 600, 650);
        auto [padUp, padDn] = makeSplitPads(c);

        padUp->cd();
        int n = (int)t.size();
        std::vector<double> tErr(n, 0.0);
        TGraphErrors* gData = new TGraphErrors(n,
            t.data(), bsa.data(), tErr.data(), bsa_err.data());
        TGraph* gMC = new TGraph(n, t.data(), mc_bsa.data());

        styleData(gData);
        styleMC(gMC);

        // Symmetric y range around 0 — use index loop, NOT range-for with pointer arithmetic
        double yabs = 0.0;
        for (int i = 0; i < (int)bsa.size(); ++i)
            yabs = std::max(yabs, std::abs(bsa[i]) + 3.0 * bsa_err[i]);
        for (double v : mc_bsa)
            yabs = std::max(yabs, std::abs(v) * 1.3);
        yabs = std::max(yabs, 0.15);
        yabs = std::min(yabs, 0.70);

        TH1F* frame = padUp->DrawFrame(tMin, -yabs, tMax, yabs);
        styleFrame(frame, "", "A_{LU}", false);

        // Zero line
        TLine* zero = new TLine(tMin, 0, tMax, 0);
        zero->SetLineStyle(2); zero->SetLineColor(kGray+1); zero->Draw();

        gMC->Draw("L SAME");
        gData->Draw("P SAME");

        TLegend* leg = new TLegend(0.55, 0.72, 0.90, 0.88);
        leg->AddEntry(gData, "Pseudo-data", "pe");
        leg->AddEntry(gMC,   "Best fit",    "l");
        leg->SetTextSize(0.045);
        leg->Draw();

        drawLabel(title, 0.18, 0.88);
        drawCLAS12Label(0.18, 0.80);

        padDn->cd();
        drawPullPanel(padDn, n, t.data(), pull.data(), tErr.data(),
                      tMin, tMax, "t' [GeV^{2}]");

        std::string fname = outDir_ + "/bsa_Q2_" + std::to_string(iQ2) + ".pdf";
        c->SaveAs(fname.c_str());
        printf("  → %s\n", fname.c_str());
        delete c;
    }

    // ═══════════════════════════════════════════════════════════════════════
    // 3. Parameter summary — fitted ± σ vs truth, pull colour-coded
    // ═══════════════════════════════════════════════════════════════════════
    void plotParamSummary(const FitResult&     result,
                          const PhysicsParams* truth) const
    {
        const bool hasTruth = (truth != nullptr);

        // Parameter rows in display order
        struct Row {
            std::string label;   // LaTeX label for axis
            double fit, err, truthVal;
        };

        // Map result arrays onto named rows
        // errors[] order: A_L(0) b_L(1) Im_LL(2) A_T(3) b_T(4) Im_TT(5) A_LT(6) A_LTi(7)
        auto fe = [&](int i) { return i < (int)result.errors.size() ? result.errors[i] : 0.0; };
        auto fv = [&](int i) { return i < (int)result.bestValues.size() ? result.bestValues[i] : 0.0; };

        std::vector<Row> rows = {
            {"A_{L}",    fv(0), fe(0), hasTruth ? truth->A_L   : fv(0)},
            {"b_{L}",    fv(1), fe(1), hasTruth ? truth->b_L   : fv(1)},
            {"Im_{LL}",  fv(2), fe(2), hasTruth ? truth->Im_LL : fv(2)},
            {"A_{T}",    fv(3), fe(3), hasTruth ? truth->A_T   : fv(3)},
            {"b_{T}",    fv(4), fe(4), hasTruth ? truth->b_T   : fv(4)},
            {"Im_{TT}",  fv(5), fe(5), hasTruth ? truth->Im_TT : fv(5)},
            {"A_{LT}",   fv(6), fe(6), hasTruth ? truth->A_LT  : fv(6)},
            {"A_{LTi}",  fv(7), fe(7), hasTruth ? truth->A_LTi : fv(7)},
        };
        int nPar = (int)rows.size();

        // ── Canvas layout: left = value comparison, right = pull bar ──────
        TCanvas* c = makeCanvas("c_param_summary", 1100, 600);
        c->Divide(2, 1, 0.005, 0.01);

        // ── Left pad: fitted value ± σ  (+ truth marker if available) ─────
        {
            c->cd(1);
            gPad->SetLeftMargin(0.22);
            gPad->SetRightMargin(0.04);
            gPad->SetTopMargin(0.10);
            gPad->SetBottomMargin(0.12);

            // Determine x range from all values + errors + truth
            double xlo = 1e30, xhi = -1e30;
            for (const auto& r : rows) {
                xlo = std::min(xlo, std::min(r.fit - 3*r.err, r.truthVal));
                xhi = std::max(xhi, std::max(r.fit + 3*r.err, r.truthVal));
            }
            double margin = (xhi - xlo) * 0.20;
            xlo -= margin; xhi += margin;

            // Frame using TH2F so we can set y-axis bin labels
            TH2F* hFrame = new TH2F("hParFrame", "",
                                    1, xlo, xhi,
                                    nPar, 0, nPar);
            hFrame->GetXaxis()->SetTitle("Parameter value");
            hFrame->GetXaxis()->SetTitleSize(0.050);
            hFrame->GetXaxis()->SetLabelSize(0.042);
            hFrame->GetYaxis()->SetLabelSize(0.055);
            hFrame->GetYaxis()->SetTickLength(0.0);
            for (int i = 0; i < nPar; ++i)
                hFrame->GetYaxis()->SetBinLabel(i+1, rows[nPar-1-i].label.c_str());
            hFrame->Draw("AXIS");

            // 1σ and 2σ bands around each truth (or fit if no truth)
            for (int i = 0; i < nPar; ++i) {
                int row = nPar - 1 - i;   // reversed for visual order
                double yc = i + 0.5;
                double ref = hasTruth ? rows[row].truthVal : rows[row].fit;
                double sig = rows[row].err;

                // ±2σ band
                TBox* b2 = new TBox(ref-2*sig, yc-0.38, ref+2*sig, yc+0.38);
                b2->SetFillColorAlpha(kYellow-9, 0.6); b2->SetLineWidth(0); b2->Draw();
                // ±1σ band
                TBox* b1 = new TBox(ref-sig, yc-0.38, ref+sig, yc+0.38);
                b1->SetFillColorAlpha(kGreen-9, 0.6); b1->SetLineWidth(0); b1->Draw();
            }

            // Vertical line at zero or truth per-parameter
            for (int i = 0; i < nPar; ++i) {
                int row = nPar - 1 - i;
                double ref = hasTruth ? rows[row].truthVal : 0.0;
                double yc = i + 0.5;
                TLine* vl = new TLine(ref, yc-0.45, ref, yc+0.45);
                vl->SetLineColor(hasTruth ? (kBlue+1) : kGray+1);
                vl->SetLineWidth(hasTruth ? 2 : 1);
                vl->SetLineStyle(hasTruth ? 1 : 2);
                vl->Draw();
            }

            // Fitted values as horizontal error bars
            for (int i = 0; i < nPar; ++i) {
                int row = nPar - 1 - i;
                double yc  = i + 0.5;
                double fit = rows[row].fit;
                double err = rows[row].err;

                // Colour by pull
                double pull = hasTruth && err > 0 ?
                    std::abs(fit - rows[row].truthVal) / err : 0.0;
                int col = pullColour(pull);

                // Central point
                TGraph* gpt = new TGraph(1, &fit, &yc);
                gpt->SetMarkerStyle(kFullCircle);
                gpt->SetMarkerSize(1.2);
                gpt->SetMarkerColor(col);
                gpt->Draw("P SAME");

                // Error bar
                TLine* el = new TLine(fit-err, yc, fit+err, yc);
                el->SetLineColor(col); el->SetLineWidth(2); el->Draw();

                // End caps
                double capH = 0.15;
                TLine* lc = new TLine(fit-err, yc-capH, fit-err, yc+capH);
                TLine* rc = new TLine(fit+err, yc-capH, fit+err, yc+capH);
                lc->SetLineColor(col); lc->SetLineWidth(2); lc->Draw();
                rc->SetLineColor(col); rc->SetLineWidth(2); rc->Draw();

                // Numerical label to the right
                char buf[64];
                snprintf(buf, sizeof(buf), "%.3g #pm %.3g", fit, err);
                TLatex* lx = new TLatex(xhi - (xhi-xlo)*0.02, yc, buf);
                lx->SetNDC(false);
                lx->SetTextSize(0.036);
                lx->SetTextAlign(31);
                lx->SetTextColor(col);
                lx->Draw();
            }

            // Legend
            TLegend* leg = new TLegend(0.25, 0.01, 0.75, 0.10);
            leg->SetNColumns(3);
            leg->SetTextSize(0.038);
            TBox* lb1 = new TBox(0,0,1,1); lb1->SetFillColorAlpha(kGreen-9,0.6);
            TBox* lb2 = new TBox(0,0,1,1); lb2->SetFillColorAlpha(kYellow-9,0.6);
            leg->AddEntry(lb1, "#pm1#sigma", "f");
            leg->AddEntry(lb2, "#pm2#sigma", "f");
            if (hasTruth) {
                TLine* tl = new TLine(0,0,1,1);
                tl->SetLineColor(kBlue+1); tl->SetLineWidth(2);
                leg->AddEntry(tl, "Truth", "l");
            }
            leg->Draw();

            drawLabel("Fit parameters", 0.25, 0.93);
        }

        // ── Right pad: pull bar chart ──────────────────────────────────────
        {
            c->cd(2);
            gPad->SetLeftMargin(0.18);
            gPad->SetRightMargin(0.08);
            gPad->SetTopMargin(0.10);
            gPad->SetBottomMargin(0.12);

            double pullMax = 3.5;
            TH2F* hPull = new TH2F("hPullFrame", "",
                                   1, -pullMax, pullMax,
                                   nPar, 0, nPar);
            hPull->GetXaxis()->SetTitle("Pull = (fit #minus truth) / #sigma_{fit}");
            hPull->GetXaxis()->SetTitleSize(0.050);
            hPull->GetXaxis()->SetLabelSize(0.042);
            hPull->GetYaxis()->SetLabelSize(0.055);
            hPull->GetYaxis()->SetTickLength(0.0);
            for (int i = 0; i < nPar; ++i)
                hPull->GetYaxis()->SetBinLabel(i+1, rows[nPar-1-i].label.c_str());
            hPull->Draw("AXIS");

            // ±1σ and ±2σ bands
            TBox* band2 = new TBox(-2, 0, 2, nPar);
            band2->SetFillColorAlpha(kYellow-9, 0.5); band2->SetLineWidth(0); band2->Draw();
            TBox* band1 = new TBox(-1, 0, 1, nPar);
            band1->SetFillColorAlpha(kGreen-9, 0.5); band1->SetLineWidth(0); band1->Draw();

            // Zero line
            TLine* zero = new TLine(0, 0, 0, nPar);
            zero->SetLineStyle(2); zero->SetLineColor(kGray+1); zero->Draw();

            if (hasTruth) {
                for (int i = 0; i < nPar; ++i) {
                    int row = nPar - 1 - i;
                    double pull = (rows[row].err > 0) ?
                        (rows[row].fit - rows[row].truthVal) / rows[row].err : 0.0;
                    double yc = i + 0.5;

                    // Clamp display but draw arrow if outside
                    bool clipped = std::abs(pull) > pullMax * 0.95;
                    double displayPull = std::max(-pullMax*0.92, std::min(pullMax*0.92, pull));

                    int col = pullColour(std::abs(pull));

                    // Horizontal bar from 0 to pull
                    double barH = 0.30;
                    TBox* bar = new TBox(std::min(0.0, displayPull), yc - barH,
                                         std::max(0.0, displayPull), yc + barH);
                    bar->SetFillColor(col);
                    bar->SetLineColor(col);
                    bar->Draw();

                    // Numerical label
                    char buf[32];
                    if (clipped) snprintf(buf, sizeof(buf), "%.1f#rightarrow", pull);
                    else         snprintf(buf, sizeof(buf), "%+.2f", pull);

                    double labelX = pull > 0 ? displayPull + 0.08 : displayPull - 0.08;
                    int    align  = pull > 0 ? 12 : 32;

                    // Keep label inside frame
                    if (labelX > pullMax * 0.85)  { labelX = pullMax  * 0.85; align = 32; }
                    if (labelX < -pullMax * 0.85)  { labelX = -pullMax * 0.85; align = 12; }

                    TLatex* lx = new TLatex(labelX, yc, buf);
                    lx->SetNDC(false);
                    lx->SetTextSize(0.042);
                    lx->SetTextAlign(align);
                    lx->SetTextColor(col);
                    lx->Draw();
                }
            } else {
                // No truth: draw a note
                TLatex* note = new TLatex(0.5, 0.5,
                    "Pull requires truth (closure test)");
                note->SetNDC(true);
                note->SetTextAlign(22);
                note->SetTextSize(0.038);
                note->SetTextColor(kGray+2);
                note->Draw();
            }

            // Legend boxes
            TLegend* leg = new TLegend(0.20, 0.01, 0.80, 0.10);
            leg->SetNColumns(2);
            leg->SetTextSize(0.038);
            TBox* lb1 = new TBox(0,0,1,1); lb1->SetFillColorAlpha(kGreen-9, 0.5);
            TBox* lb2 = new TBox(0,0,1,1); lb2->SetFillColorAlpha(kYellow-9,0.5);
            leg->AddEntry(lb1, "|pull| < 1", "f");
            leg->AddEntry(lb2, "|pull| < 2", "f");
            leg->Draw();

            drawLabel("Parameter pulls", 0.25, 0.93);
        }

        // ── Shared title at top ────────────────────────────────────────────
        c->cd(0);
        TLatex* mainTitle = new TLatex(0.5, 0.97,
            hasTruth ? "#phi SDME Closure Test  #font[52]{#minus}  CLAS12 RGA  #sqrt{s} = 10.6 GeV"
                     : "#phi SDME Fit Result  #font[52]{#minus}  CLAS12 RGA  #sqrt{s} = 10.6 GeV");
        mainTitle->SetNDC(true);
        mainTitle->SetTextAlign(23);
        mainTitle->SetTextSize(0.030);
        mainTitle->SetTextColor(kBlack);
        mainTitle->Draw();

        std::string fname = outDir_ + "/param_summary.pdf";
        c->SaveAs(fname.c_str());
        printf("  → %s\n", fname.c_str());
        delete c;
    }

    // ═══════════════════════════════════════════════════════════════════════
    // Helpers
    // ═══════════════════════════════════════════════════════════════════════

    static TCanvas* makeCanvas(const char* name, int w, int h) {
        return new TCanvas(name, name, w, h);
    }

    // Returns {upperPad, lowerPad} with shared x axis
    static std::pair<TPad*, TPad*> makeSplitPads(TCanvas* c) {
        const double split = PlotStyle::kPadSplit;
        TPad* pUp = new TPad("pUp", "", 0.0, split, 1.0, 1.0);
        TPad* pDn = new TPad("pDn", "", 0.0, 0.0,   1.0, split);

        pUp->SetBottomMargin(0.015);
        pUp->SetTopMargin(0.10);
        pUp->SetLeftMargin(0.15);
        pUp->SetRightMargin(0.05);

        pDn->SetTopMargin(0.02);
        pDn->SetBottomMargin(0.35);
        pDn->SetLeftMargin(0.15);
        pDn->SetRightMargin(0.05);

        c->cd(); pUp->Draw(); pDn->Draw();
        return {pUp, pDn};
    }

    static void styleData(TGraphErrors* g) {
        g->SetMarkerStyle(kFullCircle);
        g->SetMarkerSize(PlotStyle::kMarkerSize);
        g->SetMarkerColor(PlotStyle::kData);
        g->SetLineColor(PlotStyle::kData);
        g->SetLineWidth(1);
    }

    static void styleMC(TGraph* g) {
        g->SetLineColor(PlotStyle::kMC);
        g->SetLineWidth((int)PlotStyle::kLineWidth + 1);
    }

    static void styleFrame(TH1F* h, const char* xTitle,
                           const char* yTitle, bool isLog) {
        h->GetXaxis()->SetTitle(xTitle);
        h->GetXaxis()->SetTitleSize(0.00);   // hidden — shown on lower pad
        h->GetXaxis()->SetLabelSize(0.00);
        h->GetYaxis()->SetTitle(yTitle);
        h->GetYaxis()->SetTitleSize(0.065);
        h->GetYaxis()->SetTitleOffset(0.95);
        h->GetYaxis()->SetLabelSize(0.058);
        (void)isLog;
    }

    // Draw a pull panel in the lower pad
    static void drawPullPanel(TPad* pad, int n,
                               const double* x, const double* pullVals,
                               const double* /*xErr*/,
                               double xMin, double xMax,
                               const char* xTitle)
    {
        // ── Frame MUST be drawn first so the coordinate system is set ──────
        // Any object drawn before DrawFrame uses the default [0,1]×[0,1]
        // coordinate system and becomes invalid after DrawFrame resets it.
        double pullRange = 3.5;
        for (int i = 0; i < n; ++i)
            pullRange = std::max(pullRange, std::abs(pullVals[i]) * 1.2);
        pullRange = std::min(pullRange, 5.5);

        TH1F* hf = pad->DrawFrame(xMin, -pullRange, xMax, pullRange);
        hf->GetXaxis()->SetTitle(xTitle);
        hf->GetXaxis()->SetTitleSize(0.12);
        hf->GetXaxis()->SetTitleOffset(0.95);
        hf->GetXaxis()->SetLabelSize(0.10);
        hf->GetYaxis()->SetTitle("Pull");
        hf->GetYaxis()->SetTitleSize(0.11);
        hf->GetYaxis()->SetTitleOffset(0.55);
        hf->GetYaxis()->SetLabelSize(0.09);
        hf->GetYaxis()->SetNdivisions(505);

        // ── ±2σ / ±1σ bands — drawn AFTER frame so coordinates are valid ──
        TBox* bg2 = new TBox(xMin, -2.0, xMax,  2.0);
        bg2->SetFillColorAlpha(kYellow-9, 0.5); bg2->SetLineWidth(0); bg2->Draw();
        TBox* bg1 = new TBox(xMin, -1.0, xMax,  1.0);
        bg1->SetFillColorAlpha(kGreen-9, 0.5);  bg1->SetLineWidth(0); bg1->Draw();

        // Zero line
        TLine* zero = new TLine(xMin, 0, xMax, 0);
        zero->SetLineStyle(2); zero->SetLineColor(kGray+1); zero->Draw();

        // Redraw the frame axes on top of the boxes so tick marks are visible
        hf->Draw("AXIS SAME");

        // Pull points, colour-coded per point
        for (int i = 0; i < n; ++i) {
            double xi = x[i];
            double pi = std::max(-pullRange*0.95,
                                  std::min( pullRange*0.95, pullVals[i]));
            int col = pullColour(std::abs(pullVals[i]));

            TGraph* pt = new TGraph(1, &xi, &pi);
            pt->SetMarkerStyle(kFullCircle);
            pt->SetMarkerSize(0.9);
            pt->SetMarkerColor(col);
            pt->Draw("P SAME");

            // Arrow if clipped
            if (std::abs(pullVals[i]) > pullRange * 0.90) {
                double sign = std::copysign(1.0, pullVals[i]);
                TLine* arrow = new TLine(xi, pullRange*0.75*sign,
                                          xi, pullRange*0.92*sign);
                arrow->SetLineColor(col);
                arrow->SetLineWidth(2);
                arrow->Draw();
            }
        }
    }

    // Colour based on |pull|: green (<1) → yellow (1–2) → red (>2)
    static int pullColour(double absPull) {
        if (absPull < 1.0) return kGreen+2;
        if (absPull < 2.0) return kOrange+1;
        return kRed+1;
    }

    static void drawLabel(const char* text, double x, double y) {
        TLatex* l = new TLatex(x, y, text);
        l->SetNDC(true);
        l->SetTextSize(0.052);
        l->SetTextFont(42);
        l->Draw();
    }

    static void drawCLAS12Label(double x, double y) {
        TLatex* l = new TLatex(x, y,
            "#font[62]{CLAS12}  #font[52]{simulation}  #sqrt{s} = 10.6 GeV");
        l->SetNDC(true);
        l->SetTextSize(0.042);
        l->SetTextColor(kGray+2);
        l->Draw();
    }

    const BinDef&  bins_;
    std::string    outDir_;
};
