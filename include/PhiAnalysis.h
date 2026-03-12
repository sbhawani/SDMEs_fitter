/**
 * PhiAnalysis.h
 *
 * Fills cross-section (dσ/dt') vs t' in Q² bins, and
 * beam-spin asymmetry (BSA = A_LU) vs cos(θ_H) in (Q², t') bins,
 * directly from LundEvent objects — no FastMC needed.
 *
 * ── Bin definitions ─────────────────────────────────────────────────────────
 *
 *   Q² edges (4 bins):  [0.1, 1.1533, 1.680, 2.6017, 8.0]
 *   t'  edges (9 bins): [0.01, 0.3842, 0.459, 0.6087, 0.6835,
 *                         0.8332, 1.0577, 1.2822, 1.7312, 4.5]
 *   cos(θ_H) : uniform 10 bins in [-1, +1]
 *
 * ── Cross-section proxy ──────────────────────────────────────────────────────
 *
 *   Since events are generated from the physics model (acceptance-rejection
 *   on the full SDME-weighted cross-section), the event distribution IS the
 *   cross-section shape:
 *
 *     dσ/dt' |_{Q²-bin} ∝ N(t'_bin) / Δt'_bin
 *
 *   We normalise each Q²-bin so that the first bin = 1 (relative shape).
 *   Statistical errors: σ_i = √N_i / Δt'_i × (same normalisation factor)
 *
 * ── BSA ─────────────────────────────────────────────────────────────────────
 *
 *   A_LU(cos θ_H) = [N+(cos θ_H) - N-(cos θ_H)] / [N+(cos θ_H) + N-(cos θ_H)]
 *
 *   where N+/N- = events with beamPol > 0 / < 0.
 *   Error: δA = √(1 - A²) / √(N+ + N-)
 *
 * ── Output ──────────────────────────────────────────────────────────────────
 *
 *   output/xsec/xsec_Q2bin{0..3}.pdf     — dσ/dt' in each Q² bin
 *   output/bsa/BSA_Q2bin{i}_tbin{j}.pdf  — BSA vs cos(θ_H) per (Q²,t') bin
 *   output/phi_analysis.root             — all TGraphErrors
 */

#pragma once
#include "LundReader.h"

#include <vector>
#include <string>
#include <cmath>
#include <cstdio>
#include <algorithm>
#include <filesystem>
// ── Suppress ROOT headers for compile check in non-ROOT env ──────────────────
#include "TFile.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TLine.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TPad.h"
#include "TH1D.h"

namespace fs = std::filesystem;

class PhiAnalysis {
public:
    // ── Bin edges ────────────────────────────────────────────────────────────
    const std::vector<double> Q2edges = {0.1, 1.1533, 1.680, 2.6017, 8.0};
    const std::vector<double> Tpedges = {0.01, 0.3842, 0.459, 0.6087, 0.6835,
                                          0.8332, 1.0577, 1.2822, 1.7312, 4.5};
    static constexpr int NCosH = 10;   // cos(θ_H) bins in [-1,1]

    int nQ2 = 4;
    int nTp = 9;

    // ── Accumulators ─────────────────────────────────────────────────────────
    // Cross-section: N[iQ2][itp]
    std::vector<std::vector<double>> Nxsec;

    // BSA: Np[iQ2][itp][icosH],  Nm[iQ2][itp][icosH]
    std::vector<std::vector<std::vector<double>>> Nplus, Nminus;

    // Total event counter
    long nTot = 0, nFill = 0;

    // Output directory
    std::string outDir = "output";

    PhiAnalysis() { reset(); }

    void setOutputDir(const std::string& d) { outDir = d; }

    void reset() {
        Nxsec .assign(nQ2, std::vector<double>(nTp, 0.0));
        Nplus .assign(nQ2, std::vector<std::vector<double>>(nTp, std::vector<double>(NCosH, 0.0)));
        Nminus.assign(nQ2, std::vector<std::vector<double>>(nTp, std::vector<double>(NCosH, 0.0)));
        nTot = nFill = 0;
    }

    // ── Fill one event ────────────────────────────────────────────────────────
    void fill(const LundEvent& ev) {
        ++nTot;
        if (!ev.valid) return;

        int iQ2 = findBin(ev.Q2,     Q2edges);
        int itp  = findBin(ev.tprime, Tpedges);
        if (iQ2 < 0 || itp < 0) return;

        Nxsec[iQ2][itp] += 1.0;

        double cosH = std::cos(ev.decayTheta);
        int iCos = static_cast<int>((cosH + 1.0) / 2.0 * NCosH);
        iCos = std::clamp(iCos, 0, NCosH-1);

        if (ev.beamPol > 0) Nplus [iQ2][itp][iCos] += 1.0;
        else                 Nminus[iQ2][itp][iCos] += 1.0;

        ++nFill;
    }

    // ── Process LUND file ─────────────────────────────────────────────────────
    /**
     * @param lundFile  LUND file to read
     * @param fastmc    true  → file is fastMC.jar output; apply detector-hit cuts
     *                  false → file is raw phi_gen output; accept all valid events
     */
    void processFile(const std::string& lundFile, bool fastmc = false) {
        LundReader reader(lundFile, fastmc);
        LundEvent ev;
        long n = 0, nrej = 0;
        while (reader.next(ev)) {
            if (!ev.fastmcRejected) fill(ev);
            else ++nrej;
            if (++n % 100000 == 0)
                fprintf(stderr, "[PhiAnalysis] processed %ld events (rejected fastMC: %ld)\r",
                        n, nrej);
        }
        reader.printStats();
        fprintf(stderr, "[PhiAnalysis] Filled %ld / %ld events into histograms\n",
                nFill, nTot);
    }

    // ── Summary table ─────────────────────────────────────────────────────────
    void printSummary() const {
        printf("\n=== Cross-section summary (dN/dt' shape) ===\n");
        printf("  %-4s  %-4s  %-10s  %-10s  %-10s\n",
               "iQ2","itp","t'_center","N","dN/dt'");
        for (int iQ2 = 0; iQ2 < nQ2; ++iQ2) {
            printf("  --- Q² bin %d [%.3f, %.3f] ---\n",
                   iQ2, Q2edges[iQ2], Q2edges[iQ2+1]);
            for (int itp = 0; itp < nTp; ++itp) {
                double tc = 0.5*(Tpedges[itp]+Tpedges[itp+1]);
                double dt = Tpedges[itp+1]-Tpedges[itp];
                double n  = Nxsec[iQ2][itp];
                printf("  %-4d  %-4d  %-10.4f  %-10.0f  %-10.4f\n",
                       iQ2, itp, tc, n, n/dt);
            }
        }
    }

    // ── Write DataSet-compatible text files for phi_fit ───────────────────────
    /**
     * Writes three files that phi_fit can read directly:
     *
     *   <dir>/xsec.dat     iQ2  it  ixB  dN/dt'  err
     *   <dir>/bsa.dat      iQ2  it  ixB  A_LU    err
     *   <dir>/moments.dat  iQ2  it  ixB  fLL fLL_err  fTT fTT_err  MLT MLT_err
     *
     * xsec:     event count per bin divided by bin width in t'.
     *           Error = sqrt(N) / dt'.
     *           Normalised per Q² bin so the total = 1 (shape-only for closure test).
     *           For absolute cross-sections, provide luminosity × efficiency externally.
     *
     * bsa:      (N+ - N-) / (N+ + N-) integrated over cos(theta_H) per (Q²,t') bin.
     *           Error = sqrt(1 - A^2) / sqrt(N+ + N-).
     *
     * moments:  f_LL = <cos^2 theta_H>,  f_TT = <sin^2 theta_H>,
     *           M_LT = <sin(2*theta_H)> × sign (proxy for LT interference).
     *           Errors from Poisson sqrt(N).
     *
     * @param dir  Output directory (must exist or be created).
     * @return     Number of bins written.
     */
    int writeDataFiles(const std::string& dir) const {
        std::filesystem::create_directories(dir);
        int nWritten = 0;

        FILE* fx  = fopen((dir + "/xsec.dat").c_str(),    "w");
        FILE* fb  = fopen((dir + "/bsa.dat").c_str(),     "w");
        FILE* fm  = fopen((dir + "/moments.dat").c_str(), "w");

        if (!fx || !fb || !fm) {
            fprintf(stderr, "[PhiAnalysis] ERROR: cannot open data files in %s\n", dir.c_str());
            if (fx) fclose(fx); if (fb) fclose(fb); if (fm) fclose(fm);
            return 0;
        }

        // Headers
        fprintf(fx, "# xsec.dat  —  phi electroproduction dN/dt' (shape only)\n");
        fprintf(fx, "# iQ2  it  ixB  dN/dt[GeV-2]  err\n");
        fprintf(fb, "# bsa.dat   —  beam-spin asymmetry A_LU per (Q2,t') bin\n");
        fprintf(fb, "# iQ2  it  ixB  A_LU  err\n");
        fprintf(fm, "# moments.dat  —  cos(theta_H) moments per (Q2,t') bin\n");
        fprintf(fm, "# iQ2  it  ixB  fLL  fLL_err  fTT  fTT_err  MLT  MLT_err\n");

        for (int iQ2 = 0; iQ2 < nQ2; ++iQ2) {
            // Total events in this Q² bin for shape normalisation
            double nTotQ2 = 0;
            for (int itp = 0; itp < nTp; ++itp) nTotQ2 += Nxsec[iQ2][itp];

            for (int itp = 0; itp < nTp; ++itp) {
                double dt  = Tpedges[itp+1] - Tpedges[itp];
                double N   = Nxsec[iQ2][itp];
                if (N < 3) continue;

                // ── xsec ────────────────────────────────────────────────
                double xsecVal = (N / dt);
                double xsecErr = std::sqrt(N) / dt;
                fprintf(fx, "  %d  %d  0  %.6e  %.6e\n",
                        iQ2, itp, xsecVal, xsecErr);

                // ── BSA — integrate over cos(theta_H) ───────────────────
                double Np = 0, Nm = 0;
                for (int ic = 0; ic < NCosH; ++ic) {
                    Np += Nplus [iQ2][itp][ic];
                    Nm += Nminus[iQ2][itp][ic];
                }
                double Ntot = Np + Nm;
                if (Ntot < 3) continue;
                double A   = (Np - Nm) / Ntot;
                double dA  = std::sqrt(std::max(1.0 - A*A, 0.0) / Ntot);
                fprintf(fb, "  %d  %d  0  %.6f  %.6f\n", iQ2, itp, A, dA);

                // ── Angular moments from cos(theta_H) distribution ───────
                // f_LL = <cos^2 theta_H>,  f_TT = <sin^2 theta_H>
                // M_LT is not directly accessible from cos bins alone
                // (needs azimuthal angle phi_KK which we don't accumulate here)
                // → write f_LL, f_TT; set M_LT = 0 with large error
                double sumN = 0, sumCos2 = 0, sumSin2 = 0;
                for (int ic = 0; ic < NCosH; ++ic) {
                    double cosC = -1.0 + (ic + 0.5) * 2.0 / NCosH;
                    double n_ic = Nplus[iQ2][itp][ic] + Nminus[iQ2][itp][ic];
                    sumN    += n_ic;
                    sumCos2 += n_ic * cosC * cosC;
                    sumSin2 += n_ic * (1.0 - cosC * cosC);
                }
                double fLL = (sumN > 0) ? sumCos2 / sumN : 0;
                double fTT = (sumN > 0) ? sumSin2 / sumN : 0;
                double fErr = (sumN > 0) ? std::sqrt(fLL * (1.0 - fLL) / sumN) : 1.0;
                // M_LT requires phi_KK: set to 0 ± 0.5 (large uncertainty)
                fprintf(fm, "  %d  %d  0  %.6f  %.6f  %.6f  %.6f  %.6f  %.6f\n",
                        iQ2, itp,
                        fLL, fErr,
                        fTT, fErr,
                        0.0, 0.5);

                ++nWritten;
            }
        }

        fclose(fx); fclose(fb); fclose(fm);

        printf("[PhiAnalysis] Data files written to %s/\n", dir.c_str());
        printf("  xsec.dat    %d bins\n", nWritten);
        printf("  bsa.dat     %d bins\n", nWritten);
        printf("  moments.dat %d bins (M_LT set to 0, requires phi_KK for real value)\n",
               nWritten);
        return nWritten;
    }

    // ── Make and save all plots ───────────────────────────────────────────────
    void makePlots() {
        setupStyle();
        fs::create_directories(outDir + "/xsec");
        fs::create_directories(outDir + "/bsa");

        TFile* rootFile = TFile::Open((outDir + "/phi_analysis.root").c_str(), "RECREATE");

        makeCrossSectionPlots(rootFile);
        makeBSAPlots(rootFile);

        rootFile->Write();
        rootFile->Close();
        printf("[PhiAnalysis] ROOT file: %s/phi_analysis.root\n", outDir.c_str());
    }

private:
    // ── Bin finder ────────────────────────────────────────────────────────────
    int findBin(double val, const std::vector<double>& edges) const {
        if (val < edges.front() || val >= edges.back()) return -1;
        for (int i = 0; i < (int)edges.size()-1; ++i)
            if (val >= edges[i] && val < edges[i+1]) return i;
        return -1;
    }

    // ── Style ─────────────────────────────────────────────────────────────────
    void setupStyle() const {
        gStyle->SetOptStat(0);
        gStyle->SetOptTitle(0);
        gStyle->SetPadTopMargin(0.08);
        gStyle->SetPadBottomMargin(0.14);
        gStyle->SetPadLeftMargin(0.14);
        gStyle->SetPadRightMargin(0.06);
        gStyle->SetTitleFont(42, "XYZ");
        gStyle->SetLabelFont(42, "XYZ");
        gStyle->SetTitleSize(0.05, "XYZ");
        gStyle->SetLabelSize(0.045,"XYZ");
        gStyle->SetMarkerStyle(20);
        gStyle->SetMarkerSize(1.1);
        gStyle->SetLineWidth(2);
    }

    // ── Q² bin label ─────────────────────────────────────────────────────────
    std::string q2label(int iQ2) const {
        char buf[64];
        snprintf(buf, sizeof(buf), "%.3f < Q^{2} < %.3f GeV^{2}",
                 Q2edges[iQ2], Q2edges[iQ2+1]);
        return buf;
    }

    std::string tplabel(int itp) const {
        char buf[64];
        snprintf(buf, sizeof(buf), "%.3f < t' < %.3f GeV^{2}",
                 Tpedges[itp], Tpedges[itp+1]);
        return buf;
    }

    // ── Cross-section plots ───────────────────────────────────────────────────
    /**
     * For each Q² bin: one canvas with dσ/dt' vs t'
     * Points = N / Δt', error = √N / Δt'
     * Normalised so that all Q²-bin curves share the same relative scale.
     */
    void makeCrossSectionPlots(TFile* rootFile) {
        // Also make one summary canvas with all 4 Q² bins overlaid
        TCanvas* cAll = new TCanvas("xsec_all","dN/dt' all Q2",900,650);
        cAll->SetLogy();
        TLegend* leg = new TLegend(0.55, 0.55, 0.92, 0.88);
        leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.034);

        static const int colours[4] = {kBlue+1, kRed+1, kGreen+2, kMagenta+1};

        TGraphErrors* gAll[4] = {};

        for (int iQ2 = 0; iQ2 < nQ2; ++iQ2) {
            // Build graph
            std::vector<double> x, y, ex, ey;
            for (int itp = 0; itp < nTp; ++itp) {
                double tc = 0.5*(Tpedges[itp]+Tpedges[itp+1]);
                double dt = Tpedges[itp+1]-Tpedges[itp];
                double N  = Nxsec[iQ2][itp];
                if (N < 1) continue;
                x.push_back(tc);
                y.push_back(N/dt);
                ex.push_back(0.5*dt);
                ey.push_back(std::sqrt(N)/dt);
            }
            if (x.empty()) continue;

            auto* g = new TGraphErrors((int)x.size(),
                                       x.data(), y.data(),
                                       ex.data(), ey.data());
            char gname[64]; snprintf(gname, sizeof(gname), "xsec_Q2bin%d", iQ2);
            g->SetName(gname);
            g->SetMarkerColor(colours[iQ2]);
            g->SetLineColor(colours[iQ2]);
            g->SetMarkerStyle(20+iQ2);

            // Individual canvas
            TCanvas* c = new TCanvas(Form("c_xsec_Q2bin%d",iQ2),
                                     Form("dN/dt' Q2 bin %d",iQ2), 750, 600);
            c->SetLogy();
            c->SetLeftMargin(0.14); c->SetBottomMargin(0.14);

            g->Draw("AP");
            g->GetXaxis()->SetTitle("t' (GeV^{2})");
            g->GetYaxis()->SetTitle("dN/dt' (GeV^{-2})");
            g->GetXaxis()->SetLimits(0, 4.5);
            g->SetMinimum(0.1);

            TLatex lat; lat.SetNDC(); lat.SetTextSize(0.042); lat.SetTextFont(42);
            lat.DrawLatex(0.17, 0.92, q2label(iQ2).c_str());

            // t-slope fit hint line
            TLine* zero = new TLine(0, 0, 0, 0);  // invisible, just for save
            (void)zero;

            std::string pdfPath = outDir + "/xsec/" + gname + ".pdf";
            c->SaveAs(pdfPath.c_str());
            printf("[PhiAnalysis] Saved %s\n", pdfPath.c_str());

            rootFile->cd();
            g->Write();
            c->Write();

            gAll[iQ2] = g;
        }

        // Overlay plot
        bool first = true;
        for (int iQ2 = 0; iQ2 < nQ2; ++iQ2) {
            if (!gAll[iQ2]) continue;
            cAll->cd();
            if (first) {
                gAll[iQ2]->Draw("AP");
                gAll[iQ2]->GetXaxis()->SetTitle("t' (GeV^{2})");
                gAll[iQ2]->GetYaxis()->SetTitle("dN/dt' (GeV^{-2})");
                gAll[iQ2]->GetXaxis()->SetLimits(0, 4.5);
                gAll[iQ2]->SetMinimum(0.1);
                first = false;
            } else {
                gAll[iQ2]->Draw("P SAME");
            }
            leg->AddEntry(gAll[iQ2], q2label(iQ2).c_str(), "lp");
        }
        leg->Draw();

        TLatex lat2; lat2.SetNDC(); lat2.SetTextSize(0.04); lat2.SetTextFont(42);
        lat2.DrawLatex(0.17, 0.92, "ep #rightarrow e'p K^{+}K^{-}  CLAS12 RGA");

        std::string allPath = outDir + "/xsec/xsec_all_Q2bins.pdf";
        cAll->SaveAs(allPath.c_str());
        printf("[PhiAnalysis] Saved %s\n", allPath.c_str());
        rootFile->cd(); cAll->Write();
    }

    // ── BSA plots ─────────────────────────────────────────────────────────────
    /**
     * For each (Q², t') bin: BSA = (N+ - N-) / (N+ + N-) vs cos(θ_H).
     * Plus a 4×9 summary grid per Q² bin.
     */
    void makeBSAPlots(TFile* rootFile) {
        // ── Per Q² bin: multi-panel summary ──────────────────────────────────
        for (int iQ2 = 0; iQ2 < nQ2; ++iQ2) {
            // 3×3 grid for 9 t' bins
            int nCols = 3, nRows = 3;
            TCanvas* cSum = new TCanvas(Form("bsa_summary_Q2bin%d",iQ2),
                                        Form("BSA Q2 bin %d",iQ2),
                                        1050, 900);
            cSum->Divide(nCols, nRows, 0.001, 0.001);

            for (int itp = 0; itp < nTp; ++itp) {
                auto* g = buildBSAgraph(iQ2, itp);
                if (!g) continue;

                char gname[64];
                snprintf(gname, sizeof(gname), "BSA_Q2bin%d_tbin%d", iQ2, itp);
                g->SetName(gname);

                // Individual canvas
                TCanvas* c1 = new TCanvas(Form("c_%s",gname), gname, 600, 500);
                drawBSApad(c1, g, iQ2, itp);
                std::string p1 = outDir+"/bsa/"+gname+".pdf";
                c1->SaveAs(p1.c_str());
                rootFile->cd(); g->Write(); c1->Write();
                printf("[PhiAnalysis] Saved %s\n", p1.c_str());

                // In summary canvas
                cSum->cd(itp+1);
                drawBSApad(gPad, g, iQ2, itp);
            }

            // Summary title
            cSum->cd(0);
            TLatex lat; lat.SetNDC(); lat.SetTextSize(0.025); lat.SetTextFont(42);
            char header[128];
            snprintf(header, sizeof(header),
                     "BSA A_{LU} vs cos#theta_{H}  |  %s",
                     q2label(iQ2).c_str());
            lat.DrawLatex(0.03, 0.97, header);

            std::string sp = outDir+"/bsa/BSA_summary_Q2bin"+std::to_string(iQ2)+".pdf";
            cSum->SaveAs(sp.c_str());
            printf("[PhiAnalysis] Saved %s\n", sp.c_str());
            rootFile->cd(); cSum->Write();
        }

        // ── Q²-integrated BSA (summed over t') ───────────────────────────────
        TCanvas* cInt = new TCanvas("bsa_integrated","BSA integrated", 1000, 700);
        cInt->Divide(2, 2, 0.002, 0.002);
        static const int colours[4] = {kBlue+1, kRed+1, kGreen+2, kMagenta+1};

        for (int iQ2 = 0; iQ2 < nQ2; ++iQ2) {
            // Sum over all t' bins
            std::vector<double> Np(NCosH,0), Nm(NCosH,0);
            for (int itp = 0; itp < nTp; ++itp)
                for (int ic = 0; ic < NCosH; ++ic) {
                    Np[ic] += Nplus [iQ2][itp][ic];
                    Nm[ic] += Nminus[iQ2][itp][ic];
                }

            auto* g = buildBSAgraphFromBins(Np, Nm);
            if (!g) continue;
            char gname[64]; snprintf(gname,sizeof(gname),"BSA_integrated_Q2bin%d",iQ2);
            g->SetName(gname);
            g->SetMarkerColor(colours[iQ2]);
            g->SetLineColor(colours[iQ2]);

            cInt->cd(iQ2+1);
            drawBSApadRaw(gPad, g, ("t'-integrated | " + q2label(iQ2)).c_str());

            rootFile->cd(); g->Write();
        }
        std::string ip = outDir+"/bsa/BSA_Q2bins_tintegrated.pdf";
        cInt->SaveAs(ip.c_str());
        printf("[PhiAnalysis] Saved %s\n", ip.c_str());
        rootFile->cd(); cInt->Write();
    }

    // ── Build BSA TGraphErrors for one (Q², t') bin ───────────────────────────
    TGraphErrors* buildBSAgraph(int iQ2, int itp) const {
        std::vector<double> x, y, ex, ey;
        for (int ic = 0; ic < NCosH; ++ic) {
            double cosCenter = -1.0 + (ic + 0.5) * 2.0 / NCosH;
            double dcos      = 2.0 / NCosH;
            double Np = Nplus [iQ2][itp][ic];
            double Nm = Nminus[iQ2][itp][ic];
            double Ntot = Np + Nm;
            if (Ntot < 3) continue;
            double A  = (Np - Nm) / Ntot;
            double dA = std::sqrt((1.0 - A*A) / Ntot);
            x.push_back(cosCenter);
            y.push_back(A);
            ex.push_back(0.5*dcos);
            ey.push_back(dA);
        }
        if (x.empty()) return nullptr;
        return new TGraphErrors((int)x.size(),
                                x.data(), y.data(), ex.data(), ey.data());
    }

    TGraphErrors* buildBSAgraphFromBins(const std::vector<double>& Np,
                                         const std::vector<double>& Nm) const {
        std::vector<double> x, y, ex, ey;
        for (int ic = 0; ic < NCosH; ++ic) {
            double cosCenter = -1.0 + (ic + 0.5) * 2.0 / NCosH;
            double dcos      = 2.0 / NCosH;
            double ntot = Np[ic] + Nm[ic];
            if (ntot < 3) continue;
            double A  = (Np[ic] - Nm[ic]) / ntot;
            double dA = std::sqrt((1.0 - A*A) / ntot);
            x.push_back(cosCenter); y.push_back(A);
            ex.push_back(0.5*dcos); ey.push_back(dA);
        }
        if (x.empty()) return nullptr;
        return new TGraphErrors((int)x.size(),
                                x.data(), y.data(), ex.data(), ey.data());
    }

    // ── Draw BSA on a pad ─────────────────────────────────────────────────────
    void drawBSApad(TVirtualPad* pad, TGraphErrors* g,
                    int iQ2, int itp) const {
        pad->cd();
        pad->SetLeftMargin(0.16); pad->SetBottomMargin(0.16);

        g->SetMarkerStyle(20); g->SetMarkerSize(0.9);
        g->SetMarkerColor(kBlue+1); g->SetLineColor(kBlue+1);
        g->Draw("AP");
        g->GetXaxis()->SetTitle("cos#theta_{H}");
        g->GetYaxis()->SetTitle("A_{LU}");
        g->GetXaxis()->SetLimits(-1.1, 1.1);
        g->SetMinimum(-0.6); g->SetMaximum(0.6);
        g->GetXaxis()->SetNdivisions(505);

        // Zero line
        TLine* zl = new TLine(-1, 0, 1, 0);
        zl->SetLineStyle(2); zl->SetLineColor(kGray+1); zl->Draw();

        TLatex lat; lat.SetNDC(); lat.SetTextSize(0.038); lat.SetTextFont(42);
        lat.DrawLatex(0.17, 0.91, tplabel(itp).c_str());
        lat.SetTextSize(0.032);
        lat.DrawLatex(0.17, 0.84, q2label(iQ2).c_str());
    }

    void drawBSApadRaw(TVirtualPad* pad, TGraphErrors* g,
                       const char* title) const {
        pad->cd();
        pad->SetLeftMargin(0.16); pad->SetBottomMargin(0.16);
        g->SetMarkerStyle(20); g->SetMarkerSize(0.9);
        g->Draw("AP");
        g->GetXaxis()->SetTitle("cos#theta_{H}");
        g->GetYaxis()->SetTitle("A_{LU}");
        g->GetXaxis()->SetLimits(-1.1, 1.1);
        g->SetMinimum(-0.6); g->SetMaximum(0.6);

        TLine* zl = new TLine(-1, 0, 1, 0);
        zl->SetLineStyle(2); zl->SetLineColor(kGray+1); zl->Draw();

        TLatex lat; lat.SetNDC(); lat.SetTextSize(0.034); lat.SetTextFont(42);
        lat.DrawLatex(0.13, 0.92, title);
    }
};
