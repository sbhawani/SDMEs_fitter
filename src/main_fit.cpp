/**
 * main_fit.cpp
 *
 * phi SDME fitting pipeline — two clearly separated paths:
 *
 * ── PATH 1: Closure test ─────────────────────────────────────────────────────
 *
 *   Purpose:  Verify that the fitter recovers the true parameters used to
 *             generate the events.  No real data, no FastMC required.
 *
 *   Workflow:
 *     1. phi_gen generates events with known params  → events.lund
 *     2. EventCache::fillFromLund(events.lund)       → cache (truth-level)
 *     3. AnalysisModule evaluates model at TRUTH params → mc_truth
 *     4. DataSet::fillFromMC(mc_truth, err=5%)       → pseudo-data
 *     5. Fitter starts 10% off truth, minimises chi2
 *     6. Check: fitted params ≈ truth params (pulls ~0)
 *
 *   Usage:
 *     ./phi_fit --closure output/events.lund [--config config/default.cfg]
 *               [--output-dir output/closure]
 *
 * ── PATH 2: Real experimental data ───────────────────────────────────────────
 *
 *   Purpose:  Fit real CLAS12 cross-section and BSA measurements.
 *             Uses fastMC.jar output for acceptance correction.
 *
 *   Workflow:
 *     1. phi_gen generates MC events            → events_mc.lund
 *     2. fastMC.jar filters                     → events_fmc.lund
 *     3. EventCache::fillFromLund(fmc, fastmc=true) → cache (detector-accepted)
 *     4. Load real data: xsec.dat + bsa.dat
 *     5. Fitter minimises chi2 between real data and acceptance-corrected MC
 *
 *   Usage:
 *     ./phi_fit --data xsec.dat bsa.dat --mc-lund events_fmc.lund
 *               [--config config/default.cfg] [--output-dir output/real]
 *               [--moments moments.dat]
 *
 * ── Options ─────────────────────────────────────────────────────────────────
 *   --closure  <lund>       Path 1: closure test LUND file (truth-level)
 *   --mc-lund  <lund>       Path 2: MC LUND file (truth or fastMC-filtered)
 *   --fastmc                Mark --mc-lund as fastMC.jar output (apply hit cuts)
 *   --data  <xsec> <bsa>    Path 2: real data files
 *   --moments <mom.dat>     Optional angular moment data (both paths)
 *   --config  <cfg>         Physics parameter config file
 *   --output-dir <dir>      Where to save fit results (default: output)
 *   --stat-err <frac>       Pseudo-data stat error fraction (default: 0.05)
 *   --perturb <frac>        Starting point offset for closure test (default: 0.10)
 *   --help                  Show this message
 */

#include <cstdio>
#include <string>
#include <memory>
#include <filesystem>

#include "TRandom3.h"

#include "PhysicsParams.h"
#include "BinDef.h"
#include "EventCache.h"
#include "AnalysisModule.h"
#include "DataSet.h"
#include "Chi2Fitter.h"
#include "FitPlotter.h"

namespace fs = std::filesystem;

static void usage(const char* prog) {
    printf("\nUsage:\n");
    printf("  PATH 1 — Closure test (no real data, no FastMC):\n");
    printf("    %s --closure <events.lund> [options]\n\n", prog);
    printf("  PATH 2 — Real data fit:\n");
    printf("    %s --data <xsec.dat> <bsa.dat> --mc-lund <events_fmc.lund> [options]\n\n", prog);
    printf("Options:\n");
    printf("  --closure  <lund>       Truth-level LUND for closure test\n");
    printf("  --mc-lund  <lund>       MC LUND for model/acceptance (path 2)\n");
    printf("  --fastmc                --mc-lund is fastMC.jar output (apply detector cuts)\n");
    printf("  --data  <xsec> <bsa>    Real experimental data files\n");
    printf("  --moments <mom.dat>     Angular moment data (optional)\n");
    printf("  --config  <cfg>         Physics parameter config (default: built-in)\n");
    printf("  --output-dir <dir>      Output directory (default: output)\n");
    printf("  --stat-err <frac>       Pseudo-data stat error fraction (default: 0.05)\n");
    printf("  --perturb <frac>        Closure starting offset fraction (default: 0.10)\n");
    printf("  --help                  Show this message\n\n");
    printf("Examples:\n");
    printf("  # Closure test:\n");
    printf("  ./phi_fit --closure output/events.lund --config config/default.cfg\n\n");
    printf("  # Real data with acceptance from fastMC:\n");
    printf("  ./phi_fit --data data/xsec.dat data/bsa.dat \\\n");
    printf("            --mc-lund output/events_fmc.lund --fastmc\n\n");
}

// ── Perturb params for closure test starting point ────────────────────────────
static PhysicsParams perturb(const PhysicsParams& p, double frac) {
    PhysicsParams q = p;
    // Perturb each amplitude by ±frac, slopes and powers by ±frac/2
    q.A_L    *= (1.0 + frac);
    q.A_T    *= (1.0 - frac);
    q.Im_LL  *= (1.0 + frac * 0.8);
    q.Im_TT  *= (1.0 - frac * 0.8);
    q.A_LT   *= (1.0 + frac * 0.6);
    q.A_LTi  *= (1.0 - frac * 0.6);
    q.b_L    *= (1.0 + frac * 0.4);
    q.b_T    *= (1.0 - frac * 0.4);
    return q;
}

int main(int argc, char* argv[])
{
    printf("\n");
    printf("╔══════════════════════════════════════════════════════╗\n");
    printf("║   phi SDME Fitter  —  CLAS12 RGA phi electroproduction\n");
    printf("╚══════════════════════════════════════════════════════╝\n\n");

    if (argc < 2) { usage(argv[0]); return 1; }

    // ── Argument parsing ───────────────────────────────────────────────────────
    enum class Mode { NONE, CLOSURE, REAL_DATA } mode = Mode::NONE;

    std::string configFile;
    std::string closureLund;
    std::string mcLund;
    std::string dataXsec, dataBSA, dataMom;
    std::string outputDir = "output";
    bool   mcFastmc  = false;
    double statErr   = 0.05;
    double perturbFrac = 0.10;

    for (int i = 1; i < argc; ++i) {
        std::string a = argv[i];
        if      (a == "--help" || a == "-h")           { usage(argv[0]); return 0; }
        else if (a == "--closure"  && i+1 < argc)      { closureLund = argv[++i]; mode = Mode::CLOSURE; }
        else if (a == "--mc-lund"  && i+1 < argc)      mcLund     = argv[++i];
        else if (a == "--fastmc")                       mcFastmc   = true;
        else if (a == "--data"     && i+2 < argc)      { dataXsec = argv[++i]; dataBSA = argv[++i]; mode = Mode::REAL_DATA; }
        else if (a == "--moments"  && i+1 < argc)      dataMom    = argv[++i];
        else if (a == "--config"   && i+1 < argc)      configFile = argv[++i];
        else if (a == "--output-dir" && i+1 < argc)    outputDir  = argv[++i];
        else if (a == "--stat-err"  && i+1 < argc)     statErr    = std::stod(argv[++i]);
        else if (a == "--perturb"   && i+1 < argc)     perturbFrac = std::stod(argv[++i]);
        else { fprintf(stderr, "[phi_fit] Unknown option: %s\n", a.c_str()); usage(argv[0]); return 1; }
    }

    if (mode == Mode::NONE) {
        fprintf(stderr, "[phi_fit] ERROR: specify --closure <lund> or --data <xsec> <bsa>\n\n");
        usage(argv[0]); return 1;
    }
    if (mode == Mode::REAL_DATA && mcLund.empty()) {
        fprintf(stderr, "[phi_fit] ERROR: --data requires --mc-lund <lund>\n\n");
        usage(argv[0]); return 1;
    }

    fs::create_directories(outputDir);

    // ── 1. Load physics parameters ─────────────────────────────────────────────
    PhysicsParams params;
    if (!configFile.empty()) {
        try {
            params.load(configFile);
            printf("[phi_fit] Config: %s\n", configFile.c_str());
        } catch (const std::exception& e) {
            fprintf(stderr, "[phi_fit] ERROR loading config: %s\n", e.what());
            return 1;
        }
    } else {
        printf("[phi_fit] No config given — using built-in defaults.\n");
    }
    params.print();

    // ── 2. Define bins ──────────────────────────────────────────────────────────
    BinDef bins;
    printf("[phi_fit] Bins: Q2=%d  t'=%d  xB=%d  total=%d\n\n",
           bins.nQ2(), bins.nT(), bins.nXB(), bins.nTotal());

    // ── 3. Fill EventCache from LUND file ───────────────────────────────────────
    EventCache cache;

    if (mode == Mode::CLOSURE) {
        printf("╔══════════════════════════════════════════════╗\n");
        printf("║  PATH 1: CLOSURE TEST                        ║\n");
        printf("║  MC truth LUND + no detector cuts            ║\n");
        printf("╚══════════════════════════════════════════════╝\n\n");
        printf("[phi_fit] Reading truth LUND: %s\n", closureLund.c_str());
        int n = cache.fillFromLund(closureLund, bins, false);
        if (n == 0) {
            fprintf(stderr, "[phi_fit] ERROR: no events cached from %s\n", closureLund.c_str());
            return 1;
        }
        printf("[phi_fit] Cached %d truth events\n\n", n);
    } else {
        printf("╔══════════════════════════════════════════════╗\n");
        printf("║  PATH 2: REAL DATA FIT                       ║\n");
        printf("║  Experimental data + MC acceptance           ║\n");
        printf("╚══════════════════════════════════════════════╝\n\n");
        printf("[phi_fit] Reading MC LUND: %s  (fastmc=%s)\n",
               mcLund.c_str(), mcFastmc ? "ON" : "OFF");
        int n = cache.fillFromLund(mcLund, bins, mcFastmc);
        if (n == 0) {
            fprintf(stderr, "[phi_fit] ERROR: no events cached from %s\n", mcLund.c_str());
            if (mcFastmc)
                fprintf(stderr, "  All events rejected by detector cuts.\n"
                                "  Is this genuine fastMC.jar output?\n");
            return 1;
        }
        printf("[phi_fit] Cached %d MC events\n\n", n);
    }

    // ── 4. AnalysisModule: evaluate model at generator (truth) params ───────────
    AnalysisModule analysis(bins);
    auto mc_truth = analysis.compute(cache, params);

    printf("[phi_fit] Model at generator params:\n");
    AnalysisModule::printResults(mc_truth, bins);

    // ── Sanity check: warn if BSA is saturated ──────────────────────────────────
    {
        double bsa_max = 0.0, bsa_rms = 0.0;
        int n = 0;
        for (const auto& r : mc_truth) {
            if (r.nEvents < 5) continue;
            double a = std::abs(r.A_LU);
            bsa_max  = std::max(bsa_max, a);
            bsa_rms += r.A_LU * r.A_LU;
            ++n;
        }
        if (n > 0) bsa_rms = std::sqrt(bsa_rms / n);
        printf("[phi_fit] BSA sanity: max|A_LU|=%.3f  rms=%.3f\n", bsa_max, bsa_rms);
        if (bsa_max > 0.5) {
            fprintf(stderr,
                "\n[phi_fit] WARNING: BSA is saturated (max|A_LU|=%.2f > 0.5)!\n"
                "  The Im parameters (Im_LL, Im_TT, A_LTi) are too large.\n"
                "  The chi² surface will be flat in these directions → fit will fail.\n"
                "  Fix: reduce Im_LL, Im_TT, A_LTi in the config so that\n"
                "       |A_LU| ~ 0.05–0.20 (physical for phi at CLAS12).\n"
                "  Hint: Im_LL / (A_L + A_T) controls the BSA scale.\n"
                "        With A_L=%.1f, A_T=%.1f → Im_LL < %.1f for |A_LU| < 0.3\n\n",
                bsa_max, params.A_L, params.A_T, 0.3*(params.A_L+params.A_T));
        } else {
            printf("[phi_fit] BSA looks physical — closure test should work.\n\n");
        }
    }

    // ── 5. Build DataSet ────────────────────────────────────────────────────────
    DataSet data;
    PhysicsParams truthParams = params;

    if (mode == Mode::CLOSURE) {
        // -- Closure: pseudo-data from the MC truth model --
        printf("[phi_fit] Generating pseudo-data from MC model (stat_err=%.0f%%)...\n\n",
               statErr * 100);
        data.fillFromMC(mc_truth, statErr, /*seed=*/1234);

        // Save pseudo-data files so user can inspect them
        std::string pdDir = outputDir + "/pseudo_data";
        fs::create_directories(pdDir);
        // Write in DataSet format
        {
            FILE* fx = fopen((pdDir+"/xsec.dat").c_str(), "w");
            FILE* fb = fopen((pdDir+"/bsa.dat").c_str(),  "w");
            if (fx && fb) {
                fprintf(fx,"# pseudo-data xsec (from MC truth + %.0f%% stat smearing)\n",
                        statErr*100);
                fprintf(fb,"# pseudo-data BSA  (from MC truth + %.0f%% stat smearing)\n",
                        statErr*100);
                fprintf(fx,"# iQ2  it  ixB  dσ/dt  err\n");
                fprintf(fb,"# iQ2  it  ixB  A_LU   err\n");
                for (const auto& dp : data.xsec)
                    fprintf(fx,"  %d  %d  %d  %.6e  %.6e\n",
                            dp.iQ2,dp.it,dp.ixB,dp.value,dp.error);
                for (const auto& db : data.bsa)
                    fprintf(fb,"  %d  %d  %d  %.6f  %.6f\n",
                            db.iQ2,db.it,db.ixB,db.value,db.error);
                fclose(fx); fclose(fb);
                printf("[phi_fit] Pseudo-data saved to %s/\n\n", pdDir.c_str());
            }
        }

        // Perturb starting params for the fitter (to make the test non-trivial)
        params = perturb(truthParams, perturbFrac);
        printf("[phi_fit] Starting params perturbed by %.0f%% from truth:\n",
               perturbFrac*100);
        params.print();

    } else {
        // -- Real data: load from files --
        printf("[phi_fit] Loading real data:\n");
        printf("  xsec : %s\n", dataXsec.c_str());
        printf("  BSA  : %s\n", dataBSA.c_str());
        data.loadXsec(dataXsec);
        data.loadBSA(dataBSA);
        if (!dataMom.empty()) {
            printf("  mom  : %s\n", dataMom.c_str());
            data.loadMoments(dataMom);
        }
        printf("\n");
    }

    if (data.empty()) {
        fprintf(stderr, "[phi_fit] ERROR: DataSet is empty after loading.\n");
        return 1;
    }

    // ── 6. Evaluate MC at starting params ──────────────────────────────────────
    auto mc_start = analysis.compute(cache, params);
    printf("[phi_fit] Model at starting params:\n");
    AnalysisModule::printResults(mc_start, bins);

    // ── 7. Build fitter and run ─────────────────────────────────────────────────
    Chi2Fitter fitter(cache, data, bins, params);
    fitter.setWeights(1.0, 2.0, 1.0);  // BSA gets 2× weight

    // Fix LT amplitudes if SCHC config
    if (params.A_LT == 0.0 && params.A_LTi == 0.0) {
        printf("[phi_fit] SCHC config detected (A_LT=A_LTi=0) — fixing LT params\n");
        fitter.fixParam(6); fitter.fixParam(7);
    }

    printf("[phi_fit] chi2 at starting point: %.2f\n\n", fitter.chi2AtStart());
    printf("[phi_fit] Starting Minuit2 minimisation...\n");
    auto result = fitter.fit("Migrad");
    result.print();

    // ── 8. Evaluate MC at best-fit params ──────────────────────────────────────
    auto mc_best = analysis.compute(cache, result.bestParams);
    printf("[phi_fit] Model at best-fit params:\n");
    AnalysisModule::printResults(mc_best, bins);

    // ── 9. Print pull table ────────────────────────────────────────────────────
    printf("\n%-5s %-4s  %-12s  %-12s %-12s %-8s  %-10s %-10s %-8s\n",
           "iQ2","it'","Q2_range","data_xsec","mc_xsec","pull",
           "data_BSA","mc_BSA","pull");
    printf("%s\n", std::string(95,'-').c_str());
    for (const auto& dp : data.xsec) {
        int idx = bins.flatIndex(dp.iQ2, dp.it, dp.ixB);
        if (idx < 0 || idx >= (int)mc_best.size()) continue;
        const auto& r = mc_best[idx];
        double bsa_data=0, bsa_mc=0, bsa_err=1;
        for (const auto& db : data.bsa) {
            if (db.iQ2==dp.iQ2 && db.it==dp.it && db.ixB==dp.ixB) {
                bsa_data=db.value; bsa_err=db.error; bsa_mc=r.A_LU;
            }
        }
        double pull_x = (dp.error>0) ? (dp.value - r.dsigma_dt)/dp.error : 0;
        double pull_b = (bsa_err >0) ? (bsa_data - bsa_mc)/bsa_err       : 0;
        printf("  %2d  %2d    Q2=[%.2f,%.2f]  %10.4e %10.4e %+6.2f   %+8.4f %+8.4f %+6.2f\n",
               dp.iQ2, dp.it,
               bins.Q2bins.lo(dp.iQ2), bins.Q2bins.hi(dp.iQ2),
               dp.value, r.dsigma_dt, pull_x,
               bsa_data, bsa_mc, pull_b);
    }

    // ── 10. Closure test: compare fitted vs truth params ───────────────────────
    if (mode == Mode::CLOSURE) {
        // Gather Minuit bounds for boundary detection
        // Order matches kAL,kbL,kImLL,kAT,kbT,kImTT,kALT,kALTi
        struct ParInfo {
            const char* name;
            double truth, fit, err, lo, hi;
        };
        std::vector<ParInfo> pars = {
            {"A_L",   truthParams.A_L,   result.bestParams.A_L,   result.errors[0],   0.0,  100.0},
            {"b_L",   truthParams.b_L,   result.bestParams.b_L,   result.errors[1],   1.0,   10.0},
            {"Im_LL", truthParams.Im_LL, result.bestParams.Im_LL, result.errors[2], -60.0,   60.0},
            {"A_T",   truthParams.A_T,   result.bestParams.A_T,   result.errors[3],   0.0,  100.0},
            {"b_T",   truthParams.b_T,   result.bestParams.b_T,   result.errors[4],   1.0,   10.0},
            {"Im_TT", truthParams.Im_TT, result.bestParams.Im_TT, result.errors[5], -20.0,   20.0},
            {"A_LT",  truthParams.A_LT,  result.bestParams.A_LT,  result.errors[6], -30.0,   30.0},
            {"A_LTi", truthParams.A_LTi, result.bestParams.A_LTi, result.errors[7], -30.0,   30.0},
        };

        int nPass = 0, nWarn = 0, nFail = 0, nBound = 0;
        printf("\n╔══════════════════════════════════════════════════════════════════════╗\n");
        printf("║  CLOSURE TEST RESULTS                                                ║\n");
        printf("╠═══════════════╦═════════════╦═════════════╦════════╦════════════════╣\n");
        printf("║  Parameter    ║  Truth      ║  Fitted     ║  Pull  ║  Status        ║\n");
        printf("╠═══════════════╬═════════════╬═════════════╬════════╬════════════════╣\n");

        for (const auto& p : pars) {
            double pull = (p.err > 0) ? (p.fit - p.truth) / p.err : 0.0;
            double abspull = std::abs(pull);

            // Detect boundary hit: fitted value within 1% of bound range from edge
            double range = p.hi - p.lo;
            bool atBound = (std::abs(p.fit - p.lo) < 0.01 * range ||
                            std::abs(p.fit - p.hi) < 0.01 * range);

            const char* status;
            if (atBound) {
                status = "BOUNDARY HIT"; ++nBound; ++nFail;
            } else if (abspull < 1.0) {
                status = "PASS";         ++nPass;
            } else if (abspull < 2.0) {
                status = "warn (1-2σ)";  ++nWarn;
            } else {
                status = "FAIL (>2σ)";   ++nFail;
            }

            printf("║  %-13s ║  %10.4g  ║  %10.4g  ║ %+6.2f ║  %-14s  ║\n",
                   p.name, p.truth, p.fit, pull, status);
        }
        printf("╚═══════════════╩═════════════╩═════════════╩════════╩════════════════╝\n");
        printf("  Pull = (fitted - truth) / σ_fit\n");
        printf("  PASS: |pull|<1  |  warn: 1<|pull|<2  |  FAIL: |pull|>2 or boundary\n\n");
        printf("  Summary: %d PASS  %d warn  %d FAIL", nPass, nWarn, nFail);
        if (nBound > 0)
            printf("  (%d boundary hit%s — tighten bounds or increase truth |Im| values)",
                   nBound, nBound > 1 ? "s" : "");
        printf("\n");

        // Note on A_L/b_L anti-correlation
        double pullAL = (result.errors[0] > 0) ?
            std::abs(result.bestParams.A_L - truthParams.A_L) / result.errors[0] : 0.0;
        double pullbL = (result.errors[1] > 0) ?
            std::abs(result.bestParams.b_L - truthParams.b_L) / result.errors[1] : 0.0;
        if (pullAL > 1.0 || pullbL > 1.0) {
            printf("\n  NOTE: A_L and b_L are anti-correlated in t-slope fits.\n");
            printf("        Their individual pulls can be large even when the\n");
            printf("        dσ/dt shape is well-recovered.  Check if (A_L↑,b_L↓)\n");
            printf("        or (A_L↓,b_L↑) — opposite signs confirm anti-correlation.\n");
            double signAL = result.bestParams.A_L - truthParams.A_L;
            double signbL = result.bestParams.b_L - truthParams.b_L;
            if (signAL * signbL < 0)
                printf("        ✓ Confirmed: A_L and b_L pulls have opposite signs.\n");
        }
        printf("\n");
    }

    // ── 11. Save results ───────────────────────────────────────────────────────
    std::string resFile = outputDir + "/fit_result_params.cfg";
    result.bestParams.save(resFile);
    printf("[phi_fit] Best-fit parameters → %s\n", resFile.c_str());

    if (mode == Mode::CLOSURE) {
        std::string truthFile = outputDir + "/fit_truth_params.cfg";
        truthParams.save(truthFile);
        printf("[phi_fit] Truth parameters   → %s\n", truthFile.c_str());
    }

    // ── 12. Plots ──────────────────────────────────────────────────────────────
    printf("\n[phi_fit] Generating fit quality plots...\n");
    {
        FitPlotter plotter(bins, outputDir);
        if (mode == Mode::CLOSURE)
            plotter.plotAll(data, mc_best, result, &truthParams);
        else
            plotter.plotAll(data, mc_best, result, nullptr);
    }

    printf("[phi_fit] Converged: %s  |  chi2/ndf: %.2f\n\n",
           result.converged ? "YES" : "NO",
           result.chi2min / std::max(1, result.ndf));

    return result.converged ? 0 : 2;
}
