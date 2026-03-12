/**
 * Chi2Fitter.h
 *
 * Minuit2-based χ² minimiser for the full phi SDME fit.
 *
 * Free parameters (subset can be fixed):
 *   0  A_L    — σ_L amplitude
 *   1  b_L    — σ_L t-slope [GeV⁻²]
 *   2  Im_LL  — BSA sin(φ) from σ_LT' (LL sector)
 *   3  A_T    — σ_T amplitude
 *   4  b_T    — σ_T t-slope
 *   5  Im_TT  — BSA sin(φ) from TT interference
 *   6  A_LT   — LT interference real part
 *   7  A_LTi  — LT interference imaginary part (BSA)
 *
 * The chi2 function:
 *
 *   χ²_total = w_xsec × χ²_xsec   (dσ/dt shape in Q² bins)
 *            + w_bsa  × χ²_BSA    (A_LU in t' bins)
 *            + w_mom  × χ²_moments (f_LL, f_TT, M_LT)
 *
 * where each component is the standard sum of ((data-MC)/err)².
 *
 * The fit uses the reweighting strategy: events are pre-generated once
 * (in EventCache), and at each Minuit step only the per-event weights
 * are recomputed — no new events are generated.
 *
 * Usage:
 *   Chi2Fitter fitter(cache, data, bins);
 *   fitter.setWeights(1.0, 2.0, 1.0);   // xsec, BSA, moments
 *   fitter.fixParam(1);                 // fix b_L
 *   auto result = fitter.fit();
 *   result.bestParams.print();
 */

#pragma once
#include "PhysicsParams.h"
#include "EventData.h"
#include "BinDef.h"
#include "AnalysisModule.h"
#include "DataSet.h"
#include "ProgressBar.h"

// ROOT Minuit2 — only include when ROOT is available
#ifndef PHI_SDME_NO_ROOT
#  include "Math/Minimizer.h"
#  include "Math/Factory.h"
#  include "Math/Functor.h"
#endif

#include <vector>
#include <string>
#include <functional>
#include <cstdio>
#include <cmath>
#include <memory>

// ── Parameter descriptor ──────────────────────────────────────────────────────
struct FitParam {
    std::string name;
    double initVal;
    double stepSize;
    double lowerBound;
    double upperBound;
    bool   fixed = false;
};

// ── Fit result ────────────────────────────────────────────────────────────────
struct FitResult {
    PhysicsParams    bestParams;
    double           chi2min   = -1.0;
    int              ndf       = 0;
    bool             converged = false;

    std::vector<std::string> paramNames;
    std::vector<double>      bestValues;
    std::vector<double>      errors;

    void print() const {
        printf("\n===== Fit Result =====\n");
        printf("  chi2/ndf = %.3f / %d = %.3f\n",
               chi2min, ndf, ndf > 0 ? chi2min/ndf : -1.0);
        printf("  Converged: %s\n", converged ? "YES" : "NO");
        printf("  Parameters:\n");
        for (int i = 0; i < (int)paramNames.size(); ++i)
            printf("    %-12s = %10.4f ± %.4f\n",
                   paramNames[i].c_str(), bestValues[i], errors[i]);
        printf("======================\n\n");
    }
};

// ── Chi2Fitter ────────────────────────────────────────────────────────────────
class Chi2Fitter {
public:
    // Relative weights for each observable in chi2
    double w_xsec   = 1.0;
    double w_bsa    = 2.0;    // BSA gets extra weight — it's the most discriminating
    double w_moments = 1.0;

    // -------------------------------------------------------------------------
    Chi2Fitter(const EventCache& cache,
               const DataSet&    data,
               const BinDef&     bins,
               const PhysicsParams& initParams)
        : cache_(cache), data_(data), bins_(bins),
          analysis_(bins), initParams_(initParams),
          spinner_("Fitting", 0.15)
    {
        setupDefaultParams();
    }

    // -------------------------------------------------------------------------
    // Fix a parameter by index (it will not be varied by Minuit)
    void fixParam(int idx)   { params_[idx].fixed = true;  }
    void freeParam(int idx)  { params_[idx].fixed = false; }

    // Override initial value and step size
    void setParam(int idx, double val, double step,
                  double lo = -1e9, double hi = 1e9) {
        params_[idx].initVal    = val;
        params_[idx].stepSize   = step;
        params_[idx].lowerBound = lo;
        params_[idx].upperBound = hi;
    }

    void setWeights(double wxsec, double wbsa, double wmom) {
        w_xsec = wxsec; w_bsa = wbsa; w_moments = wmom;
    }

    // -------------------------------------------------------------------------
    // Evaluate chi2 at given parameter vector (called by Minuit at each step)
    double chi2(const double* par) const
    {
        // Unpack Minuit parameters → PhysicsParams
        PhysicsParams p = paramsFromVector(par);

        // Compute MC predictions
        auto mc = analysis_.compute(cache_, p);

        double total = 0.0;
        int    ndf   = 0;

        // ── χ²_xsec ───────────────────────────────────────────────────
        for (const auto& dp : data_.xsec) {
            int idx = bins_.flatIndex(dp.iQ2, dp.it, dp.ixB);
            if (idx < 0 || idx >= (int)mc.size()) continue;
            const auto& r = mc[idx];
            if (r.dsigma_dt <= 0 || dp.error <= 0) continue;
            double pull = (dp.value - r.dsigma_dt) / dp.error;
            total += w_xsec * pull * pull;
            ndf++;
        }

        // ── χ²_BSA ────────────────────────────────────────────────────
        for (const auto& dp : data_.bsa) {
            int idx = bins_.flatIndex(dp.iQ2, dp.it, dp.ixB);
            if (idx < 0 || idx >= (int)mc.size()) continue;
            const auto& r = mc[idx];
            if (dp.error <= 0) continue;
            double pull = (dp.value - r.A_LU) / dp.error;
            total += w_bsa * pull * pull;
            ndf++;
        }

        // ── χ²_moments ────────────────────────────────────────────────
        for (const auto& mp : data_.moments) {
            int idx = bins_.flatIndex(mp.iQ2, mp.it, mp.ixB);
            if (idx < 0 || idx >= (int)mc.size()) continue;
            const auto& r = mc[idx];

            if (mp.fLL_err > 0) {
                double pull = (mp.fLL - r.f_LL) / mp.fLL_err;
                total += w_moments * pull * pull; ndf++;
            }
            if (mp.fTT_err > 0) {
                double pull = (mp.fTT - r.f_TT) / mp.fTT_err;
                total += w_moments * pull * pull; ndf++;
            }
            if (mp.M_LT_err > 0) {
                double pull = (mp.M_LT - r.M_LT) / mp.M_LT_err;
                total += w_moments * pull * pull; ndf++;
            }
        }

        spinner_.tick(total);   // update spinner with current chi2
        return total;
    }

    // -------------------------------------------------------------------------
    // Run the full fit
    FitResult fit(const std::string& strategy = "Migrad")
    {
        // Count free parameters
        int nFree = 0;
        for (const auto& p : params_) if (!p.fixed) nFree++;
        printf("[Chi2Fitter] Running %s with %d free parameters\n",
               strategy.c_str(), nFree);
        printf("[Chi2Fitter] Cache size: %ld events\n",  (long)cache_.size());
        printf("[Chi2Fitter] Data:  %zu xsec, %zu BSA, %zu moment points\n",
               data_.xsec.size(), data_.bsa.size(), data_.moments.size());

        // Build Minuit2 minimiser
        auto minimiser = std::unique_ptr<ROOT::Math::Minimizer>(
            ROOT::Math::Factory::CreateMinimizer("Minuit2", strategy.c_str()));

        minimiser->SetMaxFunctionCalls(10000);
        minimiser->SetMaxIterations(5000);
        minimiser->SetTolerance(0.01);
        minimiser->SetPrintLevel(-1);   // silence Minuit; FitSpinner shows progress

        // Wrap chi2 as ROOT functor
        int nPar = static_cast<int>(params_.size());
        ROOT::Math::Functor fcn([this](const double* p){ return chi2(p); }, nPar);
        minimiser->SetFunction(fcn);

        // Set parameters
        for (int i = 0; i < nPar; ++i) {
            const auto& fp = params_[i];
            minimiser->SetLimitedVariable(i, fp.name, fp.initVal, fp.stepSize,
                                          fp.lowerBound, fp.upperBound);
            if (fp.fixed) minimiser->FixVariable(i);
        }

        // ── Minimise ───────────────────────────────────────────────────
        bool ok = minimiser->Minimize();
        spinner_.done(minimiser->MinValue());

        // ── Build result ───────────────────────────────────────────────
        FitResult result;
        result.converged = ok;
        result.chi2min   = minimiser->MinValue();
        result.ndf       = static_cast<int>(data_.xsec.size() +
                                            data_.bsa.size()  +
                                            3 * data_.moments.size()) - nFree;

        const double* bestPar = minimiser->X();
        const double* bestErr = minimiser->Errors();

        result.bestParams = paramsFromVector(bestPar);
        for (int i = 0; i < nPar; ++i) {
            result.paramNames.push_back(params_[i].name);
            result.bestValues.push_back(bestPar[i]);
            result.errors.push_back(bestErr[i]);
        }

        return result;
    }

    // -------------------------------------------------------------------------
    // Evaluate chi2 at the initial parameter point (sanity check before fitting)
    double chi2AtStart() const {
        std::vector<double> v(params_.size());
        for (int i = 0; i < (int)params_.size(); ++i)
            v[i] = params_[i].initVal;
        return chi2(v.data());
    }

    // -------------------------------------------------------------------------
    // Parameter names (for bookkeeping)
    std::vector<std::string> paramNames() const {
        std::vector<std::string> names;
        for (const auto& p : params_) names.push_back(p.name);
        return names;
    }

private:
    const EventCache&    cache_;
    const DataSet&       data_;
    const BinDef&        bins_;
    AnalysisModule       analysis_;
    PhysicsParams        initParams_;
    std::vector<FitParam> params_;
    mutable FitSpinner   spinner_;   // mutable: updated inside const chi2()

    // ── Minuit parameter index → field mapping ─────────────────────────
    enum ParIdx {
        kAL=0, kbL, kImLL,
        kAT,   kbT,  kImTT,
        kALT,  kALTi,
        kNPar
    };

    void setupDefaultParams() {
        params_.resize(kNPar);
        auto& p = initParams_;
        //
        // Bounds are physically motivated:
        //   A_L, A_T  — cross-section amplitudes, must be positive, O(1–100)
        //   b_L, b_T  — t-slopes in [1,10] GeV^-2 (phi electroproduction)
        //   Im_*      — BSA-driving imaginary parts.
        //               CRITICAL: bounds must be tight enough that Minuit
        //               cannot escape to a flat region of chi2.
        //               Rule of thumb: |Im| < 2 × amplitude
        //               Im_LL:  A_L ~ 30  → [-60, +60]
        //               Im_TT:  A_T ~ 10  → [-20, +20]
        //               A_LTi:  A_LT ~ 3  → [-20, +20]
        //               Wide bounds (±200) are the #1 cause of boundary hits.
        //
        params_[kAL]   = {"A_L",   p.A_L,   1.0,    0.0,  100.0};
        params_[kbL]   = {"b_L",   p.b_L,   0.1,    1.0,   10.0};
        params_[kImLL] = {"Im_LL", p.Im_LL, 0.5,  -60.0,   60.0};
        params_[kAT]   = {"A_T",   p.A_T,   0.5,    0.0,  100.0};
        params_[kbT]   = {"b_T",   p.b_T,   0.1,    1.0,   10.0};
        params_[kImTT] = {"Im_TT", p.Im_TT, 0.2,  -20.0,   20.0};
        params_[kALT]  = {"A_LT",  p.A_LT,  0.5,  -30.0,   30.0};
        params_[kALTi] = {"A_LTi", p.A_LTi, 0.5,  -30.0,   30.0};
    }

    PhysicsParams paramsFromVector(const double* par) const {
        PhysicsParams p = initParams_;   // copy fixed fields (n_L, n_T, etc.)
        p.A_L    = par[kAL];
        p.b_L    = par[kbL];
        p.Im_LL  = par[kImLL];
        p.A_T    = par[kAT];
        p.b_T    = par[kbT];
        p.Im_TT  = par[kImTT];
        p.A_LT   = par[kALT];
        p.A_LTi  = par[kALTi];
        return p;
    }
};
