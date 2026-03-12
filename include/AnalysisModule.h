/**
 * AnalysisModule.h
 *
 * Computes binned observables from the EventCache using a given PhysicsParams.
 *
 * For each (Q², t', xB) bin the module computes:
 *
 *   ── Cross-section ────────────────────────────────────────────────────────
 *   dσ/dt  =  (1/ΔΩ) × (1/L) × Σ_i w_i(params)
 *   where w_i = σ₃_i × (dsigmaT + ε_i × dsigmaL) × W_UU_i / (4π²)
 *   and the sum is acceptance-corrected by the ratio N_gen/N_reco per bin.
 *
 *   ── Beam-Spin Asymmetry ──────────────────────────────────────────────────
 *   A_LU = Σ_i [Pl_i × W_LU_i × w_i] / Σ_i [W_UU_i × w_i]
 *
 *   ── Angular moments (LL / LT / TT fractions) ─────────────────────────────
 *   f_LL = Σ_i [cos²θ × w_i] / Σ_i w_i
 *   f_TT = Σ_i [sin²θ × w_i] / Σ_i w_i
 *   M_LT = Σ_i [sin2θ × cos(φ_KK) × w_i] / Σ_i w_i   ← from W_UU^LT
 *
 * These moments are compared to data in the χ² function.
 *
 * The key physics:
 *   Under strict SCHC: f_LL + f_TT = 1, M_LT = 0
 *   A_LU driven purely by Im_LL and Im_TT
 *   Nonzero M_LT → direct evidence of L/T interference (SCHC violation)
 */

#pragma once
#include "EventData.h"
#include "BinDef.h"
#include "SDMEModel.h"
#include "PhysicsParams.h"
#include "w_kernels.hpp"
#include <vector>
#include <cmath>
#include <cassert>
#ifdef _OPENMP
#  include <omp.h>
#endif

// ── Per-bin result ────────────────────────────────────────────────────────────
struct BinResult {
    int iQ2 = -1, it = -1, ixB = -1;

    double Q2_center  = 0.0;
    double t_center   = 0.0;
    double xB_center  = 0.0;

    // Cross-section
    double dsigma_dt  = 0.0;   // [nb/GeV²]  acceptance-corrected
    double dsigmaL_dt = 0.0;   // longitudinal piece
    double dsigmaT_dt = 0.0;   // transverse piece

    // BSA
    double A_LU       = 0.0;
    double A_LU_err   = 0.0;   // statistical error estimate

    // Angular moment fractions (normalised to 1)
    double f_LL = 0.0;         // cos²θ moment  (= <cos²θ>_w)
    double f_TT = 0.0;         // sin²θ moment
    double M_LT = 0.0;         // sin2θ cosφ_KK moment (LT interference; nonzero = SCHC violation)
    double R     = 0.0;        // σ_L/σ_T  from angular distribution method (requires SCHC)
    double r04_00 = 0.0;       // SDME element r^04_00  (extracted from cos²θ shape)
    double eps_bar = 0.0;      // event-averaged ε̄ in this bin (NOT a free parameter)

    // Number of events in this bin (for error estimates)
    int    nEvents = 0;
    double sumW    = 0.0;      // sum of weights
};

// ── Analysis module ───────────────────────────────────────────────────────────
class AnalysisModule {
public:

    explicit AnalysisModule(const BinDef& bins) : bins_(bins) {}

    /**
     * Compute all bin results from the EventCache with the given PhysicsParams.
     * This is called at every Minuit step — it must be fast.
     *
     * @param cache   pre-generated accepted events
     * @param params  current trial parameters
     * @return        vector of BinResult, one per (iQ2, it, ixB) bin
     */
    std::vector<BinResult> compute(const EventCache&    cache,
                                   const PhysicsParams& params) const
    {
        int nBins = bins_.nTotal();
        std::vector<BinResult> results(nBins);

        // Initialise bin metadata
        for (int iQ2 = 0; iQ2 < bins_.nQ2(); ++iQ2)
        for (int it  = 0; it  < bins_.nT();  ++it )
        for (int ixB = 0; ixB < bins_.nXB(); ++ixB) {
            int idx = bins_.flatIndex(iQ2, it, ixB);
            results[idx].iQ2      = iQ2;
            results[idx].it       = it;
            results[idx].ixB      = ixB;
            results[idx].Q2_center = bins_.Q2bins.center(iQ2);
            results[idx].t_center  = bins_.tbins.center(it);
            results[idx].xB_center = bins_.xBbins.center(ixB);
        }

        // Accumulators
        //  [0] sum W_UU                   → cross-section numerator
        //  [1] sum Pl×W_LU×W_UU_w        → BSA numerator
        //  [2] sum W_UU                   → BSA denominator (same as [0])
        //  [3] sum cos²θ × W_UU           → f_LL numerator
        //  [4] sum sin²θ × W_UU           → f_TT numerator
        //  [5] sum sin2θ cosφ × W_UU      → M_LT numerator
        //  [6] sum W_LU (for helicity-odd cross check)
        struct Acc {
            double sumWUU     = 0.0;
            double sumWLUPol  = 0.0;
            double sumLL      = 0.0;
            double sumTT      = 0.0;
            double sumLT      = 0.0;
            double sumSigmaL  = 0.0;
            double sumSigmaT  = 0.0;
            // For proper r^04_00 extraction via angular distribution method:
            //   W(cosθ) = (3/4)[( 1 - r00) + (3r00 - 1)cos²θ]
            //   → r^04_00 = (1 - 3<cos²θ>_w) / (1 - 9<cos²θ>_w)   ... only if isotropic baseline
            //   Better: extract via weighted cos²θ moment with proper ε-per-event
            //   R = r^04_00 / (ε_bar × (1 - r^04_00))
            //   where ε_bar = Σ(ε_i × w_i) / Σ w_i   (event-averaged ε in this bin)
            double sumEpsW    = 0.0;   // Σ ε_i × w_i   (for ε̄_bin)
            int    n          = 0;
        };
        // ── Per-thread accumulator arrays ─────────────────────────────
        // Each thread gets its own copy → no false sharing, no atomics.
        // We merge them serially after the parallel section.
#ifdef _OPENMP
        int nThreads = omp_get_max_threads();
#else
        int nThreads = 1;
#endif
        // [thread][bin]
        std::vector<std::vector<Acc>> tacc(nThreads, std::vector<Acc>(nBins));

        const int nEv = static_cast<int>(cache.events.size());

#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
        for (int ei = 0; ei < nEv; ++ei) {
#ifdef _OPENMP
            int tid = omp_get_thread_num();
#else
            int tid = 0;
#endif
            const auto& ev = cache.events[ei];
            int idx = bins_.flatIndex(ev.iQ2, ev.it, ev.ixB);

            // Build helicity matrix for this event's kinematics + current params
            Wkernels::Mat4 u{};
            SDMEModel::fillMatrix(params, ev.xB, ev.Q2, ev.t, u);

            // Compute σ_T and σ_L at this kinematics
            double sigmaT, sigmaL;
            SDMEModel::getSigmaLT(params, ev.xB, ev.Q2, ev.t, sigmaT, sigmaL);

            // Full σ₃ with amplitude dependence restored
            double sigma3 = ev.sigma3_prefactor * (sigmaT + ev.eps * sigmaL);

            // Angular kernels
            auto WUU = Wkernels::UU(u, ev.eps, ev.prodPhi, ev.decayPhi);
            auto WLU = Wkernels::LU(u, ev.eps, ev.prodPhi, ev.decayPhi);

            double cosT  = std::cos(ev.decayTheta);
            double sinT  = std::sin(ev.decayTheta);
            double cos2T = cosT * cosT;
            double sin2T = sinT * sinT;

            double W_UU = cos2T             * WUU.LL
                        + M_SQRT2*cosT*sinT * WUU.LT
                        + sin2T             * WUU.TT;

            double W_LU = cos2T             * WLU.LL
                        + M_SQRT2*cosT*sinT * WLU.LT
                        + sin2T             * WLU.TT;

            double inv4pi2 = 1.0 / (4.0 * M_PI * M_PI);
            double w = sigma3 * W_UU * inv4pi2;

            Acc& a = tacc[tid][idx];
            a.sumWUU    += w;
            a.sumWLUPol += ev.beamPol * sigma3 * W_LU * inv4pi2;
            a.sumLL     += cos2T * w;
            a.sumTT     += sin2T * w;
            a.sumLT     += std::sin(2.0 * ev.decayTheta) * std::cos(ev.decayPhi) * w;
            a.sumSigmaL += ev.sigma3_prefactor * ev.eps * sigmaL * W_UU * inv4pi2;
            a.sumSigmaT += ev.sigma3_prefactor           * sigmaT * W_UU * inv4pi2;
            a.sumEpsW   += ev.eps * w;
            a.n         += 1;
        }

        // ── Merge per-thread accumulators ──────────────────────────────
        std::vector<Acc> acc(nBins);
        for (int t = 0; t < nThreads; ++t)
            for (int b = 0; b < nBins; ++b) {
                acc[b].sumWUU    += tacc[t][b].sumWUU;
                acc[b].sumWLUPol += tacc[t][b].sumWLUPol;
                acc[b].sumLL     += tacc[t][b].sumLL;
                acc[b].sumTT     += tacc[t][b].sumTT;
                acc[b].sumLT     += tacc[t][b].sumLT;
                acc[b].sumSigmaL += tacc[t][b].sumSigmaL;
                acc[b].sumSigmaT += tacc[t][b].sumSigmaT;
                acc[b].sumEpsW   += tacc[t][b].sumEpsW;
                acc[b].n         += tacc[t][b].n;
            }

        // ── Convert accumulators → physical observables ────────────────
        for (int idx = 0; idx < nBins; ++idx) {
            const Acc& a = acc[idx];
            BinResult& r = results[idx];

            r.nEvents = a.n;
            r.sumW    = a.sumWUU;

            if (a.n == 0 || a.sumWUU == 0.0) continue;

            // dσ/dt: proportional to sum of weights
            // (absolute normalisation requires luminosity + bin width;
            //  for the chi2 fit we compare shapes, so normalise by nEvents)
            r.dsigma_dt  = a.sumWUU  / a.n;
            r.dsigmaL_dt = a.sumSigmaL / a.n;
            r.dsigmaT_dt = a.sumSigmaT / a.n;

            // BSA
            r.A_LU     = a.sumWLUPol / a.sumWUU;
            r.A_LU_err = (a.n > 0) ? 1.0 / std::sqrt(a.n) : 1.0;   // stat only

            // ── Angular moments (normalised to 1) ─────────────────────
            r.f_LL = a.sumLL / a.sumWUU;
            r.f_TT = a.sumTT / a.sumWUU;
            r.M_LT = a.sumLT / a.sumWUU;

            // ── R = σ_L/σ_T  via the angular distribution method ──────
            //
            // This is NOT a Rosenbluth separation (which requires multiple
            // beam energies). At a single energy, R is extracted from the
            // shape of the K+ polar angle distribution in the phi helicity
            // frame, using the SDME element r^04_00.
            //
            // The cos²θ distribution is:
            //   W(cosθ) = (3/4π)[(1 - r^04_00) + (3·r^04_00 - 1)·cos²θ]
            //
            // Taking the weighted average:
            //   <cos²θ>_w = (1/5)(1 + 2·r^04_00)    [from integrating W]
            //   → r^04_00 = (5·<cos²θ>_w - 1) / 2
            //
            // Then with the event-averaged ε̄ (properly computed per bin):
            //   R = r^04_00 / (ε̄ · (1 − r^04_00))
            //
            // Note: this assumes SCHC holds (r^04_00 → ε̄R/(1+ε̄R)).
            // The M_LT moment above tells you how badly SCHC is violated —
            // nonzero M_LT means this R extraction has a systematic bias.
            {
                // event-averaged ε̄ in this bin (properly weighted)
                double eps_bar = a.sumEpsW / a.sumWUU;

                // r^04_00 from the cos²θ moment
                double cos2_avg = r.f_LL;   // f_LL = <cos²θ>_w already
                double r04_00   = (5.0 * cos2_avg - 1.0) / 2.0;

                // Guard: r^04_00 must be in [0, 1] physically
                r04_00 = std::max(0.0, std::min(r04_00, 0.999));

                r.eps_bar = eps_bar;
                r.r04_00  = r04_00;

                if (eps_bar > 1e-4 && (1.0 - r04_00) > 1e-6)
                    r.R = r04_00 / (eps_bar * (1.0 - r04_00));
                else
                    r.R = 0.0;
            }
        }

        return results;
    }

    // ── Print a summary table ─────────────────────────────────────────────
    static void printResults(const std::vector<BinResult>& results,
                             const BinDef& bins)
    {
        printf("\n%-4s %-4s  %-7s %-7s  %-9s  %-8s  %-6s %-6s  %-7s %-7s  %-6s\n",
               "iQ2","it","Q2","t'","dsigma/dt","A_LU",
               "r04_00","eps_bar","R=sL/sT","M_LT","nEvt");
        printf("%s\n", std::string(90,'-').c_str());
        for (const auto& r : results) {
            if (r.nEvents == 0) continue;
            printf("%-4d %-4d  %-7.3f %-7.3f  %-9.3e  %+8.4f  %-6.4f %-7.4f  %-7.3f  %+7.4f  %-6d\n",
                   r.iQ2, r.it,
                   r.Q2_center, r.t_center,
                   r.dsigma_dt,
                   r.A_LU,
                   r.r04_00, r.eps_bar, r.R,
                   r.M_LT,
                   r.nEvents);
        }
        printf("\n");
        // SCHC violation warning
        bool schc_ok = true;
        for (const auto& r : results) {
            if (r.nEvents == 0) continue;
            if (std::abs(r.M_LT) > 0.05) { schc_ok = false; break; }
        }
        if (!schc_ok)
            printf("  NOTE: |M_LT| > 0.05 in some bins → SCHC violation; R extraction is biased.\n"
                   "        The M_LT moment constrains the LT interference amplitude A_LT.\n\n");
    }

private:
    const BinDef& bins_;
};
