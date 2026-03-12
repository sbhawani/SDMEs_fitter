/**
 * generator.h
 *
 * Event generator for exclusive phi electroproduction.
 * 
 * This version replaces the hard-coded helicity matrix in the original
 * loadPhysicsInputs() with the SDMEModel which is driven by PhysicsParams,
 * allowing LL / LT / TT amplitude parameters to be varied externally.
 *
 * Output: LUND-format events  (same as original phiproduce.10.6.C)
 *         Each event contains:  e_in, p_in, e_out, p_out, phi, K+, K-
 *
 * Usage:
 *   PhysicsParams par;
 *   par.load("config/default.cfg");          // or edit struct directly
 *
 *   sim::candidate cr(par.beamEnergy, 1.02, 0.49368, 0.49368);
 *   sim::SDMEGenerator gen(&cr, par);
 *   gen.setRange(par.Q2min, par.Q2max, par.xBmin, par.xBmax);
 *
 *   double maxW = gen.scan(par.nScan);
 *   for (int i = 0; i < par.nEvents; i++) {
 *       gen.generate(maxW);
 *       sim::event ev;
 *       cr.react.getEvent(ev);
 *       if (!ev.hasNaN()) ev.show();
 *   }
 */

#pragma once

#include "reaction.h"
#include "w_kernels.hpp"
#include "cross.h"
#include "PhysicsParams.h"
#include "SDMEModel.h"

namespace sim {

// ── physicsInput (unchanged from original) ──────────────────────────────────
struct physicsInput {
    Wkernels::Mat4 u, l, s;
    double dsigmaT_dt, dsigmaL_dt;   // [nb / GeV²]
    double Pl = 0.0, SL = 0.0, ST = 0.0;
};

// ── SDMEGenerator ─────────────────────────────────────────────────────────────
class SDMEGenerator {
private:
    static constexpr double Mp   = 0.93827;
    static constexpr double alem = 1.0 / 137.035999084;

    candidate*    cand;
    PhysicsParams params;
    TRandom3      rand;

    double q2_min, q2_max;
    double xb_min, xb_max;
    double E;

    int generated, accepted;

public:
    // -------------------------------------------------------------------------
    SDMEGenerator(candidate* c, const PhysicsParams& p)
        : cand(c), params(p)
    {
        E         = cand->react.E();
        generated = 0;
        accepted  = 0;
    }

    // -------------------------------------------------------------------------
    void setRange(double q2min, double q2max, double xbmin, double xbmax) {
        q2_min = q2min; q2_max = q2max;
        xb_min = xbmin; xb_max = xbmax;
    }

    // Allow updating params at runtime (e.g. for parameter scans / chi2 fits)
    void setParams(const PhysicsParams& p) { params = p; }
    const PhysicsParams& getParams() const { return params; }

    // Seed the internal RNG (use different seeds per thread for parallel generation)
    void setSeed(unsigned seed) { rand.SetSeed(seed); }

    // -------------------------------------------------------------------------
    // Kinematic validity: Q² > 1 GeV² and above threshold
    // -------------------------------------------------------------------------
    int is_valid(double q2, double xb) {
        double upper = 2.0 * Mp * E * xb / (1.0 + Mp * xb / E);
        double lower = 2.73 * xb / (1.0 - xb);
        return (q2 > 1.0 && q2 > lower && q2 < upper) ? 1 : 0;
    }

    // -------------------------------------------------------------------------
    // Generate flat (Q²,xB) then produce full reaction kinematics
    // -------------------------------------------------------------------------
    void generate() {
        double q2 = rand.Uniform(q2_min, q2_max);
        double xb = rand.Uniform(xb_min, xb_max);
        while (is_valid(q2, xb) == 0) {
            q2 = rand.Uniform(q2_min, q2_max);
            xb = rand.Uniform(xb_min, xb_max);
        }
        cand->react.generate(q2, xb);
    }

    // -------------------------------------------------------------------------
    // 3-fold virtual photoproduction cross section prefactor
    // -------------------------------------------------------------------------
    double dsigma_3fold(double xB, double Q2, double y,
                        double eps, double sigmaT, double sigmaL) {
        double pref = alem / (2.0 * M_PI);
        double kin  = (y * y) / (1.0 - eps) * (1.0 - xB) / xB / Q2;
        return pref * kin * (sigmaT + eps * sigmaL);
    }

    // -------------------------------------------------------------------------
    // 7-fold differential cross section (UU + beam-polarisation term)
    // -------------------------------------------------------------------------
    double dsigma_7fold(double WUU, double WLU,
                        double WUL, double WLL,
                        double WUT, double WLT,
                        double Pl, double SL, double ST,
                        double sigma3) {
        double S = WUU + Pl * WLU + SL * WUL + Pl * SL * WLL
                 + ST * WUT + Pl * ST * WLT;
        return sigma3 * S / (4.0 * M_PI * M_PI);
    }

    // -------------------------------------------------------------------------
    // Update the event weight using the current PhysicsParams
    // -------------------------------------------------------------------------
    void updateWeight() {
        double t  = (cand->react.p_out - cand->react.p_in).m2();
        double Q2 = cand->react.Q2();
        double xB = cand->react.xB();

        // Build helicity matrix from current params
        physicsInput ph{};
        SDMEModel::fillMatrix(params, xB, Q2, t, ph.u);
        SDMEModel::getSigmaLT(params, xB, Q2, t, ph.dsigmaT_dt, ph.dsigmaL_dt);

        double y   = cand->react.Y();
        double eps = cand->react.Eps();

        double sigma3 = dsigma_3fold(xB, Q2, y, eps,
                                     ph.dsigmaT_dt, ph.dsigmaL_dt);

        // Compute angular kernels
        auto WUU = Wkernels::UU(ph.u, eps,
                                cand->react.prodPhi,
                                cand->react.decayPhi);
        auto WLU = Wkernels::LU(ph.u, eps,
                                cand->react.prodPhi,
                                cand->react.decayPhi);

        double cosT = std::cos(cand->react.decayTheta);
        double sinT = std::sin(cand->react.decayTheta);

        // Full angular distribution:
        //   W = cos²θ × W^LL  +  √2 cosθ sinθ × W^LT  +  sin²θ × W^TT
        double W_UU = cosT * cosT       * WUU.LL
                    + std::sqrt(2.0) * cosT * sinT * WUU.LT
                    + sinT * sinT       * WUU.TT;

        double W_LU = cosT * cosT       * WLU.LL
                    + std::sqrt(2.0) * cosT * sinT * WLU.LT
                    + sinT * sinT       * WLU.TT;

        double pol = (cand->react.rpol < 0) ? -1.0 : 1.0;

        cand->weight = dsigma_7fold(W_UU, W_LU, 0, 0, 0, 0,
                                    pol, ph.SL, ph.ST, sigma3);
    }

    // -------------------------------------------------------------------------
    // Scan phase space to find the maximum weight (for rejection sampling)
    // -------------------------------------------------------------------------
    double scan(int count) {
        double maxWeight = 0.0;
        for (int i = 0; i < count; i++) {
            generate();
            updateWeight();
            if (cand->weight > maxWeight) maxWeight = cand->weight;
        }
        return maxWeight * 1.2;   // safety margin
    }

    // -------------------------------------------------------------------------
    // Generate one accepted event using hit-or-miss rejection sampling
    // -------------------------------------------------------------------------
    void generate(double maxWeight) {
        generate();
        updateWeight();
        generated++;

        double accept = rand.Uniform(0, maxWeight);
        while (accept > cand->weight) {
            generate();
            updateWeight();
            generated++;
            accept = rand.Uniform(0, maxWeight);
        }
        accepted++;
    }

    // -------------------------------------------------------------------------
    // Print generation efficiency
    // -------------------------------------------------------------------------
    void stats() const {
        printf("\nGenerator: %d generated, %d accepted, efficiency = %.4f\n\n",
               generated, accepted,
               (generated > 0) ? static_cast<double>(accepted) / generated : 0.0);
    }

    // -------------------------------------------------------------------------
    // Return the W_LU value for the current event (useful for BSA analysis)
    // -------------------------------------------------------------------------
    double getWLU() {
        double t  = (cand->react.p_out - cand->react.p_in).m2();
        double Q2 = cand->react.Q2();
        double xB = cand->react.xB();
        double eps = cand->react.Eps();

        physicsInput ph{};
        SDMEModel::fillMatrix(params, xB, Q2, t, ph.u);

        auto WLU = Wkernels::LU(ph.u, eps,
                                cand->react.prodPhi,
                                cand->react.decayPhi);

        double cosT = std::cos(cand->react.decayTheta);
        double sinT = std::sin(cand->react.decayTheta);

        return cosT * cosT       * WLU.LL
             + std::sqrt(2.0) * cosT * sinT * WLU.LT
             + sinT * sinT       * WLU.TT;
    }

    // -------------------------------------------------------------------------
    // Return the W_UU value for the current event (useful for cross-section fits)
    // -------------------------------------------------------------------------
    double getWUU() {
        double t  = (cand->react.p_out - cand->react.p_in).m2();
        double Q2 = cand->react.Q2();
        double xB = cand->react.xB();
        double eps = cand->react.Eps();

        physicsInput ph{};
        SDMEModel::fillMatrix(params, xB, Q2, t, ph.u);

        auto WUU = Wkernels::UU(ph.u, eps,
                                cand->react.prodPhi,
                                cand->react.decayPhi);

        double cosT = std::cos(cand->react.decayTheta);
        double sinT = std::sin(cand->react.decayTheta);

        return cosT * cosT       * WUU.LL
             + std::sqrt(2.0) * cosT * sinT * WUU.LT
             + sinT * sinT       * WUU.TT;
    }
};

} // namespace sim
