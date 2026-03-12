/**
 * FastMCInterface.h
 *
 * Abstract interface for the FastMC acceptance/resolution model.
 *
 * The real FastMC (Glasgow NN+BDT, arxiv:2207.11254) takes truth-level
 * 4-momenta and returns reconstructed momenta after detector simulation.
 * Replace DummyFastMC with a wrapper around that model when available.
 *
 * Pipeline:
 *   TruthEvent → FastMCInterface::process() → RecoEvent
 *
 * The DummyFastMC implemented here:
 *   - Applies simple fiducial geometric acceptance (FD/CD angular cuts)
 *   - Smears momenta with Gaussian resolution (~0.5% for tracks)
 *   - Recomputes DIS variables and helicity-frame angles from smeared 4-momenta
 *   - This gives realistic acceptance × efficiency shape for development/testing
 *
 * To plug in the real FastMC:
 *   class RealFastMC : public FastMCInterface {
 *   public:
 *       bool process(const TruthEvent& truth, RecoEvent& reco) override {
 *           // call your NN+BDT model here
 *           // fill reco fields
 *           return reco.accepted;
 *       }
 *   };
 */

#pragma once
#include "EventData.h"
#include <cmath>
#include <random>

// ── Abstract interface ────────────────────────────────────────────────────────
class FastMCInterface {
public:
    virtual ~FastMCInterface() = default;
    virtual bool process(const TruthEvent& truth, RecoEvent& reco) = 0;
    // Returns true for GavFastMC — tells EventCache to use batch mode
    virtual bool isBatchFastMC() const { return false; }
};

// ── Dummy FastMC ──────────────────────────────────────────────────────────────
class DummyFastMC : public FastMCInterface {
public:
    // Momentum resolution (fractional sigma)
    double trackResolution  = 0.005;   // 0.5% — typical CLAS12 FD

    // Geometric acceptance (polar angle in lab frame)
    double eTheta_min  =  5.0 * M_PI / 180.0;   // e-  FD lower edge
    double eTheta_max  = 35.0 * M_PI / 180.0;   // e-  FD upper edge
    double pTheta_min  =  5.0 * M_PI / 180.0;   // p   FD lower
    double pTheta_max  = 65.0 * M_PI / 180.0;   // p   FD/CD upper (CD up to 125°)
    double kTheta_min  =  5.0 * M_PI / 180.0;   // K± FD lower
    double kTheta_max  = 35.0 * M_PI / 180.0;   // K± FD upper (no RICH → FD only)

    // Minimum momenta [GeV]
    double p_e_min   = 1.5;
    double p_p_min   = 0.3;
    double p_k_min   = 0.2;

    explicit DummyFastMC(unsigned int seed = 42) : rng_(seed) {}

    bool process(const TruthEvent& truth, RecoEvent& reco) override
    {
        reco.accepted = false;
        reco.beamPol  = truth.beamPol;

        // ── Step 1: smear 4-momenta ──────────────────────────────────────
        auto smear = [&](double px, double py, double pz, double E,
                         double& spx, double& spy, double& spz, double& sE)
        {
            double p  = std::sqrt(px*px + py*py + pz*pz);
            double sf = 1.0 + gauss() * trackResolution;
            spx = px * sf; spy = py * sf; spz = pz * sf;
            double m2 = E*E - p*p;
            sE  = std::sqrt(spx*spx + spy*spy + spz*spz + std::max(m2, 0.0));
        };

        double se_px, se_py, se_pz, se_E;
        double sp_px, sp_py, sp_pz, sp_E;
        double skp_px, skp_py, skp_pz, skp_E;
        double skm_px, skm_py, skm_pz, skm_E;

        smear(truth.e_px,  truth.e_py,  truth.e_pz,  truth.e_E,  se_px,  se_py,  se_pz,  se_E );
        smear(truth.p_px,  truth.p_py,  truth.p_pz,  truth.p_E,  sp_px,  sp_py,  sp_pz,  sp_E );
        smear(truth.kp_px, truth.kp_py, truth.kp_pz, truth.kp_E, skp_px, skp_py, skp_pz, skp_E);
        smear(truth.km_px, truth.km_py, truth.km_pz, truth.km_E, skm_px, skm_py, skm_pz, skm_E);

        // ── Step 2: geometric acceptance ────────────────────────────────
        if (!passAcceptance(se_px,  se_py,  se_pz,  p_e_min,  eTheta_min, eTheta_max)) return false;
        if (!passAcceptance(sp_px,  sp_py,  sp_pz,  p_p_min,  pTheta_min, pTheta_max)) return false;
        if (!passAcceptance(skp_px, skp_py, skp_pz, p_k_min,  kTheta_min, kTheta_max)) return false;
        if (!passAcceptance(skm_px, skm_py, skm_pz, p_k_min,  kTheta_min, kTheta_max)) return false;

        // ── Step 3: recompute DIS variables from smeared momenta ─────────
        static constexpr double Mp    = 0.93827;
        static constexpr double beamE = 10.6;

        // q = beam - e'
        double qx = -se_px, qy = -se_py, qz = beamE - se_pz;
        double qE =  beamE  - se_E;
        double Q2 = -(qE*qE - qx*qx - qy*qy - qz*qz);
        if (Q2 < 1.0) return false;

        double nu = qE;   // target at rest
        double xB = Q2 / (2.0 * Mp * nu);
        if (xB < 0.08 || xB > 0.68) return false;

        double y   = nu / beamE;
        double gg2 = 4.0 * xB * xB * Mp * Mp / Q2;
        double eps = (1.0 - y - 0.25*y*y*gg2) / (1.0 - y + 0.5*y*y + 0.25*y*y*gg2);
        if (eps < 0.0 || eps > 1.0) return false;

        // t = (p_out - p_in)²
        double dp_px = sp_px, dp_py = sp_py, dp_pz = sp_pz - 0.0;
        double dp_E  = sp_E  - Mp;
        double t = dp_E*dp_E - dp_px*dp_px - dp_py*dp_py - dp_pz*dp_pz;
        double tmin = compute_tmin(xB, Q2);
        double tprime = std::abs(t) - std::abs(tmin);
        if (tprime < 0.0 || tprime > 4.5) return false;

        // ── Step 4: recompute helicity frame angles ───────────────────────
        // (use truth angles directly for dummy — real FastMC would recompute from reco 4-vectors)
        reco.decayTheta = truth.decayTheta;
        reco.decayPhi   = truth.decayPhi;
        reco.prodPhi    = truth.prodPhi;

        // ── Step 5: fill reco event ──────────────────────────────────────
        reco.Q2      = Q2;
        reco.xB      = xB;
        reco.t       = t;
        reco.tprime  = tprime;
        reco.y       = y;
        reco.eps     = eps;
        reco.accepted = true;

        return true;
    }

private:
    std::mt19937 rng_;
    std::normal_distribution<double> gauss_{0.0, 1.0};

    double gauss() { return gauss_(rng_); }

    bool passAcceptance(double px, double py, double pz,
                        double pmin, double thmin, double thmax)
    {
        double p   = std::sqrt(px*px + py*py + pz*pz);
        if (p < pmin) return false;
        double th  = std::acos(pz / p);
        return (th >= thmin && th <= thmax);
    }
};
