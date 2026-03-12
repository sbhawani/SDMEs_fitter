/**
 * EventData.h
 *
 * Data structures passed between pipeline stages:
 *
 *   TruthEvent  →  FastMC  →  RecoEvent
 *                                 ↓
 *                           CachedEvent  (stored in EventCache for reweighting)
 *
 * The reweighting approach avoids regenerating events at every Minuit step:
 *   1. Generate once → pass through FastMC → fill EventCache
 *   2. At each minimiser step, call SDMEModel::weight(params, ev) per event
 *   3. Sum weights per bin → MC prediction for cross-section and BSA
 *
 * Angles follow the Goloskokov-Kroll / w_kernels convention:
 *   prodPhi   = Φ     (azimuthal angle of phi production plane vs lepton plane)
 *   decayPhi  = φ_KK  (K+ azimuthal angle in phi helicity frame)
 *   decayTheta = θ_H  (K+ polar angle in phi helicity frame)
 */

#pragma once
#include <vector>
#include <cmath>

// ── Raw truth-level kinematics from the generator ────────────────────────────
struct TruthEvent {
    // DIS variables
    double Q2, xB, y, eps;
    double t;           // Mandelstam t  [GeV²]  (negative, spacelike)
    double tprime;      // t' = |t| - |t_min|    [GeV²]

    // Production and decay angles
    double prodPhi;     // Φ   lepton-hadron plane angle [rad]
    double decayTheta;  // θ_H                           [rad]
    double decayPhi;    // φ_KK                          [rad]

    // Beam helicity:  +1 or -1
    double beamPol;

    // 4-momenta (lab frame, GeV)  stored as (px,py,pz,E)
    double e_px,  e_py,  e_pz,  e_E;    // scattered electron
    double p_px,  p_py,  p_pz,  p_E;    // recoil proton
    double kp_px, kp_py, kp_pz, kp_E;  // K+
    double km_px, km_py, km_pz, km_E;  // K-
};

// ── FastMC output ─────────────────────────────────────────────────────────────
struct RecoEvent {
    bool accepted = false;   // false → event killed by acceptance

    // Reconstructed DIS variables (after smearing / fiducial cuts)
    double Q2, xB, t, tprime;
    double eps, y;

    // Reconstructed angles (in phi helicity frame after reco-level boost)
    double prodPhi;
    double decayTheta;
    double decayPhi;

    double beamPol;
};

// ── Cached event — stored once, reweighted at each minimiser step ─────────────
struct CachedEvent {
    // Kinematic bin indices (set once during caching)
    int iQ2  = -1;
    int it   = -1;
    int ixB  = -1;

    // Quantities needed to evaluate W_UU and W_LU without re-running full kinematics
    double Q2, xB, t, tprime;
    double eps, y;
    double prodPhi;
    double decayTheta;
    double decayPhi;
    double beamPol;

    // Pre-computed σ₃ prefactor (depends only on xB,Q²,y,ε — not on amplitude params)
    // This is expensive to recompute, so store it.
    double sigma3_prefactor;   // [nb / GeV²]  (without dsigmaT + eps*dsigmaL)

    // t_min for this (xB, Q²) — needed for tprime
    double t_min;
};

// ── Convenience: compute t_min for given (xB, Q²) ────────────────────────────
inline double compute_tmin(double xB, double Q2)
{
    static constexpr double Mp = 0.93827;
    static constexpr double Mv = 1.020;   // phi mass
    double W2    = Mp * Mp + Q2 * (1.0 / xB - 1.0);
    double W     = std::sqrt(std::max(W2, 0.0));
    double nu    = Q2 / (2.0 * Mp * xB);
    double qstar = std::sqrt(nu * nu + Q2);
    // Mandelstam t_min (forward scattering)
    double Ephi  = (W2 + Mv * Mv - Mp * Mp) / (2.0 * W);
    double pphi  = std::sqrt(std::max(Ephi * Ephi - Mv * Mv, 0.0));
    double Ep    = (W2 - Mv * Mv + Mp * Mp) / (2.0 * W);
    double pq    = std::sqrt(qstar * qstar + Q2);   // |q*| in CM
    double cost_min = 1.0;
    double tmin = Q2 - 2.0 * (nu * Ep - pq * pphi * cost_min);
    return tmin;  // negative value
}

// EventCache class is defined in EventCache.h
