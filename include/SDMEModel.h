/**
 * SDMEModel.h
 *
 * Translates PhysicsParams → Wkernels::Mat4 helicity amplitude matrix u[ν'][ν][μ'][μ]
 *
 * The indices follow Wkernels convention:  '+' → 0, '0' → 1, '-' → 2
 *
 * The mapping between matrix elements and physics observables is:
 *
 *  ── LL sector ──────────────────────────────────────────────────────────────
 *  u[0][0][0][0]  real part   →  σ_L (longitudinal cross section)
 *                               contributes to W_UU^LL via real(U('0','0','0','0'))
 *
 *  u[0][0][0][+]  imag part   →  σ_LT' (appears in BSA)
 *                               Im(U('0','0','0','+')) × sin(φ) → W_LU^LL
 *
 *  ── TT sector ──────────────────────────────────────────────────────────────
 *  u[+][+][+][+]  real part   →  σ_T (transverse cross section)
 *                               Re(U('+','+','+','+')) → W_UU^TT leading term
 *
 *  u[+][+][0][+]  imag part   →  helicity-flip T interference → BSA
 *                               Im(U('+','+','0','+')) × sin(φ) → W_LU^TT
 *
 *  ── LT sector ──────────────────────────────────────────────────────────────
 *  u[0][+][+][+]  real + imag →  L/T interference
 *                               Both parts enter W_UU^LT (cos(φ_KK) modulation)
 *                               and W_LU^LT (sin(φ_KK) BSA contribution)
 *
 * The kinematic pre-factors sqrt(-t)/Q suppress helicity-flip terms at small t,
 * which is physically motivated: flip amplitudes vanish as t → t_min.
 *
 * Reference: Goloskokov & Kroll (Eur.Phys.J.C 53 (2008) 367);
 *            Diehl, JHEP09(2007)064 (ρ production formalism extended to φ)
 */

#pragma once
#include "w_kernels.hpp"
#include "PhysicsParams.h"
#include <cmath>

namespace SDMEModel {

/**
 * Fill the 4-index helicity amplitude matrix from PhysicsParams.
 *
 * @param p    physics parameter set
 * @param xB   Bjorken x at this phase-space point
 * @param Q2   photon virtuality  [GeV²]
 * @param t    Mandelstam t       [GeV²]  (negative in spacelike region)
 * @param u    output matrix to fill
 */
inline void fillMatrix(const PhysicsParams& p,
                       double xB, double Q2, double t,
                       Wkernels::Mat4& u)
{
    // Zero out the matrix first
    for (auto& a : u) for (auto& b : a) for (auto& c : b) for (auto& d : c)
        d = {0.0, 0.0};

    // Guard: t is negative (spacelike), use |t| where appropriate
    double mt   = std::max(-t, 0.0);          // |t|
    double sqt  = std::sqrt(mt);              // sqrt(|t|)
    double Q    = std::sqrt(std::max(Q2, 1e-6));

    // xB-dependent prefactors
    double x1       = 1.0 - xB;
    double xdepl    = std::pow(xB, p.n_L);         // for LL
    double xdept    = std::pow(xB / x1, p.n_T);    // for TT (matches original code)
    double xdeplt   = std::pow(xB, p.n_LT);        // for LT

    // t-dependent envelopes
    double tL  = std::exp(p.b_L  * t);              // t < 0, so this decays
    double tT  = std::exp((p.b_T + p.Delta_b) * t); // extra slope for T flip
    double tLT = std::exp(p.b_LT * t);

    // ── LL sector ─────────────────────────────────────────────────────────
    //
    // u[0][0][0][0]: real part → dσ_L/dt shape
    //   enters W_UU^LL as: eps * Re(U('0','0','0','0'))
    u[Wkernels::h('0')][Wkernels::h('0')]
     [Wkernels::h('0')][Wkernels::h('0')]
        = { p.A_L * tL * xdepl, 0.0 };

    // u[0][0][0][+]: imag part → drives sin(φ) in BSA (W_LU^LL)
    //   enters as: -2 sin(φ) * sqrt(ε(1-ε)) * Im(U('0','0','0','+'))
    //   pre-factor sqrt(|t|)/Q ensures helicity-flip suppression at t_min
    u[Wkernels::h('0')][Wkernels::h('0')]
     [Wkernels::h('0')][Wkernels::h('+')]
        = { 0.0, p.Im_LL * tL * xdepl * sqt / Q };

    // ── TT sector ─────────────────────────────────────────────────────────
    //
    // u[+][+][+][+]: real part → dσ_T/dt shape
    //   enters W_UU^TT as: 0.5*Re(U('+','+','+','+')) + 0.5*Re(U('-','-','+','+'))
    //   We set only the ++ term; -- is related by symmetry (both set equal below)
    u[Wkernels::h('+')][Wkernels::h('+')]
     [Wkernels::h('+')][Wkernels::h('+')]
        = { p.A_T * tT * xdept, 0.0 };

    // Mirror term: u[-][-][+][+] (needed for W_UU^TT completeness)
    // In the pure helicity-non-flip limit these are equal
    u[Wkernels::h('-')][Wkernels::h('-')]
     [Wkernels::h('+')][Wkernels::h('+')]
        = { p.A_T * tT * xdept, 0.0 };

    // u[+][+][0][+]: imag part → helicity-flip T → BSA sin(φ) from W_LU^TT
    //   enters as: -sin(φ)*sqrt(ε(1-ε)) * Im(U('+','+','0','+') + U('-','-','0','+'))
    u[Wkernels::h('+')][Wkernels::h('+')]
     [Wkernels::h('0')][Wkernels::h('+')]
        = { 0.0, p.Im_TT * tT * xdept * sqt / Q };

    // ── LT interference sector ────────────────────────────────────────────
    //
    // u[0][+][+][+]: both real and imaginary parts
    //   Real part → cos(φ_KK) modulation in W_UU^LT
    //               cos(φ+φ_KK) modulation after angle combination
    //   Imag part → sin(φ_KK) in W_LU^LT (BSA)
    //   pre-factor sqrt(|t|)/Q from helicity-flip in the L→T amplitude
    u[Wkernels::h('0')][Wkernels::h('+')]
     [Wkernels::h('+')][Wkernels::h('+')]
        = { p.A_LT  * tLT * xdeplt * sqt / Q,
            p.A_LTi * tLT * xdeplt * sqt / Q };
}

/**
 * Compute dσ_L/dt and dσ_T/dt at a given kinematic point.
 * These feed into the 3-fold prefactor σ3(xB,Q²,ε).
 */
inline void getSigmaLT(const PhysicsParams& p,
                       double xB, double Q2, double t,
                       double& sigmaT, double& sigmaL)
{
    double mt    = std::max(-t, 0.0);
    double x1    = 1.0 - xB;
    double xdepl = std::pow(xB, p.n_L);
    double xdept = std::pow(xB / x1, p.n_T);
    double tL    = std::exp(p.b_L * t);
    double tT    = std::exp((p.b_T + p.Delta_b) * t);

    sigmaL = p.A_L * tL * xdepl;
    sigmaT = p.A_T * tT * xdept;
}

} // namespace SDMEModel
