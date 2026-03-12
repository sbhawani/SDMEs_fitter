/**
 * PhysicsParams.h
 *
 * Tunable amplitude parameters for exclusive phi electroproduction.
 *
 * The full cross section decomposes into three helicity sectors:
 *
 *   dσ/dt dΩ  =  σ3(xB,Q²,ε) × (1/4π²) ×
 *                [ W_UU^LL(cosθ) + W_UU^LT(cosθ,φ_KK) + W_UU^TT(cosθ,φ_KK) ]
 *              + Pl × [ W_LU^LL + W_LU^LT + W_LU^TT ]    ← BSA terms
 *
 * Each sector is controlled by its own amplitude, t-slope and xB power:
 *
 *  LL:  Longitudinal photon → Longitudinal phi (σ_L)
 *       Drives:  cos²θ term in W_UU,  sin(φ) in BSA via Im_LL
 *
 *  TT:  Transverse photon → Transverse phi (σ_T)
 *       Drives:  sin²θ term in W_UU,  sin(φ) in BSA via Im_TT
 *
 *  LT:  L/T interference
 *       Drives:  cos(φ_KK), cos(φ+φ_KK) modulations in W_UU and BSA
 *
 * The helicity matrix elements are set in SDMEModel.h using these params.
 *
 * Parameters can be loaded from a simple key=value config file.
 */

#pragma once
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <stdexcept>
#include <cmath>

struct PhysicsParams {

    // -----------------------------------------------------------------------
    // Longitudinal sector (LL)
    // u[0][0][0][0]  →  real part → dσL/dt shape
    // u[0][0][0][+]  →  imaginary part → BSA sin(φ) from σ_LT'
    // -----------------------------------------------------------------------
    double A_L     = 30.0;   // overall L amplitude      [arb → nb/GeV²]
    double b_L     = 4.25;   // t-slope for L            [GeV⁻²]
    double n_L     = 1.0;    // xB power for L:  ~ xB^n_L
    double Im_LL   = 30.0;   // imag. part of u[0][0][0][+], controls sin(φ) BSA-LL
                              //  (enters as Im_LL * sqrt(-t)/Q)

    // -----------------------------------------------------------------------
    // Transverse sector (TT)
    // u[+][+][+][+]  →  real part → dσT/dt shape
    // u[+][+][0][+]  →  imaginary part → BSA sin(φ) from TT term
    // -----------------------------------------------------------------------
    double A_T     = 10.0;   // overall T amplitude      [arb → nb/GeV²]
    double b_T     = 4.25;   // t-slope for T            [GeV⁻²]
    double Delta_b =  1.5;   // additional slope for helicity-flip TT term [GeV⁻²]
    double n_T     = 1.0;    // xB power: (xB/(1-xB))^n_T
    double Im_TT   = -5.0;   // imag. part of u[+][+][0][+], controls sin(φ) BSA-TT
                              //  (enters as Im_TT * sqrt(-t)/Q)

    // -----------------------------------------------------------------------
    // LT interference sector
    // u[0][+][+][+]  →  drives cos(φ_KK) and sin(φ_KK) modulations
    // Both real and imaginary parts contribute to cross-section and BSA
    // -----------------------------------------------------------------------
    double A_LT    = 25.0;   // real part of u[0][+][+][+]
    double A_LTi   = 25.0;   // imaginary part of u[0][+][+][+]
    double b_LT    = 4.25;   // t-slope for LT term      [GeV⁻²]
    double n_LT    = 1.0;    // xB power for LT term

    // -----------------------------------------------------------------------
    // Beam and kinematics
    // -----------------------------------------------------------------------
    double beamEnergy = 10.6;   // GeV
    double Q2min      =  1.0;
    double Q2max      =  8.5;
    double xBmin      =  0.08;
    double xBmax      =  0.68;
    int    nEvents    =  50000;
    int    nScan      = 150000;

    // -----------------------------------------------------------------------
    // Load parameters from a key=value config file.
    // Lines starting with '#' are treated as comments.
    // -----------------------------------------------------------------------
    void load(const std::string& filename) {
        std::ifstream f(filename);
        if (!f.is_open())
            throw std::runtime_error("Cannot open config file: " + filename);

        std::string line;
        while (std::getline(f, line)) {
            // strip comments and empty lines
            auto pos = line.find('#');
            if (pos != std::string::npos) line = line.substr(0, pos);
            if (line.find_first_not_of(" \t\r\n") == std::string::npos) continue;

            std::istringstream ss(line);
            std::string key; double val;
            char eq;
            ss >> key >> eq >> val;
            if (eq != '=') continue;

            if      (key == "A_L")        A_L        = val;
            else if (key == "b_L")        b_L        = val;
            else if (key == "n_L")        n_L        = val;
            else if (key == "Im_LL")      Im_LL      = val;
            else if (key == "A_T")        A_T        = val;
            else if (key == "b_T")        b_T        = val;
            else if (key == "Delta_b")    Delta_b    = val;
            else if (key == "n_T")        n_T        = val;
            else if (key == "Im_TT")      Im_TT      = val;
            else if (key == "A_LT")       A_LT       = val;
            else if (key == "A_LTi")      A_LTi      = val;
            else if (key == "b_LT")       b_LT       = val;
            else if (key == "n_LT")       n_LT       = val;
            else if (key == "beamEnergy") beamEnergy = val;
            else if (key == "Q2min")      Q2min      = val;
            else if (key == "Q2max")      Q2max      = val;
            else if (key == "xBmin")      xBmin      = val;
            else if (key == "xBmax")      xBmax      = val;
            else if (key == "nEvents")    nEvents    = static_cast<int>(val);
            else if (key == "nScan")      nScan      = static_cast<int>(val);
            else
                std::cerr << "[PhysicsParams] WARNING: unknown key '" << key << "'\n";
        }
    }

    // -----------------------------------------------------------------------
    // Print current parameter values to stdout
    // -----------------------------------------------------------------------
    void print() const {
        printf("\n========================================================\n");
        printf("  phi electroproduction SDME generator — Physics Params \n");
        printf("========================================================\n");
        printf("  Beam energy   : %.2f GeV\n", beamEnergy);
        printf("  Q² range      : [%.2f, %.2f] GeV²\n", Q2min, Q2max);
        printf("  xB  range     : [%.3f, %.3f]\n", xBmin, xBmax);
        printf("  nEvents       : %d\n", nEvents);
        printf("--------------------------------------------------------\n");
        printf("  LL  sector:\n");
        printf("    A_L   = %8.3f   (σ_L amplitude)\n", A_L);
        printf("    b_L   = %8.3f   (t-slope, GeV⁻²)\n", b_L);
        printf("    n_L   = %8.3f   (xB power)\n", n_L);
        printf("    Im_LL = %8.3f   (BSA sin(φ) from σ_LT')\n", Im_LL);
        printf("--------------------------------------------------------\n");
        printf("  TT  sector:\n");
        printf("    A_T     = %8.3f   (σ_T amplitude)\n", A_T);
        printf("    b_T     = %8.3f   (t-slope, GeV⁻²)\n", b_T);
        printf("    Delta_b = %8.3f   (extra helicity-flip slope)\n", Delta_b);
        printf("    n_T     = %8.3f   (xB power)\n", n_T);
        printf("    Im_TT   = %8.3f   (BSA sin(φ) from TT)\n", Im_TT);
        printf("--------------------------------------------------------\n");
        printf("  LT  interference:\n");
        printf("    A_LT  = %8.3f   (real part)\n", A_LT);
        printf("    A_LTi = %8.3f   (imag part)\n", A_LTi);
        printf("    b_LT  = %8.3f   (t-slope, GeV⁻²)\n", b_LT);
        printf("    n_LT  = %8.3f   (xB power)\n", n_LT);
        printf("========================================================\n\n");
    }

    // -----------------------------------------------------------------------
    // Save current params to file (useful for logging a run's config)
    // -----------------------------------------------------------------------
    void save(const std::string& filename) const {
        std::ofstream f(filename);
        f << "# phi SDME generator config — auto-saved\n\n";
        f << "# Beam / kinematic range\n";
        f << "beamEnergy = " << beamEnergy << "\n";
        f << "Q2min      = " << Q2min      << "\n";
        f << "Q2max      = " << Q2max      << "\n";
        f << "xBmin      = " << xBmin      << "\n";
        f << "xBmax      = " << xBmax      << "\n";
        f << "nEvents    = " << nEvents    << "\n";
        f << "nScan      = " << nScan      << "\n";
        f << "\n# LL sector\n";
        f << "A_L   = " << A_L   << "\n";
        f << "b_L   = " << b_L   << "\n";
        f << "n_L   = " << n_L   << "\n";
        f << "Im_LL = " << Im_LL << "\n";
        f << "\n# TT sector\n";
        f << "A_T     = " << A_T     << "\n";
        f << "b_T     = " << b_T     << "\n";
        f << "Delta_b = " << Delta_b << "\n";
        f << "n_T     = " << n_T     << "\n";
        f << "Im_TT   = " << Im_TT   << "\n";
        f << "\n# LT interference\n";
        f << "A_LT  = " << A_LT  << "\n";
        f << "A_LTi = " << A_LTi << "\n";
        f << "b_LT  = " << b_LT  << "\n";
        f << "n_LT  = " << n_LT  << "\n";
    }
};
