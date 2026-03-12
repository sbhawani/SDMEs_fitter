/**
 * DataSet.h
 *
 * Experimental data points for the chi2 fit.
 *
 * Three types of observables are stored:
 *   1.  Cross-section points:  dσ/dt  in (Q², t') bins
 *   2.  BSA points:            A_LU   in (Q², t') bins
 *   3.  Angular moment points: f_LL, f_TT, M_LT in (Q², t') bins
 *
 * Data can be loaded from simple text files (one line per bin) or filled
 * programmatically (e.g. from a ROOT file via a TH1 reader).
 *
 * File format (space-separated, lines starting with # ignored):
 *   xsec:  iQ2  it  ixB  value  error
 *   bsa:   iQ2  it  ixB  value  error
 *   mom:   iQ2  it  ixB  fLL  fLL_err  fTT  fTT_err  MLT  MLT_err
 *
 * Example data file:
 *   # iQ2 it ixB  A_LU  stat_err
 *   0  0  0   0.102   0.045
 *   0  1  0   0.085   0.038
 *   ...
 */

#pragma once
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <cmath>

struct DataPoint {
    int iQ2 = -1, it = -1, ixB = -1;

    // Bin centres for labelling
    double Q2_center = 0.0;
    double t_center  = 0.0;
    double xB_center = 0.0;

    double value = 0.0;
    double error = 1.0;   // stat + syst combined

    // Weight in chi2:  contribution = ((value - MC) / error)²  × weight
    double chi2weight = 1.0;
};

struct MomentPoint {
    int iQ2 = -1, it = -1, ixB = -1;

    double fLL = 0.0,     fLL_err = 1.0;
    double fTT = 0.0,     fTT_err = 1.0;
    double M_LT = 0.0,    M_LT_err = 1.0;
};

struct DataSet {
    std::vector<DataPoint>  xsec;    // dσ/dt points
    std::vector<DataPoint>  bsa;     // A_LU points
    std::vector<MomentPoint> moments; // angular moment points

    // ── Load cross-section data from file ───────────────────────────────
    void loadXsec(const std::string& filename) {
        loadGeneric(filename, xsec);
        std::cout << "[DataSet] Loaded " << xsec.size()
                  << " xsec points from " << filename << "\n";
    }

    void loadBSA(const std::string& filename) {
        loadGeneric(filename, bsa);
        std::cout << "[DataSet] Loaded " << bsa.size()
                  << " BSA points from " << filename << "\n";
    }

    void loadMoments(const std::string& filename) {
        std::ifstream f(filename);
        if (!f.is_open()) {
            std::cerr << "[DataSet] Cannot open: " << filename << "\n";
            return;
        }
        std::string line;
        while (std::getline(f, line)) {
            auto pos = line.find('#');
            if (pos != std::string::npos) line = line.substr(0, pos);
            if (line.find_first_not_of(" \t") == std::string::npos) continue;
            MomentPoint mp;
            std::istringstream ss(line);
            ss >> mp.iQ2 >> mp.it >> mp.ixB
               >> mp.fLL >> mp.fLL_err
               >> mp.fTT >> mp.fTT_err
               >> mp.M_LT >> mp.M_LT_err;
            if (!ss.fail()) moments.push_back(mp);
        }
        std::cout << "[DataSet] Loaded " << moments.size()
                  << " moment points from " << filename << "\n";
    }

    // ── Fill from MC pseudo-data (for closure tests) ─────────────────────
    /**
     * Fill data set from a BinResult vector (e.g. generated with known params).
     * Adds Poisson-like statistical fluctuations scaled by fractional_stat_err.
     */
    void fillFromMC(const std::vector<struct BinResult>& mc_results,
                    double fractional_stat_err = 0.05,
                    unsigned int seed = 42)
    {
        std::mt19937 rng(seed);
        std::normal_distribution<double> gauss(0.0, 1.0);

        xsec.clear(); bsa.clear(); moments.clear();

        for (const auto& r : mc_results) {
            if (r.nEvents < 10) continue;

            double stat = fractional_stat_err;

            // Cross-section point
            DataPoint dp;
            dp.iQ2 = r.iQ2; dp.it = r.it; dp.ixB = r.ixB;
            dp.Q2_center = r.Q2_center;
            dp.t_center  = r.t_center;
            dp.xB_center = r.xB_center;
            dp.value = r.dsigma_dt * (1.0 + stat * gauss(rng));
            dp.error = r.dsigma_dt * stat;
            if (dp.error > 0) xsec.push_back(dp);

            // BSA point
            DataPoint db = dp;
            db.value = r.A_LU + r.A_LU_err * gauss(rng);
            db.error = r.A_LU_err;
            if (db.error > 0) bsa.push_back(db);

            // Moment point
            MomentPoint mp;
            mp.iQ2 = r.iQ2; mp.it = r.it; mp.ixB = r.ixB;
            mp.fLL     = r.f_LL  + stat * gauss(rng);
            mp.fLL_err = stat;
            mp.fTT     = r.f_TT  + stat * gauss(rng);
            mp.fTT_err = stat;
            mp.M_LT    = r.M_LT  + stat * gauss(rng);
            mp.M_LT_err = stat;
            moments.push_back(mp);
        }

        std::cout << "[DataSet] Filled from MC: "
                  << xsec.size() << " xsec, "
                  << bsa.size()  << " BSA, "
                  << moments.size() << " moment points\n";
    }

    bool empty() const { return xsec.empty() && bsa.empty() && moments.empty(); }

private:
    void loadGeneric(const std::string& filename,
                     std::vector<DataPoint>& pts)
    {
        std::ifstream f(filename);
        if (!f.is_open()) {
            std::cerr << "[DataSet] Cannot open: " << filename << "\n";
            return;
        }
        std::string line;
        while (std::getline(f, line)) {
            auto pos = line.find('#');
            if (pos != std::string::npos) line = line.substr(0, pos);
            if (line.find_first_not_of(" \t") == std::string::npos) continue;
            DataPoint dp;
            std::istringstream ss(line);
            ss >> dp.iQ2 >> dp.it >> dp.ixB >> dp.value >> dp.error;
            if (!ss.fail() && dp.error > 0) pts.push_back(dp);
        }
    }
};

// Include BinResult here to avoid circular dependency
#include "AnalysisModule.h"
