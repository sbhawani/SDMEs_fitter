/**
 * BinDef.h
 *
 * Bin definitions for the phi electroproduction analysis.
 *
 * Bins match the CLAS12 RGA phi analysis (B. Singh, DVEP group, March 2026):
 *   - 4 Q² bins (equal-statistics strategy)
 *   - 8 t'  bins
 *   - 1 integrated xB bin (or 3 xB bins for differential studies)
 *
 * Usage:
 *   BinDef bins;
 *   int iQ2 = bins.findQ2(2.1);    // returns 1  (second bin)
 *   int it  = bins.findT(0.55);    // returns 1  (second t' bin)
 *   bool ok = bins.inRange(2.1, 0.55, 0.3);
 */

#pragma once
#include <vector>
#include <string>
#include <cmath>
#include <cassert>
#include <stdexcept>

struct BinEdges {
    std::vector<double> edges;   // N+1 edges for N bins

    int nBins() const { return static_cast<int>(edges.size()) - 1; }

    // Returns bin index [0, nBins-1], or -1 if out of range
    int find(double val) const {
        for (int i = 0; i < nBins(); ++i)
            if (val >= edges[i] && val < edges[i+1]) return i;
        return -1;
    }

    double center(int i) const {
        assert(i >= 0 && i < nBins());
        return 0.5 * (edges[i] + edges[i+1]);
    }

    double lo(int i) const { return edges[i]; }
    double hi(int i) const { return edges[i+1]; }
    double width(int i) const { return edges[i+1] - edges[i]; }
};

// ── Full bin configuration ────────────────────────────────────────────────────
struct BinDef {

    BinEdges Q2bins;  // Q² [GeV²]
    BinEdges tbins;   // t' = |t| - |t_min|  [GeV²]
    BinEdges xBbins;  // xB (optional; default = 1 integrated bin)

    // ------------------------------------------------------------------
    // Default constructor: CLAS12 phi analysis bins
    // ------------------------------------------------------------------
    BinDef() {
        // 4 Q² bins (equal-statistics, from slide 10 of the analysis)
        Q2bins.edges = { 1.00, 1.35, 1.81, 2.63, 8.00 };

        // 8 t' bins
        tbins.edges  = { 0.00, 0.45, 0.60, 0.75, 0.90,
                         1.05, 1.27, 2.02, 4.50 };

        // 1 integrated xB bin by default
        xBbins.edges = { 0.08, 0.68 };
    }

    // ------------------------------------------------------------------
    // Convenience: enable 3 xB bins for differential studies
    // ------------------------------------------------------------------
    void enableXBbins() {
        xBbins.edges = { 0.08, 0.20, 0.35, 0.68 };
    }

    int nQ2()  const { return Q2bins.nBins(); }
    int nT()   const { return tbins.nBins();  }
    int nXB()  const { return xBbins.nBins(); }
    int nTotal() const { return nQ2() * nT() * nXB(); }

    // Returns -1 if out of range
    int findQ2(double Q2)  const { return Q2bins.find(Q2); }
    int findT (double tp)  const { return tbins.find(tp);  }
    int findXB(double xB)  const { return xBbins.find(xB); }

    // Flat bin index:  (iQ2, it, ixB) → single integer
    int flatIndex(int iQ2, int it, int ixB) const {
        return iQ2 * nT() * nXB() + it * nXB() + ixB;
    }

    // Reverse flat index → (iQ2, it, ixB)
    void unflatten(int idx, int& iQ2, int& it, int& ixB) const {
        ixB  = idx % nXB();
        it   = (idx / nXB()) % nT();
        iQ2  = idx / (nXB() * nT());
    }

    // Is this event in the analysis range?
    bool inRange(double Q2, double tp, double xB) const {
        return findQ2(Q2) >= 0 && findT(tp) >= 0 && findXB(xB) >= 0;
    }

    // Bin label for printing / ROOT histogram naming
    std::string labelQ2(int i) const {
        char buf[64];
        snprintf(buf, sizeof(buf), "Q2_%.2fto%.2f",
                 Q2bins.lo(i), Q2bins.hi(i));
        return buf;
    }
    std::string labelT(int i) const {
        char buf[64];
        snprintf(buf, sizeof(buf), "tp_%.2fto%.2f",
                 tbins.lo(i), tbins.hi(i));
        return buf;
    }
};
