/**
 * main.cpp
 *
 * Exclusive phi electroproduction event generator with tunable LL/LT/TT SDMEs.
 *
 * Usage:
 *   ./phi_gen [config_file]
 *
 * If no config file is provided, built-in defaults are used.
 * Events are written in LUND format to stdout (redirect to a file).
 *
 * Example:
 *   ./phi_gen config/default.cfg > events_default.lund
 *   ./phi_gen config/schc.cfg    > events_schc.lund        # SCHC test
 *   ./phi_gen config/noLT.cfg    > events_noLT.lund        # LT off
 *
 * The LUND event header line contains extra columns:
 *   ... Q2  xB  prodTheta  prodPhi  decayTheta  decayPhi
 * which are useful for analysis cross-checks.
 *
 * Particles per event (in order):
 *   0  beam electron   (e_in)
 *   1  target proton   (p_in)
 *   2  scattered e-    (e_out)
 *   3  recoil proton   (p_out)
 *   4  phi meson       (vec)
 *   5  K-              (decayOne)
 *   6  K+              (decayTwo)
 */

#include <cstdio>
#include <cstring>
#include <string>
#include <vector>

#include "TRandom3.h"
#include "generator.h"
#include "PhysicsParams.h"
#include "ProgressBar.h"

#ifdef _OPENMP
#  include <omp.h>
#endif

int main(int argc, char* argv[])
{
    // ── 1.  Load parameters ─────────────────────────────────────────────────
    PhysicsParams params;

    if (argc > 1) {
        std::string cfgFile = argv[1];
        try {
            params.load(cfgFile);
            fprintf(stderr, "[phi_gen] Loaded config: %s\n", cfgFile.c_str());
        } catch (const std::exception& e) {
            fprintf(stderr, "[phi_gen] ERROR: %s\n", e.what());
            return 1;
        }
    } else {
        fprintf(stderr, "[phi_gen] No config file given — using built-in defaults.\n");
        fprintf(stderr, "[phi_gen] Usage: %s [config_file]\n", argv[0]);
    }

    params.print();

    // ── 2.  Build the reaction candidate ────────────────────────────────────
    // phi → K+ K-   (PDG IDs: 333 → 321 + (-321))
    // K± mass: 0.49368 GeV
    sim::candidate cr(params.beamEnergy, 1.02, 0.49368, 0.49368);
    cr.react.setDecayIds(333, 321, -321);

    // ── 3.  Build the SDME generator ────────────────────────────────────────
    sim::SDMEGenerator gen(&cr, params);
    gen.setRange(params.Q2min, params.Q2max, params.xBmin, params.xBmax);

    // ── 4.  Phase-space scan to find max weight ──────────────────────────────
    fprintf(stderr, "[phi_gen] Scanning %d points for max weight...\n", params.nScan);
    double maxWeight = gen.scan(params.nScan);
    fprintf(stderr, "[phi_gen] Max weight = %.6e  (×1.2 safety applied)\n", maxWeight);

    // ── 5.  Generate events (parallel) ──────────────────────────────────────
    fprintf(stderr, "[phi_gen] Generating %d events...\n", params.nEvents);

#ifdef _OPENMP
    int nThreads = omp_get_max_threads();
    fprintf(stderr, "[phi_gen] Using %d OpenMP threads.\n", nThreads);
#else
    int nThreads = 1;
#endif

    // Each thread fills its own event buffer — no shared state.
    // After the parallel section we write them to stdout in thread order.
    std::vector<std::vector<sim::event>> threadEvents(nThreads);
    int chunkSize = (params.nEvents + nThreads - 1) / nThreads;

    // Progress bar — shows ev/s and ETA across all threads combined
    ProgressBar pb(params.nEvents, "Generating", 40);

#ifdef _OPENMP
#pragma omp parallel num_threads(nThreads)
#endif
    {
#ifdef _OPENMP
        int tid = omp_get_thread_num();
#else
        int tid = 0;
#endif
        sim::candidate crT(params.beamEnergy, 1.02, 0.49368, 0.49368);
        crT.react.setDecayIds(333, 321, -321);
        sim::SDMEGenerator genT(&crT, params);
        genT.setRange(params.Q2min, params.Q2max, params.xBmin, params.xBmax);
        genT.setSeed(42 + tid * 999983u);
        crT.react.rng.SetSeed(137 + tid * 999983u);

        int target = (tid < nThreads - 1) ? chunkSize
                                           : (params.nEvents - chunkSize * (nThreads-1));
        if (target < 0) target = 0;

        std::vector<sim::event>& buf = threadEvents[tid];
        buf.reserve(target + 10);

        sim::event ev;
        for (int attempts = 0; attempts < target * 10 && (int)buf.size() < target; ++attempts) {
            genT.generate(maxWeight);
            crT.react.getEvent(ev);
            if (ev.hasNaN()) continue;
            buf.push_back(ev);
            pb.tick(1);   // thread-safe atomic increment + rate-limited render
        }
    }
    pb.done();

    // ── 6.  Write all events to stdout (serialised) ───────────────────────────
    int written = 0;
    for (int t = 0; t < nThreads; ++t) {
        for (auto& ev : threadEvents[t]) {
            ev.show();
            ++written;
        }
    }

    // ── 6.  Summary ──────────────────────────────────────────────────────────
    fprintf(stderr, "[phi_gen] Written %d events.\n", written);
    gen.stats();

    return 0;
}
