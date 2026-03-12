/**
 * EventCache.h
 *
 * Fills and manages the EventCache (std::vector<CachedEvent>).
 *
 * -- Two filling strategies --------------------------------------------------
 *
 *  A)  fillFromLund(lundFile, bins, fastmc=false)          [PREFERRED]
 *        Reads events directly from a phi_gen LUND file.
 *        fastmc=false: truth-level events (closure test, no detector cuts).
 *        fastmc=true : fastMC.jar output LUND (real-data path, detector cuts
 *                      applied via LundReader hit-code filter).
 *        No generator needed, no extra scan step.
 *        This is the correct approach for both closure tests and real data fits.
 *
 *  B)  fill(params, bins, fastmc, nGen, nScan)             [legacy]
 *        Generates events internally using sim::SDMEGenerator,
 *        then passes them through a FastMCInterface.
 *        Still available but requires a working FastMC for real-data use.
 *
 * -- Closure test workflow ---------------------------------------------------
 *
 *   phi_gen config/default.cfg > output/events.lund
 *   EventCache cache;
 *   cache.fillFromLund("output/events.lund", bins);   // fastmc=false
 *   // Cache now holds truth-level events distributed per the gen params.
 *   // AnalysisModule evaluates the model at any params on these events.
 *   // DataSet::fillFromMC() produces pseudo-data.  Fit should recover gen params.
 *
 * -- Real-data workflow ------------------------------------------------------
 *
 *   // Run phi_gen + fastMC.jar on ifarm first:
 *   //   phi_gen config > events_mc.lund
 *   //   java -jar fastMC.jar events_mc.lund events_fmc.lund
 *   EventCache cache;
 *   cache.fillFromLund("events_fmc.lund", bins, true);  // fastmc=true
 *   // Cache holds detector-accepted MC events for acceptance correction.
 *   // DataSet loaded from real experimental data files.
 *   // Chi2Fitter fits real data against acceptance-corrected MC model.
 */

#pragma once
#include "EventData.h"
#include "BinDef.h"
#include "FastMCInterface.h"
#include "LundReader.h"
#include "generator.h"
#include "reaction.h"

#include <vector>
#include <string>
#include <cstdio>
#include <sys/stat.h>   // for file size
#include "ProgressBar.h"

// Forward declare GavFastMC so we can dynamic_cast without including it here
class GavFastMC;

class EventCache {
public:
    std::vector<CachedEvent> events;

    // Chunk size for batch FastMC calls (only used with GavFastMC)
    int batchSize = 10000;

    // ── Fill the cache ─────────────────────────────────────────────────────
    void fill(const PhysicsParams& params,
              const BinDef&        bins,
              FastMCInterface&     fastmc,
              int nGen  = 500000,
              int nScan = 200000)
    {
        events.clear();
        events.reserve(nGen / 5);

        // Build generator
        sim::candidate cr(params.beamEnergy, 1.02, 0.49368, 0.49368);
        cr.react.setDecayIds(333, 321, -321);
        sim::SDMEGenerator gen(&cr, params);
        gen.setRange(params.Q2min, params.Q2max, params.xBmin, params.xBmax);

        fprintf(stderr, "[EventCache] Scanning %d events for max weight...\n", nScan);
        double maxW = gen.scan(nScan);
        fprintf(stderr, "[EventCache] Max weight = %.4e\n", maxW);

        // Detect GavFastMC → use batch mode via virtual flag
        if (fastmc.isBatchFastMC()) {
            fprintf(stderr, "[EventCache] Using GavFastMC BATCH mode (chunk=%d)\n", batchSize);
            // Safe: isBatchFastMC() only returns true for GavFastMC
            GavFastMC* gavFMC = reinterpret_cast<GavFastMC*>(&fastmc);
            fillBatch(params, bins, *gavFMC, gen, cr, maxW, nGen);
        } else {
            fprintf(stderr, "[EventCache] Using sequential FastMC mode\n");
            fillSequential(params, bins, fastmc, gen, cr, maxW, nGen);
        }
    }

    int size()  const { return static_cast<int>(events.size()); }
    bool empty() const { return events.empty(); }

    // ── Fill from LUND file (preferred path) ──────────────────────────────
    /**
     * Read events from a phi_gen LUND output file and store them as
     * CachedEvents ready for AnalysisModule reweighting.
     *
     * @param lundFile   Path to the LUND file.
     * @param bins       Bin definitions for (Q², t', xB).
     * @param fastmc     false = truth-level (no detector cuts)  [closure test]
     *                   true  = fastMC.jar output; apply CLAS12 hit-code cuts
     *                           e'(7,8)  p(5-8)  K±(6,8)      [real data path]
     * @return           Number of events cached.
     */
    int fillFromLund(const std::string& lundFile,
                     const BinDef&      bins,
                     bool               fastmc = false)
    {
        events.clear();
        LundReader reader(lundFile, fastmc);
        LundEvent ev;
        long nRead = 0, nCached = 0;

        fprintf(stderr, "[EventCache] fillFromLund: %s  (fastmc=%s)\n",
                lundFile.c_str(), fastmc ? "ON" : "OFF");

        // Get file size for progress bar
        int64_t fileBytes = 0;
        {
            struct stat st;
            if (stat(lundFile.c_str(), &st) == 0) fileBytes = (int64_t)st.st_size;
        }
        ByteProgressBar bpb(fileBytes, "Reading LUND", 40);

        while (reader.next(ev)) {
            ++nRead;
            // Update byte progress ~every 10k events
            if (nRead % 10000 == 0)
                bpb.update(reader.bytesRead());

            if (!ev.valid) continue;           // rejected by detector cuts or bad kine

            // Map kinematics from LundEvent → CachedEvent
            int iQ2 = bins.findQ2(ev.Q2);
            int it   = bins.findT(ev.tprime);
            int ixB  = bins.findXB(ev.xB);
            if (iQ2 < 0 || it < 0 || ixB < 0) continue;

            // sigma3 prefactor: (α_em / 2π) × y² / (1-ε) × (1-xB) / (xB × Q²)
            double sigma3_pre = (alem / (2.0 * M_PI))
                              * (ev.y * ev.y) / (1.0 - ev.eps)
                              * (1.0 - ev.xB) / ev.xB / ev.Q2;

            CachedEvent ce;
            ce.iQ2 = iQ2;       ce.it  = it;      ce.ixB = ixB;
            ce.Q2  = ev.Q2;     ce.xB  = ev.xB;
            ce.t   = ev.t;      ce.tprime = ev.tprime;
            ce.eps = ev.eps;    ce.y   = ev.y;
            ce.prodPhi    = ev.prodPhi;
            ce.decayTheta = ev.decayTheta;
            ce.decayPhi   = ev.decayPhi;
            ce.beamPol    = ev.beamPol;
            ce.sigma3_prefactor = sigma3_pre;
            ce.t_min = compute_tmin(ev.xB, ev.Q2);
            events.push_back(ce);
            ++nCached;
        }
        bpb.done();

        reader.printStats();
        fprintf(stderr, "[EventCache] fillFromLund complete: cached %ld / %ld events\n",
                nCached, nRead);
        return static_cast<int>(nCached);
    }

private:
    static constexpr double Mp   = 0.93827;
    static constexpr double alem = 1.0 / 137.035999084;

    // ── Build TruthEvent from generator state ─────────────────────────────
    static TruthEvent makeTruth(sim::candidate& cr) {
        TruthEvent truth;
        truth.Q2         = cr.react.Q2();
        truth.xB         = cr.react.xB();
        truth.y          = cr.react.Y();
        truth.eps        = cr.react.Eps();
        truth.prodPhi    = cr.react.prodPhi;
        truth.decayTheta = cr.react.decayTheta;
        truth.decayPhi   = cr.react.decayPhi;
        truth.beamPol    = (cr.react.rpol < 0) ? -1.0 : 1.0;
        double t         = (cr.react.p_out - cr.react.p_in).m2();
        truth.t          = t;
        truth.tprime     = std::abs(t) - std::abs(compute_tmin(truth.xB, truth.Q2));
        truth.e_px  = cr.react.e_out.px();   truth.e_py  = cr.react.e_out.py();
        truth.e_pz  = cr.react.e_out.pz();   truth.e_E   = cr.react.e_out.e();
        truth.p_px  = cr.react.p_out.px();   truth.p_py  = cr.react.p_out.py();
        truth.p_pz  = cr.react.p_out.pz();   truth.p_E   = cr.react.p_out.e();
        truth.kp_px = cr.react.decayTwo.px(); truth.kp_py = cr.react.decayTwo.py();
        truth.kp_pz = cr.react.decayTwo.pz(); truth.kp_E  = cr.react.decayTwo.e();
        truth.km_px = cr.react.decayOne.px(); truth.km_py = cr.react.decayOne.py();
        truth.km_pz = cr.react.decayOne.pz(); truth.km_E  = cr.react.decayOne.e();
        return truth;
    }

    // ── Store one accepted RecoEvent into the cache ────────────────────────
    void storeReco(const RecoEvent& reco, const BinDef& bins) {
        int iQ2 = bins.findQ2(reco.Q2);
        int it  = bins.findT(reco.tprime);
        int ixB = bins.findXB(reco.xB);
        if (iQ2 < 0 || it < 0 || ixB < 0) return;

        double sigma3_pre = (alem / (2.0 * M_PI))
                          * (reco.y * reco.y) / (1.0 - reco.eps)
                          * (1.0 - reco.xB) / reco.xB / reco.Q2;

        CachedEvent ce;
        ce.iQ2 = iQ2; ce.it = it; ce.ixB = ixB;
        ce.Q2 = reco.Q2; ce.xB = reco.xB;
        ce.t = reco.t; ce.tprime = reco.tprime;
        ce.eps = reco.eps; ce.y = reco.y;
        ce.prodPhi = reco.prodPhi;
        ce.decayTheta = reco.decayTheta;
        ce.decayPhi = reco.decayPhi;
        ce.beamPol = reco.beamPol;
        ce.sigma3_prefactor = sigma3_pre;
        ce.t_min = compute_tmin(reco.xB, reco.Q2);
        events.push_back(ce);
    }

    // ── Sequential fill (DummyFastMC) ─────────────────────────────────────
    void fillSequential(const PhysicsParams&, const BinDef& bins,
                        FastMCInterface& fastmc,
                        sim::SDMEGenerator& gen, sim::candidate& cr,
                        double maxW, int nGen)
    {
        fprintf(stderr, "[EventCache] Generating %d events (sequential)...\n", nGen);
        int nAccepted = 0;
        for (int i = 0; i < nGen; ++i) {
            if (i % 100000 == 0 && i > 0)
                fprintf(stderr, "[EventCache]   %d / %d  (cached: %d)\n",
                        i, nGen, nAccepted);
            gen.generate(maxW);
            sim::event ev; cr.react.getEvent(ev);
            if (ev.hasNaN()) continue;

            TruthEvent truth = makeTruth(cr);
            RecoEvent reco;
            if (!fastmc.process(truth, reco)) continue;
            storeReco(reco, bins);
            ++nAccepted;
        }
        fprintf(stderr, "[EventCache] Done. Accepted %d / %d (%.1f%%)\n",
                nAccepted, nGen, 100.0*nAccepted/nGen);
    }

    // ── Batch fill (GavFastMC) ─────────────────────────────────────────────
    void fillBatch(const PhysicsParams&, const BinDef& bins,
                   GavFastMC& gavFMC,
                   sim::SDMEGenerator& gen, sim::candidate& cr,
                   double maxW, int nGen);
    // Defined below GavFastMC include in EventCache_impl.h — see note at bottom
};
