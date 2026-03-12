/**
 * EventCache_impl.h
 *
 * Defines EventCache::fillBatch() — must be included AFTER GavFastMC.h
 * because it calls GavFastMC::processBatch().
 *
 * Include order in main_fit.cpp:
 *   #include "EventCache.h"
 *   #include "GavFastMC.h"
 *   #include "EventCache_impl.h"   ← after both
 */

#pragma once
#include "EventCache.h"
#include "GavFastMC.h"

inline void EventCache::fillBatch(const PhysicsParams&,
                                  const BinDef& bins,
                                  GavFastMC& gavFMC,
                                  sim::SDMEGenerator& gen,
                                  sim::candidate& cr,
                                  double maxW, int nGen)
{
    fprintf(stderr, "[EventCache] Generating %d events in batches of %d via fastMC.jar...\n",
            nGen, batchSize);

    int nGenerated = 0;
    int nAccepted  = 0;

    while (nGenerated < nGen) {
        int thisChunk = std::min(batchSize, nGen - nGenerated);

        // ── Generate a batch of truth events ──────────────────────────
        std::vector<TruthEvent> truths;
        truths.reserve(thisChunk);

        for (int i = 0; i < thisChunk; ++i) {
            gen.generate(maxW);
            sim::event ev; cr.react.getEvent(ev);
            if (ev.hasNaN()) continue;
            truths.push_back(makeTruth(cr));
        }

        // ── Pass batch through fastMC.jar ──────────────────────────────
        auto recos = gavFMC.processBatch(truths);

        // ── Store accepted events ──────────────────────────────────────
        for (const auto& reco : recos) {
            if (!reco.accepted) continue;
            storeReco(reco, bins);
            ++nAccepted;
        }

        nGenerated += thisChunk;

        fprintf(stderr, "[EventCache]   %d / %d generated  |  %d cached  (%.1f%%)\n",
                nGenerated, nGen, nAccepted,
                100.0 * nAccepted / nGenerated);
    }

    fprintf(stderr, "[EventCache] Done. Accepted %d / %d (%.1f%%)\n",
            nAccepted, nGen, 100.0 * nAccepted / nGen);
}
