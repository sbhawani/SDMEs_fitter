/**
 * GavFastMC.h
 *
 * Wrapper around Gavalian's CLAS12 FastMC:
 *   java -jar /home/gavalian/Software/fastMC/fastMC.jar input.lund output.lund
 *
 * The FastMC reads a LUND file and modifies column 6 of each particle line
 * to indicate which detectors were hit:
 *
 *   0  — no detector hit                     → reject
 *   5  — DC only (all 6 superlayers)
 *   6  — DC + FTOF  (no ECAL)
 *   7  — DC + ECAL  (no FTOF)
 *   8  — DC + FTOF + ECAL
 *
 * Acceptance criteria for ep → e'p K⁺K⁻:
 *   electron :  needs ECAL for trigger/ID    → codes 7, 8
 *   proton   :  DC only acceptable (no FTOF/ECAL required in FD)  → codes 5,6,7,8
 *   K⁺       :  needs FTOF for TOF-based PID → codes 6, 8
 *   K⁻       :  same                          → codes 6, 8
 *
 * Integration strategy — BATCH mode:
 *   Because java has ~1s startup overhead, we don't call it per-event.
 *   Instead EventCache::fill() calls GavFastMC::processBatch():
 *     1.  Write N generated events to a temp LUND file
 *     2.  Call fastMC.jar on the whole file
 *     3.  Parse the output LUND — extract accepted events
 *   This is identical to the DummyFastMC interface but operates on batches.
 *
 * Usage:
 *   GavFastMC fmc("/home/gavalian/Software/fastMC/fastMC.jar");
 *   fmc.setJava("/usr/bin/java");          // optional, default = "java"
 *   fmc.setTmpDir("/scratch/phi_fit");     // optional, default = /tmp
 *
 *   // Check it works:
 *   if (!fmc.selfTest()) { std::cerr << "fastMC not available\n"; }
 *
 *   // Batch processing (called by EventCache):
 *   std::vector<TruthEvent> truths = ...; // from generator
 *   auto recos = fmc.processBatch(truths);
 *   // recos[i].accepted == true/false
 *
 * The single-event FastMCInterface::process() is also implemented for
 * compatibility with DummyFastMC — it just wraps a one-event batch call
 * (inefficient; use processBatch for production runs).
 */

#pragma once
#include "FastMCInterface.h"
#include "EventData.h"

#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <cstdio>
#include <cstdlib>
#include <stdexcept>
#include <filesystem>
#include <array>
#include <cmath>
#include <unistd.h>

namespace fs = std::filesystem;

class GavFastMC : public FastMCInterface {
public:

    // ── Configuration ────────────────────────────────────────────────────────
    std::string jarPath;                  // path to fastMC.jar
    std::string javaExe  = "java";        // java executable
    std::string tmpDir   = "/tmp";        // scratch dir for temp LUND files

    // Detector hit requirements per particle type
    // LUND column 6 value must be in the acceptance set
    std::vector<int> eAccept  = {7, 8};   // electron: needs ECAL
    std::vector<int> pAccept  = {5,6,7,8};// proton: DC sufficient
    std::vector<int> kAccept  = {6, 8};   // K±: needs FTOF for PID

    // LUND particle IDs
    static constexpr int PID_ELECTRON =  11;
    static constexpr int PID_PROTON   =  2212;
    static constexpr int PID_KP       =  321;
    static constexpr int PID_KM       = -321;

    // ── Constructor ──────────────────────────────────────────────────────────
    explicit GavFastMC(const std::string& jar,
                       const std::string& java = "java",
                       const std::string& tmp  = "/tmp")
        : jarPath(jar), javaExe(java), tmpDir(tmp) {}

    // ── Self-test: verify java + jar are accessible ───────────────────────
    bool isBatchFastMC() const override { return true; }

    bool selfTest() const {
        if (!fs::exists(jarPath)) {
            fprintf(stderr, "[GavFastMC] jar not found: %s\n", jarPath.c_str());
            return false;
        }
        std::string cmd = javaExe + " -version 2>&1";
        int rc = std::system(cmd.c_str());
        if (rc != 0) {
            fprintf(stderr, "[GavFastMC] java not found or not executable\n");
            return false;
        }
        fprintf(stderr, "[GavFastMC] OK: %s\n", jarPath.c_str());
        return true;
    }

    // ── FastMCInterface: single-event (wraps batch — use only for testing) ──
    bool process(const TruthEvent& truth, RecoEvent& reco) override {
        auto results = processBatch({truth});
        if (results.empty()) { reco.accepted = false; return false; }
        reco = results[0];
        return reco.accepted;
    }

    // ── Batch processing (main entry point) ──────────────────────────────────
    /**
     * Process a batch of truth events through fastMC.jar.
     * Returns one RecoEvent per input event (accepted or not).
     */
    std::vector<RecoEvent> processBatch(const std::vector<TruthEvent>& truths)
    {
        if (truths.empty()) return {};

        // ── 1.  Write LUND file ───────────────────────────────────────────
        std::string inFile  = tmpDir + "/phi_fastmc_in_"  + std::to_string(getpid()) + ".lund";
        std::string outFile = tmpDir + "/phi_fastmc_out_" + std::to_string(getpid()) + ".lund";

        writeLund(truths, inFile);

        // ── 2.  Run fastMC.jar ────────────────────────────────────────────
        std::string cmd = javaExe + " -jar " + jarPath
                        + " " + inFile + " " + outFile
                        + " 2>/dev/null";

        int rc = std::system(cmd.c_str());
        if (rc != 0) {
            fprintf(stderr, "[GavFastMC] fastMC.jar returned non-zero: %d\n", rc);
            std::vector<RecoEvent> out(truths.size());
            fs::remove(inFile);
            return out;
        }

        if (!fs::exists(outFile)) {
            fprintf(stderr, "[GavFastMC] output LUND not created: %s\n", outFile.c_str());
            std::vector<RecoEvent> out(truths.size());
            fs::remove(inFile);
            return out;
        }

        // ── 3.  Keep first output LUND in tmpDir for manual inspection ────
        static bool savedFirst = false;
        if (!savedFirst) {
            savedFirst = true;
            std::string keepPath = tmpDir + "/phi_fastmc_FIRST_output.lund";
            fs::copy_file(outFile, keepPath, fs::copy_options::overwrite_existing);
            fprintf(stderr, "[GavFastMC] First output LUND saved to: %s\n",
                    keepPath.c_str());
            fprintf(stderr, "  (first 15 lines):\n");
            std::ifstream dbg(keepPath);
            std::string dbgLine;
            for (int i = 0; i < 15 && std::getline(dbg, dbgLine); ++i)
                fprintf(stderr, "  %s\n", dbgLine.c_str());
            fprintf(stderr, "\n");
        }

        // ── 3.  Parse output LUND and match back to truth events ─────────
        auto recos = parseLund(outFile, truths);

        // ── 4.  Clean up temp files ───────────────────────────────────────
        fs::remove(inFile);
        fs::remove(outFile);

        return recos;
    }

    // ── Write LUND ────────────────────────────────────────────────────────────
    /**
     * Write LUND in the same format as phiproduce.C / reaction.h::getEvent().
     * This includes ALL particles: beam e⁻, target p, scattered e⁻, recoil p,
     * φ meson (intermediate), K⁻, K⁺.  FastMC expects this exact ordering.
     *
     * LUND column convention (CLAS12):
     *   1  index
     *   2  charge
     *   3  type  (1=stable/trackable, 2=intermediate/unstable)
     *   4  pid
     *   5  parent
     *   6  daughter  ← fastMC overwrites this with detector hit code
     *   7-9  px py pz [GeV]
     *   10   E  [GeV]
     *   11   mass [GeV]
     *   12-14 vx vy vz [cm]
     */
    void writeLund(const std::vector<TruthEvent>& truths,
                   const std::string& filename) const
    {
        std::ofstream f(filename);
        if (!f.is_open())
            throw std::runtime_error("[GavFastMC] Cannot write: " + filename);

        static constexpr double Me   = 0.000511;
        static constexpr double Mp   = 0.93827;
        static constexpr double Mk   = 0.49368;
        static constexpr double Mphi = 1.01946;
        static constexpr double beamE = 10.6;

        for (const auto& ev : truths) {
            // Header: 7 particles, A=1, Z=1, polT=0, polB=helicity, beamE, inter=1
            f << "7 1 1 0 " << (int)ev.beamPol << " " << beamE << " 1 7\n";

            // lambda: write one particle line
            //   idx  charge  type  pid  parent  daughter  px py pz E mass vx vy vz
            auto line = [&](int idx, int charge, int type, int pid,
                            double px, double py, double pz, double E, double mass) {
                f << std::fixed; f.precision(5);
                f << " " << idx << " " << charge << " " << type
                  << " " << pid << " 0 0"
                  << "  " << px << " " << py << " " << pz
                  << "  " << E  << "  " << mass
                  << "  0.0 0.0 0.0\n";
            };

            // Reconstruct beam e⁻ and target p from kinematics
            // (fastMC needs them to establish the reaction frame)
            line(1, -1, 0, PID_ELECTRON, 0.0, 0.0, beamE, beamE,   Me);   // beam e⁻
            line(2, +1, 0, PID_PROTON,   0.0, 0.0, 0.0,   Mp,      Mp);   // target p
            line(3, -1, 1, PID_ELECTRON, ev.e_px,  ev.e_py,  ev.e_pz,  ev.e_E,  Me);   // e'
            line(4, +1, 1, PID_PROTON,   ev.p_px,  ev.p_py,  ev.p_pz,  ev.p_E,  Mp);   // p'
            // φ meson (intermediate, not tracked)
            double phi_px = ev.kp_px+ev.km_px, phi_py = ev.kp_py+ev.km_py;
            double phi_pz = ev.kp_pz+ev.km_pz, phi_E  = ev.kp_E +ev.km_E;
            line(5,  0, 2, 333,           phi_px, phi_py, phi_pz, phi_E, Mphi);
            line(6, +1, 1, PID_KP,        ev.kp_px, ev.kp_py, ev.kp_pz, ev.kp_E, Mk);   // K⁺
            line(7, -1, 1, PID_KM,        ev.km_px, ev.km_py, ev.km_pz, ev.km_E, Mk);   // K⁻
        }
    }

    // ── Parse output LUND ─────────────────────────────────────────────────────
    /**
     * Read fastMC output and rebuild RecoEvent list.
     *
     * The output has the same 7-particle structure as input.
     * Column 6 (daughter) of FINAL-STATE particles (type==1) is now the
     * detector hit code (0,5,6,7,8).
     * We skip initial-state (type==0) and intermediate (type==2) particles.
     *
     * Robustness features:
     *   - Selects particles by PID AND type==1 to avoid duplicate PID ambiguity
     *   - Prints first-event diagnostics to stderr for debugging
     *   - Handles output file shorter than input (fastMC may skip events)
     */
    std::vector<RecoEvent> parseLund(const std::string& filename,
                                     const std::vector<TruthEvent>& truths) const
    {
        std::vector<RecoEvent> recos(truths.size());
        std::ifstream f(filename);
        if (!f.is_open()) {
            fprintf(stderr, "[GavFastMC] Cannot open output LUND: %s\n", filename.c_str());
            return recos;
        }

        std::string line;
        int evIdx = 0;
        bool printedDiag = false;

        while (std::getline(f, line) && evIdx < (int)truths.size()) {
            if (line.empty()) continue;

            // Parse header line: first token = nParticles
            std::istringstream hdr(line);
            int nPart;
            hdr >> nPart;
            if (hdr.fail() || nPart < 1 || nPart > 50) continue;

            struct ParticleInfo { int pid; int type; int detCode; double px,py,pz,E; };
            std::vector<ParticleInfo> parts;
            parts.reserve(nPart);

            for (int ip = 0; ip < nPart; ++ip) {
                if (!std::getline(f, line)) break;
                std::istringstream ps(line);
                int idx, charge, type, pid, parent, daughter;
                double px, py, pz, E, mass, vx, vy, vz;
                ps >> idx >> charge >> type >> pid >> parent >> daughter
                   >> px >> py >> pz >> E >> mass >> vx >> vy >> vz;
                if (ps.fail()) continue;
                parts.push_back({pid, type, daughter, px, py, pz, E});
            }

            // ── Diagnostics for first event ───────────────────────────────
            if (!printedDiag) {
                printedDiag = true;
                fprintf(stderr, "[GavFastMC] First output event (%d particles):\n", nPart);
                fprintf(stderr, "  %-4s %-6s %-4s %-6s %-8s\n",
                        "idx","PID","type","detCode","mom");
                for (int i = 0; i < (int)parts.size(); ++i) {
                    double p = std::sqrt(parts[i].px*parts[i].px +
                                        parts[i].py*parts[i].py +
                                        parts[i].pz*parts[i].pz);
                    fprintf(stderr, "  %-4d %-6d %-4d %-6d %.3f\n",
                            i+1, parts[i].pid, parts[i].type, parts[i].detCode, p);
                }
                // Count how many have non-zero detector codes
                int nHit = 0;
                for (auto& p : parts) if (p.detCode != 0) nHit++;
                fprintf(stderr, "  → %d / %d particles hit a detector\n\n", nHit, nPart);
            }

            // ── Map by PID, taking only FINAL-STATE (type==1) particles ──
            // This avoids the duplicate-PID ambiguity (beam e⁻ vs scattered e⁻)
            ParticleInfo *eInfo=nullptr, *pInfo=nullptr,
                         *kpInfo=nullptr, *kmInfo=nullptr;
            for (auto& p : parts) {
                if (p.type != 1) continue;   // skip beam, target, intermediate
                if (p.pid == PID_ELECTRON && !eInfo)  eInfo  = &p;
                if (p.pid == PID_PROTON   && !pInfo)  pInfo  = &p;
                if (p.pid == PID_KP       && !kpInfo) kpInfo = &p;
                if (p.pid == PID_KM       && !kmInfo) kmInfo = &p;
            }

            if (!eInfo || !pInfo || !kpInfo || !kmInfo) { evIdx++; continue; }

            // ── Check acceptance ──────────────────────────────────────────
            bool eOK  = detCodeOK(eInfo->detCode,  eAccept);
            bool pOK  = detCodeOK(pInfo->detCode,  pAccept);
            bool kpOK = detCodeOK(kpInfo->detCode, kAccept);
            bool kmOK = detCodeOK(kmInfo->detCode, kAccept);

            RecoEvent& reco = recos[evIdx];
            reco.beamPol = truths[evIdx].beamPol;

            if (!eOK || !pOK || !kpOK || !kmOK) {
                reco.accepted = false;
                evIdx++;
                continue;
            }

            // ── Recompute DIS variables ────────────────────────────────────
            static constexpr double MpM   = 0.93827;
            static constexpr double beamEv = 10.6;

            double qx = -eInfo->px, qy = -eInfo->py;
            double qz = beamEv - eInfo->pz;
            double qE = beamEv  - eInfo->E;
            double Q2 = -(qE*qE - qx*qx - qy*qy - qz*qz);
            if (Q2 < 1.0) { reco.accepted = false; evIdx++; continue; }

            double nu = qE;
            double xB = Q2 / (2.0 * MpM * nu);
            if (xB < 0.08 || xB > 0.68) { reco.accepted = false; evIdx++; continue; }

            double y   = nu / beamEv;
            double gg2 = 4.0 * xB * xB * MpM * MpM / Q2;
            double eps = (1.0 - y - 0.25*y*y*gg2) /
                         (1.0 - y + 0.5*y*y + 0.25*y*y*gg2);

            double dpE  = pInfo->E  - MpM;
            double dpt2 = pInfo->px*pInfo->px + pInfo->py*pInfo->py + pInfo->pz*pInfo->pz;
            double t    = dpE*dpE - dpt2;
            double tmin = compute_tmin(xB, Q2);
            double tprime = std::abs(t) - std::abs(tmin);
            if (tprime < 0.0 || tprime > 4.5) { reco.accepted = false; evIdx++; continue; }

            reco.Q2         = Q2;
            reco.xB         = xB;
            reco.t          = t;
            reco.tprime     = tprime;
            reco.y          = y;
            reco.eps        = eps;
            reco.prodPhi    = truths[evIdx].prodPhi;
            reco.decayTheta = truths[evIdx].decayTheta;
            reco.decayPhi   = truths[evIdx].decayPhi;
            reco.accepted   = true;

            evIdx++;
        }

        return recos;
    }

private:
    bool detCodeOK(int code, const std::vector<int>& allowed) const {
        for (int a : allowed) if (code == a) return true;
        return false;
    }
};
