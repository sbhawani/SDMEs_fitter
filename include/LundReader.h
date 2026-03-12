/**
 * LundReader.h
 *
 * Reads LUND output from phi_gen (truth) or fastMC.jar (detector-filtered).
 *
 * Two modes:
 *   TRUTH (fastmcMode=false): all kinematically valid events accepted.
 *   FASTMC (fastmcMode=true): col 6 of each particle line = CLAS12 hit code.
 *     e'  needs ECAL  -> codes {7,8}
 *     p   needs DC    -> codes {5,6,7,8}
 *     K+- needs FTOF  -> codes {6,8}
 *
 * LUND header cols:
 *   0:nPart 1:1 2:0 3:beamPol 4:0 5:11 6:beamE 7:2212 8:1 9:0.6
 *   10:Q2 11:xB 12:prodTheta 13:prodPhi 14:decayTheta 15:decayPhi
 *
 * Particle order (7 per event):
 *   0:beam e-  1:target p  2:scattered e-  3:recoil p
 *   4:phi      5:K-        6:K+
 *
 * Particle line cols:
 *   0:idx 1:charge 2:type 3:pid 4:status 5:detCode
 *   6:px 7:py 8:pz 9:E 10:mass 11:vx 12:vy 13:vz
 */

#pragma once
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>
#include <stdexcept>
#include <algorithm>
#include <cstdio>

struct LundParticle {
    int    pid     = 0;
    int    status  = 0;
    int    detCode = 0;  // 0 in truth; CLAS12 hit code in fastMC mode
    double px=0, py=0, pz=0, E=0, mass=0;
    double p() const { return std::sqrt(px*px+py*py+pz*pz); }
};

struct LundEvent {
    int    nPart      = 0;
    double beamPol    = 0;
    double Q2         = 0;
    double xB         = 0;
    double prodTheta  = 0;
    double prodPhi    = 0;
    double decayTheta = 0;  // K+ polar angle in phi helicity frame
    double decayPhi   = 0;

    double tprime  = 0;
    double t       = 0;
    double W       = 0;
    double y       = 0;
    double eps     = 0;

    LundParticle eBeam, pTarget, eOut, pOut, phi, Km, Kp;

    bool valid          = false;
    bool fastmcRejected = false;
};

class LundReader {
public:
    bool fastmcMode = false;

    // Acceptable CLAS12 hit codes per particle type
    std::vector<int> eAccept  = {7, 8};        // e': DC + ECAL
    std::vector<int> pAccept  = {5, 6, 7, 8};  // p: DC sufficient
    std::vector<int> kAccept  = {6, 8};         // K+-: DC + FTOF

    explicit LundReader(const std::string& filename, bool fastmc = false)
        : fastmcMode(fastmc), file_(filename)
    {
        if (!file_.is_open())
            throw std::runtime_error("[LundReader] Cannot open: " + filename);
        if (fastmcMode)
            fprintf(stderr, "[LundReader] FastMC mode ON -- e'(7,8)  p(5-8)  K+-(6,8)\n");
    }

    bool good()      const { return file_.good() && !file_.eof(); }
    long nRead()     const { return nRead_; }
    long nAccepted() const { return nAccepted_; }
    long nRejected() const { return nRejected_; }
    int64_t bytesRead() const {
        // tellg() returns the current position; cast to int64_t for ByteProgressBar
        return static_cast<int64_t>(const_cast<std::ifstream&>(file_).tellg());
    }

    void printStats() const {
        fprintf(stderr,
            "[LundReader] Read: %ld  Accepted: %ld  Rejected(fastMC): %ld\n",
            nRead_, nAccepted_, nRejected_);
    }

    bool next(LundEvent& ev) {
        ev = LundEvent{};
        std::string line;

        while (std::getline(file_, line)) {
            if (line.empty() || line[0] == '#') continue;

            std::istringstream hdr(line);
            hdr >> ev.nPart;
            if (hdr.fail() || ev.nPart < 1) continue;

            int d1, d2, bpid, tpid, d3;
            double tpol, bE, d4;
            hdr >> d1 >> d2 >> ev.beamPol >> tpol
                >> bpid >> bE >> tpid >> d3 >> d4
                >> ev.Q2 >> ev.xB
                >> ev.prodTheta >> ev.prodPhi
                >> ev.decayTheta >> ev.decayPhi;

            bool ok = true;
            for (int i = 0; i < ev.nPart; ++i) {
                if (!std::getline(file_, line)) { ok = false; break; }
                LundParticle part;
                std::istringstream ps(line);
                int idx, charge, type, pid, status2, detCode;
                double vx, vy, vz;
                ps >> idx >> charge >> type >> pid >> status2 >> detCode
                   >> part.px >> part.py >> part.pz >> part.E >> part.mass
                   >> vx >> vy >> vz;
                if (ps.fail()) { ok = false; break; }
                part.pid     = pid;
                part.status  = type;
                part.detCode = detCode;
                switch (i) {
                    case 0: ev.eBeam   = part; break;
                    case 1: ev.pTarget = part; break;
                    case 2: ev.eOut    = part; break;
                    case 3: ev.pOut    = part; break;
                    case 4: ev.phi     = part; break;
                    case 5: ev.Km      = part; break;
                    case 6: ev.Kp      = part; break;
                }
            }
            if (!ok) return false;

            ++nRead_;

            // FastMC detector acceptance cut
            if (fastmcMode) {
                bool eOK  = codeOK(ev.eOut.detCode, eAccept);
                bool pOK  = codeOK(ev.pOut.detCode, pAccept);
                bool kpOK = codeOK(ev.Kp.detCode,   kAccept);
                bool kmOK = codeOK(ev.Km.detCode,    kAccept);
                if (!eOK || !pOK || !kpOK || !kmOK) {
                    ev.fastmcRejected = true;
                    ++nRejected_;
                    return true;
                }
            }

            // Compute derived kinematics (same formulas as EventData.h)
            static constexpr double Mp    = 0.93827;
            static constexpr double Mv    = 1.020;
            static constexpr double beamE = 10.6;

            double dpx = ev.pOut.px, dpy = ev.pOut.py, dpz = ev.pOut.pz;
            double dE  = ev.pOut.E - Mp;
            ev.t = dE*dE - (dpx*dpx + dpy*dpy + dpz*dpz);

            double W2    = Mp*Mp + ev.Q2*(1.0/ev.xB - 1.0);
            ev.W         = (W2 > 0) ? std::sqrt(W2) : 0;
            double nu    = ev.Q2 / (2.0 * Mp * ev.xB);
            double qstar = std::sqrt(nu*nu + ev.Q2);
            double Ephi  = (W2 + Mv*Mv - Mp*Mp) / (2.0*ev.W);
            double pphi  = std::sqrt(std::max(Ephi*Ephi - Mv*Mv, 0.0));
            double Epr   = (W2 - Mv*Mv + Mp*Mp) / (2.0*ev.W);
            double pq    = std::sqrt(qstar*qstar + ev.Q2);
            double tmin  = ev.Q2 - 2.0*(nu*Epr - pq*pphi);
            ev.tprime    = std::abs(ev.t) - std::abs(tmin);

            ev.y   = nu / beamE;
            double gg   = 2.0*ev.xB*Mp / std::sqrt(ev.Q2);
            double y2g2 = ev.y*ev.y*gg*gg;
            ev.eps = (1.0 - ev.y - 0.25*y2g2) /
                     (1.0 - ev.y + 0.5*ev.y*ev.y + 0.25*y2g2);

            ev.valid = true;
            ++nAccepted_;
            return true;
        }
        return false;
    }

private:
    std::ifstream file_;
    long nRead_     = 0;
    long nAccepted_ = 0;
    long nRejected_ = 0;

    bool codeOK(int code, const std::vector<int>& allowed) const {
        return std::find(allowed.begin(), allowed.end(), code) != allowed.end();
    }
};
