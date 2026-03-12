/**
 * ProgressBar.h
 *
 * Lightweight ANSI progress display for stderr.
 * Three modes:
 *
 *   1. ProgressBar     — known-total bar  [=====>    ] 47%  12.3k ev/s  ETA 0:04
 *   2. ByteProgressBar — file-read bar    [=====>    ] 47%  234 MB / 500 MB
 *   3. FitSpinner      — Minuit spinner   ◐  call 1234  χ²=12.34  14.2 calls/s  0:32 elapsed
 *
 * All write to stderr (never stdout) so they don't pollute LUND output.
 * Thread-safe: uses std::atomic for the counter, std::mutex for the render.
 *
 * Usage — event generation:
 *   ProgressBar pb(nEvents, "Generating");
 *   pb.tick(n);          // called from any thread; renders at most 10 Hz
 *   pb.done();           // print final line with total time
 *
 * Usage — file reading:
 *   ByteProgressBar bpb(fileSizeBytes, "Reading LUND");
 *   bpb.update(bytesRead);
 *   bpb.done();
 *
 * Usage — Minuit fit:
 *   FitSpinner sp("Fitting");
 *   sp.tick(currentChi2);   // call inside chi2() callback
 *   sp.done(finalChi2);
 */

#pragma once
#include <cstdio>
#include <cstdint>
#include <cmath>
#include <atomic>
#include <mutex>
#include <chrono>
#include <string>

// ── ANSI helpers ──────────────────────────────────────────────────────────────
namespace pb_ansi {
    static const char* RESET  = "\033[0m";
    static const char* GREEN  = "\033[32m";
    static const char* CYAN   = "\033[36m";
    static const char* YELLOW = "\033[33m";
    static const char* BOLD   = "\033[1m";
    inline bool isatty_stderr() {
#ifdef _WIN32
        return false;
#else
        return ::isatty(fileno(stderr));
#endif
    }
}

// ── Shared time utilities ─────────────────────────────────────────────────────
using Clock     = std::chrono::steady_clock;
using TimePoint = std::chrono::time_point<Clock>;

static inline double elapsed_s(const TimePoint& t0) {
    return std::chrono::duration<double>(Clock::now() - t0).count();
}

static inline std::string fmt_duration(double s) {
    char buf[32];
    int h = (int)s / 3600;
    int m = ((int)s % 3600) / 60;
    int sec = (int)s % 60;
    if (h > 0) snprintf(buf, sizeof(buf), "%d:%02d:%02d", h, m, sec);
    else        snprintf(buf, sizeof(buf), "%d:%02d", m, sec);
    return buf;
}

static inline std::string fmt_rate(double r) {
    char buf[32];
    if      (r >= 1e6) snprintf(buf, sizeof(buf), "%.1fM/s", r/1e6);
    else if (r >= 1e3) snprintf(buf, sizeof(buf), "%.1fk/s", r/1e3);
    else               snprintf(buf, sizeof(buf), "%.0f/s",  r);
    return buf;
}

static inline std::string fmt_bytes(int64_t b) {
    char buf[32];
    if      (b >= (1LL<<30)) snprintf(buf, sizeof(buf), "%.1f GB", b/(double)(1LL<<30));
    else if (b >= (1LL<<20)) snprintf(buf, sizeof(buf), "%.1f MB", b/(double)(1LL<<20));
    else if (b >= (1LL<<10)) snprintf(buf, sizeof(buf), "%.1f KB", b/(double)(1LL<<10));
    else                     snprintf(buf, sizeof(buf), "%lld B", (long long)b);
    return buf;
}

// ── 1. ProgressBar — known total ─────────────────────────────────────────────
class ProgressBar {
public:
    explicit ProgressBar(int64_t total, const std::string& label = "",
                         int width = 40, double minIntervalS = 0.10)
        : total_(total), label_(label), width_(width),
          minInterval_(minIntervalS), count_(0), t0_(Clock::now()), tLast_(t0_)
    {
        color_ = pb_ansi::isatty_stderr();
        if (color_) fprintf(stderr, "\n");   // blank line before bar
    }

    // Add n to counter and re-render if enough time has passed.
    // Thread-safe — can be called from any thread.
    void tick(int64_t n = 1) {
        int64_t cur = count_.fetch_add(n, std::memory_order_relaxed) + n;
        // Rate-limit rendering to minInterval_ seconds
        auto now = Clock::now();
        {
            std::lock_guard<std::mutex> lk(mu_);
            double dt = std::chrono::duration<double>(now - tLast_).count();
            if (dt < minInterval_ && cur < total_) return;
            tLast_ = now;
        }
        render(cur, Clock::now());
    }

    // Print the final completed bar with total time.
    void done() {
        int64_t cur = count_.load();
        render(total_ > 0 ? total_ : cur, Clock::now(), /*final=*/true);
        fprintf(stderr, "\n");
    }

private:
    void render(int64_t cur, TimePoint now, bool final = false) const {
        double frac    = total_ > 0 ? std::min(1.0, (double)cur / total_) : 0.0;
        int    filled  = (int)(frac * width_);
        double dt      = elapsed_s(t0_);
        double rate    = dt > 0 ? cur / dt : 0.0;
        double etaSec  = (rate > 0 && !final) ? (total_ - cur) / rate : 0.0;

        const char* c1 = color_ ? pb_ansi::GREEN  : "";
        const char* c2 = color_ ? pb_ansi::CYAN   : "";
        const char* cr = color_ ? pb_ansi::RESET  : "";
        const char* cb = color_ ? pb_ansi::BOLD   : "";

        // Build bar string
        char bar[256] = {};
        bar[0] = '[';
        for (int i = 1; i <= width_; ++i) {
            if (i < filled)       bar[i] = '=';
            else if (i == filled) bar[i] = '>';
            else                  bar[i] = ' ';
        }
        bar[width_+1] = ']';
        bar[width_+2] = '\0';

        std::string timeStr = final ? fmt_duration(dt)
                                    : ("ETA " + fmt_duration(etaSec));

        fprintf(stderr, "\r%s%-14s%s %s%s%s %s%5.1f%%%s  %s  %s%-10s%s",
                cb, label_.c_str(), cr,
                c1, bar, cr,
                cb, frac*100.0, cr,
                fmt_rate(rate).c_str(),
                c2, timeStr.c_str(), cr);
        if (final) fprintf(stderr, "  %s✓%s", c1, cr);
        fflush(stderr);
    }

    int64_t         total_;
    std::string     label_;
    int             width_;
    double          minInterval_;
    std::atomic<int64_t> count_;
    TimePoint       t0_;
    TimePoint       tLast_;
    mutable std::mutex mu_;
    bool            color_;
};

// ── 2. ByteProgressBar — file read progress ───────────────────────────────────
class ByteProgressBar {
public:
    explicit ByteProgressBar(int64_t totalBytes, const std::string& label = "",
                             int width = 40, double minIntervalS = 0.15)
        : total_(totalBytes), label_(label), width_(width),
          minInterval_(minIntervalS), bytes_(0), t0_(Clock::now()), tLast_(t0_)
    {
        color_ = pb_ansi::isatty_stderr();
        if (color_) fprintf(stderr, "\n");
    }

    void update(int64_t bytesRead) {
        bytes_.store(bytesRead, std::memory_order_relaxed);
        auto now = Clock::now();
        {
            std::lock_guard<std::mutex> lk(mu_);
            double dt = std::chrono::duration<double>(now - tLast_).count();
            if (dt < minInterval_ && bytesRead < total_) return;
            tLast_ = now;
        }
        render(bytesRead, now);
    }

    void done() {
        render(total_ > 0 ? total_ : bytes_.load(), Clock::now(), true);
        fprintf(stderr, "\n");
    }

private:
    void render(int64_t cur, TimePoint now, bool final = false) const {
        double frac   = total_ > 0 ? std::min(1.0, (double)cur / total_) : 0.0;
        int    filled = (int)(frac * width_);
        double dt     = elapsed_s(t0_);
        double rate   = dt > 0 ? cur / dt : 0.0;   // bytes/s
        double etaSec = (rate > 0 && !final) ? (total_ - cur) / rate : 0.0;

        const char* c1 = color_ ? pb_ansi::GREEN : "";
        const char* c2 = color_ ? pb_ansi::CYAN  : "";
        const char* cr = color_ ? pb_ansi::RESET : "";
        const char* cb = color_ ? pb_ansi::BOLD  : "";

        char bar[256] = {};
        bar[0] = '[';
        for (int i = 1; i <= width_; ++i) {
            if (i < filled)       bar[i] = '=';
            else if (i == filled) bar[i] = '>';
            else                  bar[i] = ' ';
        }
        bar[width_+1] = ']';
        bar[width_+2] = '\0';

        // Format: "234 MB / 500 MB"
        std::string progress = fmt_bytes(cur) + " / " + fmt_bytes(total_);
        std::string timeStr  = final ? fmt_duration(dt) : ("ETA " + fmt_duration(etaSec));

        // Byte rate in MB/s
        char rateStr[32];
        double mbps = rate / (1<<20);
        if (mbps >= 1.0) snprintf(rateStr, sizeof(rateStr), "%.1f MB/s", mbps);
        else             snprintf(rateStr, sizeof(rateStr), "%.0f KB/s", rate/1024.0);

        fprintf(stderr, "\r%s%-14s%s %s%s%s %s%5.1f%%%s  %-18s  %s  %s%s%s",
                cb, label_.c_str(), cr,
                c1, bar, cr,
                cb, frac*100.0, cr,
                progress.c_str(),
                rateStr,
                c2, timeStr.c_str(), cr);
        if (final) fprintf(stderr, "  %s✓%s", c1, cr);
        fflush(stderr);
    }

    int64_t              total_;
    std::string          label_;
    int                  width_;
    double               minInterval_;
    std::atomic<int64_t> bytes_;
    TimePoint            t0_;
    TimePoint            tLast_;
    mutable std::mutex   mu_;
    bool                 color_;
};

// ── 3. FitSpinner — Minuit progress (unknown total) ───────────────────────────
// Shows: spinner  call N  χ²=12.34  14.2 calls/s  0:32 elapsed
class FitSpinner {
    static constexpr const char* FRAMES[] = {"◐","◓","◑","◒"};
    static constexpr int NFRAMES = 4;
public:
    explicit FitSpinner(const std::string& label = "Fitting",
                        double minIntervalS = 0.15)
        : label_(label), minInterval_(minIntervalS),
          calls_(0), chi2_(1e30), t0_(Clock::now()), tLast_(t0_), frame_(0)
    {
        color_ = pb_ansi::isatty_stderr();
        if (color_) fprintf(stderr, "\n");
    }

    // Call inside chi2() callback.  Thread-safe.
    void tick(double currentChi2) {
        int64_t n = calls_.fetch_add(1, std::memory_order_relaxed) + 1;
        // Track best chi2
        double prev = chi2_.load(std::memory_order_relaxed);
        while (currentChi2 < prev &&
               !chi2_.compare_exchange_weak(prev, currentChi2,
                                            std::memory_order_relaxed)) {}

        auto now = Clock::now();
        {
            std::lock_guard<std::mutex> lk(mu_);
            double dt = std::chrono::duration<double>(now - tLast_).count();
            if (dt < minInterval_) return;
            tLast_ = now;
            frame_ = (frame_ + 1) % NFRAMES;
        }
        render(n, currentChi2, now);
    }

    void done(double finalChi2 = -1.0) {
        int64_t n = calls_.load();
        double  c = finalChi2 >= 0 ? finalChi2 : chi2_.load();
        render(n, c, Clock::now(), true);
        fprintf(stderr, "\n");
    }

private:
    void render(int64_t n, double c2, TimePoint now, bool final = false) const {
        double dt   = elapsed_s(t0_);
        double rate = dt > 0 ? n / dt : 0.0;

        const char* c1 = color_ ? pb_ansi::YELLOW : "";
        const char* c2c= color_ ? pb_ansi::CYAN   : "";
        const char* cr = color_ ? pb_ansi::RESET  : "";
        const char* cb = color_ ? pb_ansi::BOLD   : "";

        const char* sym = final ? "✓" : FRAMES[frame_];

        fprintf(stderr,
                "\r%s%s%s  %s%-14s%s  call %-6lld  %sχ²=%-10.4f%s  %s  %s%s elapsed%s",
                c1, sym, cr,
                cb, label_.c_str(), cr,
                (long long)n,
                cb, c2, cr,
                fmt_rate(rate).c_str(),
                c2c, fmt_duration(dt).c_str(), cr);
        if (final) fprintf(stderr, "  %s✓%s", c1, cr);
        fflush(stderr);
    }

    std::string          label_;
    double               minInterval_;
    std::atomic<int64_t> calls_;
    std::atomic<double>  chi2_;
    TimePoint            t0_;
    TimePoint            tLast_;
    mutable std::mutex   mu_;
    int                  frame_;
    bool                 color_;
};
