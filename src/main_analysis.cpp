/**
 * main_analysis.cpp
 *
 * Analysis executable: dσ/dt' and A_LU(BSA) vs cos(θ_H) from phi_gen LUND output.
 *
 * Usage:
 *   ./phi_analysis <input.lund> [options]
 *
 * Options:
 *   --output-dir <dir>   Output directory (default: output)
 *   --fastmc             Input is fastMC.jar output — apply CLAS12 detector cuts:
 *                          e'(DC+ECAL)  p(DC)  K±(DC+FTOF)
 *   --help               Show this message
 *
 * Examples:
 *   # Truth-level (no FastMC):
 *   ./phi_analysis output/events.lund
 *
 *   # FastMC-filtered (input already processed by fastMC.jar):
 *   ./phi_analysis output/events_fastmc.lund --fastmc --output-dir output_fastmc
 */

#include <cstdio>
#include <string>
#include <stdexcept>

#include "LundReader.h"
#include "PhiAnalysis.h"

static void usage(const char* prog) {
    printf("Usage: %s <input.lund> [options]\n\n", prog);
    printf("Options:\n");
    printf("  --output-dir <dir>   Output directory (default: output)\n");
    printf("  --fastmc             Input is fastMC.jar output; apply detector cuts\n");
    printf("                         e'(DC+ECAL)  p(DC)  K+(DC+FTOF)  K-(DC+FTOF)\n");
    printf("  --help               Show this message\n\n");
    printf("Examples:\n");
    printf("  ./phi_analysis output/events.lund\n");
    printf("  ./phi_analysis output/events_fastmc.lund --fastmc --output-dir output_fmc\n\n");
}

int main(int argc, char* argv[])
{
    printf("\n");
    printf("╔══════════════════════════════════════════════════════╗\n");
    printf("║   phi SDME Analysis — CLAS12 RGA φ electroproduction ║\n");
    printf("╚══════════════════════════════════════════════════════╝\n\n");

    if (argc < 2) { usage(argv[0]); return 1; }

    std::string lundFile;
    std::string outputDir = "output";
    bool useFastMC = false;

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if      (arg == "--help" || arg == "-h")    { usage(argv[0]); return 0; }
        else if (arg == "--fastmc")                  useFastMC = true;
        else if (arg == "--output-dir" && i+1<argc)  outputDir = argv[++i];
        else if (arg[0] != '-')                      lundFile  = arg;
        else {
            fprintf(stderr, "[phi_analysis] Unknown option: %s\n", arg.c_str());
            usage(argv[0]); return 1;
        }
    }

    if (lundFile.empty()) {
        fprintf(stderr, "[phi_analysis] ERROR: no input LUND file given.\n\n");
        usage(argv[0]); return 1;
    }

    printf("[phi_analysis] Input LUND  : %s\n", lundFile.c_str());
    printf("[phi_analysis] Output dir  : %s\n", outputDir.c_str());
    printf("[phi_analysis] FastMC mode : %s\n\n",
           useFastMC ? "ON  (e'→ECAL, p→DC, K±→FTOF)" : "OFF (truth-level, no detector cuts)");

    PhiAnalysis ana;
    ana.setOutputDir(outputDir);

    printf("=== Q² bin edges (%d bins) ===\n", ana.nQ2);
    for (int i = 0; i <= ana.nQ2; ++i)
        printf("  edge[%d] = %.4f\n", i, ana.Q2edges[i]);
    printf("\n=== t' bin edges (%d bins) ===\n", ana.nTp);
    for (int i = 0; i <= ana.nTp; ++i)
        printf("  edge[%d] = %.4f\n", i, ana.Tpedges[i]);
    printf("\n=== cos(θ_H) : %d uniform bins in [-1, +1] ===\n\n", PhiAnalysis::NCosH);

    try {
        ana.processFile(lundFile, useFastMC);
    } catch (const std::exception& e) {
        fprintf(stderr, "[phi_analysis] ERROR: %s\n", e.what());
        return 1;
    }

    if (ana.nFill == 0) {
        fprintf(stderr, "[phi_analysis] ERROR: 0 events filled.\n");
        if (useFastMC)
            fprintf(stderr, "  All events rejected — is the LUND file genuine fastMC.jar output?\n");
        return 1;
    }

    ana.printSummary();
    printf("\n[phi_analysis] Generating plots...\n");
    ana.makePlots();
    printf("\n[phi_analysis] Done.\n\n");
    return 0;
}
