cmake_minimum_required(VERSION 3.14)
project(phi_sdme_gen CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# ── Find ROOT ────────────────────────────────────────────────────────────────
find_package(ROOT REQUIRED COMPONENTS Core RIO MathCore Minuit2 Graf Gpad Hist)
include(${ROOT_USE_FILE})

# ── Find OpenMP ──────────────────────────────────────────────────────────────
find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    message(STATUS "  OpenMP     : ${OpenMP_CXX_VERSION}  (${OpenMP_CXX_FLAGS})")
    set(OMP_LIBS OpenMP::OpenMP_CXX)
else()
    message(WARNING "  OpenMP NOT found — building single-threaded")
    set(OMP_LIBS "")
endif()

# ── Compiler flags ───────────────────────────────────────────────────────────
# -O3:       aggressive optimisation (auto-vectorise the event loop)
# -march=native: use all CPU features available on the build machine
add_compile_options(-O3 -march=native -Wall -Wextra -Wno-unused-parameter)

# ── Include directories ──────────────────────────────────────────────────────
include_directories(include)
include_directories(${ROOT_INCLUDE_DIRS})

# ── Executable 1: LUND event generator ───────────────────────────────────────
add_executable(phi_gen src/main_gen.cpp)
target_link_libraries(phi_gen ${ROOT_LIBRARIES} ${OMP_LIBS})

# ── Executable 2: Full fitting pipeline ──────────────────────────────────────
add_executable(phi_fit src/main_fit.cpp)
target_link_libraries(phi_fit ${ROOT_LIBRARIES} ${OMP_LIBS})

# ── Executable 3: LUND analysis ──────────────────────────────────────────────
add_executable(phi_analysis src/main_analysis.cpp)
target_link_libraries(phi_analysis ${ROOT_LIBRARIES} ${OMP_LIBS})

# ── Install ──────────────────────────────────────────────────────────────────
install(TARGETS phi_gen phi_fit phi_analysis DESTINATION bin)
install(DIRECTORY config/ DESTINATION config)

message(STATUS "")
message(STATUS "  phi_sdme_gen build summary")
message(STATUS "  --------------------------")
message(STATUS "  ROOT version  : ${ROOT_VERSION}")
message(STATUS "  Build type    : ${CMAKE_BUILD_TYPE}")
message(STATUS "  OpenMP        : ${OpenMP_CXX_FOUND}")
message(STATUS "  Executables:")
message(STATUS "    phi_gen      — LUND event generator  (parallel generation)")
message(STATUS "    phi_fit      — SDME fitting pipeline (parallel Minuit inner loop)")
message(STATUS "    phi_analysis — dσ/dt' and BSA plots")
message(STATUS "")
