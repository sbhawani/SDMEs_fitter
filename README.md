# phi_sdme_gen

Exclusive φ electroproduction event generator with tunable LL / LT / TT helicity amplitude parameters.

Built on top of the Argonne phi generator (Goloskokov-Kroll model base), with the helicity matrix decoupled into explicit physics parameters accessible from a config file.

---

## Physics

The 7-fold differential cross section decomposes as:

```
dσ/dQ² dxB dt d(cosθ) dφ_KK dΦ  =  σ₃(xB,Q²,ε) × (1/4π²) × W(cosθ, φ_KK, Φ; Pl)
```

where

```
W  =  cos²θ · W^LL  +  √2 cosθ sinθ · W^LT  +  sin²θ · W^TT
   +  Pl × [cos²θ · W^LL_LU  +  √2 cosθ sinθ · W^LT_LU  +  sin²θ · W^TT_LU]
```

The three angular kernels W^LL, W^LT, W^TT are computed by `w_kernels.hpp`
from the helicity amplitude matrix `u[ν'][ν][μ'][μ]`.

### Parameter → observable map

| Sector | Key param | Observable |
|--------|-----------|------------|
| LL | `A_L`, `b_L` | dσ_L/dt shape, cos²θ distribution |
| LL | `Im_LL` | sin(φ) in BSA → σ_LT' structure function |
| TT | `A_T`, `b_T` | dσ_T/dt shape, sin²θ distribution |
| TT | `Im_TT` | sin(φ) in BSA from TT interference |
| LT | `A_LT` | cos(φ_KK) and cos(φ+φ_KK) modulations in dσ |
| LT | `A_LTi` | sin(φ_KK) in BSA |

### SCHC test

Set `A_LT = A_LTi = 0` in the config file to run with strict s-channel helicity conservation.
Under SCHC, the only non-zero SDMEs are r⁰⁰₀₄ (gives σ_L/σ_T) and Im(r¹⁰₀₅) (BSA).

---

## Build

```bash
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make -j4
```

Requires: ROOT 6.x (for TRandom3)

---

## Run

```bash
# Default parameters
./phi_gen config/default.cfg > events_default.lund

# SCHC test (LT interference off)
./phi_gen config/schc.cfg > events_schc.lund

# Custom config
cp config/default.cfg config/myfit.cfg
# edit myfit.cfg ...
./phi_gen config/myfit.cfg > events_myfit.lund
```

---

## Output format

LUND format, one event per block:

```
<nParticles> 1 0 <beamPol> 0 11 <beamE> 2212 1 0.6  <Q2> <xB> <prodθ> <prodφ> <decθ> <decφ>
  1  0  11  0  0   px  py  pz  E  m  vx vy vz    # beam e-
  2  0 2212  0  0   ...                            # target p
  3  1  11  1  0   ...                            # scattered e-
  4  1 2212  1  0   ...                            # recoil p
  5  2  333  2  0   ...                            # phi
  6  1 -321  1  0   ...                            # K-
  7  1  321  1  0   ...                            # K+
```

---

## File structure

```
phi_sdme_gen/
├── CMakeLists.txt
├── README.md
├── config/
│   ├── default.cfg       # default LL+LT+TT parameters
│   └── schc.cfg          # strict SCHC (LT=0)
├── include/
│   ├── PhysicsParams.h   # tunable amplitude parameters + file I/O
│   ├── SDMEModel.h       # params → helicity matrix u[ν'][ν][μ'][μ]
│   ├── generator.h       # SDMEGenerator class (replaces original generator.h)
│   ├── w_kernels.hpp     # angular kernel calculators W^{LL,LT,TT}
│   ├── reaction.h        # reaction kinematics
│   ├── cross.h           # cross section utilities
│   ├── decay2body.h      # 2-body decay kinematics
│   └── fizika.h          # Lorentz vector library
└── src/
    └── main.cpp          # entry point
```

---

## Fitting workflow (planned)

1. Generate MC with trial parameters → pass through FastMC acceptance
2. Compare acceptance-corrected angular distributions to data
3. Extract SDMEs by moment analysis or χ² minimisation over (A_L, A_T, A_LT, Im_LL, Im_TT, A_LTi)
4. Use bootstrap replicas for uncertainty propagation
5. Feed extracted σ_L into dipole fit for gluon radius

The parameters `Im_LL`, `Im_TT`, and `A_LTi` directly control the three BSA contributions
and can be constrained simultaneously with the cross-section fits using the A_LU measurement.
# SDMEs_fitter
