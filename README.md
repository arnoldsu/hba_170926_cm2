# HBA-CM2: ACCESS-CM2 Mixed-Layer Heat Budget and Monthly Advection Decomposition (MOM direct output)

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.22908822.svg)](https://doi.org/10.5281/zenodo.22908822)

**DOI:** [10.5281/zenodo.22908822](https://doi.org/10.5281/zenodo.22908822)

## Scope

This repository calculates a native-grid ACCESS-CM2 mixed-layer heat budget
in `W m-2`.

The physical budget follows the time-varying mixed-layer framework of
Vialard and Delecluse (1998a). It additionally decomposes monthly ocean
advection into climatological-mean and monthly-anomaly components, motivated
by the decomposition used by Vialard et al. (2001).

The available ACCESS-CM2 input is monthly (`Omon`). Therefore, the advection
decomposition implemented here is a **Vialard et al. (2001)-inspired monthly
approximation**, rather than the exact within-month eddy covariance that would
require high-frequency or online model tendencies.

---

## Mixed-Layer Heat-Budget Framework

For a time-varying mixed-layer depth `h`, the mixed-layer temperature is
calculated from the thickness-weighted temperature of all model layers
intersecting the mixed layer:

```math
T_{\mathrm{MIX}}
=
\frac{1}{h}
\sum_k
T_k \Delta z_k
```

where `Delta z_k` is the thickness of model layer `k` contained within the
mixed layer.

The mixed-layer heat-content tendency is

```math
\mathrm{TEND}
=
\rho_0 C_p h
\frac{\partial T_{\mathrm{MIX}}}{\partial t}
```

and the diagnosed budget is written as

```math
\mathrm{TEND}
=
\mathrm{SHF}
+
\mathrm{AX}
+
\mathrm{AY}
+
\mathrm{AZ}
+
\mathrm{E1}
+
\mathrm{E2}
+
\mathrm{RESIDUAL}
```

where:

- `SHF` is the retained surface heat-flux forcing,
- `AX` is zonal temperature advection,
- `AY` is meridional temperature advection,
- `AZ` is reconstructed vertical temperature advection,
- `E1` is the reconstructed mixed-layer entrainment contribution,
- `E2` is an offline approximation of vertical diffusion at the mixed-layer base,
- `RESIDUAL` contains unresolved processes, approximation errors, and budget
  closure error.

The resolved three-dimensional advection is

```math
\mathrm{ADV}
=
\mathrm{AX}
+
\mathrm{AY}
+
\mathrm{AZ}
```

and the diagnosed right-hand side is

```math
\mathrm{RHS}
=
\mathrm{SHF}
+
\mathrm{ADV}
+
\mathrm{E1}
+
\mathrm{E2}
```

with

```math
\mathrm{RESIDUAL}
=
\mathrm{TEND}
-
\mathrm{RHS}
```

---

## Native Inputs and Domain

The script reads the top-level native-data links under

```text
/g/data/p66/ars599/work_budget/data
```

selected using the token

```text
bj594_piControl
```

It does **not** use `data/cm2/`.

Required fields are

```text
temp
mld
u
v
wt
sw_heat
sfc_hflux_coupler
sfc_hflux_pme
frazil_3d
```

The analysis domain is

```text
Longitude: 120E-280E
Latitude:   25S-25N
```

The actual native tracer centres extend to approximately 24.596 degrees
latitude.

Temperature is converted from the Kelvin-valued ACCESS-CM2 source to degrees
Celsius.

The mixed-layer depth is capped at

```text
50 m
```

and source temperature and velocity levels are read through

```text
60 m
```

to provide the information required immediately below the diagnosed mixed
layer.

---

## 1998a-Style Raw Budget

### Mixed-Layer Temperature and Tendency

For layer overlap thickness `Delta z_k` inside the time-varying mixed-layer
depth `h`,

```text
T_MIX = sum(T_k * Delta_z_k) / h

TEND = rho0 * Cp * h * d(T_MIX)/dt
```

or

```math
T_{\mathrm{MIX}}
=
\frac{
\sum_k T_k \Delta z_k
}{
h
}
```

and

```math
\mathrm{TEND}
=
\rho_0 C_p h
\frac{\partial T_{\mathrm{MIX}}}{\partial t}
```

The derivative uses the decoded model time coordinate converted to seconds,
rather than assuming that every month has identical duration.

---

### Surface Heat Flux

The retained surface forcing is

```text
SHF = Q_coupler + Q_pme + Q_frazil - Q_shortwave_below_MLD
```

or schematically,

```math
\mathrm{SHF}
=
Q_{\mathrm{coupler}}
+
Q_{\mathrm{pme}}
+
Q_{\mathrm{frazil}}
-
Q_{\mathrm{SW,below\,MLD}}
```

The shortwave component penetrating below the mixed-layer base is removed
because it does not directly heat the diagnosed mixed-layer volume.

---

### Horizontal Advection

Zonal temperature advection is

```math
\mathrm{AX}
=
-\rho_0 C_p
\sum_k
u_k
\frac{\partial T_k}{\partial x}
\Delta z_k
```

and meridional temperature advection is

```math
\mathrm{AY}
=
-\rho_0 C_p
\sum_k
v_k
\frac{\partial T_k}{\partial y}
\Delta z_k
```

implemented schematically as

```text
AX = -rho0 * Cp * sum[u * (dT/dx) * Delta_z]

AY = -rho0 * Cp * sum[v * (dT/dy) * Delta_z]
```

---

### Vertical Advection

The reconstructed vertical-advection contribution at the mixed-layer base is

```text
AZ = rho0 * Cp * w_b * (T_b - T_MIX)
```

or

```math
\mathrm{AZ}
=
\rho_0 C_p w_b
\left(
T_b-T_{\mathrm{MIX}}
\right)
```

where `w_b` and `T_b` represent the vertical velocity and temperature
associated with the base of the diagnosed mixed layer.

---

### Entrainment

The reconstructed entrainment contribution is

```text
E1 = rho0 * Cp * max(dh/dt, 0) * (T_b - T_MIX)
```

or

```math
\mathrm{E1}
=
\rho_0 C_p
\max
\left(
\frac{\partial h}{\partial t},
0
\right)
\left(
T_b-T_{\mathrm{MIX}}
\right)
```

Only mixed-layer deepening contributes to this reconstructed entrainment
term through `max(dh/dt, 0)`.

`E1` is an offline reconstruction and should not be interpreted as an exact
online ACCESS-CM2 model tendency.

---

### Vertical Diffusion

The reconstructed vertical-diffusion contribution at the mixed-layer base is

```text
E2 = rho0 * Cp * Kv * (dT/dz)_b
```

or

```math
\mathrm{E2}
=
\rho_0 C_p K_v
\left.
\frac{\partial T}{\partial z}
\right|_b
```

where `Kv` is a prescribed vertical diffusivity.

The default value is

```text
Kv = 1e-5 m2 s-1
```

As with `E1`, this constant-`Kv` `E2` term is an offline approximation rather
than an exact online ACCESS-CM2 vertical-mixing tendency.

---

### Raw Budget Closure

The three-dimensional resolved advection is

```text
ADV = AX + AY + AZ
```

or

```math
\mathrm{ADV}
=
\mathrm{AX}
+
\mathrm{AY}
+
\mathrm{AZ}
```

The diagnosed right-hand side is

```text
RHS = SHF + AX + AY + AZ + E1 + E2
```

or

```math
\mathrm{RHS}
=
\mathrm{SHF}
+
\mathrm{AX}
+
\mathrm{AY}
+
\mathrm{AZ}
+
\mathrm{E1}
+
\mathrm{E2}
```

and the closure residual is

```text
RESIDUAL = TEND - RHS
```

or

```math
\mathrm{RESIDUAL}
=
\mathrm{TEND}
-
\mathrm{RHS}
```

A small residual indicates better closure of the diagnosed offline budget,
but `E1` and `E2` remain reconstructed approximations rather than exact
online model tendencies.

The constants used by the calculation are

```text
rho0 = 1035 kg m-3

Cp   = 3989.24495 J kg-1 K-1

Kv   = 1e-5 m2 s-1
```

---

## 2001-Inspired Monthly Advection Decomposition

For each monthly field, the overbar represents its calendar-month climatology
over the selected calculation period and the prime represents the monthly
departure from that climatology:

```math
T
=
\overline{T}
+
T'
```

and

```math
\mathbf{u}
=
\overline{\mathbf{u}}
+
\mathbf{u}'
```

where `u` here represents the relevant velocity component.

Expanding the advective product gives four components:

```math
-\mathbf{u}\cdot\nabla T
=
-\overline{\mathbf{u}}\cdot\nabla\overline{T}
-\overline{\mathbf{u}}\cdot\nabla T'
-\mathbf{u}'\cdot\nabla\overline{T}
-\mathbf{u}'\cdot\nabla T'
```

The code labels these components according to the temperature field first
and the velocity field second:

```text
BB = -ubar   dot grad(Tbar)

AB = -ubar   dot grad(Tprime)

BA = -uprime dot grad(Tbar)

AA = -uprime dot grad(Tprime)
```

Thus,

```math
\mathrm{BB}
=
-\overline{\mathbf{u}}
\cdot
\nabla\overline{T}
```

```math
\mathrm{AB}
=
-\overline{\mathbf{u}}
\cdot
\nabla T'
```

```math
\mathrm{BA}
=
-\mathbf{u}'
\cdot
\nabla\overline{T}
```

```math
\mathrm{AA}
=
-\mathbf{u}'
\cdot
\nabla T'
```

The first suffix letter identifies the **temperature** component and the
second suffix letter identifies the **velocity** component.

The decomposition is applied separately to `AX`, `AY`, and the reconstructed
`AZ`.

The resulting total-advection components are written as

```text
ADV_BB
ADV_AB
ADV_BA
ADV_AA
```

with

```text
ADV_DECOMP_SUM =
    ADV_BB
  + ADV_AB
  + ADV_BA
  + ADV_AA
```

and

```text
ADV_DECOMP_ERROR = ADV - ADV_DECOMP_SUM
```

or

```math
\mathrm{ADV}_{\mathrm{DECOMP}}
=
\mathrm{ADV}_{BB}
+
\mathrm{ADV}_{AB}
+
\mathrm{ADV}_{BA}
+
\mathrm{ADV}_{AA}
```

with

```math
\mathrm{ADV\_DECOMP\_ERROR}
=
\mathrm{ADV}
-
\mathrm{ADV}_{\mathrm{DECOMP}}
```

A near-zero `ADV_DECOMP_ERROR` verifies the algebraic closure of the
implemented monthly decomposition.

### Important Interpretation of AA

`AA` is the product of **monthly anomalies**:

```math
\mathrm{AA}
=
-\mathbf{u}'
\cdot
\nabla T'
```

It must **not** be interpreted as the exact within-month eddy tendency
diagnosed by Vialard et al. (2001).

The exact eddy contribution requires covariance information at a temporal
resolution higher than the monthly means used here. In general,

```math
\overline{
\mathbf{u}'\cdot\nabla T'
}
```

cannot be reconstructed exactly from monthly-mean velocity and temperature
fields alone.

Therefore, this repository implements a **Vialard et al. (2001)-inspired
monthly decomposition**, not an exact reproduction of the high-frequency
eddy decomposition used in that study.

---

## Relation to Vialard and Delecluse (1998a) and Vialard et al. (2001)

Vialard and Delecluse (1998a) provide the time-varying mixed-layer
heat-budget foundation used here, with particular emphasis on salinity,
barrier-layer processes, and the western Pacific warm pool.

Vialard et al. (2001) apply a related heat-budget framework to investigate
oceanic mechanisms affecting equatorial Pacific SST during the 1997-98
El Niño and separate low-frequency advection from higher-frequency eddy
effects, including tropical instability waves.

The present ACCESS-CM2 implementation combines these ideas:

```text
1998a
  |
  +-- time-varying mixed-layer heat-budget framework
  |
  +-- mixed-layer tendency
  +-- surface forcing
  +-- horizontal advection
  +-- vertical advection
  +-- entrainment / mixing approximations

2001
  |
  +-- mean/anomaly advection decomposition
      |
      +-- BB
      +-- AB
      +-- BA
      +-- AA
```

Kelvin-wave effects are interpreted through their influence on the existing
temperature, velocity, thermocline, and advection terms. They are **not**
introduced as an additional independent term in the heat equation.

---

## Output Variables

### Raw Budget Variables

```text
T_MIX
H_EXTN
TEND
SHF
AX
AY
AZ
ADV
E1
E2
RHS
RESIDUAL
```

### Advection-Decomposition Variables

```text
AX_BB
AX_AB
AX_BA
AX_AA

AY_BB
AY_AB
AY_BA
AY_AA

AZ_BB
AZ_AB
AZ_BA
AZ_AA

ADV_BB
ADV_AB
ADV_BA
ADV_AA

ADV_DECOMP_SUM
ADV_DECOMP_ERROR
```

All heat-budget tendency variables are expressed in

```text
W m-2
```

---

## Running the Analysis

### 1. Validate Inputs

Run

```bash
python prepare_cm2_data.py
```

before submitting the heat-budget calculation.

### 2. Smoke Test

Use at least 24 months so that every calendar month has more than one sample
for constructing the monthly climatology:

```bash
qsub -v NMONTHS=24,OUTPUT=cm2_vd2001_smoke_24months.nc run_hba_cm2.pbs
```

### 3. Full Calculation

Submit the complete calculation with

```bash
qsub run_hba_cm2.pbs
```

The default output file is

```text
ACCESS_CM2_HBA_VD2001_native_25S25N_Wm2.nc
```

---

## Citation

If you use this software or methodology in your research, please cite the
archived release:

**Sullivan, Arnold. (2026).  
JK-CM2: ACCESS-CM2 Mixed-Layer Heat Budget and Monthly Advection Decomposition
(v1.0.0). Zenodo.**

**DOI:** [10.5281/zenodo.22908822](https://doi.org/10.5281/zenodo.22908822)

### BibTeX

```bibtex
@software{Sullivan_2026_HBA_CM2,
  author    = {Sullivan, Arnold},
  title     = {JK-CM2: ACCESS-CM2 Mixed-Layer Heat Budget and Monthly Advection Decomposition},
  year      = {2026},
  version   = {1.0.0},
  publisher = {Zenodo},
  doi       = {10.5281/zenodo.22908822},
  url       = {https://doi.org/10.5281/zenodo.22908822}
}
```

---

## References

Vialard, J., & Delecluse, P. (1998).
*An OGCM Study for the TOGA Decade. Part I: Role of Salinity in the Physics
of the Western Pacific Fresh Pool.*
Journal of Physical Oceanography, **28**, 1071-1088.  
https://doi.org/10.1175/1520-0485(1998)028%3C1071:AOSFTT%3E2.0.CO;2

Vialard, J., Menkes, C., Boulanger, J.-P., Delecluse, P., Guilyardi, E.,
McPhaden, M. J., & Madec, G. (2001).
*A Model Study of Oceanic Mechanisms Affecting Equatorial Pacific Sea Surface
Temperature during the 1997-98 El Niño.*
Journal of Physical Oceanography, **31**, 1649-1675.  
https://doi.org/10.1175/1520-0485(2001)031%3C1649:AMSOO%3E2.0.CO;2
