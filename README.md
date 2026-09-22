# ACCESS-CM2 mixed-layer heat budget: 1998a foundation and 2001-style decomposition

## Scope

This directory calculates a native-grid ACCESS-CM2 mixed-layer heat budget in
`W m-2`. The physical budget follows the time-varying mixed-layer framework of
Vialard and Delecluse (1998a). It additionally decomposes monthly advection
into climatological-mean and monthly-anomaly components, motivated by Vialard
et al. (2001).

The available input is monthly (`Omon`). Therefore the decomposition is a
**2001-inspired monthly approximation**, not the exact within-month eddy
covariance diagnosed from high-frequency or online model tendencies in the
2001 study.

## References and the important difference

1. Vialard, J., and P. Delecluse (1998a), *An OGCM Study for the TOGA Decade.
   Part I: Role of Salinity in the Physics of the Western Pacific Fresh Pool*,
   Journal of Physical Oceanography, 28, 1071-1088.
   DOI: https://doi.org/10.1175/1520-0485(1998)028%3C1071:AOSFTT%3E2.0.CO;2

2. Vialard, J., C. Menkes, J.-P. Boulanger, P. Delecluse, E. Guilyardi,
   M. J. McPhaden, and G. Madec (2001), *A Model Study of Oceanic Mechanisms
   Affecting Equatorial Pacific Sea Surface Temperature during the 1997-98 El
   Nino*, Journal of Physical Oceanography, 31, 1649-1675.
   DOI: https://doi.org/10.1175/1520-0485(2001)031%3C1649:AMSOO%3E2.0.CO;2

1998a provides the time-varying mixed-layer heat-budget foundation and focuses
on salinity, barrier layers, and the western Pacific warm pool. The 2001 study
uses that framework to investigate the 1997-98 El Nino and separates
low-frequency advection from eddy effects such as tropical instability waves.
Kelvin waves are interpreted through the existing budget terms; they are not a
separate additive term in the heat equation.

The exact 2001 eddy effect requires model-timestep/high-frequency tendencies
before monthly averaging. Monthly mean velocity and temperature alone cannot
recover the true within-month covariance.

## Native inputs and domain

The script reads the top-level native-data links under
`/g/data/p66/ars599/work_budget/data`, selected with token `bj594_piControl`.
It does not use `data/cm2/`.

Required fields are `temp`, `mld`, `u`, `v`, `wt`, `sw_heat`,
`sfc_hflux_coupler`, `sfc_hflux_pme`, and `frazil_3d`. The domain is 120-280E
and 25S-25N; actual native tracer centres extend to about 24.596 degrees.
Temperature is converted from the Kelvin-valued CM2 source to degrees Celsius.
MLD is capped at 50 m, and source levels are read through 60 m.

## 1998a-style raw budget

For layer overlap thickness `Delta z_k` inside the time-varying MLD `h`:

```text
T_MIX = sum(T_k Delta z_k) / h
TEND  = rho0 Cp h d(T_MIX)/dt
```

The derivative uses decoded real time in seconds. The retained surface forcing
is

```text
SHF = Q_coupler + Q_pme + Q_frazil - Q_shortwave_below_MLD
```

Horizontal advection is

```text
AX = -rho0 Cp sum[u (dT/dx) Delta z]
AY = -rho0 Cp sum[v (dT/dy) Delta z]
```

The reconstructed bottom terms are

```text
AZ = rho0 Cp w_b (T_b - T_MIX)
E1 = rho0 Cp max(dh/dt, 0) (T_b - T_MIX)
E2 = rho0 Cp Kv (dT/dz)_b
```

`E1` and constant-`Kv` `E2` are offline approximations, not exact online model
tendencies. The raw closure diagnostics are

```text
ADV      = AX + AY + AZ
RHS      = SHF + AX + AY + AZ + E1 + E2
RESIDUAL = TEND - RHS
```

Constants are `rho0=1035 kg m-3`, `Cp=3989.24495 J kg-1 K-1`, and default
`Kv=1e-5 m2 s-1`.

## 2001-inspired monthly decomposition

For each field, the bar is its calendar-month climatology over the selected
calculation period and the prime is the monthly departure:

```text
T = Tbar + Tprime
u = ubar + uprime
```

The first suffix letter is temperature and the second is velocity:

```text
BB = -ubar   dot grad(Tbar)
AB = -ubar   dot grad(Tprime)
BA = -uprime dot grad(Tbar)
AA = -uprime dot grad(Tprime)
```

The code applies this decomposition separately to `AX`, `AY`, and the
reconstructed `AZ`, and writes `ADV_BB`, `ADV_AB`, `ADV_BA`, and `ADV_AA`.
It also writes

```text
ADV_DECOMP_SUM   = ADV_BB + ADV_AB + ADV_BA + ADV_AA
ADV_DECOMP_ERROR = ADV - ADV_DECOMP_SUM
```

A near-zero decomposition error verifies algebraic closure. `AA` here is the
product of monthly anomalies. It must not be described as the exact
within-month eddy tendency of Vialard et al. (2001).

## Output variables

Raw variables: `T_MIX`, `H_EXTN`, `TEND`, `SHF`, `AX`, `AY`, `AZ`, `ADV`,
`E1`, `E2`, `RHS`, and `RESIDUAL`.

Decomposition variables: `AX_*`, `AY_*`, `AZ_*`, and `ADV_*`, where `*` is
`BB`, `AB`, `BA`, or `AA`, plus `ADV_DECOMP_SUM` and `ADV_DECOMP_ERROR`.
All budget variables use `W m-2`.

## Run

Validate inputs:

```bash
python prepare_cm2_data.py
```

Smoke test with at least 24 months so every calendar month has more than one
sample:

```bash
qsub -v NMONTHS=24,OUTPUT=cm2_vd2001_smoke_24months.nc run_hba_cm2.pbs
```

Full calculation:

```bash
qsub run_hba_cm2.pbs
```

Default output: `ACCESS_CM2_HBA_VD2001_native_25S25N_Wm2.nc`.
