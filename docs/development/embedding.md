# Constant-condition Fortran embedding

Build the archive and compiler-specific module files with:

```sh
make BUILD_NAME=embedding -j libxnet
```

The archive is `build/embedding/lib/libxnet.a`; modules are in
`build/embedding/mod`. The generated `build/embedding/xnet-config.mk`
records the compiler, module path, compile/link flags, archive, and required
external libraries. Use the same compiler and compatible precision/runtime
options for the application. The external example under `test/library/`
links through these outputs without compiling XNet sources.

The supported initial configuration is a dense serial CPU build with MPI,
OpenMP, and GPU modes off. Other builds retain their standalone programs;
the embedding procedures return `xnet_unsupported`. This does not qualify
CHIMERA/Flash-X adapters, host EOS providers, or caller-owned GPU resources.

## Three procedures

```fortran
Use xnet_api, Only: dp, xnet_controls_t, set_xnet_controls_defaults, &
  & xnet_network_t, xnet_timing_t, xnet_initialize, xnet_burn_batch, xnet_finalize

Call set_xnet_controls_defaults(controls)
controls%data_dir = 'prepared_network'
controls%nzbatchmx = 16
controls%idiag = -2
Call xnet_initialize(controls,network,ierr,message)
! Check ierr before using network or burning.
Call xnet_burn_batch(dt_s,active,zone_id,density,temperature,ye,mass_fraction, &
  & restart_dt_s,neutrino_loss,counters,timing,ierr,message,aux_x,aux_a,aux_z)
Call xnet_finalize(ierr,message)
```

Declare `controls`, `network`, and `timing` using the corresponding types;
the other arguments are ordinary caller-owned arrays. The library owns one
instance. Calls do not overlap or recurse. Initialize once, burn any number
of batches, and finalize once; finalization is also safe when already
finalized. The library keeps no association with caller arrays.

Use preprocessed network files (`iprocess=0`), produced separately by
`net_setup` when necessary. The constant-condition API has no incident
neutrino histories (`ineutrino=0`) or timestep output files (`itsout<=0`).
Diagnostic output, if requested, goes to the existing standard-output unit.
The compiled default controls need no standalone trajectory/abundance
filenames. The selected EOS can require additional inputs; the default
Starkiller provider reads `helm_table.dat` from the working directory.

On successful initialization, `network` contains copied names, integer A/Z,
binding energies, and mass excesses from the selected network. Both energy
arrays use the MeV values in `nuclear_data`. The caller owns this metadata;
it remains usable after finalization. Failed initialization does not produce
usable metadata.

## Batch arguments

The duration `dt_s` is a positive scalar in seconds. Density is fixed in
g/cm^3. Temperature is in/out K; it evolves when `iheat>0`. `ye` is in/out,
following the existing composition/EOS calculation. This path supplies two
identical thermodynamic endpoints internally and has no absolute time origin.

`active(:)` defines the submitted lane count, at most `nzbatchmx`.
Thermodynamic, restart, and neutrino-loss arrays have that length.
`mass_fraction(species,lane)` follows XNet species order and contains
dimensionless mass fractions. The result is raw: the API does not normalize
or reorder composition or compute application nuclear-energy release.
Applications retain those operations and their EOS/energy policy.

`zone_id(:,lane)` has three or four integer components for diagnostics.
Existing XNet failure filenames use the first three components. The planned
error return will let each host report its own identifiers.
Short batches and inactive holes are supported. Inactive lanes preserve
temperature, Ye, composition and restart hint, and return zero loss/counters.
Scientific values in inactive lanes are not read. An all-inactive batch
performs no burn.

Optional `aux_x(:)`, `aux_a(:)`, `aux_z(:)` must be supplied together.
They describe a fixed auxiliary mass reservoir, its mass number and charge;
it contributes to existing composition/EOS moments without becoming an
evolved network species. Positive auxiliary mass requires positive A and
`0<=Z<=A`. Zero auxiliary mass is treated as no reservoir.

`restart_dt_s(:)` is nonnegative and in/out. Zero requests fresh startup.
Positive input preserves the existing CHIMERA hint convention: the current
startup uses half the hint as an initial-step floor, rather than imposing
that timestep. Output is existing `tdel_old`, the preceding/penultimate
internal timestep; it can be zero after a single step. It is not a next-step
recommendation. The CHIMERA adapter currently supplies zero for BDF.

`neutrino_loss(:)` returns the integrated XNet weak-reaction loss in erg/g.
It is distinct from incident-neutrino capture and application neutrino-loss
corrections. `counters(5,lane)` returns the existing TS, NR, Jacobian,
derivative, and rate counters. `timing` returns per-batch elapsed, solve,
Jacobian, derivative, rate, and screening seconds; applications accumulate
their own totals.

## Status and remaining limits

`xnet_success=0`, `xnet_invalid=1`, `xnet_lifecycle_error=2`, and
`xnet_unsupported=3` describe the checked public paths. Validate configuration,
dimensions, finite active inputs and lifecycle order before a burn. Rejected
burn calls preserve in/out physical state and return zero loss/counters.
No automatic renormalization or repair is performed.

Status does not make deep failures recoverable. Existing malformed-network
readers, EOS/provider failures, allocation failures, and solver exhaustion can
still terminate execution. `full_net` can write failure inputs before
termination. A required follow-up will report exhausted-burn errors to the
host for its own fatal shutdown, without changing standalone termination.
Application retry, rollback and successful-lane continuation are not supported
or required for this interface. Future fallback to another integration method
is separate work, not a promise that a failed library instance can be reused.

Finalization releases CPU allocations owned by XNet and its selected
standalone EOS provider. It never initializes/finalizes MPI or a GPU runtime.
A future application EOS provider must release only its XNet scratch storage,
leaving the host EOS initialized. MPI, OpenMP, and caller-owned GPU embedding
remain unqualified.
