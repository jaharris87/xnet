# External serial CPU library tests

Build the production archive and network preprocessor, then compile and run
the callers through the build's exported compiler/link settings:

```sh
make BUILD_NAME=library-test -j libxnet net_setup
make -C test/library XNET_CONFIG=../../build/library-test/xnet-config.mk test
```

Use a different caller `BUILD_DIR` when checking another production
configuration, for example `BUILD_DIR=build-debug` for a debug archive.

Python 3 is required; there are no additional Python packages. The production
Starkiller EOS table must be present at its normal repository path. The tests
use the real dense serial CPU archive and selected EOS; they provide no
replacement solver, EOS, or provider routines. The first configuration is
GNU Fortran with `EOS=STARKILLER`, `MATRIX_SOLVER=dense`, and MPI, OpenMP, and
GPU disabled. Other configurations require their own evidence.

`api_consumer.F90` imports only `xnet_api`, standard intrinsic modules, and a
test fixture that contains no XNet imports. It links against `libxnet.a` and
the exported external libraries; no loose production objects are linked.
`legacy_reference.F90` is a separate process linked to the same archive. It
manually initializes the preexisting module state and calls `full_net`
directly. It does not use the public API or its setup/allocation helpers.
It starts the thermodynamic-history index at one, as the current
`t9rhofind`/`timestep` routines require for a two-point history. Starting at
zero would allow `timestep` to read element zero of a one-based history.
For fixed-temperature calls, both callers preserve the host's temperature
in kelvin rather than round-trip it through the internal GK representation.

The runner copies the public `Data_alpha` and `Data_SN160` inputs, preprocesses
those copies with `net_setup`, and stages an EOS-table symlink. All generated
data, diagnostics, and logs remain in a new directory under `build/runs/`.
Repeated runs preserve prior output. Tracked networks and references are not
modified.

The physical input is C/O at density `1e7 g/cm^3`, temperature `1.7` through
`1.85 GK`, and duration `0.1 s` per call. These conditions and integration
controls follow the maintained self-heating cases; the synthetic exact
composition avoids the rounding in historical molar-abundance files. The
first calls have `X(c12)=X(o16)=0.475` and a fixed auxiliary reservoir with
`X=0.05`, `A=56`, `Z=26`. Their input Ye includes that reservoir. A separate
call omits all auxiliary arguments and uses `X(c12)=X(o16)=0.5`, `Ye=0.5`.

Five cases exercise alpha BE with and without self-heating, alpha BDF without
self-heating, and SN160 BE/BDF with self-heating and weak-reaction neutrino
loss. Each runs a four-lane mask containing a hole, repeats the burn using
its outputs, calls a short three-lane batch, and calls an all-inactive batch.
It then omits optional auxiliaries and uses four-component diagnostic zone
IDs; the earlier calls use three-component IDs. Positive and zero restart
hints are supplied. The public caller finalizes twice, checks burns outside
the initialized lifetime, reinitializes, and repeats its first burn.

The comparison includes every species mass fraction, temperature in kelvin,
Ye, integrated neutrino loss in erg/g, the returned `tdel_old` in seconds, and
the five solver counters. The same serial archive and equivalent arithmetic
are expected to produce identical binary64 values; text output has enough
decimal digits for those values to round-trip. This compares wrapper
behavior with existing `full_net`; it is not an independent scientific
reference. No existing tolerances or reference values are changed. Timing
values are checked for finite nonnegative values, not compared across
processes. Inactive caller state must remain unchanged and its loss/counters
must be zero. Metadata names and order are compared with the staged `sunet`;
all A/Z and energy values are compared with the legacy nuclear-data arrays.

Validation checks cover invalid controls, duplicate initialization, unsupported
time-series output/preprocessing, invalid extents, incomplete auxiliary
arguments, NaN/Inf input, negative fractions/restart hints, invalid Ye, and an
auxiliary reservoir with invalid A. Recoverable validation is distinct from
existing fatal reader/EOS/solver paths; these tests do not claim recoverable
convergence failures, concurrent instances, MPI, OpenMP, GPU, or performance.
