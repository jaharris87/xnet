# XNet v9: goals and remaining work

XNet v9 should make the reaction network easier to use both as a standalone
program and from CHIMERA and Flash-X. The application-facing work should let
callers initialize the network once, submit repeated burns, obtain results,
and release XNet resources without reproducing the standalone driver or
manipulating solver workspace.

Preserve numerical behavior and CPU/GPU execution while doing this work.
Measure performance before and after the interface changes. No particular
speedup is required for the release.

[Issue #126](https://github.com/jaharris87/XNet/issues/126) coordinates the
release. This page explains the sequence; the component issues record their
implementation requirements, results, and maintainer decisions.

## Where we are

The starting CPU/GPU measurements use source and inputs from
`3825c705b26f11ea67f43d6098452a3d53379dfb`. Later development does not change
that comparison source. See [#126](https://github.com/jaharris87/XNet/issues/126)
for current implementation status.

| Work | Current result | Next step |
| --- | --- | --- |
| [#128: configuration](https://github.com/jaharris87/XNet/issues/128) | Merged: validated `xnet_controls_t`, programmatic configuration, layered standalone namelists, and migrated maintained inputs | Use the existing representation in the application interface |
| [#127: starting CPU measurements](https://github.com/jaharris87/XNet/issues/127) | Complete and accepted: serial dense CPU measurements for four network sizes on both facilities, with MA48 deferred | Use these results in #131's CPU/GPU comparison |
| [#129: application interface and library](https://github.com/jaharris87/XNet/issues/129) | Initial dense serial CPU API/archive and default-EOS equivalence tests in [PR #158](https://github.com/jaharris87/XNet/pull/158); application validation remains open | Return exhausted-burn errors for host-controlled fatal shutdown; validate CHIMERA table EOS |
| [#130: application integrations](https://github.com/jaharris87/XNet/issues/130) | Integration validation remains open | CHIMERA first, then Flash-X, with application-side adapters |
| [#131: CPU/GPU and application performance](https://github.com/jaharris87/XNet/issues/131) | Comparable GPU timings and later coupled measurements remain open | Stage A: standalone GPU comparisons now; stage B: resulting-code and repeated application burns after #129/#130 |
| [#132: release checks](https://github.com/jaharris87/XNet/issues/132) | Not yet a release candidate | Verify the final implementation, supported platforms, and documentation |

The CPU results are in the
[#127 measurement summary](https://github.com/jaharris87/XNet/issues/127#issuecomment-6371001835).
They are not a CPU/GPU comparison. The Frontier GPU smoke in
[#149](https://github.com/jaharris87/XNet/pull/149) and the Perlmutter results in
[#152](https://github.com/jaharris87/XNet/issues/152#issuecomment-5938123533)
establish execution for their named configurations and small workloads.
They do not supply the missing network-size GPU timings.

## 1. Finish the standalone CPU/GPU starting measurements

Keep the completed CPU collection in #127. Its reduced completion scope was
approved separately. The maintainer accepted its results and MA48 deferral
and authorized closure on 2026-10-23. The GPU results are not an additional
condition for that completed CPU collection.

Collect the matching GPU measurements as stage A of #131 before integrating
interface changes or making before/after performance claims. An isolated
unmerged #129 implementation may proceed while measurements run from the
pinned source below. Stage A does not depend on #129 or on coupled captures
from #130. The later coupled measurements are stage B of #131.

Use source and inputs from development commit
`3825c705b26f11ea67f43d6098452a3d53379dfb` for both facilities. This is the
[maintainer-selected source](https://github.com/jaharris87/XNet/issues/127#issuecomment-5938328100)
for the current baseline. Preserve the existing CPU results. If a later source
change becomes necessary, identify the change and determine which comparisons
need repeating; do not silently compare different source revisions.

The historical source `86e867c2a64267a674ce4fbf6a3064af39e2f4e0` remains part
of earlier records. It is not the source of the current CPU measurements.
In particular, the Perlmutter historical GPU failure and the separately
labeled adjusted-source success remain as documented in #152/#155.

The minimum remaining measurement uses the same alpha, CCSN52, SN160, and
ECSN350 inputs as the CPU runs: 1024 zones, T9=1.7, density 1.0e8 g cm^-3,
equal C12/O16 mass fractions, Ye=0.5, weak reactions and screening enabled,
self-heating off, and a 10-second integration interval. This is an early-burn
test, not equilibrium or a coupled timestep.

Run one rank on one GPU on each facility, with batch size 16, a point in the
documented `1,4,16,64,128` set. This is a starting comparison, not a claim of
optimal batch size. Use one setup/warm-up execution followed by five measured
repetitions, with the existing numerical
comparison after each timed process. Frontier uses its HIP/ROCm and OpenMP
offload build; Perlmutter uses its CUDA/OpenACC build. Confirm the installed
compiler/runtime combination and request Slurm resources before execution.

Report wall times, XNet Setup/Total/TimeStep and component timers, solver
counters, numerical results, build options, compiler/modules, launcher,
CPU affinity, and device visibility. Each standalone repetition starts a new
process and includes initialization; its internal burn timers do not measure
an already-initialized application call. Compare CPU and GPU on each machine,
and explain the compiler/library and batch-size differences. Cross-machine
results are not an isolated comparison of OpenMP and OpenACC.

MA48 was not configured for the CPU collection. A licensed-source location
was not established on either facility. The maintainer accepted sparse
measurement deferral for #127; do not infer facility-wide absence
from an unset variable or missing module. Never redistribute HSL source.

Finish stage A of #131 with compact CPU/GPU tables, raw repetitions, numerical
results, and limitations, linked from #126 and #127. Store large operational
files externally. GPU batch-size sweeps, ranks/GPU, OpenMP scaling,
self-heating studies, and optimization work
are not prerequisites for these starting measurements.

## 2. Implement the application interface from real calls

#129 owns the XNet procedures, types, and conventional Make library target.
#130 owns the application adapters and representative integration checks.
Inspect both consumers before fixing the shared interface. The maintainer
selected CHIMERA to drive the first implementation, followed by Flash-X.

Start with the existing validated configuration and trace what the caller
needs at initialization, for each burn, and at finalization. Keep species
mapping, application layouts, and energy policy in the application adapter.
Use ordinary Fortran procedures and the existing compile-time provider
practices. Keep one XNet instance and document concurrency restrictions.

The previous #129 signatures and ownership proposal are design material,
not approval to implement them unchanged. The maintainer approved a smaller
[initialize/burn/finalize direction](https://github.com/jaharris87/XNet/issues/129#issuecomment-6372794913):
constant-condition batches, temperature in K, XNet-order mass fractions,
active lanes, auxiliary composition, and the existing penultimate-timestep
meaning. Use CHIMERA table EOS for the first application integration; tests
with the standalone EOS do not establish that integration. The maintainer
[clarified the failure policy](https://github.com/jaharris87/XNet/issues/129#issuecomment-6374869197):
exhaustion after XNet's internal adaptive retries is fatal. The embedded API
must report it so the host can perform its own shutdown; standalone XNet keeps
its existing termination behavior. Application retry, rollback and
successful-lane continuation are not required. Fallback to surrogate or
alternative integration methods is possible later work, not part of this
implementation. Keep the first slice's existing fatal paths documented.
MPI/OpenMP and caller-owned GPU validation remain required before claiming
those modes.

Review detailed signatures and MPI/EOS/GPU ownership in the implementation.
Bring consequential changes to these decisions back to the maintainer;
private helper routines do not each need approval.

The first useful implementation should build a driver-free library, exercise
repeated burns through its public procedures, preserve reusable setup, and
compare the results with the existing numerical path. Continue to test the
standalone executable. Work in small PRs with independent review; a minimal
external caller is useful but does not replace a real application check.

## 3. Validate both consumers and measure repeated use

Under #130, establish a representative CPU burn sequence for each application.
Check species ordering, units, normalization, masks and partial batches,
prior timesteps, counters, self-heating, neutrino losses, and the application's
energy responsibilities where applicable. Check repeated initialization/use/
finalization and ownership. Test each GPU integration mode actually claimed;
standalone GPU success does not establish caller-owned GPU support.

If access to an application or its workload is unavailable, report the
specific limitation and ask the maintainer how it affects release scope.
Do not invent a representative workload or waive a required consumer.

Stage B of #131 owns performance of the resulting implementation: reuse #127
and stage A inputs for before/after standalone comparisons, and measure
repeated calls from reviewed consumer samples. Separate one-time setup,
per-call work, and batch evolution. Choose additional OpenMP, GPU batch,
ranks/GPU, MA48, and self-heating points
to answer concrete questions. The full documented matrix is not a Cartesian
product or an automatic requirement to run every configuration.

Use measurements to decide whether an optimization is worthwhile. A proposed
optimization, extra timing package, or provider layer needs a demonstrated
need. Keep unrelated investigations, including #138's late-time C/O question
and #139's optional HDF5 input support, separate unless the maintainer makes
one necessary for a real consumer.

## 4. Verify and document the release

#132 checks the assembled release candidate: component and numerical
regressions, applicable MPI/OpenMP behavior, configuration migration,
library/interface behavior, both integrations, and supported Frontier and
Perlmutter accelerator configurations. State unsupported combinations and
accepted limitations. Older successful runs do not automatically cover code
changed afterward.

The product documentation should explain how to build and run standalone
XNet, migrate inputs, link the library, call the public procedures, interpret
status/results, and respect ownership and platform limits. Keep measured
results and reproduction commands available alongside those claims.

The maintainer approves scientific/numerical policy, public interface and
support decisions, PR merges, and the final release. Use fresh-context review
as described in [the review playbook](review-playbook.md). Keep the release
plan in #126 and update it when a milestone completes or a decision changes;
avoid creating another hierarchy of planning issues.
