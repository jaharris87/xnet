#!/usr/bin/env python3
"""Run real archive consumers and compare with the same-revision legacy caller."""

from __future__ import annotations

import argparse
import math
import shutil
import subprocess
import tempfile
from pathlib import Path


def run(command: list[str], cwd: Path, log: Path, stdin: str | None = None) -> str:
    result = subprocess.run(
        command, cwd=cwd, input=stdin, text=True, stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT, timeout=180, check=False,
    )
    log.write_text(result.stdout)
    if result.returncode:
        raise RuntimeError(f"{command[0]} exited {result.returncode}; see {log}")
    return result.stdout


def states(output: str) -> dict[str, list[list[float | int]]]:
    result: dict[str, list[list[float | int]]] = {}
    extents: dict[str, tuple[int, int]] = {}
    current: list[list[float | int]] | None = None
    expected_species = expected_lanes = 0
    for line in output.splitlines():
        parts = line.split()
        if not parts:
            continue
        if parts[0] == "STATE":
            if len(parts) != 4 or parts[1] in result:
                raise RuntimeError(f"invalid state header: {line}")
            expected_species, expected_lanes = map(int, parts[2:])
            extents[parts[1]] = (expected_species, expected_lanes)
            current = []
            result[parts[1]] = current
        elif parts[0] == "ZONE":
            if current is None or len(parts) != 11:
                raise RuntimeError(f"invalid zone record: {line}")
            zone = int(parts[1])
            if zone != len(current) + 1:
                raise RuntimeError(f"out-of-order zone: {line}")
            current.append([*map(float, parts[2:6]), *map(int, parts[6:])])
        elif parts[0] == "X":
            if current is None or len(parts) != expected_species + 2:
                raise RuntimeError(f"invalid composition record: {line}")
            zone = int(parts[1])
            if zone != len(current) or zone > expected_lanes:
                raise RuntimeError(f"out-of-order composition: {line}")
            current[-1].extend(map(float, parts[2:]))
    for label, rows in result.items():
        species, lanes = extents[label]
        if len(rows) != lanes or any(len(row) != species + 9 for row in rows):
            raise RuntimeError(f"incomplete state {label}")
        if not rows or any(not all(math.isfinite(v) for v in row) for row in rows):
            raise RuntimeError(f"missing or nonfinite state {label}")
    return result


def check_metadata(output: str, names: list[str]) -> None:
    metadata = [line.split()[1:] for line in output.splitlines() if line.startswith("META ")]
    if [row[0] for row in metadata] != names:
        raise RuntimeError("public metadata names/order differ from sunet")
    for row in metadata:
        if len(row) != 5 or not all(math.isfinite(float(value)) for value in row[1:]):
            raise RuntimeError(f"invalid public metadata: {row}")
        a, z = map(float, row[1:3])
        if not (a > 0 and 0 <= z <= a):
            raise RuntimeError(f"invalid public A/Z: {row}")


def compare(api: dict, reference: dict, label: str) -> None:
    if label not in api or label not in reference or len(api[label]) != len(reference[label]):
        raise RuntimeError(f"missing/mismatched {label} state")
    for lane, (observed, expected) in enumerate(zip(api[label], reference[label]), start=1):
        if len(observed) != len(expected):
            raise RuntimeError(f"{label} lane {lane}: state extent differs")
        for column, (actual, legacy) in enumerate(zip(observed, expected), start=1):
            # Both paths call the same serial archive with the same data and arithmetic.
            # Full-precision decimal output round-trips each binary64 value.
            if actual != legacy:
                raise RuntimeError(
                    f"{label} lane {lane} column {column}: API {actual!r} != legacy {legacy!r}"
                )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--api", type=Path, required=True)
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument("--net-setup", type=Path, required=True)
    parser.add_argument("--work-dir", type=Path, default=Path("build/runs"))
    args = parser.parse_args()
    executables = [p.resolve(strict=True) for p in (args.api, args.reference, args.net_setup)]
    api_exe, reference_exe, setup_exe = executables
    root = Path(__file__).resolve().parents[2]
    table = root / "tools/starkiller-helmholtz/helm_table.dat"
    if not table.is_file():
        raise RuntimeError(f"missing production EOS table: {table}")
    args.work_dir.mkdir(parents=True, exist_ok=True)
    work = Path(tempfile.mkdtemp(prefix="library-", dir=args.work_dir.resolve()))
    print(f"Library consumer work: {work}", flush=True)
    networks = {}
    for network in ("alpha", "SN160"):
        data_name = f"Data_{network}"
        source = root / "test" / data_name
        destination = work / data_name
        shutil.copytree(source, destination, symlinks=False)
        run([str(setup_exe)], destination, work / f"net_setup_{network}.log",
            stdin=f"Public {network} embedding-library tests\n")
        for required in ("sunet", "netwinv", "nets3", "nets4", "match_data"):
            if not (destination / required).is_file():
                raise RuntimeError(f"net_setup did not produce {destination / required}")
        networks[network] = destination
    for network, solver, heat in (
        ("alpha", 1, 0), ("alpha", 1, 1), ("alpha", 3, 0),
        ("SN160", 1, 1), ("SN160", 3, 1),
    ):
        label = f"{network}-{'be' if solver == 1 else 'bdf'}-heat{heat}"
        case = work / label
        case.mkdir()
        data = networks[network]
        (case / data.name).symlink_to(data, target_is_directory=True)
        (case / table.name).symlink_to(table)
        command_args = [data.name, str(solver), str(heat)]
        public_output = run([str(api_exe), *command_args], case, case / "api.log")
        if "API_CHECKS_PASSED" not in public_output:
            raise RuntimeError(f"public caller did not complete checks; see {case / 'api.log'}")
        legacy_output = run([str(reference_exe), *command_args], case, case / "legacy.log")
        check_metadata(public_output, (data / "sunet").read_text().split())
        public_metadata = [line for line in public_output.splitlines() if line.startswith("META ")]
        legacy_metadata = [line for line in legacy_output.splitlines() if line.startswith("META ")]
        if public_metadata != legacy_metadata:
            raise RuntimeError(f"public A/Z or energy metadata differ from legacy nuclear data: {label}")
        public_states, legacy_states = states(public_output), states(legacy_output)
        timings = [line.split() for line in public_output.splitlines() if line.startswith("TIMING ")]
        if len(timings) != 6 or any(
            len(row) != 8 or any(not math.isfinite(float(v)) or float(v) < 0 for v in row[2:])
            for row in timings
        ):
            raise RuntimeError(f"invalid per-call timings; see {case / 'api.log'}")
        if set(public_states) != {"first", "repeat", "short", "inactive", "no_aux", "reinitialize"}:
            raise RuntimeError(f"missing public burn calls; see {case / 'api.log'}")
        for state in ("first", "repeat", "short", "inactive", "no_aux"):
            compare(public_states, legacy_states, state)
        compare(public_states, {"reinitialize": public_states["first"]}, "reinitialize")
        first = public_states["first"][0]
        repeated = public_states["repeat"][0]
        if first[9:] == repeated[9:]:
            raise RuntimeError(f"{label}: successive burns did not evolve composition")
        if heat and not first[0] > 1.7e9:
            raise RuntimeError(f"{label}: self-heating did not increase temperature")
        if network == "SN160" and not first[3] > 0:
            raise RuntimeError(f"{label}: weak-reaction neutrino loss not exercised")
        print(f"PASS {label}: five burn calls equal direct full_net; lifecycle/validation/metadata pass", flush=True)
    print("Real serial CPU library consumer checks passed", flush=True)


if __name__ == "__main__":
    main()
