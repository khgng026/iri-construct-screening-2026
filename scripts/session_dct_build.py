#!/usr/bin/env python3
"""Build pinned RecCut and check frozen graph-only fixtures, on a lead-owned slot.

Standard library at import; existing NumPy/SciPy only for the runtime numeric
quantization gate. No biological sequence, model, installation or upstream source
edits. Every attempt owns a fresh ignored directory and retains failures.
This script has been authored source-only; no compilation/fixtures run by author.
"""
from __future__ import annotations

import argparse
import datetime as dt
from decimal import Decimal, ROUND_FLOOR
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import re
import shlex
import shutil
import signal
import subprocess
import sys
import time
import uuid


ROOT = Path(__file__).resolve().parents[1]
SOURCE_SHA256 = "00f503c38f28368334adb9bbec361d2d765e401264a42b50a0b04ecda562c074"
HEADER_SHA256 = "fb0e59494bf9992252f6633d7d61ac0122299f340650bc54fabf673de7a5db69"
LICENSE_SHA256 = "dbbee091af63c74ab091bf536fe5783c0b9472254be345491ecaaf99bd24231f"
FINGERPRINT_SHA256 = "f0d3232743040e87bf82b9fb3fc218d7fdf409876a73feafdda70905e17f2564"
DESIGN_SHA256 = "7edbb283995290b704cf6d9b640dd2fa452cf73bc6b8c40a7083385fc5a0bb7e"
FLAGS = ["-std=c++17", "-O2", "-Wall", "-Wextra", "-pedantic"]


def require(value, message):
    if not value:
        raise ValueError(message)


def sha256(path):
    value = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def utcnow():
    return dt.datetime.now(dt.timezone.utc).isoformat()


def read_json(path):
    return json.loads(Path(path).read_text())


def atomic_json(path, value):
    """Replace reports only inside this attempt's already exclusively owned root."""
    temporary = path.with_name(path.name + ".tmp-" + uuid.uuid4().hex)
    with temporary.open("x") as stream:
        json.dump(value, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write("\n")
        stream.flush()
        os.fsync(stream.fileno())
    os.replace(temporary, path)


def file_record(path):
    return {"path": str(path), "sha256": sha256(path), "bytes": path.stat().st_size}


def stop_owned_group(process):
    # A compiler driver can exit before an owned descendant. Signal the group,
    # including when the leader has already exited; never target another group.
    try:
        os.killpg(process.pid, signal.SIGTERM)
    except ProcessLookupError:
        return
    try:
        process.wait(timeout=2)
    except subprocess.TimeoutExpired:
        pass
    try:
        os.killpg(process.pid, signal.SIGKILL)
    except ProcessLookupError:
        pass
    if process.poll() is None:
        process.wait(timeout=2)


def bounded_run(command, label, directory, timeout, deadline, report):
    remaining = (deadline - dt.datetime.now(dt.timezone.utc)).total_seconds()
    require(remaining > 0, "Session hard stop reached")
    timeout = min(float(timeout), remaining)
    item = {"label": label, "command": list(map(str, command)), "started_utc": utcnow(),
            "timeout_seconds": timeout, "status": "started"}
    report["commands"].append(item)
    atomic_json(directory / "report.json", report)
    stdout = directory / (label + ".stdout.txt")
    stderr = directory / (label + ".stderr.txt")
    process = None
    started = time.monotonic()
    # Locale and temp selection are scoped to owned children, never system settings.
    environment = dict(os.environ, LC_ALL="C", TMPDIR=str(directory / "tmp"))
    environment.update(report.get("command_environment_overrides", {}))
    try:
        with stdout.open("xb") as out, stderr.open("xb") as err:
            process = subprocess.Popen(item["command"], cwd=directory, env=environment,
                                       stdin=subprocess.DEVNULL, stdout=out, stderr=err,
                                       start_new_session=True)
            item["pid"] = process.pid
            try:
                item["returncode"] = process.wait(timeout=timeout)
                item["status"] = "passed" if item["returncode"] == 0 else "failed"
            except subprocess.TimeoutExpired:
                item["status"] = "timeout"
                stop_owned_group(process)
                item["returncode"] = process.returncode
    except BaseException as error:
        item["status"] = "interrupted" if isinstance(error, (SystemExit, KeyboardInterrupt)) else "error"
        item["error"] = f"{type(error).__name__}: {error}"
        raise
    finally:
        if process is not None:
            stop_owned_group(process)
        item["elapsed_seconds"] = time.monotonic() - started
        for name, path in (("stdout", stdout), ("stderr", stderr)):
            if path.exists():
                item[name] = file_record(path)
        atomic_json(directory / "report.json", report)
    require(stdout.stat().st_size <= 1024 * 1024 and stderr.stat().st_size <= 1024 * 1024,
            "Unexpectedly large command output")
    return item


def resolve_toolchain(args, directory, deadline, report):
    """Pin only this attempt's driver/SDK selection; never change developer settings."""
    requested = shutil.which(args.compiler or "clang++")
    require(requested is not None, "Clang++ unavailable; no installation attempted")
    invocation = Path(requested).absolute()
    version = bounded_run([invocation, "--version"], "compiler-version", directory, 10, deadline, report)
    require(version["status"] == "passed", "Compiler version query failed")
    version_text = Path(version["stdout"]["path"]).read_text()
    require("clang" in version_text.lower(), "Compiler is not identified as Clang")
    # Apple's /usr/bin/clang++ is a launcher, not the compiler bytes. Resolve its
    # reported installed driver while preserving the clang++ invocation name.
    installed = re.search(r"^InstalledDir: (.+)$", version_text, re.MULTILINE)
    compiler = invocation
    if sys.platform == "darwin":
        require(installed is not None, "Apple compiler did not report InstalledDir")
        compiler = Path(installed.group(1).strip()) / "clang++"
        require(compiler.is_file(), "Reported installed C++ driver is unavailable")
    report["compiler"] = {**file_record(compiler.resolve()),
                          "invocation_path": str(compiler), "version_text": version_text,
                          "requested_invocation_path": str(invocation),
                          "requested_executable": file_record(invocation.resolve())}
    if compiler != invocation:
        actual_version = bounded_run([compiler, "--version"], "actual-compiler-version",
                                     directory, 10, deadline, report)
        require(actual_version["status"] == "passed" and
                Path(actual_version["stdout"]["path"]).read_text() == version_text,
                "Reported installed compiler does not match requested launcher")

    sdk = args.sdk_root
    if sdk is None and sys.platform == "darwin":
        selection = bounded_run(["/usr/bin/xcrun", "--show-sdk-path"], "default-sdk-selection",
                                directory, 10, deadline, report)
        require(selection["status"] == "passed", "Cannot resolve the default SDK")
        sdk = Path(Path(selection["stdout"]["path"]).read_text().strip())
    toolchain_flags = []
    sdk_record = None
    if sdk is not None:
        sdk = sdk.resolve(strict=True)
        require(sdk.is_dir(), "Selected SDK is not an installed directory")
        settings_path = sdk / "SDKSettings.json"
        settings = read_json(settings_path)
        sdk_record = {"path": str(sdk), "selection": "explicit" if args.sdk_root else "xcrun-default",
                      "version": settings.get("Version"), "canonical_name": settings.get("CanonicalName"),
                      "settings": file_record(settings_path), "key_files": [],
                      "hash_scope": "SDK metadata and listed inputs, not a hash of the entire SDK tree"}
        for relative in ("System/Library/CoreServices/SystemVersion.plist", "usr/lib/libSystem.B.tbd",
                         "usr/lib/libc++.1.tbd", "usr/include/c++/v1/vector", "usr/include/stdio.h"):
            sdk_record["key_files"].append(file_record(sdk / relative))
        toolchain_flags = ["-isysroot", str(sdk)]
        report["command_environment_overrides"] = {"SDKROOT": str(sdk)}
    linker_query = bounded_run([compiler, *toolchain_flags, "-print-prog-name=ld"],
                              "linker-selection", directory, 10, deadline, report)
    require(linker_query["status"] == "passed", "Compiler linker selection query failed")
    linker_name = Path(linker_query["stdout"]["path"]).read_text().strip()
    linker_path = shutil.which(linker_name)
    require(linker_path is not None, "Selected linker is unavailable")
    linker = Path(linker_path).absolute()
    linker_version = bounded_run([linker, "-v"], "linker-version", directory, 10, deadline, report)
    require(linker_version["status"] == "passed", "Linker version query failed")
    resource_query = bounded_run([compiler, "-print-resource-dir"], "compiler-resource-directory",
                                directory, 10, deadline, report)
    require(resource_query["status"] == "passed", "Compiler resource-directory query failed")
    report["toolchain"] = {
        "compiler": report["compiler"], "sdk": sdk_record,
        "linker": {**file_record(linker.resolve()), "invocation_path": str(linker),
                   "version_text": Path(linker_version["stdout"]["path"]).read_text() +
                                   Path(linker_version["stderr"]["path"]).read_text()},
        "resource_directory": Path(resource_query["stdout"]["path"]).read_text().strip(),
        "compiler_flags": FLAGS, "toolchain_flags": toolchain_flags,
        "command_environment_overrides": report.get("command_environment_overrides", {}),
        "inherited_build_environment": {key: os.environ.get(key) for key in
            ("SDKROOT", "DEVELOPER_DIR", "TOOLCHAINS", "MACOSX_DEPLOYMENT_TARGET", "CPATH", "CPLUS_INCLUDE_PATH",
             "C_INCLUDE_PATH", "LIBRARY_PATH", "COMPILER_PATH", "GCC_EXEC_PREFIX")},
        "selection_policy": "Explicit installed paths; no fallback after failure or persistent environment changes",
    }
    identity = json.dumps(report["toolchain"], sort_keys=True, separators=(",", ":"), allow_nan=False)
    report["toolchain_identity_sha256"] = hashlib.sha256(identity.encode()).hexdigest()
    return compiler, linker, sdk, toolchain_flags


def verify_driver_plan(command, linker, sdk, directory, deadline, report):
    """Retain the driver's non-executing plan and check the recorded linker/SDK."""
    plan = bounded_run([command[0], "-###", *command[1:]], "compiler-driver-plan",
                       directory, 10, deadline, report)
    require(plan["status"] == "passed", "Compiler driver planning failed")
    text = Path(plan["stderr"]["path"]).read_text()
    commands = [shlex.split(line) for line in text.splitlines() if line.lstrip().startswith('"')]
    linker_commands = [row for row in commands if row and Path(row[0]).resolve() == linker.resolve()]
    require(len(linker_commands) == 1, "Driver plan does not use the recorded linker exactly once")
    if sdk is not None:
        row = linker_commands[0]
        require("-syslibroot" in row and Path(row[row.index("-syslibroot") + 1]).resolve() == sdk,
                "Driver plan does not use the recorded SDK for linking")
    # The plan includes this attempt's unique output path; its file hash is kept
    # separately from the stable compiler/SDK identity recorded above.
    report["driver_plan"] = {"commands": commands, "verified": True, "stderr": plan["stderr"]}


def expand_domains(domains, nodes, minimum):
    """Validate spans without sorting away original discontinuous interval order."""
    require(isinstance(domains, list) and domains, "No domains returned")
    seen, groups = set(), []
    for domain in domains:
        require(isinstance(domain, list) and domain, "Empty domain")
        members = set()
        for span in domain:
            require(isinstance(span, list) and len(span) == 2, "Invalid interval")
            begin, end = span
            require(type(begin) is int and type(end) is int and 1 <= begin <= end <= nodes,
                    "Out-of-range one-based inclusive interval")
            indices = set(range(begin, end + 1))
            require(not members.intersection(indices), "Repeated residue within a domain")
            members.update(indices)
        require(len(members) >= minimum, "Domain smaller than frozen minimum")
        require(not seen.intersection(members), "Overlapping domains")
        seen.update(members)
        groups.append(tuple(sorted(members)))
    require(seen == set(range(1, nodes + 1)), "Partition does not cover every node exactly once")
    return tuple(sorted(groups))


def parse_partition(text, name, nodes, minimum):
    parts = text.strip().split()
    require(len(parts) == 3 and parts[0] == name and parts[1].isdigit(), "Unexpected RecCut stdout")
    require(parts[2].endswith(";"), "Unterminated domain list")
    domains = []
    for domain_text in parts[2][:-1].split(";"):
        domain = []
        for interval in domain_text.split(","):
            match = re.fullmatch(r"([1-9][0-9]*)-([1-9][0-9]*)", interval)
            require(match is not None, "Malformed interval")
            domain.append([int(match[1]), int(match[2])])
        domains.append(domain)
    require(len(domains) == int(parts[1]), "Reported domain count disagrees with spans")
    return domains, expand_domains(domains, nodes, minimum)


def make_ce(fixture):
    name, nodes = fixture["id"], fixture["nodes"]
    require(re.fullmatch(r"[a-z0-9_]{1,48}", name), "Unsafe fixture identifier")
    require(type(nodes) is int and 1 <= nodes <= 128, "Fixture node cap exceeded")
    triples = []
    if "expected_domains" in fixture:
        groups = expand_domains(fixture["expected_domains"], nodes, 22)
        owner = {node - 1: group for group, members in enumerate(groups) for node in members}
        for i in range(nodes):
            for j in range(i + 5, nodes):
                if owner[i] == owner[j]:
                    triples.append((i, j, "1.000000"))
    else:
        seen = set()
        for edge in fixture["edges"]:
            i, j = edge["i"], edge["j"]
            require(type(i) is int and type(j) is int and 0 <= i < j < nodes, "Invalid probe edge")
            require((i, j) not in seen, "Repeated probe edge")
            seen.add((i, j))
            probability = format(float(edge["input_decimal"]), ".6f")
            require(probability == edge["serialized_probability"], "CE serialization differs from frozen expectation")
            triples.append((i, j, probability))
    contacts = ",".join(f"{i} {j} {p}" for i, j, p in triples)
    return f"INF   {name} {nodes}\nCON   {contacts}\n"


def independent_ce_graph(text):
    """Separate Decimal parser for the CE text actually supplied to C++ readGraph."""
    lines = text.splitlines()
    require(len(lines) == 2 and lines[0].startswith("INF ") and lines[1].startswith("CON "),
            "Fixture must contain only INF and CON")
    nodes = int(lines[0].split()[2])
    graph = [[0] * nodes for _ in range(nodes)]
    for item in lines[1][3:].strip().split(","):
        if not item:
            continue
        fields = item.split()
        require(len(fields) == 3 and re.fullmatch(r"[01]\.[0-9]{6}", fields[2]), "Malformed six-decimal contact")
        i, j = int(fields[0]), int(fields[1])
        probability = Decimal(fields[2])
        require(0 <= i < j < nodes and 0 <= probability <= 1, "Invalid contact")
        weight = int((probability * 100 + Decimal("0.5")).to_integral_value(rounding=ROUND_FLOOR))
        graph[i][j] = graph[j][i] = weight
    for i in range(nodes):
        for j in range(i + 1, min(nodes, i + 4)):
            graph[i][j] = graph[j][i] = 100
    return graph


def parse_graph(text, nodes):
    lines = text.splitlines()
    require(lines and lines[0] == str(nodes) and len(lines) == nodes + 1, "Graph probe size mismatch")
    graph = [[int(x) for x in line.split()] for line in lines[1:]]
    require(all(len(row) == nodes for row in graph), "Graph probe row size mismatch")
    return graph


def span_validation_gate(design):
    before = json.dumps(design["valid"])
    expand_domains(design["valid"], design["nodes"], design["minimum_domain_total"])
    require(json.dumps(design["valid"]) == before, "Validator changed emitted interval order")
    checks = [{"id": "valid_discontinuous_emitted_order", "passed": True}]
    for case in design["must_reject"]:
        error = None
        try:
            expand_domains(case["domains"], design["nodes"], design["minimum_domain_total"])
        except ValueError as caught:
            error = str(caught)
        checks.append({"id": case["id"], "passed": error is not None, "rejection": error})
    return {"passed": all(c["passed"] for c in checks), "checks": checks}


def constant_quantization_gate(fingerprint_path):
    """Exercise original source under errstate; no synthetic biological string."""
    previous = sys.dont_write_bytecode
    try:
        sys.dont_write_bytecode = True
        return _constant_quantization_gate(fingerprint_path)
    finally:
        sys.dont_write_bytecode = previous


def _constant_quantization_gate(fingerprint_path):
    import numpy as np
    import scipy
    require(sha256(fingerprint_path) == FINGERPRINT_SHA256, "Changed fingerprint source")
    module_name = "_iri_owned_fingerprint_probe_" + uuid.uuid4().hex
    spec = importlib.util.spec_from_file_location(module_name, fingerprint_path)
    require(spec is not None and spec.loader is not None, "Cannot load pinned fingerprint implementation")
    module = importlib.util.module_from_spec(spec)
    sys.modules[module_name] = module
    try:
        spec.loader.exec_module(module)
        checks = []
        for name in ("constant_scale", "zero_array_quantization"):
            fingerprint = module.Fingerprint()
            if name == "zero_array_quantization":
                fingerprint.embed = {15: np.zeros((22, 640)), 21: np.zeros((22, 640))}
                fingerprint.domains = ["1-22"]
            rejected = None
            try:
                with np.errstate(invalid="raise", divide="raise", over="raise"):
                    if name == "constant_scale":
                        fingerprint.scale(np.ones(480))
                    else:
                        fingerprint.quantize([3, 80, 3, 80])
            except FloatingPointError as error:
                rejected = str(error)
            checks.append({"id": name, "passed": rejected is not None and not fingerprint.quants,
                           "exception": rejected, "completed_integer_fingerprints": len(fingerprint.quants)})
        return {"passed": all(c["passed"] for c in checks), "checks": checks,
                "numpy_version": np.__version__, "scipy_version": scipy.__version__,
                "source": file_record(fingerprint_path)}
    finally:
        sys.modules.pop(module_name, None)


def run(args):
    session = args.session.resolve()
    session.relative_to((ROOT / "results").resolve())
    build_root = session / "vendor" / "builds"
    build_root.mkdir(parents=True, exist_ok=True)
    directory = build_root / ("dct-build-" + dt.datetime.now(dt.timezone.utc).strftime("%Y%m%dT%H%M%SZ") + "-" + uuid.uuid4().hex[:12])
    directory.mkdir(exist_ok=False)
    (directory / "tmp").mkdir()
    report = {"schema_version": 1, "status": "running", "started_utc": utcnow(),
              "output_root": str(directory), "script": file_record(Path(__file__).resolve()),
              "commands": [], "fixture_results": [], "build_status": "pending",
              "graph_probe_build_status": "pending",
              "fixtures": {"continuous": None, "discontinuous": None,
                           "constant_quantization": None, "span_validation": None,
                           "ce_integerization": None},
              "scope": "Native graph and numeric-array portability gates; no biological sequence or model."}
    atomic_json(directory / "report.json", report)
    exit_code = 1
    try:
        deadline = dt.datetime.fromisoformat(read_json(session / "launch.json")["hard_stop_utc"])
        require(deadline.tzinfo is not None, "Missing timezone in session deadline")
        require(dt.datetime.now(dt.timezone.utc) < deadline, "Session hard stop reached")
        if args.previous_report is not None:
            previous = args.previous_report.resolve(strict=True)
            previous.relative_to(build_root.resolve())
            previous_value = read_json(previous)
            require(previous_value.get("status") in {"failed", "interrupted"},
                    "Previous-report reference must identify a completed failed attempt")
            report["previous_failed_attempt"] = {**file_record(previous),
                "status": previous_value["status"], "error": previous_value.get("error"),
                "preservation": "Read-only reference; previous report and logs are not modified"}
        design_path = session / "inputs" / "dct_build_fixture_design.json"
        require(sha256(design_path) == DESIGN_SHA256, "Frozen fixture design checksum changed")
        design = read_json(design_path)
        require(design["version"] == "dct-native-build-fixtures-v1" and
                design["status"] == "frozen_before_compilation_or_fixture_execution", "Unfrozen fixture design")
        require(design["compiler_flags"] == FLAGS and design["repeats"] == 3, "Build contract changed")
        require([f["id"] for f in design["fixtures"]] ==
                ["single_dense_64", "two_dense_40_40", "outer_discontinuous_32_32_32"], "Fixture set changed")
        report["design"] = file_record(design_path)
        atomic_json(directory / "frozen-design.json", design)
        source = session / "vendor/DCTdomain/src/RecCut.cpp"
        header = session / "vendor/compat/bits/stdc++.h"
        license_path = session / "vendor/DCTdomain/LICENSE"
        fingerprint_path = session / "vendor/DCTdomain/src/fingerprint.py"
        for label, path, expected in (("source", source, SOURCE_SHA256),
                                      ("compat_header", header, HEADER_SHA256),
                                      ("license", license_path, LICENSE_SHA256),
                                      ("fingerprint_source", fingerprint_path, FINGERPRINT_SHA256)):
            require(sha256(path) == expected, f"Changed pinned {label}")
            report[label] = file_record(path)
        require(design["source"]["sha256"] == SOURCE_SHA256 and
                design["source"]["compat_header_sha256"] == HEADER_SHA256 and
                design["source"]["license_sha256"] == LICENSE_SHA256 and
                design["source"]["fingerprint_sha256"] == FINGERPRINT_SHA256, "Design/source mismatch")
        report["source_sha256"] = SOURCE_SHA256
        report["compatibility_header_sha256"] = HEADER_SHA256
        report["license"]["identifier"] = "Academic Free License 3.0"
        report["acquisition_manifest"] = file_record(session / "audit/dct_source_acquisition.json")
        shutil.copyfile(license_path, directory / "LICENSE.upstream")
        for name, probe in (("span_validation", lambda: span_validation_gate(design["span_validation_probe"])),
                            ("constant_quantization", lambda: constant_quantization_gate(fingerprint_path))):
            try:
                result = probe()
                report[name + "_evidence"] = result
                report["fixtures"][name] = result["passed"]
            except Exception as error:
                report[name + "_evidence"] = {"passed": False, "error": f"{type(error).__name__}: {error}"}
                report["fixtures"][name] = False
            atomic_json(directory / "report.json", report)
        compiler, linker, sdk, toolchain_flags = resolve_toolchain(args, directory, deadline, report)
        build_identity = {"script_sha256": report["script"]["sha256"],
                          "design_sha256": report["design"]["sha256"],
                          "source_sha256": SOURCE_SHA256, "compatibility_header_sha256": HEADER_SHA256,
                          "toolchain_identity_sha256": report["toolchain_identity_sha256"]}
        report["build_identity"] = build_identity
        report["build_identity_sha256"] = hashlib.sha256(json.dumps(
            build_identity, sort_keys=True, separators=(",", ":")).encode()).hexdigest()
        binary = directory / "RecCut"
        compile_command = [compiler, *FLAGS, *toolchain_flags, "-I", header.parent.parent, source, "-o", binary]
        verify_driver_plan(compile_command, linker, sdk, directory, deadline, report)
        compile_result = bounded_run(compile_command,
                                     "compile-reccut", directory, 90, deadline, report)
        report["build_status"] = compile_result["status"]
        require(compile_result["status"] == "passed" and binary.is_file(), "Native RecCut compilation failed")
        report["binary"] = file_record(binary)
        report["binary_sha256"] = report["binary"]["sha256"]
        # This separate harness calls the unchanged parser; the main binary above
        # has no preprocessor main renaming or other source transformation.
        harness = directory / "graph-probe.cpp"
        harness.write_text("#define main iri_original_reccut_main\n#include " + json.dumps(str(source)) +
            "\n#undef main\nint main(int argc, char** argv) {\n"
            "  if (argc != 2) return 2;\n  auto graph = readGraph(argv[1]);\n"
            "  std::cout << graph.size() << '\\n';\n"
            "  for (const auto& row : graph) {\n"
            "    for (int value : row) std::cout << value << ' ';\n"
            "    std::cout << '\\n';\n  }\n  return 0;\n}\n")
        report["graph_probe_source"] = file_record(harness)
        probe_binary = directory / "RecCut_graph_probe"
        probe_compile = bounded_run([compiler, *FLAGS, *toolchain_flags, "-I", header.parent.parent, harness, "-o", probe_binary],
                                    "compile-graph-probe", directory, 90, deadline, report)
        report["graph_probe_build_status"] = probe_compile["status"]
        require(probe_compile["status"] == "passed" and probe_binary.is_file(), "Graph probe compilation failed")
        report["graph_probe_binary"] = file_record(probe_binary)
        for fixture in [*design["fixtures"], design["rounding_probe"]]:
            record = {"id": fixture["id"], "status": "running", "partition_runs": []}
            report["fixture_results"].append(record)
            try:
                ce_text = make_ce(fixture)
                ce_path = directory / (fixture["id"] + ".ce")
                with ce_path.open("x") as stream:
                    stream.write(ce_text)
                record["ce"] = file_record(ce_path)
                expected_graph = independent_ce_graph(ce_text)
                atomic_json(directory / (fixture["id"] + ".expected-graph.json"), expected_graph)
                if "edges" in fixture:
                    for edge in fixture["edges"]:
                        require(expected_graph[edge["i"]][edge["j"]] == edge["expected_integer_after_adjacency"],
                                "Frozen integerization example disagrees with independent parser")
                graph_run = bounded_run([probe_binary, ce_path], fixture["id"] + "-graph", directory, 10, deadline, report)
                require(graph_run["status"] == "passed", "Native graph parser probe failed")
                actual_graph = parse_graph(Path(graph_run["stdout"]["path"]).read_text(), fixture["nodes"])
                require(actual_graph == expected_graph, "Native readGraph differs from independent CE integerization")
                record["graph_cells_checked"] = fixture["nodes"] ** 2
                record["graph_equal"] = True
                if "expected_domains" in fixture:
                    expected = expand_domains(fixture["expected_domains"], fixture["nodes"], 22)
                    record["expected_domains"] = fixture["expected_domains"]
                    for repeat in range(design["repeats"]):
                        item = {"repeat": repeat, "status": "failed"}
                        record["partition_runs"].append(item)
                        try:
                            result = bounded_run([binary, "--input", ce_path, "--name", fixture["id"]],
                                f"{fixture['id']}-partition-{repeat}", directory, 10, deadline, report)
                            require(result["status"] == "passed", "RecCut fixture process failed")
                            text = Path(result["stdout"]["path"]).read_text()
                            item["stdout_sha256"] = result["stdout"]["sha256"]
                            domains, actual = parse_partition(text, fixture["id"], fixture["nodes"], 22)
                            item["actual_domains_emitted_order"] = domains
                            item["valid_complete_disjoint_partition"] = True
                            require(actual == expected, "Partition differs from frozen intended graph groups")
                            item["status"] = "passed"
                        except Exception as error:
                            item["error"] = f"{type(error).__name__}: {error}"
                    hashes = [item.get("stdout_sha256") for item in record["partition_runs"]]
                    record["deterministic_stdout"] = all(hashes) and len(set(hashes)) == 1
                    require(record["deterministic_stdout"] and all(item["status"] == "passed" for item in record["partition_runs"]),
                            "Partition fixture failed expectation, validity or determinism gate")
                record["status"] = "passed"
            except Exception as error:
                record["status"] = "failed"
                record["error"] = f"{type(error).__name__}: {error}"
            atomic_json(directory / "report.json", report)
        require(sha256(source) == SOURCE_SHA256 and sha256(header) == HEADER_SHA256,
                "Source/header changed during build gate")
        by_id = {fixture["id"]: fixture for fixture in report["fixture_results"]}
        report["fixtures"]["continuous"] = all(by_id[name]["status"] == "passed" for name in
                                                  ("single_dense_64", "two_dense_40_40"))
        report["fixtures"]["discontinuous"] = by_id["outer_discontinuous_32_32_32"]["status"] == "passed"
        report["fixtures"]["ce_integerization"] = len(by_id) == 4 and all(f.get("graph_equal") is True for f in by_id.values())
        report["status"] = "passed" if all(v is True for v in report["fixtures"].values()) else "failed"
        exit_code = 0 if report["status"] == "passed" else 1
    except BaseException as error:
        report["status"] = "interrupted" if isinstance(error, (SystemExit, KeyboardInterrupt)) else "failed"
        report["error"] = f"{type(error).__name__}: {error}"
        exit_code = int(error.code) if isinstance(error, SystemExit) and isinstance(error.code, int) else 1
    finally:
        report["finished_utc"] = utcnow()
        report["limitations"] = "Portability/parser/fixture evidence only; no model/DCT fingerprint admission or proof of general RecCut correctness."
        atomic_json(directory / "report.json", report)
    print(json.dumps({"status": report["status"], "report": str(directory / "report.json")}))
    return exit_code


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--session", required=True, type=Path)
    parser.add_argument("--compiler", help="Explicit Clang++ executable; defaults to PATH clang++")
    parser.add_argument("--sdk-root", type=Path,
                        help="Installed macOS SDK selected only for this attempt; no automatic fallback")
    parser.add_argument("--previous-report", type=Path,
                        help="Prior failed build report retained unchanged and linked by SHA256")
    args = parser.parse_args()
    def interrupted(signum, frame):
        raise SystemExit(128 + signum)
    signal.signal(signal.SIGTERM, interrupted)
    signal.signal(signal.SIGINT, interrupted)
    return run(args)


if __name__ == "__main__":
    sys.exit(main())
