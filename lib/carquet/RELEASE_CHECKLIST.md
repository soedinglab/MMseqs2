# Checklist

Concise tick-list for maintaining carquet. Rationale and exact expectations for each item are in **[CONTRIBUTING.md](CONTRIBUTING.md)** — this file is just the actions.

## After every code change

```bash
# 1. Full test suite (Debug)
cmake -B build -DCMAKE_BUILD_TYPE=Debug -DCARQUET_BUILD_TESTS=ON && cmake --build build
( cd build && ctest --output-on-failure )

# 3. ASan + UBSan clean
cmake -B build-san -DCMAKE_BUILD_TYPE=Debug -DCARQUET_BUILD_TESTS=ON \
  -DCMAKE_C_FLAGS="-fsanitize=address,undefined -fno-sanitize-recover=all -fno-omit-frame-pointer" && cmake --build build-san
( cd build-san && ctest --output-on-failure )

# 4. Fuzz every reachable target ≥60s, sanitizers on
python3 fuzz/run_fuzzer.py all --time 60

# 7. Interop (if reader/writer/wire touched)
python3 interop/run_interop.py -v
```

- [ ] **Tests pass** — full `ctest`, all targets (not just the one you edited).
- [ ] **Unit tests added** for the new API: happy path + an error/rejection path + boundaries. Wire new `tests/test_*.c` into **both** `CMakeLists.txt` and `xmake.lua`.
- [ ] **Byte-level assertion** if you wrote a wire-format encoder (Thrift / page bytes) — assert exact bytes, not just self round-trip.
- [ ] **ASan + UBSan clean** — zero errors, zero leaks.
- [ ] **Fuzz ≥60s per reachable target**, zero crashes. Extend/add the fuzzer if the change opens a new path (`fuzz/run_fuzzer.py list` for targets); a new target goes in `fuzz/CMakeLists.txt`, the `fuzzers` list in `xmake.lua`, and `run_fuzzer.py`.
- [ ] **Docs updated** — `docs/`, public-header doc comments, CLAUDE.md, examples, and the component README that owns the area (`fuzz/README.md`, `docs/README.md`, `profiling/README.md`). Update the **main `README.md`** only for user-visible changes, keep it concise, and **link out** rather than pasting big chunks.
- [ ] **CHANGELOG bullet** in the top unreleased section (right heading; note ABI/default-bytes impact).
- [ ] xmake.lua updated
- [ ] **Interop passes** + new case added in `roundtrip_writer.c` if a new wire feature landed.

One note: codebase contains fuzz/external/parquet-testing with good and malformed apache arrow reference tests files. Might be useful for testing.

## Before tagging a release

Do the above on the final `main`, then:

```bash
# Long fuzz soak, all targets, parallel
python3 fuzz/run_fuzzer.py all --time 20 --jobs 4

git diff                                            # confirm every hit

# Interop snapshot for the record
python3 interop/run_interop.py -v

# Tag
git tag vX.Y.Z && git push origin main --tags
```

- [ ] **Re-run the "after every change" gate** on the exact commit to be tagged.
- [ ] **Fuzz soak** (hours, `--jobs N`) clean; triage/fix any crash artifact, commit useful new corpus.
- [ ] **`bump_version.py`** run; confirm `CMakeLists.txt` `project(VERSION)` and `xmake.lua` `set_version` agree. It also rewrites `version:` in `CITATION.cff` — hand-update that file's `date-released:` to the release date (the script does not touch dates).
- [ ] **CHANGELOG finalized** — unreleased bullets rolled into a dated `## vX.Y.Z` section with a one-paragraph theme + ABI/default-bytes statement.
- [ ] **Interop snapshot committed** — `interop/interop_<version>_<platform>_<date>.json`.
