-- xmake build for carquet — a pure C Apache Parquet library.
--
-- This mirrors CMakeLists.txt. Quick start:
--
--   xmake f -m release            # configure (release is the default mode)
--   xmake                         # build the library + CLI
--   xmake f --dev=y && xmake      # build tests, examples, benchmarks, interop too
--   xmake run test_core           # run a single test
--   xmake test                    # run every test target (ctest equivalent)
--
-- Enable optional pieces individually, e.g. `xmake f --tests=y --fuzz=y`.

set_project("carquet")
set_version("0.7.2")
-- gnu11, not strict c11: matches CMake's `set(CMAKE_C_STANDARD 11)` (GNU
-- extensions stay on by default). Under strict -std=c11 the POSIX feature-test
-- macros are off, so functions like strdup/posix_memalign aren't declared and
-- get implicitly typed to return int — truncating 64-bit pointers to 32-bit and
-- corrupting the heap (free of a non-allocated pointer).
set_languages("gnu11")

-- mode.release gives -O3 + NDEBUG, mode.debug gives -O0 -g. We layer the
-- carquet-specific extras (LTO, CARQUET_DEBUG) on top below.
add_rules("mode.debug", "mode.release", "mode.releasedbg")

-- CMake puts `-flto=auto` in CMAKE_C_FLAGS_RELEASE, so it applies to EVERY
-- target (library, tests, examples, CLI) — the whole program is optimized and
-- LTO'd together. Match that here at the root scope rather than per-target:
-- compiling only the library with LTO leaves a codegen mismatch that shifts
-- stack/heap layout and can surface latent UB. Library files that must stay out
-- of LTO (detect.c, the ISA-specific SIMD TUs) re-add -fno-lto per file below;
-- that flag comes last on those compiles and wins. GCC/Clang only.
if is_mode("release") and not is_plat("windows") then
    add_cflags("-flto=auto")
    add_ldflags("-flto=auto")
end

-- Compiler-independent notion of "gcc or clang" vs MSVC. CMake keys off the
-- compiler id; keying off the platform is close enough for the flag choices
-- carquet makes (Windows == MSVC-style, everything else == gcc-like).
local GCC_LIKE = not is_plat("windows")
local IS_X86   = is_arch("x86_64", "x64", "i386", "x86")
local IS_ARM   = is_arch("arm64", "arm64-v8a", "aarch64", "armv7", "armv7s", "arm", "armeabi", "armeabi-v7a")

--------------------------------------------------------------------------------
-- Build options (parallel to the CMake `option(...)` block)
--------------------------------------------------------------------------------

-- Umbrella switch: turns on tests + examples + benchmarks + interop at once.
option("dev")
    set_default(false)
    set_showmenu(true)
    set_description("Build all development targets (tests, examples, benchmarks, interop)")
option_end()

for _, o in ipairs({
    {"tests",      false, "Build tests"},
    {"examples",   false, "Build example programs"},
    {"benchmarks", false, "Build benchmark and profiling programs"},
    {"interop",    false, "Build interoperability test programs"},
    {"cli",        true,  "Build carquet CLI tool"},
    {"shared",     false, "Build shared library"},
    {"fuzz",       false, "Build fuzz targets"},
    {"arrow_cpp_benchmark", false, "Build the optional Arrow C++ comparison benchmark"},
    {"native_arch", false, "Build with -march=native (produces a HOST-ONLY binary)"},
    {"sse",  true,  "Enable SSE optimizations"},
    {"avx",  true,  "Enable AVX optimizations"},
    {"avx2", true,  "Enable AVX2 optimizations"},
    {"neon", true,  "Enable NEON optimizations"},
    {"sve",  false, "Enable SVE optimizations"},
}) do
    option(o[1])
        set_default(o[2])
        set_showmenu(true)
        set_description(o[3])
    option_end()
end

-- AVX-512 auto-disables when the compiler lacks full intrinsic support
-- (GCC < 9 misses _mm512_set_epi8), matching the CMake CheckCSourceCompiles.
option("avx512")
    set_default(true)
    set_showmenu(true)
    set_description("Enable AVX-512 optimizations")
    on_check(function (option)
        if not is_arch("x86_64", "x64", "i386", "x86") then return end
        if not option:enabled() then return end
        local ok, has = pcall(function ()
            import("lib.detect.check_cxsnippets")
            return check_cxsnippets({[[
                void carquet_avx512_probe(void) {
                    __m512i v = _mm512_set_epi8(
                        0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,
                        16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,
                        32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,
                        48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63);
                    (void)v;
                }
            ]]}, {sourcekind = "cc", includes = "immintrin.h",
                  cflags = {"-mavx512f", "-mavx512bw", "-mavx512vl"}})
        end)
        -- Only disable when the probe ran and failed; keep enabled if it errored.
        if ok and not has then
            option:enable(false)
            cprint("${yellow}Compiler lacks full AVX-512 intrinsic support — disabling AVX-512${clear}")
        end
    end)
option_end()

-- OpenMP for parallel column reading. Auto-detected (like CMake's
-- find_package(OpenMP)): the on_check snippet must compile AND link, else the
-- option disables itself. Never fetches/builds anything. User-overridable
-- (`xmake f --openmp=n`). Applied declaratively on the carquet target below —
-- doing it there (not in on_load/on_config) is what makes the flag reach the
-- shared library's own link, not just its interface.
option("brotli")
    set_default(true)
    set_showmenu(true)
    set_description("Build the BROTLI codec (libbrotli fetched from the xmake repo)")
option_end()

option("openmp")
    set_default(true)
    set_showmenu(true)
    set_description("Enable OpenMP parallel column reading (auto-disabled if unavailable)")
    on_check(function (option)
        if not option:enabled() then return end
        local ok = try { function ()
            import("lib.detect.check_cxsnippets")
            local flag = is_plat("windows") and "/openmp" or "-fopenmp"
            return check_cxsnippets({[[
                #include <omp.h>
                void carquet_omp_probe(void) {
                    #pragma omp parallel
                    { (void)omp_get_thread_num(); }
                }
            ]]}, {sourcekind = "cc", cflags = flag, ldflags = flag})
        end }
        if not ok then option:enable(false) end
    end)
option_end()

-- Flag used to compile and link OpenMP when the option is enabled.
local OPENMP_FLAG = is_plat("windows") and "/openmp" or "-fopenmp"

-- tests/examples/benchmarks/interop inherit from `dev` unless set explicitly.
local function want(feature)
    return has_config(feature) or has_config("dev")
end

--------------------------------------------------------------------------------
-- Compression dependencies (system first, auto-download fallback — like CMake's
-- find_package-then-FetchContent chain). Versions match the CMake GIT_TAGs.
--------------------------------------------------------------------------------

-- Built as self-contained static archives from xmake's package repo, at the same
-- versions as the CMake GIT_TAGs. This is portable across Linux/macOS/Windows and
-- reproducible: the binaries carry no external compression-library dependency, so
-- there are no runtime rpath/dylib issues (on macOS this also sidesteps linking a
-- stray shared libz that may sit on the default search path). First configure
-- downloads and builds these once, then caches them.
-- To link the host's system libraries instead (like CMake's find_package path),
-- drop {system = false} — xmake then prefers an installed zstd/zlib/lz4.
add_requires("zstd 1.5.6", {system = false, configs = {shared = false}})
add_requires("zlib 1.3.1", {system = false, configs = {shared = false}})
add_requires("lz4 1.10.0", {system = false, configs = {shared = false}})
if has_config("brotli") then
    add_requires("brotli 1.1.0", {system = false, configs = {shared = false}})
end

-- Parallel column reading (OpenMP) is handled by the auto-detected `openmp`
-- option above and applied on the carquet target below, mirroring CMake's
-- detect-only find_package(OpenMP): nothing is fetched or built, and a missing
-- OpenMP is silently skipped.

--------------------------------------------------------------------------------
-- The carquet library
--------------------------------------------------------------------------------

target("carquet")
    set_kind(has_config("shared") and "shared" or "static")

    -- Public API headers propagate to every dependent; src/ stays private.
    add_includedirs("include", {public = true})
    add_includedirs("src")

    add_packages("zstd", "zlib", "lz4", {public = true})
    if has_config("brotli") then
        add_packages("brotli", {public = true})
        add_defines("CARQUET_HAVE_BROTLI=1")
    end

    -- OpenMP (auto-detected via the `openmp` option). Declared here so the flag
    -- reaches carquet's own link (the shared dylib) as well as its interface;
    -- {public = true} propagates it to executables that link the static archive.
    if has_config("openmp") then
        add_cflags(OPENMP_FLAG)
        add_shflags(OPENMP_FLAG)                   -- carquet's own link when built shared (dylib)
        add_ldflags(OPENMP_FLAG, {public = true})  -- dependents (executables linking the archive)
    end

    -- When building with LTO under Homebrew LLVM, the objects/archive hold LLVM
    -- bitcode that Apple's ar/ranlib can't index — the archive is malformed and
    -- binaries crash at runtime. Point the archiver at the matching
    -- llvm-ar/llvm-ranlib, exactly like the CMakeLists does. macOS-only; no
    -- effect on Linux/Windows. (Needs on_config: target:tool() isn't available
    -- during on_load / description scope.)
    on_config(function (target)
        if not is_plat("macosx") then return end
        local cc = target:tool("cc")
        if cc and cc:find("llvm", 1, true) then
            local bindir = path.directory(cc)
            local ar, ranlib = path.join(bindir, "llvm-ar"), path.join(bindir, "llvm-ranlib")
            if os.isfile(ar) then target:set("toolset.ar", ar) end
            if os.isfile(ranlib) then target:set("toolset.ranlib", ranlib) end
        end
    end)

    -- Warnings (CMake: -Wall -Wextra -Wpedantic -Wno-unused-parameter).
    if GCC_LIKE then
        add_cflags("-Wall", "-Wextra", "-Wpedantic", "-Wno-unused-parameter")
    else
        add_cflags("/W3", "/wd4996", "/wd4244", "/wd4267")
    end

    -- Build-type extras layered onto the mode.* rules. (-flto=auto is applied
    -- globally at the root scope, matching CMAKE_C_FLAGS_RELEASE.)
    if is_mode("debug") then
        add_defines("CARQUET_DEBUG")
    end

    -- -march=native: whole-library, non-portable, host-only (see CMake note).
    if has_config("native_arch") and GCC_LIKE then
        add_cflags("-march=native")
    end

    -- Portable core sources. snappy.c and detect.c are pulled in separately
    -- below because they need per-file flags; the '|foo.c' suffix excludes them.
    add_files(
        "src/core/*.c",
        "src/thrift/*.c",
        "src/encoding/*.c",
        "src/compression/*.c|snappy.c",
        "src/simd/*.c|detect.c",
        "src/reader/*.c",
        "src/writer/*.c",
        "src/metadata/*.c",
        "src/util/*.c")

    -- Snappy uses SSSE3 (pshufb) for pattern extension on x86 (baseline since
    -- Core 2 / x86-64-v2).
    if IS_X86 and has_config("sse") and GCC_LIKE then
        add_files("src/compression/snappy.c", {cflags = {"-mssse3"}})
    else
        add_files("src/compression/snappy.c")
    end

    -- Runtime CPU feature detection must stay correct under mixed LLVM/macOS
    -- toolchains: exclude from LTO so feature flags aren't miscompiled.
    if GCC_LIKE then
        add_files("src/simd/detect.c", {cflags = {"-fno-lto"}})
    else
        add_files("src/simd/detect.c")
    end

    -- x86 SIMD. ISA-specific TUs build with -fno-lto so the code generator can't
    -- inline AVX/AVX-512 into baseline callers ahead of the dispatch check.
    if IS_X86 then
        add_defines("CARQUET_ARCH_X86")
        if has_config("sse") then
            add_defines("CARQUET_ENABLE_SSE")
            add_files("src/simd/x86/sse_ops.c",
                {cflags = GCC_LIKE and {"-msse4.2", "-fno-lto"} or {}})
        end
        if has_config("avx") then
            add_defines("CARQUET_ENABLE_AVX")
            add_files("src/simd/x86/avx_ops.c",
                {cflags = GCC_LIKE and {"-mavx", "-msse4.2", "-fno-lto"} or {"/arch:AVX"}})
        end
        if has_config("avx2") then
            add_defines("CARQUET_ENABLE_AVX2")
            add_files("src/simd/x86/avx2_ops.c",
                {cflags = GCC_LIKE and {"-mavx2", "-mbmi2", "-fno-lto"} or {"/arch:AVX2"}})
        end
        if has_config("avx512") then
            add_defines("CARQUET_ENABLE_AVX512")
            add_files("src/simd/x86/avx512_ops.c",
                {cflags = GCC_LIKE and {"-mavx512f", "-mavx512bw", "-mavx512vl", "-fno-lto"}
                                    or {"/arch:AVX512"}})
        end
    end

    -- ARM SIMD.
    if IS_ARM then
        add_defines("CARQUET_ARCH_ARM")
        if has_config("neon") then
            add_defines("CARQUET_ENABLE_NEON")
            add_files("src/simd/arm/neon_ops.c")
        end
        if has_config("sve") then
            add_defines("CARQUET_ENABLE_SVE")
            add_files("src/simd/arm/sve_ops.c",
                {cflags = GCC_LIKE and {"-march=armv8-a+sve"} or {}})
        end
    end

    -- Shared-library symbol export.
    if has_config("shared") then
        add_defines("CARQUET_BUILD_SHARED", {public = true})
        add_defines("CARQUET_BUILDING_DLL")
        if is_plat("windows") then
            add_rules("utils.symbols.export_all")
        end
    end

    -- libm for bloom_filter.c (log()); pthreads for the worker pool.
    if is_plat("linux", "macosx", "bsd") then
        add_syslinks("m", "pthread", {public = true})
    end
target_end()

--------------------------------------------------------------------------------
-- Helper for the many small executables that just link libcarquet.
--------------------------------------------------------------------------------

-- kind: "test" registers a ctest-style test; "app" is a plain executable.
local function carquet_exe(name, files, opts)
    opts = opts or {}
    target(name)
        set_kind("binary")
        set_default(false)  -- only built when its group option is on
        add_deps("carquet")
        -- carquet is a static archive; name the (static) compression packages
        -- directly on each executable so their link flags are pulled in.
        add_packages("zstd", "zlib", "lz4")
        if has_config("brotli") then add_packages("brotli") end
        add_files(files)
        if opts.basename then set_basename(opts.basename) end
        if opts.with_src then add_includedirs("src") end
        if GCC_LIKE then add_cflags("-Wno-unused-result") end
        -- Tests keep assertions live even in release — and some wrap
        -- side-effecting calls in assert() (e.g. `assert(arena_init(...) == OK)`),
        -- so if asserts are compiled out those calls VANISH and the test corrupts
        -- memory. mode.release adds -DNDEBUG; xmake emits every -U before every
        -- -D, so a trailing -UNDEBUG can't win. Strip the define outright instead.
        if opts.keep_asserts then
            -- Some tests wrap essential, side-effecting calls in assert()
            -- (e.g. `assert(carquet_schema_add_column(...) == OK)`), so if
            -- asserts are compiled out those operations never run and the test
            -- corrupts memory / reads a truncated file. xmake's mode.release
            -- rule appends `-DNDEBUG` as a raw cxflag in its on_config; a raw
            -- `-UNDEBUG` cxflag added afterward (raw cxflags keep add-order,
            -- unlike -D/-U defines which xmake regroups) re-enables asserts.
            on_config(function (target)
                local undef = target:has_tool("cc", "cl") and "/UNDEBUG" or "-UNDEBUG"
                target:add("cxflags", undef, {force = true})
            end)
        end
        -- MSVC gates <stdatomic.h> behind an opt-in flag; tests that use C11
        -- atomics directly (the library itself avoids them on MSVC) need it.
        if opts.cl_cflags then
            on_config(function (target)
                if target:has_tool("cc", "cl") then
                    target:add("cflags", opts.cl_cflags, {force = true})
                end
            end)
        end
        if opts.is_test then add_tests("default") end
    target_end()
end

--------------------------------------------------------------------------------
-- Tests
--------------------------------------------------------------------------------

if want("tests") then
    local tests = {
        "test_core", "test_thrift", "test_encodings", "test_reader",
        "test_encodings_extended", "test_compression", "test_utils",
        "test_production", "test_maturity", "test_large_schema",
        "test_edge_encodings", "test_edge_compression", "test_edge_malformed",
        "test_edge_boundaries", "test_edge_io", "test_nested",
        "test_bloom_page_index", "test_mmap", "test_advanced_api",
        "test_encoding_roundtrip", "test_writer_extensions", "test_bitunpack_wide",
        "test_page_filter", "test_append", "test_row_range", "test_dict_decode",
        "test_decimal_stats", "test_arrow_dictionary", "test_snappy_diag", "test_parallel_write",
        "test_float16", "test_geo_wkb", "test_worker_pool", "test_custom_codec",
        "test_real_world", "test_arrow_c_data", "test_batch_nested",
        "test_nested_write", "test_arrow_nested", "test_p2_conformance",
        "test_cli_format", "test_fixes_0_7_2",
    }
    -- Tests that link extra non-library sources (the CLI helpers live in the
    -- carquet_cli target, not in libcarquet).
    local test_extra_files = {
        test_cli_format = {"src/cli/commands.c"},
    }
    for _, t in ipairs(tests) do
        local exe_opts = {with_src = true, keep_asserts = true, is_test = true}
        if t == "test_worker_pool" then exe_opts.cl_cflags = "/experimental:c11atomics" end
        local files = {"tests/" .. t .. ".c"}
        for _, extra in ipairs(test_extra_files[t] or {}) do
            table.insert(files, extra)
        end
        carquet_exe(t, files, exe_opts)
    end
end

--------------------------------------------------------------------------------
-- Examples
--------------------------------------------------------------------------------

if want("examples") then
    local examples = {
        {"example_basic_write_read", "examples/basic_write_read.c"},
        {"example_data_types",       "examples/data_types.c"},
        {"example_compression",      "examples/compression_codecs.c"},
        {"example_nullable",         "examples/nullable_columns.c"},
        {"example_advanced",         "examples/advanced_features.c"},
        {"example_append",           "examples/append_rows.c"},
        {"example_page_filter",      "examples/page_filter.c"},
        {"example_nested",           "examples/nested_data.c"},
        {"example_arrow_c_data",     "examples/arrow_c_data.c"},
    }
    for _, e in ipairs(examples) do
        carquet_exe(e[1], e[2])
    end
end

--------------------------------------------------------------------------------
-- Benchmarks + profiling
--------------------------------------------------------------------------------

if want("benchmarks") then
    carquet_exe("benchmark_carquet",   "benchmark/benchmark_carquet.c")
    carquet_exe("generate_test_files", "benchmark/generate_test_files.c")

    -- Profiling targets use dirent.h / Unix profilers — skip on Windows.
    if not is_plat("windows") then
        carquet_exe("profile_core",  "benchmark/profile_core.c")
        carquet_exe("profile_read",  "profiling/profile_read.c")
        carquet_exe("profile_write", "profiling/profile_write.c")
        carquet_exe("profile_micro", "profiling/profile_micro.c", {with_src = true})
    end

    -- Optional Arrow/Parquet C++ comparison benchmark.
    if has_config("arrow_cpp_benchmark") then
        add_requires("arrow", {optional = true})
        if has_package("arrow") then
            target("benchmark_arrow_cpp")
                set_kind("binary")
                set_default(false)
                set_languages("c++20")
                add_files("benchmark/benchmark_arrow_cpp.cpp")
                add_packages("arrow")
            target_end()
        end
    end
end

--------------------------------------------------------------------------------
-- Interoperability tests (Unix only, like CMake)
--------------------------------------------------------------------------------

if want("interop") and not is_plat("windows") then
    carquet_exe("test_interop",      "interop/test_interop.c", {with_src = true})
    carquet_exe("roundtrip_writer",  "interop/roundtrip_writer.c")
end

--------------------------------------------------------------------------------
-- CLI tool (installed as `carquet`)
--------------------------------------------------------------------------------

if has_config("cli") then
    target("carquet_cli")
        set_kind("binary")
        set_basename("carquet")
        add_deps("carquet")
        add_packages("zstd", "zlib", "lz4")
        if has_config("brotli") then add_packages("brotli") end
        add_includedirs("src")
        add_files(
            "src/cli/main.c",
            "src/cli/commands.c",
            "src/cli/codegen.c",
            "src/cli/codegen_read.c",
            "src/cli/codegen_write.c")
        if GCC_LIKE then add_cflags("-Wno-unused-result") end
    target_end()
end

--------------------------------------------------------------------------------
-- Fuzz targets (clang + libFuzzer/AFL, sanitizers on). Build with e.g.
--   xmake f --fuzz=y --toolchain=clang -m debug && xmake
--------------------------------------------------------------------------------

if has_config("fuzz") then
    option("fuzz_sanitizers")
        set_default(true)
        set_showmenu(true)
        set_description("Enable sanitizers for fuzzing")
    option_end()

    option("fuzz_engine")
        set_default("libFuzzer")
        set_showmenu(true)
        set_values("libFuzzer", "AFL")
        set_description("Fuzzing engine")
    option_end()

    local san = has_config("fuzz_sanitizers")
        and {"-fsanitize=address", "-fsanitize=undefined", "-fno-omit-frame-pointer"}
        or {}
    local engine = get_config("fuzz_engine") or "libFuzzer"

    local fuzzers = {
        "fuzz_reader", "fuzz_compression", "fuzz_encodings", "fuzz_thrift",
        "fuzz_roundtrip", "fuzz_writer", "fuzz_page_filter", "fuzz_append",
        "fuzz_arrow_schema", "fuzz_arrow_c_data", "fuzz_nested_write",
        "fuzz_arrow_read",
    }
    for _, f in ipairs(fuzzers) do
        target(f)
            set_kind("binary")
            set_default(false)
            add_deps("carquet")
            add_packages("zstd", "zlib", "lz4")
            if has_config("brotli") then add_packages("brotli") end
        if has_config("brotli") then add_packages("brotli") end
            add_includedirs("src", "include")
            add_files("fuzz/" .. f .. ".c")
            if engine == "libFuzzer" then
                add_cflags("-fsanitize=fuzzer", san)
                add_ldflags("-fsanitize=fuzzer", san)
            else -- AFL
                add_defines("AFL_MAIN")
                if san and #san > 0 then
                    add_cflags(san)
                    add_ldflags(san)
                end
            end
        target_end()
    end
end

--------------------------------------------------------------------------------
-- Documentation (Doxygen). Run with `xmake docs`. Reuses the checked-in
-- Doxyfile (output goes to build/docs/html). Doxygen is not required for a
-- normal build: the task prints a clear hint if it (or graphviz) is missing.
--------------------------------------------------------------------------------

task("docs")
    on_run(function ()
        import("lib.detect.find_tool")
        local doxygen = find_tool("doxygen")
        if not doxygen then
            cprint("${red}doxygen not found${clear} — install it (e.g. `brew install doxygen graphviz`) to build the docs")
            return
        end
        if not find_tool("dot") then
            cprint("${yellow}graphviz 'dot' not found${clear} — docs will build without diagrams")
        end
        os.execv(doxygen.program, {"Doxyfile"}, {curdir = os.projectdir()})
        cprint("${green}docs generated${clear} -> build/docs/html/index.html")
    end)
    set_menu({usage = "xmake docs", description = "Generate API documentation with Doxygen (needs doxygen; graphviz for diagrams)"})
task_end()
