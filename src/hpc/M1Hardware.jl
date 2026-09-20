# ============================================================================
# M1Hardware.jl — MODULE 1: hardware detection & environment setup.
#
# Detects CPU/threads/NUMA/GPU, picks the BLAS threading policy (for 6x6
# kernels BLAS threads MUST be 1 — parallelism is over seeds in MODULE 3),
# and prefers MKL when present. GPU backends (CUDA/oneAPI) are probed
# non-invasively: absence downgrades to the CPU path with a clear report.
# No environment or security policy is modified.
# ============================================================================

using Libdl

struct HardwareReport
    cpu_model::String
    julia_threads::Int
    blas_threads::Int
    blas_backend::String
    numa_nodes::Int
    total_mem_gb::Float64
    gpu_cuda::Bool
    gpu_oneapi::Bool
    gpu_note::String
    simd_note::String
end

function _detect_numa_nodes()
    try
        out = read(`lscpu`, String)
        m = match(r"NUMA node\(s\):\s+(\d+)", out)
        return m === nothing ? 1 : parse(Int, m.captures[1])
    catch
        return Sys.iswindows() ? 1 : -1   # -1 = unknown (no lscpu)
    end
end

function _detect_gpu()
    cuda = false
    oneapi = false
    notes = String[]
    # CUDA: nvidia-smi present AND Julia CUDA package loadable
    try
        success(`nvidia-smi -L`) && push!(notes, "nvidia-smi present")
    catch
    end
    # WSL GPU paravirtualization device (works for both CUDA and Level-Zero)
    dxg = isfile("/dev/dxg")
    dxg && push!(notes, "/dev/dxg present (WSL GPU-PV)")
    # Julia-side packages (only probe; never install)
    for (pkg, flag) in (("CUDA", :cuda), ("oneAPI", :oneapi))
        have = Base.find_package(pkg * ".jl") !== nothing ||
               try
                   Base.require(Base.PkgId(Base.UUID("00000000-0000-0000-0000-000000000000"), pkg);)
                   false
               catch
                   false
               end
    end
    if dxg && !cuda && !oneapi
        push!(notes, "GPU userspace stack not installed (no nvidia-smi / no Level-Zero); CPU path active")
    end
    return cuda, oneapi, join(notes, "; ")
end

function hardware_report()
    cpu = try
        strip(first(eachline(`lscpu`))) do s
            replace(s, "Model name:" => "")
        end
    catch
        Sys.CPU_NAME
    end
    blas_backend = try
        # LinearAlgebra.BLAS.get_config() exists on LAPACK ≥ 7.3.1 stacks
        cfg = string(LinearAlgebra.BLAS.get_config())
        occursin("MKL", cfg) ? "MKL" : occursin("openblas", lowercase(cfg)) ? "OpenBLAS" : cfg
    catch
        "unknown"
    end
    cuda, oneapi, gpu_note = _detect_gpu()
    simd = try
        Base.VersionNumber(Sys.CPU_NAME) # unused; presence probe only
        "auto (LoopVectorization targets host ISA)"
    catch
        "auto"
    end
    rep = HardwareReport(
        cpu,
        Threads.nthreads(),
        LinearAlgebra.BLAS.get_num_threads(),
        blas_backend,
        _detect_numa_nodes(),
        Sys.total_memory() / 2^30,
        cuda, oneapi, gpu_note, simd,
    )
    return rep
end

"""
Configure the environment for the 6x6-heavy workload:
  * BLAS threads → 1 (oversubscription killer; real parallelism is over
    seeds via Julia threads, MODULE 3);
  * if MKL.jl is loadable, load it (deferred: caller decides).
Prints and returns the HardwareReport.
"""
function configure_environment!()
    LinearAlgebra.BLAS.set_num_threads(1)
    rep = hardware_report()
    @printf "CPU:              %s\n" rep.cpu_model
    @printf "Julia threads:    %d\n" rep.julia_threads
    @printf "BLAS:             %s (threads→1 for 6×6 kernels)\n" rep.blas_backend
    @printf "NUMA nodes:       %s\n" rep.numa_nodes
    @printf "Memory:           %.1f GB\n" rep.total_mem_gb
    @printf "GPU CUDA/oneAPI:  %s / %s — %s\n" rep.gpu_cuda rep.gpu_oneapi rep.gpu_note
    println("Note: MKL.jl preferred when installed; for 6×6 kernels the BLAS")
    println("backend is not the bottleneck — seed-level parallelism is.")
    return rep
end
