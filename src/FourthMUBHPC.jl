module FourthMUBHPC

using LinearAlgebra
using Random
using Statistics
using Printf
using Dates
using Base.Threads

# Include the static kernels and family builders
include(joinpath(@__DIR__, "StaticMUBKernels.jl"))
include(joinpath(@__DIR__, "brierley_weigert_notes.jl"))

# ================================================================
# MODULE 1: Hardware Initialization & Environment Setup
# ================================================================
module Hardware
    export detect_hardware, configure_blas

    function detect_hardware()
        cpu_count = Sys.CPU_THREADS
        thread_count = Threads.nthreads()
        gpu_available = false
        gpu_backend = "none"

        try
            if Base.find_package("CUDA") !== nothing
                gpu_available = true
                gpu_backend = "CUDA"
            end
        catch
        end

        try
            if !gpu_available && Base.find_package("oneAPI") !== nothing
                gpu_available = true
                gpu_backend = "oneAPI"
            end
        catch
        end

        return (
            cpu_cores = cpu_count,
            julia_threads = thread_count,
            gpu_available = gpu_available,
            gpu_backend = gpu_backend,
            total_mem_gb = round(Sys.total_memory() / 2^30; digits = 1),
            timestamp = string(now())
        )
    end

    function configure_blas(; target_threads::Int = 1)
        # For 6×6 kernels, BLAS threads MUST be 1 — real parallelism is
        # over seeds via Julia threads.
        target_threads = max(1, target_threads)
        BLAS.set_num_threads(target_threads)
        return BLAS.get_num_threads()
    end
end

# ================================================================
# MODULE 2: Static Manifold & Objective Function Definitions
# ================================================================
module Manifold
    export identity_basis, random_start, manifold_search,
           total_W, pair_mubness,
           build_dita_family, family_W, family_bounds,
           family_bases_karlsson_triple, family_bases_karlsson_dita

    using ..StaticMUBKernels: SMat6, SVec6, MMat6,
        pair_mubness_static, total_W_static,
        grad_total_W_static, tangent_project_static,
        unitary_retract_static, gauge_fix_static,
        random_unitary_static, fourier_matrix_static,
        identity_basis_static, manifold_descent_static

    using LinearAlgebra
    using Random

    # Re-export the Diţă builder from the parent module
    using ..FourthMUBHPC: dita_D

    const DIM = 6

    # --- Wrappers around static kernels ---

    identity_basis() = identity_basis_static()

    function pair_mubness(Bk::SMat6, Bm::SMat6)
        return pair_mubness_static(Bk, Bm)
    end

    function total_W(B0::SMat6, B1::SMat6, B2::SMat6, B3::SMat6)
        return total_W_static(B0, B1, B2, B3)
    end

    """Generate a random search seed: B0 = I, B1..B3 = random unitaries."""
    function random_start(; rng::AbstractRNG = Random.default_rng())
        B0 = identity_basis_static()
        B1 = random_unitary_static(rng)
        B2 = random_unitary_static(rng)
        B3 = random_unitary_static(rng)
        return (B0 = B0, B1 = B1, B2 = B2, B3 = B3)
    end

    """Riemannian gradient descent on U(6)^3 from a given starting point."""
    function manifold_search(B0::SMat6, B1::SMat6, B2::SMat6, B3::SMat6;
                             max_iter::Int = 400, η0::Float64 = 0.05,
                             tol::Float64 = 1e-12,
                             rng::AbstractRNG = Random.default_rng())
        return manifold_descent_static(B0, B1, B2, B3;
                                       max_iter = max_iter, η0 = η0,
                                       tol = tol, rng = rng)
    end

    # --- Family-constrained search (Karlsson & Diţă) ---

    # Karlsson A-block: the audited, load-bearing transcription
    # (from Karlsson LAA 434, 2011; arXiv:1003.4177)
    const F2_HAD = ComplexF64[1 1; 1 -1]

    function build_karlsson_A(theta::Float64, phi::Float64)
        A11 = -0.5 + im * (sqrt(3) / 2) * (cos(theta) + exp(-im * phi) * sin(theta))
        A12 = -0.5 + im * (sqrt(3) / 2) * (-cos(theta) + exp(im * phi) * sin(theta))
        return ComplexF64[A11 A12; conj(A12) -conj(A11)]
    end

    """
    Build a Karlsson K6(θ,φ,λ) complex Hadamard matrix (unnormalised, |H_ij|=1).
    Uses the Möbius-map construction.
    """
    function build_karlsson_chm(theta::Float64, phi::Float64, lambda::Float64)
        A = build_karlsson_A(theta, phi)
        B = -F2_HAD - A

        alpha_A, beta_A = A[1, 2]^2, A[1, 1]^2
        alpha_B, beta_B = B[1, 2]^2, B[1, 1]^2

        z1sq = exp(2im * lambda)
        z1 = exp(im * lambda)

        if abs(theta) < 1e-14
            # Fourier seam: algebraic resolution
            z3sq = one(ComplexF64)
            z4sq = one(ComplexF64)
            z2sq = alpha_A / beta_A
        else
            # Möbius maps
            den_A = conj(beta_A) * z1sq - conj(alpha_A)
            z3sq = abs(den_A) > 1e-12 ? (alpha_A * z1sq - beta_A) / den_A : one(ComplexF64)

            den_B = conj(beta_B) * z1sq - conj(alpha_B)
            z4sq = abs(den_B) > 1e-12 ? (alpha_B * z1sq - beta_B) / den_B : one(ComplexF64)

            # z2² from z3² = M_B(z2²)
            num2 = beta_B - z3sq * conj(alpha_B)
            den2 = alpha_B - z3sq * conj(beta_B)
            z2sq = abs(den2) > 1e-12 ? num2 / den2 : one(ComplexF64)
        end

        # Build the 6×6 CHM
        Zleft(z)  = ComplexF64[1 1; z -z]
        Zright(z) = ComplexF64[1 z; 1 -z]

        z2 = sqrt(z2sq)
        z3 = sqrt(z3sq)
        z4 = sqrt(z4sq)

        Z1 = Zleft(z1)
        Z2 = Zleft(z2)
        Z3 = Zright(z3)
        Z4 = Zright(z4)

        top = hcat(F2_HAD, Z1, Z2)
        mid = hcat(Z3, 0.5 * Z3 * A * Z1, 0.5 * Z3 * B * Z2)
        bot = hcat(Z4, 0.5 * Z4 * B * Z1, 0.5 * Z4 * A * Z2)
        H = vcat(top, mid, bot)

        return ComplexF64.(H)
    end

    """Convert an unnormalised CHM (|H_ij|=1, HH†=6I) to an orthonormal basis."""
    function chm_to_basis(H::Matrix{ComplexF64})
        return SMat6(H ./ sqrt(6))
    end

    """
    Family objective W as a function of parameters.
    mode :triple_karlsson — p = [θ₁ φ₁ λ₁ θ₂ φ₂ λ₂ θ₃ φ₃ λ₃] (9 params)
    mode :karlsson_dita   — p = [θ₁ φ₁ λ₁ x₂ θ₃ φ₃ λ₃]       (7 params)
    Returns W; penalty 72.0 for invalid parameters.
    """
    function family_W(p::Vector{Float64}; mode::Symbol = :triple_karlsson)
        try
            Bs = _family_bases(p; mode = mode)
            return total_W_static(Bs[1], Bs[2], Bs[3], Bs[4]).W
        catch
            return 72.0  # penalty: 2 × W_max
        end
    end

    function _family_bases(p::Vector{Float64}; mode::Symbol = :triple_karlsson)
        B0 = identity_basis_static()
        if mode === :triple_karlsson
            length(p) == 9 || error("triple_karlsson needs 9 parameters")
            H1 = build_karlsson_chm(p[1], p[2], p[3])
            H2 = build_karlsson_chm(p[4], p[5], p[6])
            H3 = build_karlsson_chm(p[7], p[8], p[9])
            return (B0, chm_to_basis(H1), chm_to_basis(H2), chm_to_basis(H3))
        elseif mode === :karlsson_dita
            length(p) == 7 || error("karlsson_dita needs 7 parameters")
            H1 = build_karlsson_chm(p[1], p[2], p[3])
            H2 = dita_D(p[4])
            H3 = build_karlsson_chm(p[5], p[6], p[7])
            return (B0, chm_to_basis(H1), chm_to_basis(Matrix{ComplexF64}(H2)),
                    chm_to_basis(H3))
        else
            error("unknown family mode $mode")
        end
    end

    """Box bounds for the family multi-start."""
    function family_bounds(mode::Symbol = :triple_karlsson)
        kbox = [(1e-2, π - 1e-2), (0.0, 2π), (0.0, 2π)]
        if mode === :triple_karlsson
            return reduce(vcat, [kbox for _ in 1:3])
        else
            return vcat(kbox, [(0.0, 1.0)], kbox)
        end
    end
end

# ================================================================
# MODULE 3: Multi-Threaded Optimization Kernels
# ================================================================
module Optimizer
    export run_manifold_search, run_family_search, run_anneal_search,
           run_combined_search, SearchReport, make_report

    using ..Manifold
    using ..StaticMUBKernels: SMat6, identity_basis_static, random_unitary_static
    using Random
    using Statistics
    using Base.Threads
    using LinearAlgebra
    using Printf

    struct SearchReport
        strategy::Symbol
        seeds::Int
        best_W::Float64
        near_mub_count::Int     # seeds with W - 6 < tol
        elapsed_sec::Float64
        threads::Int
        best_seed::Int
        best_params::Vector{Float64}  # empty for manifold, family params otherwise
    end

    # ------------------------------------------------------------------
    # Strategy 1: Manifold multi-start (Riemannian gradient descent)
    # ------------------------------------------------------------------

    function run_manifold_search(num_seeds::Int = 2_000;
                                 max_iter::Int = 400,
                                 tol::Float64 = 1e-6,
                                 seed0::Int = 20260915)
        nthreads = Threads.nthreads()
        BLAS.set_num_threads(1)

        local_best_W = fill(Inf, nthreads)
        local_best_seed = fill(0, nthreads)
        local_hits = fill(0, nthreads)

        t0 = time()

        Threads.@threads :static for tid in 1:nthreads
            rng = Xoshiro(seed0 * 1_000_003 + tid)
            n_each = num_seeds ÷ nthreads
            rem = num_seeds % nthreads
            my_count = n_each + (tid <= rem ? 1 : 0)

            for s in 1:my_count
                seed_id = (tid - 1) * n_each + min(tid, rem) + s - (tid <= rem ? 1 : 0)
                state = Manifold.random_start(rng = rng)
                result = Manifold.manifold_search(
                    state.B0, state.B1, state.B2, state.B3;
                    max_iter = max_iter, rng = rng)

                W = result.W
                if W < local_best_W[tid]
                    local_best_W[tid] = W
                    local_best_seed[tid] = seed_id
                end
                if abs(W - 6.0) < tol
                    local_hits[tid] += 1
                end
            end
        end

        best_idx = argmin(local_best_W)
        elapsed = time() - t0

        return SearchReport(
            :manifold,
            num_seeds,
            local_best_W[best_idx],
            sum(local_hits),
            elapsed,
            nthreads,
            local_best_seed[best_idx],
            Float64[]
        )
    end

    # ------------------------------------------------------------------
    # Strategy 2: Family-constrained BFGS multi-start
    # ------------------------------------------------------------------

    function run_family_search(num_seeds::Int = 5_000;
                               mode::Symbol = :triple_karlsson,
                               max_iter::Int = 200,
                               tol::Float64 = 1e-6,
                               seed0::Int = 20260915)
        nthreads = Threads.nthreads()
        BLAS.set_num_threads(1)

        bounds = Manifold.family_bounds(mode)
        npar = length(bounds)

        local_best_W = fill(Inf, nthreads)
        local_best_seed = fill(0, nthreads)
        local_best_p = [Float64[] for _ in 1:nthreads]
        local_hits = fill(0, nthreads)

        t0 = time()

        Threads.@threads :static for tid in 1:nthreads
            rng = Xoshiro(seed0 + tid)
            n_each = num_seeds ÷ nthreads
            rem = num_seeds % nthreads
            my_count = n_each + (tid <= rem ? 1 : 0)

            for s in 1:my_count
                seed_id = (tid - 1) * n_each + min(tid, rem) + s - (tid <= rem ? 1 : 0)
                p = [bounds[i][1] + rand(rng) * (bounds[i][2] - bounds[i][1]) for i in 1:npar]
                W, popt = _bfgs_family(p, bounds; mode = mode, max_iter = max_iter)

                if W < local_best_W[tid]
                    local_best_W[tid] = W
                    local_best_seed[tid] = seed_id
                    local_best_p[tid] = copy(popt)
                end
                if abs(W - 6.0) < tol
                    local_hits[tid] += 1
                end
            end
        end

        best_idx = argmin(local_best_W)
        elapsed = time() - t0

        return SearchReport(
            :family,
            num_seeds,
            local_best_W[best_idx],
            sum(local_hits),
            elapsed,
            nthreads,
            local_best_seed[best_idx],
            local_best_p[best_idx]
        )
    end

    """Box-constrained BFGS on the family objective (dimension ≤ 9)."""
    function _bfgs_family(p0::Vector{Float64}, bounds;
                          mode::Symbol = :triple_karlsson,
                          max_iter::Int = 200)
        n = length(p0)
        p = copy(p0)
        B = Matrix{Float64}(I, n, n)  # inverse Hessian approx
        g = zeros(n)
        gp = zeros(n)
        h = 1e-6
        f = Manifold.family_W(p; mode = mode)

        for it in 1:max_iter
            # Central finite differences
            for i in 1:n
                pi = p[i]
                step = h * max(1.0, abs(pi))
                p[i] = pi + step
                fp = Manifold.family_W(p; mode = mode)
                p[i] = pi - step
                fm = Manifold.family_W(p; mode = mode)
                p[i] = pi
                gp[i] = (fp - fm) / (2 * step)
            end
            norm(gp) < 1e-10 && break

            d = -(B * gp)

            # Backtracking line search with box projection
            alpha = 1.0
            improved = false
            for _ in 1:40
                pnew = clamp.(p + alpha .* d,
                              [b[1] for b in bounds],
                              [b[2] for b in bounds])
                fnew = Manifold.family_W(pnew; mode = mode)
                if fnew < f - 1e-12
                    svec = pnew .- p
                    yvec = gp .- g
                    sy = dot(svec, yvec)
                    if sy > 1e-14
                        Bs = B * svec
                        B .= B .- (Bs * Bs') ./ dot(svec, Bs) .+ (yvec * yvec') ./ sy
                    end
                    p .= pnew
                    f = fnew
                    g .= gp
                    improved = true
                    break
                end
                alpha *= 0.5
            end
            improved || break
        end
        return f, p
    end

    # ------------------------------------------------------------------
    # Strategy 3: Simulated annealing over family parameters
    # ------------------------------------------------------------------

    function run_anneal_search(num_seeds::Int = 2_000;
                               mode::Symbol = :triple_karlsson,
                               iters::Int = 2_000,
                               T0::Float64 = 0.5,
                               tol::Float64 = 1e-6,
                               seed0::Int = 20260915)
        nthreads = Threads.nthreads()
        BLAS.set_num_threads(1)

        bounds = Manifold.family_bounds(mode)
        npar = length(bounds)

        local_best_W = fill(Inf, nthreads)
        local_best_seed = fill(0, nthreads)
        local_best_p = [Float64[] for _ in 1:nthreads]
        local_hits = fill(0, nthreads)

        t0 = time()

        Threads.@threads :static for tid in 1:nthreads
            rng = Xoshiro(seed0 * 7 + tid)
            n_each = num_seeds ÷ nthreads
            rem = num_seeds % nthreads
            my_count = n_each + (tid <= rem ? 1 : 0)

            for s in 1:my_count
                seed_id = (tid - 1) * n_each + min(tid, rem) + s - (tid <= rem ? 1 : 0)
                p = [bounds[i][1] + rand(rng) * (bounds[i][2] - bounds[i][1]) for i in 1:npar]
                f = Manifold.family_W(p; mode = mode)
                T = T0
                acc = 0

                for it in 1:iters
                    pnew = copy(p)
                    j = rand(rng, 1:npar)
                    width = 0.05 * (bounds[j][2] - bounds[j][1]) * T / T0
                    pnew[j] = clamp(p[j] + width * (2rand(rng) - 1),
                                    bounds[j][1], bounds[j][2])
                    fnew = Manifold.family_W(pnew; mode = mode)
                    if fnew < f || rand(rng) < exp((f - fnew) / T)
                        p, f = pnew, fnew
                        acc += 1
                    end
                    if it % 50 == 0
                        rate = acc / 50
                        T *= rate > 0.48 ? 1.15 : rate < 0.32 ? 0.85 : 0.97
                        acc = 0
                        T = max(T, 1e-6)
                    end
                end

                if f < local_best_W[tid]
                    local_best_W[tid] = f
                    local_best_seed[tid] = seed_id
                    local_best_p[tid] = copy(p)
                end
                if abs(f - 6.0) < tol
                    local_hits[tid] += 1
                end
            end
        end

        best_idx = argmin(local_best_W)
        elapsed = time() - t0

        return SearchReport(
            :anneal,
            num_seeds,
            local_best_W[best_idx],
            sum(local_hits),
            elapsed,
            nthreads,
            local_best_seed[best_idx],
            local_best_p[best_idx]
        )
    end

    # ------------------------------------------------------------------
    # Strategy 4: Combined search (all three strategies)
    # ------------------------------------------------------------------

    function run_combined_search(; manifold_seeds::Int = 1_000,
                                  family_seeds::Int = 2_000,
                                  anneal_seeds::Int = 1_000,
                                  family_mode::Symbol = :triple_karlsson)
        results = SearchReport[]

        r1 = run_manifold_search(manifold_seeds)
        push!(results, r1)

        r2 = run_family_search(family_seeds; mode = family_mode)
        push!(results, r2)

        r3 = run_anneal_search(anneal_seeds; mode = family_mode)
        push!(results, r3)

        best = results[argmin([r.best_W for r in results])]

        return (strategies = results, best = best)
    end

    function make_report(report::SearchReport)
        return (
            strategy = report.strategy,
            seeds = report.seeds,
            best_W = report.best_W,
            near_mub = report.near_mub_count,
            elapsed_sec = report.elapsed_sec,
            threads = report.threads,
            best_seed = report.best_seed,
            has_params = !isempty(report.best_params)
        )
    end
end

# ================================================================
# MODULE 4: High-Precision SDP / NPA Hierarchy Solver (stub)
# ================================================================
module SDP
    export solve_npa_oracle, available_sdps, npa_hierarchy_stub

    function available_sdps()
        available = String[]
        for name in ("COSMO", "SCS", "MOSEK", "JuMP")
            try
                Base.find_package(name) !== nothing && push!(available, name)
            catch
            end
        end
        return available
    end

    function npa_hierarchy_stub(; level::Int = 2, symmetry_reduction = true)
        available = available_sdps()
        if isempty(available)
            return (
                status = :not_available,
                message = "No SDP backend detected. Install COSMO.jl, SCS.jl, or MOSEK and JuMP.jl to enable NPA level 2/3 solves.",
                level = level,
                symmetry_reduction = symmetry_reduction,
                feasible = false
            )
        end

        return (
            status = :stub,
            message = "Prepared an NPA/Lasserre-style moment-matrix model for a 4-MUB infeasibility test. The exact solve requires a dedicated high-precision backend.",
            level = level,
            symmetry_reduction = symmetry_reduction,
            feasible = false,
            available_backends = available
        )
    end

    function solve_npa_oracle(; level::Int = 2, solver::Symbol = :auto)
        if solver === :auto
            available = available_sdps()
            if isempty(available)
                return npa_hierarchy_stub(level = level)
            end
            solver = :stub
        end

        if solver === :stub
            return npa_hierarchy_stub(level = level)
        end

        return (
            status = :requested,
            message = "Requested solver $(solver) is configured externally; exact SDPs are not run in this lightweight demo build.",
            level = level,
            feasible = false
        )
    end
end

# ================================================================
# MODULE 5: Benchmarking & Diagnostics Reporting
# ================================================================
module Diagnostics
    export print_summary, print_combined_summary, benchmark_suite

    using ..Hardware
    using ..Manifold
    using ..Optimizer
    using ..SDP
    using Printf
    using Dates

    function print_summary(report::Optimizer.SearchReport)
        println("="^80)
        @printf "Fourth-MUB Search Summary  [strategy: %s]\n" report.strategy
        println("time      : ", Dates.format(now(), "yyyy-mm-dd HH:MM:SS"))
        @printf "best W    : %.10f  (target: 6.0 = MUBs)\n" report.best_W
        @printf "W − 6     : %.2e\n" (report.best_W - 6.0)
        @printf "near-MUB  : %d seeds\n" report.near_mub_count
        @printf "threads   : %d\n" report.threads
        @printf "elapsed   : %.3f s\n" report.elapsed_sec
        @printf "seeds/sec : %.0f\n" (report.seeds / max(report.elapsed_sec, 1e-6))
        if !isempty(report.best_params)
            @printf "best params: %s\n" join([@sprintf("%.6f", x) for x in report.best_params], ", ")
        end
        println("="^80)
    end

    function print_combined_summary(result)
        println("\n", "╔", "═"^78, "╗")
        println("║  COMBINED SEARCH RESULTS", " "^53, "║")
        println("╠", "═"^78, "╣")
        for r in result.strategies
            @printf "║  %-12s │ W = %-16.10f │ seeds = %-6d │ %.2f s  ║\n" r.strategy r.best_W r.seeds r.elapsed_sec
        end
        println("╠", "═"^78, "╣")
        b = result.best
        @printf "║  BEST: %-10s │ W = %.10f │ W − 6 = %.2e          ║\n" b.strategy b.best_W (b.best_W - 6.0)
        println("╚", "═"^78, "╝")
    end

    function benchmark_suite(; manifold_seeds = 500, family_seeds = 1_000,
                              anneal_seeds = 500)
        hw = Hardware.detect_hardware()
        Hardware.configure_blas(target_threads = 1)
        println("Detected hardware: ", hw)

        result = Optimizer.run_combined_search(
            manifold_seeds = manifold_seeds,
            family_seeds = family_seeds,
            anneal_seeds = anneal_seeds)

        print_combined_summary(result)

        sdp = SDP.solve_npa_oracle(level = 2)
        println("\nSDP status: ", sdp.status)
        println("SDP message: ", sdp.message)

        return (hardware = hw, result = result, sdp = sdp)
    end
end

# ================================================================
# Public entry points
# ================================================================

"""
    run_demo(; seeds=2000, strategy=:combined)

Run a quick MUB search demonstration.
  strategy ∈ {:manifold, :family, :anneal, :combined}
"""
function run_demo(; seeds::Int = 2_000, strategy::Symbol = :combined)
    hw = Hardware.detect_hardware()
    Hardware.configure_blas(target_threads = 1)
    println("Hardware summary:")
    println(hw)

    if strategy === :combined
        result = Optimizer.run_combined_search(
            manifold_seeds = seeds ÷ 3,
            family_seeds = seeds ÷ 3,
            anneal_seeds = seeds - 2 * (seeds ÷ 3))
        Diagnostics.print_combined_summary(result)
        return result
    elseif strategy === :manifold
        result = Optimizer.run_manifold_search(seeds)
        Diagnostics.print_summary(result)
        return result
    elseif strategy === :family
        result = Optimizer.run_family_search(seeds)
        Diagnostics.print_summary(result)
        return result
    elseif strategy === :anneal
        result = Optimizer.run_anneal_search(seeds)
        Diagnostics.print_summary(result)
        return result
    else
        error("Unknown strategy: $strategy. Use :manifold, :family, :anneal, or :combined")
    end
end

function startup_banner()
    hw = Hardware.detect_hardware()
    println("FourthMUBHPC initialized (laptop-optimized, StaticArrays + Riemannian descent)")
    println("CPU cores      : ", hw.cpu_cores)
    println("Julia threads  : ", hw.julia_threads)
    println("Memory         : ", hw.total_mem_gb, " GB")
    println("GPU available  : ", hw.gpu_available)
    println("GPU backend    : ", hw.gpu_backend)
    println("Objective      : MINIMIZE W → 6.0 (= 4 MUBs in d=6)")
    println("Strategies     : manifold (Riemannian), family (BFGS), anneal (SA)")
    return hw
end

end # module FourthMUBHPC
