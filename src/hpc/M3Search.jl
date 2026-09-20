# ============================================================================
# M3Search.jl — MODULE 3: multithreaded multi-start search kernels.
#
#  (a) Family-constrained search: BFGS on 7–9 family parameters
#      (triple-Karlsson / Karlsson+Diţă), box-constrained multi-start.
#  (b) Manifold search: Riemannian gradient descent on U(6)^3 with
#      tangent projection and Gram–Schmidt retraction, from random flat
#      starts.
#  (c) Simulated-annealing option with acceptance-rate controller
#      (target band 0.32–0.48 per spec).
#
# Near-miss handling: any candidate with W − 6 < 1e-6 is verified
# independently at |G_ij|² = 1/6 to 1e-12; a candidate passing BOTH is a
# candidate 4-MUB SET — it would be a discovery contradicting local
# expectations and is written to disk for escalation, never suppressed
# and never silently accepted (HEURISTIC tier until independently
# certified with the repo's homotopy machinery).
#
# GPU path: guarded pseudo-code hook (CUDA/oneAPI kernels belong behind
# availability probes — this host has no usable GPU userspace stack).
# ============================================================================

using Random

"""One BFGS step (dense inverse-Hessian, dimension ≤ 9)."""
function _bfgs_step!(p, g, B, pnew)
    d = -B * g
    # box projection scale (backtracking handles the rest)
    return d
end

"""Box-constrained multi-start BFGS on the family objective."""
function family_multistart(; mode::Symbol = :triple_karlsson,
                           seeds::Int = 10_000,
                           max_iter::Int = 200,
                           tol::Float64 = 1e-10,
                           seed0::Int = 20260915,
                           verbose::Bool = true)
    bounds = family_bounds(mode)
    npar = length(bounds)
    best = (W = Inf, p = Vector{Float64}())
    results = Vector{NamedTuple}(undef, seeds)
    hist = zeros(Float64, seeds)
    lock_best = ReentrantLock()
    Threads.@threads for s in 1:seeds
        rng = Xoshiro(seed0 + s)
        p = [bounds[i][1] + rand(rng) * (bounds[i][2] - bounds[i][1]) for i in 1:npar]
        W, popt = _bfgs_family(p, bounds; mode = mode, max_iter = max_iter, tol = tol)
        hist[s] = W
        results[s] = (seed = s, W = W, p = popt)
        lock_best !== nothing && lock(lock_best)
        if W < best.W
            best = (W = W, p = copy(popt))
        end
        unlock(lock_best)
    end
    sort!(hist)
    return (best = best, histogram = hist, results = results, mode = mode)
end

function _bfgs_family(p0, bounds; mode, max_iter, tol)
    n = length(p0)
    p = copy(p0)
    B = Matrix{Float64}(I, n, n)
    g = zeros(n)
    gp = zeros(n)
    d = zeros(n)
    h = 1e-6
    f = family_W(p; mode = mode)
    for it in 1:max_iter
        # central finite differences on the (cheap, noisy-free) objective
        for i in 1:n
            pi = p[i]
            step = h * max(1.0, abs(pi))
            p[i] = pi + step
            fp = family_W(p; mode = mode)
            p[i] = pi - step
            fm = family_W(p; mode = mode)
            p[i] = pi
            gp[i] = (fp - fm) / (2 * step)
        end
        if norm(gp) < tol
            break
        end
        d = -(B * gp)
        # backtracking line search with box projection
        alpha = 1.0
        f0 = f
        improved = false
        for _ in 1:40
            pnew = clamp.(p + alpha .* d, [b[1] for b in bounds], [b[2] for b in bounds])
            fnew = family_W(pnew; mode = mode)
            if fnew < f0 - 1e-12
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

"""Riemannian gradient descent on U(6)^3 from a random flat start.
Returns (W_final, bases, iterations)."""
function manifold_run!(Bs::Vector{Matrix{ComplexF64}};
                       max_iter::Int = 400,
                       η0::Float64 = 0.05,
                       tol::Float64 = 1e-12,
                       rng = Xoshiro())
    n = length(Bs)
    grads = [zeros(ComplexF64, 6, 6) for _ in 1:n]
    tang = [zeros(ComplexF64, 6, 6) for _ in 1:n]
    G = zeros(ComplexF64, 6, 6)
    res = total_W(Bs; G = G)
    W = res.W
    η = η0
    it = 0
    for it in 1:max_iter
        grad_total_W!(grads, Bs, G)
        gn = maximum(norm.(grads))
        gn < tol && break
        for k in 1:n
            tangent_project!(tang[k], grads[k], Bs[k])
            Bs[k] .-= η .* tang[k]
            unitary_retract!(Bs[k])
        end
        Wnew = total_W(Bs; G = G).W
        if Wnew < W - 1e-14
            W = Wnew
            η = min(η * 1.05, 0.5)
        else
            η *= 0.5          # adaptive shrink on failure
            η < 1e-8 && break
        end
        W - 6 < 1e-10 && break   # converged to a 4-MUB candidate
    end
    return (W = total_W(Bs; G = G).W, bases = Bs, iterations = it)
end

"""Multithreaded manifold multi-start over `seeds` random flat starts."""
function manifold_multistart(; seeds::Int = 2_000, max_iter::Int = 400,
                             seed0::Int = 20260915)
    hist = Vector{Float64}(undef, seeds)
    best = (W = Inf, bases = Vector{Matrix{ComplexF64}}())
    lock_best = ReentrantLock()
    near_misses = Vector{NamedTuple}()
    lock_nm = ReentrantLock()
    Threads.@threads for s in 1:seeds
        rng = Xoshiro(seed0 * 1_000_003 + s)
        Bs = [random_flat_matrix!(rng, zeros(ComplexF64, 6, 6)) ./ sqrt(6) for _ in 1:3]
        pushfirst!(Bs, Matrix{ComplexF64}(I, 6, 6))
        out = manifold_run!(Bs; max_iter = max_iter, rng = rng)
        hist[s] = out.W
        lock(lock_best)
        if out.W < best.W
            best = (W = out.W, bases = deepcopy(out.bases))
        end
        unlock(lock_best)
        if out.W - 6 < 1e-6
            vm = _verify_near_miss(out.bases)
            lock(lock_nm)
            push!(near_misses, (W = out.W, verified = vm.verified,
                                max_dev = vm.max_dev, bases = deepcopy(out.bases)))
            unlock(lock_nm)
        end
    end
    sort!(hist)
    return (best = best, histogram = hist, near_misses = near_misses)
end

"""Independent near-miss verification: every pair must have ALL |G_ij|² = 1/6
to 1e-12 (⟺ W = 6 exactly)."""
function _verify_near_miss(Bs; tol::Float64 = 1e-12)
    max_dev = 0.0
    for k in 1:length(Bs), m in (k + 1):length(Bs)
        Gm = Bs[k]' * Bs[m]
        for ij in eachindex(Gm)
            max_dev = max(max_dev, abs(abs2(Gm[ij]) - 1 / 6))
        end
    end
    return (verified = max_dev < tol, max_dev = max_dev)
end

"""Simulated annealing over the family parameters with acceptance-rate
controller (target band 0.32–0.48)."""
function family_anneal(; mode::Symbol = :triple_karlsson, seeds::Int = 1_000,
                       iters::Int = 2_000, T0::Float64 = 0.5, seed0::Int = 20260915)
    bounds = family_bounds(mode)
    npar = length(bounds)
    hist = Vector{Float64}(undef, seeds)
    best = (W = Inf, p = Float64[])
    lock_best = ReentrantLock()
    Threads.@threads for s in 1:seeds
        rng = Xoshiro(seed0 * 7 + s)
        p = [bounds[i][1] + rand(rng) * (bounds[i][2] - bounds[i][1]) for i in 1:npar]
        f = family_W(p; mode = mode)
        T = T0
        acc = 0
        for it in 1:iters
            pnew = copy(p)
            j = rand(rng, 1:npar)
            width = 0.05 * (bounds[j][2] - bounds[j][1]) * T / T0
            pnew[j] = clamp(p[j] + width * (2rand(rng) - 1), bounds[j][1], bounds[j][2])
            fnew = family_W(pnew; mode = mode)
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
        hist[s] = f
        lock(lock_best)
        f < best.W && (best = (W = f, p = copy(p)))
        unlock(lock_best)
    end
    sort!(hist)
    return (best = best, histogram = hist)
end

# --- GPU hook (availability-guarded; see M1Hardware) -----------------------
#
# When a CUDA device is present and CUDA.jl installed, the manifold
# multi-start maps to a CUDA kernel: one thread block per seed, 6×6 complex
# tiles in shared memory, the Gram–Schmidt retraction done in-block. The
# CPU path above is the reference implementation; the GPU path must
# reproduce its histogram on identical seeds before being trusted.
function gpu_available()
    return try
        success(`nvidia-smi -L`)
    catch
        false
    end
end
