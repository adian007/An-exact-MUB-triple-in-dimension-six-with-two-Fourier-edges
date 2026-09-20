# Bounded fixed-parameter precision experiment; NOT arbitrary-precision HC.
using LinearAlgebra, Random, Graphs, HomotopyContinuation, Dates, SHA
const ROOT = normpath(joinpath(@__DIR__, "..", ".."))
include(joinpath(ROOT, "src", "mub_zauner_6d_liang_chen.jl"))
const CB = Complex{BigFloat}

# Explicit current theta=0 seam convention, not a theta->0 limit.
function seam_matrix(lam::BigFloat)
    r = sqrt(BigFloat(3))/2
    a = -BigFloat(1)/2 + im*r
    A = [a conj(a); a -conj(a)]
    F = BigFloat[1 1; 1 -1]
    B = -F-A
    zl(z) = [1 1; z -z]
    zr(z) = [1 z; 1 -z]
    Z1 = zl(cis(lam)); Z2 = zl(sqrt(A[1,2]^2/A[1,1]^2))
    Z3 = zr(one(CB)); Z4 = zr(one(CB))
    [F Z1 Z2; Z3 Z3*A*Z1/2 Z3*B*Z2/2; Z4 Z4*B*Z1/2 Z4*A*Z2/2]
end
function phase_rj(t, H)
    z = [one(CB); cis.(t)]
    s = H' * z
    r = abs2.(s[1:5]) .- 6
    J = [2real(conj(s[k])*conj(H[j+1,k])*im*z[j+1]) for k=1:5,j=1:5]
    r, J, z
end
function defects(v, H)
    b = BigFloat(1)/6
    max(maximum(abs.(abs2.(v).-b)), maximum(abs.(abs2.(H'*v/sqrt(BigFloat(6))).-b)))
end
function polynomial_defect(v,H)
    z = v/v[1]; w = conj.(z)
    max(maximum(abs.(z.*w.-1)), maximum(abs.((H'*z).*(transpose(H)*w).-6)))
end
function refine(v,H,bits)
    t = angle.(CB.(v[2:6]./v[1]))
    goal = BigFloat(2)^(-div(bits,2))
    iterations = 0; status = "iteration_limit"; laststep = BigFloat(0)
    for it=1:25
        r,J,z = phase_rj(t,H); iterations=it
        if maximum(abs.(r)) < goal
            status="converged"; break
        end
        d = try J\r catch; status="singular"; break end
        all(isfinite,d) || (status="nonfinite"; break)
        accepted=false
        for back=0:10
            step = d/(BigFloat(2)^back)
            trial=t-step
            rr,_,_ = phase_rj(trial,H)
            if maximum(abs.(rr)) < maximum(abs.(r))
                t=trial; laststep=maximum(abs.(step)); accepted=true; break
            end
        end
        if !accepted; status="stagnation"; break; end
    end
    r,J,z = phase_rj(t,H)
    out=z/sqrt(BigFloat(6))
    full=polynomial_defect(out,H)
    ok=full < 10goal && defects(out,H)<goal
    condition_proxy=try opnorm(J,Inf)*opnorm(inv(J),Inf) catch; BigFloat(Inf) end
    (v=out,ok=ok,status=status,iterations=iterations,poly=full,
     mu=defects(out,H),step=laststep,condition_proxy=condition_proxy,
     displacement=maximum(abs.(out-CB.(v))))
end
function graph_stats(pool,tol)
    g=SimpleGraph(length(pool))
    for i=1:length(pool),j=i+1:length(pool)
        abs(dot(pool[i],pool[j])) < tol && add_edge!(g,i,j)
    end
    cs=maximal_cliques(g)
    isempty(cs) && return (nedge=ne(g),maxclique=0,n6=0,witness=Int[],ortho=BigFloat(NaN))
    c=cs[argmax(length.(cs))]
    ortho=length(c)<2 ? BigFloat(0) : maximum(abs(dot(pool[c[i]],pool[c[j]])) for i=1:length(c),j=i+1:length(c))
    (nedge=ne(g),maxclique=length(c),n6=count(c->length(c)==6,cs),witness=c,ortho=ortho)
end
function selftest(H,lam,bits)
    epsb=BigFloat(2)^(-div(bits,3))
    @assert maximum(abs.(H*H'-6I)) < 100eps(BigFloat)
    @assert maximum(abs.(abs.(H).-1)) < 100eps(BigFloat)
    @assert maximum(abs.(ComplexF64.(H)-build_karlsson_family(0.0,0.5,Float64(lam))))<1e-14
    t=BigFloat.(1:5)/7; r,J,_=phase_rj(t,H)
    for j=1:5
        d=zeros(BigFloat,5);d[j]=epsb
        rp,_,_=phase_rj(t+d,H); rm,_,_=phase_rj(t-d,H)
        @assert maximum(abs.((rp-rm)/(2epsb)-J[:,j])) < sqrt(epsb)
    end
    # Planted orthonormal columns check; MU filtering is tested separately.
    control=[H[:,j]/sqrt(BigFloat(6)) for j=1:6]
    @assert graph_stats(control,parse(BigFloat,"1e-40")).maxclique==6
end

function main()
    length(ARGS)>=2 || error("usage: script LAMBDA NEW_OUTPUT_DIR [SEED]")
    lamstr=ARGS[1]; outdir=abspath(ARGS[2]); seed=length(ARGS)>2 ? parse(Int,ARGS[3]) : 20260917
    isdir(outdir) && error("Refusing existing output directory: $outdir")
    mkpath(outdir)
    open(joinpath(outdir,"metadata.txt"),"w") do io
        println(io,"date=$(now()) julia=$VERSION HC=$(pkgversion(HomotopyContinuation)) seed=$seed")
        println(io,"lambda_decimal=$lamstr theta=0 phi=0.5; current resolved seam; principal square roots")
        println(io,"Float64 HC seeds; independent fixed-H BigFloat phase Newton at 256 and 512 bits; NO HP PATH TRACKING")
        println(io,"source_sha256=$(bytes2hex(sha256(read(joinpath(ROOT,"src","mub_zauner_6d_liang_chen.jl")))))")
        println(io,"script_sha256=$(bytes2hex(sha256(read(@__FILE__))))")
        println(io,"HC_source=$(pathof(HomotopyContinuation))")
    end
    H64=build_karlsson_family(0.0,0.5,parse(Float64,lamstr))
    Random.seed!(seed)
    println("START lambda=$lamstr seed=$seed");flush(stdout)
    t0=time()
    res=solve(build_numeric_pool_system(H64); seed=seed,show_progress=false,threading=false)
    println("HC seconds=$(time()-t0)");flush(stdout)
    open(joinpath(outdir,"hc.txt"),"w") do io
        show(io,MIME("text/plain"),res);println(io)
    end
    pool,stats=pool_from_raw_solutions(solutions(res),H64)
    pool=deduplicate_pool(pool)
    open(joinpath(outdir,"seeds.csv"),"w") do io
        println(io,"id,coordinate,re,im")
        for (i,v) in enumerate(pool),j=1:6
            println(io,"$i,$j,$(real(v[j])),$(imag(v[j]))")
        end
    end
    println("pool=$(length(pool)) stats=$stats");flush(stdout)
    open(joinpath(outdir,"summary.csv"),"w") do summary
        println(summary,"bits,stage,tolerance,n_pool,nedge,max_clique,n_maximal6,witness_ortho,worst_mu,chm_defect")
        for bits in (256,512)
            setprecision(bits) do
                H=seam_matrix(parse(BigFloat,lamstr));selftest(H,parse(BigFloat,lamstr),bits)
                chm=maximum(abs.(H*H'-6I))
                evalpool=[CB.(v) for v in pool]
                refined=Vector{Vector{CB}}()
                open(joinpath(outdir,"refinement_$bits.csv"),"w") do io
                    println(io,"id,accepted,status,iterations,poly,mu,last_step,condition_proxy,displacement")
                    for (i,v) in enumerate(pool)
                        q=refine(v,H,bits)
                        println(io,join((i,q.ok,q.status,q.iterations,q.poly,q.mu,q.step,q.condition_proxy,q.displacement),","))
                        if q.ok && all(maximum(abs.(q.v-w))>parse(BigFloat,"1e-30") for w in refined)
                            push!(refined,q.v)
                        end
                    end
                end
                open(joinpath(outdir,"vectors_$bits.csv"),"w") do io
                    println(io,"id,coordinate,re,im")
                    for (i,v) in enumerate(refined),j=1:6
                        println(io,"$i,$j,$(real(v[j])),$(imag(v[j]))")
                    end
                end
                for (stage,vs) in (("evaluation_only",evalpool),("newton_refined",refined))
                    mu=isempty(vs) ? BigFloat(NaN) : maximum(defects(v,H) for v in vs)
                    for tolstr in ("1e-12","1e-14","1e-16","1e-24","1e-40")
                        tol=parse(BigFloat,tolstr)
                        good=filter(v->defects(v,H)<tol,vs)
                        g=graph_stats(good,tol)
                        println(summary,join((bits,stage,tolstr,length(good),g.nedge,g.maxclique,g.n6,g.ortho,mu,chm),","));flush(summary)
                        println("bits=$bits $stage tol=$tolstr pool=$(length(good)) clique=$(g.maxclique) ortho=$(g.ortho)")
                    end
                end
            end
        end
    end
    println("DONE seconds=$(time()-t0)")
end
main()
