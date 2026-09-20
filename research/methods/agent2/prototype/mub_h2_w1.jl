# NOT RUN by Agent 2.  Small API/semantics benchmark only.
# For I and H2=[[1,1],[1,-1]], gauge z0=w0=1.
# Physical locus is w=conj(z); complexification treats z,w independently.
using HomotopyContinuation

@var z w
F = System([z*w - 1, (1 + z)*(1 + w) - 2]; variables=[z,w])
result = solve(F; show_progress=false)
roots = solutions(result)
cert = certify(F, result; show_progress=false, threading=false)

@assert length(roots) == 2
@assert all(abs(r[2] - conj(r[1])) < 1e-10 for r in roots)
@assert ncertified(cert) == 2
println((roots=roots, certified=ncertified(cert)))
