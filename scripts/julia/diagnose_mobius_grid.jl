# Diagnostic: which grid points hit singular Möbius branches?
include(joinpath(@__DIR__, "_paths.jl"))
include(joinpath(ROOT, "src", "MubSearch.jl"))
using .MubSearch
using Printf

pts = [(0.3, 0.5, 0.2), (1.0, 2.0, 0.7), (0.1, 0.1, 3.0),
       (DITA_THETA, DITA_PHI, 0.4), (DITA_THETA, DITA_PHI, pi / 2),
       (0.0, 0.5, 0.3), (0.0, 0.5, 0.0), (0.5, 1.0, 1.5)]
for (th, ph, lm) in pts
    ma = karlsson_mobius_audit(th, ph, lm)
    status = try
        H = build_karlsson_family(th, ph, lm)
        u = maximum(abs.(H * H' - 6I(6)))
        isnan(u) ? "NaN-MATRIX" : (u < 1e-10 ? "HADAMARD-OK" : "HADAMARD-BAD")
    catch e
        "THROWN"
    end
    fmt_ma = isnan(ma.id_MA_z2) ? "NaN" : @sprintf("%.4f", ma.id_MA_z2)
    fmt_mb = isnan(ma.id_MB_z2) ? "NaN" : @sprintf("%.4f", ma.id_MB_z2)
    @printf "(%.3f,%.3f,%.3f): idMAz2=%s idMBz2=%s sA1=%s sB1=%s z2sing=%s | %s\n" th ph lm fmt_ma fmt_mb ma.singular_MA_z1 ma.singular_MB_z1 ma.z2_derived_singular status
end
