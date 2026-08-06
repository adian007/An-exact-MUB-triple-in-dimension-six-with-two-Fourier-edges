# Minimal Karlsson CHM construction (no HomotopyContinuation dependency).
# Used by formalize_gauge_lemmas.jl when full pipeline cannot load.

function gauge_build_A(theta, phi)
    A11 = -0.5 + im * (sqrt(3) / 2) * (cos(theta) + exp(-im * phi) * sin(theta))
    A12 = -0.5 + im * (sqrt(3) / 2) * (-cos(theta) + exp(im * phi) * sin(theta))
    return [A11 A12; conj(A12) -conj(A11)]
end

function gauge_mobius(z, alpha, beta)
    return (alpha * z - beta) / (conj(beta) * z - conj(alpha))
end

function gauge_build_karlsson_family(theta, phi, lambda_)
    F2 = [1 1; 1 -1]
    A = gauge_build_A(theta, phi)
    B = -F2 - A
    alpha_A, beta_A = A[1, 2]^2, A[1, 1]^2
    alpha_B, beta_B = B[1, 2]^2, B[1, 1]^2
    z1 = exp(im * lambda_)
    z1sq = z1^2
    z3sq = gauge_mobius(z1sq, alpha_A, beta_A)
    z4sq = gauge_mobius(z1sq, alpha_B, beta_B)
    num = beta_B - z3sq * conj(alpha_B)
    den = alpha_B - z3sq * conj(beta_B)
    z2sq = num / den
    z2, z3, z4 = sqrt(z2sq), sqrt(z3sq), sqrt(z4sq)
    Zleft(z) = [1 1; z -z]
    Zright(z) = [1 z; 1 -z]
    Z1, Z2 = Zleft(z1), Zleft(z2)
    Z3, Z4 = Zright(z3), Zright(z4)
    top = hcat(F2, Z1, Z2)
    mid = hcat(Z3, 0.5 * Z3 * A * Z1, 0.5 * Z3 * B * Z2)
    bot = hcat(Z4, 0.5 * Z4 * B * Z1, 0.5 * Z4 * A * Z2)
    return vcat(top, mid, bot)
end

function gauge_dephase(H; tol=1e-8)
    H = copy(H)
    for j in 2:size(H, 2)
        if abs(H[1, j]) > tol
            d = conj(H[1, j]) / abs(H[1, j])
            H[:, j] .*= d
        end
    end
    for i in 2:size(H, 1)
        if abs(H[i, 1]) > tol
            d = conj(H[i, 1]) / abs(H[i, 1])
            H[i, :] .*= d
        end
    end
    return H
end

# Aliases for scripts that include full module when available
const build_A = gauge_build_A
const build_karlsson_family = gauge_build_karlsson_family
const dephase = gauge_dephase
