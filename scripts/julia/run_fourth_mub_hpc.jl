# run_fourth_mub_hpc.jl — Laptop-optimized MUB search runner.
#
# Usage:
#   julia -t auto scripts/julia/run_fourth_mub_hpc.jl
#
# The '-t auto' flag enables multithreading (uses all available CPU cores).
# Without it, Julia runs single-threaded.

include(joinpath(@__DIR__, "..", "..", "src", "FourthMUBHPC.jl"))
using .FourthMUBHPC

# --- Show hardware and configuration ---
FourthMUBHPC.startup_banner()
println()

# --- Run the combined search (all three strategies) ---
# Adjust seed counts based on your laptop's speed:
#   Quick test:   seeds = 300
#   Normal run:   seeds = 3000
#   Deep search:  seeds = 30000
seeds = parse(Int, get(ARGS, 1, "3000"))
println("Running combined search with $(seeds) seeds total...")
println()

result = FourthMUBHPC.run_demo(seeds = seeds, strategy = :combined)

# --- Individual strategy runs (uncomment to use) ---
# println("\n--- Manifold-only search ---")
# r = FourthMUBHPC.run_demo(seeds = 1000, strategy = :manifold)

# println("\n--- Family BFGS search ---")
# r = FourthMUBHPC.run_demo(seeds = 2000, strategy = :family)

# println("\n--- Simulated annealing search ---")
# r = FourthMUBHPC.run_demo(seeds = 1000, strategy = :anneal)

println("\nDone.")
