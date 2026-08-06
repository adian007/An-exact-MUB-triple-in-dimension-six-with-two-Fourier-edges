# Use baremodule to shave off a few KB from the serialized `.ji` file
baremodule FLINT_jll
using Base
using Base: UUID
import JLLWrappers

JLLWrappers.@generate_main_file_header("FLINT")
JLLWrappers.@generate_main_file("FLINT", Base.UUID("e134572f-a0d5-539d-bddf-3cad8db41a82"))
end  # module FLINT_jll
