# Use baremodule to shave off a few KB from the serialized `.ji` file
baremodule SymEngine_jll
using Base
using Base: UUID
import JLLWrappers

JLLWrappers.@generate_main_file_header("SymEngine")
JLLWrappers.@generate_main_file("SymEngine", UUID("3428059b-622b-5399-b16f-d347a77089a4"))
end  # module SymEngine_jll
