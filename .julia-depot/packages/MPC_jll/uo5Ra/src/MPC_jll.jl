# Use baremodule to shave off a few KB from the serialized `.ji` file
baremodule MPC_jll
using Base
using Base: UUID
import JLLWrappers

JLLWrappers.@generate_main_file_header("MPC")
JLLWrappers.@generate_main_file("MPC", Base.UUID("2ce0c516-f11f-5db3-98ad-e0e1048fbd70"))
end  # module MPC_jll
