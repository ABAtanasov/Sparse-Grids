include("DG_Methods.jl")

#------------------------------------------------------
# Precompute relevant 1-D gaps
#------------------------------------------------------

precomputed_vals = Dict{NTuple{2,Int},Array{Float64,1}}()

for k in 2:K_max
	for f_number in 1:k
		a = v(k,1,1,f_number,0)
		b = v(k,1,1,f_number,0.5)
		precomputed_vals[(k,f_number)] = [a, (-1)^(f_number) * b, b, (-1)^(f_number) * a]
	end
end

function gap_vals1D(k::Int, level::Int, f_number::Int)
	return (2.0)^((level-1)/2) * precomputed_vals[(k,f_number)]
end