include("DG_Functions.jl")
include("Specific_DG_Functions.jl")
include("DG_Methods.jl")

function inner_product1D(f::Function, g::Function, lvl::Int, place::Int)
    xmin = (place-1)/(1<<(pos(lvl-1)))
	xmax = (place)/(1<<(pos(lvl-1)))
	h = (x-> f(x)*g(x))
    (val, err) = hquadrature(h, xmin, xmax; reltol=1e-8, abstol=1e-8, maxevals=0)
	return val 
end
#don't need to do anything numerically here, I don't think

function symbolic_diff{T<:Real}(v::Array{T})
	n=length(v)
    k=div(n,2)
	ans = zeros(T,n)
	for i in 1:n
		if i<k
			ans[i] = i*v[i+1]
		elseif i > k && i<2k
			ans[i] = (i-k) * v[i+1]
		else
			ans[i]=0
		end
	end
	return ans
end
#working

function legendreDlegendre(f_number1::Int, f_number2::Int)
	return inner_product(Leg_coeffs[f_number1],symbolic_diff(Leg_coeffs[f_number2]))
end 


function hDh(k::Int, f_number1::Int, f_number2::Int)
	return inner_product(DG_coeffs[k][f_number1], symbolic_diff(DG_coeffs[k][f_number2]))
end

function vDv(k::Int, lvl1::Int, place1::Int, f_number1::Int, lvl2::Int, place2::Int, f_number2::Int)
	if lvl1 == lvl2
		if lvl1 == 0
			return legendreDlegendre(f_number1, f_number2)
		end
		if place1 == place2
			return hDh(k, f_number1, f_number2)*(1<<pos(lvl1-1))
		end
		return 0.0
	end
	if lvl1 < lvl2
		if lvl1 == 0
			return inner_product1D(v(k,0,1,f_number1), v(k,lvl2,place2,f_number2), lvl2, place2)
		end
		if (1<<(lvl2-lvl1))*place1 >= place2 && (1<<(lvl2-lvl1))*(place1-1) < place2
			return inner_product1D(v(k,lvl1,place1,f_number1), v(k,lvl2,place2,f_number2), lvl2, place2)
		end
		return 0.0
	return 0.0
	# if lvl1 == 0 && lvl2==0
	# 	return legendreDlegendre(f_number1,f_number2)
	# elseif lvl1 ==0
	# 	return inner_product(Leg_coeffs[f_number1],symbolic_diff(DG_coeffs[k][f_number2]))
	# elseif lvl2 ==0
	# 	return inner_product(DG_coeffs[k][f_number1],symbolic_diff(Leg_coeffs[f_number2]))
	# end
	# if place1 == place2
	# 	return hDh(k, f_number1, f_number2)
	# return 0
	# #Now we're out of the zone where we have to deal with level 0 Legendre
	# elseif lvl1 > lvl2
# 		#xmin1, xmax1 are close together, since lvl1 is finer
# 	    xmin1 = (place1-1)/(1<<(pos(lvl1-1)))
# 		xmax1 = (place1)/(1<<(pos(lvl1-1)))
# 		#xmin2, xmax2 are further apart, and need to have xmin2 < xmin1 < xmax1 < xmax2
# 	    xmin2 = (place2-1)/(1<<(pos(lvl2-1)))
# 		xmin2 = (place2)/(1<<(pos(lvl2-1)))
# 		if xmin1 < xmin2 || xmax1 > xmax2
# 			return 0.0
# 		end
# 		if xmax11s3d2
#
# 		return inner_product1D(v(k, lvl1, place1, f_number1),
# 									v(k, lvl2, place2, f_number2), lvl1, place1)
#
# 	elseif lvl2 > lvl1
#
#
# 		#return vDv(k,lvl2,place2,f_number2,lvl1,place1,f_number1)
# 		# NOOOO
# 	elseif lvl1==lvl2
# 		if place1==place2
# 			return inner_product(DG_coeffs[k][f_number1],symbolic_diff(DG_coeffs[k][f_number1]))
# 			#return inner_product1D(v(k, lvl1, place1, f_number1),
# 			#						v(k, lvl2, place2, f_number2), lvl1, place1)
# 		end
# 		return 0.0
# 	end
#
end

# function VDV(i::Int, k::Int,
# 	lvl1::NTuple{Int,D}, place1::CartesianIndex{D}, f_number1::CartesianIndex{D},
# 	lvl2::NTuple{Int,D}, place2::CartesianIndex{D}, f_number2::CartesianIndex{D})
# end

function getDerivs{D,T<:Real}(k, coeffs::Dict{CartesianIndex{D}, Array{Array{T},D}}, 
								level::Int, place::Int, f_number)
	f_numbers = ntuple(i -> k, D)
	derivs = Array(T,D)
	for f_number2 in CartesianRange(f_numbers)
		coeffs[level][place][f_number]
	end
end

function D_matrix_sparse(i::Int, n::Int, d::Int)
	D_i = Dict{CartesianIndex{D}, Array{Array{Any},D}}()
	f_numbers= ntuple(i-> k, D)
    ls = ntuple(i->(n+1),D)
	for level in CartesianRange(ls) #This really goes from 0 to l_i for each i 
        diag_level=0;
        for i in 1:D
            diag_level+=level[i]
        end
        if diag_level > n + D #If we're past the levels we care about, don't compute coeffs
            continue
        end
		ks = ntuple(i -> 1<<pos(level[i]-2), D)
		for place in CartesianRange(ks)
			for f_number in CartesianRange(f_numbers)
				#maybe literally make a dictionary applicable only when= overlap = true

			end
		end
	end
end