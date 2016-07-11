#------------------------------------------------------
# Differentiate a function along the ith component
#------------------------------------------------------
function diffOneVariable(f::Function, l::Int, x::Real)
	h = 1.0/(1<<l)
	if x< h
		x+=h
		after = f(x)
		x-=h
		before = f(x)
		return (after - before)*(1<<l)
	elseif x > 1-1/(1<<l)
		x-=h
		before = f(x)
		x+=h
		after = f(x)
		return (after - before)*(1<<l)
	end
	x+=h
	after = f(x)
	x-=2*h 
	before = f(x)
	x+=h
	return (after - before)*(1<<(l-1))
end
#working

function diff{T<:Real}(i::Int, f::Function, l::Int, x::Array{T})
	h = 1.0/(1<<l)
	if x[i]< h
		x[i]+=h
		after = f(x)
		x[i]-=h
		before = f(x)
		return (after - before)*(1<<l)
	elseif x[i] > 1-1/(1<<l)
		x[i]-=h
		before = f(x)
		x[i]+=h
		after = f(x)
		return (after - before)*(1<<l)
	end
	x[i]+=h
	after = f(x)
	x[i]-=2*h 
	before = f(x)
	x[i]+=h
	return (after - before)*(1<<(l-1))
end
#working


#Full gradient, returning an NTuple:
function grad{T<:Real}(f::Function, l::Int, x::Array{T})
	return ntuple(i-> diff(i,f,l,x))
end