include("DG_Functions.jl")

const K_max = 10;

#Going from Arrays to Polynomials

function array2poly{T<:Real}(v::Array{T},x)
    if abs(x)>1 
		return 0
	else
		n=length(v)
	    k=Int(round(n/2))
	    coeffs1 = [v[i] * x^(i-1) for i in 1:k]
	    coeffs2 = [v[i+k] * f(i-1,x) for i in 1:k]
	    return sum(coeffs1)+sum(coeffs2)
	end
end

function array2poly{T<:Real}(v::Array{T})
    return (x-> array2poly(v,x))
end

#precompute Legendre

Leg_coeffs=legendre(K_max)

function LegendreP(k,x)
    k<=K_max?0:throw(DomainError())
    return array2poly(Leg_coeffs[k+1],x)
end

function LegendreP(k)
    k<=K_max?0:throw(DomainError())
    return array2poly(Leg_coeffs[k+1])
end

#precompute DG functions

DG_coeffs=Array(Any,10)
for i in 1:K_max
    DG_coeffs[i] = DG_Basis(i)
end

function h(k,f_number,x)
    f_number<=k?0:throw(DomainError())
    return array2poly((DG_coeffs[k])[f_number],x)
end

function h(k,f_number)
    f_number<=k?0:throw(DomainError())
    return array2poly((DG_coeffs[k])[f_number])
end

# Need to ask Erik some questions about precomputation