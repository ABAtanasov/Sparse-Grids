include("../Specific_DG_Functions.jl")
include("../DG_Derivative.jl")

function leg(f_number::Int, x::Real)
    return sqrt(2.0)*LegendreP(f_number, 2*x-1) 
end

function basis(level::Int, place::Int, f_number::Int, x::Real)
    return leg(f_number, (1<<level)*x - (place-1)) * (2.0)^(level/2)
end

function basis(level::Int, place::Int, f_number::Int)
    return x->basis(level, place, f_number, x)
end

function dleg(f_number::Int, x::Real)
    return sqrt(2.0)*dLegendreP(f_number, 2*x-1) *2
end

function dbasis(level::Int, place::Int, f_number::Int, x::Real)
    return dleg(f_number, (1<<level)*x - (place-1)) * (2.0)^(level/2) * (1<<level)
end

function dbasis(level::Int, place::Int, f_number::Int)
    return x->dbasis(level, place, f_number, x)
end



function get_coeffs(level::Int, k::Int, f::Function)
    coeffs = Array(Float64,(1<<level, k+1))
    for place in 1:(1<<level)
        for f_number in 0:k
            coeffs[place,f_number+1] = hquadrature(x->(basis(level,place,f_number,x)*f(x)),
            (place-1)/(1<<level), place/(1<<level), abstol=1.0e-10)[1]
        end
    end
    return coeffs
end

function get_vcoeffs(level::Int, k::Int, f::Function)
    vcoeffs = Array(Float64,(1<<level)*(k+1))
    i=1
    for place in 1:(1<<level)
        for f_number in 0:k
            vcoeffs[i] = hquadrature(x->(basis(level,place,f_number,x)*f(x)),
                            (place-1)/(1<<level), place/(1<<level), abstol=1.0e-10)[1]
            i+=1
        end
    end
    return vcoeffs
end



function reconstruct_coeffs(coeffs::Array{Float64,2}, x::Real)
    value = 0.0;
    level = Int(round(log2(length(coeffs[:,1]))))
    k = length(coeffs[1,:])-1
    for place in 1:(1<<level)
        for f_number in 0:k
            value += coeffs[place,f_number+1]*basis(level,place,f_number,x)
        end
    end
    return value
end

function reconstruct_vcoeffs(level::Int, k::Int, vcoeffs::Array{Float64,1}, x::Real)
    value = 0.0
    i=1
    for place in 1:(1<<level)
        for f_number in 0:k
            value += vcoeffs[i]*basis(level,place,f_number,x)
            i+=1
        end
    end
    return value
end

function legvDv(level, place1, f_number1, place2, f_number2)
    if place1 == place2
        return hquadrature(x->(basis(level,place1,f_number1,x)*dbasis(level,place2,f_number2,x)),
            (place1-1)/(1<<level), place1/(1<<level), abstol=1.0e-10)[1]
    end
    return 0.0
end

function D_matrix(level::Int, k::Int)
    mat = spzeros((1<<level)*(k+1), (1<<level)*(k+1))
    i=1
    for place1 in 1:(1<<level)
        for f_number1 in 0:k
            j=1
            for place2 in 1:(1<<level)
                for f_number2 in 0:k
                    ans = legvDv(level, place1, f_number1, place2, f_number2)
                    if abs(ans)<1.0e-15
                        j+=1
                        continue
                    end
                    mat[i,j]=ans
                    j+=1
                end
            end
            i+=1
        end
    end
    return mat
end


function legfDeltav(level, f::Function, place2, f_number2)
    point1 = (place2-1)/(1<<level)
    point2 = (place2)/(1<<level)
    
    mean1 = 0.5* (f(point1-5.0e-16)+f(point1+5.0e-16))
    mean2 = 0.5* (f(point2-5.0e-16)+f(point2+5.0e-16))
    
    val1 = basis(level,place2,f_number2,point1+5.0e-16)
    val2 = basis(level,place2,f_number2,point2-5.0e-16)
    
    return mean2 * val2 - mean1 * val1
end



function legvDeltav(level, place1, f_number1, place2, f_number2)
    point1 = (place2-1)/(1<<level)
    point2 = (place2)/(1<<level)
    
    mean1 = 0.5* (basis(level,place1,f_number1,point1-5.0e-16)+basis(level,place1,f_number1,point1+5.0e-16))
    mean2 = 0.5* (basis(level,place1,f_number1,point2-5.0e-16)+basis(level,place1,f_number1,point2+5.0e-16))
    
    if place2==(1<<level)
        mean2 = basis(level,place1,f_number1,point2-5.0e-16)
    end
    if place2==1
        mean1 = basis(level,place1,f_number1,point1+5.0e-16)
    end
    
    val1 = basis(level,place2,f_number2,point1+5.0e-16)
    val2 = basis(level,place2,f_number2,point2-5.0e-16)
    
    return mean2 * val2 - mean1 * val1
end


function Delta_matrix(level::Int, k::Int)
    mat = spzeros((1<<level)*(k+1), (1<<level)*(k+1))
    i=1
    for place1 in 1:(1<<level)
        for f_number1 in 0:k
            j=1
            for place2 in 1:(1<<level)
                for f_number2 in 0:k
                    ans = legvDeltav(level, place1, f_number1, place2, f_number2)
                    if abs(ans)<1.0e-15
                        j+=1
                        continue
                    end
                    mat[i,j]=ans
                    j+=1
                end
            end
            i+=1
        end
    end
    return mat
end