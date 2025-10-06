using Plots
using LaTeXStrings
using Polynomials
using PrettyTables

"""
    simple_iteration( g, x1; N=100, tol=1e-10  ) 

Computes the simple iteration ``x_{n+1} = g(x_n)``. Inputs ``g`` (a function) and ``x1`` (an initial starting point). Optional Inputs: ``N`` (max number of iterations) and ``tol`` (tolerance with the stopping criteria ``|f(x_n)| < tol``). Outputs ``x``, the vector (x1, x2, ...) 
"""
function simple_iteration( g, x1; N=100, tol=1e-10 )
    println( "x[n+1] = g(x[n]) with x_1 = ", x1 )
    x = [ x1 ]
    for n in 2:N
        push!( x, g(x[n-1]) )
        if (abs(g(x[end]) - x[end]) < tol)
            break
        elseif (x[end] == Inf)
            @warn "simple iteration diverges to Inf";
            break
        elseif (x[end] == -Inf)
            @warn "simple iteration diverges to -Inf";
            break
        end 
    end
    return orderOfConvergenceF(x, x->(g(x) - x) )
end

"""
Computes the iteration ``x_{n+1} = x_n -  f(x_n) / f'(x_n)``.

## Input
`f`   - function \n
`df`  - derivative of ``f`` \n
`x1`  - initial starting point

## Optional Inputs:    
`N`   - max number of iterations \n
`tol` - tolerance 

## Output
`x`   - vector ``(x_1, x_2, ...)`` 
"""
function Newton( f, f_prime, x1; N=100, tol=1e-10)
    x = [x1]
    r = 0;
    for n in 2:N
        push!( x, x[n-1] - f(x[n-1])/f_prime(x[n-1]) )
        r = abs(f(x[end]));

        if (r < tol)
            break;
        end
    end
    return orderOfConvergenceF( x, f )
end

"""
    orderOfConvergence( x, ξ; α=0 )

Computes the approximate order of convergence of the method

## Input:              
`x`   - vector of iterates \n
`ξ`   - zero that you are trying to approximate 

## Optional Input: 
`α`   - order of convergence

## Output: 
Returns a table containing the errors ``| x_n - ξ |``, order of convergence ``α ≈ ( \\log|x_{n+1}-ξ| ) / ( \\log|x_n-ξ| )``, and asymptotic error constant ``μ``,
"""
function orderOfConvergence( x, ξ; α=0 )
    err = @. abs( x - ξ )
    logerr = @. log10( err )
    ratios = [NaN; [logerr[i+1] / logerr[i] for i in 1:length(logerr)-1]]
    if (α == 0) 
        α = ratios[end]
    end
    αr = round(α)
    mu = [NaN; [err[i+1] / err[i]^α for i in 1:length(err)-1]]
    pretty_table( 
        [1:length(x) x err ratios mu];
        column_labels = [
           ["iteration", "sequence", "abs. error", "ratio", "μ (α = $αr)"],
           [ "n", "x[n]", "e[n]=|x[n]-ξ|", "log e[n]/log e[n-1]", "e[n]/e[n-1]^α"   ]
       ]
    )
end 

"""
    orderOfConvergenceF( x, f )

Computes the approximate order of convergence of the method

## Input:              
`x`   - vector of iterates \n
`f`   - function you are trying to find roots of 

## Optional Input: 
`α`   - order of convergence

## Output: 
Returns a table containing the residuals ``| f(x_n) |`` and computational order of convergence ``α ≈ ( \\log|f(x_{n+1})| ) / ( \\log|f(x_n)| )``
"""
function orderOfConvergenceF( x, f )
    err = @. abs( f( x ) )
    logerr = @. log10( err )
    ratios = [NaN; [logerr[i+1] / logerr[i] for i in 1:length(logerr)-1]]
    pretty_table( 
        [1:length(x) x err ratios];
        column_labels = [
           ["iteration", "sequence", "residual", "approximate α"],
           [ "n", "x[n]", "r[n]=|f(x[n])|", "log r[n]/log r[n-1]"]
       ]
    )
end 

"""
    ChebyshevNodes( n; type=2) 

Returns the ``n+1`` Chebyshev nodes ``X = \\{x_0, ..., x_n\\}`` of type type (default being type II) 
"""
function ChebyshevNodes( n; type=2) 
    if (type == 1)
        return @. cos( (2*(0:n)+1)*pi/(2*(n+1)) )
    else
        return @. cos( (pi/n)*(0:n) )
    end
end

"""
    ChebyshevInterpolant( f, n ) 

Returns the polynomial of degree ``n`` interpolating ``f`` at Chebyshev nodes using the Barycentric formula
"""
function ChebyshevInterpolant( f, n ) 
    x = ChebyshevNodes(n)
    y = f.(x)

    λ = zeros(n+1)
    λ[1] = 1/2
    λ[n+1] = (-1)^n/2
    for j in 1:(n-1)
        λ[j+1] = (-1)^j
    end

    return t -> sum( @. λ * y / (t - x) ) / sum( @. λ / (t - x) )
end


"""
    Runge( x; a=25 ) 

Returns the Runge function
"""
function Runge( x; a=25 )
    return 1 / ( 1 + a*x^2 )
end










println("✓ file included! \n
using: Plots, LaTeXStrings, Polynomials, PrettyTables \n
Functions included: 
    simple_iteration, 
    Newton, 
    orderOfConvergence, 
    ChebyshevNodes \n
Use @doc <<function>> for help")