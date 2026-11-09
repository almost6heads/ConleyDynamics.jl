export ouroboros_vf
export vdp_vf, J_vdp_vf
export vdp_vf_trev, J_vdp_vf_trev
export vdp_vf_ode!
export vf_curv

# Ouroboros Counter Example Vector Field
function ouroboros_vf(x::Vector{Float64})
    
    if norm(x - [2.0, 3.5]) < 0.1
        return [0.0, 1.0]
    elseif norm(x - [2.25, 3.25]) < 0.1
        return [1.0, 1.0]
    elseif norm(x - [2.5, 3.0]) < 0.1
        return [1.0, 0.0]
    elseif norm(x - [2.75, 2.75]) < 0.1
        return [1.0, 1.0]
    elseif norm(x - [2.75, 2.25]) < 0.1
        return [1.0,-1.0]
    elseif norm(x - [2.5, 2.0]) < 0.1
        return [1.0, 0.0]
    elseif norm(x - [2.25, 1.75]) < 0.1
        return [1.0,-1.0]
    elseif norm(x - [2.0, 1.5]) < 0.1
        return [0.0,-1.0]
    elseif norm(x - [1.75, 1.75]) < 0.1
        return [1.0, 1.0]
    elseif norm(x - [1.5, 2.0]) < 0.1
        return [1.0, 0.0]
    elseif norm(x - [1.25, 2.25]) < 0.1
        return [1.0, 1.0]
    elseif norm(x - [1.25, 2.75]) < 0.1
        return [1.0,-1.0]
    elseif norm(x - [0.25, 2.25]) < 0.1
        return [1.0, 1.0]
    elseif norm(x - [0.75, 2.25]) < 0.1
        return [1.0,-1.0]
    else
        return [2.1, 2.5] - x
    end
end

# Van Der Pol Vector Field
function vdp_vf(v::AbstractArray)
    x, y = v
    mu = 1.0

    xd = y
    yd = mu*(1 - x^2)*y - x

    return [xd, yd]
end

# Jacobian of Van Der Pol Vector Field
function J_vdp_vf(v::AbstractArray)
    x, y = v
    mu = 1

    Jxx = 0
    Jxy = 1
    Jyx = -1 - mu*2*y*x
    Jyy = mu*(1 - x^2)

    return [Jxx, Jxy, Jyx, Jyy]
end

# Time-reversed Van Der Pol Vector Field
function vdp_vf_trev(v::AbstractArray)
    x, y = v
    mu = 1.0

    xd = -y
    yd = -mu*(1 - x^2)*y + x

    return [xd, yd]
end

# Jacobian of Time-reversed Van Der Pol Vector Field
function J_vdp_vf_trev(v::AbstractArray)
    x, y = v
    mu = 1

    Jxx = 0
    Jxy = -1
    Jyx = 1 + mu*2*y*x
    Jyy = mu*(x^2 - 1)

    return [Jxx, Jxy, Jyx, Jyy]
end

# Calculate 2D vector field curvature
function vf_curv(v::AbstractArray, vf::Function, Jvf::Function)
    
    fx, fy = vf(v)
    Jxx, Jxy, Jyx, Jyy = Jvf(v)

    K_num = Jyx*fx^2 - Jxy*fy^2 + fx*fy*(Jyy - Jxx)
    K_den = sqrt(fx^2 + fy^2)^3

    return abs(K_num/K_den)
end

# Van Der Pol vector field for DifferentialEquations ODE solver
function vdp_vf_ode!(du,u,p,t)
    x, y = u
    mu = 1.0

    du[1] = y
    du[2] = mu*(1 - x^2)*y - x
end