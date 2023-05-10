using Random, Distributions

using VortexDynamics
using ComponentArrays

name = "vortex-arc"


function vortex_group( N, center, γ, σ_xy, σ_γ)
    x = rand( Normal( center[1], σ_xy ), N )
    y = rand( Normal( center[2], σ_xy ), N )
    γ = rand( Normal( γ, σ_γ ), N ) ./ N # divide by N1 to maintain ~ constant total circulation
    return x, y, γ
end

N  = 10
angles = collect(LinRange(0,π/2,10))
u0 = ComponentArray(
x = cos.(angles),
y = sin.(angles)
)
γ = ones(N) ./ N

# 2. Simulate the differential equation
using OrdinaryDiffEq

t0 = 0.
tf = 20.
Δ = 1e-1;

diffeq = ODEProblem( vortex_biot_savart!, u0, (t0,tf),γ )
sol = solve( diffeq, Tsit5(); saveat=Δ )

# extract just a subset of timesteps (not really necessary)
subselect = 1:100
us = @view sol.u[subselect]
ts = @view sol.t[subselect]



# extract horizontal and vertical positions of vortices for each snapshot and merge (hcat) them into a matrix
ux = hcat(getindex.(us,:x)...)
uy = hcat(getindex.(us,:y)...)

# 3. Evaluate velocity and vorticity fields on a grid

# get velocity and vorticity fields
grid_x = LinRange(-1.2,1.2,101)
grid_y = LinRange(-1.2,1.2,101)
vxs, vys, Ωs = getfields( us, ts; px=grid_x, py=grid_y, γ, ν=1e-3 );

# 4. Store to numpy files
using NPZ
using Dates

uniquelabel = string(today())*"-"*randstring(2)
npzwrite("results/$name-tracks-" * uniquelabel *".npz", 
    ux = ux,
    uy = uy,
    t = ts,
    gamma = γ)


npzwrite("results/$name-fields-"* uniquelabel *".npz", 
    grid_x = collect(grid_x),
    grid_y = collect(grid_y),
    t = ts,
    gamma = γ,
    Vx = cat(vxs...; dims=3),
    Vy = cat(vys...; dims=3),
    Omega = cat(Ωs...; dims=3))


# visualize first snapshot of vorticity and the vortex track overlaid
using CairoMakie
using Makie
ax = Makie.heatmap(grid_x, grid_y, Ωs[end], colormap=:balance )

for k = 1:length(γ)
    lines!( getindex.( getindex.(us, :x), k ),
            getindex.( getindex.(us, :y), k ) )
end
ax
