using Random, Distributions

using VortexDynamics
using ComponentArrays

name = "vortexpair-groups"

σ = 0.1

# 1. Define initial configuration of vortices.
N1 = 3 
x1 = rand( Normal( 0.0, σ ), N1 )
y1 = rand( Normal( 0.5, σ ), N1 )
γ1 = rand( Normal( 1.0, σ ), N1 ) ./ N1 # divide by N1 to maintain ~ constant total circulation

N2 = 6
x2 = rand( Normal( 0.0, σ ), N2 )
y2 = rand( Normal( 1.0, σ ), N2 )
γ2 = rand( Normal( -0.5, σ ), N2 ) ./ N2 # divide by N2 to maintain ~ constant total circulation

u0 = ComponentArray(
x = [x1; x2],
y = [y1; y2]
)
γ = [γ1; γ2]
nVortices = N1 + N2

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
ax = Makie.heatmap(grid_x, grid_y, Ωs[1], colormap=:balance )

for k = 1:length(γ)
    lines!( getindex.( getindex.(us, :x), k ),
            getindex.( getindex.(us, :y), k ) )
end
ax
