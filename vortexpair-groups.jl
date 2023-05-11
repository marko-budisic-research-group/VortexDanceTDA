using Random, Distributions

using VortexDynamics
using ComponentArrays

name = "vortexpair-groups"

result_dir="/Users/nipuni/Desktop/Research/VortexDanceTDA/Results-Nipuni/vortex-pair-groups/"
pixel_size=[16, 101]
result_folder = ["t100-16x16","t100-101x101"]
noICs = 30


σ = 0.1

function vortex_group( N, center, γ, σ_xy, σ_γ)
    x = rand( Normal( center[1], σ_xy ), N )
    y = rand( Normal( center[2], σ_xy ), N )
    γ = rand( Normal( γ, σ_γ ), N ) ./ N # divide by N1 to maintain ~ constant total circulation
    return x, y, γ
end

########################################################
function randomSetSum(target_sum, num_elements)
    # Generate random numbers
    numbers = rand(num_elements)
    
    # Calculate the sum of the numbers
    sum_numbers = sum(numbers)
    
    # Scale the numbers to match the target sum
    scaled_numbers = (target_sum / sum_numbers) * numbers
    
    # Round the scaled numbers to integers
    rounded_numbers = round.(Int, scaled_numbers)
    
    # Calculate the adjustment needed to match the target sum
    adjustment = target_sum - sum(rounded_numbers)
    
    # Adjust one of the numbers randomly to match the target sum
    index = rand(1:num_elements)
    rounded_numbers[index] += adjustment
    
    return rounded_numbers
end

########################################################


 for i = 1:2
    for j = 1:noICs
        # 1. Define initial configuration of vortices.
        c = round.((2 .* rand(Float64, 4) .- 1)./ 1.5, digits=1) 
        groups = 2
        noVortices = randomSetSum(11, groups)

        x1, y1, γ1 = vortex_group( noVortices[1], (c[1], c[2]), 1.0, 0.1, 0.1 )
        x2, y2, γ2 = vortex_group( noVortices[2], (c[3], c[4]), -0.5, 0.1, 0.1 )

        u0 = ComponentArray(
        x = [x1; x2],
        y = [y1; y2]
        )
        γ = [γ1; γ2]

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
        grid_x = LinRange(-1.2,1.2,pixel_size[i])
        grid_y = LinRange(-1.2,1.2,pixel_size[i])
        vxs, vys, Ωs = getfields( us, ts; px=grid_x, py=grid_y, γ, ν=1e-3 );

        # 4. Store to numpy files
        using NPZ
        using Dates

        uniquelabel = string(today())*"-"*randstring(2)
        npzwrite(result_dir *"vortex-tracks/"*result_folder[i]*"/$name-tracks-" * uniquelabel *".npz", 
            ux = ux,
            uy = uy,
            t = ts,
            gamma = γ)


        npzwrite(result_dir *"vortex-fields/"*result_folder[i]*"/$name-fields-"* uniquelabel *".npz", 
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

        if i == 2
            record(ax, result_dir *"$name-"* uniquelabel *".mp4", 1:length(Ωs); framerate = 60) do i
                Makie.heatmap!(grid_x, grid_y, Ωs[i], colormap=:balance )
            end;
        end

    end
end