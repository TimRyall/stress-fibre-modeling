using Plots, Optim, DataFrames, CSV
using Plots.PlotMeasures
using Distributions: Poisson, sample
using StatsBase: mean
include("model-parameters.jl")
include("model.jl")
include("plot-simulation.jl")

# Orginal default simulation 
# Simulate a stress fibre with specified parameters, choosing whether to visualise the sim and/or save it as a gif (plotting is required to write)
function main(params::Parameters,PLOTSIM::Bool=false,WRITESIM::Bool=false,filename::String="",set_fps::Int64=100)
    # Set up simulation
    X, Y, contractile_force, params = reset_model(params)

    # Evolve simulation
    for t in 1:params.T
        println(t)
        # Store current simulation status
        contractile_force[t] = params.k * (X[t][params.focal_adhesions[2]] - X[t][params.focal_adhesions[1]] - (params.B - params.A)) # Calculate current contractile force between focal adhesions
        
        #=
        if PLOTSIM 
            plotsimu = plot_sim(X[t],Y,t,params,WRITESIM);
            avg_cf = round((sum(contractile_force[1:t])/t); digits = 2)
            if contractile_force[t] <= 0
                plotcf = plot(1:t,contractile_force[1:t],xlabel="Time",ylabel="Contractile Force",legend=false,c=:green, title="Average Contractile Force: $(avg_cf)");
            else
                plotcf = plot(1:t,contractile_force[1:t],xlabel="Time",ylabel="Contractile Force",legend=false,c=:red, title="Average Contractile Force: $(avg_cf)");
            end
            plot(plotcf, plotsimu, layout = (1, 2), legend = false)
            plot!(size=(1200,400))
            plot!(margin = 10mm)
        WRITESIM && (params = save_sim(params))
        end
        =#

        avg_cf = round((sum(contractile_force[1:t])/t); digits = 2)

        PLOTSIM && plot_sim(X[t],Y,t,params,WRITESIM)
        WRITESIM && (params = save_sim(params))
        # Store current cf graph
        PLOTSIM && contractile_force[t] <= 0  && plot(1:t,contractile_force[1:t],xlabel="Time",ylabel="Contractile Force",legend=false,c=:green,title="Average Contractile Force: $(avg_cf)")
        PLOTSIM && contractile_force[t] > 0 && plot(1:t,contractile_force[1:t],xlabel="Time",ylabel="Contractile Force",legend=false,c=:red,title="Average Contractile Force: $(avg_cf)")
        WRITESIM && (params = save_sim_cf(params))

        # Perform minimisation to advance simulation to next time step
        od = OnceDifferentiable(x -> E(x, X[t], O(X[t][params.filaments],params.L), Oᵩ(X[t][params.filaments],X[t][params.focal_adhesions],params.L,params.l), Y, params), X[t]; autodiff=:forward)
        X[t+1] = Optim.minimizer(optimize(od, X[t], LBFGS()))

        # Perform random filament and motor turnover
        X[t+1] = turnover(X[t+1],params)
        # Update attached filaments to motors
        Y = update_af(X[t+1],Y,params)
    end

    # Output results
    print("fps:")
    print(set_fps)
    WRITESIM && gif(params.anim, "sim.gif", fps=set_fps)
    WRITESIM && gif(params.anim_cf, "cf.gif", fps=set_fps)
    return contractile_force[1:end-1]
end;


# Simulation where the fiber is cut down the middle deleting all filaments in the cut range
function main_cut1(params::Parameters,PLOTSIM::Bool=false,WRITESIM::Bool=false,filename::String="",set_fps::Int64=5,cuttimestep::Int64=100)
    # Set up simulation
    X, Y, contractile_force, params = reset_model(params)

    # Evolve simulation
    for t in 1:params.T
        # Store current simulation status
        contractile_force[t] = params.k * (X[t][params.focal_adhesions[2]] - X[t][params.focal_adhesions[1]] - (params.B - params.A)) # Calculate current contractile force between focal adhesions
        PLOTSIM && plot_sim(X[t],Y,t,params,WRITESIM)
        WRITESIM && (params = save_sim(params))

        # Perform minimisation to advance simulation to next time step
        od = OnceDifferentiable(x -> E(x, X[t], O(X[t][params.filaments],params.L), Oᵩ(X[t][params.filaments],X[t][params.focal_adhesions],params.L,params.l), Y, params), X[t]; autodiff=:forward)
        X[t+1] = Optim.minimizer(optimize(od, X[t], LBFGS()))

        # Perform random filament and motor turnover
        X[t+1] = turnover(X[t+1],params)
        # Update attached filaments to motors
        Y = update_af(X[t+1],Y,params)

        # Perform cut down the middle of the fiber if at correct cut time
        cut_position = (params.B + params.A)/2 # Cut position is in the middle of the fiber
        if t == cuttimestep-1
            X[t+1], params = cut_v1(X[t+1],params,cut_position)
        end
    end

    # Output results
    WRITESIM && gif(params.anim, "$(filename).gif", fps=set_fps)
    return contractile_force[1:end-1]
end;

# Simulation where the fiber is cut down the middle severing all filaments in the cut range
# i.e for each filament severed it will create two 'new' filaments from its halves
function main_cut2(params::Parameters,PLOTSIM::Bool=false,WRITESIM::Bool=false,filename::String="",set_fps::Int64=5,cuttimestep::Int64=100)
    # Set up simulation
    X, Y, contractile_force, params = reset_model(params)
    distance = Vector{Float64}(undef,params.T+1) # Distance between 2 severed halves
    distance2 = Vector{Float64}(undef,params.T+1) # Distance between 2 severed halves

    # Evolve simulation
    println("Time Step:")
    for t in 1:params.T
        println(t)
        # Store current simulation status
        contractile_force[t] = params.k * (X[t][params.focal_adhesions[2]] - X[t][params.focal_adhesions[1]] - (params.B - params.A)) # Calculate current contractile force between focal adhesions
        # Distance between 2 severed halves
        lh_total = 0.0
        lh_n = 0
        rh_total = 0.0
        rh_n = 0
        for fil in params.filaments
            if (params.A + params.B)/2 > X[t][fil] # Filament in left half
                lh_total += X[t][fil] 
                lh_n += 1
            else # Filament in right half
                rh_total += X[t][fil]
                rh_n += 1
            end
        end
        distance[t] = (rh_total/rh_n) - (lh_total/lh_n)  
        
        lh = params.A
        rh = params.B
        for fil in params.filaments
            if (params.A + params.B)/2 > X[t][fil] # Filament in left half
                lh = max(lh, X[t][fil] + params.L[fil]/2)
            else # Filament in right half
                rh = min(rh, X[t][fil] - params.L[fil]/2)
            end
        end
        distance2[t] = rh - lh

        PLOTSIM && plot_sim(X[t],Y,t,params,WRITESIM)
        WRITESIM && (params = save_sim(params))

        # Perform minimisation to advance simulation to next time step
        od = OnceDifferentiable(x -> E(x, X[t], O(X[t][params.filaments],params.L), Oᵩ(X[t][params.filaments],X[t][params.focal_adhesions],params.L,params.l), Y, params), X[t]; autodiff=:forward)
        X[t+1] = Optim.minimizer(optimize(od, X[t], LBFGS()))

        # Perform random filament and motor turnover
        X[t+1] = turnover(X[t+1],params)
        # Update attached filaments to motors
        Y = update_af(X[t+1],Y,params)

        # Perform cut down the middle of the fiber if at correct cut time
        cut_position = (params.B + params.A)/2 # Cut position is in the middle of the fiber
        cut_width = (params.Lmin)/2
        if t == cuttimestep-1
            X[t+1], params = cut_v2(X[t+1],params,cut_position,cut_width)
        end
    end

    # Output results
    WRITESIM && gif(params.anim, "$(filename).gif", fps=set_fps)
    return contractile_force[1:end-1], distance[1:end-1], distance2[1:end-1]
end;