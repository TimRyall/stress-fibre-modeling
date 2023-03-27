using Plots, RCall

function plot_sim(x,y,t,p,write::Bool=false)
    padding = p.Lmin
    X_span = (p.A-padding,p.B+padding) # Range for filament movement
    # Set up domain
    plot_curr = plot(xlims=X_span, ylims=(0,findmax(p.filament_vertical_pos)[1]+1), title="Time = $t"); # For square: size=(800,800)
    # Display filaments
    a = x[p.filaments] - p.L./2 # Start of filaments
    b = x[p.filaments] + p.L./2 # End of filaments
    line_colour = reshape([pol > 0 ? :blue : :red for pol in p.P], 1, length(p.P)) # Polarity of filaments
    plot!([a b]', [p.filament_vertical_pos p.filament_vertical_pos]', legend=false, lc=line_colour)
    # Display motors and their attachments
    for m in p.motors
        af = y[m-p.N]
        if isempty(af)
            scatter!([x[m]], [0], markercolor=:black, markerstroke=0, markersize=3)
        else
            # Change from filament index to vertical position
            af_new = zeros(2)
            af_new[1] = p.filament_vertical_pos[af[1]]
            af_new[2] = p.filament_vertical_pos[af[2]]
            scatter!([x[m]], [mean(af_new)], markercolor=:black, markerstroke=0, markersize=3)
            plot!([x[m];x[m]], [af_new[1];af_new[end]], legend=false, lc=:black)
        end
    end
    # Display focal adhesions
    p.N % 2 == 0 ? height = (findmax(p.filament_vertical_pos)[1]+1)/2 : height = findmax(p.filament_vertical_pos)[1]/2 # Ensure focal adhesions are not overlapping a filament
    plot!([x[p.focal_adhesions[1]]-p.l/2 x[p.focal_adhesions[2]]-p.l/2;x[p.focal_adhesions[1]]+p.l/2 x[p.focal_adhesions[2]]+p.l/2], [height height; height height], legend=false, lc=:black, linewidth=2)
    if !write
        plot!(show=true)
        sleep(0.1)
    end
    return plot_curr
end

function save_sim(p)
    frame(p.anim)
    return p
end

function save_sim_cf(p)
    frame(p.anim_cf)
    return p
end

# Constructing LOESS Regression plot using RCall package, and LOESS Regression Package in R
function loess_plot(x_vec::Vector{Float64},y_vec::Vector{Float64}, param::String, filename::String)
    gr(); plot(); # Load GR plotting back-end and clear previous plots
    @rput y_vec x_vec param filename
    R"""
    library(spatialEco)
    loess_model = loess.ci(y_vec, x_vec, p=0.99, plot = FALSE, span=0.6)
    pdf(file = filename)
    plot(x_vec, y_vec, xlab="Time Step", ylab="Distance", main="LOESS Regression")
    lines(x_vec, loess_model$loess, col="red")
    lines(x_vec, loess_model$uci, col="red", lty=2)
    lines(x_vec, loess_model$lci, col="red", lty=2)
    polygon(x = c(x_vec, rev(x_vec)),
            y = c(loess_model$lci,
                rev(loess_model$uci)),
            col =  adjustcolor("red", alpha.f = 0.10), border = NA)
    dev.off()
    """
end;

# Constructing LOESS Regression log plot using RCall package, and LOESS Regression Package in R
function loess_plot_log(x_vec::Vector{Float64},y_vec::Vector{Float64}, param::String, filename::String)
    gr(); plot(); # Load GR plotting back-end and clear previous plots
    @rput y_vec x_vec param filename
    R"""
    library(spatialEco)
    loess_model = loess.ci(y_vec, x_vec, p=0.99, plot = FALSE, span=0.6)
    pdf(file = filename)
    plot(x_vec, y_vec, log="x", xlab=param, ylab="Contractile Force", main="LOESS Regression")
    lines(x_vec, loess_model$loess, col="red")
    lines(x_vec, loess_model$uci, col="red", lty=2)
    lines(x_vec, loess_model$lci, col="red", lty=2)
    polygon(x = c(x_vec, rev(x_vec)),
            y = c(loess_model$lci,
                rev(loess_model$uci)),
            col =  adjustcolor("red", alpha.f = 0.10), border = NA)
    dev.off()
    """
end;
