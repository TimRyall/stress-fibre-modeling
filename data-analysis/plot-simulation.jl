using Plots, RCall

function plot_sim(x,y,t,p,write::Bool=false)
    padding = p.L
    X_span = (p.A-padding,p.B+padding) # Range for filament movement
    # Set up domain
    plot(xlims=X_span, ylims=(0,p.N+1), title="Time = $t") # For square: size=(800,800)
    # Display filaments
    a = x[p.filaments] .- p.L/2 # Start of filaments
    b = x[p.filaments] .+ p.L/2 # End of filaments
    line_colour = reshape([pol > 0 ? :blue : :red for pol in p.P], 1, length(p.P)) # Polarity of filaments
    plot!([a b]', [p.filaments p.filaments]', legend=false, lc=line_colour)
    # Display motors and their attachments
    for m in p.motors
        af = y[m-p.N]
        if isempty(af)
            scatter!([x[m]], [0], markercolor=:black, markerstroke=0, markersize=3)
        else
            scatter!([x[m]], [mean(af)], markercolor=:black, markerstroke=0, markersize=3)
            plot!([x[m];x[m]], [af[1];af[end]], legend=false, lc=:black)
        end
    end
    # Display focal adhesions
    p.N % 2 == 0 ? height = (p.N+1)/2 : height = p.N/2 # Ensure focal adhesions are not overlapping a filament
    plot!([x[p.focal_adhesions[1]]-p.l/2 x[p.focal_adhesions[2]]-p.l/2;x[p.focal_adhesions[1]]+p.l/2 x[p.focal_adhesions[2]]+p.l/2], [height height; height height], legend=false, lc=:black, linewidth=2)
    if !write
        plot!(show=true)
        sleep(0.1)
    end
end

function save_sim(p)
    frame(p.anim)
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
    plot(x_vec, y_vec, xlab=param, ylab="Contractile Force", main="LOESS Regression")
    lim <- par("usr")
    rect(lim[1]-1, 0, lim[2]+1, lim[4]+1, border = adjustcolor("gainsboro", alpha.f = 0.20), col = adjustcolor("gainsboro", alpha.f = 0.20))
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

function loess_plot_2(x_vec::Vector{Float64},y_vec::Vector{Float64}, param::String, filename::String)
    gr(); plot(); # Load GR plotting back-end and clear previous plots
    @rput y_vec x_vec param filename
    R"""
    library(spatialEco)
    loess_model = loess.ci(y_vec, x_vec, p=0.99, plot = FALSE, span=0.6)
    pdf(file = filename)
    plot(x_vec, y_vec, xlab=param, ylab="Contractile Force", main="LOESS Regression")
    lim <- par("usr")
    rect(lim[1]-1, 0, lim[2]+1, lim[4]+1, border = adjustcolor("gainsboro", alpha.f = 0.20), col = adjustcolor("gainsboro", alpha.f = 0.20))
    rect(lim[1]-1, lim[3]-1, 0.25, lim[4]+1, border = adjustcolor("palegreen", alpha.f = 0.20), col = adjustcolor("palegreen", alpha.f = 0.20))
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

function loess_plot_3(x_vec::Vector{Float64},y_vec::Vector{Float64}, param::String, filename::String)
    gr(); plot(); # Load GR plotting back-end and clear previous plots
    @rput y_vec x_vec param filename
    R"""
    library(spatialEco)
    loess_model = loess.ci(y_vec, x_vec, p=0.99, plot = FALSE, span=0.6)
    pdf(file = filename)
    plot(x_vec, y_vec, xlab=param, ylab="Contractile Force", main="Model fit to Simuation Data")
    lim <- par("usr")
    rect(lim[1]-1, 0, lim[2]+1, lim[4]+1, border = adjustcolor("gainsboro", alpha.f = 0.20), col = adjustcolor("gainsboro", alpha.f = 0.20))
    rect(lim[1]-1, lim[3]-1, 0.25, lim[4]+1, border = adjustcolor("palegreen", alpha.f = 0.20), col = adjustcolor("palegreen", alpha.f = 0.20))
    x = seq(0,1,0.01)
    fit <- nls(y_vec ~ (-(2/3)*x_vec + (1/6)) * ((1-x_vec))^2 * A ,start=list(A=2))
    result <- summary(fit)$coefficients
    K = result[1]
    y = (-(2/3)*x + (1/6)) * ((1-x))^2 * K
    #curve(fit, col="black", add = TRUE)
    lines(x,y, col="black")
    lines(x_vec, loess_model$uci, col="red", lty=2)
    lines(x_vec, loess_model$lci, col="red", lty=2)
    polygon(x = c(x_vec, rev(x_vec)),
            y = c(loess_model$lci,
                rev(loess_model$uci)),
            col =  adjustcolor("red", alpha.f = 0.10), border = NA)
    x_mod = seq(0,1,0.025)
    y_mod = (-(2/3)*x_mod + (1/6)) * ((1-x_mod))^2 * 2
    r_2_value = cor(y_mod, y_vec)^2
    mse = mean((y_vec - y_mod)^2)
    rp = vector('expression',2)
    rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
        list(MYVALUE = format(r_2_value,dig= 4)))[2]
    rp[2] = substitute(expression(italic(MSE) == MYOTHERVALUE), 
        list(MYOTHERVALUE = format(mse, digits = 2)))[2]
    legend('topright', legend = rp, bty = 'n')
    dev.off()
    print(r_2_value)
    print(mse)
    print("test")  
    print("K:")
    print(K) 
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