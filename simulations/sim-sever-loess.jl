include("src/main.jl")
num_time_steps = 3000
p = Parameters(δ=1.0,β=0.0,T=num_time_steps,Lmin=1.0,Lsd=0.5, N=100, M=50)
# Visualise simulation with reference parameters
# main_cut2(p,true,true,"testing3",25,100);

# Plot contractile force of stress fibres with early motor drop-offs and motor turnover
# & Plot avgerage distance bewteen the two severed halves of the fiber

cf, dt, dt2 = main_cut2(p,false,false,"",25,100)
distance_vec = dt2[100:100:num_time_steps]
distance_vec_avg = dt[100:100:num_time_steps]
time_vec = collect(100:100:num_time_steps)
time_vec2 = convert(Array{Float64,1}, time_vec)

loess_plot(time_vec2,distance_vec,"Test","example-figures/endpoints-loess-end.pdf")
loess_plot(time_vec2,distance_vec_avg,"Test","example-figures/endpoints-loess-avg.pdf")
