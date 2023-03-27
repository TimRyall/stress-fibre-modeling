include("src/main.jl")
num_time_steps = 250
p = Parameters(δ=0,β=0.0,T=num_time_steps,Lmin=2,Lsd=0, N=40, M=20)
# Visualise simulation with reference parameters
main_cut2(p,true,true,"testing3",25,100);

# Plot contractile force of stress fibres with early motor drop-offs and motor turnover
# & Plot avgerage distance bewteen the two severed halves of the fiber
conf = []
dist = []
dist2 = []
print("Graphs produced... ")
for i in 1:3
    cf, dt, dt2 = main_cut2(p,false,false,"",25,100)
    push!(conf,cf)
    push!(dist,dt[100:end])
    push!(dist2,dt2[100:end]) 
    print(i)
    print(" ")
end
plot(1:num_time_steps,conf,xlabel="Time",ylabel="Contractile Force",legend=false)
savefig("example-figures/ref-cf.pdf")
plot(100:num_time_steps,dist,xlabel="Time",ylabel="Distance",legend=false, title="Avg distance between severed fiber half filaments")
savefig("example-figures/ref-recoil-average.pdf")
plot(100:num_time_steps,dist2,xlabel="Time",ylabel="Distance",legend=false,title="Distance between severed fiber half endpoints")
savefig("example-figures/ref-recoil-endpoints.pdf")

loess_plot(param_vec,mean_cf_vec,"Motor Drop-off Distance","example-figures/endpoints-loess.pdf")


plot(100:num_time_steps,sum(dist)./length(dist),xlabel="Time",ylabel="Distance",legend=false, title="Mean distance over the trials")
savefig("example-figures/ref-recoil-average-mean.pdf")
plot(100:num_time_steps,sum(dist2)./length(dist2),xlabel="Time",ylabel="Distance",legend=false,title="Mean distance over the trials")
savefig("example-figures/ref-recoil-endpoints-mean.pdf")

