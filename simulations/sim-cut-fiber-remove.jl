include("src/main.jl")

# Visualise simulation with reference parameters
main_cut1(Parameters(δ=1,β=0.01),true,false,"",25,100);

# Plot contractile force of stress fibres with early motor drop-offs and motor turnover
conf = []
for i in 1:1
    cf = main_cut1(Parameters(δ=1,β=0.01),false,false,"",100)
    push!(conf,cf)
end
plot(1:250,conf,xlabel="Time",ylabel="Contractile Force",legend=false)
savefig("example-figures/ref-cf.pdf")