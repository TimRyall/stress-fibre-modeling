include("src/main.jl")

# Visualise simulation with reference parameters
# main(Parameters(δ=1,β=0.01,Lmin=2,Lsd=0),true,true,"testing2",25);

# Plot contractile force of stress fibres with early motor drop-offs and motor turnover
p = Parameters(δ=0,β=0.0,T=1000,Lmin=2,Lsd=0,N=30,M=15)

conf = []
for i in 1:1
    cf = main(p,true,true,"testing2", 50)
    push!(conf,cf)
end
plot(1:(p.T),conf,xlabel="Time",ylabel="Contractile Force",legend=false)
savefig("example-figures/ref-cf.pdf")