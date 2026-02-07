using LinearAlgebra

domain = range(0,2*π, step=0.01)
#print(parentmodule(plot))
plot(domain, [cos.(domain), sin.(domain)], title="Sin & Cos", colors=["red", "blue"], label=["Cos" "Sin"])