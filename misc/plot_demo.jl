using Revise, Plots
using LaTeXStrings
pgfplotsx()
x = -2π:0.01:2π
y = sin.(x)
plot(x, y, framestyle=:origin, draw_arrow=true)
plot!(grid=false)
plot!(axis=true)
xlabel!(L"x")
ylabel!(L"y")
xticks!([0])
yticks!([0])