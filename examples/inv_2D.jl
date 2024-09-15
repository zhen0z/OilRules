## A simple 2D example for permeability inversion

using DrWatson
@quickactivate "OilRules-example"

using OilRules
using LinearAlgebra
using PyPlot
using Flux
using LineSearches

sim_name = "2D-K-inv"
exp_name = "channel"

mkpath(datadir())
mkpath(plotsdir())

## grid size
n = (30, 1, 15)
d = (30.0, 30.0, 30.0)

## permeability
K0 = 20 * md * ones(n)
K = deepcopy(K0)
K[:,:,8:10] .*= 6.0

ϕ = 0.25
model0 = jutulModel(n, d, ϕ, K1to3(K0))
model = jutulModel(n, d, ϕ, K1to3(K))

## simulation time steppings
tstep = 40 * ones(50)
tot_time = sum(tstep)

## injection & production
inj_loc = (3, 1, 9) .* d
prod_loc = (28, 1, 9) .* d
irate = 5e-3
q = jutulForce(irate, [inj_loc, prod_loc])

## set up modeling operator
S = jutulModeling(model, tstep)

## simulation
mesh = CartesianMesh(model)
T(x) = log.(KtoTrans(mesh, K1to3(x)))

K_init = deepcopy(K0)

@time state0 = S(T(K0), q)
@time state = S(T(K), q)

state_init = deepcopy(state0)
f(K) = .5 * norm(S(T(K),q)[1:length(tstep)*prod(n)]-state[1:length(tstep)*prod(n)])^2
ls = BackTracking(order=3, iterations=10)

lower, upper = 0.9*minimum(K), 1.1*maximum(K)
prj(x) = max.(min.(x,upper),lower)
# Main loop
niterations = 50
fhistory = zeros(niterations)

for j=1:niterations

    fval = f(K0)
    g = gradient(()->f(K0), Flux.params(K0))[K0]
    p = -g/norm(g, Inf) * norm(K_init, Inf)
    
    println("Inversion iteration no: ",j,"; function value: ",fval)
    fhistory[j] = fval

    # linesearch
    function f_(α)
        misfit = f(prj(K0 .+ α * p))
        @show α, misfit
        return misfit
    end

    step, fval = ls(f_, 1.0, fval, dot(g, p))

    # Update model and bound projection
    global K0 = prj(K0 .+ step .* p)
end

## plotting
state0 = S(T(K0), q)

fig_name = @strdict n d ϕ tstep irate niterations lower upper inj_loc prod_loc

fig=figure(figsize=(20,12));
subplot(1,3,1);
imshow(K_init[:,1,:]', vmin=minimum(K), vmax=maximum(K)); colorbar(); title("initial permeability")
subplot(1,3,2);
imshow(K0[:,1,:]', vmin=minimum(K), vmax=maximum(K)); colorbar(); title("inverted permeability")
subplot(1,3,3);
imshow(K[:,1,:]', vmin=minimum(K), vmax=maximum(K)); colorbar(); title("true permeability")
tight_layout()
safesave(joinpath(plotsdir(sim_name, exp_name), savename(fig_name; digits=6)*"_K.png"), fig);
close(fig)

fig=figure(figsize=(20,12));
subplot(1,3,1);
imshow(reshape(state_init.states[end].Saturations, n[1], n[end])', vmin=0, vmax=1); colorbar(); title("initial saturation")
subplot(1,3,2);
imshow(reshape(state0.states[end].Saturations, n[1], n[end])', vmin=0, vmax=1); colorbar(); title("inverted saturation")
subplot(1,3,3);
imshow(reshape(state.states[end].Saturations, n[1], n[end])', vmin=0, vmax=1); colorbar(); title("true saturation")
tight_layout()
safesave(joinpath(plotsdir(sim_name, exp_name), savename(fig_name; digits=6)*"_saturation.png"), fig);
close(fig)

fig=figure(figsize=(20,12));
subplot(1,3,1);
imshow(reshape(state_init.states[end].Pressure, n[1], n[end])', vmin=minimum(state.states[end].Pressure), vmax=maximum(state.states[end].Pressure)); colorbar(); title("initial pressure")
subplot(1,3,2);
imshow(reshape(state0.states[end].Pressure, n[1], n[end])', vmin=minimum(state.states[end].Pressure), vmax=maximum(state.states[end].Pressure)); colorbar(); title("inverted pressure")
subplot(1,3,3);
imshow(reshape(state.states[end].Pressure, n[1], n[end])', vmin=minimum(state.states[end].Pressure), vmax=maximum(state.states[end].Pressure)); colorbar(); title("true pressure")
tight_layout()
safesave(joinpath(plotsdir(sim_name, exp_name), savename(fig_name; digits=6)*"_pressure.png"), fig);
close(fig)

fig=figure(figsize=(20,12));
plot(fhistory); title("loss")
tight_layout()
safesave(joinpath(plotsdir(sim_name, exp_name), savename(fig_name; digits=6)*"_loss.png"), fig);
close(fig)

x0 = T(K0)
x_init = T(K_init)
x = T(K)

## plotting
trans_x = reshape(x[1:(n[1]-1)*n[end]], n[1]-1, n[end])
trans_z = reshape(x[(n[1]-1)*n[end]+1:end], n[1], n[end]-1)

trans_x0 = reshape(x0[1:(n[1]-1)*n[end]], n[1]-1, n[end])
trans_z0 = reshape(x0[(n[1]-1)*n[end]+1:end], n[1], n[end]-1)

trans_x_init = reshape(x_init[1:(n[1]-1)*n[end]], n[1]-1, n[end])
trans_z_init = reshape(x_init[(n[1]-1)*n[end]+1:end], n[1], n[end]-1)

fig=figure(figsize=(20,12));
subplot(1,3,1);
imshow(trans_x_init', vmin=minimum(trans_x), vmax=maximum(trans_x)); colorbar(); title("initial transmissibility x")
subplot(1,3,2);
imshow(trans_x0', vmin=minimum(trans_x), vmax=maximum(trans_x)); colorbar(); title("inverted transmissibility x")
subplot(1,3,3);
imshow(trans_x', vmin=minimum(trans_x), vmax=maximum(trans_x)); colorbar(); title("true transmissibility x")
tight_layout()
safesave(joinpath(plotsdir(sim_name, exp_name), savename(fig_name; digits=6)*"_transx.png"), fig);
close(fig)

fig=figure(figsize=(20,12));
subplot(1,3,1);
imshow(trans_z_init', vmin=minimum(trans_z), vmax=maximum(trans_z)); colorbar(); title("initial transmissibility z")
subplot(1,3,2);
imshow(trans_z0', vmin=minimum(trans_z), vmax=maximum(trans_z)); colorbar(); title("inverted transmissibility z")
subplot(1,3,3);
imshow(trans_z', vmin=minimum(trans_z), vmax=maximum(trans_z)); colorbar(); title("true transmissibility z")
tight_layout()
safesave(joinpath(plotsdir(sim_name, exp_name), savename(fig_name; digits=6)*"_transz.png"), fig);
close(fig)

