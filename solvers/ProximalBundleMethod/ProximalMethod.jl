abstract type AbstractMethod end

struct Bundle
	ref	# constraint/variable reference
	y	# evaluation point
	fy	# evaluation value
	g	# subgradient
end

#=
	This structure contains necessary information for bundle methods.
=#
mutable struct Model{T<:AbstractMethod}
	n::Int64		# dimension in the solution space
	m::JuMP.Model	# Bundle model
	k::Int64 		# iteration counter
	maxiter::Int64	# maximum number of iterations

	y::Vector{Float64}		# current iterate
	fy::Float64	# objective values at the current iterate
	g::Array{Float64,1}		# subgradients

	# user-defined function to evaluate f
	# and return the value and its subgradients
	evaluate_f
	# History of bundles
	history::Dict{Tuple{Int64,Int64},Bundle}
    tol::Float64
	solver

	# Placeholder for extended structures
	ext

	function Model{T}(n::Int64, func,tol) where {T<:AbstractMethod}
		bundle = new()
		bundle.n = n
		bundle.m = JuMP.Model(with_optimizer(OSQP.Optimizer,verbose=0))
		bundle.k = 0
		bundle.maxiter = 500000
		bundle.y = zeros(n)
        bundle.tol=tol
		bundle.fy = 0
		bundle.g = zeros(0)
		bundle.evaluate_f = func
		bundle.history = Dict{Tuple{Int64,Int64},Bundle}()

		# initialize bundle model
		initialize!(bundle)

		return bundle
	end
end

function runb(bundle::Model{<:AbstractMethod})

	add_initial_bundles!(bundle)
	bundle.k += 1

	while true
		status = solve_bundle_model(bundle)
		if status != MOI.OPTIMAL
			println("TERMINATION: Invalid status from bundle model.")
			break
		end
		display_info!(bundle)
		if termination_test(bundle)
			break
		end
		evaluate_functions!(bundle)
		manage_bundles!(bundle)
		update_iteration!(bundle)
	end
end

#=
	Abstract functions are defined below.
	For each method, the functions may be implemented.
=#

const AbstractModel = Model{AbstractMethod}

function initialize!(bundle::AbstractModel)
end

function add_initial_bundles!(bundle::AbstractModel)
end

function solve_bundle_model(bundle::AbstractModel)
	return :Optimal
end

function termination_test(bundle::AbstractModel)
	return true
end

function evaluate_functions!(bundle::AbstractModel)
end

function manage_bundles!(bundle::AbstractModel)
end

function update_iteration!(bundle::AbstractModel)
end

function display_info!(bundle::AbstractModel)
end

getsolution(bundle::AbstractModel)::Array{Float64,1} = Array{Float64,1}(undef,bundle.n)
getobjectivevalue(bundle::AbstractModel)::Float64 = NaN

#=
 	Add the implementations of bundle methods
=#





#=
	Implementation of Proximal Bundle Method.

	The implementation is based on
	Krzysztof C. Kiwiel, "Proximity control in bundle methods for convex nondifferentiable minimization"
	Mathematical Programming 46(1-3), 1990
=#

#=
TODO: purge_cut() is not functioning correctly when the model is modified by user.
So, this is disabled.
=#
const disable_purge_cuts = true

abstract type ProximalMethod <: AbstractMethod end

# Algorithm-specific structure
mutable struct ProximalModelExt
	# Algorithm-specific parameters
	u::Float64
	u_min::Float64
	M_g::Int64
	ϵ_float::Float64	# tolerance for floating point comparison
	ϵ_s::Float64
	ϵ_v::Float64
	m_L::Float64
	m_R::Float64

	x0::Array{Float64,1}	# current best solution (at iteration k)
	x1::Array{Float64,1}	# new best solution (at iteration k+1)
	fx0::Float64	# current best objective values
	fx1::Float64	# new best objective values (at iteration k+1)
	d::Array{Float64,1}
	v::Float64
	i::Int64
	α::Float64

	function ProximalModelExt(n::Int64)
		ext = new()
		ext.u = 1e3#0.1
		ext.u_min = 1e-2#1.0e-2
		ext.M_g = 1e+6
		ext.ϵ_float = 1.0e-8
		ext.ϵ_s = 1e-6
		ext.ϵ_v = Inf
		ext.m_L = 1e-4#1.0e-4
		ext.m_R = 0.5#0.5
		ext.x0 = zeros(n)
		ext.x1 = zeros(n)
		ext.fx0 = 0
		ext.fx1 = 0
		ext.d = zeros(n)
		ext.v = 0
		ext.i = 0
		ext.α = 0
		return ext
	end
end

const ProximalModel = Model{ProximalMethod}

function initialize!(bundle::ProximalModel)
	# Attach the extended structure
	bundle.ext = ProximalModelExt(bundle.n)
    bundle.ext.ϵ_s=bundle.tol
	# create the initial bundle model



	@variable(bundle.m,x[1:bundle.n])

	@variable(bundle.m, θ)
	@objective(bundle.m, Min,θ+ 0.5 * bundle.ext.u * sum((x[i] - bundle.ext.x0[i])^2 for i in 1:bundle.n))
end

function add_initial_bundles!(bundle::ProximalModel)

    bundle.fy, bundle.g=bundle.evaluate_f(bundle.y)
    bundle.ext.fx0 = copy(bundle.fy)


	# add bundles

    add_cut(bundle, bundle.g, bundle.fy, bundle.y)

end

function solve_bundle_model(bundle::ProximalModel)
	# solve the bundle model
    optimize!(bundle.m)
    status = termination_status(bundle.m)



	if status == MOI.OPTIMAL
		# variable references
		x = getindex(bundle.m, :x)
		θ = getindex(bundle.m, :θ)

		# get solutions

        bundle.y = JuMP.value.(x)
        bundle.ext.d = bundle.y - bundle.ext.x0


        bundle.ext.v = JuMP.value(θ) - bundle.ext.fx0
        return status
	end


end

function termination_test(bundle::Model{<:ProximalMethod})
	if bundle.ext.v >= -bundle.ext.ϵ_s * (1 + abs(bundle.ext.fx0))
		println("TERMINATION: Optimal: v = ", bundle.ext.v)
		return true
	end
	if bundle.k >= bundle.maxiter
		println("TERMINATION: Maximum number of iterations reached.")
		return true
	end
	return false
end

function evaluate_functions!(bundle::Model{<:ProximalMethod})
	# evaluation function f
	bundle.fy, bundle.g = bundle.evaluate_f(bundle.y)

	# descent test
	descent_test(bundle)
end

function manage_bundles!(bundle::ProximalModel)
	if disable_purge_cuts
		ncuts_purged = purge_cuts(bundle)
	end

	# add cuts

    gd= bundle.g' * bundle.ext.d
    bundle.ext.α = bundle.ext.fx0 - (bundle.fy - gd)
    if -bundle.ext.α + gd > bundle.ext.v + bundle.ext.ϵ_float
        add_cut(bundle, bundle.g, bundle.fy, bundle.y)
    end

end

function update_iteration!(bundle::ProximalModel)
	# update u
	updated = update_weight(bundle)

	# Update objective function
	if updated
		update_objective!(bundle)
	end

	bundle.k += 1

	bundle.ext.x0 = copy(bundle.ext.x1)
	bundle.ext.fx0 = copy(bundle.ext.fx1)
end

getsolution(bundle::Model{<:ProximalMethod})::Array{Float64,1} = bundle.ext.x0
getobjectivevalue(bundle::Model{<:ProximalMethod})::Float64 = bundle.ext.fx0

function descent_test(bundle::Model{<:ProximalMethod})
	if bundle.fy <= bundle.ext.fx0 + bundle.ext.m_L * bundle.ext.v
		bundle.ext.x1 = copy(bundle.y)
		bundle.ext.fx1 = copy(bundle.fy)
	else
		bundle.ext.x1 = copy(bundle.ext.x0)
		bundle.ext.fx1 = copy(bundle.ext.fx0)
	end
end

function add_cut(bundle::ProximalModel, g::Array{Float64,1}, fy::Float64, y::Array{Float64,1}; store_cuts = true)
	x = getindex(bundle.m, :x)
	θ = getindex(bundle.m, :θ)
	constr = @constraint(bundle.m, fy + sum(g[i] * (x[i] - y[i]) for i in 1:bundle.n) <= θ)
	if !disable_purge_cuts && tore_cuts
		bundle.history[bundle.k] = Bundle(constr, deepcopy(y), fy, g)
	end
end

function purge_cuts(bundle::ProximalModel)
	ncuts = length(bundle.history)
	ncuts_to_purge = ncuts - bundle.ext.M_g
	cuts_to_purge = Tuple{Int64,Int64}[]
	if ncuts_to_purge > 0
		for (refkey,hist) in bundle.history
			if getdual(hist.ref) < -1.0e-3#8
				push!(cuts_to_purge, refkey)
				ncuts_to_purge -= 1
			end
			if ncuts_to_purge <= 0
				break
			end
		end
	end

	if length(cuts_to_purge) > 0
		for refkey in cuts_to_purge
			delete!(bundle.history, refkey)
		end

		solver = bundle.m.solver
		bundle.m = Model(solver=solver)
        @variable(bundle.m,x[i=1:bundle.n])

        @variable(bundle.m, θ)
		@objective(bundle.m, Min,θ + 0.5 * bundle.ext.u * sum((x[i] - bundle.ext.x1[i])^2 for i in 1:bundle.n))

		for (k,h) in bundle.history
			add_cut(bundle, h.g, h.fy, h.y, k[1], store_cuts = false)
			delete!(bundle.history, k)
		end
	end

	return length(cuts_to_purge)
end

function update_weight(bundle::Model{<:ProximalMethod})
	updated = false

	# update weight u
	u = bundle.ext.u
	if bundle.fy <= bundle.ext.fx0 + bundle.ext.m_L * bundle.ext.v
		if bundle.ext.i > 0 && bundle.fy <= bundle.ext.fx0 + bundle.ext.m_R * bundle.ext.v
			u = 2 * bundle.ext.u * (1 - (bundle.fy - bundle.ext.fx0) / bundle.ext.v)
		elseif bundle.ext.i > 3
			u = bundle.ext.u / 2
		end

		newu = max(u, bundle.ext.u/10, bundle.ext.u_min)
		bundle.ext.ϵ_v = max(bundle.ext.ϵ_v, -2*bundle.ext.v)
		bundle.ext.i = max(bundle.ext.i+1,1)
		if newu != bundle.ext.u
			bundle.ext.i = 1
		end
		updated = true
	else
		p = Compat.norm(-bundle.ext.u .* bundle.ext.d, 1)
		α_tilde = -p^2 / bundle.ext.u - bundle.ext.v

		bundle.ext.ϵ_v = min(bundle.ext.ϵ_v, p + α_tilde)
		if bundle.ext.α > max(bundle.ext.ϵ_v, -10*bundle.ext.v) && bundle.ext.i < -3
			u = 2 * bundle.ext.u * (1 - (bundle.fy - bundle.ext.fx0) / bundle.ext.v)
		end
		newu = min(u, 10*bundle.ext.u)
		bundle.ext.i = min(bundle.ext.i-1,-1)
		if newu != bundle.ext.u
			bundle.ext.i = -1
		end
	end
	# newu = 1.0e-6
	if newu != bundle.ext.u
		bundle.ext.u = newu
		updated = true
	end

	return updated
end

function update_objective!(bundle::ProximalModel)
	x = getindex(bundle.m, :x)
	θ = getindex(bundle.m, :θ)
	@objective(bundle.m, Min, θ+ 0.5 * bundle.ext.u * sum((x[i] - bundle.ext.x1[i])^2 for i in 1:bundle.n))
end


function display_info!(bundle::Model{<:ProximalMethod})
	Compat.Printf.@printf("Iter %d: fx0 %e, fx1 %e, fy %e, v %e, u %e, i %d\n",
	bundle.k, bundle.ext.fx0, bundle.ext.fx1, bundle.fy, bundle.ext.v, bundle.ext.u, bundle.ext.i)
end
