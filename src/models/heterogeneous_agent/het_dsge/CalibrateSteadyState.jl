# this should be run after running GovDebt.jl
# targeted moments:
varlinc_target  = 0.7  # var(log(annual income))
vardlinc_target = 0.23 # var(one year changes in log(annual income))
mpc_target      = 0.16
pc0_target      = 0.10 # fraction with zero bond holdings
# parameters: pLH, pHL, zlo, sH/sL. We normalize so s and z both have mean 1.

using Distributions
#using PyPlot
using Roots
using JLD

function mollifier(z::AbstractFloat,ehi::AbstractFloat,elo::AbstractFloat)
    # mollifier function
    In = 0.443993816237631
    if z<ehi && z>elo
        temp = -1.0 + 2.0*(z-elo)/(ehi-elo)
        out = (2.0/(ehi-elo))*exp(-1.0/(1.0-temp*temp))/In
    else
        out = 0.0
    end

    return out
end

ni = 10000 # number of individuals in the sample we use for simulating skill process
uz = rand(ni,8) # uniform random variable used to get transitory shock z; need 8 quarters to get 1 year change
us = rand(ni,8) # uniform random variable used to get persistent markov chain shock s

# we will create a single z sample which has zlo=0, zhi = 2. We will then rescale it as we change zlo.
nz = 1000
zgrid = collect(linspace(0.,2.,nz))
zprob = zeros(nz)
for i=1:nz
	zprob[i] = 2*mollifier(zgrid[i],2.,0.)/nz
end
zprob = zprob/sum(zprob)
zcdf = cumsum(zprob)

zave = 0.5*zgrid[1:nz-1]+0.5*zgrid[2:nz]

function zsample(uz::Array{Float64,2}, zgrid::AbstractArray, zcdf::AbstractArray, ni::Int, nz::Int)
	zs = zeros(ni,8)
	for i=1:ni
		for t=1:8
			for iz=1:nz-1
				if uz[i,t]>zcdf[iz] && uz[i,t] <= zcdf[iz+1]
					zs[i,t] = zave[iz]
				end
			end
		end
	end
	return zs
end
tic()
zs = zsample(uz,zgrid,zcdf,ni,nz)
toc()

function ssample(us::Array{Float64,2}, P::Array{Float64,2}, πss::AbstractArray, ni::Int)
	shist = ones(Int,ni,8)
	for i=1:ni
		if us[i,1] > πss[1]
			shist[i,1] = 2
		end
		for t=2:8
			if us[i,t] > P[shist[i,t-1],1]
				shist[i,t] = 2
			end
		end
	end
	return shist
end

function ln_annual_inc(zhist::Array{Float64,2}, us::Array{Float64,2}, zlo::AbstractFloat, P::Array{Float64,2}, πss::AbstractArray, sgrid::AbstractArray, ni::Int)
  	s_inds = ssample(us,P,πss,ni)
	linc1 = zeros(ni)
	linc2 = zeros(ni)
	for i=1:ni
		inc1 = 0.
		for t=1:4
#            @show zhist[i, t]-1
#            @show zlo
#            @show sgrid[1]
#            @show sgrid[2]
#            @show s_inds[i,t]-1
			zshock = 1. + (1. - zlo)*(zhist[i,t]-1.)
			sshock = (sgrid[1] + (sgrid[2] - sgrid[1])*(s_inds[i,t] - 1.))
#            @show zshock
#            @show sshock
			inc1 += zshock*sshock
		end
		inc2 = 0.
		for t=5:8
			zshock = 1. + (1. - zlo)*(zhist[i,t]-1.)
			sshock = (sgrid[1] + (sgrid[2] - sgrid[1])*(s_inds[i,t] - 1.))
			inc2 += zshock*sshock
		end
		linc1[i] += log(inc1)
     	linc2[i] += log(inc2)
	end
	return (linc1, linc2)
end

function skill_moments(sH_by_sL::Real, zlo::Real, pLH::AbstractFloat, pHL::AbstractFloat, zhist::Array{Float64,2}, us::Array{Float64,2}, ni::Int = 10000)
	πL = pHL/(pLH+pHL)
	πss = [πL;1-πL]
	P = [[1.-pLH pLH];[pHL 1.-pHL]]
	slo = 1./(πL+(1-πL)*sH_by_sL)
	shi = sH_by_sL*slo
	sgrid = [slo; shi]
    @show zlo
	(linc1, linc2) = ln_annual_inc(zs,us,zlo,P,πss,sgrid,ni)
	return (var(linc1), var(linc2 - linc1))
end

tic()
(varlinc, vardlinc) = skill_moments(1.2, 0.5, 0.1, 0.1, zs, us)
toc()


function best_fit(pLH::AbstractFloat, pHL::AbstractFloat, varlinc_target::AbstractFloat, vardlinc_target::AbstractFloat, sH_by_sL_grid::AbstractArray, zlo_grid::AbstractArray,
				  n1::Int, n2::Int, zhist::Array{Float64,2}, us::Array{Float64,2}, ni::Int = 10000)
	n = n1*n2
	dist = zeros(n)
	varlinc = zeros(n)
	vardlinc = zeros(n)
	for i1=1:n1
		for i2=1:n2
			i = n2*(i1-1)+i2
			(varlinc[i],vardlinc[i]) = skill_moments(sH_by_sL_grid[i1], zlo_grid[i2], pLH, pHL, zs, us)
			dist[i] = abs(varlinc[i] - varlinc_target) + abs(vardlinc[i] - vardlinc_target)
		end
	end
	(mindist,imin) = findmin(dist)
	i2 = mod(imin-1,n2)+1
	i1 = Int((imin-i2)/n2 + 1)
	return (sH_by_sL_grid[i1], zlo_grid[i2], varlinc[imin], vardlinc[imin])
end

np = 100
sH_by_sL_grid = linspace(3.,9.,np)
zlo_grid      = linspace(0.+eps(), 0.8-eps(), np)
# tic()
# varlinc_target  = 0.7
# vardlinc_target = 0.23
# (sH_by_sL, zlo, varlinc, vardlinc) = best_fit(0.02,0.03,varlinc_target,vardlinc_target,sH_by_sL_grid, zlo_grid,np,np,zs,us)
# toc()

# out = [varlinc varlinc_target vardlinc vardlinc_target sH_by_sL zlo]

"""
THIS ONE GETS CALLED
"""
function best_fit(pLH::AbstractFloat, pHL::AbstractFloat, varlinc_target::AbstractFloat, vardlinc_target::AbstractFloat, sH_by_sL_grid::AbstractArray, zlo_grid::AbstractArray,
				  n1::Int, n2::Int, zhist::Array{Float64,2}, us::Array{Float64,2}, T::Int = 8, ni::Int = 10000)
	n = n1*n2
	dist = zeros(n)
	varlinc = zeros(n)
	vardlinc = zeros(n)
	for i1=1:n1
		for i2=1:n2
			i = n2*(i1-1)+i2
			(varlinc[i], vardlinc[i]) = skill_moments(sH_by_sL_grid[i1], zlo_grid[i2], pLH, pHL, zs, us)

			dist[i] = abs(varlinc[i] - varlinc_target) + abs(vardlinc[i] - vardlinc_target)
		end
	end
    print(dist)
	(mindist,imin) = findmin(dist)
	i2 = mod(imin-1,n2)+1
	i1 = Int((imin-i2)/n2 + 1)
	return (sH_by_sL_grid[i1], zlo_grid[i2], varlinc[imin], vardlinc[imin])
end

loss(x::Vector{Float64}, target::Vector{Float64}) = sum(abs.(x-target))
function grad(f, x::Matrix{Float64}, h::Float64 = 0.001)
    N = length(x)
    ∇fx = Matrix{Float64}(undef, N, 2)
    fx = [i for i in f(x)]
    for i=1:N
        x[i]  += h
        ∇fx[i, :] = ([i for i in f(x)] - fx) / h
        x[i]  -= h
    end
    return ∇fx, fx
end
function best_fit(pLH::AbstractFloat, pHL::AbstractFloat, varlinc_target::AbstractFloat, vardlinc_target::AbstractFloat, sH_by_sL_lo::AbstractFloat, zlo_lo::AbstractFloat,
                  zlo_hi::AbstractFloat, us::Array{Float64,2}, max_iter::Int = 50, initial_guess::Vector{Float64} = [6.3, 0.03], α::AbstractFloat=0.001, η::AbstractFloat = 0.1)
	min_varlinc  = 0.0
	min_vardlinc = 0.0
    min_dist     = 20.0

    α0 = α
    w  = initial_guess #'; initial_guess']
#    @show w
    iter_w_no_improv = 0

    @show w[2]
    skill_moments_f(x) = skill_moments(x[1], x[2], pLH, pHL, zs, us)

    while iter_w_no_improv < max_iter
        ∇fx, fx = grad(skill_moments_f, w)
        target = [varlinc_target; vardlinc_target]
        dist = loss(fx, target)
        if dist < min_dist
            min_dist = dist
            min_varlinc, min_vardlinc = fx[1], fx[2]
            iter_w_no_improv = 0
            α0 = α
        else
            iter_w_no_improv += 1
            α0 *= 0.9
        end
        #is the solution just to return the derivative from \grad fx?
#        @show -η*∇fx
        w += -η*∇fx

        if w[1] < sH_by_sL_lo
            w[1] = sH_by_sL_lo
        end
        if w[2] < zlo_lo
            w[2] = zlo_lo
        elseif w[2] > zlo_hi
            w[2] = zlo_hi
        end
    end
    sH_by_sL_argmin = w[1]
    zlo_argmin = w[2]
    print(min_dist)

	return (sH_by_sL_argmin, zlo_argmin, min_varlinc, min_vardlinc)
end

function ave_mpc(m::AbstractArray, c::AbstractArray, agrid::AbstractArray, aswts::AbstractArray, na::Int, ns::Int)
	mpc = 0.
	for iss=1:ns
		i = na*(iss-1)+1
		mpc += aswts[i]*m[i]*(c[i+1] - c[i])/(agrid[2] - agrid[1])
		for ia=2:na
			i = na*(iss-1)+ia
			mpc += aswts[i]*m[i]*(c[i] - c[i-1])/(agrid[ia] - agrid[ia-1])
		end
	end
	return mpc
end

function frac_zero(m::AbstractArray, c::AbstractArray, agrid::AbstractArray, aswts::AbstractArray, ns::Int)
	return sum(aswts.*m.*(c.==repeat(agrid,ns)))
end

nLH = 2
nHL = 2

pLH_grid = linspace(0.01,0.1,nLH)
pHL_grid = linspace(0.01,0.1,nHL)
nn = nLH*nHL
sH_by_sL_grid2 = zeros(nn)
zlo_grid2 = zeros(nn)
varlinc_grid = zeros(nn)
vardlinc_grid = zeros(nn)
mpc_grid  = zeros(nn)
pc0_grid  = zeros(nn)

# RECA
sH_by_sL_lo = 3.0
sH_by_sL_hi = 9.0
zlo_lo      = 0.+eps()
zlo_hi      =  0.8-eps()
@show "1=============================================================="
for iLH = 1:nLH
	for iHL = 1:nHL
		i = nHL*(iLH-1)+iHL
        @time out = best_fit(pLH_grid[iLH], pHL_grid[iHL], varlinc_target, vardlinc_target, sH_by_sL_lo, zlo_lo, zlo_hi, us)
@show "2=============================================================="
		@time (sH_by_sL_grid2[i], zlo_grid2[i], varlinc_grid[i], vardlinc_grid[i]) = best_fit(pLH_grid[iLH],pHL_grid[iHL],varlinc_target,vardlinc_target,sH_by_sL_grid, zlo_grid,np,np,zs,us)
		(f, sgrid, swts) = persistent_skill_process(sH_by_sL_grid2[i], pLH_grid[iLH], pHL_grid[iHL], ns)
		(agrid, awts) = cash_grid(sgrid, ω, H, r, η, γ, T, zlo, na)
		aswts = kron(swts, awts)
		qfunction(x) = mollifier(x, 2. - zlo_grid[i], zlo_grid[i])
		println("parameters:")
		println([sH_by_sL_grid2[i] zlo_grid2[i] pLH_grid[iLH] pHL_grid[iHL]])
		(ell, c, m, β, report) = findss(na, ns, βlo, βhi, R, ω, H, η, T, γ, bg, qfunction, agrid, sgrid,aswts, Win, f)
		mpc_grid[i] = ave_mpc(m,c,agrid,aswts,na,ns)
		pc0_grid[i] = frac_zero(m,c,agrid,aswts,ns)
	end
end
