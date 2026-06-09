#!/usr/bin/env julia
using Arrow
using DataFrames
using MixedModels
using GLM
using CSV
using ArgParse
using Random
using PooledArrays
using LinearAlgebra # Required for BLAS settings

function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--arrow"
            help = "Path to input Arrow file"
            required = true
        "--nSims"
            help = "Number of simulations"
            arg_type = Int
            default = 1000
        "--out"
            help = "Output CSV"
            required = true
    end
    return parse_args(s)
end

function check_convergence(m)
    success_statuses = (:SUCCESS, :FTOL_REACHED, :XTOL_REACHED, :STOPVAL_REACHED)
    return m.optsum.returnvalue in success_statuses
end

function main()
    BLAS.set_num_threads(1)

    args = parse_commandline()

    if !isfile(args["arrow"])
        exit(1)
    end

    df = DataFrame(Arrow.Table(args["arrow"]))
    df.pop = PooledArray(df.pop)
    df.trt = PooledArray(df.trt)
    
    results = DataFrame(
        chrom = String[], pos = Int[], 
        LRT_chisq = Float64[], PB_p_val = Float64[], 
        singular = Bool[], converged = Bool[], 
        sim_fail_rate = Float64[], error = String[]
    )
    
    results_lock = ReentrantLock()

    gdf = groupby(df, [:chrom, :pos])
    
    all_keys = collect(keys(gdf))
    
    Threads.@threads for key in all_keys
        
        rng = MersenneTwister(1234 + Threads.threadid())

        subdf = gdf[key]
        
        if nrow(subdf) == 0 || length(unique(subdf.pop)) < 2
            lock(results_lock) do
                push!(results, (string(key.chrom), key.pos, NaN, NaN, false, false, NaN, "NotEnoughData"))
            end
            continue
        end

        wts = Float64.(subdf.total)
        f_null = @formula(success / total ~ gen + trt + (1 | pop))
        f_full = @formula(success / total ~ gen * trt + (1 | pop))
        
        try
            m0 = fit(MixedModel, f_null, subdf, Binomial(), wts=wts; fast=true)
            m1 = fit(MixedModel, f_full, subdf, Binomial(), wts=wts; fast=true)
            
            is_sing = issingular(m0) || issingular(m1)
            is_conv = check_convergence(m0) && check_convergence(m1)
            
            lrt_obs = 2 * (loglikelihood(m1) - loglikelihood(m0))
            if lrt_obs < 0 lrt_obs = 0.0 end
            
            if !is_conv
                lock(results_lock) do
                    push!(results, (string(key.chrom), key.pos, lrt_obs, NaN, is_sing, false, NaN, "ObservedFitFailed"))
                end
                continue
            end

            m0_frozen = deepcopy(m0) 

            m0.optsum.optimizer = :LN_NELDERMEAD
            m0.optsum.maxfeval = 10000 
            m0.optsum.maxtime = 2.0 
            
            m1.optsum.optimizer = :LN_NELDERMEAD
            m1.optsum.maxfeval = 10000
            m1.optsum.maxtime = 2.0 

            better_count = 0
            valid_sims = 0
            consecutive_failures = 0
            aborted = false
            
            for i in 1:args["nSims"]
                
                y_sim = simulate(rng, m0_frozen)
                
                try
                    refit!(m0, y_sim; fast=true)
                    if !check_convergence(m0)
                        consecutive_failures += 1
                        if consecutive_failures > 50 && valid_sims == 0
                             aborted = true
                             break
                        end
                        continue 
                    end
                    ll0 = loglikelihood(m0)
                    
                   refit!(m1, y_sim; fast=true)
                    if !check_convergence(m1)
                        consecutive_failures += 1
                        continue
                    end
                    ll1 = loglikelihood(m1)
                    
                    consecutive_failures = 0
                    
                    lrt_sim = 2 * (ll1 - ll0)
                    if lrt_sim < 0 lrt_sim = 0.0 end
                    
                    if lrt_sim >= lrt_obs
                        better_count += 1
                    end
                    valid_sims += 1
                    
                catch e
                    consecutive_failures += 1
                    continue 
                end
            end
            
            lock(results_lock) do
                if aborted
                    push!(results, (string(key.chrom), key.pos, lrt_obs, NaN, is_sing, true, 1.0, "Aborted_TooManyFailures"))
                else
                    fail_rate = 1.0 - (valid_sims / args["nSims"])
                    if valid_sims > (args["nSims"] * 0.5) 
                        pval = (better_count + 1) / (valid_sims + 1)
                        push!(results, (string(key.chrom), key.pos, lrt_obs, pval, is_sing, true, fail_rate, "OK"))
                    else
                        push!(results, (string(key.chrom), key.pos, lrt_obs, NaN, is_sing, true, fail_rate, "HighSimFailure"))
                    end
                end
            end

        catch e
            err_msg = replace(string(e), "\n" => " ")[1:min(end, 200)]
            lock(results_lock) do
                push!(results, (string(key.chrom), key.pos, NaN, NaN, false, false, NaN, err_msg))
            end
        end
    end
    
    CSV.write(args["out"], results)
    println("Done.")
end

main()