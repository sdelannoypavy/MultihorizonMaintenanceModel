
using Dates
using CSV

include("optimization.jl")

fichier_csv = "computing_time.csv"

Random.seed!(1)


# Pfixed parameters
P0 = 100
Nb_scen_prod = 1
Quota = 1

res = []

function test(fail_prob_list, S_list)
    # warning: Ns is number of scenarios for each jumb probability, not total number of scenarios

    for p in fail_prob_list
        res_T = []
        
        for S in S_list
            local t_start = now()
            
            num_vars, num_constrs = V(p,S)

            local t_end = now()
            local t_elapsed = floor(Int,(t_end - t_start)/Second(1))
            push!(res_T, t_elapsed)
            CSV.write(fichier_csv, [(; param1=p, param3=S, param3=num_vars, param4=num_constrs)], append=true)
        end
        push!(res, res_T)

    end

end

# testing parameters
fail_prob_list = [0.0005]
S_list = [1,3]

test(fail_prob_list, S_list)