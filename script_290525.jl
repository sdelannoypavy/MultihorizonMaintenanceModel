
using Dates
using CSV

include("optimization.jl")

fichier_csv = "computing_time.csv"

# Pfixed parameters
P0 = 100
Nb_scen_prod = 1
Quota = 1

res = []

function test(fail_prob_list, S_list, nb_state_list)
    # warning: Ns is number of scenarios for each jumb probability, not total number of scenarios

    for p in fail_prob_list
        res_T = []
        
        for S in S_list

            for nb_state in nb_state_list
                local t_start = now()
                
                num_vars, num_constrs = V(p,S,nb_state)

                local t_end = now()
                local t_elapsed = floor(Int,(t_end - t_start)/Second(1))
                push!(res_T, t_elapsed)
                CSV.write(fichier_csv, [(; param1=p, param2=S, param3=nb_state, param4=num_vars, param5=num_constrs, param6=t_elapsed)], append=true)
            end

        end
        push!(res, res_T)

    end

end

# testing parameters
fail_prob_list = [0.05]
S_list = [1,5,10]
nb_state_list = [2]

test(fail_prob_list, S_list,nb_state_list)