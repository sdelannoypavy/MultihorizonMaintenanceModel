using Gurobi
using JuMP  
using Distributions
using Random
using StatsBase

function V(s,q_i,S,h,product,nb_state,P_w,P_m,Q,Vals)

    # Components have the same degradation model

    # Random initial state (might not be representative of the reality...)
    x_0 = zeros(nb_state)      
    x_0[s] = 1 

    T = 60

    # maintenance duration is the same for all components
    d = max(s - 1,1)

    model = Model(Gurobi.Optimizer)

    @variable(model, x[1:(S+1),1:(T+1),1:nb_state] >= 0)

    @variable(model, c[1:(T+1)] >= 0)
    @variable(model, p_evac[1:S,1:T,1:nb_state] >= 0)

    @variable(model, m[1:T], Bin)
    @variable(model, u[1:(S+1),1:T], Bin)
    @variable(model, ong_m[1:(S+1),1:T], Bin) #ongoing maintenance
    @variable(model, q[1:(S+1),1:T] >= 0)
    @variable(model, k[1:T], Bin)
    @variable(model, q_f, Int)
    @variable(model, δ[0:Q], Bin)    

    #@objective(model, Min, sum(c[t] for t in 1:T) + sum(x[S+1,T+1,c,i]*V_table[i] for c in 1:nb_comp, i in 1:nb_state))
    @objective(model, Min, sum(c[t] for t in 1:(T+1)))


    @constraint(model, [s in 1:(S+1), t in 2:(T+1), i in 1:nb_state], x[s,t,i] == (1-u[s,t-1])*sum(P_w[i,j]*x[s,t-1,j] for j in 1:nb_state) + u[s,t-1]*sum(P_m[i,j]*x[s,t-1,j] for j in 1:nb_state))

    @constraint(model, [s in 1:(S+1), i in 1:nb_state], q[s,1] == d*m[1])
    @constraint(model, [s in 1:(S+1), t in 1:(T-1)], q[s,t+1] == q[s,t] - u[s,t] + d*m[t+1])
    @constraint(model, [s in 1:S, t in 1:(T)], u[s,t] + 1 - h[s,t] >= (1/T)*q[s,t])
    @constraint(model, [t in 1:(T)], u[S+1,t] >= (1/T)*q[S+1,t])
    @constraint(model, [s in 1:(S+1), t in 1:(T)], u[s,t] <= q[s,t])
    @constraint(model, [s in 1:S, t in 1:(T)], u[s,t] <= h[s,t])
    @constraint(model, [s in 1:(S+1), t in 1:(T)], ong_m[s,t] >= (1/T)*q[s,t])

    @constraint(model, [s in 1:(S+1), i in 1:nb_state], x[s,1,i] == x_0[i])

    # proba of curtailement for failure = proba of union x[c] = n^c
    @constraint(model, [t in 1:T], c[t] >= (1/S)*sum((1 - ong_m[s,t])*x[s,t,nb_state]*product[s,t] + ong_m[s,t]*(1-k[t])*product[s,t] for s in 1:S))

    @constraint(model, q_f == q_i - sum(k[t] for t in 1:T)) 
    @constraint(model, q_f >= 0)   
    M = 1e6            
    @constraint(model, [j in 0:Q], q_f - j <=  M * (1 - δ[j]))
    @constraint(model, [j in 0:Q], q_f - j >= -M * (1 - δ[j]))

    @constraint(model, [j in 0:Q], c[T+1] + M * (1 - δ[j]) >= (1/S)*sum(sum(x[s,T+1, i] * Vals[i, j+1] for i in 1:nb_state) for s in 1:S))


    optimize!(model)

    status = termination_status(model)
    println("Statut de l'optimisation: $status")

    cost = objective_value(model)

    #numvars = MOI.get(model, Gurobi.ModelAttribute("NumVars"))
    #numcons = MOI.get(model, Gurobi.ModelAttribute("NumConstrs"))


    m = [value(m[t]) for t in 1:T]
    k = [value(m[t]) for t in 1:T]
    #u = [value(u[s,t]) for  s in 1:S, t in 1:T]
    #q = [value(q[s,t]) for  s in 1:S, t in 1:T]
    #ong_m = [value(ong_m[s,t]) for  s in 1:S, t in 1:T]
    #x = [value(x[s,t,nb_state]) for  s in 1:S, t in 1:(T+1)]
    #c = [value(c[t]) for t in 1:T]
    #sum_m = sum(value(m[t]) for t in 1:T)

    if status == MOI.OPTIMAL
        return(cost, m,k)
    else
        println("Aucune solution optimale trouvée.")
    end    

end


function Bellman(years,Q,S,h,product,nb_state,P_w,P_m)

    Tmax = 6*years

    Vals_old = zeros(nb_state, Q+1)

    # Stockage des valeurs à chaque pas de temps (optionnel)
    Vals_all = Array{Float64, 3}(undef, nb_state, Q+1, Tmax)

    # Stockage des politiques optimales à chaque pas de temps
    Policies_m = Array{Float64}(undef, nb_state, Q+1, Tmax, 60)
    Policies_k = Array{Float64}(undef, nb_state, Q+1, Tmax, 60)

    for t in 1:Tmax
        Vals_new = zeros(nb_state, Q+1)
        for s in 1:nb_state
            for q in 0:Q
                h_t = h[t, :, :]  # taille (S, T)
                product_t = product[t, :, :]
                val, m_opt, k_opt = V(s,q,S,h_t,product_t,nb_state,P_w,P_m,Q,Vals_old)
                # V doit retourner (valeur, action optimale)
                Vals_new[s, q+1] = val
                Policies_m[s,q+1,Tmax - t + 1,:] = m_opt
                Policies_k[s,q+1,Tmax - t + 1,:] = k_opt
            end
        end
        Vals_all[:, :, t] = Vals_new
        Vals_old = Vals_new
    end

    return Vals_all, Policies_m, Policies_k

end

function simulate_period(Policies_m_t,state,h,nb_state)

    states = [state]

    d = max(state - 1,1)

    u = [0 for i in 1:60]
    q = [0 for i in 1:60]

    q[1] == d*Policies_m_t[1]

    last_state = state

    if (q[1]>0)&&(h[1]==1)
        u[1] == 1
    end

    for t in 1:60-1
        q[t+1] == q[t] - u[t] + d*Policies_m_t[t+1]
        if (q[t+1]>0)&&(h[t+1]==1)
            u[t+1] == 1
        end

        probas = [0.0 for i in 1:nb_state]

        if u[t+1] == 0
            probas = P_w[last_state, :]
        else 
            probas = P_m[last_state, :]
        end
        last_state = sample(1:nb_state, Weights(probas))

        push!(states, last_states)
    end

    return states

end 


function simulate(Policies_m,Policies_k,h,nb_state, years,Q)

    state_list = []

    last_state = 1
    last_q = Q

    for t in 1:(6*years)

        Policies_k_t = Policies_k[last_state,last_q,t,:]
        sum_k = sum(Policies_k_t )

        Policies_m_t = Policies_m[last_state,last_q,t,:]
        h_t = h[t, :, :]

        new_states = simulate_period(Policies_m_t,last_state,h_t,nb_state)

        last_q -= sum_k
        last_state = new_states[end]
        push!(state_list,new_states)

    end

    return state_list

end


#script

#parameters

p = 5.591277606933184e-5
nb_state = 12

Q = 1

S = 3

# wind park always produces the maximum amount of electricity
product_month = [100 for s in 1:S, t in 1:60]

# random wave height, accessible with proba 0.9
d = Bernoulli(0.9)
h_month = [Int(rand(d)) for s in 1:S, t in 1:60]

P_w = zeros(nb_state, nb_state)

for i in 1:(nb_state - 1)
    P_w[i,i] = 1.0 - p 
    P_w[i+1,i] = p
end

P_w[nb_state,nb_state] = 1.0


P_m = zeros(nb_state, nb_state)

for i in 1:(nb_state)
    P_m[1,i] = 1.0 
end




years = 1

h = Array{Float64}(undef, 6*years, S, 60)
for t in 1:(6*years)
    h[t,:,:] = h_month
end

product = Array{Float64}(undef, 6*years, S, 60)
for t in 1:(6*years)
    product[t,:,:] = product_month
end


#Vals_all, Policies_m, Policies_k = Bellman(years,Q,S,h,product,nb_state,P_w,P_m)
#state_list = simulate(Policies_m,Policies_k,h,nb_state, years,Q)