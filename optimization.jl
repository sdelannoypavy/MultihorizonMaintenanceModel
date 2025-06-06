
using Gurobi
using JuMP  
using Distributions

function V(p,S,nb_state)

    V_table = [0.0 for _ in 1:nb_state]

    #p = 0.0005

    #P_w = [[1.0 - p, 0.0],
    #       [p, 1.0]] 

    P_w = zeros(nb_state, nb_state)

    for i in 1:(nb_state - 2)
        P_w[i,i] = 1.0 - 2*p 
        P_w[i+1,i] = p
        P_w[nb_state,i] = p
    end

    P_w[nb_state - 1,nb_state - 1] = 1.0 - p 
    P_w[nb_state, nb_state - 1] = p

    P_w[nb_state,nb_state] = 1.0

    #P_m = [[1.0, 1.0],
    #       [0.0, 0.0]] 

    P_m = zeros(nb_state, nb_state)

    for i in 1:(nb_state)
        P_m[1,i] = 1.0 
    end

    x_0 = [1.0; zeros(nb_state - 1)]

    C_x = [100.0 - (i-1)*100/(nb_state-1) for i in 1:nb_state]

    #nb_state = 2

    T = 30
    #S = 3

    d_c = 5

    prod = [100 for s in 1:S, i in 1:30]

    d = Bernoulli(0.9)
    h = [Int(rand(d)) for s in 1:S, i in 1:30]

    K = 2

    model = Model(Gurobi.Optimizer)

    @variable(model, x[1:(S+1),1:(T+1),1:nb_state] >= 0)

    @variable(model, c[1:T] >= 0)
    @variable(model, p_evac[1:S,1:T,1:nb_state] >= 0)

    @variable(model, m[1:T], Bin)
    @variable(model, u[1:(S+1),1:T], Bin)
    @variable(model, ong_m[1:(S+1),1:T], Bin) #ongoing maintenance
    @variable(model, q[1:(S+1),1:T] >= 0)
    @variable(model, k[1:T], Bin)

    @objective(model, Min, sum(c[t] for t in 1:T) + sum(x[S+1,T+1,i]*V_table[i] for i in 1:nb_state))


    @constraint(model, [s in 1:(S+1), t in 2:(T+1), i in 1:nb_state], x[s,t,i] == (1-u[s,t-1])*sum(P_w[i,j]*x[s,t-1,j] for j in 1:nb_state) + u[s,t-1]*sum(P_m[i,j]*x[s,t-1,j] for j in 1:nb_state))

    @constraint(model, [s in 1:(S+1), i in 1:nb_state], q[s,1] == d_c*m[1])
    @constraint(model, [s in 1:(S+1), t in 1:(T-1)], q[s,t+1] == q[s,t] - u[s,t] + d_c*m[t+1])
    @constraint(model, [s in 1:S, t in 1:(T)], u[s,t] + 1 - h[s,t] >= (1/T)*q[s,t])
    @constraint(model, [t in 1:(T)], u[S+1,t] >= (1/T)*q[S+1,t])
    @constraint(model, [s in 1:(S+1), t in 1:(T)], u[s,t] <= q[s,t])
    @constraint(model, [s in 1:S, t in 1:(T)], u[s,t] <= h[s,t])
    @constraint(model, [s in 1:(S+1), t in 1:(T)], ong_m[s,t] >= (1/T)*q[s,t])

    @constraint(model, [s in 1:(S+1), i in 1:nb_state], x[s,1,i] == x_0[i])

    @constraint(model, sum(k[t] for t in 1:T) <= K)

    @constraint(model, [s in 1:S, t in 1:T, i in 1:nb_state], p_evac[s,t,i] >= prod[s,t] - C_x[i])
    @constraint(model, [t in 1:T], c[t] >= (1/S)*sum((1 - u[s,t])*sum(x[s,t,i]*p_evac[s,t,1] for i in 1:nb_state) + ong_m[s,t]*(1-k[t])*prod[s,t] for s in 1:S))



    optimize!(model)

    status = termination_status(model)
    println("Statut de l'optimisation: $status")

    cost = objective_value(model)

    numvars = MOI.get(model, Gurobi.ModelAttribute("NumVars"))
    numcons = MOI.get(model, Gurobi.ModelAttribute("NumConstrs"))


    if status == MOI.OPTIMAL
        #m = [value(m[t]) for  t in 1:T]
        #u = [value(u[s,t]) for  s in 1:S, t in 1:T]
        #q = [value(q[s,t]) for  s in 1:S, t in 1:T]
        #ong_m = [value(ong_m[s,t]) for  s in 1:S, t in 1:T]
        #x = [value(x[s,t,nb_states]) for  s in 1:S, t in 1:(T+1)]

        return(numvars,numcons)
    else
        println("Aucune solution optimale trouvée.")
    end

end
