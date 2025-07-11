MTBF = 588
# using Polynomials

psim = (12/MTBF)/365

function proba(p,i,j)
    # j > i
    return binomial(1680-i+1,j-i)*(p^(j-i))*((1-p)^(1680-j+1))
end 

function proba_ii(p,i)
    p_ii = 1
    for j in (i+1):12
        p_ii -= proba(p,i,j)
    end
    return p_ii
end

function compute_MTBF(p)

    # T = vector of time before absorption
    T = [0.0 for i in 1:12]

    for i in 11:-1:1
        T[i] = (1/(1 - proba_ii(p,i)))*(1 + sum(proba(p,i,j)*T[j] for j in (i+1):12))
    end

    return T[1]

end



P = vcat([100.0 for i in 1:11], [0.0])

#binomial(1677, 8) overflows
# proba(1,3) is 10 times smaller than proba(1,2)
# proba(1,2) is very close to proba(11,12)

#CCl: we take p(i,i+1) = psim for each i, p(i,j) = 0 if j > i+1