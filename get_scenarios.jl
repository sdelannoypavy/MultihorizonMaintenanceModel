using DelimitedFiles

scenarios_h = Array{Float64}(undef, 10, 6, 62)
scenarios_p = Array{Float64}(undef, 10, 6, 62)

for j in 1:6
    # Lire la matrice 10x62 depuis le fichier
    slice_h = readdlm("scenarios_h$(j).txt",Int)
    # Remettre dans le bon endroit
    scenarios_h[:, j, :] = slice_h
    for j in 1:6
    # Lire la matrice 10x62 depuis le fichier
    slice_p = readdlm("scenarios_p$(j).txt")
    # Remettre dans le bon endroit
    scenarios_p[:, j, :] = slice_p
end
end

h = Array{Float64}(undef, 30*6, 10, 62)
p = Array{Float64}(undef, 30*6, 10, 62)

for T in 1:(30*6)

    for s in 1:10
        type = mod(T,6)
        if mod(T,6) == 0
            type = 6
        end
        h[T,s,:] = scenarios_h[s,type,:] # we use the same scenarios for each jan/fev strategic period
        p[T,s,:] = scenarios_p[s,type,:]
    end

end