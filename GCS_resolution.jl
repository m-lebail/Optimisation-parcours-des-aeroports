using JuMP
using Gurobi
using Random
using Plots


Random.seed!(1234)

include("lecture_distances.jl")

file = "instances_aerodromes-20230918/instance_20_1.txt"
#Get instance data
n,d,f,Amin,Nr,R,regions,coords,D = readInstance(file)

#D = init_arc_matrix(coords) 

# #Is t in R_i with i in {1,k} ?
# i_t = 0
# for (key,value) in regions
#     if f in value
#         i_t = key
#         break
#     end
# end

model = JuMP.Model(Gurobi.Optimizer)

#variables
#@variable(model,x[i in 1:n,j in 1:n],Bin)

#On se contente de la relaxation continue du problème au départ, car on cherche seulement
#à déterminer les contraintes de soustours à ajouter pour obtenir notre solution
@variable(model,x[i in 1:n,j in 1:n] >= 0)

#variables bianires, y[i] à 1 si l'aéroport i est visité
#@variable(model,y[i in 1:n],Bin)
@variable(model,y[i in 1:n] >= 0)

#fonction objective
@objective(model,Min, sum( sum(D[i,j]*x[i,j] for j in 1:n ) for i in 1:n )  )


### CONTRAINTES

#contrainte sur la distance maximale R
@constraint(model,[i in 1:n, j in 1:n], D[i,j]*x[i,j] <= R)

#contrainte pour nous assurer que y[i] à 1 si un arc arrive sur i
@constraint(model,[i in 1:n], n*y[i] >= sum(x[j,i] for j in 1:n))

#contrainte pour nous assurer que y[i] à 0 si aucun arc n'arrive sur i
@constraint(model,[i in 1:n], sum(x[j,i] for j in 1:n) >= y[i] )

#contrainte qui nous assure de visiter Amin aéroport
@constraint(model, sum(y[i] for i in 1:n) - y[d] >= Amin-1)

#contrainte de conservation du flux
for i in 1:n
    if i == d
        @constraint(model,sum(x[i,j] for j in 1:n) - sum(x[j,i] for j in 1:n) == 1 )
    elseif i == f
        @constraint(model,sum(x[i,j] for j in 1:n) - sum(x[j,i] for j in 1:n) == -1 )
    else
        @constraint(model,sum(x[i,j] for j in 1:n) - sum(x[j,i] for j in 1:n) == 0 )
    end
end

#contrainte qui nous assure de parcourir toutes les régions, on s'assure d'avoir de l'activité dans chaque région
#Pour cela on s'assure que pour chaque région, si on regarde les arcs sortants et entrants de tous les noeuds
#qui compose la région, au moins un est activé
@constraint(model,[k in 1:Nr], sum(x[ip,j] for ip in regions[k] for j in 1:n) + sum(x[i,jp] for jp in regions[k] for i in 1:n) >= 1 )

#Pour éviter le surplace
@constraint(model,[i in 1:n],x[i,i]==0)


#resolution
JuMP.optimize!(model)

#affichage des résultats
obj_value_chemin = JuMP.objective_value(model)
println("##############################################################################################")
println("Valeur objectif pour le chemin entre les aéroports : ", obj_value_chemin)
println("##############################################################################################")

#Solution
x_star=JuMP.value.(x)
y_star = JuMP.value.(y)


#println("x_star : ",x_star)
#println("y_star :" ,y_star)

#for n_iter in 1:100 

##Deuxième problème d'optimisation pour la violation des CONTRAINTES

model_contrainte = JuMP.Model(Gurobi.Optimizer)

@variable(model_contrainte,a[i in 1:n],Bin)
@variable(model_contrainte,b[j in 1:n],Bin)
@variable(model_contrainte,c[i in 1:n,j in 1:n] >= 0)


#fonction objective
@objective(model_contrainte,Max, sum( sum(b[i] * x_star[i,j] for j in 1:n ) for i in 1:n ) - sum( x_star[i,j] * c[i,j] for i in 1:n for j in 1:n ) )


##CONTRAINTES
@constraint(model_contrainte,a[d]==0)
@constraint(model_contrainte,a[f]==0)

@constraint(model_contrainte,sum(a[i] for i in 1:n) >= 2 )
@constraint(model_contrainte,sum(a[i] for i in 1:n) <= n-2)

@constraint(model_contrainte,sum(b[i] for i in 1:n) == 1)
@constraint(model_contrainte, [i in 1:n],b[i] <= a[i])


@constraint(model_contrainte,[i in 1:n, j in 1:n],c[i,j] <= a[i])
@constraint(model_contrainte,[i in 1:n, j in 1:n],c[i,j] <= 1-a[j])

@constraint(model_contrainte,[i in 1:n, j in 1:n],c[i,j] >= a[i]-a[j])


#resolution
JuMP.optimize!(model_contrainte)

#affichage des résultats
obj_value_soustours = JuMP.objective_value(model_contrainte)
println("##############################################################################################")
println("Valeur objectif pour la violation de la contrainte des sous-tours : ", obj_value_soustours)
println("##############################################################################################")

#Solution
S_star=JuMP.value.(a)
k_star = JuMP.value.(b)


#println("S : ",S_star)
#println("k : ",k_star)

S_indices = []
not_S_indices = []

iter = 1
while obj_value_soustours>0 && iter < 100
    #On récupère les numéros des noeuds qui composent S
    global S_indices = findall(x -> x == 1.0, S_star)
    global not_S_indices = findall(x -> x == 0.0, S_star)
   
    #On ajoute la contrainte sur le sous-ensemble S que l'on a identifié
    @constraint(model,[k in S_indices], sum( x[i,j] for i in S_indices for j in not_S_indices) >= sum( x[k,j] for j in 1:n ) )
    
    #resolution
    JuMP.optimize!(model)

    #affichage des résultats
    global obj_value_chemin = JuMP.objective_value(model)
    println("##############################################################################################")
    println("Valeur objectif pour le chemin entre les aéroports : ", obj_value_chemin)
    println("##############################################################################################")


    #Solution
    global x_star=JuMP.value.(x)
    global y_star = JuMP.value.(y)

    # println("x_star : ",x_star)
    # println("y_star :" ,y_star)

    #On met à jour la fonction objectif de notre sous-problème
    @objective(model_contrainte,Max, sum( sum(b[i] * x_star[i,j] for j in 1:n ) for i in 1:n ) - sum( x_star[i,j] * c[i,j] for i in 1:n for j in 1:n ) )

    #resolution
    JuMP.optimize!(model_contrainte)

    #affichage des résultats
    global obj_value_soustours = JuMP.objective_value(model_contrainte)
    println("##############################################################################################")
    println("Valeur objectif pour la violation de la contrainte des sous-tours : ", obj_value_soustours)
    println("##############################################################################################")

    #Solution
    global S_star=JuMP.value.(a)
    global k_star = JuMP.value.(b)


    # println("S : ",S_star)
    # println("k : ",k_star)
    global iter = iter + 1

end

println("")
println("####################################################################")
println("Nombre de contraintes soustrours ajoutées : ", iter-1)
println("####################################################################")
println("")

#On peut maintenant résoudre finalement notre problème de départ à variables entières
for i in 1:n, j in 1:n
    set_binary(x[i, j])
end

for i in 1:n
    set_binary(y[i])
end

#resolution
JuMP.optimize!(model)

#affichage des résultats
obj_value_chemin = JuMP.objective_value(model)
#affichage des résultats
println("##############################################################################################")
println("Valeur objectif pour le chemin entre les aéroports : ", obj_value_chemin)
println("##############################################################################################")


x_star=JuMP.value.(x)
y_star = JuMP.value.(y)


indice=d
print(d)
noeuds=1
chemin = []
while indice!=f
    for j in 1:n
        if x_star[indice,j]==1
            print("->",j)
            global indice=j
            global noeuds+=1
            push!(chemin,indice)
        end
    end
end 
verif=sum(x_star[i,j] for i in 1:n, j in 1:n)
println("\nNombre d'arcs activés : ", verif)
println("Nombres de noeuds visités : ",noeuds)
println("Nombre d'aéroports sur lesquels on atterit : ",1 + sum(y_star[k] for k in 1:n))



x = coords[:,1]
y = coords[:,2]

#arc_color = :red

# # Define the markers for nodes d and f
# d_marker = :diamond
# f_marker = :square



# Number of distinct colors you want
num_colors = length(regions)

# Define a color space, you can choose any color space you like
color_space = LCHab

# Create an array to store the generated colors
palette = ColorTypes.Color[]

# Generate distinct colors from a perceptually uniform color map
color_map = distinguishable_colors(num_colors, colorant"blue")
    
for i in 1:num_colors
    push!(palette, color_map[i])
end

# Create a color palette with unique colors for each list in L
#palette = [:red, :green, :blue, :orange, :purple, :yellow, :cyan, :magenta]

def = :grays

plt = scatter(x, y, label="Aeroports",legend=false,color=:grays)


for (i, points) in enumerate(regions)
    # Assign a unique color to this list of points
    color = palette[i]
    #println(points[2])
    # Scatter plot for the points in this list with the assigned color
    #println([x[p] for p in points[2]], [y[p] for p in points[2]])
    scatter!([x[p] for p in points[2]], [y[p] for p in points[2]], label="List $(i)", markercolor=color)

end

# Plot the arcs based on x_ij matrix
for i = 1:n
    for j = 1:n
        if x_star[i, j] == 1
            plot!([x[i], x[j]], [y[i], y[j]], label="Arc $(i) to $(j)") #linecolor=arc_color)
        end
    end
end

#Start and end points
annotate!([(x[d], y[d], text("O", 20, :yellow))])
annotate!([(x[f], y[f], text("O", 20, :rainbow))])

xlabel!("X-axis")
ylabel!("Y-axis")
title!("Carte des aéroports")


display(plt)


