using JuMP
using Gurobi
using Random
using ColorSchemes
using Plots
using ColorTypes

Random.seed!(1234)

include("lecture_distances.jl")
include("creation_arcs.jl")
file = "instances_aerodromes-20230918/instance_70_1.txt"
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
@variable(model,x[i in 1:n,j in 1:n],Bin)


## On peut remplacer y par z comme ils jouent ici le même rôle
#variables bianires, y[i] à 1 si l'aéroport i est visité
#@variable(model,y[i in 1:n],Bin)

#variable correspondant à l'"auxiliary flow" q du "Single-flow formulation" (SF)
@variable(model,q[i in 1:n, j in 1:n] >= 0)

@variable(model,z[k in setdiff(1:n, [d])],Bin)

#fonction objective
@objective(model,Min, sum( sum(D[i,j]*x[i,j] for j in 1:n ) for i in 1:n )  )


### CONTRAINTES

#contrainte sur la distance maximale R
@constraint(model,[i in 1:n, j in 1:n], D[i,j]*x[i,j] <= R)

#contrainte pour nous assurer que y[i] à 1 si un arc arrive sur i
#@constraint(model,[i in 1:n], n*y[i] >= sum(x[j,i] for j in 1:n))

#contrainte pour nous assurer que y[i] à 0 si aucun arc n'arrive sur i
#@constraint(model,[i in 1:n], sum(x[j,i] for j in 1:n) >= y[i] )

#contrainte qui nous assure de visiter Amin aéroport
@constraint(model, sum(z[k] for k in setdiff(1:n, [d])) >= Amin-1)


#contrainte de conservation du flux
for i in 1:n
    if i == d
        @constraint(model,sum(x[i,j] - x[j,i] for j in 1:n) == 1 )
    elseif i == f
        @constraint(model,sum(x[i,j] - x[j,i] for j in 1:n) == -1 )
    else
        @constraint(model,sum(x[i,j] - x[j,i] for j in 1:n) == 0 )
    end
end

#contrainte qui nous assure de parcourir toutes les régions, on s'assure d'avoir de l'activité dans chaque région
#Pour cela on s'assure que pour chaque région, si on regarde les arcs sortants et entrants de tous les noeuds
#qui compose la région, au moins un est activé
@constraint(model,[k in 1:Nr], sum(sum(x[r,j] + x[j,r] for j in 1:n) for r in regions[k]) >= 1 )


#Pour éviter le surplace
@constraint(model,[i in 1:n],x[i,i]==0)

############ Additional constraints corresponding to the Single-flow formulation #############

#Si le flot est positif sur l'arc ij si et seulement si celui-ci appartient au chemin s-t
@constraint(model,[i in 1:n, j in 1:n],q[i,j] <= (n-1)*x[i,j])

#Le flot quittant s est égal au nombre de noeuds qui seront atteints (au sens entrant) par le chemin s-t
@constraint(model,sum(q[d,j] for j in 1:n) == sum(z[k] for k in setdiff(1:n, [d])))

#L'équilibre entre le flot sortant et entrant sur chaque noeud est obtenu grâce à z
@constraint(model,[k in setdiff(1:n, [d])], sum(q[i,k] for i in 1:n) - sum(q[k,j] for j in 1:n) == z[k])

#z marque les noeuds visités par le chemin (c'est y !?)
@constraint(model,[k in setdiff(1:n, [d])], sum(x[i,k] for i in 1:n) == z[k])

#resolution
JuMP.optimize!(model)

#affichage des résultats
obj_value = JuMP.objective_value(model)
println("Objective value : ", obj_value)

#Solution
x_star=JuMP.value.(x)
z_star = JuMP.value.(z)

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
println("Nombre d'aéroports visités : ",1 + sum(z_star[k] for k in setdiff(1:n, [d])))

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