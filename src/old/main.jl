include("hamiltonian_path.jl")
include("visualization.jl")

using .HamiltonianPath: PATHS, paths_to_file 
using .Visualization: contact_map, plot_cube_ctc_map, x, y, z, lx, ly, lz


# Finding all paths 
println("We find $length(PATHS) unique paths.")
println(first(PATHS))
println(typeof(PATHS))

# Saving paths to a file
paths_to_file("paths.csv")
println("Paths saved to paths.csv")

# Finding the contact map of a given path
println("The first path is : \n", first(PATHS))
println("Its contact map is : \n", contact_map(first(PATHS)))

# Visualizing the first path
println("Visualizing the first path")
liste_contact = contact_map(first(PATHS))
display(plot_cube_ctc_map(x, y, z, lx, ly, lz, liste_contact))