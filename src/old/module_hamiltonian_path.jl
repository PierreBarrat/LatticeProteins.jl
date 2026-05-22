module HamiltonianPath
    using StaticArrays
    using DataStructures
    using CSV
    using DataFrames

    # ------------------------------ FUNCTIONS ------------------------------
    #ID (1-27) to Coordinates (0-2, 0-2, 0-2)
    function id_to_coords(id)
        z = div(id - 1, 9) #id - 1 because Julia indices start at 1
        rem = mod(id - 1, 9) #remainder of the entire division
        y = div(rem, 3)
        x = mod(rem, 3)
        return (x, y, z)
    end

    #Coord to ID
    function coords_to_id(x, y, z)
        return x + 3y + 9z + 1
    end

    #Generating the table of 24 rotations (only det = 1)
    #Reflections are excluded when counting chiral forms separately
    function generate_rotation_table()
        #permutation table : table[id, r] = new_id -> under rotation r, node id becomes new_id
        table = zeros(Int, 27, 48) #27 -> each node ; 24 -> each rotation
        
        #6 possibles transposition : [1,2,3] -> (x,y,z)  ;  [1,3,2] -> (x,z,y)  ;  [2,1,3] -> (y,x,z)  etc.
        perms = [[1,2,3], [1,3,2], [2,1,3], [2,3,1], [3,1,2], [3,2,1]]
        
        #Each combination corresponds to a change of sign on the axes : total 2^3 = 8 sign combination
        signs = [[s1, s2, s3] for s1 in [-1,1], s2 in [-1,1], s3 in [-1,1]]
        
        #so we got 6*8 = 48 rotations
        #we only keep 24 rotations because we want to keep mirror inversions (so mirror are not considered as rotations)

        rot_idx = 1
        for p in perms, s in signs
            #Calculating the determinant to retain only the physical rotations
            
            #A transposition is even if the number of permutation is even, determinant of even=1, else -1
            p_parity = (p == [1,2,3] || p == [2,3,1] || p == [3,1,2]) ? 1 : -1
            
            #The transformation is a permutation of the axes followed by sign reversals
            #Its determinant is det(P)*det(D) = parity(permutation) × (s1*s2*s3)
            #This comes from the general property: det(AB) = det(A)*det(B)
            #Since the sign matrix is ​​diagonal, its determinant is the product of its diagonal
            det = p_parity * (s[1] * s[2] * s[3])

            for id in 1:27 #rotation applied to all points
                c = id_to_coords(id) .- 1 #Centering on 0 for rotation, (0,1,2) becomes (-1,0,1)
                new_c = (c[p[1]]*s[1], c[p[2]]*s[2], c[p[3]]*s[3]) #application of the permutation
                table[id, rot_idx] = coords_to_id(new_c[1]+1, new_c[2]+1, new_c[3]+1)
            end
            rot_idx += 1 #next rotation
        end
        return table
    end

    #Canonical formatting function
    function get_canonical(path)
        best_path = copy(path)
        #similar(path) create an object of same type and same lenght but empty
        current_variant = similar(path) #to temporarily store an applied rotation
        rev_variant = similar(path) #to temporarily store the reversed version of the path
        
        for s in 1:48 #each rotation
            #application of rotation
            for i in 1:27 #each node
                current_variant[i] = ROT_TABLE[path[i], s] #new id after rotation
            end
            
            #We are looking for the minimal form among all the rotations
            if current_variant < best_path #lexicographical comparison ; [1,3,9] < [1,4,2] because 1==1 then 3<4 
                best_path .= current_variant #copy element by element (.)
            end
        end
        return best_path
    end

    #Connectivity graph
    function build_graph()
        adj = [Int[] for _ in 1:27] #init of 27 empty vectors
        for id in 1:27
            x, y, z = id_to_coords(id)
            #we consider the 6 directions in the cube (+x,-x,+y,-y,+z,-z)
            for (dx, dy, dz) in [(1,0,0), (-1,0,0), (0,1,0), (0,-1,0), (0,0,1), (0,0,-1)]
                nx, ny, nz = x+dx, y+dy, z+dz
                if 0 <= nx < 3 && 0 <= ny < 3 && 0 <= nz < 3 #if the neighbor exist
                    push!(adj[id], coords_to_id(nx, ny, nz)) #add the neighbor to the vector
                end
            end
        end
        return adj
    end

    #Backtracking with filtering by Set
    function find_all_unique_conformations()
        #the Set guarantees uniqueness
        unique_paths = OrderedSet{Vector{Int8}}()
        visited = zeros(Bool, 27)
        current_path = zeros(Int8, 27)

        #Recursive function to explore all Hamiltonian paths
        function backtrack(node, depth)
            visited[node] = true #indicates that this node is explored
            current_path[depth] = node #we add the node to the current path
            
            if depth == 27
                #we transform the path into its unique reference version
                push!(unique_paths, get_canonical(current_path)) #if already in Set, not added
            else
                for neighbor in ADJ[node]
                    if !visited[neighbor] #if this neighbor has not been visited
                        backtrack(neighbor, depth + 1) #backtrack with this new node
                    end
                end
            end
            #after exploring all paths starting from this node, it is freed up for other paths
            visited[node] = false 
        end
        println("Starting backtracking")
        println("Calculation of all conformations")
        
        println("Step 1: Exploration from the corner (node 1)")
        backtrack(1, 1) 
        
        println("Step 2: Exploration from the center of the face (node 5)")
        backtrack(5, 1) 
        
        return unique_paths, length(unique_paths)
    end

    # Function for saving paths 
    function paths_to_file(filename)
        """
        Save the paths to a CSV file, each path is a row and nodes are separated by commas.
        """
        open(filename, "w") do file
            for path in PATHS
                println(file, join(path, ",")) 
            end
        end
    end


    #function for keeping only a part of the data (to test quickly for example)
    #percentage needs to be superior to 0 and max 100
    function reduce_data(paths, percentage=10)
        @assert 0< percentage <=100 "Value of percentage needs to be superior to 0 and max 100"
        step=100.0 / percentage
        result=OrderedSet{Vector{Int8}}()
        next=1.0

        for i in 1:length(paths)
            if i>= next
                push!(result, paths[i])
                next += step
            end
        end
        return result
    end

    # ------------------------------ Variables ------------------------------
    # cube 3*3*3
    const N = 3
    const TOTAL_NODES = N^3
    const ROT_TABLE = generate_rotation_table()
    const ADJ = build_graph()
    paths, _ = find_all_unique_conformations()
    const PATHS = paths


    # ------------------------------ Exports ------------------------------
    #Functions
    export id_to_coords, coords_to_id, generate_rotation_table, get_canonical, build_graph, find_all_unique_conformations, reduce_data, paths_to_file
    #Global variables
    export ADJ, PATHS, ROT_TABLE

end #module HamiltonianPath