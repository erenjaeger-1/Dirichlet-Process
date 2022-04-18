# # Dirichlet Process #

# The code has been commented, but please read below for a full explanation 
# of the output. 

# Iterations must be increased when we increase the number of mutations 
# This is so we can identify more unique configurations
Dirichlet.Process <- function(alpha, mutations, iterations = 1000){
  run.time = Sys.time()   # Start time 
  configurations = vector(mode = "list", length = iterations)   # Initialise a list to hold the configurations
  for(i in 1:iterations){                                       
    if(mutations == 1){ config.n1 = t(as.matrix(c(1, 1)))       # If n = 1 
    config.n1 = cbind(config.n1, 1)       # There is only 1 configuration so the probability is 1
    colnames(config.n1) = c("Subclone", "Size", "Probability")  # Name columns  
    return(config.n1)  # Return the configuration
    }else{             # If n > 1 #
      category = 1     # Initialise a category with 1 mutation
      for(j in 2:mutations){                # Start from the 2nd mutation because we have allocated 1st. 
        if(alpha/(alpha + j-1) > runif(1)){ # Mutation into new category with prob alpha/(alpha+j-1)
          category = c(category, j)
        }else{  # Else randomly assign the mutation to an existing category 
          category = c(category, sample(unique(category), 1))
        }
      }
     configurations[[i]] = as.data.frame(category)         # Store configurations as df's in the list
     if(i%%500==0) print(paste(i, "Iterations completed")) # Printout to indicate progess so that we know the program is running 
    } 
  }
  
  # Now we calculate the probability of each configuration:
  unique.configs = unique(configurations) # Retrieve the unique configurations 
  print(paste(length(unique.configs), "Unique Configurations"))  # With enough iterations the number of unique
                                                                 # configurations will follow the Bell Numbers
  matrix = matrix(nrow = length(unique.configs), ncol = iterations) # Initialise matrix to hold probabilites
  for(i in 1:length(unique.configs)){    # Count number of times a particular config is
    for(j in 1:iterations){              # repeated and store in matrix
      matrix[i,j] = cbind(sum(identical(unique.configs[[i]], configurations[[j]])))
    }  # Each row corresponds to one unique configuration and the sum of the 
  }    # row is the number of times it is replicated  
  prob = as.matrix(apply(matrix, 1, function(x){sum(x) / iterations}))  # Prob = The sum of each row / number of configurations 
  unique.configs$Probabilities = prob 
  max.prob = unique.configs$Probabilities[prob == max(prob)]  # Select the max probability 
  max.config = unique.configs[which(unique.configs$Probabilities == max(unique.configs$Probabilities))] # Select the corresponding configuration(s) 
  
  # Clean up the output. Information regarding the most likely configuration(s)
  # will be displayed in two data frames:
  max.summary = lapply(max.config, function(x) as.data.frame(table(x)))     # Summarise the max configuration in the form of a table
  max.summary = lapply(max.summary, setNames, c("Subclone", "Size"))        # Name columns
  max.config = lapply(max.config, function(x) cbind(1:mutations, x))        # Extra details of max config. This df tells us the 
                                                                            # category that each mutation was assigned to
  max.config = lapply(max.config, setNames, c("Mutation", "Into Subclone")) # Name columns of this dataframe
  run.time = Sys.time() - run.time                           # Calculate the runtime
  return(list(max.config, max.summary, max.prob, run.time))  # Return results in a list
}



# Run with alpha = 1.5, mutations = 3 and iterations = 1000 

run.1 = Dirichlet.Process(alpha = 1.5, mutations = 3, iterations = 1000)
run.1

# Run with alpha = 1.5, mutations = 10 and iterations = 10,000

run.2 = Dirichlet.Process(alpha = 1.5, mutations = 10, iterations = 10000)
run.2 

# Run with alpha = 1.5, mutations = 30 and iterations = 10,000 
# Note: this takes < 6 minutes to run. The second loop seems to be very slow. 

run.3 = Dirichlet.Process(alpha = 1.5, mutations = 30, iterations = 10000)
run.3 

# Run when mutations = 1

run.4 = Dirichlet.Process(alpha = 1.5, mutations = 1, iterations = 1000)
run.4 


 
#################################
#   OUTPUT WHEN MUTATIONS > 1   #
#################################


### First element of the list: ###

# Dataframe(s) refers to the most probable clustering configuration(s).
# Specifically it tells us the subclone that each mutations was assigned to. 
# The number of rows are therefore equal to the number of mutations. 

# We do this because some unique configurations will look identical when 
# represented in the format of a table. An example of this is given by 2 of the 
# configurations that result from 4 mutations. Consider the configuration 1,2,1,2 and 1,2,2,1
# where each number is a mutation and their value refers to the subclone they were assigned
# Specifically, the first sequence tells us that the 1st and 3rd mutations were assigned to the
# same subclone. In contrast, the second sequence tells us that the 1st and 4th mutations were 
# assigned to the same subclone. This information is lost if we convert to a table.

# Example:
# These are all the configurations possible when we have 4 mutations and 2 categories
# a and b look identical when they have been converted to tables. 

a = data.frame(x = c(1,2,2,1))  
b = data.frame(x = c(1,2,1,2))
c = data.frame(x = c(1,1,3,3))
as.data.frame(table(a))  
as.data.frame(table(b))
as.data.frame(table(c))
a
b
# Not possible to distinguish between a and b even though the mutations in each  
# subclone are different. 


### Second element of the list: ###

# Dataframe(s) summarises the dataframe(s) above in the form of tables. 
# The number of rows is equal to the number of subclones. The number of each
# subclone refers to the mutation that first gave rise to it. For example, subclone 7
# is named this because mutation 7 was first assigned to it. The size column then refers
# to the number of mutations that were assigned to this subclone. We include this extra 
# dataframe(s) because it is easier to see the output when we have more mutations. 


### Third element of the list: ###

# Probability of the most likely clustering solution(s).


### Fourth element of the list: ###

# Runtime 















