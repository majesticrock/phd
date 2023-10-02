from itertools import product

def iterate_containers(*containers):
    # Generate combinations of indices for each container
    index_combinations = product(*[range(len(container)) for container in containers])
    
    for indices in index_combinations:
        values = [container[index] for container, index in zip(containers, indices)]
        yield values
        
def total_size(*containers):
    size = 1
    for container in containers:
        size *= len(container)
    return size
    
## Example usage
#import numpy as np
#
#Ts = np.array([1, 2, 3])
#Us = np.array([4, 5])
#Vs = np.array([6, 7, 8])
#
#for T, U, V in iterate_containers(Ts, Us, Vs):
#    # Do something with T, U, V
#    print(f"T: {T}, U: {U}, V: {V}")
#
#print(total_size(Ts, Us, Vs))