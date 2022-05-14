import numpy as np
import random

def initial_pi_state(pi_matrix):
    pi_index = get_random_interval(pi_matrix)
    return(pi_index)

def get_new_symbol_index(state, B_matrix):
    symbol_index = get_random_interval(B_matrix[state])
    return(symbol_index)

def get_next_state(state, A_matrix):
    new_state = get_random_interval(A_matrix[state][:])
    return (new_state)

# Make a list of intervals up to a 1000 representing the probabilities in the list
# For instance p = [0.1,0.1,0.8] gives the intervals [0,100,200,1000]
# Randomly choose a number between 1 and 1000
# Return the index of the interval within which this random number falls

def get_random_interval(list_of_probabilities):
    number_of_intervals = len(list_of_probabilities)
    factor = 1000
    interval_list = [0]
    for value in list_of_probabilities:
        interval_list.append(int(interval_list[-1] + factor * value))
    random.seed()
    random_int = random.randint(1, factor)
    new_index = np.searchsorted(interval_list, random_int, side='left') - 1
    new_index = new_index if new_index < number_of_intervals else number_of_intervals-1
    return(new_index)


pi_matrix = [0.5,0.5]
A_matrix = [[0.95,0.05],[0.1,0.9]]
B_matrix = [[1/6,1/6,1/6,1/6,1/6,1/6],[0.1,0.1,0.1,0.1,0.1,0.5]]
list_of_symbols = [1,2,3,4,5,6]
list_of_states = ['F','L']

# initialize
symbol_index_pattern = []
state_index = []
state = initial_pi_state(pi_matrix)
symbol_index = get_new_symbol_index(state, B_matrix)
symbol_index_pattern.append(symbol_index)
state_index.append(state)
for i in range(100):
    state = get_next_state(state, A_matrix)
    state_index.append(state)
    symbol_index = get_new_symbol_index(state, B_matrix)
    symbol_index_pattern.append(symbol_index)
for i in range(100):
    print(list_of_states[state_index[i]],end="")
print()
for i in range(100):
    print(symbol_index_pattern[i],end="")
print()
for i in range(100):
    print(list_of_symbols[symbol_index_pattern[i]],end="")

