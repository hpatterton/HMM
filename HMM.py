import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)

def gamma(forward_matrix, backward_matrix, number_of_states, number_of_observation):
    gamma = np.zeros((number_of_states,number_of_observation),dtype=float)
    product_alpha_beta = np.zeros((number_of_states,number_of_observation),dtype=float)
    product_alpha_beta = forward_matrix * backward_matrix
    for t in range(number_of_observation):
        sum = np.sum(product_alpha_beta[:,t])
        for i in range(number_of_states):
            gamma[i,t] = product_alpha_beta[i,t]/sum
    return(gamma)

def forward(pi_matrix, a_matrix, b_matrix, pattern_list):
    number_of_states = len(a_matrix)
    length = len(pattern_list)
    alpha_matrix = np.zeros(number_of_states, dtype=float)
    temp_alpha_matrix = np.zeros(number_of_states, dtype=float)
    alpha_matrix = np.copy(pi_matrix)
    alpha_results = np.zeros((number_of_states, length), dtype=float)
    for t in range(length):
        for j in range(number_of_states):
            if (t == 0):
                temp_alpha_matrix[j] = alpha_matrix[j] * b_matrix[j, pattern_list[t]]
            else:
                temp_alpha_matrix[j] = np.dot(alpha_matrix, a_matrix[:, j]) * b_matrix[j, pattern_list[t]]
            alpha_results[j, t] = temp_alpha_matrix[j]
        alpha_matrix = np.copy(temp_alpha_matrix)
    return (alpha_results)

def forward_probability(pi_matrix, a_matrix, b_matrix, pattern_list):
    alpha_results = forward(pi_matrix, a_matrix, b_matrix, pattern_list)
    # Currently not returning 'probability_of_hmm' from this function. You can if you are interested in it.
    probability_of_hmm = np.sum(alpha_results[:,len(pattern_list)-1])
    return probability_of_hmm

def backward(pi_matrix, a_matrix, b_matrix, pattern_list):
    number_of_different_symbols = len(b_matrix[0])
    number_of_states = len(a_matrix)
    length = len(pattern_list)
    beta_matrix = np.ones((number_of_states, 1), dtype=float) # these are the beta(t=T) = 1 values
    temp_beta_matrix = np.zeros((number_of_states, 1), dtype=float)
    beta_results = np.ones((number_of_states, length), dtype=float)
    # Start from T-1. Beta(T-1) = a_ij*beta(T)*b_j(k)
    for i in range(length-2, -1, -1):  # start
        # Python returns a row vector from b_matrix[:,x], so we have to reshape it to a column vector
        # otherwise the multiplication does not multiply corresponding elements
        # We first multiply the beta(t+1) values with corresponding b_j(k) values
        beta_matrix = beta_matrix * b_matrix[:, pattern_list[i + 1]].reshape(number_of_states,1)
        for j in range(number_of_states):
            # Now we take the dot product of the a_ij row with the beta_matrix column
            temp_beta_matrix[j] = np.dot(a_matrix[j, :], beta_matrix[:])
            beta_results[j, i] = temp_beta_matrix[j]
        beta_matrix = np.copy(temp_beta_matrix)
    return (beta_results)

def backward_probability(pi_matrix, a_matrix, b_matrix, pattern_list):
    beta_results = backward(pi_matrix, a_matrix, b_matrix, pattern_list)
    # When done with t=1, multiply each corresponding beta_j(t=1) with b_j(k) and pi_j
    # These are all row vectors, and thus the same shape.  The sum gives beta_zero(*)
    # beta_zero(*) is the probability that the HMM lambda (A,B,pi) generated the
    # observed sequence of symbols
    beta_zero = np.sum(pi_matrix[:]*beta_results[:,0]*b_matrix[:,pattern_list[0]])
    return(beta_zero)

def viterbi(pi_matrix, a_matrix, b_matrix, pattern_list):
    number_of_states = len(a_matrix)
    length = len(pattern_list)
    delta_matrix = np.zeros((number_of_states, length), dtype=float)
    temp_delta_matrix = np.zeros(number_of_states, dtype=float)
    psi_matrix = np.zeros((number_of_states, length), dtype=int)
    path_matrix = np.zeros((length), dtype=int)
    for t in range(length):
        for current_state in range(number_of_states):
            for previous_state in range(number_of_states):
                if (t == 0):  # handle t=1 use pi_matrix
                    temp_delta_matrix[previous_state] = pi_matrix[previous_state] * b_matrix[
                        previous_state, pattern_list[t]]
                else:
                    temp_delta_matrix[previous_state] = delta_matrix[previous_state, t - 1] * a_matrix[
                        previous_state, current_state] * b_matrix[current_state, pattern_list[t]]
            if (t == 0):
                delta_matrix[current_state, t] = temp_delta_matrix[current_state]
            else:
                delta_matrix[current_state, t] = np.max(temp_delta_matrix)
            psi_matrix[current_state, t] = np.argmax(temp_delta_matrix)  # returns the index of the maximum
    # starting at the max value in the last column, trace the route back
    path_matrix[length - 1] = np.argmax(delta_matrix[:, length - 1])
    for position in range(length - 1, 0, -1):
        path_matrix[position - 1] = psi_matrix[path_matrix[position], position]
    return (path_matrix)

def xi(forward_matrix, backward_matrix, a_matrix, b_matrix, number_of_states, number_of_observations):
    # The transition_i_to_j matrix is 3D. i x j x t
    transition_i_to_j = np.zeros((number_of_observations, number_of_states, number_of_states), dtype=float)
    #
    sum_of_transitions = np.zeros((number_of_observations), dtype=float)
    # The xi matrix is 3D. i x j x t
    xi_matrix = np.zeros((number_of_observations, number_of_states, number_of_states), dtype=float)
    for t in range(number_of_observations-1):
        for i in range(number_of_states):
            for j in range(number_of_states):
                transition_i_to_j[t,i,j] = forward_matrix[i,t]*a_matrix[i,j]*b_matrix[j,pattern_list[t+1]]*backward_matrix[j,t+1]
    # Now add up all the transition frequencies at each time t
    for t in range(number_of_observations-1):
        sum = 0
        for i in range(number_of_states):
            for j in range(number_of_states):
                sum += transition_i_to_j[t,i,j,]
        sum_of_transitions[t] = sum
    # Normalise with sum and write to xi_matrix
    for t in range(number_of_observations-1):
        for i in range(number_of_states):
            for j in range(number_of_states):
                xi_matrix[t,i,j] = transition_i_to_j[t,i,j]/sum_of_transitions[t]
    return(xi_matrix)

def new_A_matrix(xi_matrix, gamma_matrix, number_of_states, number_of_symbols):
    new_A_matrix = np.zeros((number_of_states,number_of_states), dtype=float)
    for i in range(number_of_states):
        number_of_times_in_state_i = np.sum(gamma_matrix[i,:])
        for j in range(number_of_states):
            number_of_transitions_from_i_to_j = 0
            for t in range(number_of_symbols-1):
                number_of_transitions_from_i_to_j += xi_matrix[t,i,j]
            new_A_matrix[i,j] = number_of_transitions_from_i_to_j/number_of_times_in_state_i
    return(new_A_matrix)

def new_pi_matrix(gamma_matrix):
    new_pi = gamma_matrix[:,0]
    return(new_pi)

def get_symbol_pattern_indexes(series_of_symbols, symbol):
    begin = 0
    index = 0
    list_of_indexes = []
    while (index != -1 and begin < len(series_of_symbols)):
        try:
            index = series_of_symbols[begin:].index(symbol)
        except ValueError:
            index = -1
        if (index >= 0):
            list_of_indexes.append(index + begin)
            begin = begin + index + 1
    return(list_of_indexes)

def new_B_matrix(gamma_matrix, number_of_states, pattern_of_symbols):
    set_of_symbols = sorted(set(pattern_of_symbols)) # returns a list where unique symbols are sorted
    number_of_different_symbols = len(set_of_symbols)
    new_B_matrix = np.zeros((number_of_states,number_of_different_symbols),dtype=float)
    for symbol_index in range(number_of_different_symbols):
        symbol_indexes_in_pattern = get_symbol_pattern_indexes(pattern_of_symbols, set_of_symbols[symbol_index])
        for state in range(number_of_states):
            symbol_observed_in_state = 0
            for t in symbol_indexes_in_pattern:
                symbol_observed_in_state += gamma_matrix[state,t]
            number_of_times_in_state = np.sum(gamma_matrix[state,:])
            new_B_matrix[state,symbol_index] = symbol_observed_in_state/number_of_times_in_state
    return(new_B_matrix)

def baum_welch(pi_matrix, a_matrix, b_matrix, pattern_list):
    number_of_states = len(pi_matrix)
    number_of_observations = len(pattern_list)
    forward_result = forward(pi_matrix, a_matrix, b_matrix, pattern_list)
    backward_result = backward(pi_matrix, a_matrix, b_matrix, pattern_list)
    gamma_result = gamma(forward_result, backward_result, number_of_states, number_of_observations)
    xi_result = xi(forward_result, backward_result, a_matrix, b_matrix, number_of_states, number_of_observations)
    new_pi = new_pi_matrix(gamma_result)
    new_A = new_A_matrix(xi_result, gamma_result, number_of_states, number_of_observations)
    new_B = new_B_matrix(gamma_result, number_of_states, pattern_list)
    return(new_pi, new_A, new_B)


pi_matrix = np.array([0.4,0.3,0.3],float)
A_matrix = np.array([[0.1,0.4,0.5],[0.3,0.4,0.3],[0.2,0.3,0.5]],float)
B_matrix = np.array([[0.5,0.5],[0.5,0.5],[0.5,0.5]],float)
pattern_list = [0,0,0,0,0,1,1,1,1,1]
number_of_states = len(pi_matrix)
number_of_observations = len(pattern_list)


# hmm_probabilities =[]
# x_values = [0]
# hmm_probability=forward_probability(pi_matrix, A_matrix, B_matrix, pattern_list)
# hmm_probabilities.append(hmm_probability)
# print('Start','p(O|lambda=)',hmm_probability)

# for i in range(1000):
#     pi_matrix, A_matrix, B_matrix = baum_welch(pi_matrix, A_matrix, B_matrix, pattern_list)
#     hmm_probability = forward_probability(pi_matrix, A_matrix, B_matrix, pattern_list)
#     hmm_probabilities.append(hmm_probability)
#     x_values.append(x_values[-1]+1)
#     print(i,'p(O|lambda=)',hmm_probabilit
# print(pi_matrix,'\n', A_matrix,'\n', B_matrix)

forward_result = forward(pi_matrix, A_matrix, B_matrix, pattern_list)
print("Forward result",'\n',forward_result)
backward_result = backward(pi_matrix, A_matrix, B_matrix, pattern_list)
print("Backward result",'\n',backward_result)
viterbi_result =  viterbi(pi_matrix, A_matrix, B_matrix, pattern_list)
print("Viterbi result",'\n',viterbi_result)
gamma_result = gamma(forward_result, backward_result, number_of_states, number_of_observations)
print("Gamma result",'\n',gamma_result)
xi_result = xi(forward_result, backward_result, A_matrix, B_matrix, number_of_states, number_of_observations)
print("Xi result",'\n',xi_result)
# A_matrix_new = new_A_matrix(xi_result, gamma_result, number_of_states, number_of_observations)
# B_matrix_new = new_B_matrix(gamma_result, number_of_states, pattern_list)
# new_pi = new_pi_matrix(gamma_result)


# ********** Plot the result
# fig,ax = plt.subplots()
# linewidth = 1.0
# #print(data_array[:,0])
# #for i in range(1,len(hmm_probabilities)):
# ax.plot(x_values, hmm_probabilities, lw = linewidth, color = "black")
# ax.set_xlabel('Iteration', labelpad=5, color='k', fontname='Arial', fontsize=16)
# ax.set_ylabel(r'p(O|$\lambda$)', labelpad=5, color='k', fontname='Arial', fontsize=16)
# ax.set_title('',fontname='Arial', fontsize=18, pad=10, loc='center')
# #ax.set_aspect(45/15)
# ax.tick_params(direction='in', pad=10,bottom=True, top=True, left=True, right=True)
# #ax.set_ticks_position('both')
# #plt.yticks(size='12')
# #ax.xaxis.set_major_locator(MultipleLocator(10))
# #plt.xticks([0,5,10,15,20], ['0','5','10','15','20'],size='12')
# plt.tight_layout()
# plt.show()