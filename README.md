# generate_HMM_pattern
A python program to generate an HMM pattern from a supplied pi-, A- and B-matrix and symbol list
Specify a pi, A and B matrix.  For example:

pi_matrix = [0.5,0.5]
A_matrix = [[0.95,0.05],[0.1,0.9]]
B_matrix = [[1/6,1/6,1/6,1/6,1/6,1/6],[0.1,0.1,0.1,0.1,0.1,0.5]]
list_of_symbols = [1,2,3,4,5,6]
list_of_states = ['F','L']

This is the example from Durbin et al. Biological Sequence Analysis. Cambridge University Press. 1998, where a caniso player used a fair (F, state=0) or loaded (L, state = 1) dice.  
He/she picks a dice randomly (pi_matrix = [0.5,0.5]) and chooses to stay with the fair dice 95% or the time, or to stay with the loaded dice 90% of the time, i.e.,
A_matrix = [[0.95,0.05],[0.1,0.9]].
The fair dice has an exqual probability of showing any of the 6 fases, whereas the loaded dice has a 50% chance to display "6", and a 0.1 chance of displaying any of the other 
five faces. B_matrix = [[1/6,1/6,1/6,1/6,1/6,1/6],[0.1,0.1,0.1,0.1,0.1,0.5]].

generate_HMM_pattern uses these matrixes to generate a series of symbols based on the HMM parameters.  You can obviously specify your own symbol list, number of states,
and pi-, A- and B-matrixes.
