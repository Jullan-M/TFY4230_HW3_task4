__author__ = 'Jullan'
# -*- coding: utf-8 -*-
#Made by Jullan
import numpy as np
from sympy.utilities.iterables import permute_signs
import timeit
from matplotlib import pyplot as plt

KB = 1.38064852e-23
J_PAIR = np.array([1,1])
B = 1
BETA = 1

#   Generates the possible 2^n spin-configurations
#   for the n spins omega_m,i on the n chains at a
#   given i.
def generate_spin_conf(n):
    sample_tup = tuple(np.ones(n))
    return np.asarray(list(permute_signs(sample_tup)))

def generate_spin_conf_slow(n):
    #   Slow, but intuitive implementation
    sigma = np.ones((2**n, n))
    for l in range(2**n):
        for k in range(n):
            sigma[l][k] = (-1)**(l//(2**k))
    return sigma

def generate_P_matrix(sigma_arr, J, B, k, n, beta):
    arr_len = len(sigma_arr)
    P_arr = np.zeros((arr_len, arr_len))
    for l in range(arr_len):
        for m in range(arr_len):
            #   These _sum-variables are sums over the chains, k.
            site_sum = np.sum(sigma_arr[l]*sigma_arr[m])
            chain_sum = np.sum(sigma_arr[l]*np.roll(sigma_arr[l], -1))
            b_sum = np.sum(sigma_arr[l] + sigma_arr[m])
            P_arr[l, m] = np.exp( beta * (J[0] * site_sum + J[1] * chain_sum + B/2 * b_sum))
    return P_arr

class Ising_2D:
    def __init__(self, n, beta = BETA, kB = KB, J = J_PAIR, B = B):
        #   Various constants/variables
        self.beta = beta
        self.kB = kB
        self.J = J
        self.B = B

        #   Length of arrays, chains, sites.
        self.n = n
        self.P_len = 2**n

        #   Arrays/matrices
        self.sigma_arr = generate_spin_conf(n)
        self.P_matrix = np.zeros((self.P_len, self.P_len))
        self.lambdas = np.zeros(self.P_len)

    def generate_P_matrix(self):
        for l in range(self.P_len):
            for m in range(self.P_len):
                #   These _sum-variables are sums over the chains, k.
                site_sum = np.sum(self.sigma_arr[l] * self.sigma_arr[m])
                chain_sum = np.sum(self.sigma_arr[l] * np.roll(self.sigma_arr[l], -1))
                b_sum = np.sum(self.sigma_arr[l] + self.sigma_arr[m])

                self.P_matrix[l, m] = np.exp(self.beta * (self.J[0] * site_sum + self.J[1] * chain_sum + self.B / 2 * b_sum))










#   Oppgave 4, a)
sigma_arr = generate_spin_conf(5)
print(sigma_arr)
print(len(sigma_arr))

#   Oppgave 4, b)
is1 = Ising_2D(1)
is1.generate_P_matrix()
print(is1.P_matrix)




