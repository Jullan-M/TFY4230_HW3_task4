__author__ = 'Jullan'
# -*- coding: utf-8 -*-
#Made by Jullan
import numpy as np
from sympy.utilities.iterables import permute_signs
import timeit
from matplotlib import pyplot as plt


#   Generates the possible 2^n spin-confirgurations
#   for the n spins omega_m,i on the n chains at a
#   given i.
kB = 1.38064852e-23
J = 1
B = 10


def generate_spin_conf(n):
    sample_tup = tuple(np.ones(n))
    return np.asarray(list(permute_signs(sample_tup)))

def generate_spin_conf_slow(n):
    #   Slow implementation    #
    sigma = np.ones((2**n, n))
    for l in range(2**n):
        for k in range(n):
            sigma[l][k] = (-1)**(l//(2**k))
    return sigma

def generate_P_matrix(sigma_arr, J, B):
    arr_len = len(sigma_arr)
    P_arr = np.zeros((arr_len, arr_len))
    for l in range(arr_len):
        P_arr[l] = np.exp(- sigma_arr[l]*sigma_arr[l] - sigma_arr[l] * sigma_arr[])


#   Oppgave 4, a)
sigma_arr = generate_spin_conf(5)
print(sigma_arr)
print(len(sigma_arr))





