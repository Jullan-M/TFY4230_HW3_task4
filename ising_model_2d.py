__author__ = 'Jullan'
# -*- coding: utf-8 -*-
# Made by Jullan
import numpy as np
from sympy.utilities.iterables import permute_signs
import timeit
from matplotlib import pyplot as plt


def get_cmap(n, name='hsv'):
    '''Returns a function that maps each index in 0, 1, ..., n-1 to a distinct
    RGB color; the keyword argument name must be a standard mpl colormap name.'''
    return plt.cm.get_cmap(name, n)


KB = 1.38064852e-23
J_PAIR = np.array([1, 2])
B = 5
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
    #   Old, unused function
    arr_len = len(sigma_arr)
    P_arr = np.zeros((arr_len, arr_len))
    for l in range(arr_len):
        for m in range(arr_len):
            #   These _sum-variables are sums over the chains, k.
            site_sum = np.sum(sigma_arr[l]*sigma_arr[m])
            chain_sum = np.sum(sigma_arr[l]*np.roll(sigma_arr[l], -1))
            b_sum = np.sum(sigma_arr[l] + sigma_arr[m])
            P_arr[l, m] = np.exp(
                beta * (J[0] * site_sum + J[1] * chain_sum + B/2 * b_sum))
    return P_arr


class Ising_2D:
    def __init__(self, n, beta=BETA, kB=KB, J=J_PAIR, B=B):
        #   Various constants/variables
        self.beta = beta
        self.kB = kB
        self.J = J
        self.B = B
        if (n == 2):
            self.B = B/2
        elif (n == 1):
            self.B = 0

        #   Length of arrays, chains, sites.
        self.n = n
        self.P_len = 2**n

        #   Arrays/matrices
        self.sigma_arr = generate_spin_conf(self.n)
        self.P_matrix_nopow = np.empty((self.P_len, self.P_len))
        self.P_matrix = np.empty_like(self.P_matrix_nopow)
        self.generate_P_matrix()

    def update_sigma_arr(self):
        self.sigma_arr = generate_spin_conf(self.n)

    def generate_P_matrix(self):
        for l in range(self.P_len):
            for m in range(self.P_len):
                #   These _sum-variables are sums over the chains, k.
                site_sum = np.dot(self.sigma_arr[l], self.sigma_arr[m])
                chain_sum = np.dot(
                    self.sigma_arr[l], np.roll(self.sigma_arr[l], -1))
                b_sum = np.sum(self.sigma_arr[l] + self.sigma_arr[m])

                self.P_matrix_nopow[l, m] = np.exp(
                    self.J[0] * site_sum + self.J[1] * chain_sum + self.B / 2 * b_sum)
        self.P_matrix = self.P_matrix_nopow ** self.beta

    def refresh_P_matrix_beta(self):
        self.P_matrix = self.P_matrix_nopow**self.beta

    def find_lambdas(self):
        return np.linalg.eigvals(self.P_matrix)


def lambda_beta(ising, betas):
    lambdas_beta = np.empty((len(betas), 2**ising.n))
    for j, beta in enumerate(betas):
        #print("hey!,", i)
        ising.beta = beta
        ising.refresh_P_matrix_beta()
        lambdas_beta[j] = np.real(ising.find_lambdas())
    return np.sort(lambdas_beta)


def lambda_betaB(ising, Bs, betas):
    lambdas_betaB = np.empty((len(Bs), len(betas)))
    for i, B in enumerate(Bs):
        ising.B = B
        ising.generate_P_matrix()
        lambdas_betaB[i] = lambda_beta(ising, betas)[:, -1]
    return lambdas_betaB


def magn_betaB(ising, lambdas_betaB):
    return np.gradient(np.log(lambdas_betaB), Bs, axis=0) / betas / ising.n


if (__name__ == "__main__"):

    #   Oppgave 4, a)
    sigma_arr = generate_spin_conf(5)
    print(sigma_arr)
    print(len(sigma_arr))

    #   Oppgave 4, b)
    is1 = Ising_2D(2)
    print(is1.P_matrix)

    #   Oppgave 4, c)

    #   Setup
    n_max = 10
    betas = np.linspace(0.005, 1.000, 101)

    for i in range(1, n_max + 1):
        ising = Ising_2D(i)
        print("n=", ising.n)
        cmap = get_cmap(ising.P_len)
        lambdas_arr = lambda_beta(ising, betas)
        plt.figure()
        plt.title(r"$n$ = " + str(ising.n) + "; $J_{\parallel}$,$J_{\\bot}$ = " + str(
            ising.J[0]) + "," + str(ising.J[1]) + "; $B$ = " + str(ising.B))
        for i in range(2**ising.n):
            plt.semilogy(betas, lambdas_arr[:, i], color=cmap(i))
        plt.xlabel(r"$\beta$", fontsize=16)
        plt.ylabel(r"$\lambda$-values", fontsize=16)
        plt.grid()
        plt.savefig(f'eigvals_{ising.n:.1g}.pdf')
        plt.show()

    #   Oppgave 4, d)
    betas = np.linspace(0.005, 1.000, 101)
    Bs = np.linspace(-1, 1, 201)
    is1 = Ising_2D(5)

    lambdas = lambda_betaB(is1, Bs, betas)
    magns = magn_betaB(is1, lambdas)
    bB, Bb = np.meshgrid(betas, Bs)

    fig, ax = plt.subplots()
    cs = ax.contourf(bB, Bb, magns, levels=np.linspace(-1, 1, 101))
    cbar = fig.colorbar(cs, ticks=(-1, 0, 1))
    cbar.set_label(r'$m/n$')
    plt.xlabel(r'$\beta$')
    plt.ylabel(r'$H$')
    plt.savefig('magnetization.pdf')

    #   Oppgave 4, e)
    enthalpy = np.gradient(- np.log(lambdas), betas, axis=1) / is1.n
    spec_heat = - betas**2 * np.gradient(enthalpy, betas, axis=1)

    fig, ax = plt.subplots()
    cs = ax.contourf(bB, Bb, enthalpy, 100)
    cbar = fig.colorbar(cs)
    cbar.set_label(r'$H_e(\beta)/N$')
    plt.xlabel(r'$\beta$')
    plt.ylabel(r'$B$')
    plt.savefig('enthalpy.pdf')

    fig, ax = plt.subplots()
    cs = ax.contourf(bB, Bb, spec_heat, 100)
    cbar = fig.colorbar(cs)
    cbar.set_label(r'$C_B(\beta)/Nnk_B$')
    plt.xlabel(r'$\beta$')
    plt.ylabel(r'$B$')
    plt.savefig('specific_heat.pdf')

    plt.show()
