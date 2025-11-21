import numpy as np
import matplotlib.pyplot as plt
import sympy as sp
import control as ctrl

# -------------------------------------------------------
#  Funções auxiliares
# -------------------------------------------------------


def eps_from_ap(a_pass):
    """
    Calcula epsilon a partir do ripple de passagem em dB.
    """
    return np.sqrt(10**(0.1 * a_pass) - 1)


def cheb_order(a_pass, a_stop, w_pass, w_stop):
    """
    Calcula a ordem mínima do filtro Chebyshev Tipo I.
    - a_pass e a_stop em dB
    - w_pass e w_stop em rad/s (ou equivalentes, pois é normalizado)
    """
    eps = eps_from_ap(a_pass)

    num = np.acosh(np.sqrt((10**(0.1*a_stop) - 1)/(10**(0.1*a_pass) - 1)))
    den = np.acosh(w_stop / w_pass)

    n = abs(np.ceil(num / den))
    n = int(n)
    return n, eps


def chebyshev_poly(n, x):
    """
    Calcula o polinômio de Chebyshev de 1a espécie:
    - para |x| <= 1:  Cn(x) = cos(n * arccos(x))
    - para |x|  > 1:  Cn(x) = cosh(n * arccosh(x))
    """
    x = np.atleast_1d(x)

    C = np.zeros_like(x, dtype=float)

    # Região |x| <= 1
    mask1 = np.abs(x) <= 1
    C[mask1] = np.cos(n * np.arccos(x[mask1]))

    # Região |x| > 1
    mask2 = np.abs(x) > 1
    C[mask2] = np.cosh(n * np.acosh(x[mask2]))

    return C if len(C) > 1 else C[0]

# -------------------------------------------------------
#  Cálculo dos polos do Chebyshev Tipo I CORRIGIDO
# -------------------------------------------------------


def cheb_poles(n, eps, w_c=1.0):  # CORREÇÃO: adicionado w_c com valor default
    """
    Calcula os polos do Chebyshev Tipo I com escalonamento de frequência.
    """
    beta = np.arcsinh(1/eps) / n  # CORREÇÃO: usar arcsinh em vez de acosh
    poles = []

    for k in range(1, n+1):
        theta = (2*k - 1) * np.pi / (2*n)
        sigma = -np.sinh(beta) * np.sin(theta)
        omega = np.cosh(beta) * np.cos(theta)
        # Escalonar para a frequência de corte desejada
        pole = w_c * (sigma + 1j*omega)
        poles.append(pole)

    return np.array(poles)

# -------------------------------------------------------
#  Função de transferência do Chebyshev CORRIGIDA
# -------------------------------------------------------


def cheb_transfer(n, eps, w_c):
    poles = cheb_poles(n, eps, w_c)  # Agora funciona com 3 argumentos

    # Constante do ganho DC
    if n % 2 == 0:  # n par
        G0 = 1 / np.sqrt(1 + eps**2)
    else:  # n ímpar
        G0 = 1

    # Construir polinômio do denominador
    den_poly = np.poly(poles)
    den = np.real(den_poly)

    # Ajustar o numerador para ganho DC correto
    num = [G0 * den[0]]

    return num, den, poles

# -------------------------------------------------------
#  PRINT FUNÇÃO DE TRANSFERÊNCIA
# -------------------------------------------------------
