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
    # Teste de erro
    assert a_pass < a_stop, "\n===================================\nERRO, a_pass deve ser menor que a_stop\n==================================="
    assert w_pass < w_stop, "Para Chebyshev PB, w_pass < w_stop"

    eps = eps_from_ap(a_pass)

    num = np.acosh(np.sqrt((10**(0.1*a_stop) - 1)/(10**(0.1*a_pass) - 1)))
    den = np.acosh(w_stop / w_pass)

    n = abs(np.ceil(num / den))
    n = int(n)
    return n, eps


def cheb_order_pa(a_pass, a_stop, w_pass, w_stop):
    """
    Ordem de Chebyshev Tipo I para passa-alta.
    """
    # Transformação passa-alta → passa-baixa equivalente
    # LP equivalents
    wp_lp = 1
    ws_lp = w_pass / w_stop   # precisa ser > 1

    assert ws_lp > 1, "\n===================================\nERRO: para Chebyshev HP, w_pass deve ser MAIOR que w_stop \n==================================="

    eps = eps_from_ap(a_pass)

    num = np.acosh(np.sqrt((10**(0.1*a_stop)-1)/(10**(0.1*a_pass)-1)))
    den = np.acosh(ws_lp)

    n = int(np.ceil(num / den))
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
#  Cálculo dos polos do Chebyshev Tipo I P
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
#  Função de transferência do Chebyshev PB
# -------------------------------------------------------


def cheb_transfer(n, eps, w_c):
    poles = cheb_poles(n, eps, w_c)

    # Constante do ganho DC
    if n % 2 == 0:
        G0 = 1 / np.sqrt(1 + eps**2)
    else:
        G0 = 1

    # Construir polinômio do denominador (já está na frequência correta)
    den_poly = np.poly(poles)
    den = np.real(den_poly)

    # Ajustar o numerador para ganho DC correto
    num = [G0 * den[0]]

    return num, den, poles

# -------------------------------------------------------
#  Função de transferência do Chebyshev PA
# -------------------------------------------------------


def cheb_transfer_pa(n, eps, w_c):
    """
    Função de transferência Chebyshev Tipo I passa-alta.
    """
    # Polos do protótipo PB normalizado (wc = 1)
    poles_lp = cheb_poles(n, eps, 1.0)

    # Transformação LP -> HP: polos HP = w_c / polos_LP
    poles_hp = w_c / poles_lp

    # Denominador do filtro PA
    den = np.real(np.poly(poles_hp))

    # Ganho em alta frequência
    if n % 2 == 0:
        Ginf = 1.0 / np.sqrt(1.0 + eps**2)
    else:
        Ginf = 1.0

    # Numerador: Ginf * s^n
    num = [Ginf] + [0] * n

    return num, den, poles_hp