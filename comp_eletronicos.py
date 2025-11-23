from Chebychev import cheb_transfer, cheb_order, cheb_order_pa
import numpy as np
import control as ctrl
import sympy as sp


import sys
sys.path.append(
    r"C:\Users\danie\OneDrive\Área de Trabalho\Facul\PDS\P2\codigos\design_filtro")

# =====================================
# SEPARANDO EM FILTROS DE MENOR ORDEM
# ======================================

def n_filtros_PB(a_pass, a_stop, w_pass, w_stop):

    n, trash = cheb_order(a_pass, a_stop, w_pass, w_stop)
    k_fillters_sec_order = 0
    k_fillters_frt_order = 0

    if n % 2 == 1:
        k_filters_fst_order = 1
        k_filters_sec_order = (n - 1) // 2
    else:
        k_filters_fst_order = 0
        k_filters_sec_order = n // 2

    return k_filters_fst_order, k_filters_sec_order

def n_filtros_PA(a_pass, a_stop, w_pass, w_stop):

    n, trash = cheb_order_pa(a_pass, a_stop, w_pass, w_stop)
    k_fillters_sec_order = 0
    k_fillters_frt_order = 0

    if n % 2 == 1:
        k_filters_fst_order = 1
        k_filters_sec_order = (n - 1) // 2
    else:
        k_filters_fst_order = 0
        k_filters_sec_order = n // 2

    return k_filters_fst_order, k_filters_sec_order

def separa_func(num, den, k_filters_fst_order, k_filters_sec_order):
    """
    Decompõe H(s) em k1 filtros de 1ª ordem e k2 de 2ª ordem.
    """
    H_total = ctrl.TransferFunction(num, den)

    # Polos da função completa
    poles = np.roots(den)

    filtros = []
    usados = np.zeros(len(poles), dtype=bool)

    # ---- 1ª ORDEM ---------------------------------------------------
    for i in range(len(poles)):
        if usados[i]:
            continue

        p = poles[i]

        # Polo real → 1ª ordem
        if np.isclose(p.imag, 0):
            den1 = [1, -p.real]
            num1 = [1]
            F = ctrl.TransferFunction(num1, den1)
            filtros.append(F)
            usados[i] = True

            if len([f for f in filtros if f.den[0][0].size == 2]) == k_filters_fst_order:
                break

    # ---- 2ª ORDEM ---------------------------------------------------
    for i in range(len(poles)):
        if usados[i]:
            continue

        p = poles[i]

        # Procurar conjugado
        j = np.where(
            np.isclose(poles.real,  p.real) &
            np.isclose(poles.imag, -p.imag) &
            (~usados)
        )[0]

        if len(j) == 0:
            continue

        j = j[0]

        # coeficientes do polinômio s² - 2α s + (α² + β²)
        alpha = p.real
        beta = p.imag
        den2 = [1, -2*alpha, alpha*alpha + beta*beta]
        num2 = [1]

        F = ctrl.TransferFunction(num2, den2)
        filtros.append(F)

        usados[i] = True
        usados[j] = True

        if len([f for f in filtros if f.den[0][0].size == 3]) == k_filters_sec_order:
            break

    return filtros

