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
    Retorna SEMPRE uma lista plana de TransferFunctions.
    """
    # dividindo o ganho entre as fts menores
    k_filters = k_filters_fst_order + k_filters_sec_order
    gain_total = num[0]   # numerador original
    pedaco_ganho = gain_total ** (1.0 / k_filters)

    H_total = ctrl.TransferFunction(num, den)

    # Polos da função completa
    poles = np.roots(den)

    filtros = []
    usados = np.zeros(len(poles), dtype=bool)

    # ------------------------------------------------------------
    # 1ª ORDEM → polos reais
    # ------------------------------------------------------------
    for i, p in enumerate(poles):
        if usados[i]:
            continue

        if np.isclose(p.imag, 0, atol=1e-9):  # polo real
            den1 = [1, -p.real]
            # separando entre partes menores
            num1 = [pedaco_ganho]

            F = ctrl.TransferFunction(num1, den1)

            filtros.append(F)
            usados[i] = True

            # conta quantos já foram montados
            count_1st = sum(len(f.den[0][0]) == 2 for f in filtros)
            if count_1st == k_filters_fst_order:
                break

    # ------------------------------------------------------------
    # 2ª ORDEM → polos complexos conjugados
    # ------------------------------------------------------------
    for i, p in enumerate(poles):
        if usados[i]:
            continue

        # procura o conjugado
        conj_idx = None
        for j in range(len(poles)):
            if i != j and not usados[j]:
                if np.isclose(poles[j].real,  p.real, atol=1e-9) and \
                   np.isclose(poles[j].imag, -p.imag, atol=1e-9):
                    conj_idx = j
                    break

        if conj_idx is None:
            continue  # não achou conjugado → pula

        # se chegou aqui → temos par conjugado
        alpha = p.real
        beta = p.imag
        den2 = [1, -2*alpha, alpha*alpha + beta*beta]
        num2 = [pedaco_ganho]

        F = ctrl.TransferFunction(num2, den2)

        filtros.append(F)
        usados[i] = True
        usados[conj_idx] = True

        # conta quantos já foram montados
        count_2nd = sum(len(f.den[0][0]) == 3 for f in filtros)
        if count_2nd == k_filters_sec_order:
            break

    return filtros
