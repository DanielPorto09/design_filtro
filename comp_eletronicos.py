from Chebychev import cheb_transfer, cheb_order
import sys
sys.path.append(
    r"C:\Users\danie\OneDrive\Área de Trabalho\Facul\PDS\P2\codigos\design_filtro")


# SEPARANDO EM FILTROS DE MENOR ORDEM

def n_filtros(a_pass, a_stop, w_pass, w_stop):

    n, trash = cheb_order(a_pass, a_stop, w_pass, w_stop)
    k_fillters_sec_order = 0
    k_fillters_frt_order = 0

    if n % 2 == 1:
        k_filters_fst_order = 1
        k_filters_sec_order = (n - 1) // 2
    else:
        k_filters_sec_order = n // 2

    return k_filters_fst_order, k_filters_sec_order
