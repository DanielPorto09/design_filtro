# FUNÇÕES E BIBLIOTECAS PUBLICAS E AUTORAIS UTILIZADAS
import sys
from Chebychev import cheb_transfer, cheb_order, eps_from_ap, cheb_poles, cheb_transfer_pa, cheb_order_pa
from comp_eletronicos import n_filtros_PA, n_filtros_PB, separa_func
from graphs import anlize_em_frequencia, PrintaPlota_lista_de_filtros
import numpy as np
import control as ctrl

# Variaveis globais
choice = int(input("Escolha 1 (PB) ou 2 (PA): "))
filtros = []

# Caminho para encontrar pasta
sys.path.append(
    r"C:\Users\danie\OneDrive\Área de Trabalho\Facul\PDS\P2\codigos\design_filtro")


def print_transfer(num, den, poles):
    print("\n===============================")
    print(" Função de Transferência H(s) ")
    print("===============================\n")
    print("Numerador:")
    print(num, "\n")
    print("Denominador ( em Polos):")
    for i, p in enumerate(poles, start=1):
        print(f"Polo {i}:  {p.real:+.5f}  {p.imag:+.5f}j")


if choice == 1:

    print("== Filtro Chebyshev Passa-Baixa == ")
    # Dados do Filtro
    f_pass = float(input("Add frequencia passante (Hz):"))
    f_stop = float(input("Add frequencia de Corte (Hz):"))
    a_pass = float(input("Atenuação de banda passante(dB):"))
    a_stop = float(input("Atenuação de banda de corte(dB):"))

    # Normalização
    w_pass = 2*np.pi*f_pass
    w_stop = 2*np.pi*f_stop

    # Cálculo da ordem e de ε
    n, eps = cheb_order(a_pass, a_stop, w_pass, w_stop)

    print("Ordem mínima:", n)
    print("Epsilon:", eps)

    num, den, poles = cheb_transfer(n, eps, w_pass)
    print_transfer(num, den, poles)

    H = ctrl.TransferFunction(num, den)
    # print("\nH(T) NA FORMA PADRÃO", H)

    k_filters_fst_order, k_filters_sec_order = n_filtros_PB(a_pass, a_stop, w_pass, w_stop)
    filtros = separa_func(num, den, k_filters_fst_order, k_filters_sec_order)

elif choice == 2:

    print("== Filtro Chebyshev Passa-Alta == ")
    # Dados do Filtro
    f_pass = float(input("Add frequencia passante (Hz):"))
    f_stop = float(input("Add frequencia de Corte (Hz):"))
    a_pass = float(input("Atenuação de banda passante(dB):"))
    a_stop = float(input("Atenuação de banda de corte(dB):"))

    # Normalização
    w_pass = 2*np.pi*f_pass
    w_stop = 2*np.pi*f_stop

    # Cálculo da ordem e de ε
    n, eps = cheb_order_pa(a_pass, a_stop, w_pass, w_stop)

    print("Ordem mínima:", n)
    print("Epsilon:", eps)

    num, den, poles = cheb_transfer_pa(n, eps, w_pass)
    print_transfer(num, den, poles)

    H = ctrl.TransferFunction(num, den)
    # print("\nH(T) NA FORMA PADRÃO", H)

    k_filters_fst_order, k_filters_sec_order = n_filtros_PA(a_pass, a_stop, w_pass, w_stop)
    filtros = separa_func(num, den, k_filters_fst_order, k_filters_sec_order)

# PLOTA AS COISAS -> FUNÇÕES DE GRAPHS
# anlize_em_frequencia(H, f_pass, f_stop, a_pass, a_stop, w_pass, w_stop)
# PrintaPlota_lista_de_filtros(filtros)
