# FUNÇÕES E BIBLIOTECAS PUBLICAS E AUTORAIS UTILIZADAS
import sys
from Chebychev import cheb_transfer, cheb_order, cheb_transfer_pa, cheb_order_pa
from tratamento_TF import n_filtros_PA, n_filtros_PB, separa_func, add_zeros_em_FT
from comp_eletronicos_pb import process_filter_list
from comp_eletronicos_pa import process_filter_list_pa
from graphs import anlize_em_frequencia, PrintaPlota_lista_de_filtros, print_resultados
import numpy as np
import control as ctrl
import sympy as sp

# Variaveis globais
choice = int(input("Escolha 1 (PB) ou 2 (PA): "))
filtros = []

# Caminho para encontrar pasta
sys.path.append(
    r"C:\Users\danie\OneDrive\Área de Trabalho\Facul\PDS\P2\codigos\design_filtro")


# teste de saida dos filtros
def extract_from_any(f):
    # Extrai num e den de vários formatos possíveis.

    # Caso seja TransferFunction
    if isinstance(f, ctrl.TransferFunction):
        num = f.num[0][0]
        den = f.den[0][0]
        return list(map(float, num)), list(map(float, den))

    # Caso seja dicionário
    if isinstance(f, dict):
        num = f["num"]
        den = f["den"]
        return list(map(float, num)), list(map(float, den))

    # Caso contrário: não reconhecido
    raise TypeError("Formato de filtro não suportado.")


def print_filtros(filtros):

    for idx, filtro in enumerate(filtros):
        if not isinstance(filtro, (dict, ctrl.TransferFunction)):
            print(f"ERRO: tipo inesperado dentro de filtros -> {type(filtro)}")
            print(f"Conteúdo:", filtro)
            continue

        print("\n==============================")
        print(f"  Filtro {idx+1}")
        print("==============================")

        num, den = extract_from_any(filtro)

        print("\nNumerador:")
        for i, c in enumerate(num):
            print(f"  a{i}: {c}")

        print("\nDenominador:")
        for i, c in enumerate(den):
            print(f"  b{i}: {c}")

        print("==============================\n")


def print_transfer_poles(num, den, poles):
    print("\n===============================")
    print(" Função de Transferência H(s) ")
    print("===============================\n")
    print("Numerador:")
    print(num, "\n")
    print("Denominador ( em Polos):")
    for i, p in enumerate(poles, start=1):
        print(f"Polo {i}:  {p.real:+.5f}  {p.imag:+.5f}j")


if choice == 1:
    # CALCULO DA FT FILTRO(S) ANALOGICO
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
    print_transfer_poles(num, den, poles)

    H = ctrl.TransferFunction(num, den)
    print("\nH(T) NA FORMA PADRÃO", H)

    k_filters_fst_order, k_filters_sec_order = n_filtros_PB(
        a_pass, a_stop, w_pass, w_stop)
    filtros = separa_func(num, den, k_filters_fst_order, k_filters_sec_order)
# TESTE
    print_filtros(filtros)
# CALCULO DOS VALORES DOS COMPONENTES ELETRÔNICOS
    resultados = process_filter_list(filtros)

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
    print_transfer_poles(num, den, poles)

    H = ctrl.TransferFunction(num, den)
    print("\nH(T) NA FORMA PADRÃO", H)

    k_filters_fst_order, k_filters_sec_order = n_filtros_PA(
        a_pass, a_stop, w_pass, w_stop)
    filtros = separa_func(num, den, k_filters_fst_order, k_filters_sec_order)
    filtros = add_zeros_em_FT(filtros)

    # CALCULO DOS VALORES DOS COMPONENTES ELETRÔNICOS
    resultados = process_filter_list_pa(filtros)


# PLOTA AS COISAS -> FUNÇÕES DE GRAPHS
anlize_em_frequencia(H, f_pass, f_stop, a_pass, a_stop, w_pass, w_stop)
PrintaPlota_lista_de_filtros(filtros)
# print_resultados(resultados)
