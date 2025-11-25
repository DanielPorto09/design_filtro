import math
from comp_eletronicos_pb import solve_biquad_svf

import sys
sys.path.append(
    r"C:\Users\danie\OneDrive\Área de Trabalho\Facul\PDS\P2\codigos\design_filtro")


def solve_first_order_pa(num, den):
    """
    Resolve o sistema para filtros de 1ª ordem.
    H(s) = a*s / (b*s + c)

    Modelo:
        a = R2/R1
        b = 1
        c = 1/C*R1
    """
    print("\n--- Filtro de Primeira Ordem ---")

    R1 = float(input("Escolha um valor para R1 [Ohms]: "))

    a = float(num[0])     # a = R2/R1
    b = float(den[0])     # b = 1
    c = float(den[1])     # c = 1/C*R1

    R2 = a * R1
    C1 = 1 / (c*R1)

    return {
        "tipo": 1,
        "R1": R1,
        "R2": R2,
        "C1": C1,
        "C2": 0,
        "Ra": 0,
        "Rb": 0,
        "K": a
    }

# adapta process_filter_list para usar solve_biquad_auto (substitui chamada anterior)


def process_filter_list_pa(filtros):
    resultados = []

    for f in filtros:
        num = f.num[0][0]
        den = f.den[0][0]
        ordem = len(den) - 1

        if ordem == 1:
            # usa sua função existente
            data = solve_first_order_pa(num, den)
        elif ordem == 2:
            # chama o projetista automático (MFB first, SK fallback)
            a_val = float(num[0])
            c_val = float(den[1])
            d_val = float(den[2])
            data = solve_biquad_svf(a_val, c_val, d_val)

        else:
            raise ValueError("Filtro não é 1ª nem 2ª ordem!")

        resultados.append(data)

    return resultados
