from math import sqrt
from Chebychev import cheb_transfer, cheb_order, cheb_order_pa
import numpy as np
import control as ctrl
import sympy as sp


import sys
sys.path.append(
    r"C:\Users\danie\OneDrive\Área de Trabalho\Facul\PDS\P2\codigos\design_filtro")

# ===================================
# CALCULO DOS COMPONENTES DOS FILTROS
# ===================================


def solve_first_order(num, den):
    """
    Resolve o sistema para filtros de 1ª ordem.
    H(s) = a / (b*s + c)

    Modelo:
        a = R2/R1
        b = R2*C1
        c = 1  (normalizado)
    """
    print("\n--- Filtro de Primeira Ordem ---")

    R1 = float(input("Escolha um valor para R1 [Ohms]: "))

    a = float(num[0])     # a = R2/R1
    b = float(den[0])     # b = R2*C1
    c = float(den[1])     # normalmente = 1

    R2 = a * R1
    C1 = b / R2

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

# utilitário: verificar se resistores estão em faixa prática


def _resistors_reasonable(Rvals, R_min=1.0, R_max=1e7):
    return all((R_min <= rv <= R_max) for rv in Rvals)

# === MFB: equações corretas (retornam dict) ===


def _solve_mfb_given_caps(a_val, c_val, d_val, C1, C2, chute=[10e3, 10e3, 10e3]):
    """
    Resolve MFB para R1, R2, R3 dado C1,C2 e coeficientes a,c,d.
    Retorna dict com R1,R2,R3 ou lança Exception se não convergir.
    """
    R1, R2, R3 = sp.symbols("R1 R2 R3", positive=True)

    eq1 = d_val - 1/(R2*R3*C1*C2)  # wn^2
    eq2 = c_val - (1/C1)*(1/R2 + 1/R3 + 1/R1)  # wn/Q
    eq3 = a_val - (R3/R1)  # ganho DC

    try:
        sol = sp.nsolve([eq1, eq2, eq3], [R1, R2, R3],
                        chute, tol=1e-6, maxsteps=60)
    except Exception as e:
        raise

    R1v = float(sol[0])
    R2v = float(sol[1])
    R3v = float(sol[2])

    return {"R1": R1v, "R2": R2v, "R3": R3v}

# === SALLEN-KEY: equações (non-inverting SK with gain K) ===


def _solve_sallen_key(a_val, c_val, d_val, C, chute=[10e3, 10e3, 1.0]):
    """
    Resolve Sallen-Key não-inversor (com ganho K a determinar),
    assumindo C1 = C2 = C.
    Variáveis: R1, R2, K
    Equações (denominator comparadas com s^2 + c s + d):
       1) d = 1/(R1*R2*C*C)
       2) c = (C*(R1 + R2) + C*R1*(1 - K)) / (R1*R2*C*C)
          => simplifica para numerator /(R1*R2*C*C)
       3) a = K   (DC gain)
    """
    R1, R2, K = sp.symbols("R1 R2 K", positive=True)

    eq1 = d_val - 1/(R1*R2*C*C)
    # termo s: C*(R1 + R2) + C*R1*(1 - K)  -> divide por (R1*R2*C*C)
    eq2 = c_val - (C*(R1 + R2) + C*R1*(1 - K)) / (R1*R2*C*C)
    eq3 = a_val - K

    try:
        sol = sp.nsolve([eq1, eq2, eq3], [R1, R2, K],
                        chute, tol=1e-6, maxsteps=80)
    except Exception as e:
        raise

    return {"R1": float(sol[0]), "R2": float(sol[1]), "K": float(sol[2])}

# === função principal que tenta MFB e faz fallback SK ===


def solve_biquad_svf(a_val, c_val, d_val):
    import math

    # frequência natural
    wn = math.sqrt(d_val)

    # Q
    Q = wn / c_val

    # escolha automática de C
    C = 10e-9   # 10 nF padrão

    # escolha de R integrador
    R_int = 1/(wn*C)

    # ganho → Rg/Rf
    K = a_val

    # retorno no formato
    return {
        "tipo": 2,
        "topologia": "SVF",
        "R1": R_int,
        "R2": R_int,
        "C1": C,
        "C2": C,
        "Ra": K,
        "Rb": 0,
        "K": K
    }

# adapta process_filter_list para usar solve_biquad_auto (substitui chamada anterior)


def process_filter_list(filtros):
    resultados = []

    for f in filtros:
        num = f.num[0][0]
        den = f.den[0][0]
        ordem = len(den) - 1

        if ordem == 1:
            # usa sua função existente
            data = solve_first_order(num, den)
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
