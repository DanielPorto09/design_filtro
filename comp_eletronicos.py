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
    k_filters_sec_order = 0
    k_filters_frt_order = 0

    if n % 2 == 1:
        k_filters_fst_order = 1
        k_filters_sec_order = (n - 1) // 2
    else:
        k_filters_fst_order = 0
        k_filters_sec_order = n // 2

    return k_filters_fst_order, k_filters_sec_order


def separar_ordens(poles):
    polos_real = []
    polos_complex = []
    usados = np.zeros(len(poles), dtype=bool)

    for i, p in enumerate(poles):
        if usados[i]:
            continue

        if np.isclose(p.imag, 0):
            polos_real.append(p.real)
            usados[i] = True
        else:
            conj_index = np.where(
                np.isclose(poles.real, p.real) &
                np.isclose(poles.imag, -p.imag) &
                (~usados)
            )[0]
            if len(conj_index) > 0:
                j = conj_index[0]
                polos_complex.append((p, poles[j]))
                usados[i] = usados[j] = True

    return polos_real, polos_complex


def gerar_secoes_normalizadas(polos_real, polos_complex, w_c):
    secoes = []

    # 1ª ORDEM
    for p in polos_real:
        a0 = 1 / (-p / w_c)
        b0 = 1  # ganho normalizado
        den = [1, a0]
        num = [b0]
        H = ctrl.TransferFunction(num, den)
        secoes.append(H)

    # 2ª ORDEM
    for p, pc in polos_complex:
        alpha = p.real
        beta = p.imag

        w0 = np.sqrt(alpha**2 + beta**2)
        Q = w0 / (-2 * alpha)

        a1 = 1/(Q)
        a0 = 1/(1)
        den = [1, a1, a0]
        num = [1]
        H = ctrl.TransferFunction(num, den)
        secoes.append(H)

    return secoes

def separa_func(num, den, k_filters_fst_order, k_filters_sec_order, w_c, K=1.0):
    poles = np.roots(den)
    polos_real, polos_complex = separar_ordens(poles)
    secoes = gerar_secoes_normalizadas(polos_real, polos_complex, w_c)
    return secoes

# ==================================================================
# ENCONTRANDO COEFICIENTES DOS COMPONENTES DE CADA FILTRO DA CASCATA
# ==================================================================

R_LIMIT_DEFAULT = 1e6


def solve_fst_PB(b0, a0, K, C1, r_limit=R_LIMIT_DEFAULT):
    """
    Resolve seção de 1ª ordem .
    A TF esperada é H(s) = b0 / (a0*s + 1) (ou ajustada conforme normalização).

    R1 = 1 / (a0 * C1)
    """
    resultado = {"ordem": 1, "C1": C1, "K": K}

    if abs(a0) < 1e-30:
        resultado.update(
            {"error": "a0 quase zero, não é possível resolver R1."})
        return resultado

    R1 = 1.0 / (a0 * C1)

    if not (0 < R1 < r_limit):
        resultado.update(
            {"error": f"R1 encontrado fora do intervalo físico (R1={R1:.6g} Ω). Limite = {r_limit}"})
        return resultado

    resultado["R1"] = float(R1)

    if abs(b0 - K*a0) > max(1e-9, 1e-6 * abs(b0)):
        resultado["warning"] = f"Inconsistência: b0 != K*a0 (b0={b0:.6g}, K*a0={(K*a0):.6g})."
    else:
        resultado["ok"] = True

    return resultado


def solve_scd_PB(b0, a1, a0, K, C1, C2, r_limit=R_LIMIT_DEFAULT):
    """
     solver para seção de 2ª ordem baseado na topologia Sallen-Key passa-baixa.

    Usa SymPy para resolver as equações:
      (1) 1/(R1*R2*C1*C2) = a0
      (2) 1/(R1*C1) + 1/(R2*C1) + (1-K)/(R2*C2) = a1
      (3) K/(R1*R2*C1*C2) = b0

    Retorna dicionário com R1,R2 ou erro.
    """
    resultado = {"ordem": 2, "C1": C1, "C2": C2, "K": K}

    if abs(a0) < 1e-30:
        resultado.update(
            {"error": "a0 quase zero; impossível resolver (divisão por zero)."})
        return resultado

    if abs(b0 - K*a0) > max(1e-9, 1e-6 * abs(b0)):
        resultado["warning"] = f"b0 != K*a0 (b0={b0:.6g}, K*a0={(K*a0):.6g}) — resultado pode ser inválido."

    R1, R2 = sp.symbols('R1 R2', positive=True, real=True)

    eq1 = sp.Eq(1/(R1 * R2 * C1 * C2), a0)
    eq2 = sp.Eq((1/(R1 * C1)) + (1/(R2 * C1)) + ((1 - K)/(R2 * C2)), a1)
    eq3 = sp.Eq(K/(R1 * R2 * C1 * C2), b0)

    sols = []
    try:
        raw = sp.solve([eq1, eq2, eq3], [R1, R2], dict=True)
        if raw:
            for s in raw:
                try:
                    r1_val = float(sp.N(s[R1]))
                    r2_val = float(sp.N(s[R2]))
                except Exception:
                    continue
                sols.append((r1_val, r2_val))
    except Exception:
        raw = None

    reais_positivas = []
    for (r1_val, r2_val) in sols:
        if np.isreal(r1_val) and np.isreal(r2_val) and r1_val > 0 and r2_val > 0 and r1_val < r_limit and r2_val < r_limit:
            reais_positivas.append({"R1": float(r1_val), "R2": float(r2_val)})

    if not reais_positivas:
        resultado.update(
            {"error": "Nenhuma solução real positiva encontrada.", "raw_solutions": sols})
        return resultado

    # escolha heurística: menor resistor máximo
    escolha = min(reais_positivas, key=lambda s: max(s["R1"], s["R2"]))
    resultado.update(escolha)
    resultado["ok"] = True
    resultado["raw_solutions"] = reais_positivas

    return resultado


def processar_filtros_PB(lista_Hn):

    print("=== Parâmetros globais (serão usados para todas as seções) ===")
    Ra = float(input("Digite Ra (ohms): "))
    Rb = float(input("Digite Rb (ohms, diferente de 0): "))

    assert Rb != 0, "Rb não pode ser zero."

    K = 1.0 + (Ra / Rb)
    print(f"K calculado (original) = {K}")

    C1_global = float(
        input("Digite valor de C1 (Farad) — usado nas seções 1ª e 2ª ordem: "))
    C2_global = float(
        input("Digite valor de C2 (Farad) — usado apenas nas seções 2ª ordem: "))

    resultados = []

    for i, H in enumerate(lista_Hn, start=1):
        num = np.array(H.num[0][0], dtype=float)
        den = np.array(H.den[0][0], dtype=float)
        deg = len(den) - 1

        print(f"\n=== Processando filtro {i} (grau denom = {deg}) ===")

        # CORREÇÃO: Escalar numerador para compatibilidade b0 = K*a0
        if deg >= 1:
            a0 = den[-1]  # termo constante do denominador
            b0_atual = num[0] if len(num) > 0 else 0
            b0_desejado = K * a0
            
            if abs(b0_atual - b0_desejado) > 1e-6 and abs(b0_atual) > 1e-10:
                # Escalar todo o numerador
                fator = b0_desejado / b0_atual
                num = num * fator
                print(f"  Numerador escalado por fator {fator:.6f} para compatibilidade")

        if deg == 2:
            # den = [1, a1, a0]  (normalizado para monic)
            b0 = float(num[0]) if len(num) >= 1 else 0.0
            a1 = float(den[1])
            a0 = float(den[2])

            res = solve_scd_PB(b0, a1, a0, K, C1_global, C2_global)
            res.update({"filtro": i})
            resultados.append(res)

        elif deg == 1:
            # 1ª ordem: den = [1, a0]
            b0 = float(num[0]) if len(num) >= 1 else 0.0
            a0 = float(den[1])

            res = solve_fst_PB(b0, a0, K, C1_global)
            res.update({"filtro": i})
            resultados.append(res)

        else:
            resultados.append(
                {"filtro": i, "error": f"Denominador com grau {deg} não suportado por este solver."})

    return resultados