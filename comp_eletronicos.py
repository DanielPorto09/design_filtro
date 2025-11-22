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


# ==================================================================
# ENCONTRANDO COEFICIENTES DOS COMPONENTES DE CADA FILTRO DA CASCATA
# ==================================================================

def solve_fst_PB(b0, a0, K, C1):

    resultado = {"ordem": 1, "C1": C1, "K": K}
   # nao da pra a0 ser 0, se nao da M.
    if abs(a0) < 1e-30:
        resultado.update(
            {"error": "a0 quase zero, não é possível resolver R1."})
        return resultado

    # Cálculo do resistor
    R1 = 1.0 / (a0 * C1)
    R1 = float(R1)

    # --- Limite máximo para R1 ---
    limite = 1e6  # 1 megaohm
    if not (0 < R1 < limite):
        resultado.update({
            "error": (
                f"R1 encontrado fora do intervalo físico (R1={R1:.3g} Ω). "
                f"Limite máximo = {limite} Ω."
            )
        })
        return resultado

    resultado["R1"] = R1

    # Verificação de consistência b0 ≈ K*a0
    if abs(b0 - K*a0) > max(1e-9, 1e-6 * abs(b0)):
        resultado["warning"] = (
            f"Inconsistência: b0 != K*a0 (b0={b0:.6g}, K*a0={(K*a0):.6g})."
        )
    else:
        resultado["ok"] = True

    return resultado


def solve_scd_PB(b0, a1, a0, K, C1, C2):
    """
    Resolve R1, R2 para seção 2ª ordem.
    Retorna dict com R1, R2, C1, C2, K ou None se sem solução física.
    """

    resultado = {"ordem": 2, "C1": C1, "C2": C2, "K": K}

    # checagem básica: a0 não nulo
    if abs(a0) < 1e-30:
        resultado.update(
            {"error": "a0 quase zero; impossível resolver (divisão por zero)."})
        return resultado

    # checar consistência b0 == K*a0 (não exigimos, apenas avisamos)
    if abs(b0 - K*a0) > max(1e-9, 1e-6 * abs(b0)):
        resultado["warning"] = f"b0 != K*a0 (b0={b0:.6g}, K*a0={(K*a0):.6g}) — resultado pode ser inválido."

    # variáveis simbólicas
    R1, R2 = sp.symbols('R1 R2', positive=True, real=True)

    eq1 = sp.Eq(b0, K / (R1 * R2 * C1 * C2))
    eq2 = sp.Eq(a1, (1/(R1*C1)) + (1/(R2*C1)) + ((1-K)/(R2*C2)))

    sols = sp.solve([eq1, eq2], [R1, R2], dict=True)

    # filtrar soluções reais e positivas
    reais_positivas = []
    for s in sols:
        try:
            r1_val = float(s[R1])
            r2_val = float(s[R2])
        except Exception:
            continue

        limite = 1e6  # limite superior para resistores

        if (np.isreal(r1_val) and np.isreal(r2_val) and r1_val > 0 and r2_val > 0 and r1_val < limite and r2_val < limite):
            reais_positivas.append({"R1": r1_val, "R2": r2_val})

    if not reais_positivas:
        resultado.update(
            {"error": "Nenhuma solução real positiva encontrada.", "raw_solutions": sols})
        return resultado

    # escolher a primeira solução física (poderíamos aplicar heurísticas aqui)
    escolha = reais_positivas[0]
    resultado.update(escolha)
    resultado["ok"] = True

    return resultado


def processar_filtros_PB(lista_Hn):

    # INPUTS
    print("=== Parâmetros globais (serão usados para todas as seções) ===")
    Ra = float(input("Digite Ra (ohms): "))
    Rb = float(input("Digite Rb (ohms, diferente de 0): "))

    # verificações obs: fisicamente não da pra ter K negativo, por isso Ra < Rb
    assert Rb != 0, "Rb não pode ser zero."
    assert Rb != Ra, "Rb não pode ser igual a Ra (isso faria K = 0)."
    assert Ra < Rb, "Ra deve ser menor que Rb para garantir 0 < K < 1."

    K = abs(1.0 - (Ra / Rb))
    print(f"K calculado = {K:.6g}")

    C1_global = float(
        input("Digite valor de C1 (Farad) — usado nas seções 1ª e 2ª ordem: "))
    C2_global = float(
        input("Digite valor de C2 (Farad) — usado apenas nas seções 2ª ordem: "))

    resultados = []

    for i, H in enumerate(lista_Hn, start=1):
        # extrair num/den (vetores numpy)
        num = np.array(H.num[0][0], dtype=float)
        den = np.array(H.den[0][0], dtype=float)

        print(f"\n=== Processando filtro {i} (grau denom = {len(den)-1}) ===")

        if len(den) == 3:
            # 2ª ordem
            b0 = float(num[0]) if len(num) >= 1 else 0.0
            a1 = float(den[1])
            a0 = float(den[2])

            res = solve_scd_PB(b0, a1, a0, K, C1_global, C2_global)
            res.update({"filtro": i})
            resultados.append(res)

        elif len(den) == 2:
            # 1ª ordem
            b0 = float(num[0]) if len(num) >= 1 else 0.0
            a0 = float(den[1])

            res = solve_fst_PB(b0, a0, K, C1_global)
            res.update({"filtro": i})
            resultados.append(res)

        else:
            resultados.append({
                "filtro": i,
                "error": f"Denominador com grau {len(den)-1} não suportado por este solver."
            })

    return resultados
