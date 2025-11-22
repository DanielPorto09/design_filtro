"""
design_sallen_key.py

Arquivo final que integra seu fluxo antigo com o solver Sallen-Key passa-baixa
compatível com seções Chebyshev (1ª e 2ª ordem).

Características:
- Mantém nomenclatura em português (solve_fst_PB, solve_scd_PB, processar_filtros_PB)
- Usa K dependente de Ra e Rb exatamente como antes: K = abs(1 - Ra/Rb)
- Apresenta aviso com o K padrão Sallen-Key (K_sk = 1 + Rb/Ra) — recomendado
- Aceita C1 e C2 arbitrários
- Limite opcional de resistores (r_limit, padrão 1e6)
- Entrada: lista de control.TransferFunction (lista_Hn), mesma estrutura que você já usava
- Saída: lista de dicionários por seção (formato compatível com seus resultados anteriores)

OBS: imagem da topologia (enviada pelo usuário) em:
/mnt/data/e039135a-2608-45ef-8960-2e53af5177e8.png

"""

import sys
import numpy as np
import control as ctrl
import sympy as sp


# -----------------------------
# Funções auxiliares antigas (mantidas)
# -----------------------------

def n_filtros_PB(a_pass, a_stop, w_pass, w_stop):
    n, trash = cheb_order(a_pass, a_stop, w_pass, w_stop)
    if n % 2 == 1:
        k_filters_fst_order = 1
        k_filters_sec_order = (n - 1) // 2
    else:
        k_filters_fst_order = 0
        k_filters_sec_order = n // 2
    return k_filters_fst_order, k_filters_sec_order


def separa_func(num, den, k_filters_fst_order, k_filters_sec_order):
    """Decompõe H(s) em seções de 1ª e 2ª ordem (retorna lista de control.TransferFunction)
    Mantive sua implementação principal (compatível com seu código anterior).
    """
    H_total = ctrl.TransferFunction(num, den)
    poles = np.roots(den)
    filtros = []
    usados = np.zeros(len(poles), dtype=bool)

    # 1ª ordem: polos reais
    for i in range(len(poles)):
        if usados[i]:
            continue
        p = poles[i]
        if np.isclose(p.imag, 0):
            den1 = [1, -p.real]
            num1 = [1]
            F = ctrl.TransferFunction(num1, den1)
            filtros.append(F)
            usados[i] = True
            if len([f for f in filtros if f.den[0][0].size == 2]) == k_filters_fst_order:
                break

    # 2ª ordem: pares complexos conjugados
    for i in range(len(poles)):
        if usados[i]:
            continue
        p = poles[i]
        j = np.where(
            np.isclose(poles.real,  p.real) &
            np.isclose(poles.imag, -p.imag) &
            (~usados)
        )[0]
        if len(j) == 0:
            continue
        j = j[0]
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


# -----------------------------
# Solvers atualizados (Sallen-Key)
# -----------------------------

R_LIMIT_DEFAULT = 1e6


def solve_fst_PB(b0, a0, K, C1, r_limit=R_LIMIT_DEFAULT):
    """
    Resolve seção de 1ª ordem (compatível com seu formato anterior).
    A TF esperada é H(s) = b0 / (a0*s + 1) (ou ajustada conforme normalização).

    R1 = 1 / (a0 * C1)
    """
    resultado = {"ordem": 1, "C1": C1, "K": K}

    if abs(a0) < 1e-30:
        resultado.update({"error": "a0 quase zero, não é possível resolver R1."})
        return resultado

    R1 = 1.0 / (a0 * C1)

    if not (0 < R1 < r_limit):
        resultado.update({"error": f"R1 encontrado fora do intervalo físico (R1={R1:.6g} Ω). Limite = {r_limit}"})
        return resultado

    resultado["R1"] = float(R1)

    if abs(b0 - K*a0) > max(1e-9, 1e-6 * abs(b0)):
        resultado["warning"] = f"Inconsistência: b0 != K*a0 (b0={b0:.6g}, K*a0={(K*a0):.6g})."
    else:
        resultado["ok"] = True

    return resultado


def solve_scd_PB(b0, a1, a0, K, C1, C2, r_limit=R_LIMIT_DEFAULT):
    """
    Novo solver para seção de 2ª ordem baseado na topologia Sallen-Key passa-baixa.

    Usa SymPy para resolver as equações:
      (1) 1/(R1*R2*C1*C2) = a0
      (2) 1/(R1*C1) + 1/(R2*C1) + (1-K)/(R2*C2) = a1
      (3) K/(R1*R2*C1*C2) = b0

    Retorna dicionário com R1,R2 ou erro.
    """
    resultado = {"ordem": 2, "C1": C1, "C2": C2, "K": K}

    if abs(a0) < 1e-30:
        resultado.update({"error": "a0 quase zero; impossível resolver (divisão por zero)."})
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

    # filtrar soluções reais e positivas e dentro do limite
    reais_positivas = []
    for (r1_val, r2_val) in sols:
        if np.isreal(r1_val) and np.isreal(r2_val) and r1_val > 0 and r2_val > 0 and r1_val < r_limit and r2_val < r_limit:
            reais_positivas.append({"R1": float(r1_val), "R2": float(r2_val)})

    if not reais_positivas:
        resultado.update({"error": "Nenhuma solução real positiva encontrada.", "raw_solutions": sols})
        return resultado

    # escolha heurística: menor resistor máximo
    escolha = min(reais_positivas, key=lambda s: max(s["R1"], s["R2"]))
    resultado.update(escolha)
    resultado["ok"] = True
    resultado["raw_solutions"] = reais_positivas

    return resultado


# -----------------------------
# Processamento principal (mantém API antiga)
# -----------------------------

def processar_filtros_PB(lista_Hn):
    """
    Interativo: pede Ra,Rb,C1,C2 e processa cada seção (lista de control.TransferFunction)
    Retorna lista de resultados (dicionários) — formato semelhante ao seu código anterior.
    """
    print("=== Parâmetros globais (serão usados para todas as seções) ===")
    Ra = float(input("Digite Ra (ohms): "))
    Rb = float(input("Digite Rb (ohms, diferente de 0): "))

    assert Rb != 0, "Rb não pode ser zero."

    # Mantemos cálculo original de K, mas mostramos também K padrão SK
    K = abs(1.0 - (Ra / Rb))
    K_sk = 1.0 + (Rb / Ra) if Ra != 0 else None
    print(f"K calculado (original) = {K:.6g}")
    if K_sk is not None:
        print(f"K padrão Sallen-Key (recomendado) = {K_sk:.6g}")

    C1_global = float(input("Digite valor de C1 (Farad) — usado nas seções 1ª e 2ª ordem: "))
    C2_global = float(input("Digite valor de C2 (Farad) — usado apenas nas seções 2ª ordem: "))

    resultados = []

    for i, H in enumerate(lista_Hn, start=1):
        # extrair num/den (vetores numpy) usando control.TransferFunction API
        num = np.array(H.num[0][0], dtype=float)
        den = np.array(H.den[0][0], dtype=float)
        deg = len(den) - 1

        print(f"\n=== Processando filtro {i} (grau denom = {deg}) ===")

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
            resultados.append({"filtro": i, "error": f"Denominador com grau {deg} não suportado por este solver."})

    return resultados


