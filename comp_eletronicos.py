# comp_eletronicos_corrected.py - funções de separação e dessíntese Sallen-Key
from Chebychev import cheb_transfer, cheb_order, cheb_order_pa
import numpy as np
import control as ctrl
import sympy as sp


def n_filtros_PB(a_pass, a_stop, w_pass, w_stop):
    n, trash = cheb_order(a_pass, a_stop, w_pass, w_stop)
    if n % 2 == 1:
        k_filters_fst_order = 1
        k_filters_sec_order = (n - 1) // 2
    else:
        k_filters_fst_order = 0
        k_filters_sec_order = n // 2
    return k_filters_fst_order, k_filters_sec_order


def n_filtros_PA(a_pass, a_stop, w_pass, w_stop):
    n, trash = cheb_order_pa(a_pass, a_stop, w_pass, w_stop)
    if n % 2 == 1:
        k_filters_fst_order = 1
        k_filters_sec_order = (n - 1) // 2
    else:
        k_filters_fst_order = 0
        k_filters_sec_order = n // 2
    return k_filters_fst_order, k_filters_sec_order


def separar_ordens(poles, tol_real=1e-8, tol_pair=1e-6):
    poles = np.array(poles, dtype=complex)
    usados = np.zeros(len(poles), dtype=bool)
    polos_1ord = []
    polos_2ord = []

    for i, p in enumerate(poles):
        if usados[i]:
            continue
        if np.isclose(p.imag, 0.0, atol=tol_real, rtol=0.0):
            polos_1ord.append(p)
            usados[i] = True
            continue
        conj_target = np.conjugate(p)
        diffs = np.abs(poles - conj_target)
        cand_idx = [j for j in range(len(poles)) if (
            not usados[j]) and j != i and diffs[j] <= tol_pair]
        if cand_idx:
            j = cand_idx[0]
            if p.imag > 0:
                polos_2ord.append((p, poles[j]))
            else:
                polos_2ord.append((poles[j], p))
            usados[i] = True
            usados[j] = True
        else:
            polos_1ord.append(p)
            usados[i] = True
    return polos_1ord, polos_2ord


def gerar_secoes_normalizadas(polos_real, polos_complex):
    """Gera TransferFunctions monic (numerador = 1) para cada seção a partir dos polos."""
    secoes = []

    # 1ª ORDEM
    for p in polos_real:
        den = np.real_if_close(np.poly([p]))
        den = np.array(den, dtype=float)
        if abs(den[0]) < 1e-30:
            continue
        den = den / den[0]
        num = np.array([1.0])
        H = ctrl.TransferFunction(num, den)
        secoes.append(H)

    # 2ª ORDEM
    for p, pc in polos_complex:
        den = np.real_if_close(np.poly([p, pc]))
        den = np.array(den, dtype=float)
        if abs(den[0]) < 1e-30:
            continue
        den = den / den[0]
        num = np.array([1.0])
        H = ctrl.TransferFunction(num, den)
        secoes.append(H)

    return secoes


def separa_func(num, den, k_filters_fst_order, k_filters_sec_order, w_c):
    """Decompõe H(s) em seções monic 1ª e 2ª ordem e retorna lista de TransferFunction."""
    poles = np.roots(den)
    polos_real, polos_complex = separar_ordens(poles)

    if k_filters_fst_order is not None and len(polos_real) != k_filters_fst_order:
        print(
            f"Warning: esperado {k_filters_fst_order} polos 1ª ordem, encontrado {len(polos_real)}.")
    if k_filters_sec_order is not None and len(polos_complex) != k_filters_sec_order:
        print(
            f"Warning: esperado {k_filters_sec_order} polos 2ª ordem, encontrado {len(polos_complex)}.")

    secoes = gerar_secoes_normalizadas(polos_real, polos_complex)
    return secoes


R_LIMIT_DEFAULT = 1e6


def solve_fst_PB(b0, a0, K, C1, r_limit=R_LIMIT_DEFAULT):
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

    escolha = min(reais_positivas, key=lambda s: max(s["R1"], s["R2"]))
    resultado.update(escolha)
    resultado["ok"] = True
    resultado["raw_solutions"] = reais_positivas
    return resultado


def processar_filtros_PB(lista_Hn, K_total=1.0):
    print("=== Parâmetros globais (serão usados para todas as seções) ===")
    Ra = float(input("Digite Ra (ohms): "))
    Rb = float(input("Digite Rb (ohms, diferente de 0): "))

    assert Ra != 0.0, "Ra não pode ser zero para configuração não-inversora."
    assert Rb != 0.0, "Rb não pode ser zero."

    K_sk = 1.0 + (Rb / Ra)
    print(f"K do bloco Sallen-Key (não-inversor) = {K_sk:.6g}")

    K_effective = float(K_total) * float(K_sk)
    print(f"K efetivo (K_total * K_sk) = {K_effective:.6g}")

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

        if deg >= 1:
            a0 = den[-1]
            b0_atual = num[0] if len(num) > 0 else 0.0
            b0_desejado = K_effective * a0
            if abs(b0_atual) < 1e-12:
                num = np.array([b0_desejado])
                print(
                    f"  Numerador definido para {b0_desejado:.6g} (b0_atual ~ 0)")
            else:
                fator = b0_desejado / b0_atual
                if not np.isfinite(fator) or abs(fator) > 1e8:
                    print(
                        f"  Atenção: fator de escala muito grande ({fator:.6g}). Verifique C1/C2 ou w_c.")
                else:
                    if abs(fator - 1.0) > 1e-9:
                        num = num * fator
                        print(
                            f"  Numerador escalado por fator {fator:.6f} para compatibilidade")

        if deg == 2:
            b0 = float(num[0]) if len(num) >= 1 else 0.0
            a1 = float(den[1])
            a0 = float(den[2])
            res = solve_scd_PB(b0, a1, a0, K_effective, C1_global, C2_global)
            res.update({"filtro": i})
            resultados.append(res)

        elif deg == 1:
            b0 = float(num[0]) if len(num) >= 1 else 0.0
            a0 = float(den[1])
            res = solve_fst_PB(b0, a0, K_effective, C1_global)
            res.update({"filtro": i})
            resultados.append(res)

        else:
            resultados.append(
                {"filtro": i, "error": f"Denominador com grau {deg} não suportado."})

    return resultados
