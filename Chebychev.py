import numpy as np
import matplotlib.pyplot as plt
import sympy as sp
import control as ctrl

# -------------------------------------------------------
#  Funções auxiliares
# -------------------------------------------------------


def eps_from_ap(a_pass):
    """
    Calcula epsilon a partir do ripple de passagem em dB.
    """
    return np.sqrt(10**(0.1 * a_pass) - 1)


def cheb_order(a_pass, a_stop, w_pass, w_stop):
    """
    Calcula a ordem mínima do filtro Chebyshev Tipo I.
    - a_pass e a_stop em dB
    - w_pass e w_stop em rad/s (ou equivalentes, pois é normalizado)
    """
    # Teste de erro
    assert a_pass < a_stop, "\n===================================\nERRO, a_pass deve ser menor que a_stop\n==================================="
    assert w_pass < w_stop, "Para Chebyshev PB, w_pass < w_stop"

    eps = eps_from_ap(a_pass)

    num = np.acosh(np.sqrt((10**(0.1*a_stop) - 1)/(10**(0.1*a_pass) - 1)))
    den = np.acosh(w_stop / w_pass)

    n = abs(np.ceil(num / den))
    n = int(n)
    return n, eps


def cheb_order_pa(a_pass, a_stop, w_pass, w_stop):
    """
    Ordem de Chebyshev Tipo I para passa-alta.
    """
    # Transformação passa-alta → passa-baixa equivalente
    # LP equivalents
    wp_lp = 1
    ws_lp = w_pass / w_stop   # precisa ser > 1

    assert ws_lp > 1, "\n===================================\nERRO: para Chebyshev HP, w_pass deve ser MAIOR que w_stop \n==================================="

    eps = eps_from_ap(a_pass)

    num = np.acosh(np.sqrt((10**(0.1*a_stop)-1)/(10**(0.1*a_pass)-1)))
    den = np.acosh(ws_lp)

    n = int(np.ceil(num / den))
    return n, eps


def chebyshev_poly(n, x):
    """
    Calcula o polinômio de Chebyshev de 1a espécie:
    - para |x| <= 1:  Cn(x) = cos(n * arccos(x))
    - para |x|  > 1:  Cn(x) = cosh(n * arccosh(x))
    """
    x = np.atleast_1d(x)

    C = np.zeros_like(x, dtype=float)

    # Região |x| <= 1
    mask1 = np.abs(x) <= 1
    C[mask1] = np.cos(n * np.arccos(x[mask1]))

    # Região |x| > 1
    mask2 = np.abs(x) > 1
    C[mask2] = np.cosh(n * np.acosh(x[mask2]))

    return C if len(C) > 1 else C[0]


def cheb_poles(n, eps, w_c=1.0): 
    """
    Calcula os polos do Chebyshev Tipo I com escalonamento de frequência.
    """
    beta = np.arcsinh(1/eps) / n  # CORREÇÃO: usar arcsinh em vez de acosh
    poles = []

    for k in range(1, n+1):
        theta = (2*k - 1) * np.pi / (2*n)
        sigma = -np.sinh(beta) * np.sin(theta)
        omega = np.cosh(beta) * np.cos(theta)
        # Escalonar para a frequência de corte desejada
        pole = w_c * (sigma + 1j*omega)
        poles.append(pole)

    return np.array(poles)


def cheb_transfer(n, eps, w_c):
    """
    Gera H(s) Chebyshev Tipo I passa-baixa conforme fórmulas analíticas.
    Corrige os casos de n par e n ímpar.
    """

    # ----------------------------
    # 1) Calcula polos normalizados
    # ----------------------------
    poles = cheb_poles(n, eps, w_c)

    # ----------------------------
    # 2) Agrupa polos em reais e complexos conjugados
    # ----------------------------
    real_poles = []
    complex_pairs = []

    usados = np.zeros(len(poles), dtype=bool)
    for i, p in enumerate(poles):
        if usados[i]:
            continue

        if abs(p.imag) < 1e-12:
            real_poles.append(p.real)
            usados[i] = True
        else:
            # procura conjugado
            for j in range(i+1, len(poles)):
                if usados[j]:
                    continue
                if abs(poles[j].real - p.real) < 1e-9 and \
                   abs(poles[j].imag + p.imag) < 1e-9:
                    complex_pairs.append((p, poles[j]))
                    usados[i] = usados[j] = True
                    break

    # ----------------------------
    # 3) Montagem do denominador
    # ----------------------------
    den = [1.0]

    if n % 2 == 1:
        # termo linear extra: (s + sinh(D))
        beta = np.arcsinh(1/eps) / n
        D = np.sinh(beta)
        den = np.convolve(den, [1, D*w_c])  # escalonado
    else:
        D = None

    # adiciona pares quadráticos
    for p, pc in complex_pairs:
        alpha = p.real
        beta = p.imag
        den = np.convolve(den, [1, -2*alpha, alpha*alpha + beta*beta])

    # adiciona polos reais ímpares (se houver)
    if n % 2 == 0:
        for pr in real_poles:
            den = np.convolve(den, [1, -pr])

    den = np.real_if_close(den)

    # ----------------------------
    # 4) Montagem do numerador
    # ----------------------------
    if n % 2 == 0:
        # n par: fórmula 10^(0.05*ap) * prod(B2m)
        G0 = 1 / np.sqrt(1 + eps*eps)
        num = [G0 * den[0]]
    else:
        # n ímpar: sinh(D) * prod(B2m)
        beta = np.arcsinh(1/eps) / n
        D = np.sinh(beta)

        num = [D * den[0]]

    return num, den.tolist(), poles


def cheb_transfer_pa(n, eps, w_c):
    """
    Chebyshev Tipo I – Passa-Alta (PA)
    Corrigido para tratar n par e n ímpar conforme a topologia correta.
    """
    # ======================================================
    # 1) Polos passa-baixa normalizados (wc = 1)
    # ======================================================
    poles_lp = cheb_poles(n, eps, 1.0)   

    # ======================================================
    # 2) Transformação LP → HP   (p_hp = wc / p_lp)
    # ======================================================
    poles_hp = np.array([w_c / p for p in poles_lp])

    # ======================================================
    # 3) Agrupar polos hp em reais e pares complexos
    # ======================================================
    usados = np.zeros(len(poles_hp), dtype=bool)
    real_poles = []
    complex_pairs = []

    for i, p in enumerate(poles_hp):
        if usados[i]:
            continue

        if abs(p.imag) < 1e-12:
            real_poles.append(p.real)
            usados[i] = True
        else:
            for j in range(i+1, len(poles_hp)):
                if usados[j]:
                    continue

                if abs(p.real - poles_hp[j].real) < 1e-9 and \
                   abs(p.imag + poles_hp[j].imag) < 1e-9:

                    complex_pairs.append((p, poles_hp[j]))
                    usados[i] = usados[j] = True
                    break

    # ======================================================
    # 4) Construção do denominador HP
    # ======================================================
    den = [1.0]

    # ----- termo linear extra (n ímpar)
    if n % 2 == 1:
        beta = np.arcsinh(1/eps) / n
        D = np.sinh(beta)
        # termo HP correspondente a (s + D) do PB:
        # (s + D) → (s + w_c*D)
        den = np.convolve(den, [1.0, -w_c * D])
    else:
        D = None   # só para manter compatível

    # ----- pares complexos transformados
    for p, pc in complex_pairs:
        alpha = p.real
        beta = p.imag
        den = np.convolve(den, [1.0, -2*alpha, alpha*alpha + beta*beta])

    # ----- polos reais restantes (se existirem)
    for pr in real_poles:
        den = np.convolve(den, [1.0, -pr])

    den = np.real_if_close(den).tolist()

    # ======================================================
    # 5) Numerador HP — obrigatório s^n
    # ======================================================
    if n % 2 == 0:
        # ganho DC equivalente ao PB espelhado
        K = 1.0 / np.sqrt(1 + eps**2)
    else:
        beta = np.arcsinh(1/eps) / n
        D = np.sinh(beta)
        K = D

    # numerador HP = K * s^n
    num = [K] + [0]*n

    return num, den, poles_hp
