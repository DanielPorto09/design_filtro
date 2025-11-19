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
    eps = eps_from_ap(a_pass)

    num = np.acosh(np.sqrt((10**(0.1*a_stop) - 1)/(10**(0.1*a_pass) - 1)))
    den = np.acosh(w_stop / w_pass)

    n = abs(np.ceil(num / den))
    n = int(n)
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

# -------------------------------------------------------
#  Cálculo dos polos do Chebyshev Tipo I CORRIGIDO
# -------------------------------------------------------


def cheb_poles(n, eps, w_c=1.0):  # CORREÇÃO: adicionado w_c com valor default
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

# -------------------------------------------------------
#  Função de transferência do Chebyshev CORRIGIDA
# -------------------------------------------------------


def cheb_transfer(n, eps, w_c):
    poles = cheb_poles(n, eps, w_c)  # Agora funciona com 3 argumentos

    # Constante do ganho DC
    if n % 2 == 0:  # n par
        G0 = 1 / np.sqrt(1 + eps**2)
    else:  # n ímpar
        G0 = 1

    # Construir polinômio do denominador
    den_poly = np.poly(poles)
    den = np.real(den_poly)

    # Ajustar o numerador para ganho DC correto
    num = [G0 * den[0]]

    return num, den, poles

# -------------------------------------------------------
#  PRINT FUNÇÃO DE TRANSFERÊNCIA
# -------------------------------------------------------


def print_transfer(num, den, poles):
    print("\n===============================")
    print(" Função de Transferência H(s) ")
    print("===============================\n")

    print("Numerador:")
    print(num, "\n")

    print("Denominador ( em Polos):")
    for i, p in enumerate(poles, start=1):
        print(f"Polo {i}:  {p.real:+.5f}  {p.imag:+.5f}j")


"""    print("\nCoeficientes do polinômio do denominador:")
    for i, c in enumerate(den):
        print(f"s^{len(den)-i-1}: {c}")
    print("\n")
"""
# -------------------------------------------------------
#  Aplicação do Chebychev do filtro
# -------------------------------------------------------


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

# CORREÇÃO: Usar w_pass como frequência de corte, não w_stop
num, den, poles = cheb_transfer(n, eps, w_pass)
print_transfer(num, den, poles)

H = ctrl.TransferFunction(num, den)
print("\nH(T) NA FORMA PADRÃO", H)

print("\nVerificação da resposta:")

# CORREÇÃO: Usar freqresp em vez de bode para obter valores numéricos
w_test = np.array([1e-6, w_pass, w_stop])  # DC, banda passante, banda de corte
freq_names = ['DC', 'Frequência de passagem', 'Frequência de corte']

for w, name in zip(w_test, freq_names):
    # Usar freqresp para obter a resposta em frequência em pontos específicos
    response = H(1j * w)
    mag = np.abs(response)
    freq_hz = w / (2*np.pi)
    print(f"{name} ({freq_hz:.1f} Hz): {20*np.log10(mag):.2f} dB")

# -------------------------------------------------------
#  DIAGRAMAS CORRIGIDOS
# -------------------------------------------------------

plt.figure(figsize=(13, 9))

# Root Locus
plt.subplot(2, 2, 1)
ctrl.rlocus(H)
plt.title("Root Locus")
plt.grid(True)

# Diagrama de Bode
plt.subplot(2, 2, 2)
# CORREÇÃO: usar plot=False (minúsculo)
mag_bode, phase_bode, omega_bode = ctrl.bode(H, dB=True, plot=False)
plt.semilogx(omega_bode, 20*np.log10(mag_bode))
plt.title("Diagrama de Bode - Magnitude")
plt.xlabel('Frequência (rad/s)')
plt.ylabel('Magnitude (dB)')
plt.grid(True)

# Diagrama de Fase
plt.subplot(2, 2, 3)
plt.semilogx(omega_bode, phase_bode)
plt.title("Diagrama de Bode - Fase")
plt.xlabel('Frequência (rad/s)')
plt.ylabel('Fase (graus)')
plt.grid(True)

# Resposta em frequência com pontos críticos (em Hz)
plt.subplot(2, 2, 4)
w_plot = np.logspace(np.log10(100), np.log10(10000), 1000)  # 100 a 10000 rad/s
mag_plot, phase_plot, omega_plot = ctrl.bode(H, w_plot, plot=False)

plt.semilogx(omega_plot/(2*np.pi), 20*np.log10(mag_plot), 'b-', linewidth=2)
plt.axvline(f_pass, color='r', linestyle='--', label=f'f_pass = {f_pass} Hz')
plt.axvline(f_stop, color='g', linestyle='--', label=f'f_stop = {f_stop} Hz')
plt.axhline(-a_pass, color='r', linestyle=':', label=f'Ripple = {a_pass} dB')
plt.axhline(-a_stop, color='g', linestyle=':',
            label=f'Atenuação = {a_stop} dB')

# Marcar os pontos de verificação
for w, name in zip(w_test, freq_names):
    response = H(1j * w)
    mag_point = np.abs(response)
    freq_hz = w / (2*np.pi)
    plt.plot(freq_hz, 20*np.log10(mag_point), 'ro', markersize=8)

plt.xlabel('Frequência (Hz)')
plt.ylabel('Magnitude (dB)')
plt.title('Resposta em Frequência - Filtro Chebyshev')
plt.legend()
plt.grid(True)
plt.ylim(-a_stop-5, 5)  # Ajustar limites para melhor visualização

plt.tight_layout()
plt.show()
