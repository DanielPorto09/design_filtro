# FUNÇÕES E BIBLIOTECAS PUBLICAS E AUTORAIS UTILIZADAS
from Chebychev import cheb_transfer, cheb_order, eps_from_ap, cheb_poles, cheb_transfer_pa, cheb_order_pa
import numpy as np
import matplotlib.pyplot as plt
import control as ctrl

# Caminho para encontrar pasta
import sys
sys.path.append(
    r"C:\Users\danie\OneDrive\Área de Trabalho\Facul\PDS\P2\codigos\design_filtro")


choice = int(input("Escolha 1 (PB) ou 2 (PA): "))


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
    print("\nH(T) NA FORMA PADRÃO", H)

    print("\nVerificação da resposta:")


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

    H = ctrl.TransferFunction(num, den)
    print("\nH(T) NA FORMA PADRÃO", H)


# -------------------------------------------------------
#  Rad/s -> Hz para melhor visualização de alguns trechos
# -------------------------------------------------------

w_test = np.array([1e-6, w_pass, w_stop])  # DC, banda passante, banda de corte
freq_names = ['DC', 'Frequência de passagem', 'Frequência de corte']

for w, name in zip(w_test, freq_names):
    # Usar freqresp para obter a resposta em frequência em pontos específicos
    response = H(1j * w)
    mag = np.abs(response)
    freq_hz = w / (2*np.pi)
    print(f"{name} ({freq_hz:.1f} Hz): {20*np.log10(mag):.2f} dB")

# -------------------------------------------------------
#  DIAGRAMAS DE ANALISE DO MÉTODO
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
plt.ylim(-a_stop-5, 5)

plt.tight_layout()
plt.show()
