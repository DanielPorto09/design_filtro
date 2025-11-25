import matplotlib.pyplot as plt
import control as ctrl
import numpy as np

# -------------------------------------------------------
#  DIAGRAMAS DE ANALISE DA FT ORIGINAL
# -------------------------------------------------------

def anlize_em_frequencia(H, f_pass, f_stop, a_pass, a_stop, w_pass, w_stop):
    """
    Recebe parametros do filtro e plota root-locus, diagrama de bode e frequencia 
    de corte e passagem
    """

    # DC, banda passante, banda de corte
    w_test = np.array([1e-6, w_pass, w_stop])
    freq_names = ['DC', 'Frequência de passagem', 'Frequência de corte']

    for w, name in zip(w_test, freq_names):
        # Usar freqresp para obter a resposta em frequência em pontos específicos
        response = H(1j * w)
        mag = np.abs(response)
        freq_hz = w / (2*np.pi)
        # print(f"{name} ({freq_hz:.1f} Hz): {20*np.log10(mag):.2f} dB")

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
    w_plot = np.logspace(np.log10(100), np.log10(10000),
                         1000)  # 100 a 10000 rad/s
    mag_plot, phase_plot, omega_plot = ctrl.bode(H, w_plot, plot=False)

    plt.semilogx(omega_plot/(2*np.pi), 20 *
                 np.log10(mag_plot), 'b-', linewidth=2)
    plt.axvline(f_pass, color='r', linestyle='--',
                label=f'f_pass = {f_pass} Hz')
    plt.axvline(f_stop, color='g', linestyle='--',
                label=f'f_stop = {f_stop} Hz')
    plt.axhline(-a_pass, color='r', linestyle=':',
                label=f'Ripple = {a_pass} dB')
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


# ==================================================
# EXIBE FUNÇÕES SEPARADAS PRA CADA FILTRO DA CASCATA
# ==================================================

def PrintaPlota_lista_de_filtros(filtros):
    # PRINTA  CADA FILTRO
    print("=== Lista de filtros encontrados ===")
    for i, f in enumerate(filtros):
        print(f"\nFiltro {i+1}:")
        print(f)

# PLOTAR DIAGRAMA DE BODE DE CADA FILTRO
    for i, f in enumerate(filtros):
        mag, phase, omega = ctrl.bode(
            f, dB=True, Hz=False, deg=True, plot=False)

        plt.figure(figsize=(8, 6))
        ctrl.bode(f, dB=True, Hz=True, deg=True, label=f"Filtro {i+1}")
        plt.suptitle(f"Bode do Filtro {i+1}")
        plt.show()
    

# ================================
# EXIBE OS VALORES DOS COMPONETNES
# ================================

def print_resultados(resultados):
    """
    Recebe a lista de dicionários 'resultados' e imprime cada filtro formatado.
    """
    print("\n=== RESULTADOS DOS FILTROS ===")

    for i, data in enumerate(resultados):
        print(f"\n--- Filtro {i+1} ---")
        print(f"Tipo: {data['tipo']}ª ordem")
        print(f"R1 = {data['R1']}")
        print(f"R2 = {data['R2']}")
        print(f"C1 = {data['C1']}")
        print(f"C2 = {data['C2']}")
        print(f"Ra = {data['Ra']}")
        print(f"Rb = {data['Rb']}")
        print(f"K  = {data['K']}")

