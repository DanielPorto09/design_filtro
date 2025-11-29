========================
Antes de rodar o projeto 
========================
  1. Instale as seguintes bibliotecas Python:

    pip install numpy
    pip install matplotlib
    pip install sympy
    pip install control

  2. Mude o path para o lugar no seu computador nos Arquivos:

    comp_eletronicos_pa.py
    comp_eletronicos_pb.py
    main.py
    tratamento_TF.py

    path atual:
    sys.path.append(r"C:\Users\danie\OneDrive\Área de Trabalho\Facul\PDS\P2\codigos\design_filtro")


===================================
Arquivos do Projeto e Suas Funções
===================================

  📄 Chebychev.py

    Implementa TODAS as funções matemáticas do filtro Chebyshev, inclui:

      cálculo de ripple (ε)
      cálculo da ordem mínima
      cálculo dos polos
      montagem dos numeradores e denominadores
      funções para Passa-Baixa (PB) e Passa-Alta (PA)

  📄 comp_eletronicos_pb.py

    Implementa funções responsáveis por calcular componentes eletrônicos reais (R, C) para filtros passa-baixa
    e decompõe H(s) para ordem ordem inferiores, inclui:

      Projeto para 1ª ordem (R1, R2, C1),
      Projeto automático para 2ª ordem usando SVF,
      Função process_filter_list() 

  📄 comp_eletronicos_pa.py

    Semelhante ao arquivo anterior, mas para filtros passa-alta, inclui:

      Projeto de 1ª ordem PA,
      Projeto automático SVF para 2ª ordem,
      process_filter_list_pa()

  📄 graphs.py

    Contém todas as funções de visualização, inclui:

      Diagrama de Bode,
      Root Locus,
      Curvas de magnitude para verificação das especificações,
      Impressão de cada filtro da cascata,
      Impressão dos valores de componentes calculados.

  📄 tratamento_TF.py

    Arquivo responsável por manipular a função de transferência completa, incluindo:

      cálculo de quantos filtros de 1ª e 2ª ordem são necessários,
      decomposição da FT em seções menores (cascateamento),
      criação dos filtros individuais (TFs),
      adição de zeros nos filtros PA 

  📄 main.py

    Arquivo principal (ponto de entrada).

    Fluxo do código:

    1.Pergunta ao usuário PB ou PA;
    2.Pede as especificações (frequências, atenuações);
    3.Calcula a FT Chebyshev;
    4.Decompõe em filtros menores;
    5.Calcula os componentes eletrônicos;
    6.Plota gráficos;
    7.Exibe relatórios.