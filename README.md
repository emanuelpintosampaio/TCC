## Como executar no Google Colab

Você pode rodar este projeto diretamente no Google Colab, sem instalar nada localmente.

### Passos:

1° Acesse: [https://colab.research.google.com/](https://colab.research.google.com/)

2° Clique em **"Novo notebook"**

3° Clique em **"Conectar"**

4° Cole o código abaixo em uma célula e execute:

```python
import os
import time
from IPython.display import Image, display, clear_output

if not os.path.exists('IC'):
    !git clone https://github.com/offemanuel/IC.git
%cd IC

!pip install -q numpy matplotlib

# Menu 
while True:
    print("ESCOLHA O RAIO DA GOTA PARA A SIMULAÇÃO:")
    print("1 - 30 µm")
    print("2 - 500 µm")
    
    opcao_raio = input("\nDigite o número da opção desejada: ").strip()
    
    if opcao_raio == "1":
        raio = "30"
        print(f"\nRaio selecionado: {raio} µm")
    elif opcao_raio == "2":
        raio = "500"
        print(f"\nRaio selecionado: {raio} µm")
    else:
        print("Opção de raio inválida! Tente novamente.\n")
        continue
    
    time.sleep(1)
    print("1 - Sistema de EDO completo")
    print("2 - Simplificações exponenciais")
    print("3 - Diferença: Sistema de EDO Completo vs Simplificações Exponenciais")
    
    opcao_principal = input("\nDigite o número da opção desejada: ").strip()
    
    if opcao_principal == "1":
        print(f"\n Sistema de EDO Completo | Raio: {raio} µm")
        print("1 - Subpassos")
        #print("3 - Passo adaptativo (PID)")
        print("2 - Passo fixo\n")
        
        opcao_secundaria = input("Digite a opção para EDO Completo: ").strip()
        if opcao_secundaria == "1":
            !python IC2_completo_RK3_subpassos_{raio}.py
        elif opcao_secundaria == "3":
            !python IC2_completo_RK3_passo_adaptativo_PID_{raio}.py
        elif opcao_secundaria == "2":
            !python IC2_completo_RK3_dt_fixo_{raio}.py
        else:
            print("Opção inválida!")
    
    elif opcao_principal == "2":
        print(f"\n Simplificações Exponenciais | Raio: {raio} µm")
        print("1 - Subpassos")
        #print("3 - Passo adaptativo (PID)")
        print("2 - Passo fixo\n")
        
        opcao_secundaria = input("Digite a opção para Simplificações: ").strip()
        if opcao_secundaria == "1":
            !python IC2_simplificacao_RK3_subpassos_{raio}.py
        elif opcao_secundaria == "3":
            !python IC2_simplificacao_RK3_passo_adaptativo_PID_{raio}.py
        elif opcao_secundaria == "2":
            !python IC2_simplificacao_RK3_dt_fixo_{raio}.py
        else:
            print("Opção inválida!")
    
    elif opcao_principal == "3":
        print(f"\n Executando Comparação (Completo vs Simplificações) | Raio: {raio} µm")
        !python IC2_RK3_all_{raio}.py
    
    else:
        print("Opção do menu principal inválida!")
        continue  

    print("BUSCANDO GRÁFICOS GERADOS...")
    
    arquivos = [
        f"grafico_subpassos_completo_{raio}.png",
        f"grafico_subpassos_simplificacao_{raio}.png",
        f"grafico_dt_adaptativo_completo_{raio}.png",
        f"grafico_dt_adaptativo_simplificacao_{raio}.png",
        f"grafico_dt_fixo_completo_{raio}.png",
        f"grafico_dt_fixo_simplificacao_{raio}.png",
        f"grafico_all_{raio}.png"
    ]

    graficos_exibidos = 0
    for arq in arquivos:
        if os.path.exists(arq):
            print(f"Exibindo: {arq}")
            display(Image(arq))
            graficos_exibidos += 1

    if graficos_exibidos == 0:
        print("Nenhum gráfico encontrado para exibir.")
    else:
        time.sleep(3)

    continuar = input("\nDeseja fazer outra simulação? (s/n): ").strip().lower()
    if continuar != 's':
        print("\nEncerrando o programa...")
        break
    
```

---
