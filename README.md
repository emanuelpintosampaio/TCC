## Como executar no Google Colab

Você pode rodar este projeto diretamente no Google Colab, sem instalar nada localmente.

### Passos:

1. Acesse: [https://colab.research.google.com/](https://colab.research.google.com/)
2. Clique em **"Novo notebook"**
3. Clique em **"Conectar"**
4. Cole o código abaixo em uma célula e execute:

```python
import os
from IPython.display import Image, display

if not os.path.exists('IC'):
    !git clone https://github.com/offemanuel/IC.git
%cd IC

!pip install -q numpy matplotlib

# Menu Principal
print("\n")
print("1 - Sistema de EDO completo")
print("2 - Simplificações exponenciais")
print("3 - Diferença: Sistema de EDO Completo vs Simplificações Exponenciais")

opcao_principal = input("Digite o número da opção desejada: ")

if opcao_principal == "1":
    print("\n Sistema de EDO Completo")
    print("1 - Subcycling")
    print("2 - Passo adaptativo (PID)")
    print("3 - Passo fixo\n")
    
    opcao_secundaria = input("Digite a opção para EDO Completo: ")
    if opcao_secundaria == "1":
        !python IC2_completo_RK3_subcycling.py
    elif opcao_secundaria == "2":
        !python IC2_completo_RK3_passo_adaptativo_PID.py
    elif opcao_secundaria == "3":
        !python IC2_completo_RK3_dt_fixo.py
    else:
        print("Opção inválida!")

elif opcao_principal == "2":
    print("\n Simplificações Exponenciais")
    print("1 - Subcycling")
    print("2 - Passo adaptativo (PID)")
    print("3 - Passo fixo\n")
    
    opcao_secundaria = input("Digite a opção para Simplificações: ")
    if opcao_secundaria == "1":
        !python IC2_simplificacao_RK3_subcycling.py
    elif opcao_secundaria == "2":
        !python IC2_simplificacao_RK3_passo_adaptativo_PID.py
    elif opcao_secundaria == "3":
        !python IC2_simplificacao_RK3_dt_fixo.py
    else:
        print("Opção inválida!")

elif opcao_principal == "3":
    print("\n Executando Comparação (Completo vs Simplificações)")
    !python IC2_RK3_all.py

else:
    print("Opção do menu principal inválida!")

# Gráficos
arquivos = ["grafico_subcycling.png","grafico_adaptativo.png","grafico_fixo.png","grafico_completo.png","grafico_all.png"]

for arq in arquivos:
    if os.path.exists(arq):
        print(f"Exibindo: {arq}")
        display(Image(arq))
```

---
