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

# Menu de escolha
print("--- Simulador de Integração Numérica ---")
print("1 - Subcycling")
print("2 - Passo adaptativo (PID)")
print("3 - Passo fixo")
print("4 - Passo fixo + Passo adaptativo + Subcycling")
print("5 - RODAR TUDO (Sequencial)\n")

opcao = input("Digite o número da opção: ")

if opcao == "1":
    !python IC2_completo_RK3_subcycling.py
elif opcao == "2":
    !python IC2_completo_RK3_passo_adaptativo_PID.py
elif opcao == "3":
    !python IC2_completo_RK3_dt_fixo.py
elif opcao == "4":
    !python IC2_completo_RK3_all.py
elif opcao == "5":
    !python IC2_completo_RK3_subcycling.py
    !python IC2_completo_RK3_passo_adaptativo_PID.py
    !python IC2_completo_RK3_dt_fixo.py
    !python IC2_completo_RK3_all.py
else:
    print("Opção inválida!")

print("\n--- Resultados Gerados ---")
arquivos = [
    "grafico_subcycling.png",
    "grafico_adaptativo.png",
    "grafico_fixo.png",
    "grafico_completo.png"
]

for arq in arquivos:
    if os.path.exists(arq):
        print(f"Exibindo: {arq}")
        display(Image(arq))
```

---
