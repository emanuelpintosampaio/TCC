## Como executar no Google Colab

Você pode rodar este projeto diretamente no Google Colab, sem instalar nada localmente.

### Passos:

1. Acesse: [https://colab.research.google.com/](https://colab.research.google.com/)
2. Clique em **"Novo notebook"**
3. Clique em **"Conectar"**
4. Cole o código abaixo em uma célula e execute:

```python
# Clonar o repositório
!git clone https://github.com/offemanuel/IC.git
%cd IC

# Instalar dependências
!pip install numpy matplotlib

# Configurar gráficos
%matplotlib inline

# Menu de escolha
print("Escolha qual código deseja executar:\n")
print("1 - Subcycling")
print("2 - Passo adaptativo (PID)")
print("3 - Passo fixo")
print("4 - Passo fixo + Subcycling")

opcao = input("Digite o número da opção: ")

if opcao == "1":
    !python IC2_completo_RK3_subcycling.py
elif opcao == "2":
    !python IC2_completo_RK3_passo_adaptativo_PID.py
elif opcao == "3":
    !python IC2_completo_RK3_dt_fixo.py
elif opcao == "4":
    !python IC2_completo_RK3_all.py
else:
    print("Opção inválida!")

# Exibir gráfico
from IPython.display import Image, display
import os

if os.path.exists("grafico.png"):
    display(Image("grafico.png"))
else:
    print("Gráfico não encontrado.")
```

---
