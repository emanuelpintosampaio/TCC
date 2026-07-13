## Como executar o menu das simulações no Google Colab

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

# Clona o repositório na primeira execução e acessa a pasta dos códigos.
if not os.path.exists('TCC'):
    !git clone https://github.com/emanuelpintosampaio/TCC.git
    %cd TCC/CODE
    !pip install -q numpy matplotlib
elif os.path.basename(os.getcwd()) != 'CODE':
    if os.path.exists('TCC/CODE'):
        %cd TCC/CODE
    elif os.path.exists('CODE'):
        %cd CODE

# Associação entre cada simulação, o script correspondente e o gráfico gerado.
SIMULACOES = {
    "SCONTROL": {
        "descricao": "Sistema Completo | SSPRK3 | Passo Fixo",
        "script": "SCONTROL_{raio}.py",
        "grafico": "SCONTROL_{raio}.png",
    },
    "ECONTROL": {
        "descricao": "Sistema Completo | Euler Explícito | Passo Fixo",
        "script": "ECONTROL_{raio}.py",
        "grafico": "ECONTROL_{raio}.png",
    },
    "E1": {
        "descricao": "Sistema Completo | Euler Explícito | Subpassos",
        "script": "E1_{raio}.py",
        "grafico": "E1_{raio}.png",
    },
    "E2": {
        "descricao": "Sistema Completo | Euler Explícito | PID",
        "script": "E2_{raio}.py",
        "grafico": "E2_{raio}.png",
    },
    "E3": {
        "descricao": "Sistema Simplificado | Euler Explícito | Passo Fixo",
        "script": "E3_{raio}.py",
        "grafico": "E3_{raio}.png",
    },
    "E4": {
        "descricao": "Sistema Simplificado | Euler Explícito | Subpassos",
        "script": "E4_{raio}.py",
        "grafico": "E4_{raio}.png",
    },
    "E5": {
        "descricao": "Sistema Simplificado | Euler Explícito | PID",
        "script": "E5_{raio}.py",
        "grafico": "E5_{raio}.png",
    }
}

OPCOES_RAIO = {"1": "30", "2": "500"}

# Menu principal.
while True:
    clear_output(wait=True)

    print("=" * 60)
    print(" ESCOLHA O RAIO DA GOTA PARA A SIMULAÇÃO".center(60))
    print("=" * 60)

    for k, v in OPCOES_RAIO.items():
        print(f" [{k}] - {v} µm")

    print("-" * 60)

    time.sleep(0.5)
    opcao_raio = input("Opção de Raio (1 ou 2): ").strip()

    if opcao_raio not in OPCOES_RAIO:
        print("\nOpção de raio inválida. Tente novamente.")
        time.sleep(1.5)
        continue

    raio = OPCOES_RAIO[opcao_raio]

    print("\n" + "=" * 60)
    print(" ESCOLHA O TIPO DE SIMULAÇÃO".center(60))
    print("=" * 60)

    opcoes_sim = list(SIMULACOES.keys())

    for i, nome in enumerate(opcoes_sim, 1):
        print(f" [{i}] - {nome:<8} : {SIMULACOES[nome]['descricao']}")

    opcao_todas = len(opcoes_sim) + 1
    print(f" [{opcao_todas}] - Comparação : Executar todas as simulações")
    print("-" * 60)

    time.sleep(0.5)
    opcao_sim = input("Opção de Simulação: ").strip()

    # Lista dos gráficos que serão exibidos ao final da execução.
    graficos_para_exibir = []

    if opcao_sim.isdigit() and 1 <= int(opcao_sim) <= len(opcoes_sim):
        nome_sim = opcoes_sim[int(opcao_sim) - 1]
        sim = SIMULACOES[nome_sim]

        script = sim["script"].format(raio=raio)
        grafico = sim["grafico"].format(raio=raio)

        print(f"\nExecutando Simulação {nome_sim} | Raio: {raio} µm")
        print(f"Arquivo: {script}")

        os.system(f"python {script}")
        graficos_para_exibir.append(grafico)

    elif opcao_sim == str(opcao_todas):
        print(f"\nExecutando comparação completa | Raio: {raio} µm")

        script_all = f"all_{raio}.py"
        print(f"Arquivo: {script_all}")

        os.system(f"python {script_all}")

        # O script all gera dois gráficos comparativos.
        graficos_para_exibir.append(f"all_fixo_sub_{raio}.png")
        graficos_para_exibir.append(f"all_pid_{raio}.png")

    else:
        print("\nOpção inválida. Tente novamente.")
        time.sleep(1.5)
        continue

    # Aguarda a escrita dos arquivos antes de exibi-los.
    print("\nBuscando gráficos gerados...")
    time.sleep(0.5)

    exibidos = 0
    for arq in graficos_para_exibir:
        if os.path.exists(arq):
            print(f"Exibindo: {arq}")
            display(Image(arq))
            exibidos += 1

    if exibidos == 0:
        print("\nNenhum gráfico encontrado. Verifique se o script foi executado corretamente.")

    time.sleep(0.5)

    if input("\nDeseja fazer outra simulação? (s/n): ").strip().lower() != "s":
        print("\nEncerrando o programa.")
        break
    
```

---
