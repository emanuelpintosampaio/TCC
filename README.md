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

SIMULACOES = {
    "ECONTROL": {
        "descricao": "EDO completo | passo fixo",
        "script":    "IC2_completo_RK3_dt_fixo_{raio}.py",
        "grafico":   "grafico_dt_fixo_completo_{raio}.png",
    },
    "E1": {
        "descricao": "EDO completo | controle por subpassos",
        "script":    "IC2_completo_RK3_subpassos_{raio}.py",
        "grafico":   "grafico_subpassos_completo_{raio}.png",
    },
    "E2": {
        "descricao": "EDO simples | passo fixo",
        "script":    "IC2_simplificacao_RK3_dt_fixo_{raio}.py",
        "grafico":   "grafico_dt_fixo_simplificacao_{raio}.png",
    },
    "E3": {
        "descricao": "EDO simples | controle por subpassos",
        "script":    "IC2_simplificacao_RK3_subpassos_{raio}.py",
        "grafico":   "grafico_subpassos_simplificacao_{raio}.png",
    }}

OPCOES_RAIO = {"1": "30", "2": "500"}

while True:
    # Seleção de raio
    print("\nESCOLHA O RAIO DA GOTA PARA A SIMULAÇÃO")
    for k, v in OPCOES_RAIO.items():
        print(f"{k} - {v} µm{' ' * (41 - len(v))}")

    opcao_raio = input("Opção: ").strip()
    if opcao_raio not in OPCOES_RAIO:
        print("Opção de raio inválida! Tente novamente.")
        continue
    raio = OPCOES_RAIO[opcao_raio]
    time.sleep(0.5)

    # Seleção da simulação 
    print("\nESCOLHA O TIPO DE SIMULAÇÃO")
    opcoes_sim = list(SIMULACOES.keys())
    for i, nome in enumerate(opcoes_sim, 1):
        desc = SIMULACOES[nome]["descricao"]
        linha = f"{i} - {nome}: {desc}"
        print(f"{linha:<58}")
    print(f"{len(opcoes_sim)+1} - Comparação: todas as simulações{' ' * 25}")

    opcao_sim = input("Opção: ").strip()

    graficos_para_exibir = []

    if opcao_sim.isdigit() and 1 <= int(opcao_sim) <= len(opcoes_sim):
        nome_sim = opcoes_sim[int(opcao_sim) - 1]
        sim = SIMULACOES[nome_sim]
        script  = sim["script"].format(raio=raio)
        grafico = sim["grafico"].format(raio=raio)

        print(f"\nExecutando Simulação {nome_sim} — {sim['descricao']} | Raio: {raio} µm\n")
        os.system(f"python {script}")
        graficos_para_exibir.append(grafico)

    elif opcao_sim == str(len(opcoes_sim) + 1):
        print(f"\nExecutando comparação completa | Raio: {raio} µm\n")
        os.system(f"python IC2_RK3_all_{raio}.py")
        graficos_para_exibir = [sim["grafico"].format(raio=raio) for sim in SIMULACOES.values()] + [f"grafico_all_{raio}.png"]

    else:
        print("Opção inválida! Tente novamente.")
        continue

    # Gráficos 
    print("\nBUSCANDO GRÁFICOS GERADOS...")
    exibidos = 0
    for arq in graficos_para_exibir:
        if os.path.exists(arq):
            print(f"Exibindo: {arq}")
            display(Image(arq))
            exibidos += 1
    if exibidos == 0:
        print("Nenhum gráfico encontrado para exibir.")
    else:
        time.sleep(3)

    # Continuar? 
    if input("\nDeseja fazer outra simulação? (s/n): ").strip().lower() != 's':
        print("\nEncerrando o programa...")
        break
    
```

---
