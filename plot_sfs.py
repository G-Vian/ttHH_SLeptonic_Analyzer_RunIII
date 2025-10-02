import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import os
import sys

def make_sf_histograms(process_name):
    """
    Lê o arquivo .csv com os dados de SF e gera histogramas em formato .jpeg.
    """
    print(f"--- Iniciando processo para a amostra: {process_name} ---")

    # 1. Define os caminhos de entrada e saída
    base_dir = f"/eos/user/g/gvian/job/{process_name}_SFs"
    input_file = os.path.join(base_dir, f"{process_name}_sf_data.csv")
    output_dir_plots = os.path.join(base_dir, "histograms_jpeg")

    # Verifica se o arquivo de entrada existe
    if not os.path.exists(input_file):
        print(f"\nAVISO: Arquivo de entrada não encontrado em '{input_file}'")
        print(f"Pulando o processo '{process_name}'.\n")
        return # Pula para o próximo processo em vez de encerrar

    # Cria o diretório de saída para os plots se ele não existir
    if not os.path.exists(output_dir_plots):
        os.makedirs(output_dir_plots)
        print(f"Diretório de saída para plots criado em: '{output_dir_plots}'")

    # 2. Carrega os dados usando pandas
    print(f"Lendo dados de '{input_file}'...")
    df = pd.read_csv(input_file)
    print(f"Dados carregados com sucesso: {len(df)} léptons encontrados.")

    # 3. Separa os dados em elétrons e múons
    electrons = df[df['lep_is_ele'] == 1]
    muons = df[df['lep_is_ele'] == 0]
    print(f"Elétrons: {len(electrons)} | Múons: {len(muons)}")

    # Define os SFs a serem plotados
    electron_sfs = {'sf_trigger': 'Electron Trigger SF', 'sf_reco': 'Electron Reco SF', 'sf_id': 'Electron ID SF'}
    muon_sfs = {'sf_trigger': 'Muon Trigger SF', 'sf_id': 'Muon ID SF', 'sf_iso': 'Muon Iso SF'}
    
    # 4. Gera e salva os histogramas
    plot_range = (0.9, 1.1)
    
    print("\nGerando histogramas para Elétrons...")
    for sf_name, title in electron_sfs.items():
        if electrons.empty: continue
        
        plt.figure(figsize=(10, 7))
        electrons[sf_name].hist(bins=100, range=plot_range)
        plt.title(f"{process_name} - {title}", fontsize=16)
        plt.xlabel("Scale Factor Value", fontsize=12)
        plt.ylabel("Events", fontsize=12)
        plt.grid(True)
        
        output_filename = os.path.join(output_dir_plots, f"electron_{sf_name}.jpeg")
        plt.savefig(output_filename)
        plt.close()
        print(f" -> Salvo: {output_filename}")

    print("\nGerando histogramas para Múons...")
    for sf_name, title in muon_sfs.items():
        if muons.empty: continue

        plt.figure(figsize=(10, 7))
        muons[sf_name].hist(bins=100, range=plot_range)
        plt.title(f"{process_name} - {title}", fontsize=16)
        plt.xlabel("Scale Factor Value", fontsize=12)
        plt.ylabel("Events", fontsize=12)
        plt.grid(True)

        output_filename = os.path.join(output_dir_plots, f"muon_{sf_name}.jpeg")
        plt.savefig(output_filename)
        plt.close()
        print(f" -> Salvo: {output_filename}")
        
    print(f"\n--- Processo '{process_name}' concluído com sucesso! ---")


# ===================================================================
# PONTO DE ENTRADA DO SCRIPT - LÓGICA DE LINHA DE COMANDO MODIFICADA
# ===================================================================
if __name__ == "__main__":
    # Pega a lista de nomes de processos dos argumentos do terminal
    # sys.argv[0] é o nome do script, então pegamos do 1 em diante
    process_names = sys.argv[1:]

    # Verifica se pelo menos um nome de processo foi fornecido
    if not process_names:
        print("\nERRO: Nenhum nome de processo foi fornecido.")
        print("Uso: python plot_sfs.py <nome_do_processo_1> [nome_do_processo_2] ...")
        sys.exit(1)

    # Itera sobre cada nome de processo fornecido e gera os plots
    for name in process_names:
        make_sf_histograms(process_name=name)
        print("-" * 60) # Adiciona um separador para clareza entre os processos
