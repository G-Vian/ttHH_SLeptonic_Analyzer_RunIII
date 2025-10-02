# -*- coding: utf-8 -*-

import matplotlib
matplotlib.use('Agg') # Essencial para rodar em servidores sem interface grafica

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import os
import sys

def make_sf_scatter_plots(process_name):
    """
    Le o arquivo .csv com os dados de SF e gera graficos de dispersao (scatter plots)
    dos SFs vs. pT e vs. eta, salvando-os como .jpeg.
    """
    print(f"--- Iniciando processo para a amostra: {process_name} ---")

    # 1. Define os caminhos de entrada e saida
    base_dir_input = f"/eos/user/g/gvian/job/{process_name}_SFs"
    input_file = os.path.join(base_dir_input, f"{process_name}_sf_data.csv")
    
    # Salva os plots em um diretorio local para facil acesso
    output_dir_plots = f"scatter_plots_{process_name}_jpeg"

    # Verifica se o arquivo de entrada existe
    if not os.path.exists(input_file):
        print(f"\nAVISO: Arquivo de entrada nao encontrado em '{input_file}'")
        print(f"Pulando o processo '{process_name}'.\n")
        return

    # Cria o diretorio de saida para os plots (localmente) se ele nao existir
    if not os.path.exists(output_dir_plots):
        os.makedirs(output_dir_plots)
        print(f"Diretorio de saida para plots criado em: '{output_dir_plots}' (no diretorio atual)")

    # 2. Carrega os dados
    print(f"Lendo dados de '{input_file}'...")
    df = pd.read_csv(input_file)
    print(f"Dados carregados com sucesso: {len(df)} leptons encontrados.")

    # 3. Separa os dados
    electrons = df[df['lep_is_ele'] == 1]
    muons = df[df['lep_is_ele'] == 0]
    print(f"Eletrons: {len(electrons)} | Muons: {len(muons)}")

    # Define os SFs a serem plotados
    electron_sfs = {'sf_trigger': 'Electron Trigger SF', 'sf_reco': 'Electron Reco SF', 'sf_id': 'Electron ID SF'}
    muon_sfs = {'sf_trigger': 'Muon Trigger SF', 'sf_id': 'Muon ID SF', 'sf_iso': 'Muon Iso SF'}
    
    # 4. Gera e salva os graficos de dispersao
    
    print("\nGerando scatter plots para Eletrons...")
    for sf_name, title in electron_sfs.items():
        if electrons.empty: continue
        
        # --- Plot SF vs. pT ---
        fig, ax = plt.subplots(figsize=(12, 8))
        ax.scatter(electrons['lep_pt'], electrons[sf_name], alpha=0.5, s=5) # s=5 define o tamanho do ponto
        ax.set_title(f"{process_name} - {title} vs. pT", fontsize=16)
        ax.set_xlabel("Electron pT [GeV]", fontsize=12)
        ax.set_ylabel("Scale Factor Value", fontsize=12)
        ax.grid(True, which="both", linestyle='--')
        ax.set_xscale('log')
        ax.get_xaxis().set_major_formatter(mticker.ScalarFormatter())
        ax.set_xticks([10, 20, 30, 50, 75, 100, 200, 500])
        
        output_filename_pt = os.path.join(output_dir_plots, f"electron_{sf_name}_vs_pt.jpeg")
        fig.savefig(output_filename_pt)
        plt.close(fig)
        print(f" -> Salvo: {output_filename_pt}")

        # --- Plot SF vs. Eta ---
        fig, ax = plt.subplots(figsize=(12, 8))
        ax.scatter(electrons['lep_eta'], electrons[sf_name], alpha=0.5, s=5)
        ax.set_title(f"{process_name} - {title} vs. Eta", fontsize=16)
        ax.set_xlabel("Electron Eta", fontsize=12)
        ax.set_ylabel("Scale Factor Value", fontsize=12)
        ax.grid(True)
        ax.set_xlim(-2.5, 2.5) # Limita a faixa do eta para melhor visualizacao
        
        output_filename_eta = os.path.join(output_dir_plots, f"electron_{sf_name}_vs_eta.jpeg")
        fig.savefig(output_filename_eta)
        plt.close(fig)
        print(f" -> Salvo: {output_filename_eta}")


    print("\nGerando scatter plots para Muons...")
    for sf_name, title in muon_sfs.items():
        if muons.empty: continue

        # --- Plot SF vs. pT ---
        fig, ax = plt.subplots(figsize=(12, 8))
        ax.scatter(muons['lep_pt'], muons[sf_name], alpha=0.5, s=5)
        ax.set_title(f"{process_name} - {title} vs. pT", fontsize=16)
        ax.set_xlabel("Muon pT [GeV]", fontsize=12)
        ax.set_ylabel("Scale Factor Value", fontsize=12)
        ax.grid(True, which="both", linestyle='--')
        ax.set_xscale('log')
        ax.get_xaxis().set_major_formatter(mticker.ScalarFormatter())
        ax.set_xticks([10, 20, 30, 50, 100, 200, 500, 1000])

        output_filename_pt = os.path.join(output_dir_plots, f"muon_{sf_name}_vs_pt.jpeg")
        fig.savefig(output_filename_pt)
        plt.close(fig)
        print(f" -> Salvo: {output_filename_pt}")

        # --- Plot SF vs. Eta ---
        fig, ax = plt.subplots(figsize=(12, 8))
        ax.scatter(muons['lep_eta'], muons[sf_name], alpha=0.5, s=5)
        ax.set_title(f"{process_name} - {title} vs. Eta", fontsize=16)
        ax.set_xlabel("Muon Eta", fontsize=12)
        ax.set_ylabel("Scale Factor Value", fontsize=12)
        ax.grid(True)
        ax.set_xlim(-2.4, 2.4)

        output_filename_eta = os.path.join(output_dir_plots, f"muon_{sf_name}_vs_eta.jpeg")
        fig.savefig(output_filename_eta)
        plt.close(fig)
        print(f" -> Salvo: {output_filename_eta}")
        
    print(f"\n--- Processo '{process_name}' concluido com sucesso! ---")


# ===============================================
# PONTO DE ENTRADA DO SCRIPT
# ===============================================
if __name__ == "__main__":
    process_names = sys.argv[1:]

    if not process_names:
        print("\nERRO: Nenhum nome de processo foi fornecido.")
        print("Uso: python3 plot_sfs.py <nome_do_processo_1> [nome_do_processo_2] ...")
        sys.exit(1)

    for name in process_names:
        make_sf_scatter_plots(process_name=name)
        print("-" * 60)
##
#     How to use : python3 plot_sfs.py <process name>
##
##    To work first time may have to use : 
##
##    python3 -m pip uninstall pandas
##    python3 -m pip install --ignore-installed --no-cache-dir pandas
