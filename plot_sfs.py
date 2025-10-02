# -*- coding: utf-8 -*-

import matplotlib
matplotlib.use('Agg') # Essencial para rodar em servidores sem interface grafica

import pandas as pd
import numpy as np # Importamos numpy para definir os bins facilmente
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import os
import sys

def make_sf_plots(process_name):
    """
    Le o arquivo .csv com os dados de SF, gera graficos de dispersao
    e graficos separados para a media do SF por bin.
    """
    print(f"--- Iniciando processo para a amostra: {process_name} ---")

    # 1. Define os caminhos
    base_dir_input = f"/eos/user/g/gvian/job/{process_name}_SFs"
    input_file = os.path.join(base_dir_input, f"{process_name}_sf_data.csv")
    output_dir_plots = f"sf_plots_{process_name}_jpeg"

    if not os.path.exists(input_file):
        print(f"\nAVISO: Arquivo de entrada nao encontrado em '{input_file}'")
        print(f"Pulando o processo '{process_name}'.\n")
        return

    if not os.path.exists(output_dir_plots):
        os.makedirs(output_dir_plots)
        print(f"Diretorio de saida para plots criado em: '{output_dir_plots}' (no diretorio atual)")

    # 2. Carrega e prepara os dados
    print(f"Lendo dados de '{input_file}'...")
    df = pd.read_csv(input_file)
    print(f"Dados carregados com sucesso: {len(df)} leptons encontrados.")

    electrons = df[df['lep_is_ele'] == 1]
    muons = df[df['lep_is_ele'] == 0]
    print(f"Eletrons: {len(electrons)} | Muons: {len(muons)}")

    electron_sfs = {'sf_trigger': 'Electron Trigger SF', 'sf_reco': 'Electron Reco SF', 'sf_id': 'Electron ID SF'}
    muon_sfs = {'sf_trigger': 'Muon Trigger SF', 'sf_id': 'Muon ID SF', 'sf_iso': 'Muon Iso SF'}
    
    # 3. Define os bins com escala "normal" (linear)
    pt_bins = np.linspace(0, 500, 25) # 25 bins lineares de 0 a 500 GeV
    eta_bins = np.linspace(-2.5, 2.5, 25) # 25 bins lineares de -2.5 a 2.5
    
    # 4. Gera e salva os graficos
    
    print("\nGerando plots para Eletrons...")
    for sf_name, title in electron_sfs.items():
        if electrons.empty: continue
        
        # --- PLOT 1: Scatter SF vs. pT ---
        fig, ax = plt.subplots(figsize=(12, 8))
        ax.scatter(electrons['lep_pt'], electrons[sf_name], alpha=0.2, s=5)
        ax.set_title(f"{process_name} - {title} vs. pT (Dispersao)", fontsize=16)
        ax.set_xlabel("Electron pT [GeV]", fontsize=12)
        ax.set_ylabel("Scale Factor Value", fontsize=12)
        ax.grid(True, linestyle='--')
        ax.set_xlim(pt_bins[0], pt_bins[-1])
        
        output_filename = os.path.join(output_dir_plots, f"electron_{sf_name}_vs_pt_scatter.jpeg")
        fig.savefig(output_filename)
        plt.close(fig)
        print(f" -> Salvo: {output_filename}")

        # --- PLOT 2: Media de SF vs. pT ---
        electrons['pt_bin'] = pd.cut(electrons['lep_pt'], bins=pt_bins)
        stats = electrons.groupby('pt_bin')[sf_name].agg(['mean', 'std']).dropna()
        bin_centers = [b.mid for b in stats.index]
        
        fig, ax = plt.subplots(figsize=(12, 8))
        ax.errorbar(bin_centers, stats['mean'], yerr=stats['std'], fmt='o-', color='red', capsize=5)
        ax.set_title(f"{process_name} - {title} vs. pT (Media por Bin)", fontsize=16)
        ax.set_xlabel("Electron pT [GeV]", fontsize=12)
        ax.set_ylabel("Mean Scale Factor", fontsize=12)
        ax.grid(True, linestyle='--')
        ax.set_xlim(pt_bins[0], pt_bins[-1])

        output_filename = os.path.join(output_dir_plots, f"electron_{sf_name}_vs_pt_mean.jpeg")
        fig.savefig(output_filename)
        plt.close(fig)
        print(f" -> Salvo: {output_filename}")

        # --- PLOT 3: Scatter SF vs. Eta ---
        fig, ax = plt.subplots(figsize=(12, 8))
        ax.scatter(electrons['lep_eta'], electrons[sf_name], alpha=0.2, s=5)
        ax.set_title(f"{process_name} - {title} vs. Eta (Dispersao)", fontsize=16)
        ax.set_xlabel("Electron Eta", fontsize=12)
        ax.set_ylabel("Scale Factor Value", fontsize=12)
        ax.grid(True, linestyle='--')
        ax.set_xlim(eta_bins[0], eta_bins[-1])
        
        output_filename = os.path.join(output_dir_plots, f"electron_{sf_name}_vs_eta_scatter.jpeg")
        fig.savefig(output_filename)
        plt.close(fig)
        print(f" -> Salvo: {output_filename}")

        # --- PLOT 4: Media de SF vs. Eta ---
        electrons['eta_bin'] = pd.cut(electrons['lep_eta'], bins=eta_bins)
        stats = electrons.groupby('eta_bin')[sf_name].agg(['mean', 'std']).dropna()
        bin_centers = [b.mid for b in stats.index]

        fig, ax = plt.subplots(figsize=(12, 8))
        ax.errorbar(bin_centers, stats['mean'], yerr=stats['std'], fmt='o-', color='red', capsize=5)
        ax.set_title(f"{process_name} - {title} vs. Eta (Media por Bin)", fontsize=16)
        ax.set_xlabel("Electron Eta", fontsize=12)
        ax.set_ylabel("Mean Scale Factor", fontsize=12)
        ax.grid(True, linestyle='--')
        ax.set_xlim(eta_bins[0], eta_bins[-1])

        output_filename = os.path.join(output_dir_plots, f"electron_{sf_name}_vs_eta_mean.jpeg")
        fig.savefig(output_filename)
        plt.close(fig)
        print(f" -> Salvo: {output_filename}")

    print("\nGerando plots para Muons...")
    for sf_name, title in muon_sfs.items():
        if muons.empty: continue
        
        # --- PLOT 1: Scatter SF vs. pT ---
        fig, ax = plt.subplots(figsize=(12, 8))
        ax.scatter(muons['lep_pt'], muons[sf_name], alpha=0.2, s=5)
        ax.set_title(f"{process_name} - {title} vs. pT (Dispersao)", fontsize=16)
        ax.set_xlabel("Muon pT [GeV]", fontsize=12)
        ax.set_ylabel("Scale Factor Value", fontsize=12)
        ax.grid(True, linestyle='--')
        ax.set_xlim(pt_bins[0], pt_bins[-1])
        
        output_filename = os.path.join(output_dir_plots, f"muon_{sf_name}_vs_pt_scatter.jpeg")
        fig.savefig(output_filename)
        plt.close(fig)
        print(f" -> Salvo: {output_filename}")

        # --- PLOT 2: Media de SF vs. pT ---
        muons['pt_bin'] = pd.cut(muons['lep_pt'], bins=pt_bins)
        stats = muons.groupby('pt_bin')[sf_name].agg(['mean', 'std']).dropna()
        bin_centers = [b.mid for b in stats.index]
        
        fig, ax = plt.subplots(figsize=(12, 8))
        ax.errorbar(bin_centers, stats['mean'], yerr=stats['std'], fmt='o-', color='red', capsize=5)
        ax.set_title(f"{process_name} - {title} vs. pT (Media por Bin)", fontsize=16)
        ax.set_xlabel("Muon pT [GeV]", fontsize=12)
        ax.set_ylabel("Mean Scale Factor", fontsize=12)
        ax.grid(True, linestyle='--')
        ax.set_xlim(pt_bins[0], pt_bins[-1])

        output_filename = os.path.join(output_dir_plots, f"muon_{sf_name}_vs_pt_mean.jpeg")
        fig.savefig(output_filename)
        plt.close(fig)
        print(f" -> Salvo: {output_filename}")

        # --- PLOT 3: Scatter SF vs. Eta ---
        fig, ax = plt.subplots(figsize=(12, 8))
        ax.scatter(muons['lep_eta'], muons[sf_name], alpha=0.2, s=5)
        ax.set_title(f"{process_name} - {title} vs. Eta (Dispersao)", fontsize=16)
        ax.set_xlabel("Muon Eta", fontsize=12)
        ax.set_ylabel("Scale Factor Value", fontsize=12)
        ax.grid(True, linestyle='--')
        ax.set_xlim(eta_bins[0], eta_bins[-1])
        
        output_filename = os.path.join(output_dir_plots, f"muon_{sf_name}_vs_eta_scatter.jpeg")
        fig.savefig(output_filename)
        plt.close(fig)
        print(f" -> Salvo: {output_filename}")

        # --- PLOT 4: Media de SF vs. Eta ---
        muons['eta_bin'] = pd.cut(muons['lep_eta'], bins=eta_bins)
        stats = muons.groupby('eta_bin')[sf_name].agg(['mean', 'std']).dropna()
        bin_centers = [b.mid for b in stats.index]

        fig, ax = plt.subplots(figsize=(12, 8))
        ax.errorbar(bin_centers, stats['mean'], yerr=stats['std'], fmt='o-', color='red', capsize=5)
        ax.set_title(f"{process_name} - {title} vs. Eta (Media por Bin)", fontsize=16)
        ax.set_xlabel("Muon Eta", fontsize=12)
        ax.set_ylabel("Mean Scale Factor", fontsize=12)
        ax.grid(True, linestyle='--')
        ax.set_xlim(eta_bins[0], eta_bins[-1])

        output_filename = os.path.join(output_dir_plots, f"muon_{sf_name}_vs_eta_mean.jpeg")
        fig.savefig(output_filename)
        plt.close(fig)
        print(f" -> Salvo: {output_filename}")
        
    print(f"\n--- Processo '{process_name}' concluido com sucesso! ---")


if __name__ == "__main__":
    process_names = sys.argv[1:]
    if not process_names:
        print("\nERRO: Nenhum nome de processo foi fornecido.")
        print("Uso: python3 plot_sfs.py <nome_do_processo_1> [nome_do_processo_2] ...")
        sys.exit(1)

    for name in process_names:
        make_sf_plots(process_name=name)
        print("-" * 60)
##
#     How to use : python3 plot_sfs.py <process name>
##
##    To work first time may have to use : 
##
##    python3 -m pip uninstall pandas
##    python3 -m pip install --ignore-installed --no-cache-dir pandas
