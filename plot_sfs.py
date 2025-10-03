# -*- coding: utf-8 -*-

import matplotlib
matplotlib.use('Agg')

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import os
import sys

def create_plots(df, lepton_type, sf_name, title, x_var, bins, output_dir, process_name):
    # ... (esta função auxiliar não precisa de alterações) ...
    if df.empty: return

    # --- PLOT 1: Scatter Plot ---
    fig, ax = plt.subplots(figsize=(12, 8))
    ax.scatter(df[x_var], df[sf_name], alpha=0.5, s=5, color='blue', label='Individual Events')
    ax.set_title(f"{process_name} - {title} vs. {x_var.replace('lep_', '').upper()} (Scatter)", fontsize=16)
    ax.set_xlabel(f"{lepton_type} {x_var.replace('lep_', '')}", fontsize=12)
    ax.set_ylabel("Scale Factor Value", fontsize=12)
    ax.grid(True, linestyle='--')
    ax.text(0.0, 1.05, "CMS Private Work", transform=ax.transAxes, fontsize=12, fontweight='bold', va='bottom', ha='left')
    ax.text(1.0, 1.05, r"2024 year, 108.96 fb$^{-1}$ (13.6 TeV)", transform=ax.transAxes, fontsize=12, ha='right', va='bottom')
    if x_var == 'lep_pt':
        ax.set_xlim(bins[0], bins[-1])
    else: # eta
        ax.set_xlim(bins[0], bins[-1])
    ax.legend(loc='lower right')
    output_filename = os.path.join(output_dir, f"{lepton_type.lower()}_{sf_name}_vs_{x_var.replace('lep_', '')}_scatter.jpeg")
    fig.savefig(output_filename, bbox_inches='tight')
    plt.close(fig)
    print(f" -> Salvo: {output_filename}")

    # --- PLOT 2: Mean Profile Plot ---
    df.loc[:, 'bin'] = pd.cut(df[x_var], bins=bins)
    stats = df.groupby('bin')[sf_name].agg(['mean', 'std']).dropna()
    bin_centers = [b.mid for b in stats.index]
    fig, ax = plt.subplots(figsize=(12, 8))
    ax.errorbar(bin_centers, stats['mean'], yerr=stats['std'], fmt='o', color='red', capsize=5, label='Mean +/- Std. Dev.')
    ax.set_title(f"{process_name} - {title} vs. {x_var.replace('lep_', '').upper()} (Mean)", fontsize=16)
    ax.set_xlabel(f"{lepton_type} {x_var.replace('lep_', '')}", fontsize=12)
    ax.set_ylabel("Mean Scale Factor", fontsize=12)
    ax.grid(True, linestyle='--')
    ax.legend(loc='lower right')
    ax.text(0.0, 1.05, "CMS Private Work", transform=ax.transAxes, fontsize=12, fontweight='bold', va='bottom', ha='left')
    ax.text(1.0, 1.05, r"2024 year, 108.96 fb$^{-1}$ (13.6 TeV)", transform=ax.transAxes, fontsize=12, ha='right', va='bottom')
    if x_var == 'lep_pt':
        ax.set_xlim(bins[0], bins[-1])
    else: # eta
        ax.set_xlim(bins[0], bins[-1])
    output_filename = os.path.join(output_dir, f"{lepton_type.lower()}_{sf_name}_vs_{x_var.replace('lep_', '')}_mean.jpeg")
    fig.savefig(output_filename, bbox_inches='tight')
    plt.close(fig)
    print(f" -> Salvo: {output_filename}")


def make_sf_plots(process_name):
    print(f"--- Iniciando processo para a amostra: {process_name} ---")

    base_dir_input = f"/eos/user/g/gvian/job/{process_name}_SFs"
    input_file = os.path.join(base_dir_input, f"{process_name}_sf_data.csv")
    output_dir_plots = f"sf_plots_{process_name}_jpeg"

    if not os.path.exists(input_file):
        print(f"\nAVISO: Arquivo de entrada nao encontrado em '{input_file}'")
        return

    if not os.path.exists(output_dir_plots):
        os.makedirs(output_dir_plots)
        print(f"Diretorio de saida para plots criado em: '{output_dir_plots}' (no diretorio atual)")

    print(f"Lendo dados de '{input_file}'...")
    try:
        df = pd.read_csv(input_file, error_bad_lines=False, warn_bad_lines=True)
    except Exception as e:
        print(f"Ocorreu um erro ao ler o arquivo CSV: {e}")
        return

    # ===================================================================
    # --- NOVO PASSO: Limpeza e Conversão de Dados ---
    # ===================================================================
    sf_columns = ['sf_trigger', 'sf_reco', 'sf_id', 'sf_iso']
    for col in sf_columns:
        # Força a coluna a ser numérica. Qualquer valor que não puder ser convertido
        # se tornará 'NaN' (Not a Number) e será ignorado nos cálculos.
        df[col] = pd.to_numeric(df[col], errors='coerce')

    # Remove qualquer linha onde as colunas essenciais tenham valores inválidos (NaN)
    df.dropna(subset=sf_columns, inplace=True)
    # ===================================================================

    print(f"Dados carregados e limpos com sucesso: {len(df)} leptons validos encontrados.")

    electrons = df[df['lep_is_ele'] == 1].copy()
    muons = df[df['lep_is_ele'] == 0].copy()
    print(f"Eletrons: {len(electrons)} | Muons: {len(muons)}")

    electron_sfs = {'sf_trigger': 'Electron Trigger SF', 'sf_reco': 'Electron Reco SF', 'sf_id': 'Electron ID SF'}
    muon_sfs = {'sf_trigger': 'Muon Trigger SF', 'sf_id': 'Muon ID SF', 'sf_iso': 'Muon Iso SF'}
    
    pt_bins = np.linspace(0, 500, 25)
    eta_bins = np.linspace(-2.5, 2.5, 25)
    
    print("\nGerando plots para Eletrons...")
    for sf_name, title in electron_sfs.items():
        create_plots(electrons, "Electron", sf_name, title, 'lep_pt', pt_bins, output_dir_plots, process_name)
        create_plots(electrons, "Electron", sf_name, title, 'lep_eta', eta_bins, output_dir_plots, process_name)

    print("\nGerando plots para Muons...")
    for sf_name, title in muon_sfs.items():
        create_plots(muons, "Muon", sf_name, title, 'lep_pt', pt_bins, output_dir_plots, process_name)
        create_plots(muons, "Muon", sf_name, title, 'lep_eta', eta_bins, output_dir_plots, process_name)
        
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
