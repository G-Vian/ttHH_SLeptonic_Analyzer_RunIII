import os
import subprocess

# valores fixos
exe = "./ttHHanalyzer_trigger"
folder = "filelistTrigger/TTH"
xsec = "0.002892"
year = "2024"
dataType = "MC"
sampleName = "ttH_MC_V15"

# ------------------------------------------------
# OPÇÃO 1: Ler todos os arquivos .txt da pasta (padrão)
txt_files = sorted([f for f in os.listdir(folder) if f.endswith(".txt")])

# ------------------------------------------------
# OPÇÃO 2: Ler apenas arquivos dentro de um intervalo específico
# Basta descomentar este bloco e comentar o bloco da opção 1 acima
#
# start_index = 18
# end_index = 22  # inclusive
# txt_files = [f"full_TTH_24_15_{i}.txt" for i in range(start_index, end_index + 1)]
# ------------------------------------------------

if not txt_files:
    print(f"Nenhum .txt encontrado em {folder}")
    exit(1)

for i, txt_file in enumerate(txt_files):
    input_txt = os.path.join(folder, txt_file)
    
    # saída pode ser pelo índice ou pelo nome do arquivo de entrada
    # output_root = f"ttH_V15_{i}.root"  # opção 1: usar índice
    output_root = os.path.splitext(txt_file)[0] + ".root"  # opção 2: manter nome do .txt
    
    cmd = [
        exe,
        input_txt,
        output_root,
        xsec,
        year,
        dataType,
        sampleName
    ]
    
    print("Executando:", " ".join(cmd))
    subprocess.run(cmd)

