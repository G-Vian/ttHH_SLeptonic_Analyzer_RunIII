import os
import subprocess

# ------------------------------
# Configurações fixas
exe = "./ttHHanalyzer_trigger"
folder = "filelistTrigger/TTH"   # pasta com os arquivos .txt

# Defina o nome do processo aqui
process_name = "ttH_MC_V15"

# Pasta de saída no EOS
output_folder = f"/eos/user/g/gvian/{process_name}"

xsec = "0.002892"
year = "2024"
dataType = "MC"
sampleName = process_name
# ------------------------------

# Cria a pasta de saída se não existir
os.makedirs(output_folder, exist_ok=True)

# ------------------------------------------------
# OPÇÃO 1: Ler todos os arquivos .txt da pasta
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
    
    # saída usando o nome do arquivo .txt
    output_root = os.path.join(output_folder, os.path.splitext(txt_file)[0] + ".root")
    
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
