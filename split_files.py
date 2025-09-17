# salvar como split_files.py
# rodar com: python3 split_files.py <nome_do_arquivo_txt>

import os
import sys

# Pasta onde estão os arquivos de entrada
base_dir = "filelistTrigger"

if len(sys.argv) < 2:
    print("Uso: python3 split_files.py <nome_do_arquivo_txt>")
    sys.exit(1)

# Nome do arquivo de entrada (ex.: dataset.txt)
input_filename = sys.argv[1]
input_path = os.path.join(base_dir, input_filename)

# Garante que o arquivo existe
if not os.path.isfile(input_path):
    print(f"Arquivo {input_path} não encontrado!")
    sys.exit(1)

# Nome da pasta de saída = nome do arquivo sem extensão
output_dir = os.path.join(base_dir, os.path.splitext(input_filename)[0])
os.makedirs(output_dir, exist_ok=True)

# Lê as linhas do arquivo de entrada
with open(input_path, "r") as f:
    lines = [line.strip() for line in f if line.strip()]

# Cria um txt para cada linha dentro da subpasta
for i, line in enumerate(lines):
    filename = os.path.join(output_dir, f"{os.path.splitext(input_filename)[0]}_{i}.txt")
    with open(filename, "w") as f:
        f.write(line + "\n")
    print(f"Criado {filename}")
