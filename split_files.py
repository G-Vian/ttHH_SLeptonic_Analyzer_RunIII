# salvar como split_files.py
# rodar com: python3 split_files.py <nome_da_pasta_de_saida>

import os
import sys

# Arquivo de entrada com as linhas
input_file = "filelistTrigger"

# Pega o nome da pasta de saída a partir do argumento da linha de comando
if len(sys.argv) < 2:
    print("Uso: python3 split_files.py <nome_da_pasta_de_saida>")
    sys.exit(1)

output_dir = sys.argv[1]

# Cria a pasta se não existir
os.makedirs(output_dir, exist_ok=True)

# Lê todas as linhas do arquivo de entrada
with open(input_file, "r") as f:
    lines = [line.strip() for line in f if line.strip()]

# Cria um txt para cada linha
for i, line in enumerate(lines):
    filename = os.path.join(output_dir, f"full_TTH_24_15_{i}.txt")
    with open(filename, "w") as f:
        f.write(line + "\n")
    print(f"Criado {filename}")
