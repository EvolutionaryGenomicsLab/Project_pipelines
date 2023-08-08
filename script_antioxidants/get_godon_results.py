import json # Para ler os resultados
import os # Para poder fazer a lista dos arquivos na pasta
import numpy as np # Para ver quais sÃ­tios tem a probabilidade posterior de estar na categoria de selecao positiva > 0.95
import itertools

print("Gene\tLRT\tSite")

for arquivo_json in os.listdir("/dados/giovanna/analises_IC/teste_Giovanna/results/GODON/output"):
        if arquivo_json.endswith(".json"):
                arquivo = open(arquivo_json, "r")
                output_json = json.load(arquivo)
                lnL_null = output_json['H0']['maxLnL']
                lnL_alte = output_json['H1']['maxLnL']
                LRT = 2 * (lnL_alte - lnL_null)
                if LRT >= 3.84:
                        print(f"{arquivo_json[:-5]}\t{LRT:.2f}\t{','.join(str((index+1,item)[0]) for (index,item) in enumerate(output_json['H1']['final']['sitePosteriorBEB']) if item >= 0.95)}") 
