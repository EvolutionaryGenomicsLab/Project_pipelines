#!/bin/bash

lnl=$(for file in $(cat ~giovanna/analises_IC/teste_Giovanna/data/lista_genes.txt); do echo $'\n'lnL ${file}: && grep -H "lnL" ${file}/${file}_* | sed 's/:/ /; s/  / /g; s/)://' | awk ' -F " " {print $1" "$4" "$3}' | sort -k 2 -r; done)
omega=$(for file in $(cat ~giovanna/analises_IC/teste_Giovanna/data/lista_genes.txt); do echo $'\n'omegas ${file}: && grep -H "(dN/dS)" ${file}/${file}_* | sed 's/w (dN\/dS) for branches:  / / ; s/omega (dN\/dS) =  / /' | awk '{print length, $0}' | sort ; done)
sites=$(for file in $(cat ~giovanna/analises_IC/teste_Giovanna/data/lista_genes.txt); do grep -A 50 "(BEB)" ${file}/${file}_modelA_* | awk '/([0-9]\s[A-Z]\s[0-1]\.[0-9]|[0-9]\s(-)\s[0-1]\.[0-9])/ {print $0}';done)

echo "${lnl}\n" > results_lnl_PAML.txt
echo "${omega}\n" > results_omega_PAML.txt
echo "${sites}"\n > results_BSM.txt
