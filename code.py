import re
import os
from cigar import Cigar

idem = []
poli_ssr = []

# Itera sobre os arquivos no diretório atual
for bed in os.listdir(os.getcwd()):
    # Verifica se o nome do arquivo segue um padrão específico
    if re.match(".*\w+-\w+\.bed", bed):
        # Extrai informações do nome do arquivo usando expressões regulares
        # para obter informações sobre início, fim, nome do indivíduo e locus do SSR
        inicio_ssr = re.findall(r'.*_(\d+)', bed)
        inicio_ssr = int(str(inicio_ssr[0]))
        fim_ssr = re.findall(r'.*-(\d+).bed', bed)
        fim_ssr = int(str(fim_ssr[0]))
        nome_indiv = re.findall(r'^(\w+_\w+).', bed)
        locus_ssr = re.findall(r'(scaffold\d+(?:_size\d+)?)', bed)

        # Abre o arquivo para leitura
        f = open(bed, "r")
        linha = f.readline()
        conta_read = {}
        # Conta a ocorrência de cada read no arquivo
        while True:
            linha = f.readline()
            if linha == "":
                break
            linha = linha.replace('n', '')
            coluna = linha.split('\t')
            read = coluna[3]
            if read in conta_read.keys():
                conta_read[read] += 1
            else:
                conta_read[read] = 1
        f.close()

        # Abre novamente o arquivo para mais operações
        f = open(bed, 'r')
        linha = f.readline()
        conta_padrao = {}
        posicao_inicial_read = {}
        while True:
            linha = f.readline()
            if linha == '':
                break
            linha = linha.replace('\n', '')
            coluna = linha.split('\t')
            read = coluna[3]
            if conta_read[read] == 1:
                continue
            st_cigar = coluna[6]
            inicio_read = coluna[1]
            if read not in posicao_inicial_read:
                posicao_inicial_read[read] = inicio_read

            chave_final = posicao_inicial_read[read] + '_' + st_cigar
            if chave_final in conta_padrao.keys():
                conta_padrao[chave_final] += 1
            else:
                conta_padrao[chave_final] = 1
            start_read = chave_final.split('_')
            start_read = start_read[0]
            start_read = int(start_read)
        f.close()

        # Verifica padrões específicos no arquivo para identificar reads polimórficos
        for chave_final, valor in conta_padrao.items():
            if valor % 2 == 0 and valor >= 6:
                cigar_1 = chave_final.split('_')[1]
                c = Cigar(cigar_1)
                parser = c.items()
                contador = start_read
                fim_do_indel = 0
                for i in parser:
                    posicao_cigar = i[0]
                    tipo_cigar = i[1]
                    for cigar in tipo_cigar:
                        if cigar == 'M' and posicao_cigar == 127:
                            idem.append(bed + "_" + "127M" + "(" + str(conta_padrao[chave_final] / 2) + ")")
                        if cigar in ('M', 'H', 'S'):
                            contador = contador + posicao_cigar
                        if cigar in ('I', 'D'):
                            fim_do_indel = contador + posicao_cigar
                            if inicio_ssr <= contador <= fim_ssr or inicio_ssr <= fim_do_indel <= fim_ssr:
                                poli_ssr.append(bed + "_" + cigar_1 + "(" + str(conta_padrao[chave_final] / 2) + ")")

# Escreve os resultados em arquivos de texto
with open('indiv_polimorf.txt', 'w') as arquivo:
    for s in set(poli_ssr):
        arquivo.write("%s\n" % s)
with open('indiv_idem.txt', 'w') as arquivo:
    for s in set(idem):
        arquivo.write('%s\n' % s)
