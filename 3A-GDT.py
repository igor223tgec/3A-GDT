#!/usr/bin/env python3

from multiprocessing import cpu_count
from shutil import copyfile
import re
import sys
import os
import argparse
from argparse import RawTextHelpFormatter
import subprocess
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

version = "1.20"

num_threads = cpu_count()

help = """
3A-DGT - Triple-A DGT version """+version+""" - 08 nov 2023
The All-Against-All Distance Graphic Tool
(c) 2023. Igor Custódio dos Santos & Arthur Gruber

Usage: 3A-DGT -i <input_file> -m <method> -o <output_directory>

Mandatory parameters:
-i input_file <file name>       Input file (FASTA sequence file)
-m method <pa|msa|both>         Alignment method:
                                  pa - pairwise
                                  msa - multiple sequence
                                  both - both methods
-o output <directory name>      Output directory to store files

Additional parameters:
-h --help                       This help page
-c color                        Color palette from matplotlib package for graphs (default=RdBu)
-l lines <yes|no>               Use cell dividing lines on heatmap (default=yes)
-n names* <id|taxon|all>        Names for labels:
                                  id - YP_003289293.1 (default)
                                  name - Drosophila totivirus 1
                                  all - YP_003289293.1 - Drosophila totivirus 1
-mafft "mafft parameters"	Only "--maxiterate" and "--thread" parameters for MAFFT run.
-iqtree "iqtree parameters"	Only "-m", "-bb" and "-nt" parameters for IQTREE run.
-version			Show the program's version.

*Only valid for NCBI naming format:
YP_003289293.1 RNA-dependent RNA polymerase [Drosophila totivirus 1]
For other naming formats, all terms are used
"""

parser = argparse.ArgumentParser(add_help=False, formatter_class=RawTextHelpFormatter)
parser.add_argument('-i')
parser.add_argument('-m')
parser.add_argument('-o')
parser.add_argument('-c', default='RdBu')
parser.add_argument('-l', default='yes')
parser.add_argument('-n', default="id")
parser.add_argument('-h', '--help', action='store_true')
parser.add_argument('-v', '--version', action='store_true')
parser.add_argument('-mafft', default='--maxiterate 1000 --thread '+str(round(num_threads/2)))
parser.add_argument('-iqtree', default='-m Q.pfam+F+R6 -bb 1000 -nt '+str(round(num_threads/2)), type=str)
args = parser.parse_args()

def mandatory_param_check(args):
  not_specified = []
  if args.i == None:
    not_specified.append("-i [input_file]")
  if args.m == None:
    not_specified.append("-m [method]")
  if args.o == None:
    not_specified.append("-o [output]")

  if len(not_specified) != 0:
    if len(not_specified) == 1:
      not_specified = not_specified[0]
    else:
      not_specified = ", ".join(not_specified)
    print("""Error: The command is missing mandatory parameters: """+not_specified+""".
Use the command template below adding any optional parameters you want.

3A-GDT -i [input_file] -m [method] -o [output]""")
    sys.exit()


def check_fasta(fasta):
  if os.path.isfile:
    with open(fasta, 'r') as file:
      lines = file.readlines()
    if not lines:
      print("Error: Input file is empty.")
      sys.exit()
    if not lines[0].startswith(">"):
      print("Error: Input file is not in FASTA format.")
      sys.exit()
  else:
    print("Error: Input file not found in the working directory.")
    sys.exit()
  type_header = "ncbi"
  for line in lines:
    if line.startswith(">") and not re.match(r'>.* .* \[.*\]', line):
      type_header= "not_ncbi"
      break

  is_protein = None  # Inicializamos como None para determinar o tipo na primeira sequência.
  for line in lines:
    line = line.strip()
    if line.startswith(">"):
      # Se a linha começa com ">", é um cabeçalho e não contém sequências.
      continue
    elif any(base in "RDEILKMrdeilkm" for base in line.upper()):
      is_protein = True
      break
    else:
      # Se não encontramos caracteres não nucleotídicos até agora, assumimos
      # que a sequência é de nucleotídeos.
      is_protein = False

  if is_protein:
    return type_header, "protein"
  elif is_protein is not None:
    return type_header, "nucleotide"


def check_param(args):
  if args.l not in ['yes', 'no']:
    print("Error: Line parameter (-l) not recognized.")
    sys.exit()

  color = args.c
  if not sns.color_palette(color):
    print("Error: Color parameter (-c) not recognized.")
    sys.exit()

  if args.n not in ['id', 'name', 'all']:
    print("Error: Label parameter (-n) not recognized.")
    sys.exit()

def isInt(value):
  try:
    int(value)
    return True
  except:
    return False

def check_mafft(param_mafft):
  num_threads = cpu_count()
  mafft_list = param_mafft.split(' ')
  if len(mafft_list) % 2 != 0:
    print("ERROR: MAFFT parameters not accepted. It has an odd number of parameters.")
    sys.exit()
  for i in range(0, len(mafft_list), 2):
    if mafft_list[i] == "--maxiterate":
      if not isInt(mafft_list[i+1]):
        print("ERROR: MAFFT parameter --maxiterate is not integer")
        sys.exit()
    elif mafft_list[i] == "--thread":
      if not isInt(mafft_list[i+1]):
        print("ERROR: MAFFT parameter --thread is not integer")
        sys.exit()
      else:
        if int(mafft_list[i+1]) > num_threads:
          print("Number of processors ("+mafft_list[i+1]+") exceeds the total CPUs available on the server. Using "+str(round(num_threads))+" CPUs.")
          mafft_list[i+1] = str(round(num_threads))
    else:
      print("ERROR: MAFFT parameter '"+mafft_list[i]+"' not accepted by 3A-DGT. Please correct your syntax.")
      sys.exit()
  if len(mafft_list) == 2:
    if "--maxiterate" not in mafft_list:
      mafft_list.append("--maxiterate")
      mafft_list.append("1000")
    if "--thread" not in mafft_list:
      mafft_list.append("--thread")
      mafft_list.append(str(round(num_threads/2)))

  mafft_cmd = " ".join(mafft_list)

  return mafft_cmd


def check_iqtree(param_iqtree):
  num_threads = cpu_count()
  iqtree_list = param_iqtree.split(' ')
  if len(iqtree_list) % 2 != 0:
    print("ERROR: IQ-TREE parameters not accepted. It has an odd number of parameters.")
    sys.exit()
  for i in range(0, len(iqtree_list), 2):
    if iqtree_list[i] == "-bb":
      if not isInt(iqtree_list[i+1]):
        print("ERROR: IQ-TREE parameter -bb is not integer.")
        sys.exit()
      if isInt(iqtree_list[i+1]) and int(iqtree_list[i+1]) < 1000:
        print("ERROR: IQ-TREE parameter -bb must be greater than 1000.")
        sys.exit()
    elif iqtree_list[i] == "-nt":
      if not isInt(iqtree_list[i+1]):
        print("ERROR: IQ-TREE parameter -nt is not integer")
        sys.exit()
      else:
        if int(iqtree_list[i+1]) > num_threads:
          print("Number of processors ("+iqtree_list[i+1]+") exceeds the total CPUs available on the server. Using "+str(round(num_threads))+" CPUs.")
          iqtree_list[i+1] = str(round(num_threads))
    elif iqtree_list[i] == "-m":
      continue
    else:
      print("ERROR: IQ-TREE parameter '"+iqtree_list[i]+"' not accepted by 3A-DGT. Please correct your syntax.")
      sys.exit()

  if len(iqtree_list) == 2:
    if "-bb" not in iqtree_list:
      iqtree_list.append("-bb")
      iqtree_list.append("1000")
    if "-nt" not in iqtree_list:
      iqtree_list.append("-nt")
      iqtree_list.append(str(round(num_threads/2)))
    if "-m" not in iqtree_list:
      iqtree_list.append("-m")
      iqtree_list.append("Q.pfam+F+R6")

  iqtree_cmd = " ".join(iqtree_list)
  return iqtree_cmd


def where_to_save_graphics(output):
  if os.path.isdir(output+"_dir"):
    print("Output directory alrealdy exists: "+output+"_dir")
    i = 2
    if os.path.isdir(output+"_dir/graphics_dir"):
      while True:
        if os.path.isdir(output+"_dir/graphics_dir_"+str(i)):
          i +=1
          continue
        else:
          save_graphics_dir = output+"_dir/graphics_dir_"+str(i)
          break
    else:
      save_graphics_dir = output+"_dir/graphics_dir"
  else:
    os.mkdir(output+"_dir")
    save_graphics_dir = output+"_dir/graphics_dir"

  return save_graphics_dir

def entry_file(args, type_header):
  if type_header == 'not_ncbi':
    with open(args.i, 'r') as file:
      with open(args.o+"_dir/renamed_"+args.i, 'w') as file2:
        for line in file:
          if line.startswith(">"):
            line = line.replace('(', '')
            line = line.replace(')', '')
            line = line.replace(']', '')
            line = line.replace('[', '')
            line = line.replace(' ', "_")
          file2.write(line)
    enter = args.o+"_dir/"+args.o+"_concatened_headers.fasta"

  elif type_header == "ncbi":
    enter = args.i

  return enter

def correct_label(data, args):
  if args.n == 'id':
    return data
  elif args.n == "name" or args.n == "all":
    new_data = []
    with open(args.i, "r") as read_fasta:
      for line in read_fasta:
        if line.startswith(">"):
          new_data.append([line[1:].split(' ')[0], line.split("[")[1].rsplit("]")[0]])
    for o in data:
      for t in o:
        for n in new_data:
          if t == n[0]:
            if args.n == "name":
              data[data.index(o)][o.index(t)] = n[1]
            elif args.n == "all":
              data[data.index(o)][o.index(t)] = n[0]+" - "+n[1]
    return data

def run_needle(fasta_file, output_dir, type_header, type_sequence):
  save_needle_dir = output_dir+"_dir/needle_dir"
  if os.path.isdir(save_needle_dir) == False:
    os.mkdir(save_needle_dir)
  sequencias = []  # Inicializa uma lista para armazenar as sequências
  sequencia_atual = []  # Inicializa uma lista temporária para armazenar a sequência atual

  with open(fasta_file, 'r') as arquivo:
    for linha in arquivo:
      linha = linha.strip()
      if linha.startswith(">"):
        # Se a linha é um cabeçalho, isso indica o início de uma nova sequência
        if sequencia_atual:
          sequencias.append(sequencia_atual)  # Adiciona a sequência à lista
        sequencia_atual = [linha]  # Inicia uma nova sequência com o cabeçalho
      else:
        # Se a linha não é um cabeçalho, ela faz parte da sequência atual
        sequencia_atual.append(linha)

    if sequencia_atual:
      sequencias.append(sequencia_atual)  # Adiciona a última sequência

  needle_file = save_needle_dir+"/"+output_dir+".needle"
  with open(needle_file, 'w') as open_needle:
    open_needle.write("")

  for seq in sequencias:
    in_seq = seq[0][1:].split(' ')[0]
    in_seq_file = seq[0][1:].split(' ')[0]+'.fasta'
    out_seq_file = seq[0][1:].split(' ')[0]+'.needle'
    with open(in_seq_file, 'w') as arquivo:
      arquivo.write('\n'.join(seq))
    with open('comparated_sequences.fasta', 'w') as arquivo:
      for l in sequencias[sequencias.index(seq):]:
        arquivo.write('\n'.join(l)+'\n')
    try:
      if type_sequence == 'protein':
        subprocess.run(["needle", "-sprotein1", in_seq_file, "-sprotein2",
                       "comparated_sequences.fasta", "-gapopen", "10.0",
                       "-gapextend", "0.5", "-datafile", "EBLOSUM62",
                       "-outfile", out_seq_file])
      if type_sequence == 'nucleotide':
        subprocess.run(["needle", "-snucleotide1", in_seq_file, "-snucleotide2",
                       "comparated_sequences.fasta", "-gapopen", "10.0",
                       "-gapextend", "0.5", "-datafile", "EDNAFULL",
                       "-outfile", out_seq_file])
    except:
      print("Error at running Needle.")
      sys.exit
    else:
      print("Running Needle with "+in_seq+"...")
    os.remove(in_seq_file)
    os.remove('comparated_sequences.fasta')
    with open(out_seq_file, 'r') as open_moment_needle:
      with open(needle_file, 'a') as open_needle:
        for linha in open_moment_needle:
          open_needle.write(linha)
    os.remove(out_seq_file)

def data_from_needle(file, type_sequence):
  data = []
  data_ident = []
  data_sim = []
  with open(file, 'r') as open_needle:
    pair_sim = []
    pair_ident = []
    for needle_line in open_needle:
      if needle_line.startswith("# 1:"):
        pair_sim.append(needle_line[5:-1])
        pair_ident.append(needle_line[5:-1])
      elif needle_line.startswith("# 2:"):
        pair_sim.append(needle_line[5:-1])
        pair_ident.append(needle_line[5:-1])
      elif needle_line.startswith("# Identity:"):
        pair_ident.append(round(float(needle_line.split("(")[1].split("%")[0])/100, 4))
      if type_sequence == "protein":
        if needle_line.startswith("# Similarity:"):
          pair_sim.append(round(float(needle_line.split("(")[1].split("%")[0])/100, 4))
      else:
        continue
      if len(pair_sim) == 3:
        data_sim.append(pair_sim)
        pair_sim = []
      if len(pair_ident) == 3:
        data_ident.append(pair_ident)
        pair_ident = []
    data.append(data_ident)
    if type_sequence == "protein":
      data.append(data_sim)
  return data

def save_matrix(perc_list, save_path):
  for triples in perc_list:
    if perc_list.index(triples) == 0:
      matrix_name = save_path+"_ident_matrix.csv"
    if perc_list.index(triples) == 1:
      matrix_name = save_path+"_simil_matrix.csv"
    order = []
    for l in triples:
      if l[0] not in order:
        order.append(l[0])
    matriz = []
    for k in range(len(order)):
      linha = [None] * len(order)
      matriz.append(linha)

    for element in triples:
      i = order.index(element[0])
      j = order.index(element[1])
      if i < j:
        h = i
        i = j
        j = h
      matriz[j][i] = matriz[i][j] = element[2]

    table = {}
    table[''] = order
    for l in range(0, len(matriz)):
      table[order[l]] = matriz[l]
    final_table = pd.DataFrame(table)
    final_table.set_index('', inplace=True)

    final_csv = final_table.to_csv(index = True)

    with open(matrix_name, "w") as file:
      file.write(final_csv)

def run_mafft(fasta_file, output, cmd_mafft):
  save_align_dir = output+"_dir/mafft_dir"
  if os.path.isdir(save_align_dir) == False:
    os.mkdir(save_align_dir)
  real_path_fasta = os.path.realpath(fasta_file)
  real_path_output = os.path.realpath(save_align_dir+"/"+output+".align")
  error_out = os.path.realpath(save_align_dir+"/error_"+output+"_mafft")
  cmd = 'mafft '+cmd_mafft+' '+real_path_fasta+' > '+real_path_output+' 2> '+error_out
  subprocess.call(cmd, shell=True)

def run_iqtree(align_file, output, cmd_iqtree):
  save_iqtree_dir = output+"_dir/iqtree_dir"
  if os.path.isdir(save_iqtree_dir) == False:
    os.mkdir(save_iqtree_dir)
  copyfile(align_file, save_iqtree_dir+"/"+output)
  real_path_align = os.path.realpath(save_iqtree_dir+"/"+output)
  error_out = os.path.realpath(save_iqtree_dir+"/error1_"+output+"_iqtree")
  error_out2 = os.path.realpath(save_iqtree_dir+"/error2_"+output+"_iqtree")
  cmd = "iqtree2 -s "+real_path_align+" "+cmd_iqtree+" 1>"+error_out+" 2>"+error_out2
  try:
    subprocess.call(cmd, shell=True)
  except Exception as err:
    print("IQTREE reported an error:\n\n"+err)
  else:
    print("Running IQ-TREE...")
  os.remove(save_iqtree_dir+"/"+output)

def data_from_mldist(mldist_file):
  data = []
  with open(mldist_file, 'r') as arquivo:
    lines = arquivo.readlines()
    for line in lines[1:]:
      for distance in line.split()[1:]:
        data.append([line.split()[0], lines[line.split()[1:].index(distance)+1].split()[0], round(float(distance), 4)])
  return data

def create_clustermap(all_data, args):
  for data in all_data:
    print("Generating file: "+data)

    if args.n == "name" or args.n == "all":
      new_data = []
      with open(args.i, "r") as read_fasta:
        for line in read_fasta:
          if line.startswith(">"):
            new_data.append([line[0][1:].split(' ')[0], line.split("[")[1].rsplit("]")[0]])
      for o in all_data[data]:
        for t in o:
          for n in new_data:
            if t == n[0]:
              if args.n == "name":
                all_data[data][all_data[data].index[o]][o.index(t)] = n[1]
              elif args.n == "all":
                all_data[data][all_data[data].index[o]][o.index(t)] = n[0]+" - "+n[1]

    order = []
    for l in all_data[data]:
      if l[0] not in order:
        order.append(l[0])

    # Gera nome do arquivo final.
    jpg_file_name = data + '.jpg'
    svg_file_name = data + '.svg'

    # Cria matriz vazia.
    matriz = []
    for k in range(len(order)):
        linha = [None] * len(order)
        matriz.append(linha)

    # Posiciona os valores de similaridade na matriz de acordo com seus respectivos códigos, guiados pela ordem dos códigos do arquivo de referência.
    for element in all_data[data]:
        i = order.index(element[0])
        j = order.index(element[1])
        if i < j:
            h = i
            i = j
            j = h
        matriz[j][i] = matriz[i][j] = element[2]

    table = {}
    table[''] = order
    for l in range(0, len(matriz)):
      table[order[l]] = matriz[l]

    final_table = pd.DataFrame(table)
    final_table.set_index('', inplace=True)

    color = args.c
    if sns.color_palette(color):
      color_msa = color
      if color.endswith("_r"):
        color_pa = color[:-2]
      else:
        color_pa = color+"_r"
    else:
      print("Error: Invalid color (-c).")
      sys.exit()

    if args.l == "yes":
      linew = 0.15
    elif args.l == "no":
      linew = 0

    if jpg_file_name.endswith("_mldist.jpg"):
      clustermap = sns.clustermap(final_table, method='complete', cmap=color_msa,
                   xticklabels=True, yticklabels = True, linewidths=linew,
                   vmin=0, vmax=9, cbar_pos=(.02, .32, .03, .2),
                   figsize=(25, 25))

    elif jpg_file_name.endswith("_simil.jpg") or jpg_file_name.endswith("_ident.jpg"):
      clustermap = sns.clustermap(final_table, cmap=color_pa,
                   xticklabels=True, yticklabels = True,
                   cbar_pos=(.02, .32, .03, .2), linewidths=linew,
                   vmin=0, vmax=1,
                   figsize=(25, 25))

    clustermap.ax_row_dendrogram.set_visible(False)

    plt.savefig(jpg_file_name)
    plt.savefig(svg_file_name, format='svg')

if __name__ == '__main__':
  if not len(sys.argv)>1:
    print(help)
  elif args.help == True:
    print(help)
  elif args.version == True:
    print("""
3A-DGT - Triple-A DGT version """+version+""" - 08 nov 2023
The All-Against-All Distance Graphic Tool
(c) 2023. Igor Custódio dos Santos & Arthur Gruber
""")
  else:
    mandatory_param_check(args)
    check_param(args)
    cmd_iqtree = check_iqtree(args.iqtree)
    cmd_mafft = check_mafft(args.mafft)
    type_header, type_sequence = check_fasta(args.i)
    if type_header == 'not_ncbi':
      args.n = 'id'
    if args.m not in ['pa', 'msa', 'both']:
      print("Error: Method (-m) not recognized.")
      sys.exit()
    save_graphics_dir = where_to_save_graphics(args.o)
    save_align_dir = args.o+"_dir/mafft_dir"
    save_needle_dir = args.o+"_dir/needle_dir"
    save_iqtree_dir = args.o+"_dir/iqtree_dir"

    args.i = entry_file(args, type_header)

    all_data = {}
    if args.m == "pa" or args.m == "both":
      if os.path.isfile(save_needle_dir+"/"+args.o+".needle"):
        print("Needle file already exists. Skipping Needle step.")
      else:
        run_needle(args.i, args.o, type_header, type_sequence)
      data = data_from_needle(save_needle_dir+"/"+args.o+".needle", type_sequence)
      save_matrix(data, save_needle_dir+"/"+args.o)
      if len(data) == 1:
        data[0] = correct_label(data[0], args)
        all_data[save_graphics_dir+"/"+args.o+"_ident"] = data[0]
      else:
        data[0] = correct_label(data[0], args)
        all_data[save_graphics_dir+"/"+args.o+"_ident"] = data[0]
        data[1] = correct_label(data[1], args)
        all_data[save_graphics_dir+"/"+args.o+"_simil"] = data[1]

    if args.m == "msa" or args.m == "both":
      if os.path.isfile(save_align_dir+"/"+args.o+".align"):
        print("Mafft aligned file already exists. Skipping Mafft step.")
      else:
        run_mafft(args.i, args.o, cmd_mafft)
      if os.path.isfile(save_iqtree_dir+"/"+args.o+".mldist"):
        print("IQ-TREE maximum-likelihood distance matrix file already exists. Skipping IQ-TREE step.")
      else:
        run_iqtree(save_align_dir+"/"+args.o+".align", args.o, cmd_iqtree)

      data = data_from_mldist(save_iqtree_dir+"/"+args.o+".mldist")
      data = correct_label(data, args)
      all_data[save_graphics_dir+"/"+args.o+"_mldist"] = data

    os.mkdir(save_graphics_dir)

    create_clustermap(all_data, args)


