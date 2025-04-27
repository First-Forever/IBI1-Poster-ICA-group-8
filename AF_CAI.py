import pandas as pd
import re

# check if the input is right
# 只检查了是否有正确的碱基对，没有检查长度
def input_check(sequence):
    check = re.findall('[A|U|C|G]', sequence)
    if len(check) != len(sequence):
        return False
    else:
        return True 

# Choose the species you want to check
# 还没定好哪些物种
species = input('Choose the species you want to check:')

# change tsv to csv
tsv = open(str(species) + "_RSCU.tsv", 'r')
tsv2 = open(str(species) + '_RSCU_2.tsv', 'w')
counter = 0
for line in tsv:
    counter += 1
    if counter >= 8: # 这个标题去除条件可能还要再改改
        if counter >= 9: # 这个判定条件也改一下
            #print(type(line))
            line = re.sub('T', 'U', line) # sub T to U
            #print(line)
        tsv2.write(line)
tsv.close()
tsv2.close()

RSCU_table = pd.read_table(str(species) + '_RSCU_2.tsv', sep = '\t')
RSCU_table.to_csv(str(species) + '_RSCU.csv', index = False)
RSCU = pd.read_csv(str(species) + '_RSCU.csv')
# print(human_RSCU.head())
# print(human_RSCU[human_RSCU['CODON'] == 'GCT']['RSCU'])

# Create a dictionary that stores RSCU_max for each AA
RSCU_max_dic = {
    }
for i in RSCU['Amino acid']:
    AA_RSCU = RSCU[RSCU['Amino acid'] == i]['RSCU'].max()
    # print(type(AA_RSCU))
    if i not in RSCU_max_dic:
        RSCU_max_dic[i] = AA_RSCU
    else:
        continue
# print(RSCU_max_dic['Thr])

CAI = 1.000
## 这里先假定是一段标准的mRNA（AUG开头，终止密码子结尾），认为sequence是一段字符串
## 我还假定了有一个统计密码子总数的量L
L = 3.000 # The total number of codons in the sequence
frequency = {}
def CAI_calculator(sequence, CAI, L):
    for i in range(0, len(sequence),3):
        codon = sequence[i:i+3]
        # search RSCU for every codon
        RSCU_k = RSCU[RSCU['CODON'] == codon]['RSCU'].tolist()
        # search AA of that codon
        AA = RSCU[RSCU['CODON'] == codon]['Amino acid'].tolist()
        # search RSCU_max for that AA
        RSCU_max = RSCU_max_dic[AA[0]]
        #print(AA, RSCU_k, RSCU_max)
        # calculate CAI
        CAI *= RSCU_k[0] / RSCU_max
        #print(CAI)
    CAI = CAI ** (1/L)
    round(CAI, 3)
    return CAI
# 为了测试假定了sequence
sequence = 'AUGAGUUUG'
if input_check(sequence):
    print(CAI_calculator(sequence, CAI, L))
else:
    print('Wrong input')
