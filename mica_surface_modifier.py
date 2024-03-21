import numpy as np 
import  random
import argparse

parser = argparse.ArgumentParser()

parser.add_argument("-f", "--file", help="Mica_input_pdb_file")

parser.add_argument("-o", "--output", help="Mica_edited_output_file")

args = parser.parse_args()


original = args.file

modified = args.output 

f = open(original, "r") 



#Storing MG and K positions to replace

K = { }

MG = { }


for line in f:
    if line[13:16].strip() == "K1" and int(line[23:27].strip()) < 460:
        K[line[6:11]] = line[31:55]
    if line[13:16].strip() == "MG":
       MG[line[6:11]] = line[31:55]


#Take 115 potassium randomly from the dictionary K

Pot = [ ]
for k, v in K.items():
    Pot.append(k)

Pot1 = random.sample(Pot, 115)



#Coordinates of 115 Magnesium 

MG1 = [ ]

for k1, v1 in MG.items():
    MG1.append(k1)

MG2 = random.sample(MG1, 115)


print (MG2)

#Now store the coordinates of the randomly selected atoms from random.sample

ran_coor_MG = [MG[i] for i in MG2]
ran_coor_K = [K[j] for j in Pot1]



#Perform swapping 

K_rep_w_MG = { }

MG_rep_w_K =  { }


for i in range(0, len(MG2)):
    K_rep_w_MG[Pot1[i]] = ran_coor_MG[i]
    MG_rep_w_K[MG2[i]] = ran_coor_K[i]


#We are doing repeating operation we need a function for calculating the differences between MG and its associated water atoms


def cal_diff(coor_mg, coor_x):
    diff = [(float(coor_mg[0]) -float(coor_x[0])), ( float(coor_mg[1]) - float(coor_x[1])), (float(coor_mg[2]) - float(coor_x[2]))]
    return diff
       



#This is a function taking Mg and diff and spit out the new coordinates
def new_c(cor1,cor2):
    add_1 = [(float(cor1[0]) + float(cor2[0])), (float(cor1[1]) + float(cor2[1])), (float(cor1[2]) + float(cor2[2]))]
    return add_1 

f5 = open(original, "r")

f6 = f5.readlines()
def store_atoms(x):
    wat_atm = { }
    for i in range(x+1, x+19):
        for line in f6:
            if int(line[6:12].strip()) == x:
               coor_mg = [line[31:39].strip(), line[39:46].strip(), line[46:55].strip()]
            if int(line[6:12].strip()) == i:
               wat_atm[str(i)] = cal_diff(coor_mg, [line[31:39].strip(), line[39:46].strip(), line[46:55].strip()])
    return wat_atm 

#What this function will return a dictonary with atom number associated with MG residues and having differences
    
#Store the Magnesium atoms, calculate the differences 


#Function to raise the Magnesium above the surface and use these new Mg parameters to match with 

#Take the new coordiantes of Magnesium and add 1.50 along the z axis and we need write them to the file REPLACED.pdb

#Take the Mg atom in the file

#raise

def raise1(coor):
    diff_1 = [(float(coor[0])), (float(coor[1])), (float(coor[2]) - 1.50)]
    return diff_1
    

Mg_n = { }
for mg, cor in MG_rep_w_K.items():
    Mg_n[mg] =  str('{:7.3f}'.format(raise1(cor.split())[0]))  + str('{:8.3f}'.format(raise1(cor.split())[1]))  +   str('{:8.3f}'.format(raise1(cor.split())[2])+ " " )

new_coord = { }
new1 = { }
#Generate coordinates for water_atoms associate with magnesium

for mg1, cor1 in Mg_n.items():
    t = store_atoms(int(mg1))
    for t1, cor2 in t.items():
        new1[t1] = new_c(cor1.split(), cor2) # Calling function to store the modified coordinates of water atoms 

    

new_atoms = { }
for k7, v7 in new1.items():
    new_atoms[k7] =  str('{:7.3f}'.format(v7[0]))  + str('{:8.3f}'.format(v7[1]))  +   str('{:8.3f}'.format(v7[2]))+ " " 



#Generate the new coordinates for sampled magnesium along with water associated with them


K_rep_w_MG.update(Mg_n)
f.close()
K_rep_w_MG.update(new_atoms)
f1 = open(original, "r")
with open(modified, "w") as fp:
     for line1 in f1:
         if line1[6:11] in K_rep_w_MG:
             fp.writelines(line1[0:31] + K_rep_w_MG[line1[6:11]] + line1[55:78] + '\n')
         else:
             fp.writelines(line1)
