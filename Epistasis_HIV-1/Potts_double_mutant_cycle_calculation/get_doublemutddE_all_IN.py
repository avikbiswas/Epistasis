#!/usr/bin/env python3
from scipy import *
import numpy as np
import os,sys
from io import StringIO

def alphabet_compare1(alphabet1):
    if(alphabet1=='A'):
        num=0
    elif(alphabet1=='B'):
        num=4
    elif(alphabet1=='C'):
        num=8
    elif(alphabet1=='D'):
        num=12
    return num
    
    #return "ABCDEFGH".index(alphabet1)*8

def alphabet_compare2(alphabet2):
    if(alphabet2=='A'):
        num=0
    elif(alphabet2=='B'):
        num=1
    elif(alphabet2=='C'):
        num=2
    elif(alphabet2=='D'):
        num=3
    return num

def other_alphabet(alphabet):
    return list("ABCD".replace(alphabet, ''))


J=np.load('../data/models/J_IN.npy')
Q_squared=J.shape[1]
Q=int(np.sqrt(Q_squared)+0.2)
L=int((1+np.sqrt(1+8*J.shape[0]))/2+0.2)
print(L,Q)

map_index=np.arange(0,int(J.shape[0]))
J_mapped=np.zeros([L,L])
J_mapped[np.triu_indices(L,1)]=map_index
J_mapped=J_mapped+J_mapped.T
print(J_mapped.shape)


conseq=np.loadtxt('../data/IN/in.nl43.reduce4.seq',dtype='str',ndmin=1)
con=list(conseq[0])
con2=[]
#con2=alphabet_compare2(con)
print(len(con))

wild_list=[]
pos_list=[]
mut_list=[]
reduced_mut_list=[]

my_list1=[]
my_list2=[]
mut1_list=[]
mut2_list=[]
my_pos1_list=[]
my_pos2_list=[]
mut2_alpha_list=[]

f=open('in.nl43.ddE.all.txt',"w")


for c in con:
    con2.append(alphabet_compare2(c))

for i in range(0,len(con)):
    posi=i
    J_get=J_mapped[posi,:]
    gamma=con[posi]
    alpha_list=other_alphabet(con[posi])
    dE_tot2=0
    for alpha in alpha_list:
        dE_tot2=0
        for j,co in enumerate(J_get):
            J_ij=J[int(co)]
        
            if(posi<j):
                E1=J_ij[(con2[posi]*4)+con2[j]]
                E2=J_ij[alphabet_compare1(alpha)+con2[j]]
                dE=E1-E2
                #dE=E2-E1
        
            elif(posi==j):
                continue

            elif(posi>j):
                E1=J_ij[(con2[j]*4)+con2[posi]]
                E2=J_ij[(con2[j]*4)+alphabet_compare2(alpha)]
                dE=E1-E2
                #dE=E2-E1
        
            dE_tot2=dE_tot2+dE
        
        ttdl=gamma+str(posi+1)+alpha
        mut2_alpha_list.append(alpha)
        mut2_list.append(ttdl)

        my_pos2_list.append(i)
        my_list2.append(dE_tot2)
#print my_pos1_list
#print my_pos2_list
#print mut2_list
#print mut1_list

for i,dE1 in enumerate(my_list2):
    posi=my_pos2_list[i]
    
    for j,dE2 in enumerate(my_list2):
        if (j>i):
            posj=my_pos2_list[j]

            J_ij=J[int(J_mapped[posi,posj])]

            if(posi==posj):
                continue

            elif(posi<posj):
                m=con2[posi]*4
                n=con2[posj]
                m_prime=alphabet_compare1(mut2_alpha_list[i])
                n_prime=alphabet_compare2(mut2_alpha_list[j])

            elif(posi>posj):
                m=con2[posj]*4
                n=con2[posi]
                m_prime=alphabet_compare1(mut2_alpha_list[j])
                n_prime=alphabet_compare2(mut2_alpha_list[i])



            dE_sglmut1=dE1
            dE_sglmut2=dE2
        
            deldelE=J_ij[m+n_prime]+J_ij[m_prime+n]-J_ij[m+n]-J_ij[m_prime+n_prime]
 
            #deldelE=J_ij[m+n]+J_ij[m_prime+n_prime]-J_ij[m+n_prime]-J_ij[m_prime+n]

            print(f"{mut2_list[i]}-{mut2_list[j]}", deldelE, file=f)

f.close()      

