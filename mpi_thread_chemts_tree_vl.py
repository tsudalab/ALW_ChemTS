from subprocess import Popen, PIPE
from math import *
import random
import numpy as np
import random as pr
from copy import deepcopy
import itertools
import time
import math
import os
import shutil

import tensorflow as tf

import argparse
import subprocess
from load_model import loaded_model
from keras.preprocessing import sequence
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import Descriptors
from rdkit.Chem import MolFromSmiles, MolToSmiles
from rdkit.Chem import Crippen
import sys
from threading import Thread, Lock, RLock
import threading
from Queue import Queue
from mpi4py import MPI
from RDKitText import tansfersdf
from SDF2GauInput import GauTDDFT_ForDFT
from GaussianRunPack import GaussianDFTRun

import sascorer
import pickle
import gzip
import networkx as nx
from rdkit.Chem import rdmolops


smiles_max_len = 81 # zinc dataset
state_length = 64


class chemical:

    def __init__(self):

        self.position=['&']
    def Clone(self):

        st = chemical()
        st.position= self.position[:]
        return st

    def SelectPosition(self,m):

        self.position.append(m)

    def Getatom(self):
        return [i for i in range(self.num_atom)]

class Node:

    def __init__(self, position = None, parent = None, state = None, nodelock=threading.Lock()):
        self.position = position
        self.parentNode = parent
        self.childNodes = []
        self.child=None
        self.wins = 0
        self.re_max = 0
        self.visits = 0
        self.depth=0
        self.expanded=[]
        self.nodeadded=[]
        self.random_node=[]
        self.all_posible=[]
        self.generate_smile=[]
        self.node_index=[]
        self.valid_smile=[]
        self.new_compound=[]
        self.nodelock=nodelock
        self.ucb=[]
        self.core_id=[]
        self.virtual_loss=0
        self.num_thread_visited=0
        self.all_probs=[]


    def Selectnode(self, ts_strategy, search_parameter, alpha):
        #self.nodelock.acquire()
 
        ucb=[]
        ntv_list=[]
        base_list=[]
        bias_list=[]
        max_list=[]
        for i in range(len(self.childNodes)):
            C = search_parameter
            cNodei = self.childNodes[i]
            if ts_strategy == 'uct': 
                ucb.append(alpha*(cNodei.wins)/(0.0001+cNodei.visits+cNodei.num_thread_visited)+
                           (1-alpha)*cNodei.re_max/(1+cNodei.num_thread_visited)+
                           C*sqrt(2*log(self.visits+self.num_thread_visited)/(0.0001+cNodei.visits+cNodei.num_thread_visited)))
            elif ts_strategy == 'puct':
                prob=self.all_probs[i]
                ucb.append(alpha*(cNodei.wins)/(0.001+cNodei.visits+cNodei.num_thread_visited)+
                           (1-alpha)*cNodei.re_max/(1+cNodei.num_thread_visited)+
                           C*(np.tanh(2*prob-1)+1)/2*sqrt((self.visits+self.num_thread_visited))/(1+cNodei.visits+cNodei.num_thread_visited))
            ntv_list.append(cNodei.num_thread_visited)
            base_list.append(alpha*(cNodei.wins)/(0.001+cNodei.visits+cNodei.num_thread_visited)+(1-alpha)*cNodei.re_max/(1+cNodei.num_thread_visited))
            bias_list.append(ucb[-1] - base_list[-1])
            max_list.append(cNodei.re_max)

        m = np.amax(ucb)
        indices = np.nonzero(ucb == m)[0]
        ind=pr.choice(indices)
        s=self.childNodes[ind]
        return s

    def Addnode(self, m):

        self.nodeadded.remove(m)
        n = Node(position = m, parent = self)
        self.childNodes.append(n)
        return n



    def Update(self, result, add_vis_count = 1):
        self.visits += add_vis_count
        self.wins += result
        if self.re_max < result:
            self.re_max = result

    def delete_virtual_loss(self):
        self.num_thread_visited += -1
        self.virtual_loss=0

    def expanded_node1(self, model, state, val):


        all_nodes=[]

        end="\n"
        position=[]
        position.extend(state)
        total_generated=[]
        new_compound=[]
        get_int_old=[]
        for j in range(len(position)):
            get_int_old.append(val.index(position[j]))

        get_int=get_int_old
        x=np.reshape(get_int,(1,len(get_int)))
        x_pad= sequence.pad_sequences(x, maxlen=82, dtype='int32', padding='post', truncating='pre', value=0.)  #zinc 250,000 
        ex_time=time.time()

        for i in range(1):
            global graph
            with graph.as_default():
                predictions=model.predict(x_pad)
                preds=np.asarray(predictions[0][len(get_int)-1]).astype('float64')
                preds = np.log(preds) / 1.0
                preds = np.exp(preds) / np.sum(np.exp(preds))
                next_probas=np.argsort(preds)[-5:]
                next_probas=list(next_probas)

        if 0 in next_probas:
            next_probas.remove(0)
        all_nodes=next_probas

        self.expanded=all_nodes
        exfi_time=time.time()-ex_time


    def expanded_node(self, model,state,val):
        all_nodes=[]

        end="\n"
        position=[]
        position.extend(state)
        total_generated=[]
        new_compound=[]
        get_int_old=[]
        for j in range(len(position)):
            get_int_old.append(val.index(position[j]))

        get_int=get_int_old
        x=np.reshape(get_int,(1,len(get_int)))
        x_pad= sequence.pad_sequences(x, maxlen=smiles_max_len, dtype='int32',
            padding='post', truncating='pre', value=0.)
        for i in range(60):
            global graph
            with graph.as_default():
                predictions=model.predict(x_pad)
                preds=np.asarray(predictions[0][len(get_int)-1]).astype('float64')
                preds = np.log(preds) / 1.0
                preds = np.exp(preds) / np.sum(np.exp(preds))
                next_probas = np.random.multinomial(1, preds, 1)
                next_int=np.argmax(next_probas)
                all_nodes.append(next_int)

        all_nodes=list(set(all_nodes))
    
        self.expanded=all_nodes


    def expanded_node_puct(self, model,state,val):
        all_nodes=[]

        end="\n"
        position=[]
        position.extend(state)
        total_generated=[]
        new_compound=[]
        get_int_old=[]
        for j in range(len(position)):
            get_int_old.append(val.index(position[j]))

        get_int=get_int_old
        x=np.reshape(get_int,(1,len(get_int)))
        x_pad= sequence.pad_sequences(x, maxlen=smiles_max_len, dtype='int32',padding='post', truncating='pre', value=0.)
        for i in range(1):
            global graph
            with graph.as_default():
                predictions=model.predict(x_pad)
                preds=np.asarray(predictions[0][len(get_int)-1]).astype('float64')
                preds = np.log(preds) / 1.0
                preds = np.exp(preds) / np.sum(np.exp(preds))
                next_probas = np.random.multinomial(1, preds, 1)
                next_int=np.argmax(next_probas)
        
        ordered_preds = np.sort(preds)[::-1]
        ordered_index = np.argsort(preds)[::-1]
        cut_index = 0
        p_sum = 0
        for i in range(len(ordered_preds)):
            p_sum += ordered_preds[i]
            if p_sum > 0.99:
                cut_index = i+1
                break
        all_nodes = ordered_index[:cut_index]
        all_probs = ordered_preds[:cut_index]
        self.expanded=all_nodes
        self.all_probs=all_probs



    def node_to_add(self, all_nodes,val):
        added_nodes=[]
        for i in range(len(all_nodes)):
            added_nodes.append(val[all_nodes[i]])

        self.nodeadded=added_nodes


    def random_node_to_add(self, all_nodes,val):
        added_nodes=[]
        for i in range(len(all_nodes)):
            added_nodes.append(val[all_nodes[i]])

        self.random_node=added_nodes










"""Define some functions used for RNN"""



def chem_kn_simulation(model,state,val,added_nodes):
    all_posible=[]

    end="\n"

    position=[]
    position.extend(state)
    position.append(added_nodes)
    total_generated=[]
    new_compound=[]
    get_int_old=[]
    for j in range(len(position)):
        get_int_old.append(val.index(position[j]))

    get_int=get_int_old

    x=np.reshape(get_int,(1,len(get_int)))
    x_pad= sequence.pad_sequences(x, maxlen=smiles_max_len, dtype='int32',
        padding='post', truncating='pre', value=0.)
    while not get_int[-1] == val.index(end):
        predictions=model.predict(x_pad)
        preds=np.asarray(predictions[0][len(get_int)-1]).astype('float64')
        preds = np.log(preds) / 1.0
        preds = np.exp(preds) / np.sum(np.exp(preds))
        next_probas = np.random.multinomial(1, preds, 1)
        next_int=np.argmax(next_probas)
        a=predictions[0][len(get_int)-1]
        next_int_test=sorted(range(len(a)), key=lambda i: a[i])[-10:]
        get_int.append(next_int)
        x=np.reshape(get_int,(1,len(get_int)))
        x_pad = sequence.pad_sequences(x, maxlen=smiles_max_len, dtype='int32',
            padding='post', truncating='pre', value=0.)
        if len(get_int)>state_length:
            break
    total_generated.append(get_int)
    all_posible.extend(total_generated)


    return all_posible




def predict_smile(all_posible,val):
    new_compound=[]
    for i in range(len(all_posible)):
        total_generated=all_posible[i]

        generate_smile=[]

        for j in range(len(total_generated)-1):
            generate_smile.append(val[total_generated[j]])
        generate_smile.remove("&")
        new_compound.append(generate_smile)

    return new_compound


def make_input_smile(generate_smile):
    new_compound=[]
    for i in range(len(generate_smile)):
        middle=[]
        for j in range(len(generate_smile[i])):
            middle.append(generate_smile[i][j])
        com=''.join(middle)
        new_compound.append(com)
    return new_compound




def ChemTS_run(rootnode,result_queue,lock,chem_model,ts_strategy,search_parameter,num_simulations, gau_parallel,simulation_time, output_file,alpha,objective,num_rollout,charge_check,SA_score_check):
    """----------------------------------------------------------------------"""
    """----------------------------------------------------------------------"""
    global maxnum
    global gau_file_index
    global ind_mol
    start_time=time.time()
    while time.time()-start_time<simulation_time:
        node = rootnode
        state=['&']
        """selection step"""
        node_pool=[]
        lock.acquire()

        while len(node.expanded)>0 and node.nodeadded==[] and len(node.childNodes)==len(node.expanded):
            node = node.Selectnode(ts_strategy, search_parameter, alpha)
            state.append(node.position)
        depth.append(len(state))

        """this if condition makes sure the tree not exceed the maximum depth"""
        if len(state)>state_length:
            re=-10


            while node != None:
                node.Update(re)
                node = node.parentNode
            lock.release()
        else:
            """expansion step"""
            m = None
            if node.expanded==[]:
                if ts_strategy == 'uct':
                    node.expanded_node(chem_model,state,val)
                elif ts_strategy == 'puct':
                    node.expanded_node_puct(chem_model,state,val)
                node.node_to_add(node.expanded,val)
                node.random_node_to_add(node.expanded,val)
                
                if node.nodeadded!=[]:
                    
                    m=node.nodeadded[0]
            else:
                if node.nodeadded!=[]:
                    m=node.nodeadded[0]

            if m == None:
                m = val[random.choice(node.expanded)]
                print('randomly selected')
            else:
                if m != '\n':
                    node = node.Addnode(m)
                else:
                    node.nodeadded.remove(m)
                    lock.release()   
                    continue

            lock.release()
	   
	    """simulation step"""
            for ro in range(num_rollout):
                                
                lock.acquire()
                """add virtual loss"""
                node_tmp = node
                while node_tmp != None:
                    node_tmp.num_thread_visited+=1
                    node_tmp = node_tmp.parentNode
                print 'rootnode.num_thread_visited', rootnode.num_thread_visited
                lock.release()

                lock.acquire()
                maxnum+=1
                ind_mol+=1
                print('free_core_id_prev', len(free_core_id),'use_core_id', len(use_core_id))

                dest_core=random.choice(free_core_id)
                use_core_id.append(dest_core)
                free_core_id.remove(dest_core)
                print('dest_core', dest_core)
                try:
                    comm.send([state,m,ind_mol], dest=dest_core, tag=START)
                    lock.release()
                except:
                    print('comm.send failed')
                    free_core_id.append(dest_core)
                    use_core_id.remove(dest_core)
                    lock.acquire()
                    """backpropation step"""
                    while node!= None:
                        node.Update(0, add_vis_count = 0)
                        node.delete_virtual_loss()
                        node = node.parentNode
                    lock.release()
                    
                    continue

                try:
                    data = comm.recv(source=dest_core, tag=MPI.ANY_TAG, status=status)   

                    lock.acquire()
                    free_core_id.append(data[2])
                    use_core_id.remove(data[2])
                    print('data[2]', data[2], 'dest_core', dest_core)
                    lock.release()
                except:
                    print('comm.recv failed.')
                    lock.acquire()
                    free_core_id.append(dest_core)
                    use_core_id.remove(dest_core)
                    
                    """backpropation step"""
                    while node!= None:
                        node.Update(0, add_vis_count = 0)
                        node.delete_virtual_loss()
                        node = node.parentNode
                    lock.release()
                    
                    continue
                print('free_core_id', free_core_id)

                tag = status.Get_tag()
                if tag == DONE:
                    lock.acquire()
                    all_compounds.append(data[1])
                    lock.release()
                    if data[0]!=-1000 and data[3] >= 0:
                    

                        if objective == 'WL_IN':
                            if data[3] > 0.1:                                                                                                          
                                re = (np.tanh(0.003*(data[0]-400)) + 1)/2                                                                     
                            else:                                                                                                                      
                                re= (1/(-np.log10(data[3]))) * (np.tanh(0.003*(data[0]-400)) + 1)/2
                        elif objective == 'HL':
                            #HOMO/LUMO                                                                                                                            
                            re = 1 - data[5]/10.
                        elif objective == 'WL':
                            re = (np.tanh(0.003*(data[0]-400)) + 1)/2

                        if SA_score_check:
                            data[11]
                            re = re*((-np.tanh(data[11]-4)+1)/2)
                        lock.acquire()
                        wave_compounds.append(data[1])
                        wave.append(data[0])
                        deen_list.append(data[4])
                        uv_intensity_list.append(data[3])
                        gap_list.append(data[5])
                        wl_list_list.append(data[6])
                        intensity_list_list.append(data[7])
                        reward_list.append(re)
                        index_list.append(data[8])
                        mol_weight_list.append(data[9])
                        logP_list.append(data[10])
                        SA_score_list.append(data[11])
                        depth_list.append(data[12])

                        with open('/home/terayama/csvcom_.csv','wb') as file: 
                            for line1 in wave_compounds:
                                file.write(str(line1))
                                file.write('\n')
                        with open('/home/terayama/csvwave_.csv','wb') as file:
                            for line2 in wave:
                                file.write(str(line2))
                                file.write('\n')
                    
                        with open('/home/terayama/'+output_file,'wb') as file:
                        
                            file.write('#Search strategy, '+ts_strategy+', search_parameter, '+str(search_parameter)+', alpha, '+str(alpha)+', objective,'+str(objective)+', parallel simulations, '+str(num_simulations)+', gaussian parallel, '+str(gau_parallel)+', simulation_time (h), '+str(simulation_time/3600)+', num_rollout, '+str(num_rollout)+'charge_check, '+str(charge_check)+'SA_score_check, '+str(SA_score_check)+'\n')
                            file.write('#Compound, index, wavelength, uv_intensity, reward,  deen, gap, mol_weight, logP, SA_score, TS depth,  wavelength_list, intensity_list \n')
                            for i in range(len(wave_compounds)):
                                file.write(str(wave_compounds[i])+', ')
                                file.write(str(index_list[i])+', ')
                                file.write(str(wave[i])+', ')
                                file.write(str(uv_intensity_list[i])+', ')
                                file.write(str(reward_list[i])+', ')
                                file.write(str(deen_list[i])+', ')
                                file.write(str(gap_list[i])+', ')
                                file.write(str(mol_weight_list[i])+', ')
                                file.write(str(logP_list[i])+', ')
                                file.write(str(SA_score_list[i])+', ')
                                file.write(str(depth_list[i])+', ')
                                for wl_i in wl_list_list[i]:
                                    file.write(str(wl_i)+', ')
                                for int_i in intensity_list_list[i]:
                                    file.write(str(int_i)+', ')
                                file.write('\n')
                            
                        lock.release()
                    if data[0]==-1000:
                        re=0
                    if data[3]<0:
                        re=0

                lock.acquire()
                if re == None:
                    re = 0
                """backpropation step"""
                while node!= None:
                    #print "node.parentNode:",node.parentNode
                    if re <= 0:
                        node.Update(re, add_vis_count = 0)
                    else:
                        node.Update(re, add_vis_count = 1)
                    node.delete_virtual_loss()
                    node = node.parentNode
                lock.release()


    result_queue.put([all_compounds,wave_compounds,depth,wave,maxnum,uv_intensity_list,deen_list,gap_list,reward_list,index_list,mol_weight_list,logP_list,SA_score_list,depth_list])

def charge_check(mol):
    print 'charge_checking'
    standard_valence_list = [0, 1, 0, 1, 2, 3, 4, 3, 2, 1, 0, 1, 2, 3, 4, 3, 2, 1, 0, 1, 2]
    check = True
    for atom in mol.GetAtoms():
        if standard_valence_list[atom.GetAtomicNum()] != atom.GetExplicitValence():
            check = False
            break
    return check

def gaussion_workers(chem_model,val,gau_parallel,charge_check):
    while True:
        simulation_time=time.time()
        task = comm.recv(source=0, tag=MPI.ANY_TAG, status=status)
        tag = status.Get_tag()
        if tag==START:
            state=task[0]
            m=task[1]
            ind=task[2]
            all_posible=chem_kn_simulation(chem_model,state,val,m)
            generate_smile=predict_smile(all_posible,val)
            new_compound=make_input_smile(generate_smile)
            score=[]
            kao=[]
            intensity = -1000000
            deen = 1000000
            gap = 1000000
            mol_weight = 0
            SA_score = 10
            logP = 0
            dp = len(state)
            intensity_list = []
            wl_list = []
            standard_valence_list = [0, 1, 0, 1, 2, 3, 4, 3, 2, 1, 0, 1, 2, 3, 4, 3, 2, 1, 0, 1, 2]

            try:
                m = Chem.MolFromSmiles(str(new_compound[0]))
                mol_weight = Descriptors.MolWt(m)
                logP = Crippen.MolLogP(m)
                SA_score = sascorer.calculateScore(m)
                m_H = Chem.AddHs(m)
                standard_valence_list = [0, 1, 0, 1, 2, 3, 4, 3, 2, 1, 0, 1, 2, 3, 4, 3, 2, 1, 0, 1, 2]
                ccheck = True
                if charge_check:
                    for atom in m_H.GetAtoms():
                        if standard_valence_list[atom.GetAtomicNum()] != atom.GetExplicitValence():
                            ccheck = False
                            break

                if not ccheck:
                    m = None

            except:
                m=None


            if m!=None:
                try:
                    stable=tansfersdf(str(new_compound[0]),ind)
                except:
                    stable = -1
                if stable==1.0:
                    try:
			SDFinput = 'CheckMolopt'+str(ind)+'.sdf'
                        calc_sdf = GaussianDFTRun('B3LYP', '3-21G*', gau_parallel, 'OPT energy deen nmr uv homolumo', SDFinput, 0)
                        outdic = calc_sdf.run_gaussian()
                        wavelength = outdic['uv'][0]
                        
                        if os.path.isfile('CheckMol'+str(ind)+'.sdf'):
                            shutil.move('CheckMol'+str(ind)+'.sdf', 'dft_result')
                        if os.path.isfile('CheckMolopt'+str(ind)+'.sdf'):
                            shutil.move('CheckMolopt'+str(ind)+'.sdf', 'dft_result')
                        if os.path.isfile('CheckMolopt'+str(ind)+'.com'):
                            shutil.move('CheckMolopt'+str(ind)+'.com', 'dft_result')
                        if os.path.isfile('CheckMolopt'+str(ind)+'.log'):
                            shutil.move('CheckMolopt'+str(ind)+'.log', 'dft_result')
                        if os.path.isfile('CheckMolopt'+str(ind)+'.chk'):
                            os.remove('CheckMolopt'+str(ind)+'.chk')
                    except:
                        wavelength=None
                        if os.path.isfile('CheckMolopt'+str(ind)+'.sdf'):
                            os.remove('CheckMolopt'+str(ind)+'.sdf')
                        if os.path.isfile('CheckMol'+str(ind)+'.sdf'):
                            os.remove('CheckMol'+str(ind)+'.sdf')
                else:
                    wavelength=None
                    if os.path.isfile('CheckMolopt'+str(ind)+'.sdf'):
                        os.remove('CheckMolopt'+str(ind)+'.sdf')
                    if os.path.isfile('CheckMol'+str(ind)+'.sdf'):
                        os.remove('CheckMol'+str(ind)+'.sdf')

                if wavelength!=None and wavelength!=[]:
                    wavenum=wavelength[0]
                    intensity=outdic['uv'][1][0]
                    deen=outdic['deen']
                    gap=outdic['gap']
                    wl_list = outdic['uv'][0]
                    intensity_list = outdic['uv'][1]
                else:
                    wavenum=-1000
            else:
                wavenum=-1000
            score.append(wavenum)
            score.append(new_compound[0])
            score.append(rank)
            score.append(intensity)
            score.append(deen)
            score.append(gap)
            score.append(wl_list)
            score.append(intensity_list)
            score.append(ind)
            score.append(mol_weight)
            score.append(logP)
            score.append(SA_score)
            score.append(dp)

            comm.send(score, dest=0, tag=DONE)
            simulation_fi_time=time.time()-simulation_time
            print("simulation_fi_time:",simulation_fi_time)
        if tag==EXIT:
            MPI.Abort(MPI.COMM_WORLD)

    comm.send([-1000,'',0,0,0,0,[],[],ind,0,0,0,0], dest=0, tag=EXIT)



if __name__ == "__main__":
    comm=MPI.COMM_WORLD
    size=comm.size
    rank=comm.rank
    status=MPI.Status()
    READY, START, DONE, EXIT = 0, 1, 2, 3

    val=['\n', '&', 'C', '(', ')', 'c', '1', '2', 'o', '=', 'O', 'N', '3', 'F', '[C@@H]', 'n', '-', '#', '/', '[nH]', 'Br', '[C@H]', 'Cl', '[C@]', '[C@@]', '\\', '4', '5', '6', '7', 'I']

    chem_model=loaded_model()
    graph = tf.get_default_graph()
    chemical_state = chemical()

    ts_strategy = 'puct' #'uct', 'puct' 
    search_parameter = 2 #If ts_strategy=='uct', 0 < search_parameter < 1. If ts_strategy=='puct', default value is 5 (AlphaGo). 
    num_simulations = 2048 # core - 1, max: 2560 (skylake)
    gau_parallel = 1
    num_rollout = 3
    simulation_time = 3600*120 # 3600*24 # max: 168h
    alpha = 1 # alph*mean + (1 - alpha)*max + bais
    objective = 'WL' # 'WL_IT', 'HL', 'WL'
    charge_check = True # True or False
    SA_score_check = True # True or False
    output_file = 'csvresult_FP_'+ts_strategy+'_C'+str(search_parameter)+'_alpha'+str(alpha)+'_obj'+objective+'_para'+str(num_simulations)+'_time'+str(simulation_time/3600)+'h_rollout'+str(num_rollout)+'_CC'+str(charge_check)+'_SA'+str(SA_score_check)+'_1101.csv'

    thread_pool=[]
    lock=Lock()
    gau_file_index=0

    """initialization of the chemical trees and grammar trees"""
    root=['&']
    rootnode = Node(position= root)
    maxnum=0
    ind_mol=0
    reward_dis=[]
    all_compounds=[]
    wave_compounds=[]
    wave=[]
    deen_list = []
    gap_list = []
    uv_intensity_list = []
    wl_list_list = []
    intensity_list_list = []
    reward_list = []
    index_list = []
    mol_weight_list = []
    logP_list = []
    SA_score_list = []
    depth_list = []
    depth=[]
    result=[]
    result_queue=Queue()
    free_core_id=range(1,num_simulations+1)
    use_core_id = []

    if rank==0:
        for thread_id in range(num_simulations):
            thread_best = Thread(target=ChemTS_run,args=(rootnode,result_queue,lock,chem_model,ts_strategy,search_parameter,num_simulations, gau_parallel,simulation_time,output_file,alpha,objective,num_rollout,charge_check,SA_score_check))
            thread_pool.append(thread_best)

        for i in range(num_simulations):
            thread_pool[i].start()

        for i in range(num_simulations):
            thread_pool[i].join()
        for i in range(num_simulations):
            result.append(result_queue.get())
        comm.Abort()
        for i in range(len(free_core_id)):
            comm.send(None, dest=i+1, tag=EXIT)
    else:
	    gaussion_workers(chem_model,val,gau_parallel, charge_check)
