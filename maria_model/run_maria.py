#####################import#################################
#basic
from __future__ import print_function
import numpy as np
import cPickle as pickle
#keras related
import keras
from keras.models import model_from_json
print(keras.__version__)
#ml and math related
from sklearn.metrics import roc_auc_score
from random import shuffle
import random
from scipy.stats import percentileofscore
# %matplotlib inline
# import matplotlib
# import matplotlib.pyplot as plt

####################Get relative path####################################
path_additional = "/".join(__file__.split("/")[:-1])

#####################path and dictionary directory########################
path_dict = path_additional+'/supporting_file/'
print('loading data from '+path_dict)
# path_save = '/cstor/stanford/rbaltman/users/bchen45/mcl_data/model_weight/'
# path_encoding = '/home/stanford/rbaltman/users/bchen45/code/python_general/encoding_dict/'
# path_final_model = '/cstor/stanford/rbaltman/users/bchen45/mcl_data/final_model/'
# path_mcl2 = '/home/stanford/rbaltman/users/bchen45/data/MCL_data/'
# path_list = '/cstor/stanford/rbaltman/users/bchen45/mcl_data/supporting_data/'
# 
# one_gene_path = '/home/stanford/rbaltman/users/bchen45/data/protein_general/human_proteinome_oneline.str'


dict_aa = pickle.load(open(path_dict+'aa_21_sparse_encoding.dict','r'))
dict_aa['-'] = [0]*19 + [-1,0]

onegenestr = pickle.load(open(path_dict+'human_proteinome_oneline.str'
,'r'))
dict_gene2seq = pickle.load(open(path_dict+'gene2seq.dict','r'))
cleavage_default = 'AAAAAAAKP--' #0.39 for probability
dict_name = 'NetMHC_DR.dict'
dict_dr = pickle.load(open(path_dict+dict_name,'r'))
dict_mcl_rna_blood = pickle.load(open(path_dict+'dict_mcl_rna_blood.dict','r'))

print('Start loading cached binding score distribution')
dict_binding_distr = '' #pickle.load(open(path_dict+'binding_distr.dict','r')) 
print('Loaded cached binding score distribution')


################## Parameters ######################################
fix_len = 19
len_mhc = 19
print('Each DR allele presented by a '+str(fix_len)+'-AA long pseudosequence')
MAXLEN = 26
print('The maximum length allowed for MARIA MHC-II is '+str(MAXLEN))
batch0 = 1024*16
char_len = 21
len_char = 21


#####################model and weight #################################
#function to import keras NN model
def import_model(path_model, model_name0,weight_name0):
    model_name0 = path_model+ model_name0
    weight_name0 = path_model + weight_name0
    model0 = model_from_json(open(model_name0).read())
    model0.load_weights(weight_name0)
    return model0

#path 
path_model = path_additional+'/models/'

#binding model
#path_save = '/cstor/stanford/rbaltman/users/bchen45/mcl_data/model_weight/'
model1 = 'rnn_merge_all_netmhc_data_n64_sparse_d0.45_dronly.json'
weight1 = 'rnn_merge_all_netmhc_data_n64_sparse_d0.45_dronly.h5'
model_binding = import_model(path_model,model1,weight1)
model_binding.compile(loss='binary_crossentropy', optimizer='adam')

#cleavage model
model1 = 'nn_classifier_cleavage_n32_updown6_d0.5.json'
weight1 = 'nn_classifier_cleavage_n32_updown6_d0.5.weight'
model_cleavage = import_model(path_model,model1,weight1)
model_cleavage.compile(loss='binary_crossentropy', optimizer='adam')

#merge model MARIA
model1 = 'rnn_withiedb_4feature_merge_d0.4_934auc.json'
weight1 = 'rnn_withiedb_4feature_merge_d0.4_934auc.h5'
model_merge = import_model(path_model,model1,weight1)
model_merge.compile(loss='categorical_crossentropy', optimizer='adam')
#model_merge.summary()


################ NN related Encoding functions ################ 
def encoding_y(list0,class0=2):
    list_out = []
    for x in list0:
        vec0 = [0]*class0
        vec0[int(x)] = 1
        list_out.append(vec0)
    return list_out

def encoding_line(str0, max_len, char_len, dict_aa, classn = 2):
    coded0 = np.zeros((max_len,char_len))
    for i,char0 in enumerate(str0):
        coded0[i,:] = dict_aa[char0] 
    #print(str0)
    #print(coded0)
    return coded0

def encoding(matrix0, input0, len0,dict_aa, char_len):
    for i, sentence in enumerate(input0):
        matrix0[i] = encoding_line(sentence, len0, char_len,dict_aa)
    return matrix0

def encoding_data(list0,MAXLEN=12,dict_aa=dict_aa,char_len=21):
    #encoding   
    X_0_m = np.zeros((len(list0), MAXLEN, char_len))
    X_encoded = encoding(X_0_m,list0, MAXLEN, dict_aa, char_len)
    return X_encoded

def encoding_data(list0,MAXLEN=MAXLEN,dict_aa=dict_aa,char_len=char_len ):
    #encoding   
    X_0_m = np.zeros((len(list0), MAXLEN, char_len))
    X_encoded = encoding(X_0_m,list0, MAXLEN, dict_aa, char_len)
    return X_encoded

def encoding_fixed(list0,MAXLEN=MAXLEN,dict_aa=dict_aa,char_len=char_len,add_placeh=0):
    X_0_m = np.zeros((len(list0), MAXLEN, char_len))
    X_encoded = encoding(X_0_m,list0,MAXLEN, dict_aa,char_len)
    X_encoded = X_encoded.reshape((X_encoded.shape[0],MAXLEN*char_len))
    if add_placeh > 0:
        add_block = np.ones((X_encoded.shape[0],add_placeh))
        print(X_encoded.shape)
        X_encoded = np.concatenate((X_encoded,add_block),axis=1)
        print(X_encoded.shape)
    return np.array(X_encoded)

def encode_apair(data_pair,dict_aa=dict_aa,len_mhc=fix_len,MAXLEN=MAXLEN,char_len=char_len):
    x_mhc = encoding_fixed(data_pair[0],len_mhc,dict_aa,char_len,add_placeh=0)
    x_seq = encoding_data(data_pair[1],MAXLEN,dict_aa,char_len)
    return [x_mhc, x_seq]


def get_binding_rank(list_predicted,pep_list,mhc0_pep,dict_distr):
    list_out = []
    for score0, pep0 in zip(list_predicted,pep_list):
        len0 = len(pep0)
        rank0 = percentileofscore(dict_distr[mhc0_pep][len0],score0)/float(100)
        list_out.append(rank0)
    return list_out

def run_MARIA_binding(model0,mhc_pep,pep_list,dict_aa=dict_aa,dict_dr=dict_dr,dict_distr='',
              len_mhc=fix_len,MAXLEN=MAXLEN,char_len=char_len,batch0=batch0,
                      translate_dr=True,output_ranking=False):
    if translate_dr:
        mhc_pep = dict_dr[mhc_pep]
    data0 = [[mhc_pep]*len(pep_list),pep_list]
    #print(data0)
    x = encode_apair(data0,dict_aa,len_mhc,MAXLEN,char_len)
    #print(x)
    list_predicted = model0.predict(x,verbose=0,batch_size=batch0)
    if output_ranking:
        list_predicted = get_binding_rank(list_predicted,pep_list,mhc_pep,dict_distr)
    return list_predicted

#gene expresion normalization
def convert_gn2exp(list_gn,dict_rna,scaling=10,def_val=5):
    list_out = []
    for gn0 in list_gn:
        if gn0 in dict_rna:
            list_out.append(dict_rna[gn0]+0.001)
#             if dict_rna[gn0] < 0.001:
#                 print(dict_rna[gn0])
        else:
            list_out.append(def_val)
    list_out = np.log10(list_out)/float(def_val)
    return list_out

#cleavage related functions
def get_cleave_ingene(frag0,seq0,flanking=6):
    out = []
    seq0 = '------'+seq0+'------'
    index0 = seq0.find(frag0)
    #print(index0)
    if index0 > -1:
        part1 = seq0[index0-flanking:index0]
        part2 = seq0[index0+len(frag0):index0+len(frag0)+flanking]
        #print(part1)
        #print(part2)
        out = part1 + part2
    return out

def get_cleave_inall(frag0,seq0,flanking=6):
    out = []
    index0 = seq0.find(frag0)
    #print(index0)
    if index0 > -1:
        part1 = seq0[index0-flanking:index0]
        part2 = seq0[index0+len(frag0):index0+len(frag0)+flanking]
        #print(part1)
        #print(part2)
        out = part1 + part2
    return out

def get_cleave_scores(model_cleavage,pep_list,gn_list,oneline=onegenestr, 
                      flanking=6, cleavage_default = 'AAAAAAAKP--',char_len=21,
                      dict_aa=dict_aa):
    cle_list = []
    n_unfound = 0
    #get cleavage sites
    for pep0, gn0 in zip(pep_list,gn_list):
        if gn0 in dict_gene2seq:
            seq0 = dict_gene2seq[gn0]
            out = get_cleave_ingene(pep0,seq0,flanking=flanking)
        if out == []:
            out = get_cleave_inall(pep0,oneline,flanking=flanking)
        if out == []:
            out = cleavage_default
            n_unfound += 1
        cle_list.append(out)
    #run model on cleavage sites
    x_predict = encoding_data(cle_list,MAXLEN=flanking*2,char_len=char_len,
                              dict_aa=dict_aa).reshape(len(cle_list),flanking*2*char_len)
    scores = model_cleavage.predict(x_predict)[:,0]
    print(n_unfound)
    return scores



#################### Running MARIA ####################
def get_binding_rank(list_predicted,pep_list,mhc_list,dict_distr):
    list_out = []
    if len(mhc_list[0]) == 1 :
        mhc_list = [mhc_list]*len(pep_list)
    for score0, pep0,mhc0_pep in zip(list_predicted,pep_list,mhc_list):
        len0 = max(len(pep0),8)
        rank0 = percentileofscore(dict_distr[mhc0_pep][len0],score0)/float(100)
        list_out.append(rank0)
    return list_out

#this version is for generating training too
#in this version, mhc_list is previously translated and has the same size of pep_list
#output encoding sequence matrix too
def run_MARIA_binding_mhclist(model0,mhc_list,pep_list,dict_aa=dict_aa,dict_dr=dict_dr,dict_distr='',
              len_mhc=fix_len,MAXLEN=MAXLEN,char_len=char_len,batch0=batch0,
                      output_ranking=False):
    data0 = [mhc_list,pep_list]
    #print(data0)
    x = encode_apair(data0,dict_aa,len_mhc,MAXLEN,char_len)
    #print(x)
    list_predicted = model0.predict(x,verbose=0,batch_size=batch0)[:,0]
    #print(list_predicted)
    if output_ranking:
        list_predicted = get_binding_rank(list_predicted,pep_list,mhc_list,dict_distr)
    return list_predicted, x[1]

def add_binding_to_dataset(train_list,model_binding,model_cleavage,ranking=True):
    mhc0_list = []
    mhc1_list = []
    seq_list = train_list[1]
    gene_list = []
    
    for fixed0 in train_list[0]:
        allele0 = fixed0[0]
        allele1 = fixed0[1]
        gene_list.append(fixed0[2])
        mhc0_list.append(allele0)
        mhc1_list.append(allele1)
        
    score0, _ = run_MARIA_binding_mhclist(model_binding,mhc0_list,seq_list,dict_distr=dict_binding_distr,
                 output_ranking=ranking)
    score1, seq_encoded_list = run_MARIA_binding_mhclist(model_binding,mhc1_list,seq_list,dict_distr=dict_binding_distr,
                 output_ranking=ranking)
    score_exp = convert_gn2exp(gene_list,dict_mcl_rna_blood)
    cleavage_score = get_cleave_scores(model_cleavage,seq_list,gene_list)
    #print(len(score_exp))
    score_list = np.array([score0,score1,score_exp,cleavage_score]).T
    print(len(score_list),len(seq_encoded_list))
    train_new = [score_list,seq_encoded_list]
    return train_new

##checking the model on training data
def check_model(model_merge, x_train_pos, x_train_neg,batch0 = 1024,
                record_file=1,cut_off0=0.5,average0='macro'):
    list_predicted_pos = model_merge.predict(x_train_pos,verbose=1,batch_size=batch0)[:,1]
    #list_predicted_pos1 = model_merge.predict(x_train_pos1,verbose=1,batch_size=batch0)[:,1]
    list_predicted_neg = model_merge.predict(x_train_neg,verbose=1,batch_size=batch0)[:,1]
    #np.sum(list_predicted_neg < 0.5)/float(len(list_predicted_neg))
    #list_predict_pos = [max(x1,x2) for x1,x2 in zip(list_predicted_pos0,list_predicted_pos1)]
    #list_predict_scores = list_predict_pos + list_predicted_neg
    #list_true = [1]*len(list_predict_pos) + [0] * list_predicted_neg
    #print(np.sum(list_predict_pos >= 0.5)/float(len(list_predict_pos)))
    print('Sensitivity:')
    sens0 = np.sum(np.array(list_predicted_pos) >= cut_off0)/float(len(list_predicted_pos))
    print(sens0)
    print('Specificity:')
    spec0 = 1-np.sum(np.array(list_predicted_neg) >= cut_off0)/float(len(list_predicted_neg))
    print(spec0)
    if record_file != 1:
        record_file.write('Sensitivity='+str(sens0)+' Specificity='+str(spec0)+' ')
    #plt.scatter(list_predicted_pos0,list_predicted_pos1)
    auc_model = cal_auc_pos_neg(list_predicted_pos,list_predicted_neg,average0=average0)
    return auc_model

def make_mhc_list(mhc_list,pep_list):
    mhc_list = mhc_list.replace(' ','')
    mhc_list = mhc_list.split(',')
    if len(mhc_list) == 1:
        if mhc_list[0] in dict_dr:
            dr_seq = dict_dr[mhc_list[0]]
            out = [[dr_seq,dr_seq]]*len(pep_list)
        else:
            out = 'DR allele not recognized.'
    elif len(mhc_list) == 2:
        if mhc_list[0] in dict_dr and mhc_list[1] in dict_dr:
            dr_seq0 = dict_dr[mhc_list[0]]
            dr_seq1 = dict_dr[mhc_list[1]]
            out = [[dr_seq0,dr_seq1]]*len(pep_list)
        else:
            out = 'DR alleles not recognized.'
    else:
        out = "Please give one DR allele or two DR alleles separated by a comma only (e.g. HLA-DRB1*07:01,HLA-DRB1*01:01)"
    return out

def predict_with4(pep_list,mhc_list,gn_list,
                  model_merge=model_merge,model_cleavage=model_cleavage,
                  model_binding=model_binding,ranking=False):
    #make the data into the right format
    data0 = [[],[]]
    #print(pep_list)
    for mhc0,gn0,pep0 in zip(mhc_list,gn_list,pep_list):
        data0[0].append([mhc0[0],mhc0[1],gn0])
        data0[1].append(pep0)
    data0_encoded = add_binding_to_dataset(data0,model_binding,model_cleavage,ranking=ranking)
    scores = model_merge.predict(data0_encoded)[:,1]
    return scores

def slide_gene(str0,len0=15):
    list_out = []
    for i in range(0,len(str0)-len0):
        list_out.append(str0[i:i+len0])
    return list_out

def scan_gene(gn0,len0=15):
    out0 = 'Unknown gene'
    if gn0 in dict_gene2seq:
        pep_full = dict_gene2seq[gn0]
        if gn0 in dict_mcl_rna_blood:
            print(dict_mcl_rna_blood[gn0])
        pep_list = slide_gene(pep_full,len0=len0)
        gn_list = [gn0]*len(pep_list)
        out0 = [gn_list,pep_list]

    return out0


        
        