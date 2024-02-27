import itertools as it
import scipy as sc
import copy as cp
import numpy as np
from multiprocessing import Process,Manager

FILTER_P_THRES=0.01
FILTER_TOPN_PERF_GENE=5
MAX_PROCESS_NUM=71

def gen_all_comb(max_len,res,head_need_comb):
    res.extend(list(it.combinations(head_need_comb,max_len)))
def cal_all_gene(offset,gene_list,all_comb,cli_num,label,gene,res):
    for g in range(offset,len(gene_list),MAX_PROCESS_NUM):
        print('process id,gene id',offset,g)
        k_gene=gene_list[g]
        v_gene=gene[k_gene]
        f=open('res/'+str(k_gene)+'.txt','w')
        f.write('gene_name\tcomb_item_num\tcomb_item\tspearman coe\tspearman pvalue\n')
        all_src=[]
        tmp_data=np.zeros(cli_num)
        for c in all_comb:
            for i in c:
                tmp_data=tmp_data+np.array(label[i])
            coe,pvalue=sc.stats.spearmanr(v_gene,tmp_data)
            all_src.append((coe,pvalue)+c)
        print('start rank')
        rank=[]
        for src in all_src:
            if src[1]<FILTER_P_THRES:
                rank.append(src)
        rank.sort(key=lambda x:abs(x[0]),reverse=True)
        for i_top in range(min(FILTER_TOPN_PERF_GENE,len(rank))):
            tmp_comb=rank[i_top]
            f.write(k_gene+'\t'+str(len(tmp_comb)-2)+'\t"'+str(tmp_comb[2:])+'"\t'+str(tmp_comb[0])+'\t'+str(tmp_comb[1])+'\n')
        f.close()
if __name__=='__main__':
    FNAME='240203NEWinput.txt'
    f=open(FNAME,'r')
    fo=open('res_'+FNAME,'w')
    ls=f.readlines()
    f.close()
    ls=[l[:-1] if l[-1]=='\n' else l for l in ls]

    FIRST_LABEL_NAME='ACR Malar Rash'

    head=ls[0].split('\t')
    data_list=ls[1:]
    data_list=[l.split('\t') for l in data_list]

    cli_num=len(ls)-1

    first_label_idx=head.index(FIRST_LABEL_NAME)
    head_need_comb=head[first_label_idx:]

    data_titled=[]

    gene={}
    label={}

    for g in range(1,first_label_idx):
        gene[head[g]]=[]
    for g in range(first_label_idx,len(head)):
        label[head[g]]=[]
    for l in data_list:
        for g in range(1,first_label_idx):
            gene[head[g]].append(float(l[g]))
        for g in range(first_label_idx,len(head)):
            label[head[g]].append(int(l[g]))
    print('gen all comb,cal spear man')

    print('gene num:',len(gene))
    print('label num:'+str(len(head_need_comb)))

    res=Manager().list()

    pros=[]
    for k in range(1,len(head_need_comb)+1):
        pros.append(Process(target=gen_all_comb,args=(k,res,head_need_comb)))

    print('start: '+str(len(pros))+' process to gen comb')
    for p in pros:
        p.start()
    for p in pros:
        p.join()
    print('done: '+str(len(pros))+' process to gen comb')

    #merge all comb
    all_comb=tuple(res)
    print('done: merge all comb num:'+str(len(all_comb)))

    gene_list=list(gene.keys())
    pros=[]
    for process_offset in range(MAX_PROCESS_NUM):
       pros.append(Process(target=cal_all_gene,args=(process_offset,gene_list,all_comb,cli_num,label,gene,res)))

    print('start: '+str(len(pros))+' process to cal gene')
    for p in pros:
        p.start()
    for p in pros:
        p.join()
    print('done: '+str(len(pros))+' process to gene')
    
    fo.write('gene_name\tcomb_item_num\tcomb_item\tspearman coe\tspearman pvalue\n')
    for k,v in gene.items():
        tmp_f=open('res/'+k+'.txt','r')
        ls=tmp_f.readlines()
        tmp_f.close()
        if len(ls)>=2 and ls[1]!='':
            for l in ls[1:]:
                fo.write(l)
    print('merge res done')
