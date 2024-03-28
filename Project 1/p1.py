from numpy import * 
import numpy as np
from Bio import SeqIO
import sys
import configparser
# dp cell
class  cell:
  def __init__(self, D, S,I):
    self.D = D
    self.S = S
    self.I= I
# for alignment storing
class alignment:
  def __init__(self, s1_c, dir,s2_c):
    self.s1_c = s1_c
    self.dir = dir
    self.s2_c= s2_c


def similarity(a,b):
  return match_score if a==b else mismatch_score

# reading config file
def read_config(filename):
    config_file = configparser.ConfigParser()
    config_file.read(filename)
    global match_score, mismatch_score, h, g
    h=int(config_file['params']['h'])
    g=int(config_file['params']['g'])
    match_score=int(config_file['params']['match'])
    mismatch_score=int(config_file['params']['mismatch'])
    
# loading sequence
def load_sequences(filename):
    sequences= list(SeqIO.parse(filename, "fasta"))
    global s1,s2,m,n
    s1=sequences[0].seq
    s2=sequences[1].seq
    if len(sequences[0].seq)>len(sequences[1].seq):
      s1=sequences[0].seq
      s2=sequences[1].seq
    else:
      s1=sequences[1].seq
      s2=sequences[0].seq


    m,n=len(s1),len(s2)
    if n>5000:
      m,n=5000,5000
    print("s1 length:",m)
    print("s1 length:",n)


# global alignment algorithm
def global_alignment():
  print('Scroes: match= ',match_score,'  Mismatch: ',mismatch_score,' h: ',h,' g:',g,file=f_global)
  print('',file=f_global)
  print('Sequence 1 = "s1", length = ',m, 'characters',file=f_global)
  print('Sequence 2 = "s2", length = ',n, 'characters',file=f_global)
  dp=[]  
  for i in range(m+1):
      col =[]
      for j in range(n+1):
        if i==0 and j==0:
          col.append(cell(0,0,0))
        elif i==0 and j!=0:
          col.append(cell(-99999,-99999,h+(j*g)))
        elif i!=0 and j==0:
          col.append(cell((h+(i*g)),-99999,-99999))
        else:
          col.append(cell(0,0,0))  
      
      dp.append(col)
  print(shape(dp))
  for i in range(m+1):
    for j in range(n+1):
      if i==0 or j==0:
        continue
      dp[i][j].D=max(int(dp[i-1][j].D)+g, int(dp[i-1][j].S)+h+g, int(dp[i-1][j].I)+h+g)
      dp[i][j].I=max(int(dp[i][j-1].D)+h+g, int(dp[i][j-1].S)+h+g, int(dp[i][j-1].I)+g)
      dp[i][j].S=max(dp[i-1][j-1].S, dp[i-1][j-1].D, dp[i-1][j-1].I)+similarity(s1[i-1],s2[j-1]) 

  print('',file=f_global)
  print('Global Optimal Score: ',max(dp[m][n].S, dp[m][n].D, dp[m][n].I ),file=f_global)
  return dp

# retrace algorithm for global alignment
def retrace_global(dp):
  cur_i,cur_j = m,n
  dir=np.argmax([dp[m][n].S,dp[m][n].D,dp[m][n].I]) # 0: subsitution, 1: deletion 2: insertion
  matctes=0
  mismatches=0
  opening_gap=0
  gap_extension=0
  alignment_store=[]

  while cur_i != 0  and cur_j != 0:
    if dir==1: # deleteion case
      if dp[cur_i][cur_j].D==dp[cur_i-1][cur_j].D+g:
        dir=1
      elif  dp[cur_i][cur_j].D==dp[cur_i-1][cur_j].S+h+g:
        dir=0
        opening_gap=opening_gap+1
      else:
        dir=2
        opening_gap=opening_gap+1
      
      alignment_store.append(alignment(s1[cur_i-1],' ','-'))
      cur_i=cur_i-1
      gap_extension=gap_extension+1  
    elif dir==0: # sub case
      if s1[cur_i-1]==s2[cur_j-1]:
        matctes=matctes+1
        pen=match_score
        alignment_store.append(alignment(s1[cur_i-1],'|',s2[cur_j-1]))
      else:
        mismatches=mismatches+1
        pen=mismatch_score
        alignment_store.append(alignment(s1[cur_i-1],' ',s2[cur_j-1]))
      if dp[cur_i][cur_j].S==dp[cur_i-1][cur_j-1].D +pen:
        dir=1
      elif dp[cur_i][cur_j].S==dp[cur_i-1][cur_j-1].S+pen:
        dir=0
      else:
        dir=2
      cur_i=cur_i-1
      cur_j=cur_j-1
    else: # insertion case
      if dp[cur_i][cur_j].I==dp[cur_i][cur_j-1].D +g+h:
        dir=1
        opening_gap=opening_gap+1
      elif dp[cur_i][cur_j].I==dp[cur_i][cur_j-1].S +g+h:
        dir=0
        opening_gap=opening_gap+1
      else:
        dir=2

      alignment_store.append(alignment('-',' ',s2[cur_j-1]))
      cur_j=cur_j-1
      gap_extension=gap_extension+1  
      
    
    if cur_i==m and cur_j==n:
      continue

  print('',file=f_global)
  print('Matchs: ',matctes,'MisMatchs: ',mismatches,'Opening gaps: ',opening_gap, 'gap_extension: ',gap_extension,file=f_global)
  print('Identities = {0}/{1} ({2}%), Gaps: {3}/{4} ({5}%)'.format(matctes,m,(matctes*100/m),gap_extension,m,((gap_extension)*100/m)),file=f_global)
  score_op=matctes*match_score+ mismatches*mismatch_score+ opening_gap*h +gap_extension*g
  return alignment_store

# printing of global alignment

def print_output_global(alignment_result):
  print('',file=f_global)
  print('*********************',file=f_global)
  print('global alignment printing',file=f_global)
  print('*********************',file=f_global)
  print('Length of allignment store:',shape(alignment_result)[0])
  len=60
  start_index,end_index,counter_s1,counter_s2, miss_s1,miss_s2= shape(alignment_result)[0]-len, shape(alignment_result)[0], 1, 1,0,0
  counter_s1,counter_s2, miss_s1,miss_s2=  1, 1,0,0
  while end_index!=0:
    miss_s1=0
    miss_s2=0
    a=alignment_result[start_index:end_index]
    print('{0: <3}{1: <7}'.format('s1',counter_s1),end=" ",file=f_global)
    for i in reversed(a):
      print(i.s1_c, end=' ',file=f_global)
      if i.s1_c=='-':
        miss_s1=miss_s1+1
    counter_s1=(counter_s1+len)- miss_s1
    if counter_s1>m:
      counter_s1=m+1
    print(' ',counter_s1-1,file=f_global)

    print('{0: <3}{1: <7}'.format('',''),end=" ",file=f_global)
    for i in reversed(a):
      print(i.dir, end=' ',file=f_global)
    print('',file=f_global)

    print('{0: <3}{1: <7}'.format('s2',counter_s2),end=" ",file=f_global)
    for i in reversed(a):
      print(i.s2_c, end=' ',file=f_global)
      if i.s2_c=='-':
        miss_s2=miss_s2+1
    counter_s2=(counter_s2+len)- miss_s2
    if counter_s2>n:
      counter_s2=n+1

  
    print(' ',counter_s2-1,'\n',file=f_global) 
    
    end_index=start_index
    start_index=start_index-len
    if start_index<0:
      start_index=0
    

# local alignment algorithm
def local_alignment():
  print('Scroes: match= ',match_score,'  Mismatch: ',mismatch_score,' h: ',h,' g:',g,file=f_local)
  print('',file=f_local)
  print('Sequence 1 = "s1", length = ',m, 'characters',file=f_local)
  print('Sequence 2 = "s2", length = ',n, 'characters',file=f_local)
  dp_local=[]
  for i in range(m+1):
      col =[]
      for j in range(n+1):
        col.append(cell(0,0,0))
      
      dp_local.append(col)

  cur_max=-1
  i_max=-1
  j_max=-1
  for i in range(m+1):
    for j in range(n+1):
      if i==0 or j==0:
        continue
      dp_local[i][j].D=max(dp_local[i-1][j].D+g, dp_local[i-1][j].S+h+g, dp_local[i-1][j].I+h+g,0)
      dp_local[i][j].I=max(dp_local[i][j-1].D+h+g, dp_local[i][j-1].S+h+g, dp_local[i][j-1].I+g,0)
      dp_local[i][j].S=max(dp_local[i-1][j-1].S, dp_local[i-1][j-1].D, dp_local[i-1][j-1].I,0)+similarity(s1[i-1],s2[j-1])
      max_cell=max(dp_local[i][j].D,dp_local[i][j].I,dp_local[i][j].S)

      if max_cell>cur_max:
        i_max,j_max=i,j
        cur_max=max_cell

  print('',file=f_local)
  print('Local Optimal score:',max(dp_local[i_max][j_max].D,dp_local[i_max][j_max].I,dp_local[i_max][j_max].S),file=f_local)
  return dp_local,i_max,j_max

# retrace of local alignment
def retrace_local(dp_local,i_max,j_max):
  cur_i,cur_j = i_max,j_max
  dir=np.argmax([dp_local[cur_i][cur_j].S,dp_local[cur_i][cur_j].D,dp_local[cur_i][cur_j].I])
  matctes=0
  mismatches=0
  opening_gap=0
  prev=0
  gap_extension=0
  alignment_store=[]
  local_alignment_len=0
  
  while True:
    # local alignment breaking condition
    if dir==0:
      if dp_local[cur_i][cur_j].S==0:
        break
    elif dir==1:
      if dp_local[cur_i][cur_j].D==0:
        break
    elif dir==2:
      if dp_local[cur_i][cur_j].I==0:
        break

    if dir==1: # deleteion case
      if dp_local[cur_i][cur_j].D==dp_local[cur_i-1][cur_j].D+g:
        dir=1
      elif  dp_local[cur_i][cur_j].D==dp_local[cur_i-1][cur_j].S+h+g:
        dir=0
        opening_gap=opening_gap+1
      else:
        dir=2
        opening_gap=opening_gap+1
      
      alignment_store.append(alignment(s1[cur_i-1],' ','-'))
      cur_i=cur_i-1
      gap_extension=gap_extension+1  
    elif dir==0: # sub case
      if s1[cur_i-1]==s2[cur_j-1]:
        matctes=matctes+1
        pen=match_score
        alignment_store.append(alignment(s1[cur_i-1],'|',s2[cur_j-1]))
      else:
        mismatches=mismatches+1
        pen=mismatch_score
        alignment_store.append(alignment(s1[cur_i-1],' ',s2[cur_j-1]))
      if dp_local[cur_i][cur_j].S==dp_local[cur_i-1][cur_j-1].D +pen:
        dir=1
      elif dp_local[cur_i][cur_j].S==dp_local[cur_i-1][cur_j-1].S+pen:
        dir=0
      else:
        dir=2
      cur_i=cur_i-1
      cur_j=cur_j-1
    else: # insertion case
      if dp_local[cur_i][cur_j].I==dp_local[cur_i][cur_j-1].D +g+h:
        dir=1
        opening_gap=opening_gap+1
      elif dp_local[cur_i][cur_j].I==dp_local[cur_i][cur_j-1].S +g+h:
        dir=0
        opening_gap=opening_gap+1
      else:
        dir=2

      alignment_store.append(alignment('-',' ',s2[cur_j-1]))
      cur_j=cur_j-1
      gap_extension=gap_extension+1

    local_alignment_len=local_alignment_len+1  
      
    
    if cur_i==m and cur_j==n:
      continue

  print('',file=f_local)
  print('Matchs: ',matctes,'MisMatchs: ',mismatches,'Opening gaps: ',opening_gap, 'gap_extension: ',gap_extension,file=f_local)
  print('Identities = {0}/{1} ({2}%), Gaps: {3}/{4} ({5}%)'.format(matctes,local_alignment_len,(matctes*100/local_alignment_len),gap_extension,local_alignment_len,((gap_extension)*100/local_alignment_len)),file=f_local)
  score_op=matctes*match_score+ mismatches*mismatch_score+ opening_gap*h +gap_extension*g
  return alignment_store,local_alignment_len

# alignment_result_local, local_length=retrace_local(dp_local,i_max,j_max)

def print_output_local(alignment_result_local,i_max,j_max,local_length):
  print('',file=f_local)
  print('*********************',file=f_local)
  print('Local alignment printing',file=f_local)
  print('*********************',file=f_local)
  len=60
  start_index= local_length-60
  end_index= local_length
  counter_s1,counter_s2, miss_s1,miss_s2= 1, 1,i_max,j_max
  print('Local Length: ',shape(alignment_result_local))
 

  while end_index>0:
    miss_s1,miss_s2=0,0
    a=alignment_result_local[start_index:end_index]
    print('{0: <7}'.format('s1'),end=" ",file=f_local)
    for i in reversed(a):
      print(i.s1_c, end=' ',file=f_local)

    print('',file=f_local)
    print('{0: <7}'.format(' '),end=" ",file=f_local)
    for i in reversed(a):
      print(i.dir, end=' ',file=f_local)
    print('',file=f_local)
    print('{0: <7}'.format('s2'),end=" ",file=f_local)
    for i in reversed(a):
      print(i.s2_c, end=' ',file=f_local)
    print('\n',file=f_local)

        
    end_index=start_index
    if(end_index==0):
      break
    start_index=start_index-len
    if start_index<0:
      start_index=0
    

if __name__ == "__main__":
    args= list(sys.argv)
    fasta_file= args[1]
    algorithm=int(args[2])
    config_file=args[3]
    read_config(config_file) # reading config file
    load_sequences(fasta_file) # loading sequences
    if algorithm==0:
      global f_global      
      f_global= open("global_output.txt", "w") 
      dp_global=global_alignment()
      alignment_result=retrace_global(dp_global) # getting alignment result
      print_output_global(alignment_result) # printing alignment result
    elif algorithm==1: # local alignment
      global f_local
      f_local= open("local_output.txt", "w")
      dp_local,i_max,j_max=local_alignment() 
      alignment_result_local, local_length=retrace_local(dp_local,i_max,j_max) # getting alignment result
      print_output_local(alignment_result_local,i_max,j_max,local_length) # printing alignment result
    



