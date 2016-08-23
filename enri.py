

#import stuff
import os
import random
import math
from scipy import stats as scstats
from matplotlib import pyplot as plt
from tabulate import tabulate as tabulate
import matplotlib.cm as cm
import sys
import numpy as np
import scipy.stats as stats

#### data preprocessing ####

class Enri():
  def __init__(self):
    self.pydir = os.getcwd()
    self.paths = open('%s/infiles/paths.txt' % self.pydir).read().splitlines()
    self.dirs = self.make_dirs()
    self.data, self.headers = self.parse_train_txt()
    self.pocket1, self.pocket1_headers = self.get_pocket2(self.data, self.headers)


  def make_dirs(self):
    ''' make output directories  '''
    pydir = os.getcwd()
    curdir = os.listdir(pydir)
    dirs = ['figures', 'tables', 'outfiles']
    path_dirs = []
    for dir in dirs:
      path_dir = pydir  + '/%s' % dir
      path_dirs.append(path_dir)
      if dir not in curdir:
        os.system('mkdir %s' % path_dir)
    return path_dirs

  def find_files(self, filename):
    ''' find files, approximate matching '''
    file_paths = []
    for path in self.paths:
      for root, dirs, files in os.walk(path):
        for file in files:
          if filename in file:
            file_path = os.path.join(root,file)
            file_paths.append(file_path)
    return file_paths

  def find_files_from_path(self, filename, path):
    ''' find files, approximate matching '''
    file_paths = []
    for root, dirs, files in os.walk(path):
      for file in files:
        if filename in file:
          file_path = os.path.join(root,file)
          file_paths.append(file_path)
    return file_paths
 
  
  def find_files_pydir(self, filename):
    ''' find files in pydir and its subdir. approximate matching '''
    filepaths = []
    for root, dirs, files in os.walk(self.pydir): 
      for file in files:
        if filename in file:
          filepath = os.path.join(root, file)
          filepaths.append(filepath)
    print filepaths
    return filepaths


  def rename_files(self,filename):
    ''' make pdbs dir, rename pdb files, dump to pdbs dir '''
    pdbdir = '/Users/rahmadakbar/uc/enrichments/pdbs'
    os.system('mkdir %s' % pdbdir)
    files = self.find_files(filename)
    print len(files)
    for file in files:
      parts = file.split('/')
      if 'high' in  file:
        new_name = 'h_' + parts[-1]
      else:
        new_name = 'l_' +  parts[-1]
      new_file = pdbdir + '/' + new_name
      command = 'cp %s %s' % (file, new_file)
      os.system(command)
  
  def pdb2desc(self, filename):
    ''' use dogsitescorer to harvest pockets and descriptors from pdb files '''
    files = self.find_files(filename)
    print len(files)
    dogout_dir = '/'.join(files[0].split('/')[:-1]) + '/dogout'
    os.system('mkdir %s' % dogout_dir)
    exception_path = dogout_dir + '/exception.txt'
    exception_file = open(exception_path, 'w')
    for file in files[:]:
      print file
      try:
        contents = open(file).read().split('TER')
        #print contents
        ligand = contents[-1].splitlines()[1][17:20]
        protname = file.split('/')[-1].split('.')[0]
        command = 'dogsite -p %s -o %s/%s -r %s -w3 -d -i' % (file, dogout_dir, protname, ligand)
        os.system(command)
      except Exception:
        exception_file.write(file + '\n')
        #pass
    exception_file.close()


  def pdb2desc_from_path(self, pdbpath):
    ''' use dogsitescorer to harvest pockets and descriptors from pdb files '''
    files = self.find_files_from_path('pdb', pdbpath)
    print len(files)
    dogout_dir = '/'.join(files[0].split('/')[:-1]) + '/dogout'
    os.system('mkdir %s' % dogout_dir)
    exception_path = dogout_dir + '/exception.txt'
    exception_file = open(exception_path, 'w')
    for file in files[:]:
      print file
      try:
        contents = open(file).read().split('TER')
        #print contents
        ligand = contents[-1].splitlines()[1][17:20]
        protname = file.split('/')[-1].split('.')[0]
        command = 'dogsite -p %s -o %s/%s -r %s -w3 -d -i' % (file, dogout_dir, protname, ligand)
        os.system(command)
      except Exception:
        exception_file.write(file + '\n')
        #pass
    exception_file.close()

  def name2firstcol(self, filename):
    ''' add protein name to the 1st column of desc.txt files
        output .edt files '''
    files = self.find_files(filename)
    for file in files[:]:
      data = ''
      names = self.path2names(file)
      protname = names[1][:-4]
      contents = open(file).readlines()
      for line in contents[:]:
        if line.startswith('name'):
          data += line
        else:
          newline = protname + line
          data += newline
      outpath = names[0] + '/' + names[1] + '.edt'
      outfile = open(outpath, 'w') 
      outfile.write(data)


  def name2firstcol_from_path(self, pdbpath):
    ''' add protein name to the 1st column of desc.txt files
        output .edt files '''
    files = self.find_files_from_path('desc.txt', pdbpath)
    for file in files[:]:
      data = ''
      names = self.path2names(file)
      protname = names[1][:-4]
      contents = open(file).readlines()
      for line in contents[:]:
        if line.startswith('name'):
          data += line
        else:
          newline = protname + line
          data += newline
      outpath = names[0] + '/' + names[1] + '.edt'
      outfile = open(outpath, 'w') 
      outfile.write(data)

  def merge_edt(self, filename='.edt'):
    ''' merge all .edt files to a single file '''
    files = self.find_files(filename)
    names = self.path2names(files[0])
    path = names[0]
    merged_content = open(files[0]).read()
    for file in files[1:]:
      contents = open(file).readlines()[1:] # skip header
      for line in contents:
        merged_content += line
    outfile = open('%s/desc_merged.txt' % path, 'w')
    outfile.write(merged_content)
    outfile.close()
    print len(open('%s/desc_merged.txt' % path).readlines())


  def merge_edt_from_path(self, pdbpath):
    ''' merge all .edt files to a single file '''
    files = self.find_files_from_path('.edt', pdbpath)
    names = self.path2names(files[0])
    path = names[0]
    merged_content = open(files[0]).read()
    for file in files[1:]:
      contents = open(file).readlines()[1:] # skip header
      for line in contents:
        merged_content += line
    outfile = open('%s/desc_merged.txt' % path, 'w')
    outfile.write(merged_content)
    outfile.close()
    print len(open('%s/desc_merged.txt' % path).readlines())

    
  def parse_edt(self, filename = 'desc_merged.txt'):
    ''' parse desc_merged.txt, output features and headers'''
    file = open(self.find_files(filename)[0]).read().splitlines()
    unwanted_index = [3]
    data = [line.split('\t') for line in file]
    data2 = [[x for i,x in enumerate(item) if i not in unwanted_index] for item in data[:]]
    headers = data2[0][:]
    data2 = [[item[0]]+[float(x) for x in item[1:]] for item in data2[1:]]
    print len(data2)
    print headers
    print data2[0]
    return data2, headers


  def parse_desc_merged_txt(self, filepath):
    ''' parse *desc.txt file given absolute path '''
    file = open(filepath).read().splitlines()
    unwanted_index = [3,63]
    unwanted_headers = ['lig_name', 'UNK']
    data = [line.split('\t') for line in file]
    data2 = [[x for i,x in enumerate(item) if i not in unwanted_index] for item in data[:]]
    headers = data[0][:]
    contents = data[1:]
    float_contents = []
    float_headers = [header for header in headers if header not in unwanted_headers]
    for content in contents:
      protname = content[0]
      descs = content[1:]
      float_content = [protname]
      for i in range(1, len(headers)):
        header = headers[i]
        if header not in unwanted_headers:
          value = float(content[i])
          float_content.append(value)
      float_contents.append(float_content)
    return float_contents, float_headers
    #data2 = [[item[0]]+[float(x) for x in item[1:]] for item in data2[1:]]
    #return data2, headers

  def parse_train_txt(self):
    ''' parse tain.txt file, outputs features and headers'''
    #train_file = self.pydir + '/outfiles/hdstat_data.tsv'
    #train_file = self.pydir + '/infiles/train2.txt'
    train_file = self.pydir + '/infiles/train3.txt'
    file = open(train_file).read().splitlines()
    unwanted_index = [3] #previously was [3]
    data = [line.split('\t') for line in file]
    data2 = [[x for i,x in enumerate(item) if i not in unwanted_index] for item in data[:]]
    headers = data2[0][:]
    data2 = [[item[0]]+[float(x) for x in item[1:]] for item in data2[1:]]
    return data2, headers

  def parse_train_txt2(self, filepath):
    ''' parse tain.txt file from an absoulute path, outputs features and headers'''
    train_file = filepath
    file = open(train_file).read().splitlines()
    unwanted_index = [3]
    data = [line.split('\t') for line in file]
    data2 = [[x for i,x in enumerate(item) if i not in unwanted_index] for item in data[:]]
    headers = data2[0][:]
    data2 = [[item[0]]+[float(x) for x in item[1:]] for item in data2[1:]]
    return data2, headers


  def parse_train_txts(self):
    ''' parse train.txt files. returns a list of features and headers  '''
    train_files = []
    traintxt_path = self.pydir + '/infiles'
    for root, dirs, files in  os.walk(traintxt_path):
      for file in files:
        if  'train' in file:
          train_file = os.path.join(root, file)
          train_files.append(train_file)
    train_datas = []
    for train_file in train_files:
      names = self.path2names(train_file)
      name = names[1]
      train_data = self.parse_train_txt2(train_file)
      train_data.append(name)
      train_datas.append(train_data)
    return train_datas

  def plot_train_datas(self):
    ''' plots histogram of the features from each train data '''
    train_datas = self.parse_train_txts()
    for train_data in train_datas:
      data, headers, name = train_data
      clean_data, clean_headers = self.clean_data(data, headers)
      self.plot_hist(clean_data, clean_headers, name)

  def get_pocket(self, pocket, filename = 'desc_merged.txt'):
    ''' returns a specific pocket from the file. for instance, setting pocket =
      P_1 returns only pocket P_1   '''
    data, headers = self.parse_edt(filename)
    pockets = []
    for item in data:
      name = item[0]
      if pocket == name[-3:]:
        pockets.append(item)
    return pockets, headers

  def get_pocket2(self, data, headers, pocket = 'P_1'):
    ''' returns a specific pocket from the file. for instance, setting pocket =
      P_1 returns only pocket P_1   '''
    pockets = []
    for item in data:
      name = item[0]
      if pocket == name[-3:]:
        pockets.append(item)
    return pockets, headers
 
  
  def split_con_nom(self, filename = 'desc_merged.txt'):
    data, headers= self.get_pocket('P_1', filename)
    con_index, nom_index = [0], [0]
    for i,item in enumerate(data[100][1:]):
      rem = item%1
      if rem != 0:
        con_index.append(i+1)
      else:
        nom_index.append(i+1)
    con_data, nom_data = [], []
    for item in data:
      name = item[0]
      con_vals, nom_vals = [], []
      con_vals.append(name)
      nom_vals.append(name)
      for i,val in enumerate(item[1:]):
        index = i+1
        if index in con_index:
          con_vals.append(val)
        else:
          nom_vals.append(val)
      con_data.append(con_vals)
      nom_data.append(nom_vals)
    con_headers, nom_headers = [headers[0]], [headers[0]]
    for i in range(len(headers[1:])):
      index = i+1
      if index in con_index:
        con_headers.append(headers[index])
      else:
        nom_headers.append(headers[index])
    return [con_data, con_headers, 'con_data'],[nom_data, nom_headers, 'nom_data']


  def output_data_files(self, filename):
    ''' output the splitted continous and nominal data onto two filse
      continous.txt and nominal.txt '''
    datas = self.split_con_nom(filename)
    for data in  datas:
      content = ''
      name = data[2]
      outname = self.dirs[2] + '/%s.txt' %name
      outfile = open(outname, 'w')
      content   +=  ','.join(data[1]) + '\n'
      for item in data[0]:
        content += ','.join(str(x) for x in item) + '\n'
      outfile.write(content)
      outfile.close()
     
  def split_l_h(self,list):
    ''' split high (h) and low (l). the list format should follow that of
      parse_edt output'''
    l = []
    h = []
    for item in list:
      name = item[0] 
      if name[0] == 'h':
        h.append(item)
      else:
        l.append(item)
    return l, h

  def plot_con_hist(self, filename='desc_merged.txt'):
    data = self.split_con_nom(filename)
    cons =  data[0]
    con_data = cons[0]
    con_headers = cons[1]
    l,h = self.split_l_h(con_data)
    for i in range(1,len(con_headers)):
      l_xs = [x[i] for x in l]
      h_xs = [x[i] for x in h]
      header = con_headers[i].replace('/', '_')
      plot_name = self.dirs[0] + '/%s_hist.pdf' % header
      plt.figure()
      plt.hist(l_xs, bins = 20, alpha = 0.5, label = 'l', edgecolor = 'none')
      plt.hist(h_xs, bins = 20, alpha = 0.5, label = 'h', edgecolor = 'none')
      plt.legend()
      plt.savefig(plot_name)
    print 'figures are saved to %s' % self.dirs[0]


  def plot_nom_hist(self, filename='desc_merged.txt'):
    data = self.split_con_nom(filename)
    cons =  data[0]
    nom_data = cons[0]
    nom_headers = cons[1]
    l,h = self.split_l_h(nom_data)
    for i in range(1,len(nom_headers)):
      l_xs = [x[i] for x in l]
      h_xs = [x[i] for x in h]
      header = nom_headers[i].replace('/', '_')
      plot_name = self.dirs[0] + '/%s_hist.pdf' % header
      plt.figure()
      plt.hist(l_xs, bins = 20, alpha = 0.5, label = 'l', edgecolor = 'none')
      plt.hist(h_xs, bins = 20, alpha = 0.5, label = 'h', edgecolor = 'none')
      plt.legend()
      plt.savefig(plot_name)
    print 'figures are saved to %s' % self.dirs[0]

  
  def plot_hist(self, data, headers, figname):
    ''' outputs histogram of each features '''
    l,h = self.split_l_h(data)
    for i in range(1,len(headers)):
      l_xs = [x[i] for x in l]
      h_xs = [x[i] for x in h]
      header = headers[i].replace('/', '_')
      plotpath = self.dirs[0] + '/%s_%s_hist.pdf' % (header,figname)
      plt.figure()
      plt.hist(l_xs, bins = 20, alpha = 0.5, label = 'l', edgecolor = 'none')
      plt.hist(h_xs, bins = 20, alpha = 0.5, label = 'h', edgecolor = 'none')
      plt.legend()
      plt.savefig(plotpath)
      plt.close()
      #break
    print 'figures are saved to %s' % self.dirs[0]


  def plot_hist_nolabel(self, data, headers, outdir):
    ''' outputs histogram of each features '''
    for i in range(1,len(headers)):
      xs = [x[i] for x in data]
      mean, sd = round(self.mean(xs),3), round(self.sd(xs),3)
      plotlabel = 'mean %s, sd %s' % (mean, sd)
      header = headers[i].replace('/', '_')
      plotpath = outdir + '/%s_hist.pdf' % header
      plt.figure()
      plt.hist(xs, bins = 50, alpha = 0.5, edgecolor = 'none', label=plotlabel)
      plt.legend(frameon = False)
      plt.title(headers[i])
      plt.savefig(plotpath)
      plt.close()
      #break
    print 'figures are saved to %s' % outdir


  def rows2columns(self, data, headers):
    ''' transform data from parse_edt rows list to a columns list '''
    coldata = [[] for x in headers]
    for item in data:
      for i, val in enumerate(item):
        coldata[i].append(val)  
    return coldata, headers

  def columns2rows(self, coldata, headers):
    ''' transforms columns to rows '''
    rows = len(coldata[0])
    rowdata = [[] for x in range(rows)]
    for  items in coldata:
      for i,item in enumerate(items):
        rowdata[i].append(item)
    return rowdata, headers


  def get_nonzero_features(self, data, headers):
    ''' get non zero features from data, headers list. drop UNK and lig_cov'''
    data, headers = self.rows2columns(data, headers)
    new_headers = []
    new_data = []
    unwanted_vars = ['UNK', 'lig_cov']
    for i, item in enumerate(data[:]):
      var_name = headers[i]
      if i == 0:
        new_data.append(item)
        new_headers.append(headers[i])
      elif i > 0 and sum(item) > 0 and var_name not in unwanted_vars:
        new_data.append(item) 
        new_headers.append(headers[i])
    rowdata, new_headers = self.columns2rows(new_data, new_headers)
    return rowdata, new_headers

  def clean_data(self, data, headers):
    ''' returns clean data, grabs only pocket P_1 and removes zero features '''
    pockets, headers = self.get_pocket2(data, headers)
    nz_data, nz_headers = self.get_nonzero_features(pockets, headers)
    return nz_data, nz_headers    

  def clean_data_continuous(self, data, headers, sample_index):
    ''' returns clean data, grabs only continuous data for pocket P_1 and removes zero features '''
    pockets, headers = self.get_pocket2(data, headers)
    nz_data, nz_headers = self.get_nonzero_features(pockets, headers)
    con_data, con_headers = self.continuous_data(nz_data, nz_headers, sample_index)
    #print con_headers
    return con_data, con_headers    


  def path2names(self, filepath):
    ''' returns file name with and without extention for a given path '''
    parts = filepath.split('/')
    path = '/'.join(parts[:-1])
    name_ex = parts[-1]
    name = name_ex.split('.')[0]
    return path, name, name_ex

   

### adaptive synthetic sampling ###
  def euclidian_distance(self, p, q):
    ''' get euclidian distance from two points. points are a vector of n
    dimensions'''
    d = math.sqrt(sum([(p_i - q_i)**2 for p_i, q_i in zip(p, q) ]))
    return d
  
  def vector_diff(self, p,q):
    ''' get the difference between two vectors '''
    dif = [p_i - q_i for p_i, q_i in zip(p,q)]
    return  dif
  
  def vector_sum(self, p,q):
    ''' get the sum between two vectors '''
    sum = [p_i + q_i for p_i, q_i in zip(p,q)]
    return sum
  
  def scalar_mul(self, p, s):
    ''' returns scalar multiplication of vector p and scalar s '''
    mul = [p_i*s for p_i in p]
    return mul

  def scalar_mul_random(self, p):
    ''' returns scalar multiplication of vector p and scalar s '''
    s = random.random()
    mul = [p_i*s for p_i in p]
    return mul

  def generate_synthetic_data(self, x, p):
    ''' generate a synthetic data x_new. p is a randomly choosen nearest neighbor of x'''
    diff = self.vector_diff(p,x)
    random_diff = self.scalar_mul_random(diff)
    x_new = self.vector_sum(x, random_diff)
    return x_new


  def get_knn(self, pockets, headers,  k = 10):
    ''' returns and write-out  a dictionary containing k nearest neighbors (knn) '''
    knns = {}
    for pocket in pockets:
      neighbors = []
      pocket_name = pocket[0]
      for pocket2 in pockets:
        if pocket != pocket2: 
          distance = self.euclidian_distance(pocket[1:], pocket2[1:])
          pocket2_name = pocket2[0] 
          neighbors.append((pocket2_name, distance))
      knn = sorted(neighbors, key = lambda neighbor: neighbor[1])[:k] 
      knns[pocket_name] = knn
    return knns, pockets, headers


  def get_minority(self, knns):
    ''' returns  the minority class and its neighbors '''
    mins =  {}
    label = 'h_'
    for item in knns:
      if item[:2] == label:
        mins[item] = knns[item]
    return mins

  def gamma(self, mins):
    ''' returns a dictionary containing gamma ratio for each data point in the
      minority class mins. delta_i is number of neigbors belonging to the
      majority class, k total number of neigbors. z is normalizing factor so
      that sum of all gammas are 1'''
    k = len(mins.items()[0][1])
    gammas = {}
    for item in mins:
      delta_i = 0.0
      neighbors = mins[item]
      for neighbor in neighbors:
        name = neighbor[0]
        if name[:2] != item[:2]:
          delta_i = delta_i + 1
      gamma = delta_i/k
      gammas[item] = gamma
    norm_gammas = {}
    z = sum(gammas.values())
    for item in gammas:
      gamma = gammas[item]
      norm_gamma = gamma/z
      norm_gammas[item] = norm_gamma
    return norm_gammas

  def synthetic_data_number(self, knns, beta = 0.65):
    ''' returns the number of synteti data needed to balance the classes for
    each minority class member. g is the difference between the two classes '''
    mins = self.get_minority(knns)
    norm_gammas = self.gamma(mins)
    g = (len(knns)- len(mins)) * beta
    gs = {}
    for item in mins:
      gamma_i = norm_gammas[item]
      g_i = gamma_i * g
      g_i =  int(round(g_i))
      gs[item] = g_i
    return gs
  

  def adasyn(self, data, headers, k = 10,  beta = 0.65):
    ''' genearate new samples for each member of the minority class by using
    its gamma distribution (norm_gammas) to decide the number of needed
    synthetic sample(gs). this is known as adapative synthetic sampling. now
    use only neighbor with class h to gernerate the synthetic samples '''
    knns, pockets, headers = self.get_knn(data, headers, k)
    gs = self.synthetic_data_number(knns, beta)
    pocket_dict = {}
    for pocket in pockets:
      name = pocket[0]
      feats = pocket[1:]
      pocket_dict[name] = feats
    syn_pockets = {}
    for item in gs:
      g = gs[item]
      name_parts = item.split('_')
      for i in range(g):
        name = [name_parts[0]] + ['syn%s' % str(i)] + name_parts[1:]
        name = '_'.join(name)
        neighbors = knns[item]
        neighbors2 = [neighbor for neighbor in neighbors if neighbor[0][0] == name[0]]
        if len(neighbors2) > 1:
          random_neighbor = random.choice(neighbors2)[0]
        else:
          random_neighbor = random.choice(neighbors)[0]
        x = pocket_dict[item] 
        p = pocket_dict[random_neighbor]
        new_x = self.generate_synthetic_data(x,p)
        syn_pockets[name] = new_x
    merged_pockets = pocket_dict.copy()
    merged_pockets.update(syn_pockets)
#    outname = self.dirs[2] + '/adasyn_k%s_beta%s.txt' % (str(k), str(beta))
#    outfile = open(outname, 'w')
#    outfile.write(','.join(headers) + '\n')
#    for item in merged_pockets:
#      content = [item] + [str(x) for x in merged_pockets[item]]
#      content = ','.join(content) + '\n'
#      outfile.write(content)
#    outfile.close()
    adasyn_data = [[item] + merged_pockets[item] for item in merged_pockets]
    return adasyn_data, headers


  def plot_adasyn(self, k, beta): 
    ''' plot adaptive synthetic samples for train data '''
    adasyn_data ,adasyn_headers = self.adasyn(self.pocket1, self.pocket1_headers, k, beta)
    print self.pocket1
    sys.exit()
    nz_data, nz_headers = self.get_nonzero_features(adasyn_data, adasyn_headers)
    plotname = 'adasyn_k%s_beta%s' % (k, beta)
    self.plot_hist(nz_data, nz_headers, plotname)
#### naive bayes #####

#  def continuous_index(self, filename):
#    ''' returns continius variables indices and thei headers '''
#    data, headers = self.get_pocket('p_1', filename)
#    data, headers = self.get_nonzero_features2(data, headers)
#    con_index = [0] 
#    con_headers =[headers[0]] 
#    sample = data[120] 
#    for i in range(1,len(sample)):
#      item = sample[i]
#      rem = item%1
#      if rem != 0:
#        con_index.append(i)
#        con_headers.append(headers[i])
#    return con_index, con_headers

  def find_sample_index(self):
    ''' find sample index to be used to fetch continuous data '''
    train_datas= self.parse_train_txts()
    con_counts_collections = [] 
    for train_data in train_datas:
      data, headers, name = train_data
      data, headers = self.clean_data(data, headers)
      con_counts = []
      for i, items in enumerate(data):
        con_count = 0
        for item in items[1:]:
          rem = item%1
          if rem != 0:
            con_count += 1
        con_counts.append([con_count, i])
      con_counts_collections.append(con_counts)
    max_indices = []
    for item in con_counts_collections:
      con_counts = sorted(item)
      max = con_counts[-1][0]
      max_index = [x[1] for x in con_counts if x[0] == max]
      max_indices.append(max_index)
    intersect = set(max_indices[0]).intersection(max_indices[1])
    print list(intersect)[0]
    return list(intersect)[0]


  def continuous_data(self, data, headers, sample_index):
    ''' returns continuous data. use indices from continuos_index '''
    con_index = [0] 
    con_headers =[headers[0]] 
    sample = data[sample_index] 
    for i in range(1,len(sample)):
      item = sample[i]
      rem = item%1
      if rem != 0:
        con_index.append(i)
        con_headers.append(headers[i])
    con_data = []
    for item in data:
      con = [item[i] for i in con_index]
      con_data.append(con)
    return con_data, con_headers

#  def nominal_index(self, filename):
#    ''' returns nominal indices and their headers '''
#    data, headers = self.get_pocket('p_1', filename)
#    data, headers = self.get_nonzero_features2(data, headers)
#    nom_index = [0] 
#    nom_headers =[headers[0]] 
#    sample = data[120] 
#    for i in range(1,len(sample)):
#      item = sample[i]
#      rem = item%1
#      if rem == 0:
#        nom_index.append(i)
#        nom_headers.append(headers[i])
#    return nom_index, nom_headers

  def nominal_data(self, data, headers, sampel_index):
    ''' returns nominal data '''
    nom_index = [0] 
    nom_headers =[headers[0]] 
    sample = data[sampel_index] 
    for i in range(1,len(sample)):
      item = sample[i]
      rem = item%1
      if rem == 0:
        nom_index.append(i)
        nom_headers.append(headers[i])
    nom_data = []
    for item in data:
      nom = [item[i] for i in nom_index]
      nom_data.append(nom)
    return nom_data, nom_headers



  def mean(self, list):
    ''' compute sample mean '''
    mean = sum(list)/float(len(list))
    return mean

  def demean(self, list, mean):
    ''' substarct mean from teh samples '''
    demean = [x_i - mean for x_i in list]
    return demean
  
  def sum_squared(self, list):
    ''' return sum of squred of a collection '''
    ss = sum(x_i**2 for x_i in list)
    return ss

  def sd(self, list):
    ''' return standard deviation of a colection '''
    mu = self.mean(list)
    n = float(len(list))
    demean = self.demean(list, mu) 
    ss = self.sum_squared(demean)
    sd = math.sqrt(ss/n)
    return sd

  def p_normal_cdf(self, x, mu, sigma): 
    ''' return cumulative probabilty from a cumulative distribution function. assumes the
    distribution is normal '''
    half_sigma = sigma/2.0
    x_up = x + sigma
    x_down = x - sigma
    p_up  = (1 + math.erf((x_up-mu)/sigma/math.sqrt(2)))/float(2)
    p_down  = (1 + math.erf((x_down-mu)/sigma/math.sqrt(2)))/float(2)
    p = p_up - p_down
#    p  = (1 + math.erf((x-mu)/sigma/math.sqrt(2)))/float(2)
    return p

  def get_params(self, list):
    ''' returns mu and sigma for a collection '''
    mu = self.mean(list)
    sd = self.sd(list)
    if mu == 0:
      mu = 1e-14
    if sd == 0:
      sd = 1e-14
    return mu, sd

  def get_params_dict(self, data, headers):
    ''' retuns parameters dictionary for a given data '''
    data, headers = self.rows2columns(data, headers)
    params  = {}
    params[headers[0]] = data[0]
    for i in range(1, len(headers)):
      header = headers[i]
      values = data[i]
      mu, sigma = self.get_params(values)
      params[header] = mu, sigma
    return params

  def get_kde_dict(self, data, headers):
    ''' get a gaussian kernel density estimator for each feature in the data '''
    data, headers = self.rows2columns(data, headers)
    kdes = {}
    kdes[headers[0]] = data[0]
    for i in range(1, len(headers)):
      header = headers[i]
      values = data[i]
      kde = scstats.kde.gaussian_kde(values)
      minval = min(values)
      kdes[header] = kde,minval
    return kdes

  def split_train_test(self, data, fraction = 0.1):
    '''  split data to train and test set according to the given percentage'''
    n = int(round(fraction*len(data)))
    random.shuffle(data)
    test = data[:n]
    train = data[n:]
    return train, test
  
  def join_p(self, list):
    ''' return a joint probability from a collection '''
    join = 1.0
    for item in list:
      join = join*item
    return join
  
  def random_subset_validation(self, data, headers, fraction = 0.5):
    ''' random cross validate the data. by default, splits 10% for test, 90% for training '''
    l, h = self.split_l_h(data)
    train_l, test_l = self.split_train_test(l, fraction)
    train_h, test_h = self.split_train_test(h, fraction)
    test = test_l + test_h
    train = train_l + train_h
    params_l = self.get_params_dict(train_l, headers)
    params_h = self.get_params_dict(train_h, headers)
    prior_l, prior_h = float(len(train_l))/len(train), float(len(train_h))/len(train)
    predictions = []
    for items in test:
      probs_l = [prior_l]
      probs_h = [prior_h]
      sample_name = items[0]
      for i in range(1, len(items)):
        var_name = headers[i]
        mu_l, sigma_l = params_l[var_name]
        mu_h, sigma_h = params_h[var_name]
        x = items[i]
        prob_l = self.p_normal_cdf(x, mu_l, sigma_l)
        prob_h = self.p_normal_cdf(x, mu_h, sigma_h)
        probs_l.append(prob_l)
        probs_h.append(prob_h)
      join_l = self.join_p(probs_l)
      join_h = self.join_p(probs_h)
      p_ratio = join_h/join_l
      sample_class = sample_name[0]
      if p_ratio > 1:
        predicted_class = 'h'
      else:
        predicted_class = 'l'
      prediction = [sample_name, sample_class, predicted_class, p_ratio]
      predictions.append(prediction)
    return predictions

  def confusion_matrix(self, predictions):
    ''' returns confusion matrix for the given predictions. input follows that
    of random_subset_validation format '''
    classes = self.split_l_h(predictions)
    c_matrix = [[] for x in range(len(classes))] 
    for i, items in enumerate(classes):
      label = items[0][0][0]
      c_matrix[i].append(label)
      match = 0
      notmatch = 0 
      for item in items:
        predicted_label = item[2]
        if label == predicted_label:
          match += 1
        else:
          notmatch += 1
      c_matrix[i].append(match)
      c_matrix[i].append(notmatch)
    #print c_matrix
    return c_matrix 

  def kfold_split(self, data, k = 10):
    ''' split data to kfold '''
    fraction = 1.0/k
    portion = int(fraction * len(data))
    splits = []
    for i in range(0, len(data), portion):
      split = data[i:i+portion]
      splits.append(split)
    return splits

  
  def kfold_cross_validation(self, data, headers, k = 10):
    ''' do k fold validation. default is 10 '''
    l, h = self.split_l_h(data)
    splits_l = self.kfold_split(l, k)
    splits_h = self.kfold_split(h, k)
    c_matrices = []
    for i in range(len(splits_l)):
      test_l, test_h = splits_l[i],splits_h[i]
      train_l = []
      train_h = []
      for i2 in range(len(splits_l)):
        if i2 != i:
          train_l += splits_l[i2]
          train_h += splits_h[i2]
      train = train_l + train_h
      test = test_l + test_h
      test = test_l + test_h
      train = train_l + train_h
      params_l = self.get_params_dict(train_l, headers)
      params_h = self.get_params_dict(train_h, headers)
      prior_l, prior_h = float(len(train_l))/len(train), float(len(train_h))/len(train)
      predictions = []
      for items in test:
        probs_l = [prior_l]
        probs_h = [prior_h]
        sample_name = items[0]
        for i in range(1, len(items)):
          var_name = headers[i]
          mu_l, sigma_l = params_l[var_name]
          mu_h, sigma_h = params_h[var_name]
          x = items[i]
          prob_l = self.p_normal_cdf(x, mu_l, sigma_l)
          prob_h = self.p_normal_cdf(x, mu_h, sigma_h)
          probs_l.append(prob_l)
          probs_h.append(prob_h)
        join_l = self.join_p(probs_l)
        join_h = self.join_p(probs_h)
        p_ratio = join_h/join_l
        sample_class = sample_name[0]
        if p_ratio > 1:
          predicted_class = 'h'
        else:
          predicted_class = 'l'
        prediction = [sample_name, sample_class, predicted_class, p_ratio]
        predictions.append(prediction)
      c_matrix = self.confusion_matrix(predictions)
      c_matrices.append(c_matrix)
    ave_c_matrix = [[] for x in range(len(c_matrices[0]))]
    for c_matrix in c_matrices:
      for i in range(len(ave_c_matrix)):
        if len(ave_c_matrix[i]) == 0:
          ave_c_matrix[i].append(c_matrix[i][0])
          ave_c_matrix[i].append(c_matrix[i][1])
          ave_c_matrix[i].append(c_matrix[i][2])
        else:
          ave_c_matrix[i][1] += c_matrix[i][1]
          ave_c_matrix[i][2] += c_matrix[i][2]
    for item in ave_c_matrix:
      for i in range(1, len(ave_c_matrix)+1):
        item[i] = item[i]/float(k)
    return ave_c_matrix


  def kfold_cross_validation2(self, data, headers, k = 10):
    ''' do k fold validation. default is 10. now understands discrpancy between
    the k splits, some times they are different due to rounding discrepancies '''
    l, h = self.split_l_h(data)
    splits_l = self.kfold_split(l, k)
    splits_h = self.kfold_split(h, k)
    if len(splits_h) < len(splits_l):
      splits_l = splits_l[:-1]
    elif len(splits_h) > len(splits_l):
      splits_h = splits_h[:-1]
    c_matrices = []
    for i in range(len(splits_l)):
      test_l, test_h = splits_l[i],splits_h[i]
      train_l = []
      train_h = []
      for i2 in range(len(splits_l)):
        if i2 != i:
          train_l += splits_l[i2]
          train_h += splits_h[i2]
      train = train_l + train_h
      test = test_l + test_h
      #test = test_l + test_h
      #train = train_l + train_h
      params_l = self.get_params_dict(train_l, headers)
      params_h = self.get_params_dict(train_h, headers)
      prior_l, prior_h = float(len(train_l))/len(train), float(len(train_h))/len(train)
      predictions = []
      for items in test:
        probs_l = [prior_l]
        probs_h = [prior_h]
        sample_name = items[0]
        for i in range(1, len(items)):
          var_name = headers[i]
          mu_l, sigma_l = params_l[var_name]
          mu_h, sigma_h = params_h[var_name]
          #print var_name, mu_l, sigma_l, mu_h, sigma_h
          x = items[i]
          prob_l = self.p_normal_cdf(x, mu_l, sigma_l)
          prob_h = self.p_normal_cdf(x, mu_h, sigma_h)
          probs_l.append(prob_l)
          probs_h.append(prob_h)
        join_l = self.join_p(probs_l)
        join_h = self.join_p(probs_h)
        p_ratio = join_h/join_l
        sample_class = sample_name[0]
        if p_ratio > 1:
          predicted_class = 'h'
        else:
          predicted_class = 'l'
        prediction = [sample_name, sample_class, predicted_class, p_ratio]
        predictions.append(prediction)
      c_matrix = self.confusion_matrix(predictions)
      c_matrices.append(c_matrix)
    ave_c_matrix = [[] for x in range(len(c_matrices[0]))]
    for c_matrix in c_matrices:
      for i in range(len(ave_c_matrix)):
        if len(ave_c_matrix[i]) == 0:
          ave_c_matrix[i].append(c_matrix[i][0])
          ave_c_matrix[i].append(c_matrix[i][1])
          ave_c_matrix[i].append(c_matrix[i][2])
        else:
          ave_c_matrix[i][1] += c_matrix[i][1]
          ave_c_matrix[i][2] += c_matrix[i][2]
    for item in ave_c_matrix:
      for i in range(1, len(ave_c_matrix)+1):
        item[i] = item[i]/float(k)
    return ave_c_matrix


  def kfold_cross_validation3(self, data, headers, k = 10):
    ''' do k fold validation. default is 10. Uses gaussian kde to estimate
    densitiesnow understands discrpancy between
    the k splits, some times they are different due to rounding discrepancies '''
    l, h = self.split_l_h(data)
    splits_l = self.kfold_split(l, k)
    splits_h = self.kfold_split(h, k)
    if len(splits_h) < len(splits_l):
      splits_l = splits_l[:-1]
    elif len(splits_h) > len(splits_l):
      splits_h = splits_h[:-1]
    c_matrices = []
    for i in range(len(splits_l)):
      test_l, test_h = splits_l[i],splits_h[i]
      train_l = []
      train_h = []
      for i2 in range(len(splits_l)):
        if i2 != i:
          train_l += splits_l[i2]
          train_h += splits_h[i2]
      train = train_l + train_h
      test = test_l + test_h
      #test = test_l + test_h
      #train = train_l + train_h
      kdes_l = self.get_kde_dict(train_l, headers)
      kdes_h = self.get_kde_dict(train_h, headers)
#      prior_l, prior_h = float(len(train_l))/len(train), float(len(train_h))/len(train)
      prior_l, prior_h = 0.5, 0.5
      predictions = []
      for items in test:
        probs_l = [prior_l]
        probs_h = [prior_h]
        sample_name = items[0]
        for i in range(1, len(items)):
          var_name = headers[i]
          #mu_l, sigma_l = params_l[var_name]
          #mu_h, sigma_h = params_h[var_name]
          #print var_name, mu_l, sigma_l, mu_h, sigma_h
          kde_l, min_l = kdes_l[var_name]
          kde_h, min_h = kdes_h[var_name]
          x = items[i]
          prob_l = kde_l.evaluate(x)
          prob_h = kde_h.evaluate(x)
#          prob_l = kde_l.integrate_box_1d(min_l, x)
#          prob_h = kde_h.integrate_box_1d(min_h, x)
          probs_l.append(prob_l)
          probs_h.append(prob_h)
        join_l = self.join_p(probs_l)
        join_h = self.join_p(probs_h)
        p_ratio = join_h/join_l
        sample_class = sample_name[0]
        if p_ratio > 1:
          predicted_class = 'h'
        else:
          predicted_class = 'l'
        prediction = [sample_name, sample_class, predicted_class, p_ratio]
        predictions.append(prediction)
      c_matrix = self.confusion_matrix(predictions)
      c_matrices.append(c_matrix)
    ave_c_matrix = [[] for x in range(len(c_matrices[0]))]
    for c_matrix in c_matrices:
      for i in range(len(ave_c_matrix)):
        if len(ave_c_matrix[i]) == 0:
          ave_c_matrix[i].append(c_matrix[i][0])
          ave_c_matrix[i].append(c_matrix[i][1])
          ave_c_matrix[i].append(c_matrix[i][2])
        else:
          ave_c_matrix[i][1] += c_matrix[i][1]
          ave_c_matrix[i][2] += c_matrix[i][2]
    for item in ave_c_matrix:
      for i in range(1, len(ave_c_matrix)+1):
        item[i] = item[i]/float(k)
    return ave_c_matrix



  def tp_fp(self, c_matrix):
    ''' returns true positive false positive rates from c_matrix. assume true=ih'''
    rates = []
    for item in c_matrix:
      rate = float(item[1])/sum(item[1:])
      rates.append(rate)
    return rates

  def tp_fp2(self, c_matrix):
    ''' returns true positive false positive rates from c_matrix. assume true=ih'''
    rates = []
    tp = c_matrix[1][1]
    fn = c_matrix[1][2]
    tn = c_matrix[0][1]
    fp = c_matrix[0][2] 
    tpr = round(tp/(tp+fn),3)
    fpr = round(fp/(fp+tn),3)
    rates.append(fpr)
    rates.append(tpr)
    return rates


  def roc_beta(self, data, headers, beta_low = 0.2, beta_high = 50, skip = 0.2):
    ''' returns tp fp rates of betas range from beta_low to beta_h '''
    n = int((beta_high-beta_low)/skip)
    params = [round(beta_low + i*skip, 2) for i in range(n)]
    rates = []
    for param in params:
      adasyn_data, adasyn_headers = self.adasyn(data, headers, k = 10 , beta = param)
      c_matrix = self.kfold_cross_validation2(adasyn_data, headers)
      rate = self.tp_fp2(c_matrix)
      rate.append(param)
      rates.append(rate)
    return rates 

  def write_roc_beta(self, rates, filename = 'roc_betas'):
    contents = ''
    rates = sorted(rates)
    for rate in rates:
      content = ','.join(str(round(x, 4)) for x in rate)  + '\n'
      contents += content
    outname = self.dirs[2] + '/%s.txt' % filename
    outfile = open(outname, 'w')
    outfile.write(contents)
    outfile.close()

  def plot_roc(self, rates, figname):
    ''' plot tpr as a function of fpr. this is known as roc curve'''
    fprs = [x[0] for x in rates] 
    tprs = [x[1] for x in rates]
    figpath = self.dirs[0]
    figname = figpath + '/%s.pdf' % figname
    plt.figure()
    plt.scatter(fprs, tprs, alpha = 0.5, edgecolor = 'none', color = 'blue')
    diags = [0.1*x for x in range(11)]
    plt.plot(diags, diags, alpha = 0.5, color = 'orange')
    plt.axis([0, 1, 0,1])
    plt.xlabel('False positive rate')
    plt.ylabel('True positive rate')
    plt.savefig(figname) 
    plt.close()  

  def plot_train_rocs(self):
    ''' plots receiver operating characteristics (roc) for each train data in
    infiles. writes rates to outfiles as well'''
    train_datas = self.parse_train_txts()
    for train_data in train_datas:
      data, headers, name = train_data
      clean_data, clean_headers = self.clean_data_continuous(data, headers,250)
      beta_low, beta_high = 0.2, 50
      rates = self.roc_beta(clean_data, clean_headers, beta_low, beta_high)
      filename = 'roc_beta%s_%s' % (beta_high, name)
      self.write_roc_beta(rates, filename)
      self.plot_roc(rates, filename)

  def parse_rates_file(self, rates_file):
    ''' parse ratesfile deposited in outfiles dir '''
    rates = open(rates_file).read().splitlines()
    rates = [x.split(',') for x in rates]
    fprs = [x[0] for x in rates]
    tprs = [x[1] for x in rates]
    return fprs, tprs
   
  def plot_train_rocs2(self):
    ''' plots train roc from rates files in outfiles dir '''
    rates_files = []
    outfiles_path = self.dirs[2]
    for root, dirs, files in os.walk(outfiles_path):
      for file in files:
        if 'roc_beta50' in file:
          rates_file  = os.path.join(root, file)
          rates_files.append(rates_file)
    plt.figure()
    diags = [0.1*x for x in range(11)]
    plt.plot(diags, diags, alpha = 0.5, color = 'orange')
    colors = ['b', 'g','r','c','m','y']
    color_skip = 1.0/len(rates_files)
    color_i = [i*color_skip for i in range(len(rates_files))]
    color_map = cm.rainbow(color_i)
    for i, rates_file in enumerate(rates_files[:]):
      names = self.path2names(rates_file)
      filename = names[1]
      fprs, tprs = self.parse_rates_file(rates_file)
      plt.scatter(fprs, tprs, alpha = 0.9, edgecolor = 'none', label = filename, color = color_map[i])
    plt.axis([0, 1, 0, 1])
    plt.xlabel('False positive rate')
    plt.ylabel('True positive rate')
    plt.legend(loc = 'lower right', frameon = False, fontsize = 10)
    figpath = self.dirs[0] + '/roc_trains.pdf'
    plt.savefig(figpath)
    #plt.show()    

  def predict(self, unseens, unseens_headers, beta):
    ''' returns prediction for unseeen data, params are trained on continuos
    data which means the unseens need to be continuos as well'''
    con_data, con_headers = self.clean_data_continuous(self.data, self.headers, 250)
    train_data, train_headers = self.adasyn(con_data, con_headers, beta = beta)
    train_l, train_h = self.split_l_h(train_data)
    params_l = self.get_params_dict(train_l, con_headers)
    params_h = self.get_params_dict(train_h, con_headers)
    prior_l, prior_h = float(len(train_l))/len(train_data), float(len(train_h))/len(train_data)
    #prior_l, prior_h = 0.5, 0.5
    predictions = []
    for items in unseens:
      probs_l = [prior_l]
      probs_h = [prior_h]
      sample_name = items[0]
      for i in range(1, len(items)):
        var_name = unseens_headers[i]
        mu_l, sigma_l = params_l[var_name]
        mu_h, sigma_h = params_h[var_name]
        x = items[i]
        prob_l = self.p_normal_cdf(x, mu_l, sigma_l)
        prob_h = self.p_normal_cdf(x, mu_h, sigma_h)
        probs_l.append(prob_l)
        probs_h.append(prob_h)
      join_l = self.join_p(probs_l)
      join_h = self.join_p(probs_h)
      p_ratio = join_h/join_l
      if p_ratio > 1:
        predicted_class = 'h'
      else:
        predicted_class = 'l'
      prediction = [sample_name, predicted_class, round(p_ratio, 4)]
      predictions.append(prediction)
    predictions  = sorted(predictions, key = lambda item: item[-1], reverse = True)
    return predictions


  def file2predictions(self, filepath, beta, ranker):
    ''' returns predictions for a given file. only consider pocket P_1'''
    con_data, con_headers = self.clean_data_continuous(self.data, self.headers,250)
    train_data, train_headers = self.adasyn(con_data, con_headers, beta = beta)
    train_l, train_h = self.split_l_h(train_data)
    params_l = self.get_params_dict(train_l, con_headers)
    params_h = self.get_params_dict(train_h, con_headers)
    prior_l, prior_h = float(len(train_l))/len(train_data), float(len(train_h))/len(train_data)
    #prior_l, prior_h = 0.5, 0.5
    unseens, unseens_headers = self.parse_desc_merged_txt(filepath)
    unseens, unseens_headers = self.get_pocket2(unseens, unseens_headers)
    predictions = []
    for items in unseens:
      probs_l = [prior_l]
      probs_h = [prior_h]
      sample_name = items[0]
      for i in range(1, len(items)):
        var_name = unseens_headers[i]
        x = items[i]
        if var_name in params_l and x != 0:
          mu_l, sigma_l = params_l[var_name]
          mu_h, sigma_h = params_h[var_name]
          #x = items[i]
          prob_l = self.p_normal_cdf(x, mu_l, sigma_l)
          prob_h = self.p_normal_cdf(x, mu_h, sigma_h)
          probs_l.append(prob_l)
          probs_h.append(prob_h)
      join_l = self.join_p(probs_l)
      join_h = self.join_p(probs_h)
      p_ratio = round(join_h/join_l,4)
      wp_ratio = round(p_ratio * join_h, 4)
      if join_h == 0:
        join_h = 1e-14
      iwp_ratio = round(p_ratio/join_h, 4)
      if p_ratio > 1:
        predicted_class = 'h'
      else:
        predicted_class = 'l'
      #prediction = [sample_name, predicted_class, p_ratio, wp_ratio, iwp_ratio]
      prediction = [sample_name, predicted_class, p_ratio, wp_ratio]
      predictions.append(prediction)
    if ranker == 'wp':
      predictions  = sorted(predictions, key = lambda item: item[-1], reverse = True)
    else:
      predictions  = sorted(predictions, key = lambda item: item[-2], reverse = True)
    return predictions, ranker


  def file2predictions2(self, filepath, beta, ranker):
    ''' Uses kde, returns predictions for a given file. only consider pocket P_1'''
    con_data, con_headers = self.clean_data_continuous(self.data, self.headers,250)
    train_data, train_headers = self.adasyn(con_data, con_headers, beta = beta)
    train_l, train_h = self.split_l_h(train_data)
    kdes_l = self.get_kde_dict(train_l, con_headers)
    kdes_h = self.get_kde_dict(train_h, con_headers)
    prior_l, prior_h = float(len(train_l))/len(train_data), float(len(train_h))/len(train_data)
    #prior_l, prior_h = 0.5, 0.5
    unseens, unseens_headers = self.parse_desc_merged_txt(filepath)
    unseens, unseens_headers = self.get_pocket2(unseens, unseens_headers)
    predictions = []
    for items in unseens:
      probs_l = [prior_l]
      probs_h = [prior_h]
      sample_name = items[0]
      for i in range(1, len(items)):
        var_name = unseens_headers[i]
        if var_name in kdes_l:
          x = items[i]
          kde_l = kdes_l[var_name][0]
          kde_h = kdes_h[var_name][0]
          prob_l = kde_l(x)
          prob_h = kde_h(x)
          probs_l.append(prob_l)
          probs_h.append(prob_h)
      join_l = self.join_p(probs_l)
      join_h = self.join_p(probs_h)
      p_ratio = round(join_h/join_l,4)
      wp_ratio = round(p_ratio * join_h, 4)
      if join_h == 0:
        join_h = 1e-14
      iwp_ratio = round(p_ratio/join_h, 4)
      if p_ratio > 1:
        predicted_class = 'h'
      else:
        predicted_class = 'l'
      #prediction = [sample_name, predicted_class, p_ratio, wp_ratio, iwp_ratio]
      prediction = [sample_name, predicted_class, p_ratio, wp_ratio]
      predictions.append(prediction)
    if ranker == 'wp':
      predictions  = sorted(predictions, key = lambda item: item[-1], reverse = True)
    else:
      predictions  = sorted(predictions, key = lambda item: item[-2], reverse = True)
    return predictions, ranker


  def top_n_predicted(self, predictions, n):
    ''' writes outfile for top n predictions from predictions list '''
    top_n = predictions[:n]
    contents = ''
    contents += 'name,predicted_class,p_ratio \n'
    for prediction in top_n:
      content = ','.join(str(item) for item in prediction) + '\n'
      contents += content
    outdir = self.dirs[2]
    outname = outdir + '/top%s_predicted.txt' % n
    outfile = open(outname, 'w')
    outfile.write(contents)
    outfile.close()

  
  def file2top_predicted(self, filepath, n, beta, ranker):
    ''' writes outfile for top n preditions from a given file '''
    predictions, ranker = self.file2predictions(filepath, beta, ranker)
    names = self.path2names(filepath)
    top_n = predictions[:n]
    contents = ''
    #contents += 'name,predicted_class,p_ratio,wp_ratio,iwp_ratio\n'
    contents += 'name,predicted_class,p_ratio,wp_ratio\n'
    for prediction in top_n:
      content = ','.join(str(item) for item in prediction) + '\n'
      contents += content
    outdir = self.dirs[2]
    outname = outdir + '/%s_top%s_predicted_%s_%s.txt' % (names[1], n, ranker, beta)
    outfile = open(outname, 'w')
    outfile.write(contents)
    outfile.close()


  def file2top_predicted2(self, filepath, n, beta, ranker, outdir):
    ''' writes outfile for top n preditions from a given file '''
    predictions, ranker = self.file2predictions(filepath, beta, ranker)
    names = self.path2names(filepath)
    top_n = predictions[:n]
    contents = ''
    #contents += 'name,predicted_class,p_ratio,wp_ratio,iwp_ratio\n'
    contents += 'name,predicted_class,p_ratio,wp_ratio\n'
    for prediction in top_n:
      content = ','.join(str(item) for item in prediction) + '\n'
      contents += content
    #outdir = self.dirs[2]
    outname = outdir + '/%s_top%s_predicted_%s_%s.txt' % (names[1], n, ranker, beta)
    outfile = open(outname, 'w')
    outfile.write(contents)
    outfile.close()

### publication stuff ##

  def plot_hist_mult(self, data, headers, name):
    ''' plots multiple histograms on a single plot '''
    l_data, h_data = self.split_l_h(data)
    cl_datas = self.rows2columns(l_data, headers)
    ch_datas = self.rows2columns(h_data, headers)
    cl_data = cl_datas[0]
    ch_data = ch_datas[0]
    for i in range(1, len(headers))[:]:
      fignum = i
      rownum, colnum = 5,3
      plt.subplot(rownum, colnum,fignum) 
      plt.hist(cl_data[i], bins = 50, alpha = 0.5, edgecolor = 'none', label='l')
      plt.hist(ch_data[i], bins = 50, alpha = 0.5, edgecolor = 'none', label='h')
      plt.tick_params(axis = 'both', labelsize  = 6) 
      #plt.legend()
      plt.locator_params(axis='x', nbins=4)
      plt.locator_params(axis='y', nbins=4)
      title = headers[i]
      plt.title(title, fontsize=10)
      plt.xlabel('values',fontsize=6)
      plt.ylabel('frequency',fontsize=6)
    plt.tight_layout()
    figname = name + '_mult_hist.pdf'
    figpath = self.dirs[0] + '/' + figname
    plt.savefig(figpath)
    plt.close()

  def tab_cmatrix(self,cmatrix, name):
    '''tables a confussion matrix '''
    tp,fp = cmatrix[1][1], cmatrix[1][2]
    tn,fn = cmatrix[0][1], cmatrix[0][2]
    headers = ['Predicted high', 'Predicted low']
    cmatrix2 = [['True high',tp,fp],['True low',fn,tn]]
    tablecontent = tabulate(cmatrix2, headers, tablefmt='latex')
    tablename = name + '_cmatrix.tex'
    tablepath = self.dirs[1] + '/' + tablename
    tableout = open(tablepath, 'w')
    tableout.write(tablecontent)

  def tab_fp_tp(self,rates, name):
    '''tables tpr and fpr rates, rates are in the format [fpr, tpr] '''
    headers = ['FPR', 'TPR']
    rates = [rates, ['','']]
    tablecontent = tabulate(rates, headers, tablefmt='latex')
    tablename = name + '_fp_tp.tex'
    tablepath = self.dirs[1] + '/' + tablename
    tableout = open(tablepath, 'w')
    tableout.write(tablecontent)

  def plot_train3_roc_rates(self):
    ''' returns paths for train3 fpr and tpr (rates) files''' 
    rates_paths = self.find_files_from_path('roc_beta50_train3',  self.dirs[2])
    fp_tp_list = []
    plotnames = []
    for rates_path in rates_paths:
      names = self.path2names(rates_path)
      filenames = names[-1].split('_')
      plotname = filenames[-1][:-4]
      rates = self.parse_rates_file(rates_path)
      fp_tp_list.append(rates)
      plotnames.append(plotname)
    color_numbers = len(plotnames)
    color_skip = 1.0/color_numbers 
    color_indices = [i*color_skip for i in range(color_numbers)]
    color_map = cm.rainbow(color_indices)
    for i,plotname in enumerate(plotnames):
      fp, tp = fp_tp_list[i]
      color_index = color_map[i]
      #plt.scatter(fp, tp, alpha = 0.5, label=plotname, edgecolor='none', color=color_index)
      plt.scatter(fp, tp, alpha = 0.75, label=plotname, edgecolor='none', color=color_index)
    plt.plot([0,1],[0,1], color = 'orange')
    plt.axis([0,1,0,1])
    plt.legend(frameon=False, loc='lower right', fontsize=10)
    plt.xlabel('FPR')
    plt.ylabel('TPR')
    figname = 'roc_mult_scatter2.pdf'
    figpath = self.dirs[0] + '/' + figname
    plt.savefig(figpath)


  def ranker_test(self): 
    ''' returns rankers distribution for the given test betas'''
    test_betas = [0.65, 6.5] 
    rankers = ['p', 'wp'] 
    iters = range(100)
    ranker_predictions = {} 
    for ranker in rankers:
      beta_results = []
      means = []
      sds = []
      for beta in test_betas:
        beta_result = []
        for i in iters:
          #predictions, ranker = self.file2predictions('/Users/rahmadakbar/uc/enri/enri_rc9/infiles/train3.txt',beta,ranker)
          predictions, ranker = self.file2predictions('/Users/rahmadakbar/uc/enri/enri_rc9/infiles/hddata.tsv',beta,ranker)
          top10 = predictions[:10]
          hcount = 0.0
          for item in top10:
            if item[0][0] == 'h':
              hcount += 1
          fraction = hcount/10
          beta_result.append(fraction)
        mean = self.mean(beta_result)
        sd = self.sd(beta_result)
        beta_results.append(beta_result)
        means.append(mean)
        sds.append(sd)
      ranker_predictions[ranker] = beta_results, means, sds, test_betas
    for i, ranker in enumerate(ranker_predictions):
      fignum = i
      rownum, colnum = 1,2
      plt.subplot(rownum,colnum,fignum)
      fractions, means, sds, betas = ranker_predictions[ranker]
      label0 = 'beta %s, mean  %s, sd %s' % (betas[0], round(means[0],3), round(sds[0],3))
      label1 = 'beta %s, mean  %s, sd %s' % (betas[1], round(means[1],3), round(sds[1],3))
      plt.axis([0,100,0,100])
      x1,x2 = np.array(fractions[0])*100, np.array(fractions[1])*100
      plt.hist(x1, bins=25, edgecolor='none', alpha=0.5, label=label0)
      plt.hist(x2, bins=25, edgecolor='none', alpha=0.5, label=label1)
      plt.legend(frameon=False, fontsize = 10)
      plt.title(ranker)
      plt.xlabel('percent correct in top 10')
      plt.ylabel('frequency')
    plt.tight_layout()
    #plotname = 'rankers_hist_train3.pdf'
    plotname = 'rankers_hist_train3_hd.pdf'
    plotpath = self.dirs[0] + '/' + plotname
    plt.savefig(plotpath)
    plt.close()
    os.system('open %s' % plotpath)


  def ranker_test_kde(self): 
    ''' returns rankers distribution for the given test betas'''
    test_betas = [0.65, 6.5] 
    rankers = ['p', 'wp'] 
    iters = range(100)
    ranker_predictions = {} 
    for ranker in rankers:
      beta_results = []
      means = []
      sds = []
      for beta in test_betas:
        beta_result = []
        for i in iters:
          predictions, ranker = self.file2predictions2('/Users/rahmadakbar/uc/enri/enri_rc9/infiles/train3.txt',beta,ranker)
          top10 = predictions[:10]
          hcount = 0.0
          for item in top10:
            if item[0][0] == 'h':
              hcount += 1
          fraction = hcount/10
          beta_result.append(fraction)
        mean = self.mean(beta_result)
        sd = self.sd(beta_result)
        beta_results.append(beta_result)
        means.append(mean)
        sds.append(sd)
      ranker_predictions[ranker] = beta_results, means, sds, test_betas
    for i, ranker in enumerate(ranker_predictions):
      fignum = i
      rownum, colnum = 1,2
      plt.subplot(rownum,colnum,fignum)
      fractions, means, sds, betas = ranker_predictions[ranker]
      label0 = 'beta %s, mean  %s, sd %s' % (betas[0], round(means[0],3), round(sds[0],3))
      label1 = 'beta %s, mean  %s, sd %s' % (betas[1], round(means[1],3), round(sds[1],3))
      plt.axis([0.2,1,0,100])
      plt.hist(fractions[0], bins=30, edgecolor='none', alpha=0.5, label=label0)
      plt.hist(fractions[1], bins=30, edgecolor='none', alpha=0.5, label=label1)
      plt.legend(frameon=False, fontsize = 10)
      plt.title(ranker)
      plt.xlabel('% correct in top 10')
    plt.tight_layout()
    plotname = 'rankers_train3kde_hist.pdf'
    plotpath = self.dirs[0] + '/' + plotname
    plt.savefig(plotpath)
    plt.close()

  def parse_tsv(self, filepath):
    ''' parses tsv file from a given path '''
    contents = open(filepath).read().splitlines()
    headers = contents[0].split('\t')
    newcontents = []
    for item in contents[1:]:
      newcontent = item.split('\t')
      newcontents.append(newcontent)
    return newcontents, headers

  def tab_ef(self, filepath):
    ''' tabulates ef_table.tsv'''
    data, headers = self.parse_tsv(filepath)
    table = tabulate(data, headers, tablefmt='latex')
    path, name, name_ext = self.path2names(filepath)
    tablepath = path + '/' + '%s.tex' % name
    tablefile = open(tablepath, 'w') 
    tablefile.write(table)
    tablefile.close()


 
  def pub_stuff(self):
    ''' do stuff for the publication'''
    clean_data, clean_headers = self.clean_data_continuous(self.pocket1, self.pocket1_headers, 250)
    self.plot_hist_mult(clean_data, clean_headers, 'train3_2')
#    predictions = self.predict(clean_data, clean_headers,0.5)
#    orig_data, orig_headers = self.adasyn(clean_data, clean_headers,10,0)
#    cmatrix1 = self.kfold_cross_validation2(orig_data, orig_headers)
#    balance_data, balance_headers = self.adasyn(clean_data, clean_headers,10,0.65)
#    cmatrix2 = self.kfold_cross_validation2(balance_data, balance_headers)
#    rates1, rates2 = self.tp_fp2(cmatrix1), self.tp_fp2(cmatrix2)
#    print rates1
#    self.tab_cmatrix(cmatrix1, 'train3_orig')
#    self.tab_cmatrix(cmatrix2, 'train3_balanced')
#    self.tab_fp_tp(rates1, 'train3_orig')
#    self.tab_fp_tp(rates2, 'train3_balanced')
#    self.plot_hist_mult(balance_data, balance_headers, 'train3_balanced2')
#    self.plot_train3_roc_rates()
#    self.ranker_test()
#    l, h = self.split_l_h(self.pocket1)
#    print len(l), len(h)
#    print len(self.pocket1)
#    self.tab_ef('/Users/rahmadakbar/uc/enri/enri_write/tables/ef_mmm.tsv')


#### vh feedbacks ####


  def plot_cds(self, cds1, cds2, header):
    ''' plots cumulative density values from gkde '''
    plt.plot([item[0] for item in cds1], [item[1] for item in cds1], label = 'l', color = 'blue')
    plt.plot([item[0] for item in cds2], [item[1] for item in cds2], label = 'h', color = 'green')
    plt.title(header, fontsize = 10)
    plt.locator_params(nbins=6)
    plt.tick_params(labelsize = 6)

  def kstest_twosamples(self, sample1, sample2, header):
    ''' ks- test 2 samples ,returns data points and the corresponding cdf'''
    min1, min2 = min(sample1), min(sample2)
    gkde1 = scstats.gaussian_kde(sample1) #pdf estimates 
    gkde2 = scstats.gaussian_kde(sample2) #pdf estimates 
    #self.plot_cfs(cf1,cf2, header)
    Ds = []
    cd1s = [] #cumulative densities
    cd2s = [] #cumulative densities
    if len(sample1) >= len(sample2):
      for item in sample1:
        cd1 =  gkde1.integrate_box_1d(min1,item)
        cd2 = gkde2.integrate_box_1d(min2,item)
        D = abs(cd1-cd2)
        Ds.append((D,item))
        cd1s.append((item, cd1))
        cd2s.append((item, cd2))
    else:
      for item in sample2:
        cd1 =  gkde1.integrate_box_1d(min1,item)
        cd2 = gkde2.integrate_box_1d(min2,item)
        D = abs(cd1-cd2)
        Ds.append((D,item))
        cd1s.append((item, cd1))
        cd2s.append((item, cd2))
    #self.plot_cds(cd1s,cd2s, header)
    maxD = max(Ds)
    return  maxD,Ds, cd1s, cd2s

  def plot_cd_mult(self, data, headers, name):
    ''' plots multiple density plots on a single plot '''
    l, h = self.split_l_h(data)
    l_cols, l_headers = self.rows2columns(l, headers)
    h_cols, h_headers = self.rows2columns(h, headers)
    nrows, ncols = 5,3
    sum_maxD = 0
    for i in range(2, len(l_cols)):
      maxD, cd1s, cd2s =  self.kstest_twosamples(l_cols[i], h_cols[i], headers[i])
      print maxD,  headers[i]
      sum_maxD += maxD[0]
      plot_number = i
      cd1s, cd2s = sorted(cd1s), sorted(cd2s)
      plt.subplot(nrows, ncols, plot_number)
      self.plot_cds(cd1s, cd2s, headers[i])
    print sum_maxD
    plt.tight_layout()
    figname = name + '_mult_plot.pdf'
    figpath = self.dirs[0] + '/' + figname
    plt.savefig(figpath)
    plt.close()

  def tabulate_maxD(self, data, headers, name):
    ''' tabulate cumulative densities '''
    l, h = self.split_l_h(data)
    l_cols, l_headers = self.rows2columns(l, headers)
    h_cols, h_headers = self.rows2columns(h, headers)
    maxDs = []
    for i in range(1, len(l_cols)):
      maxD,Ds, cd1s, cd2s = self.kstest_twosamples(l_cols[i], h_cols[i], headers[i])
      featurename = headers[i]
      Dvals = [item[0] for item in Ds]
      D, x = maxD[0], maxD[1]
      Ds = [item[0] for item in Ds]
      mu,sd = np.mean(Ds), np.std(Ds)
      #p = self.p_normal_cdf(x,0,sd)
      p1  = (1 + math.erf((x)/sd/math.sqrt(2)))/float(2)-1e-14
      p1 = (1-p1)
      p  = (1 + math.erf((mu)/sd/math.sqrt(2)))/float(2)-1e-14
      pval = format(1-p,'.3g')
      p2 =  stats.ttest_1samp(Dvals,0)
      print pval, p1, p2
      D = round(D,3)
      maxDs.append([featurename, D, pval])
    sys.exit()
    maxDs = sorted(maxDs, key=lambda item: item[1],reverse = True)
    table_headers = ['Feature', 'D', 'Pval']
    tablecontent = tabulate(maxDs, table_headers, tablefmt='latex')
    tablename = name + '_Dstatistics.tex'
    tablepath = self.dirs[1] + '/' + tablename
    tableout = open(tablepath, 'w')
    tableout.write(tablecontent)
    
  def hdstats_data(self, data, headers):
    ''' selects only high D statisitics features, returns only selected
    featuress as new data '''
    wanted_headers = ['name', 'poc_cov', 'ell c/a', 'hydrophobicity', 'simpleScore', 'drugScore']
    #wanted_headers = ['name', 'surface', 'lid', 'hull','ell c/a', 'hydrophobicity', 'simpleScore', 'drugScore']
    wanted_indices = []
    for i, header in enumerate(headers):
      if header in wanted_headers:
        wanted_indices.append(i)
    hdstats_data = []
    for sample in data:
      #new_item = [item[i] for i,item in enumerate(item) if i in wanted_indices]
      new_sample = []
      for i,item in enumerate(sample):
        if i in wanted_indices:
          new_sample.append(sample[i])
      hdstats_data.append(new_sample)
    return hdstats_data, wanted_headers 

  def write_tsv(self, data, headers, name):
    ''' writes a tab separated values file from the given data and headers'''
    content = ''
    content += '\t'.join(headers) + '\n'
    for item in data:
      new_item = '\t'.join([str(x) for x in item]) + '\n'
      content += new_item
    outname = name + '.tsv'
    outpath = self.dirs[2] + '/' + outname
    outfile = open(outpath, 'w')
    outfile.write(content)
    outfile.close()

  def vh_feedbacks(self):
    ''' do stuff from vh feedbacks'''
    clean_data, clean_headers = self.clean_data_continuous(self.pocket1, self.pocket1_headers, 250)
    l, h = self.split_l_h(clean_data)
    l_cols, l_headers = self.rows2columns(l, clean_headers)
    h_cols, h_headers = self.rows2columns(h, clean_headers)
#    print self.kstest_twosamples(l_cols[13], h_cols[13], l_headers[13])[0] 
#    self.plot_cd_mult(clean_data, clean_headers, 'cds_train3')
#    balance_data, balance_headers = self.adasyn(clean_data, clean_headers,10,0.65)
#    self.plot_cd_mult(balance_data, balance_headers, 'cds_train3_balanced')
#    self.tabulate_maxD(clean_data, clean_headers, 'train3')
#    self.tabulate_maxD(balance_data, balance_headers, 'train3_balance')
#    hdstats_data, hdstats_headers = self.hdstats_data(clean_data, clean_headers)
#    hdstatsbalance_data, hdstatsbalance_headers = self.adasyn(hdstats_data, hdstats_headers,10,0.65)
#    cmatrix_hdstats = self.kfold_cross_validation3(hdstats_data, hdstats_headers)
#    cmatrix_train3kde = self.kfold_cross_validation3(clean_data, clean_headers)
#    cmatrix_hdstatsbalance = self.kfold_cross_validation3(hdstatsbalance_data, hdstatsbalance_headers)
#    self.tab_cmatrix(cmatrix_hdstats, 'train3_hdstats_kde')
#    self.tab_cmatrix(cmatrix_hdstatsbalance, 'train3_hdstatsbalance_kde')
#    rates_hdstats, rates_hdstatsbalance = self.tp_fp2(cmatrix_hdstats), self.tp_fp2(cmatrix_hdstatsbalance)
#    self.tab_fp_tp(rates_hdstats, 'hdstats_kde')
#    self.tab_fp_tp(rates_hdstatsbalance, 'hdstatsbalance_kde')
#    rates_train3kde = self.tp_fp2(cmatrix_train3kde)
#    self.tab_fp_tp(rates_train3kde, 'train3kde')
#    cmatrix_train3balancekde = self.kfold_cross_validation3(balance_data, balance_headers)
#    rates_train3balancekde = self.tp_fp2(cmatrix_train3balancekde)
#    self.tab_fp_tp(rates_train3balancekde, 'train3balancekde')
#    self.write_tsv(hdstats_data, hdstats_headers, 'hdstat_data')
#    self.ranker_test() 
#    self.get_kde_dict(clean_data, clean_headers)

### miscellaneous functions ###

  def select_adescriptor(self, filepath, pocket, variable):
    '''selects and outputs feature according to the given arguments '''
    variable_index = ''
    data, headers = self.parse_train_txt2(filepath)
    coldata, colheaders = self.rows2columns(data,headers)
    select_data = [coldata[0]]
    select_headers = [colheaders[0]]
    for i,header in enumerate(colheaders):
      if variable == header:
        select_data.append(coldata[i])
        select_headers.append(header)
        break
    path, name, name_ext = self.path2names(filepath)
    outpath = path + '/' + name + '_' + variable + '.txt'
    outfile = open(outpath, 'w')
    content = select_headers[0] + '\t' + select_headers[1] + '\n'
    for i, item in enumerate(select_data[0]):
      itempocket = '_'.join(item.split('_')[-2:])
      if pocket == itempocket:
        content += item + '\t' + str(select_data[1][i]) + '\n'
    outfile.write(content)
    outfile.close
        

# vh feedbacks2#

  def write_latextable(self,tcontent,theaders, tname):
    '''
    Writes table in latex format to tables dir
    '''
    tpath = self.dirs[1] + '/%s.tex' % tname
    filecontent = tabulate(tcontent,theaders, tablefmt = 'latex')
    tablefile = open(tpath,'w')
    tablefile.write(filecontent)
    tablefile.close()


  def check_traindata(self):
    '''
    Checks train data, how many proteins, how many conformations in each class.
    '''
    trd, trh = self.parse_train_txt()
    cd, ch = self.clean_data_continuous(trd, trh,250)
    pdb = {}
    for datum in cd:
      name = datum[0]
      nameparts = name.split('_')
      #pdbid = '_'.join(nameparts[1:3])
      pdbid = nameparts[1]
      label = nameparts[0]
      if pdbid not in pdb:
        pdb[pdbid] = [0,0] # init dict, values [h,l]
      if label == 'h':
        pdb[pdbid][0] += 1
      elif label == 'l':
        pdb[pdbid][1] += 1
    total_l = 0
    total_h = 0
    theaders = ['pdbid','high', 'low']
    tcontents = [] 
    for pdbid in pdb:
      values = pdb[pdbid]
      total_l += values[1]
      total_h += values[0]
      tcontent = [pdbid] + values
      tcontents.append(tcontent)
    tcontents.append(['sum', total_h, total_l])
    filename = 'train_h_l'
    print pdb
    print len(pdb)
    print total_l + total_h
    self.write_latextable(tcontents, theaders,filename)

  def feedbacks2(self):
    '''    
    '''
    self.check_traindata()

#end vh feedbacs2#


#vh feedbacs3#


  def feedbacks3_wf(self):
    '''
    Workflow for feedback 3
    '''
    data, headers = self.clean_data_continuous(self.pocket1, self.pocket1_headers, 250)
    bdata, bheaders = self.adasyn(data, headers,10,0.65)
    #self.plot_hist_mult(data, headers,'train3')
    #self.plot_hist_mult(bdata, bheaders,'train3_balanced')
    #hddata, hdheaders = self.hdstats_data(data, headers)
    #self.tabulate_maxD(data, headers, 'train3')
    self.tabulate_maxD(bdata, bheaders, 'train3_balanced')
    #self.ranker_test()
    #self.ranker_test() #hddata
#end vh feedbacs3#


####### usage examples and development notes ######
    
e = Enri()
#print e.find_files('.pdb')
#e.rename_files('.pdb')
#e.pdb2desc('.pdb')
#e.name2firstcol('desc.txt')
#e.merge_edt()
#e.file2top_predicted('/Users/rahmadakbar/uc/enri/desc_merged/3l3x_desc_merged.txt',10,6.5,'wp')
#e.pub_stuff()
#e.vh_feedbacks()
#e.plot_adasyn(10,0.65)
#e.feedbacks2()
e.feedbacks3_wf() 
