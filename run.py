import numpy as np
from sample_parser import sample_parser
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import nimfa
import Tkinter
import plotly.plotly as py
import plotly.graph_objs as go
import csv
import sys
import os
import random
import math

def make_rand(X,Y,var):
      #generates a random matric from a normal distribution with dimensions X,Y and variance var
      #######became obsolete later in the implementation######
      out_mat = []
      std = math.sqrt(var)
      for x in range(X):
            temp = []
            for y in range(Y):
                  # i = random.randint(0, X-1)
                  # j = random.randint(0, Y-1)

                  val = np.random.normal(0, std)

                  # val = data[i,j]
                  # val = data
                  temp.append(val)
            out_mat.append(temp)
      return np.array(out_mat)


def get_diff(a, b):
      #gets the difference between two matricies

      ###### make gamma distribution??


      sub = np.subtract(a, b)
      dist = np.linalg.norm(sub)

      return dist


def find_lim(data,cutoff):
      #finds which rank crosses the difference cutoff
      i = 1

      print str(cutoff) + ' cutoff'
      print 'rank, difference:'
      while(True):

            print "factoring with " + str(i)
            nmf = nimfa.Nmf(data, rank = i, seed="nndsvd", update='euclidean')
            nmf_fit = nmf()
            b = nmf_fit.basis()
            co = nmf_fit.coef()
            print "done facoring"


            comp = np.dot(b,co)
            diff = get_diff(data,comp)
            print i, diff
            if (diff < cutoff):
                  return i, b, co

            i+=1


def fraction(samp, basis, coef):
      # this is the additive part of the program below
      return (basis.sum())*coef / samp.sum()



def main(tissue_name, row_num):
      # main method has a lot of parts, many of which were commented out later on to save time


      look = '_'.join(tissue_name)
      writepath = 'data/' + look + '.csv'

      writepath2 = 'data/' + tissue_name[0] + '.csv'

      if os.path.exists(writepath2):
            print 'type 1'

            sp0 = sample_parser(writepath2, 1, file_map = 'tissues.txt', name = [tissue_name[0]])
            # sp0.data = sp0.data

      else:
            print 'type2'

            sp0 = sample_parser('data.gct', 2, file_map = 'tissues.txt', name = [tissue_name[0]])

            if sp0.is_valid():
                  sp0.preprocess()
                  sp0.data = sp0.get_relevant_tissues()
                  sp0.header = sp0.data[1,:]
                  sp0.save_file()


            else:

                  print "failure"
                  return -1

      sp0.data = sp0.get_relevant_genes()

      variances = sp0.get_var()

      min_variance = min(variances)
      max_variance = max(variances)
      avg_variance = sum(variances)/len(variances)
      variances.sort()

      med_variance = variances[len(variances)/2]
      first_per = variances[len(variances)/100]
      seven5_per = variances[len(variances)*75/100]

      eight5 = variances[len(variances)*85/100]
      nine6 = variances[len(variances)*96/100]
      nine7 = variances[len(variances)*97/100]
      nine8 = variances[len(variances)*98/100]
      ninenine_per = variances[len(variances)*99/100]



      print str(min_variance) + ' is min variance'
      print str(max_variance) + ' is max variance'
      print str(avg_variance) + ' is avg variance'
      print str(med_variance) + ' is med variance'
      print str(first_per) + ' is first percentile'
      print str(seven5_per) + ' is 75th variance'
      print str(eight5) + ' is 85th variance'
      print str(nine6) + ' is 96th percentile'
      print str(nine7) + ' is 97th variance'
      print str(nine8) + ' is 98th percentile'
      print str(ninenine_per) + ' is 99th percentile'


      if os.path.exists(writepath):
            print 'type 1'
            sp = sample_parser(writepath, 1, file_map = 'tissues.txt', name = tissue_name)


      else:
            print 'type2'
            sp = sample_parser('data.gct', 2, file_map = 'tissues.txt', name = tissue_name)
            if sp.is_valid():
                  sp.preprocess()
                  sp.data = sp.get_relevant_tissues()
                  sp0.header = sp0.data[1,:]
                  sp.save_file()

            else:
                  print "failure"
                  return -1

      sp.data = sp.get_relevant_genes()
      print "done parsing"




      tissue = sp.data
      tissuex = np.array(tissue)
      tissue = (tissuex[1:, 2:]).astype(float)

      print "getting difference"
      size = tissue.shape

      ######
      #prints out the norms of different randomly generated normal matricies
      ######

      # rand_mat = make_rand(size[0], size[1], variance)
      rand_mat = np.random.normal(0, math.sqrt(min_variance), (size[0], size[1]))
      #### do using built in for speed
      norm_rand = np.linalg.norm(rand_mat)
      print str(norm_rand) + ' norm of random min'


      rand_mat2 = np.random.normal(0, math.sqrt(max_variance), (size[0], size[1]))
      #### do using built in for speed
      norm_rand2 = np.linalg.norm(rand_mat2)
      print str(norm_rand2) + ' norm of random max'


      rand_mat3 = np.random.normal(0, math.sqrt(avg_variance), (size[0], size[1]))
      #### do using built in for speed
      norm_rand3 = np.linalg.norm(rand_mat3)
      print str(norm_rand3) + ' norm of random avg'

      rand_mat4 = np.random.normal(0, math.sqrt(med_variance), (size[0], size[1]))
      #### do using built in for speed
      norm_rand4 = np.linalg.norm(rand_mat4)
      print str(norm_rand4) + ' norm of random med'


      rand_mat5 = np.random.normal(0, math.sqrt(first_per), (size[0], size[1]))
      #### do using built in for speed
      norm_rand5 = np.linalg.norm(rand_mat5)
      print str(norm_rand5) + ' norm of first percentile'


      rand_mat6 = np.random.normal(0, math.sqrt(seven5_per), (size[0], size[1]))
      #### do using built in for speed
      norm_rand6 = np.linalg.norm(rand_mat6)
      print str(norm_rand6) + ' norm of 75 percentile'

      rand_mat11 = np.random.normal(0, math.sqrt(eight5), (size[0], size[1]))
      #### do using built in for speed
      norm_rand11 = np.linalg.norm(rand_mat11)
      print str(norm_rand11) + ' norm of 85 percentile'

      rand_mat8 = np.random.normal(0, math.sqrt(nine6), (size[0], size[1]))
      #### do using built in for speed
      norm_rand8 = np.linalg.norm(rand_mat8)
      print str(norm_rand8) + ' norm of 96th percentile'


      rand_mat9 = np.random.normal(0, math.sqrt(nine7), (size[0], size[1]))
      #### do using built in for speed
      norm_rand9 = np.linalg.norm(rand_mat9)
      print str(norm_rand9) + ' norm of 97th percentile'


      rand_mat10 = np.random.normal(0, math.sqrt(nine8), (size[0], size[1]))
      #### do using built in for speed
      norm_rand10 = np.linalg.norm(rand_mat10)
      print str(norm_rand10) + ' norm of 98th percentile'


      rand_mat7 = np.random.normal(0, math.sqrt(ninenine_per), (size[0], size[1]))
      #### do using built in for speed
      norm_rand7 = np.linalg.norm(rand_mat7)
      print str(norm_rand7) + ' norm of 99 percentile'



      #######
      # this is where the limit rank is found
      #######

      print "finding limit"
      rank, factor, load = find_lim(tissue, norm_rand10)



      ###### once i began hard coding rank to save time #######
      # rank = 4

      # print "factoring with " + str(rank)
      # nmf = nimfa.Nmf(tissue, rank = rank, seed="nndsvd", update='euclidean')
      # nmf_fit = nmf()
      # factor = nmf_fit.basis()
      # load = nmf_fit.coef()
      # print "done facoring"




##############
      # rank = 4

      out_tis = np.array(['name', 1, 2, 3, 4]).reshape(rank+1,1)

      big_data = [[['name'], ['1'], ['2'], ['3'], ['4']]] # this is column 1
      for i in range(len(tissue_name)):
            col = []
            small_data = []
            temp = (sp.get_relevant_tissues(tissue_name[i])[1:,2:]).astype(float)

            col.append(tissue_name[i])
            small_data.append([tissue_name[i]])

            for j in range(rank):
                  num = 0
                  total = 0
                  d = []
                  for k in range(temp.shape[1]):

                        d.append(fraction(temp[:, k], factor[:, j], load[j, k]))
                        total += fraction(temp[:, k], factor[:, j], load[j, k])
                        num += 1
                  small_data.append(d)
                  col.append(total/num)
            big_data.append(small_data)
            column = np.array(col).reshape(rank+1, 1)

            out_tis = np.concatenate((out_tis, column), axis = 1)

      outarray = np.asarray(out_tis)

      ##### saves data so i only needed to do it once

      with open('data/output.csv', 'w+') as f:
            writer = csv.writer(f, delimiter = ',',  quotechar = '"')
            writer.writerows(outarray)
      print 'done'

################




      print 'histogram making'


      # x = 1
      ###### this is the basis element being analyzed
      y = row_num
      # x2 = 3
      # y2 = 1

      rows = []
      for x in range(1, 5+1):


            row = go.Histogram(
                x=big_data[x][y],
                opacity=0.75,
                name=big_data[x][0][0],

            )
            rows.append(row)



      layout = go.Layout(
          title='Basis Element ' + big_data[0][y][0],
          xaxis=dict(
              title='Value'
          ),
          yaxis=dict(
              title='Count'
          ),
          barmode='overlay',
          bargap=0.25,
          bargroupgap=0.3
      )


      fig = go.Figure(data=rows, layout=layout)
      plot_url = py.plot(fig, filename='overlaid-histogram')







if __name__ == '__main__':

      main(sys.argv[1:-1], sys.argv[-1])


