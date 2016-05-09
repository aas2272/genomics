import time
import numpy as np
import csv
import os
import statistics
import math


class sample_parser():

      def __init__(self, file_data, type, file_map, name = None):
            #constructor
            self.file_name = file_data #this is actual data
            self.tissue_name = name#this is tissue name, should be a list even if onl one element
            if type == 2:#meaning that it already saw this tissue set
                  self.gene_dict = self.parse_genes(file_map)#mapping of samples to types
                  self.types = self.gene_dict.keys()#different types of tissues available
                  self.data = None#data will be obtained soon
                  self.header = None

            else:
                  self.data, self.header = self.get_data()#gets the data
                  self.gene_dict = self.parse_genes(file_map)
                  self.types = name



      def get_samples(self):
            return self.types

      def is_valid(self): #checks if all tissues inputted are valid tissues
            for x in self.tissue_name:
                  if x not in self.types:
                        return False

            return True

      def preprocess(self):#preprocess is a stand alone method that gets the data if need be
            self.data, self.header = self.get_data()

      def save_file(self):#saves the file in the format of tissues passed in separated by '_'
            look = '_'.join(self.tissue_name)
            with open('data/'+look+'.csv', 'w+') as f:
                  writer = csv.writer(f, delimiter = '\t',  quotechar = '"')
                  writer.writerows(self.data)

      def get_data(self):
            #reads in the data and turns it into a numpy array

            with open(self.file_name,'r') as f:
                  data_iter = csv.reader(f, delimiter = '\t', quotechar = '"')
                  data = [data for data in data_iter]

            header = data[0]
            print "making numpy:"

            my_data = np.array(data)

            print "done with numpy"

            return my_data, header

      def parse_genes(self, file_in):
            #parses through the mapping and makes a dictionary of tissue/samples


            f = open(file_in, 'r')

            tissues_as_key = {}
            for line in f:
                  whole_line = line.split()

                  if len(whole_line) < 2:
                        continue

                  value = whole_line[0]
                  key = ''
                  for x in whole_line[1:]:
                        key += x + '_'

                  try:
                        tissues_as_key[key[:-1]].append(value)

                  except:
                        tissues_as_key[key[:-1]] = []
                        tissues_as_key[key[:-1]].append(value)

            f.close()
            print "got tissues"
            return tissues_as_key




      def get_relevant_tissues(self, tissue = 'tis'):
            #gets only the tissues that are needed for this matrix, uses this list to only keep the samples specified
            if tissue == 'tis':
                  tissue = self.tissue_name
            else:
                  tissue = [tissue]
            tissue_index = [0,1]
            tissue_name = []
            for z in tissue:
                  for a in self.gene_dict[z]:
                        tissue_name.append(a)


            for x in range(2, len(self.header)):#first 2 columns are only name, description
                  if self.header[x] in tissue_name:
                        tissue_index.append(x)

            data = self.data[:, tissue_index]
            print "got cols"

            return data


      def get_relevant_genes(self):
            #keeps only genes where 70% of the data is above .1
            cols = self.data
            relevant_tissues = [cols[0]]

            for gene in cols[1:]:


                  ### this is for relevent genes
                  total = 0
                  for x in gene[2:]:
                        if float(x) > .1:
                              total += 1
                  if (total / float(len(gene) ))  > .7:
                        relevant_tissues.append(gene)

            print "have relevent genes"
            data = np.array(relevant_tissues)
            return data

      def get_var(self):
            #gets the minimum variance of the rows
            min_var = float('inf')
            max_var = 0.0
            sum_var = 0.0

            num_rows = 0

            variances = []
            for x in self.data[1:,2:]:
                  num_rows += 1
                  data = x.astype(float)
                  # print data
                  var = statistics.variance(data)
                  variances.append(var)
                  if var < min_var:
                        min_var = var

                  if var > max_var:
                        max_var = var
                  sum_var += var
            return variances

