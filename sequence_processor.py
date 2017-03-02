
# coding: utf-8

# In[1]:

#This is the script to remove non-coding RNA genes and incomplete sequences.
from Bio import SeqIO
from Bio.Alphabet import IUPAC

def conding_length_filter (file_name):
    valid_records = [] #Create a new blank list for those valid sequences
    for rec in SeqIO.parse(file_name, "fasta",IUPAC.unambiguous_dna): #Read a .fasta file and parse one sequence each time
        if  len(rec)==151 and rec.seq[50:53] == "ATG": # Check if the 51-53th are ATG and if the total length is 151 nt
            rec.description="valid_"+' '.join(rec.description.split()[2:-6]) 
    # Make minor changes to the description line by adding "valid" and removing unnecessary information. 
            valid_records.append(rec) # Put the valid sequence and its description to the list "valid_records"
    SeqIO.write(valid_records, "valid_"+file_name, "fasta") # At last output the list of valid sequences as a new .fasta file


# In[ ]:




# In[ ]:

#This is the script to divide sequences into two groups based on the promoter elements.
from Bio import SeqIO
from Bio.Alphabet import IUPAC
import re

def Inr_m5_sorter (file_name):
    records_Inr = [] #Create a new blank list for sequences with Inr promoter
    records_m5 = [] #Create a new blank list for sequences with m5 promoter
    for rec in SeqIO.parse(file_name, "fasta",IUPAC.unambiguous_dna): #Read a .fasta file and parse one sequence each time
        if  re.search ('[TCA]CA[TCA][TA]', str(rec.seq[20:51])): 
    # Slice the sequence within 30 nt upstream the start condon,convert it to a string and check if there's an Inr element. 
            rec.description="Inr_"+rec.description[12:] # Label "Inr_" in the beginning of the description line
            records_Inr.append(rec) # Add the sequence with Inr promoter and its description to the list "records_Inr"
        elif  re.search ('CCTTT', str(rec.seq[30:51])): 
    # Slice the sequence within 20 nt upstream the start condon,convert it to a string and check if there's an m5 element.
            rec.description="m5_"+rec.description[12:] # Label "m5_" in the beginning of the description line
            records_m5.append(rec) # Add the sequence with m5 promoter and its description to the list "records_m5"
    SeqIO.write(records_Inr, "Inr_"+file_name, "fasta") # At last output the list of sequences with Inr element as a new .fasta file
    SeqIO.write(records_m5, "m5_"+file_name, "fasta") # At last output the list of sequences with m5 element as a new .fasta file


# In[ ]:

#This is the script to trim the sequence with Inr promoter to its transcription start site (TSS)
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.Alphabet import IUPAC
import re

def Inr_TSS (file_name):
    trimmed_records_Inr = [] #Create a new blank list for trimmed sequences with Inr promoter
    for rec in SeqIO.parse(file_name, "fasta",IUPAC.unambiguous_dna): #Read a .fasta file and parse one sequence each time
        rec.seq=Seq(re.search(r"[TCA]C(A[TCA][TA][ATCG]{0,25}ATG[ATCG]{98}$)",str(rec.seq)).group(1),IUPAC.unambiguous_dna) 
    #The adenine of ‘TCANWY’ in an Inr element is the TSS. Remove the 5' part and obtain the putative transcript.
    #Remember that a sequence is more than just a string. When you use a regular expression, you must turn the sequence 
    #into a string first, but before you save it in a .fasta file, you must change it back to a sequence again! 
        rec.description="trimmed_"+rec.description[12:]+" "+str(len(rec))+"nt" 
    # Label "trimmed" in the beginning of the description line and add length in the end.
        trimmed_records_Inr.append(rec) # Add the trimmed sequence with Inr promoter and its description to the list "trimmed_records_Inr"
    SeqIO.write(trimmed_records_Inr, "trimmed_"+file_name, "fasta") # At last output the above list as a new .fasta file
    


# In[ ]:

#This is the script to trim the sequence with m5 promoter to its transcription start site (TSS)
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.Alphabet import IUPAC
import re

def m5_TSS (file_name):
    trimmed_records_m5 = [] #Create a new blank list for trimmed sequences with m5 promoter
    for rec in SeqIO.parse(file_name, "fasta",IUPAC.unambiguous_dna): #Read a .fasta file and parse one sequence each time
        rec.seq=Seq(re.search(r"C(CTTT[ATCG]{0,15}ATG[ATCG]{98}$)",str(rec.seq)).group(1),IUPAC.unambiguous_dna) 
    #The second cytosine of ‘CCTTT’ in an m5 element is the TSS, Remove the 5' part and obtain the putative transcript.
        rec.description="trimmed_"+rec.description[12:]+" "+str(len(rec))+"nt" 
    # Label "trimmed" in the beginning of the description line and add length in the end.
        trimmed_records_m5.append(rec) # Add the trimmed sequence and its description with m5 promoter to the list "records_m5"
    SeqIO.write(trimmed_records_m5, "trimmed_"+file_name, "fasta") # At last output the list as a new .fasta file


# In[ ]:

#This is the script to shuffle each trimmed sequence in 1000 different ways
import random
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio.Alphabet import IUPAC

def shuffle_1000 (file_name):
    shuffled_records = [] #Create a new blank list for the shuffled sequences.
    for rec in SeqIO.parse(file_name, "fasta",IUPAC.unambiguous_dna): #Read a .fasta file and parse one sequence each time
        for i in range(1000): #Create a variable i keeping the index of each shuffled sequence, run the loop for 1000 times.
            rec_list = list(rec.seq) #Convert the sequence to a list
            random.shuffle(rec_list) #Shuffle the elements in the list
            shuffled_rec = SeqRecord(Seq("".join(rec_list), rec.seq.alphabet),id="Shuffled_%i" % (i+1),description="Based on %s" % rec.id)
    #Join the shuffled list to a string, convert the string to a sequence format and add description to the .fasta sequence.
            shuffled_records.append(shuffled_rec) #Add the shuffled sequence and its description to the list "shuffled_records"
    SeqIO.write(shuffled_records, "shuffled_"+file_name, "fasta") # At last output the list as a new .fasta file


# In[ ]:

#This is the script to produce a series of sliding windows with the width of 30 bp for each sequence
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio.Alphabet import IUPAC

from EdwardBio import 

def window_30nt (file_name):
    window_records=[] #Create a new blank list for the windows of every sequence
    for rec in SeqIO.parse(file_name, "fasta",IUPAC.unambiguous_dna): #Read a .fasta file and parse one sequence each time
        for i in range(len(rec)-29): #Create a variable i keeping the index of each window, run the loop for n times based on the length of the sequence.
            index=101-len(rec)+i #The index of each window is decided by its relative location to the start codon "ATG"
            window_rec= SeqRecord(Seq(str(rec.seq)[i:i+30], rec.seq.alphabet),id="Window_%i" % index,description="The %i window of %s" %(index,rec.id))
    #Convert the sequence to a string, slice the string to 30 bp windows based on the index and add description to each window. 
            window_records.append(window_rec) #Add each window and its description to the list "window_records"
    SeqIO.write(window_records, "window_"+file_name, "fasta") # At last output the list of windows as a new .fasta file


# In[ ]:

import pandas as pd
from pandas import DataFrame
folding_energy = pd.read_csv('energy_output.csv',index_col=0,header=None)
#Read a .csv file. "index_col=0" makes the first column the index for each row. "header=None" says there's no colomn index.
column1=DataFrame(folding_energy.mean(axis=1), columns=['Mean'])
#Calculate the mean value of each row, add index "Mean" on top of this column and save it to a variable as a DataFrame
column2=DataFrame(folding_energy.std(axis=1), columns=['STD'])
#Calculate the mean value of each row, add index "STD" on top of this column and save it to a variable as a DataFrame
column=pd.concat([column1, column2], axis=1)
#Concatenate the above two columns together
column.to_csv('energy_statistaics.csv')
#write the result as a new .csv file.

