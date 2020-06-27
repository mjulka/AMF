#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 14 12:36:41 2019

@author: Yulia
"""

#list of all genera from ncbi taxanomy database
#gb query  Glomeromycotina [subtree] AND genus [rank] 
from Bio import Entrez
from Bio import SeqIO


#how much data in NCBI to establish retmax
Entrez.email = "XXXX@gmail.com" #always say NCBI who you are, use real e-mail!
handle = Entrez.egquery(term="Glomeromycetes[subtree] AND genus [rank]")
record = Entrez.read(handle)
for row in record["eGQueryResult"]:
#    print(row["DbName"], row["Count"]) #for all data
    if row["DbName"]=="taxonomy":
        genera_count = (row["Count"])
        print(genera_count)
#search all genera in the high taxon
Entrez.email = "XXXX@gmail.com" #use real e-mail 
genus_handle = Entrez.esearch(db="Taxonomy", term="Glomeromycetes [subtree] AND genus [rank]", retmax=genera_count)
genus_record = Entrez.read(genus_handle)
genus_list = genus_record["IdList"]
print(genus_list) #list of txid
taxa_number=len(genus_list) #how much elements in genus_list
print(taxa_number)
#for each genera in the genus_list
count=0
while count < taxa_number:  
    print(count)
    #get genus name to name of files
    Entrez.email = "XXXX@gmail.com" 
    genus = Entrez.efetch(db="Taxonomy", id=genus_list[count], retmode="xml")    
    genus_info = Entrez.read(genus)
    genus_name=genus_info[0]['ScientificName']
    print(genus_name)
    #search nucleotides
    query=('txid'+genus_list[count]+'[Organism] AND (internal transcribed spacer 1 [Title] OR ITS1 [Title]) AND 5.8S [Title] AND (internal transcribed spacer 2 [Title] OR ITS2 [Title])')
    #how much entries in nucleotide database
    Entrez.email = "XXXX@gmail.com"
    handle = Entrez.egquery(term=query)
    record = Entrez.read(handle)
    for row in record["eGQueryResult"]:
#        print(row["DbName"], row["Count"]) #for all data
        if row["DbName"]=="nuccore":
            nuc_count = (row["Count"])
#            print(nuc_count)
    #take accession numbers
    if nuc_count != 0:
        Entrez.email = "XXXX@gmail.com"
        handle= Entrez.esearch(db="nucleotide", term=query, retmax=nuc_count)
        record_sqv = Entrez.read(handle)
#    print(record_sqv["Count"]) # same as nuc_count
        idlist = ",".join(record_sqv["IdList"])
       #download sequences
       #for genbank format
        Entrez.email = "XXXX@gmail.com"
        handle_gb = Entrez.efetch(db="nuccore", id=idlist, rettype="gb", retmode="text")
        text = handle_gb.read() #what will be in the file
        out_handle = open(genus_name +"_" + genus_list[count]+".gb", "w")
        out_handle.write(text)
        out_handle.close()
        #filter unidentified sqv and write in the GB file  
        names_test = SeqIO.parse (genus_name +"_" + genus_list[count]+".gb", "genbank") 
        full_name=[] #filter records with full taxonomy identification
        for seq_record in names_test:
            seq_id = seq_record.id
            species = seq_record.annotations["organism"]
            if "sp." in species:
                continue
            else:
                sequence=(seq_record.seq)
                full_name.append(seq_record)
        SeqIO.write(full_name, genus_name +"species.gb", 'genbank')
#filter bad sequenced with ambiguous nucleotides, resuts in genbank and fasta formats
        record = SeqIO.parse (genus_name +"species.gb", "genbank")
        gut = []
        bad = []
        for seq_record in record:
            seq_id = seq_record.id
            sqv=(seq_record.seq)
            unread = sqv.count('N')
            ambig2=sqv.count('Y') + sqv.count('R') + sqv.count('S') + sqv.count('W') + sqv.count('K') + sqv.count('M')
            ambig3=sqv.count('B') + sqv.count('D') + sqv.count('H') + sqv.count('V')
            ambiguous=unread+ambig2+ambig3
            if ambiguous > 3:
                bad.append(seq_record)
            else:
                gut.append(seq_record)
            SeqIO.write(bad, genus_name + '_ambiguous_nucleotide.gb', 'genbank')           
            SeqIO.write(gut, genus_name + '_filtered.gb', 'genbank') 
            SeqIO.write(gut, genus_name + '_filtered.fasta', 'fasta') 
    count+=1






       
