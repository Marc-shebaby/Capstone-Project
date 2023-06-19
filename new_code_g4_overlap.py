
import csv
import os, sys
import shutil

import pandas as pd
import re
import numpy as np
from G4Hunter import main,Soft
import requests
from Bio import Seq,Entrez
import requests

#!/bin/py

'''#######################################################'''

'''#######################################################
get_result(position of snp in cDna, x=distance to go backward or forward, g=gene name, index specifies which snp of the gene is being computed) function computes the genomic distance of a snp. 
'''
def get_result(pos_of_mut,x,g,index): 
   check=0 # used as a flag to indicate that the start codon was found in the file
   distance=0 # distance is incremented when start codon is found
   exon=0 # exon is incremented for each capital letter that is found after the start codon
   dist_start_codon=0
   CDS_loc=gene_dic[g][-1] # position of the start codon
   CDS=0
   for m,l in enumerate (li):
    for i,k in enumerate (l[0:]):
            if(exon<pos_of_mut):
                if k.isupper() and check==0 and CDS<CDS_loc:
                    CDS+=1
                    dist_start_codon+=1
                elif k.islower() and check==0:
                    dist_start_codon+=1
                elif k.isupper and check!=0:
                    exon+=1
                    distance+=1
                elif k.islower() and check!=0:
                    distance+=1
                elif CDS==CDS_loc:
                    dist_start_codon=-1
                    exon=exon+2
                    distance=distance+2
                    check=1

            
            
            if exon==pos_of_mut:
                if x!='':
                   
                   destination=x 
                   destination=int(destination)
                   if destination <0:
                       destination=abs(destination)
                       distance=distance-destination
                       use=i
                       
                       if destination<use:
                            if use-destination+3 >=len(l):
                              print("*The sequence that contains the SNP is: "+ l[use-destination-2]+l[use-destination-1]+ "{!r}".format(l[use-destination]))
                            else:
                              print("*The sequence that contains the SNP is: "+ l[use-destination-2]+l[use-destination-1]+ "{!r}".format(l[use-destination])+l[use-destination+1]+l[use-destination+2]+l[use-destination+3])
                     
                       elif destination>use:
                           back_l=m
                           destination=destination-(use+1) #subtracting the index reached from destination
                           back_l=back_l-1
                      
                           while(destination> len(l)):
                               destination=destination-len(l)
                               back_l=back_l-1
                           get_seq=li[back_l].strip()
                           gt_ln=len(get_seq)-destination
                           if gt_ln+3>=len(get_seq):
                                print("*The sequence that contains the SNP is: "+get_seq[gt_ln-3]+get_seq[gt_ln-2]+"\"%s\""%get_seq[gt_ln-1]+get_seq[gt_ln])
                           else:
                                print("*The sequence that contains the SNP is: "+get_seq[gt_ln-3]+get_seq[gt_ln-2]+"\"%s\""%get_seq[gt_ln-1]+get_seq[gt_ln]+get_seq[gt_ln+1]+get_seq[gt_ln+2])
                         
                       elif destination==use:
                           back_l=m
                           back_l=back_l-1
                           get_seq=li[back_l].strip()
                           print("*The sequence that contains the SNP is: "+get_seq[len(get_seq)-4]+get_seq[len(get_seq)-3]+get_seq[len(get_seq)-2]+"\"%s\""%get_seq[len(get_seq)-1])
                   elif destination>0:
                       distance=distance+destination
                       use=i
                       if destination<use:
                            if use+destination+3 >= len(l):
                                print("*The sequence that contains the SNP is: "+ l[use+destination-3]+l[use+destination-2]+ "{!r}".format(l[use+destination-1]))
                            else:
                                print("*The sequence that contains the SNP is: "+ l[use+destination-3]+l[use+destination-2]+ "{!r}".format(l[use+destination-1])+l[use+destination]+l[use+destination+1]+l[use+destination+2])
                       
                       elif destination>use:
                           for_w=m
                           destination=destination-(len(l)-use) # len(l)-(use): subtracting the index reached from the size of the current line to obtain the number of bases that will be traveled after moving forward on the line.
                           for_w=for_w+1
                           while(destination>len(l)):
                               destination=destination-len(l)
                               for_w=for_w+1
                           get_seq=li[for_w].strip()                 
                           gt_ln=destination                     
                           if gt_ln+3>=len(get_seq):
                               print("*The sequence that contains the SNP is: "+"\"%s\""%get_seq[gt_ln]+get_seq[gt_ln+1]+get_seq[gt_ln]+get_seq[gt_ln+2]+get_seq[gt_ln+3]) 
                           else:
                               print("*The sequence that contains the SNP is: "+get_seq[gt_ln-2]+get_seq[gt_ln-1]+"\"%s\""%get_seq[gt_ln]+get_seq[gt_ln+1]+get_seq[gt_ln+2])
                       elif destination==use:
                           for_w=m
                           for_w=back_l+1
                           get_seq=li[for_w].strip()
                           print("*The sequence that contains the SNP is: "+get_seq[len(get_seq)-4],get_seq[len(get_seq)-3],get_seq[len(get_seq)-2]+"\"%s\""%get_seq[len(get_seq-1)])
                else:
                   if i+1>=len(l):
                        print("*The sequence that contains the SNP is: ",l[i-4],l[i-3],l[i-2], l[i-1],"\"%s\""%l[i])
                   else:
                        print("*The sequence that contains the snp is: ",l[i-4],l[i-3],l[i-2], l[i-1],"\"%s\""%l[i],l[i+1])
                break
    if exon == pos_of_mut:
               break
   # outside the for loops
   dist_start_codon=str(dist_start_codon)
   if index==0: # for the first SNP of the gene
       
       gene_dic[g+"_"+dist_start_codon]=gene_dic[g] # change the key to include the distnace of the start codon
       del(gene_dic[g])
   gene_dic[g+"_"+dist_start_codon][index]=gene_dic[g+"_"+dist_start_codon][index]+(distance,) # the distance is added to the tuple corresponding to the SNP

   return distance
'''#######################################################'''



'''#######################################################
#Getting distance for G4 sequences of each gene
'''
def read_g(s,file,df): # s= number of lines in the g4 file, file= path of the file, df=dataframe
    
    g4=file
    if os.path.exists(g4):
        

        file_g=open(g4,"r")
        
        cypg_reader = csv.reader(file_g, delimiter='\t')
        store={} #stores the end position of each g4 sequence as key and their start position as their values
        namee="" #name of the gene
        matching_key=""
        
        for i,size in enumerate (cypg_reader):
            #for n,e in enumerate (size):
            
                if ">" in size[0]: 
                    if len(store)!=0:
                        
                        namee=namee.strip() # in case of white spaces
                        # extracting the keys of the dictionary as a numpy array of strings
                        keys = np.array(list(gene_dic.keys()))
                        
                        # checking if 'name' is in any of the keys using boolean indexing
                        if np.any(np.char.find(keys, namee) >= 0):
                            matching_key = keys[np.char.find(keys, namee) >= 0][0]
                          
                            start_codon=matching_key.split('_')[1] #getting the distance of start codon for the gene
                            get_best(start_codon,matching_key,store,df,namee)
                            store={}
                        
                    get_name=size[0] 
                    namee=get_name[1:] #name of the gene
                    
                    print("\n"," The G4 sequences for ",namee,":")
                
                 
                elif "Start" not in size[0]: #to skip the line with column names
                    
                    store[int(size[1])]=int(size[0]) # end is the key and start is the value
                    
                    
                
                if i==s: #to call get_best for the last gene
                    namee=namee.strip()
                    keys = np.array(list(gene_dic.keys()))
                
                # checking if 'name' is in any of the keys using boolean indexing
                    if np.any(np.char.find(keys, namee) >= 0):
                        matching_key = keys[np.char.find(keys, namee) >= 0][0]
                        start_codon=matching_key.split('_')[1] #getting the distance of start codon for the gene
                        get_best(start_codon,matching_key,store,df,namee)
            
        file_g.close()
        
def get_best(start_codon,k,store,df,gene_name): # start_codon: distance of the start codon, k:key, store: dictionary ,df: dataframe and gene_name
    start_codon=int(start_codon)
   
    distance=0
    
    computed_dis=[]
    dis_dic={} # distance between G sequence and start codon as key, value is start of the g sequence
      
    for G in store: # obtaining the distance of G4 complexes relative to the start codon and storing them in a dictionary
       
        if G< start_codon: # if the G sequence is found upstream from the start codon
            distance=-(start_codon-G)
            dis_dic[distance]=store[G]
           
        else:
            distance=G-start_codon
            dis_dic[distance]=store[G]
           
 
    df["distance_btw_snp"] = None  #new column in dataframe to store distance between G4 and SNP
    
     # computing distance between G4 complexes and snp distances from previous code to extract the closest distance.
    for t in gene_dic[k]: # looping through the SNPS of the indexed gene
        computed_dis=[]
        
        for g in dis_dic:  
            if g<t[3]:  
                if g<0 and t[3]>0: #case if g4_sequence is upstream from the start codon
                        distance=t[3]+abs(g)
                        computed_dis.append(distance)
                        df.loc[df['Start'] == dis_dic[g],'distance_btw_snp']=distance
                     
                elif g<0 and t[3]<0:  #case if g4_sequence and SNP are upstream from the start codon
                        distance=abs(g)-abs(t[3])
                        computed_dis.append(distance)
                        df.loc[df['Start'] == dis_dic[g],'distance_btw_snp']=distance
                     
                
                else:
                    distance=t[3]-g
                    computed_dis.append(distance)
                    df.loc[df['Start'] == dis_dic[g],'distance_btw_snp']=distance
                
                    
            else:
                if g>0:
                    g4_start=dis_dic[g]-start_codon #computing the distance from the start codon to the start of the g4 sequence
                    if t[3]<0:
                        distance=abs(t[3])+g4_start
                        computed_dis.append(distance)
                        df.loc[df['Start'] == dis_dic[g],'distance_btw_snp']=distance
                      
                
                    else:
                        if t[3]>g4_start:
                            distance=0
                            computed_dis.append(distance)
                            df.loc[df['Start'] == dis_dic[g],'distance_btw_snp']=distance
                           
                        else:  
                            distance=g4_start-t[3]
                            computed_dis.append(distance)
                            df.loc[df['Start'] == dis_dic[g],'distance_btw_snp']=distance
                          
                elif t[3]<0 and g<0:
                    g4_start=start_codon-dis_dic[g]

                    distance=abs(t[3])-abs(g4_start)
                    computed_dis.append(distance)
                    df.loc[df['Start'] == dis_dic[g],'distance_btw_snp']=distance
                 
        
        best=min(computed_dis)
        if t[1]=="":
            if best !=0:
               print("-The G4 sequence that has the closest distance of ",best," with the SNP (c",t[0].strip(),t[2].strip(),") is: \n",df.loc[(df['Gene']==gene_name) & (df['distance_btw_snp']==best)].iloc[:,1:].to_string(index=False))
            else:
               print("-This G4 sequence overlaps with the SNP (c",t[0].strip(),t[2].strip(),"): \n",df.loc[(df['Gene']==gene_name) & (df['distance_btw_snp']==best)].iloc[:,1:].to_string(index=False))
        else:
           if best !=0:
               print("-The G4 sequence that has the closest distance of ",best," with the SNP (c.",t[0].strip(),t[1].strip(),t[2].strip(),") is: \n",df.loc[(df['Gene']==gene_name) & (df['distance_btw_snp']==best)].iloc[:,1:].to_string(index=False))
           else:
               print("-This G4 sequence overlaps with the SNP (c.",t[0].strip(),t[1].strip(),t[2].strip(),"): \n",df.loc[(df['Gene']==gene_name) & (df['distance_btw_snp']==best)].iloc[:,1:].to_string(index=False))
    '''#######################################################
    End of get_best() 
    #######################################################'''
# In this part of the code the UCSC Genomic browser API is used
# three functions to retrieve the sequences of genes


def get_gene(gene_name):
    url = f"https://api.genome.ucsc.edu/getData/track?track=knownGene&genome=hg38"
    response = requests.get(url)
    data = response.json()
    
    if 'error' in data:
        print("Error:", data['error'])
        return None
    
    
    gene_info = data['knownGene']# obtaining a list of dictionaries

    df= pd.DataFrame(gene_info)
    filter_gene=df[(df["geneName"]==gene_name)&(df["rank"]==1)].iloc[0]

    chrom=filter_gene['chrom'] #get chromosome
    strand=filter_gene['strand'] # + or - strand
    start=filter_gene['chromStart']
    end=filter_gene['chromEnd']
    # get transcript name
    transcript=filter_gene['name']
    exon_loc=filter_gene['chromStarts'].split(',') #convert the value of string into a list
    exon_sizes=filter_gene['blockSizes'].split(',')
  

    return chrom,strand,start,end,exon_loc,exon_sizes,transcript


def get_coding_strand(genome,start,end,chrom,strand,exon_loc,exon_sizes,gene_name,path):
    if strand=="-":
        url = f"https://api.genome.ucsc.edu/getData/sequence?genome={genome};start={start};end={end+1000};chrom={chrom}"
        response = requests.get(url)
        data = response.json()
        seq_object = Seq.Seq(data['dna']) #getting the sequence of the positive strand
        final_seq=fix_bases(seq_object,exon_loc,exon_sizes,strand)
        final_seq="".join(final_seq)
        final_seq=Seq.Seq(final_seq)
        reverse_complement = final_seq.reverse_complement() #getting the negative strand (coding)
        reverse_complement="".join(reverse_complement)
        result=""
        for i in range(0, len(reverse_complement), 50):
            result += reverse_complement[i:i+50] + "\n"
        with open( f"{path}/{gene_name}.fasta", "w") as file:
            gene_=">"+gene_name
            file.write(gene_+"\n")
            file.write(result)
    else:
         url = f"https://api.genome.ucsc.edu/getData/sequence?genome={genome};start={start-1000};end={end};chrom={chrom}"
         response = requests.get(url)
         data = response.json()
         seq_object = Seq.Seq(data['dna']) #getting the sequence of the positive strand
         final_seq=fix_bases(seq_object,exon_loc,exon_sizes,strand)
         final_seq="".join(final_seq)
         result=""
         for i in range(0, len(final_seq), 50):
            result += final_seq[i:i+50] + "\n"
         with open( f"{path}/{gene_name}.fasta", "w") as file:
            gene_=">"+gene_name
            file.write(gene_+"\n")
            file.write(result)
   
def fix_bases(seq,exon_loc,exon_sizes,strand):
    seq=seq.lower()
    seq_1000_before=0
    rest_seq=0
    if strand=="+":
        seq_1000_before=seq[:1000] # the first 1000 base pairs upstream the promoter are kept as lower case
        rest_seq=list(seq[1000:])
      
    else:
        rest_seq=list(seq)

    for i,c in enumerate(exon_loc):
         if i<len(exon_loc)-1:
            rest_seq[int(c):int(c)+int(exon_sizes[i])] = "".join(rest_seq[int(c):int(c)+int(exon_sizes[i])]).upper()#exons should be in capital letters
    if seq_1000_before!=0:
          seq=seq_1000_before+"".join(rest_seq)
    else:
        seq="".join(rest_seq)

    return seq

# fetch_cds_from_ncbi uses Entrez API to retrieve the position of the CDS (start codon)
def fetch_cds_from_ncbi(transcript_id,gene):
    Entrez.email = 'marcshababy02@gmail.com'  # Provide your email address
    handle = Entrez.esearch(db='nucleotide',term=transcript_id)

    record = Entrez.read(handle)
    print(handle)
    print(record["IdList"])
    handle.close()

    fetch_handle = Entrez.efetch(db='nucleotide',id=record["IdList"][0],rettype="gb")
    record=fetch_handle.read()
    indx=record.find("CDS")
    CDS=record[indx:indx+25]
    CDS=CDS.replace(" ","")
    match = re.search(r'\d+',CDS)
    fetch_handle.close()
    first_number=""
    if match:
        first_number = match.group()
        first_number=first_number.strip() # remove any spaces
        first_number=int(first_number)
    return first_number

# End of functions


'''#######################################################
#Calling G4 hunter Tool
'''
if __name__ == "__main__":
    try:
    
        inputrepository, outputrepository, window, score = main(sys.argv[1:]) # sys.argv is a list in Python that contains the command-line arguments passed to the script, and sys.argv[1:] slices the list to exclude the first argument (the script name). The resulting values returned by main() are being unpacked into the variables inputfile, outputfile, window, and score.
      

        fname=inputrepository.split("/")[-1]
        name=fname.split(".")
    except ValueError:
        print ('\033[1m' +"\n \t Oops! invalide parameters  \n" +'\033[0;0m')
        print ("--------------------------------------------------------------------\n")
        sys.exit()
    except UnboundLocalError:
        print ('\033[1m' +"\n \t Oops! invalide parameters  \n" +'\033[0;0m')
        print ("--------------------------------------------------------------------\n")
        sys.exit()
    
    '''
This part below of the code reads the table.csv to get the snps for every gene. A gene dictionary is created (gene_dic), where the key is the name of the gene and the value is a list of tuple such that every tuple is a snp.
name of gene:[(snp, distance to go backwards or forwards ?,variation)]
'''
    gene_dic={}
    table=input("Enter folder path for table.csv:") 
    with open(table, newline='') as csvfile:
        # Create a CSV reader object
       genes=[]
       
       csvreader=csv.reader(csvfile, delimiter=',')
       # Iterate over each row in the CSV file
       for row in csvreader:
        # Process the row data
            gene_snps= filter(lambda cell: cell.strip() != '', row) # every row contains gene name and all its SNPs
 
            gene_snps=list(gene_snps) # converts filter object to list
            genes.append(gene_snps[0])
       
            for index, snp in enumerate(gene_snps): #loop over every row to extract the snp positions.
                if index!=0:
                    cDNA_distance=snp
                    pattern = r'(\d+)([-|+]\d+)?(.*)'
                    match = re.search(pattern, cDNA_distance)
                    if match:
                        group1 = match.group(1)
                        group2 = match.group(2) if match.group(2) is not None else ''
                        group3=match.group(3)
                        if gene_snps[0] not in gene_dic:
                            gene_dic[gene_snps[0]]=[]
                            gene_dic[gene_snps[0]].append((group1,group2,group3))
                        else:
                            gene_dic[gene_snps[0]].append((group1,group2,group3))
           
                        '''#######################################################'''
        
# calling the api function #
    for gene_name in genes:
        chrom,strand,start,end,exon_loc,exon_sizes,transcript=get_gene(gene_name)
        get_coding_strand('hg38',start,end,chrom,strand,exon_loc,exon_sizes,gene_name,inputrepository)
        gene_dic[gene_name].append(transcript)
        cds=fetch_cds_from_ncbi(transcript)
        gene_dic[gene_name].append(cds)

    '''#######################################################'''

    '''
This part below of the code creates a directory where the results for G4 hunter will be stored.
'''
    OPF= os.listdir(outputrepository)
    
    flag=False
    for dir in OPF:
        DIR="Results"
        if dir== DIR:
            print ("true",DIR)
            flag=True
    if flag==True:
        shutil.rmtree(outputrepository+"\\"+DIR+"\\")
        os.makedirs(outputrepository+"\\"+DIR+"\\", mode=0o777)        #
        print ('\033[1m' +"\n \t Re-evaluation of G-quadruplex propensity with G4Hunter " +'\033[0;0m')
        print ("\n#####################################")
        print ("#    New Results directory Created  #")
        print ("#####################################\n")
    else:
        DIR="Results"
        os.makedirs(outputrepository+"\\"+DIR+"\\", mode=0o777)        #
        print ("\n########################################################################")
        print ("#                            Results directory Created                 #")
        print("########################################################################\n")
    
    
    files= os.listdir(inputrepository)
    for file in files:
  # Construct the full file path
        file_path = os.path.join(inputrepository, file)
        fname=file_path.split("\\")[-1]
        filefasta=fname.split(".")
        filein=open(file_path,"r")

        print ("\n Input file:", '\033[1m' + filefasta[0]+'\033[0;0m')
        # directory of output files
        namee = os.path.basename(filefasta[0])
        Res1file= open (outputrepository+"/"+DIR+"/"+"-G4_sequences"+".txt", "a")
        Res2file= open (outputrepository+"/"+DIR+"/"+"-G4_Merged.txt", "a")
    #=========================================
    
    
        soft1=Soft()
        ScoreListe, DNASeq, NumListe, HeaderListe=soft1.GFinder(filein, window)

        for i in range(len(DNASeq)):
            G4Seq=soft1.GetG4(DNASeq[i],Res1file, ScoreListe[i], float(score), int(window),filefasta[0],len(NumListe[i]))# HeaderListe[i], len(NumListe[i]))
            if (len(G4Seq)>0):
                MSCORE=soft1.WriteSeq(DNASeq[i],Res2file,ScoreListe[i], G4Seq, filefasta[0], int(window), len(NumListe[i]))
        filein.close()

        print ("\n Results files and Score Figure are created in:   ")
        print (outputrepository,"/",DIR,"/","\n ")


        Res1file.close()
        Res2file.close() 
    '''#######################################################
    end of tool

    '''######################################################

# Iterate over the list of files

    results=[]
    files= os.listdir(inputrepository)
    os.chdir(inputrepository)
    
    
    for file in files:
        
        new_extension = ".txt"
        
    # Split the file name and extension
        file_name, extension = os.path.splitext(file)


    # Concatenate the file name with the new extension
        new_filename = file_name + new_extension
    
        

    # Rename the file
        os.rename(file, new_filename)


  # Construct the full file path
        file_path = os.path.join(inputrepository, new_filename)  
      
        if os.path.isfile(file_path):
            gene=new_filename.replace('.txt','')
            print("\n-- SNP location Results for the ",gene, " gene: ")
           
     

            for index, value in enumerate  (gene_dic[gene]):
                if index < len(gene_dic[gene])-1: # last index in the value is not a SNP
            
                    fasta_file=open(file_path,"r")
                    skip_header=fasta_file.readline()
                    li=fasta_file.readlines()
                    li = [line.replace('\n', '') for line in li]
                    if value[1]=="":
                        print("\n** snp c.",value[0].strip(),value[2].strip(),": ")
                    else:
                        print("\n** snp c.",value[0].strip(),value[1].strip(),value[2].strip(),":")
              
                    print("Distance between the snp and start codon is:",get_result(int(value[0]),value[1],gene,index),"base pairs")
          
                    fasta_file.close()
        else:  
         print("path doesnt exist") 
    
    # converting back to fasta
    files= os.listdir(inputrepository)
    os.chdir(inputrepository)
        
    for filename in files:
        fasta=".fasta"
        file_name, extension = os.path.splitext(filename)

    # Concatenate the file name with the new extension
        new_filename = file_name + fasta
    
        
        # rename the file
        os.rename(filename, new_filename)

    # call read_g
    G_overlap=outputrepository+"/Results/-G4_Merged.txt"

    if os.path.abspath(G_overlap):
        absolute_path = os.path.abspath(G_overlap)
        data = []
        number_line=0
        current_gene = None
        with open(G_overlap) as f: #to convert the file of the G4Hunter to a dataframe
            for line in f:
                number_line=number_line+1
                if line.startswith('>'):
                    current_gene = line.strip().lstrip('>')
                elif current_gene:
                     fields = line.strip().split('\t')
                     if 'Start' not in fields[0]:
                        start = int(fields[0])
                        end = int(fields[1])
                        sequence = fields[2]
                        length = int(fields[3])
                        score = float(fields[4])
                        data.append((current_gene, start, end, sequence, length, score))

        df = pd.DataFrame(data, columns=['Gene', 'Start', 'End', 'Sequence', 'Length', 'Score'])
        
        number_line=number_line-1
        read_g(number_line,absolute_path,df)
       
    else:
        print("Outputrepository not found")
