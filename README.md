# Mapping of SNP with G4-quadruplexes
The G4Hunter code for python can be found in an independent repository at https://github.com/AnimaTardeb/G4Hunter.git

- [Requirments](#requirments)
- [Data Entry](#data-entry)
- [Instructions](#instructions)
- [Procedure](#procedure)
## **Requirments** 
- Python 3.9 recommended  
- Biopython  
- Matplotlib  
- Numpy 
## **Data Entry**
- Fasta files with 1000 nucleobase pairs before the promoter are accepted  
- Fasta files have to include capital letters for bases in the exon regions and small letters for bases in the inton regions for the code to function properly
- The location of the selected SNPs should be in the cDNA

## **Instructions**
1. Installation
> sudo pip install matplotlib  
>sudo pip install numpy
2. Make sure each fasta file represents one whole gene and that all of the fasta files are in the same directory.
3. Add the desired SNPs of the genes in a csv file as the following
> ![Table](https://i.postimg.cc/yNTtgbyF/table.png)
4. Open the terminal and make sure that the current directory includes the code  
5. Run the code by the following command
> python <"code_name.py"> -i <"inputrepository"> -o <"outputrepository"> -w <"window size"> -s <"score threshold">  
> inputrepository: Folder path that contains the fatsa files  
> outputrepository: Folder path to store the results of the G4Hunter tool  
6. Enter the direct path to the csv file  
>Enter folder path for table.csv: #input here
<br>

# **Procedure**

```python
#!/bin/py
import csv
import os, sys
import shutil
import re
import numpy as np
from G4Hunter import main,Soft # imports the function of the G4Hunter tool
```
# Functions 
 
>Functions are called inside the conditional block statement: if __name__ == "__main__": 

## get_result()
The function **get_result()** takes 4 parameters and computes the genomic distance of a particular SNP of a gene from the start codon.  
- It loops in the fasta file and counts the capital bases after the start codon as **exon**, all the bases after the start codon as **distance** and all the bases before the start codon as **dist_start_codon**.
- Once the number of exon reaches the distance of the SNP in the cDNA. Then the function checks if there is a particular distance upstream or downstream of the reached position that need to be traveled and recomputes the total distance accordingly.

```python
'''get_result(position of snp in cDna, x=distance to go backward or forward, g=gene name, index specifies which snp of the gene is being computed) function computes the genomic distance of a snp. '''

def get_result(pos_of_mut,x,g,index): 
   check=0 # used as a flag to indicate that the start codon was found in the file
   distance=0 # distance starts to increment with every base after the start codon is found
   exon=0 # exon is incremented for each capital letter that is found after the start codon
   dist_start_codon=0 # distance to reach the start codon from the beginning of the file
   for m,l in enumerate (li): 
       for i,k in enumerate (l[0:50]): # loops through every base pair in the file
           
           if (k.isupper() and exon <pos_of_mut): # exon is incremented for every capital base found after the start codon until it reaches the distance of the SNP in the cDNA (stored in csv file)
               if ("ATG"  in l) and check==0: #finds the start codon
                   if i==49:
                       check=check+1
                   if l.find("ATG")<i:# if Capital letter is found after the start codon but still in the same line, then distance and exon are incremented
                       exon=exon+1 
                       distance=distance+1
                   if l.find("ATG")>i:# Bases before the start codon are computed seperately in the variable dist_start_codon
                       dist_start_codon=dist_start_codon+1
                elif check >0:
                   exon=exon+1
                   distance=distance+1
               elif ('ATG' not in l and check==0):   # bases before the start codon  
                   dist_start_codon=dist_start_codon+1

           elif (k.islower() and check ==0): #small letters => intronic region
               if ("ATG"  in l):
                   if i==49:
                       check=check+1
                   if l.find("ATG")<i: # small letters after the start codon; only the distance is incremeted
                       distance=distance+1
                   if l.find("ATG")>i: #small letters before the start codon having the start codon and the letters in the same line
                       dist_start_codon=dist_start_codon+1
                       
               elif ('ATG' not in l and check==0): #small letters before the start codon and are not in the same line
                   dist_start_codon=dist_start_codon+1
                        
           elif (exon< pos_of_mut and k in lowers and check>0):# Finding a small letter after detecting the start codon
               distance=distance+1
           if exon==pos_of_mut: # The distance of the coding region is found; is equivalent to the SNP location only when the next condition is false.
               if x!='': # If there is a particular distance to move backward or forward, this means that the SNP is found in the intronic region
                   destination=x 
                   destination=int(destination)
                   if destination <0: # Checking whether the next step is moving backward to find the SNP starting from the position reached in the coding region
                       destination=abs(destination)
                       distance=distance-destination  #total distance
                       use=i-1
                       
                       if destination<use: # if the index is still pointing to the same line reached after going backward
                           print("*The sequence that contains the snp is: "+ l[use-destination-2]+l[use-destination-1]+ "{!r}".format(l[use-destination])+l[use-destination+1]+l[use-destination+2]+l[use-destination+3])
                       elif destination>use: # if going backwards will cause the index to go previous lines
                           back_l=m # index of file
                           destination=destination-(use-1) 
                           back_l=back_l-1 #current line index -1
                      
                           while(destination>50): #while the distance of going backward is greater than moving 1 whole line
                               destination=destination-50 # going 50 bases backward is equivalent to moving 1 line upstream
                               back_l=back_l-1 # current line index -1
                           get_seq=li[back_l].strip()
                           gt_ln=len(get_seq)-destination 
                           print("*The sequence that conatains the snp is: "+get_seq[gt_ln-2]+get_seq[gt_ln-1]+get_seq[gt_ln]+"\"%s\""%get_seq[gt_ln+1]+get_seq[gt_ln+2]+get_seq[gt_ln+3]) # SNP location 
                       elif destination==use: # if going backward is equivalent to moving exactly one line backward from the index position
                           back_l=m
                           back_l=back_l-1
                           get_seq=li[back_l].strip()
                           print("*The sequence that conatains the snp is: "+get_seq[len(get_seq)-4]+get_seq[len(get_seq)-3]+get_seq[len(get_seq)-2]+"\"%s\""%get_seq[len(get_seq)-1])

                   elif destination>0: # Checking whether the next step is moving forward to find the SNP starting from the position reached in the coding region
                       distance=distance+destination # total distance
                       use=i-1
                       if destination<use:# if the index is still pointing to the same line reached after going forward
                           print("*The sequence that contains the snp is: "+ l[use+destination-2]+l[use+destination-1]+ "{!r}".format(l[use+destination])+l[use+destination+1]+l[use+destination+2]+l[use-destination+3])
                       elif destination>use: # if going forward will cause the index to go previous lines
                           for_w=m
                           destination=destination-(use-1)   
                           for_w=for_w+1 #current line index +1
                           while(destination>50): #while the distance of going forward is greater than moving 1 whole line
                               destination=destination-50 # going 50 bases forward is equivalent to moving 1 line downstream
                               for_w=for_w+1 #current line index +1
                           get_seq=li[for_w].strip()                 
                           gt_ln=len(get_seq)+destination                     
                           print("*The sequence that contains the snp is: "+get_seq[gt_ln-2]+get_seq[gt_ln-1]+get_seq[gt_ln]+"\"%s\""%get_seq[gt_ln+1]+get_seq[gt_ln+2]+get_seq[gt_ln+3])
                       elif destination==use: # if going forward is equivalent to moving exactly one line from the index position
                           for_w=m
                           for_w=back_l+1
                           get_seq=li[for_w].strip()
                           print("*The sequence that contains the snp is: "+get_seq[len(get_seq)-4],get_seq[len(get_seq)-3],get_seq[len(get_seq)-2]+"\"%s\""%get_seq[len(get_seq)-1])
               else:
                   print("*The sequence that contains the snp is: ",l[i-4],l[i-3],l[i-2],"\"%s\""% l[i-1],l[i],l[i+1])
               break
       if exon == pos_of_mut:
               break

   dist_start_codon=str(dist_start_codon)
   if index==0:
       gene_dic[g+"_"+dist_start_codon]=gene_dic[g]
       del(gene_dic[g])
   gene_dic[g+"_"+dist_start_codon][index]=gene_dic[g+"_"+dist_start_codon][index]+(distance,) # in here the key value would be a combination of the gene's name and the distance of its start codon (example: key=CYP_1200)

   return distance #total distance returned
#######################################################
```
## read_g()
**read_g** function takes 3 parameters: the size of the g4 file created by the G4HUnter tool, the name of the file and a dataframe **df** that contains the G4 sequences
The following steps are taken by this function:
- Reads the G4 file
- Stores the start and the end of all the predicted G4 sequences of a particular gene in a dictionary called **store** such that the key is the end position and the value is the start position 
- Calls **get_best** function that returns the closest positions of the G4 sequence to the SNPs positions of a gene  
- reassigns **store** as an empty dictionary to use it again for the next gene 
- repeats the steps for all the genes

```python
#######################################################
#Getting distance for G4 sequences of each gene

def read_g(s,file,df): #s=size of the g4 file f=filename
    g4=file
    if os.path.exists(g4): # check if the outputrepository is available
        file_g=open(g4,"r") # open the g4 file created by the G4Hunter tool
        cypg_reader = csv.reader(file_g, delimiter='\t')
        store={} #stores the end position of each g4 sequence as key and their start position as their values
        namee=""
        matching_key=""
        
        for i,size in enumerate (cypg_reader):
           
            
                if ">" in size[0]: # check if a new gene is reached
                    if len(store)!=0: # the dictionary is still empty for teh first gene
                        namee=namee.strip()
                        # extracting the keys of the dictionary as a numpy array of strings
                        keys = np.array(list(gene_dic.keys()))
                        if np.any(np.char.find(keys, namee) >= 0):  # checking if 'name' is in any of the keys using boolean indexing
                            matching_key = keys[np.char.find(keys, namee) >= 0][0]
                          
                            start_codon=matching_key.split('_')[1] #getting the distance of start codon for the gene
                            get_best(start_codon,matching_key,store,df,namee)
                            store={}
                        
                        
                    get_name=size[0]
                    namee=get_name[1:]
                    
                    print("\n",namee,":")
                
                elif "Start" not in size[0]: # "Start" not in e to skip the header
                    
                    store[int(size[1])]=int(size[0])     # size is a row that describes one G4 sequence. Hence size[1]=End position and e= start position
                if i==s:
                    namee=namee.strip()
                    keys = np.array(list(gene_dic.keys()))
                
                    # checking if 'name' is in any of the keys using boolean indexing
                    if np.any(np.char.find(keys, namee) >= 0):
                        matching_key = keys[np.char.find(keys, namee) >= 0][0]
                        start_codon=matching_key.split('_')[1] #getting the distance of start codon for the gene
                        get_best(start_codon,matching_key,store,df,namee)
            
        file_g.close()
```
## get_best()
**get_best()** function takes 5 parameters: distance of the start codon, key of the gene_dic **k**, a dictionary (**store**) that contains all the G4 sequences of a gene,a dataframe **df** used to index the chosen G4 sequence, and the gene name **gene_name**
- At first, the function checks whether the G4 sequence is either upstream or downstream of the start codon as it iterates over the G4 sequences
- Since the G4 file contains the position of the sequences not relative to the start codon, it should not be compared to the SNP locations yet. Thus, the distance has to be computed accordingly by substracting the **dist_start_codon** from the start position of the G4 sequence and adding the value in a new dictionary **dis_dic**
- Then, the distance between the SNP and the G4 sequence can be computed by either substracting the end position of the G4 sequence from the SNP position if the SNP is found downstream of the G4 sequence, or by substracting the SNP position from the start position of the g4 sequence if the SNP is found upstream of the G4 sequence.
- The computed distances between the SNP and the G4 sequences are stored in a list called **computed_dis**
- Finally the best distance which is the the least distance is extracted from **computed_dis**

```python     
def get_best(dist_start_codon,k,store,df,gene_name):
    dist_start_codon=int(dist_start_codon)
    
    distance=0
    
    computed_dis=[]
    dis_dic={}
      
    for G in store: # obtaining the distance of G4 complexes relative to the start codon and storing them in a dictionary (start= second column value in G4Hunter file)
       
        if G< dist_start_codon: # if the G sequence is found upstream from the start codon
            distance=-(dist_start_codon-G)
            dis_dic[distance]=store[G]
            
        else: #else it is found upstream
            distance=G-dist_start_codon
            dis_dic[distance]=store[G]
           
    temp_dictionary={}

    
    
    for t in gene_dic[k]: # computing distance between G4 complexes and snp distances to extract the closest distance.
        computed_dis=[]
        for g in dis_dic:  
            if g<t[3]:  #case where SNP position is downstream of the g4 sequence (end position)
                if g<0 and t[3]>0: #case if g4_sequence is upstream from the start codon and SNP is not
                        distance=t[3]+abs(g) # adding the distances together because G4 is upstream of the start codon and SNP is not
                        computed_dis.append(distance)
                        temp_dictionary[distance]=dis_dic[g]
                elif g<0 and t[3]<0:  #case if g4_sequence and SNP are upstream from the start codon
                        distance=abs(g)-abs(t[3])
                        computed_dis.append(distance)
                        temp_dictionary[distance]=dis_dic[g]
                
                else: # none is upstream of the start codon
                    distance=t[3]-g # SNP position - End position of g4
                    computed_dis.append(distance)
                    temp_dictionary[distance]=dis_dic[g]
                    
                    
            else: # SNP is upstream of the end position of g4 sequence => start position of g4 is used
                if g>0:
                    g4_start=dis_dic[g]-dist_start_codon # start position => computing the distance from the start codon to the start of the g4 sequence
                    if t[3]<0: # SNP position is upstream of the start codon
                        distance=abs(t[3])+g4_start  # adding the distances together because SNP is upstream of the start codon and g4 is not
                        computed_dis.app
                        computed_dis.append(distance)
                        temp_dictionary[distance]=dis_dic[g]
                
                    else:# SNP and G4 are downstream the start codon
                        if t[3]>g4_start: # SNP is upstream of the start position but downstream of the end position of g4 sequence => SNP and g4 sequnece overlap distance=0
                            distance=0
                            computed_dis.append(distance)
                            temp_dictionary[distance]=dis_dic[g]
                        else: # SNP is upstream of the start position of g4 sequence
                            distance=g4_start-t[3]  #Start position of g4 - SNP position
                            computed_dis.append(distance)
                            temp_dictionary[distance]=dis_dic[g]
                elif t[3]<0 and g<0: # both are upstream of the start codon
                    g4_start=dist_start_codon-dis_dic[g]

                    distance=abs(t[3])-abs(g4_start)
                    computed_dis.append(distance)
                    temp_dictionary[distance]=dis_dic[g]
        
        best=min(computed_dis) # best computed distance is extracted
        if t[1]=="":
            if best !=0:
                print("-The G4 sequence that has the closest distance of ",best," with the SNP (c.",t[0].strip(),t[2].strip(),") is: \n",df.loc[(df['Gene']==gene_name) & (df['Start']==temp_dictionary[best])].iloc[:,1:].to_string(index=False))
            else:
                print("-This G4 sequence overlaps with the SNP (c.",t[0].strip(),t[2].strip(),"): \n",df.loc[(df['Gene']==gene_name) & (df['Start']==temp_dictionary[best])].iloc[:,1:].to_string(index=False))
        else:
           if best !=0:
               print("-The G4 sequence that has the closest distance of ",best," with the SNP (c.",t[0].strip(),t[1].strip(),t[2].strip(),") is: \n",df.loc[(df['Gene']==gene_name) & (df['Start']==temp_dictionary[best])].iloc[:,1:].to_string(index=False))
           else:
               print("-This G4 sequence overlaps with the SNP (c.",t[0].strip(),t[1].strip(),t[2].strip(),"): \n",df.loc[(df['Gene']==gene_name) & (df['Start']==temp_dictionary[best])].iloc[:,1:].to_string(index=False))
    #######################################################
    End of get_best() 
    #######################################################
```

The script starts running from this block
```python
#######################################################
if __name__ == "__main__": 
    try:
     inputrepository, outputrepository , window, score = main(sys.argv[1:]) # sys.argv is a list in Python that contains the command-line arguments passed to the script, and sys.argv[1:] slices the list to exclude the first argument (the script name). The resulting values returned by main() are being unpacked into the variables inputrepository, outputrepository, window size, and score threshold.
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
```
This Part of the code reads the csv file provided by the user and creates a dictionary called **gene_dic** where the keys are initialy the gene names alone and the values are a list of tuple, such that every tuple represents one SNP for the gene.
### gene_dic representation  
> 'CYP2D6': [('506', '-1', 'G>A')]=> *CYP2D6 has one SNP as c.506-1G>A*

```python
'''This part below of the code reads the table.csv to get the snps for every gene. A gene dictionary is created (gene_dic), where the key is the name of the gene and the value is a list of tuple such that every tuple is a snp.
name of gene:[(snp, distance to go backwards or forwards ?,variation)]'''

    gene_dic={}
    table=input("Enter folder path for table.csv:") 
    with open(table, newline='') as csvfile:
        # Create a CSV reader object
        csvreader = csv.reader(csvfile, delimiter=',')

    # Iterate over each row in the CSV file
        for row in csvreader:
        # Process the row data
            gene_snps= filter(lambda cell: cell.strip() != '', row) # every row contains gene name and all its SNPs
 
            gene_snps=list(gene_snps) # converts filter object to list
       
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
``` 
This part of the code is taken from the G4Hunter tool. It creates the G4 file for the genes in the outputdirectory in a folder called Results
``` python 
#This part below of the code creates a directory where the results for G4 hunter will be stored.

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
        
        namee = os.path.basename(filefasta[0])
        Res1file= open (outputrepository +"/"+DIR+"/"+"-G4_sequences"+".txt", "a")
        Res2file= open (outputrepository +"/"+DIR+"/"+"-G4_Merged.txt", "a")
    #=========================================
        soft1=Soft()
        ScoreListe, DNASeq, NumListe, HeaderListe=soft1.GFinder(filein, window)

        for i in range(len(DNASeq)):
            G4Seq=soft1.GetG4(DNASeq[i],Res1file, ScoreListe[i], float(score), int(window),filefasta[0],len(NumListe[i]))
            if (len(G4Seq)>0):
                MSCORE=soft1.WriteSeq(DNASeq[i],Res2file,ScoreListe[i], G4Seq, filefasta[0], int(window), len(NumListe[i]))
        filein.close()
        print ("\n Results files and Score Figure are created in:   ")
        print (outputrepository,"/",DIR,"/","\n ")


        Res1file.close()
        Res2file.close() 
    #######################################################
```
This part of the code renames the fasta files in the inputdirectory from .fasta to .txt in order to be able to read the files
```python
#######################################################
#the code below converts the fasta files into txt files and calls get_result()

# Iterate over the list of files

    results=[]
    user_input=inputrepository
    
    files= os.listdir(user_input)
    os.chdir(user_input)
    
    
    for file in files:
        
        new_extension = ".txt"
        
    # Split the file name and extension
        file_name, extension = os.path.splitext(file)


    # Concatenate the file name with the new extension
        new_filename = file_name + new_extension
    
        

    # Rename the file
        os.rename(file, new_filename)
```
In this part of the code, as the fasta files are being iterated the gene name is extracted and is used as a key in **gene_dic** to further iterate over the values (SNP representations). Then, **get_result** function is called to compute the distance of this particular SNP of the gene.  
The new **gene_dic** after calling the **get_result** would be as the following  
> {CYP2D6_1019: [('506', '-1', 'G>A', 1847)]} *the name of the gene is proceeded by the **dist_start_codon** and the total distance of the SNP from the start codon in the genomic file is added in the tuple* 

Once the the **get_result** is used on all the files, they will be renamed again to .fasta instead of .txt  

```python

  # Construct the full file path
        file_path = os.path.join(user_input, new_filename)   
      
        if os.path.isfile(file_path):
            gene=new_filename.replace('.txt','')
            print("\n--Results for the ",gene, " gene: ")
     

            for index, value in enumerate  (gene_dic[gene]):
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

    files= os.listdir(user_input)
    os.chdir(user_input)
        
    for filename in files:
        fasta=".fasta"
        file_name, extension = os.path.splitext(filename)

    # Concatenate the file name with the new extension
        new_filename = file_name + fasta
    
        # rename the file
        os.rename(filename, new_filename)
```
This part of the code opens the G4 file created by G4Hunter and creates a dataframe for thcalls **read_g** function
```python
    # call read_g
    G_overlap=outputrepository+"/Results/-G4_Merged.txt"

    if os.path.abspath(G_overlap):
        absolute_path = os.path.abspath(G_overlap)
        data = []
        number_line=0
        current_gene = None
        with open(G_overlap) as f:
            for line in f:
                number_line=number_line+1
                if line.startswith('>'):
                    current_gene = line.strip().lstrip('>')
                elif current_gene:
                     fields = line.strip().split('\t')
                     if len(fields) <= 6 and 'Start' not in fields[0]:
                        start = int(fields[0])
                        end = int(fields[1])
                        sequence = fields[2]
                        length = int(fields[3])
                        score = float(fields[4])
                        data.append((current_gene, start, end, sequence, length, score))

        df = pd.DataFrame(data, columns=['Gene', 'Start', 'End', 'Sequence', 'Length', 'Score'])
        
        number_line=number_line-1
        read_g(number_line,absolute_path,df) 
        ```
