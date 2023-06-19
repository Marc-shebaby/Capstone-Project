'''import pandas as pd
table=input("Enter folder path for table.csv:") 
with open(table, newline='') as csvfile:
        # Create a CSV reader object
       
        df = pd.read_csv(csvfile)
        gene_names = df.iloc[:, 0].tolist() #first column
        print(gene_names)'''
'''from pyensembl import EnsemblRelease

# Load the Ensembl release
ensembl = EnsemblRelease(98)

# Specify the gene name
gene_name = "CYP3A5"

# Get the gene information
gene = ensembl.genes_by_name(gene_name)

# Retrieve the start and end coordinates
start = gene[0].start
end = gene[0].end

print("Start:", start)
print("End:", end)'''
'''from pyensembl import Genome

# Specify the genome assembly version
genome = Genome("GRCh38")

# Retrieve gene start and end coordinates
gene_name = "CYP3A5"
gene = genome.genes_by_name(gene_name)[0]
start = gene.start
end = gene.end

# Print the coordinates
print("Start:", start)
print("End:", end)'''
'''from ucsc.api import Hub, Genome, Track, TrackSchema, Chromosome, Sequence

chromosomes = Chromosome.get(genome='hg38', track='knownGene')
trackk = Track.find('hg38','knownGene') 
genome = Genome.find('hg38')
track = genome.findTrack('knownGene')


chromosomeFragment = track.trackData(genome='hg38', chrom='7')
print(chromosomeFragment[0].__dict__)'''


'''from Bio import Entrez,SeqIO,Seq

def get_gene_info(gene_name):
    Entrez.email = 'marcshababy02@gmail.com'  # Provide your email address for Entrez API
    
    # Use Entrez to search for the gene
    handle = Entrez.esearch(db='gene', term=gene_name)
    record = Entrez.read(handle)
   

    # Get the gene ID from the search results
    gene_ids = record['IdList']
    gene_id=None
    for genes_id in gene_ids:
        fetch_handle = Entrez.efetch(db='gene', id=genes_id, retmode='xml')
        gene_records = Entrez.read(fetch_handle)
        if gene_records[0]['Entrezgene_gene']['Gene-ref']['Gene-ref_locus'] == gene_name:
            gene_id=genes_id
    fetch_handle.close()
    handle.close()
    Entrez.email = 'marcshababy02@gmail.com'
   
    

    handle.close()
    # Use Entrez to retrieve the gene record

    handle = Entrez.efetch(db='gene', id=gene_id, retmode='xml')
    record = Entrez.read(handle)
    
   
    # Extract the chromosome, start, end, and strand information
    chromosome = record[0]['Entrezgene_locus'][0]['Gene-commentary_products'][0]['Gene-commentary_accession']
    start = record[0]['Entrezgene_locus'][0]['Gene-commentary_seqs'][0]['Seq-loc_int']['Seq-interval']['Seq-interval_from']
    end = record[0]['Entrezgene_locus'][0]['Gene-commentary_seqs'][0]['Seq-loc_int']['Seq-interval']['Seq-interval_to']
    print(chromosome,start,end)
    handle.close()
    return chromosome, start, end

# Example usage
gene_name = 'CYP3A5'
get_gene_info(gene_name)
#chromosome, start, end = get_gene_info(gene_name)
#print(f"The gene {gene_name} is located on chromosome {chromosome}.")
#print(f"The start position is {start}, the end position is {end}")
def try_me():
    c=1
    e=3
    f="me"
    return c,e,f
chr, strt, end =try_me()
'''

import requests
from Bio import Seq
import numpy as np
import pandas as pd

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
    
    exon_loc=filter_gene['chromStarts'].split(',') #convert the value of string into a list
    exon_sizes=filter_gene['blockSizes'].split(',')
  

    return chrom,strand,start,end,exon_loc,exon_sizes


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
gene_name='CYP2C8'
chrom,strand,start,end,exon_loc,exon_sizes=get_gene(gene_name)
get_coding_strand('hg38',start,end,chrom,strand,exon_loc,exon_sizes,gene_name,'C:/Users/DELL/Desktop/cyp/test3')


'''
from Bio import Entrez
import requests
import json
import re
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
    print(CDS)
    match = re.search(r'\d+',CDS)
    if match:
        first_number = match.group()
        print(first_number)
    fetch_handle.close()
   
    return 1
    e=Entrez.read(fetch_handle)
    json_data=json.dumps
    with open( f"C:/Users/DELL/Desktop/Trynuc.json", "w") as file:
        
        json.dump(e, file)
    
    handle.close()
    def get_CDS(enst_id):
    ensembl_rest_url = f"https://rest.ensembl.org/lookup/id/{enst_id}?expand=1"
    response = requests.get(ensembl_rest_url, headers={"Content-Type": "application/json"})
    data=response.json()
    print(data)
# Example usage
transcript_id = 'NM_000770'  # Replace with your desired transcript ID
gne="CYP2C8"
fetch_cds_from_ncbi(transcript_id,gne)'''

