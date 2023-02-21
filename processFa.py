import gzip
from Bio import SeqIO

def parse_fasta(fasta_file):
    """_reads in a zipped fasta file and returns a dictioary of record.id:record.seq
    Input: 
        file: path to gzipped fasta file
    Returns:
        dictionary: a dictioary of seqRecords objects, indexed by record.id
    """
    seq_dict = {}
    with gzip.open(fasta_file, "rt") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            seq_dict[record.id] = record
    return seq_dict
