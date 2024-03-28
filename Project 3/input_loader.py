from Bio import SeqIO
def load_alphabet(file):
    f = open(file, "r")
    alphabet=f.read()
    alphabet = alphabet.replace(' ','')
    alphabet=alphabet.strip()
    return alphabet

# loading sequences
def load_sequences(filename):
    sequences= list(SeqIO.parse(filename, "fasta"))
    seq=sequences[0].seq

    global sequence_len
    sequence_len=len(seq)
    #return seq[:1000]+'$'
    return seq