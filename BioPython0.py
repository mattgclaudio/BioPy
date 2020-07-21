from Bio.Seq import Seq, MutableSeq
from Bio import SeqIO
from Bio.Alphabet import IUPAC
from Bio.Alphabet import generic_dna
from Bio.Alphabet import *
from Bio.Data import CodonTable
import pprint as pprint

# Transcription and Translation
coding_strand = Seq("ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG", IUPAC.unambiguous_dna)
template_dna = coding_strand.reverse_complement()
messenger_rna = coding_strand.transcribe()
# print(messenger_rna)
# print(messenger_rna.translate(table=2, to_stop=True))

# If the complete CDS i.e full set of encoded dna no stop codons has a start and stop all clean etc
# then some bacteria might have stop codons which differ from mammalian samples
# so specify table= X e.g. gene.translate(table="Bacterial" to_stop=True) and

gene = Seq("GTGAAAAAGATGCAATCTATCGTACTCGCACTTTCCCTGGTTCTGGTCGCTCCCATGGCA" +
           "GCACAGGCTGCGGAAATTACGTTAGTCCCGTCAGTAAAATTACAGATAGGCGATCGTGAT" +
           "TATGAATGGCGAGGCAATCGCTGGCACCTACACGGACCGCCGCCACCGCCGCGCCACCAT" +
           "AAGAAAGCTCCTCATGATCATCACGGCGGTCATGGTCCAGGCAAACATCACCGCTAA",
           generic_dna)

# print(gene.translate(table="Bacterial", to_stop=True))

# In the bacterial genetic code GTG is a valid start codon, and while it does normally encode Valine,
# if used as a start codon it should be translated as methionine.
# This happens if you tell Biopython your sequence is a complete CDS:

gene.translate(table="Bacterial", cds=True)

# In addition to telling Biopython to translate an alternative start codon as methionine,
# using this option also makes sure your sequence really is a valid CDS (you’ll get an exception if not).


# What is this notation below?
# my_seq = Seq("GATCG", IUPAC.unambiguous_dna)
# for index, letter in enumerate(my_seq):
#    print("%i %s" % (index, letter))

standard_table = CodonTable.unambiguous_dna_by_name["Standard"]
# or by_id[1]

mito_table = CodonTable.unambiguous_dna_by_name["Vertebrate Mitochondrial"]
# or by_id[2]
# standard_table.stop_codons
# print(standard_table.stop_codons)

# compare sequences as strings

j1 = Seq("ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG", IUPAC.unambiguous_dna)
j2 = Seq("ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATTG", IUPAC.unambiguous_dna)

str(j1) == str(j2)

# You cant index and change regular sequences, but if you make them mutable you can do basically anything you want

j3 = j2.tomutable()
j3[0::1]
j3.pop()

# or do this
mutable_seq = MutableSeq("GCCATTGTAATGGGCCGCTGAAAGGGTGCCCGA", IUPAC.unambiguous_dna)

# print(mutable_seq.remove("G"))

