"""
In this chapter we’ll discuss in more detail the Bio.SeqIO module, which was briefly introduced in Chapter 2 and
also used in Chapter 4. This aims to provide a simple interface for working with assorted sequence file formats in a
uniform way. See also the Bio.SeqIO wiki page (http://biopython.org/wiki/SeqIO), and the built in documentation (also
online):

This deal with SeqRecord objects, they contain a Seq a Seq object (see Chapter 3) plus annotation like an identifier
and description.

Big Datasets can have prohibitive overhead, consider a low-level SimpleFastaParser or FastqGeneralIterator which
return just a tuple of strings for each record (see Section 5.6).

The workhorse function Bio.SeqIO.parse() is used to read in sequence data as SeqRecord objects.

This function expects two arguments:

    Path to dataset
    The first argument is a handle to read the data from, or a filename. A handle is typically a file opened for
    reading, but could be the output from a command line program, or data downloaded from the internet (see Section
    5.3). See Section 24.1 for more about handles.

    e.g. "fasta" "genbank"
    The second argument is a lower case string specifying sequence format
    See http://biopython.org/wiki/SeqIO for a full listing
    of supported formats. <---- There are a whole lot.

    There is an optional argument alphabet to specify the alphabet to be used. This is useful for file formats like FASTA
    where otherwise Bio.SeqIO will default to a generic alphabet.

The Bio.SeqIO.parse() function returns an iterator which gives SeqRecord objects.
 Iterators are typically used in a for loop as shown below.
"""
from Bio import SeqIO

'''
This will print out the id number, dna sequence and length of the record for every SeqRecord in the array.

for seq_record in SeqIO.parse("ls_orchid.fasta", "fasta"):
    print(seq_record.id)
    print(repr(seq_record.seq))
    print(len(seq_record))


for seq_record in SeqIO.parse("ls_orchid.gbk", "genbank"):
    print(seq_record.id)
    print(repr(seq_record.seq))
    print(len(seq_record))
'''

jackal = SeqIO.parse("ls_orchid.gbk", "genbank")

# Take off the first SeqRecord object to test on
pop = next(jackal)

# Hold the functions which are "not special" i.e. ones in everyday use to reference
# they are ['annotations', 'dbxrefs'(this one holds any database cross references:), 'description',
# 'features', 'format', 'id', 'letter_annotations', 'lower', 'name', 'reverse_complement', 'seq', 'translate', 'upper']
# Features could have their own wiki page, refer in future.
# for ref in pop.features:
#     print(ref)

# Prints SeqRecord as a string
# print(pop.format("fasta"))
print(pop.annotations["source"])

aret = []

for seqrec in SeqIO.parse("ls_orchid.gbk", "genbank"):
    aret.append(seqrec.annotations["source"])
# print(aret)


# or the above operation can be done with a list annotation
all_species = [
    seq_record.annotations["organism"] for seq_record in SeqIO.parse("ls_orchid.gbk", "genbank")
]
# print(all_species)

# So here we are using the function of enumerate which will tag each consecutive entry in the series of SeqRecord
# objects with an index (for index) and the contents of the SeqRecord as record (for index, record)
# So supposedly after this I can access the attributes of the SeqRecord object with the dot notation...
# Would not let me do that just iterating with the SeRecIter.next.attribute

"""
SAMPLE record.description = P.haynaldianum 5.8S rRNA gene and ITS1 and ITS2 DNA

A DNA barcode is a short piece of DNA sequence used for species determination and discovery. 
The internal transcribed spacer (ITS/ITS2) region has been proposed as the standard DNA barcode
 for fungi and seed plants and has been widely used in DNA barcoding analyses for other biological groups,
  for example algae, protists and animals. The ITS region consists of both ITS1 and ITS2 regions.
  From this: https://pubmed.ncbi.nlm.nih.gov/25187125/
  
  
for index, record in enumerate(SeqIO.parse("ls_orchid.gbk", "genbank")):
    print("index %i, Description = %s, length %i, with %i features"
          % (index, record.description, len(record.seq), len(record.features)))

"""
