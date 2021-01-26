# BioPy: Exploration of bioinformatics libraries in Python, now with Tensorflow.

The Biopython module provides a function set for parsing DNA, RNA, Protein sequences from file, and performing a wide range of manipulations on a given sequence.
> [Biopython](https://biopython.org/) is a set of freely available tools for biological computation written in Python by an international team of developers.

Biopy represents a sequence as a `SeqRecord` [object](https://biopython.org/wiki/SeqRecord).
Relevant properties:
>`[..., 'annotations', 'dbxrefs', 'description', 'features', 'format', 'id', 'letter_annotations', 'name', 'seq']`

BioPy has a number of modules
- Bio.AlignIO - alignment input/output 
- Bio.PopGen - population genetics
- Bio.PDB - structural bioinformatics
- Biopython’s BioSQL interface

For reading and writing purely DNA/RNA sequences, [SeqIO](https://biopython.org/wiki/SeqIO) is only class needed however. 

So I pulled the virus sequences/vaccines ( not gathered with the interfaces BioPy can query, sadly) and wrote some Python to reformat them in a manner which 
a Tensorflow model could interpret. 
`[batch, timestep, feature]` 

The stable model does work, although I switched to using the catagorical crossentropy loss function (full disclosure the math explanations for this, adam, etc.
are beyond me, I just pick different ones for effect) and it starts at  ~4 and is then flirting with .6 after running just the three virus/vaccine pairs
I initially gathered. This seems suspicious. 

>  _In bioinformatics, a sequence alignment is a way of arranging the sequences of DNA, RNA, or protein
> to identify regions of similarity that may be a consequence of functional, structural, or evolutionary relationships between the sequences._
