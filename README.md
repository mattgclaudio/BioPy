# BioPy: Exploration of bioinformatics libraries in Python, now with Tensorflow.

The Biopython module provides a function set for parsing DNA, RNA, Protein sequences from file, and performing a wide range of manipulations on a given sequence.
> [Biopython](https://biopython.org/) is a set of freely available tools for biological computation written in Python by an international team of developers.
Biopy represents a sequence as a `SeqRecord` [object](https://biopython.org/wiki/SeqRecord). Relevant properties:
>`[..., 'annotations', 'dbxrefs', 'description', 'features', 'format', 'id', 'letter_annotations', 'name', 'seq']`

BioPy has a number of other modules
- Bio.AlignIO - alignment input/output 
- Bio.PopGen - population genetics
- Bio.PDB - structural bioinformatics
- Biopython’s BioSQL interface
For reading and writing purely DNA/RNA sequences, [SeqIO](https://biopython.org/wiki/SeqIO) is the way to go.



> In bioinformatics, a sequence alignment is a way of arranging the sequences of DNA, RNA, or protein
>to identify regions of similarity that may be a consequence of functional, structural, or evolutionary relationships between the sequences.
