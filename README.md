# feast-etc
place for some of my scripts to go

ones that work right now:

FEAST Fasta Editing and Subsampling Tool

        It is useful for modifying and subsampling .fasta files based on seqid and taxonomic information.
        ___________________________
        Commands:
        Simplify/SimplifyKeep: Turns NCBI Blast output into nicer .fasta format. Also works with "6 sseqid stitle sseq"      
                               blast+ outfmt.
        Remove: removes sequences from .fasta if the seqid contains X.
        PFAM: removes sequences from .fasta if they do not contain pfam conserved domain X.
        AppendRanks: adds taxonomic data from NCBI Taxonomy of rank X (phylum, class, etc) to seqID.
        Separate: splits .fasta file into two files. one has all seqs where seqIDs contains X, other has the remaining seqs.
        Extract: creates a new .fasta file containing only sequences whose seqID contains X.
        SubSample: creates a new fasta file containing N sequences per unique str of rank X. To be used after AppendRanks.
        Summarize: summarizes what taxa are in a .fasta file by rank X. eg, 5 Alphas, 2 Deltas.
        Compare: Compares taxa across two files. Useful for finding errors, or seeing what taxa is missing a gene or two.
        MultiDataSubSample: As in SubSample, but considers multiple datasets and optimizes for shared species (for future 
                            concatenation purposes).
        Merge: creates a new .fasta file containing all sequences from all X input files.
FISH Fasta Id Shortening Helper

        It is used to modify, fix, and/or shorten sequence IDs in fasta files and trees.
        ________________
        Commands:
        Length: shortens seqIDs to a max length
        Duplicates: adds numbers to identical seqIDs to differentiate them
        FixID: removes weird characters like ';_+= from seqIDs
        FixAA: removes any non-standard characters from AA sequence
        Common: common name-shortening tricks to reduce length eg Bacteria->Bac
        Manual: allows manual name-shortening
        TenChars: makes seqids ten characters long and referencable. >HomoSapi01
        Bayes: makes seqids that won't error in MrBayes
        Piece: given appended bar-seperated taxonomy, shorten by depth
        WriteFasta: writes your changes to a new fasta
        WriteNewick: takes a newick, replaces original names (from fasta or info) in it with your changes.
        WriteInfo: saves original, changed pair to be read later (if you need to change back ever, do it)
        
Useful things:

        Fasta class.
                unfinished, but good for coding small things with. keeps original seqIDs and seqs, allows modifications, can 
                write new fasta.
        MakeSpeciesTrees.py
                Currently working on, will probably be somewhat specific to the engaging cluster @ mit.

