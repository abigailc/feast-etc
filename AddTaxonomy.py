#!/usr/bin/python
#last edit abigailc@ACTAEON
#purpose: adds taxonomy information to .fasta file sequence IDs
#purpose2: smaller and more customizable than implementation within FISH_2.py

###CLASS
class Fasta:
    def __init__(self, name):
        #this is the to-be-modified version of sequence IDs and sequence-Data
        # ALWAYS keep IDS and SEQS the same length. id[1] should ALWAYS correspond to seq[1].
        self.name = name
        self.ids = []
        self.seqs = []
        # these are the original SEQids and Sequences. They should never be modified after generation in gen_original_lists or blast_to_fasta
        self.original_ids = []
        self.original_seqs = []
        #obsolete
        self.species_names = []
        #list of gi / accession numbers
        self.numbers = []
        #list of taxid numbers
        self.taxid = []
        #list of dictionaries of taxonomy rank information from NCBI
        self.taxonomy = []
    #loads your file as IDS and SEQUENCES
    def gen_original_lists(self):
        fastaname = self.name
        with open(fastaname) as fastafile:
            for line in fastafile:
                if "\n" == line:
                    pass
                if ">" in line:
                    #write the previous AA seq
                    try:
                        AAseq=AAseq.strip()
                        self.seqs.append(AAseq)
                        self.original_seqs.append(AAseq)
                    except:
                        pass
                        #initialize a new AAseq
                    AAseq = ""
                    #format the seqID
                    newline = line.strip()
                    newline = line.strip(">")
                    #write the seqID
                    self.ids.append(newline.strip())
                    self.original_ids.append(newline.strip())
                else:
                    AAseq = AAseq+line
                    AAseq=AAseq.strip()
            #catch the last AAseq pass
            self.seqs.append(AAseq)
            self.original_seqs.append(AAseq)
        print("Initial sequence and ID lists created. Contains "+str(len(self.ids))+" sequences")
    #gets number (gi num or acc. num)
    def gen_numbers(self):
        print("Finding identification numbers...")
        for item in self.ids:
            number = re.sub("(.*)(\|)(.*)","\\3", item)
            self.numbers.append(number)
    #gets taxID
    def SetTaxID(self):
        print("acquiring taxids")
        self.taxid = []
        for item in self.numbers:
            #we are erroring and i don't know why.
 
            GItoTAXID = "xmllint --xpath '/GBSet/GBSeq/GBSeq_feature-table/GBFeature/GBFeature_quals/GBQualifier[GBQualifier_name=\"db_xref\"]/GBQualifier_value/text()' \"http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&id="+item+"&retmode=xml\""
            
            try:
                futuretaxid = subprocess.check_output(GItoTAXID, shell=True)
            except:
                print("Error with "+item)
                futuretaxid = "NA"
            taxid = re.sub("(taxon:)([0-9]*)(.*)", "\\2", futuretaxid)
            self.taxid.append(taxid)
    def GetTaxonomy(self):
        print("getting taxonomic information... this might take a while")
        total = len(self.taxid)
        current = 0
        self.taxonomy = []
        if self.taxid == []:
            print("You need to generate taxids first.. lets try")
            self.SetTaxID()
            total = len(self.taxid)
        for item in self.taxid:
            current += 1
            done = float(current)/float(total)#will be, like, 0.2
            percent = done*100 #will be, like, 20
            if percent.is_integer() is True:
                if int(percent)%10 == 0:
                    print(str(int(percent))+" percent done")
            taxid = item
            ranklist = "superkingdom kingdom phylum class order family genus"
            ranklist = ranklist.split()
            taxdict = {}
            for r in ranklist:
                TAXIDtoRANKNAME = "xmllint --xpath '/TaxaSet/Taxon/LineageEx/Taxon[Rank=\""+r+"\"]/ScientificName/text()' \"http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=taxonomy&id=" + taxid + "\""
                try:
                    rankname = subprocess.check_output(TAXIDtoRANKNAME, shell=True)
                except:
                    rankname = "NA"
                    rankname = re.sub(" ", "_", rankname)
                taxdict[r]=rankname
            self.taxonomy.append(taxdict)
          
    #this should modify the seqIDS from Species_name|ACC##### to Euk|Meta|Chord|Species_name|ACC#####
    def AppendTaxonomy(self, ranklist = "NA"):
        for item in self.ids:
            index = self.ids.index(item)
            rankdict = self.taxonomy[index]
            if ranklist == "NA":
                newitem = rankdict["superkingdom"]+"|"+rankdict["kingdom"]+"|"+rankdict["phylum"]+"|"+rankdict["class"]+"|"+rankdict["order"]+"|"+rankdict["family"]+"|"+rankdict["genus"]+"|"+item
                self.ids[index] = newitem
            else:
                a = ""
                for rank in ranklist:
                    a = a + rankdict[rank]+"|"
                a = a + item
                self.ids[index] = a
    def gen_new_fasta(self, new_fasta_name):
        #this should print the changed seqids and changed AA sequences to file.
        newfasta = new_fasta_name
        # print(len(self.original_ids))
        # print(len(self.ids))
        # print(len(self.original_seqs))
        # print(len(self.seqs))
        with open (newfasta, "w") as new:
            for i in range(len(self.original_ids)):
                new.write(">"+self.ids[i].strip()+"\n")
                # print(i)      #
                #unclear if this needs a "\n" after it... check.#TODO
                                #print(self.seqs)
                                #print(type(self.seqs[i]))
                new.write(self.seqs[i]+"\n")
        print("Finished, your new fasta file is located at "+newfasta)
        #done
    def gen_species_lists(self):
        for item in self.ids:
            taxon = re.sub("([^_]*)([A-Z][a-z]*_[a-z]*)(.*)", "\\2", item)
            if "#" in taxon:
                print ("TAXON error in gen_species_lists():" + taxon)
            self.species_names.append(taxon)
    def PrintTaxonInfo(self):
        self.gen_species_lists()
        try:
            a=self.name.split(".")
            taxonfile = a[0]+"_taxinfo.txt"
        except:
            taxonfile = self.name+"_taxinfo.txt"
        with open(taxonfile, "w") as new:
            for item in self.ids:
                index = self.ids.index(item)
                seqid = item
                taxid = self.taxid[index]
                sciname = self.species_names[index]
                ranks = "superkingdom kingdom phylum class order family genus"
                ranklist = ranks.split()
                taxonomy = ""
                for r in ranklist:
                    try:
                        tax = self.taxonomy[index][r]
                    except:
                        pass
                    taxonomy = taxonomy+";"+tax
                new.write(item+"\t"+sciname+"\t"+taxid+"\t"+taxonomy)
                    
            

#parser

if __name__ == "__main__":

    print("Running in terminal")
    import subprocess
    import sys
    import argparse
    import os
    import re
    parser = argparse.ArgumentParser(description="All")
    parser.add_argument("directory", nargs='?', default=os.getcwd(), type=str, help="type name of directory to run in (where .nex resides)")
    parser.add_argument("-r", "--ranks", action = "store", default = False, help="give specific ranks to append, if you dont want them all. eg -r phylum")
    parser.add_argument("-f", "--fasta", action = "store", default = False, help="give a .fasta file")
    parser.add_argument("-i", "--info", action = "store_true", default = False, help="toggle saves taxonomy information to new file")


    
    args = parser.parse_args()
    #change dir if given
    try:
        os.chdir(args.directory)
    except:
        print ("didn't change dir")
    #run the thing
    MyFasta = Fasta(args.fasta)
    MyFasta.gen_original_lists()
    MyFasta.gen_numbers()
    MyFasta.SetTaxID()
    MyFasta.GetTaxonomy()
    if " " in args.ranks:
        ranklist = args.ranks.split()
    else:
        ranklist = [args.ranks]
    if args.ranks:
         MyFasta.AppendTaxonomy(ranklist)
    else:
         MyFasta.AppendTaxonomy()
    #get a name for the new fasta
    try:
        base, ext = args.fasta.split(".")
        newfastaname = base+"_Taxo.fasta"
    except:
        newfastaname = args.fasta+"_Taxo.fasta"
    MyFasta.gen_new_fasta(newfastaname)
    if args.info is True:
        MyFasta.PrintTaxonInfo()
    print("done, your new file it at: "+newfastaname)
