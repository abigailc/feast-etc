# #!/usr/bin/python

# last edit abigailc@Actaeon on october 24 2016

### usage $ python Concatenate_2.py /Users/abigailc/Documents/Test/ -f "Asco1.fasta Asco2.fasta Asco3.fasta Asco4.fasta" -o AscoCC.fasta


#this program will open n aligned files, and concatenate them.
#sequences must either have exactly identicle tip names across files
#or be in shortened format " >Taxonomy|Info|Whatever|Species_name|##########" and use the -s flag

###########start##############
def master(fasta_list, output_name, species_mode = False):
    #generate your fasta objects and lists of identifiers
    list_spec_lists, list_id_lists, list_of_fasta_objects = make_fasta_info(fasta_list)
    #either correlate them my species_name or by seqid
    if species_mode is False:
        corr_ids = correlate_ids(list_id_lists)
    elif species_mode is True:
        corr_ids = correlate_ids(list_spec_lists)
    final = Concatenate(corr_ids, list_of_fasta_objects, output_name)
    print("Finished!")
    print("Your new concat file is located at: "+final)
    CCfas = Fasta(final)
    CCfas.gen_original_lists(final)
    print("There are: "+str(CCfas.number_seqs())+" sequences")
    print("And "+str(CCfas.number_of_sites())+" sites per sequence")

#this generates an object of class Fasta for each input .fasta file, and generates the species_lists and id lists.
def make_fasta_info(fasta_list):
    list_of_species_lists = []
    list_of_idlists = []
    list_of_fas = []
    for f in fasta_list:
        g = Fasta(f)
        g.gen_original_lists(f)
        g.gen_species_lists()
        species_list = g.species_list
        ids_list = g.ids
        list_of_idlists.append(ids_list)
        list_of_species_lists.append(species_list)
        list_of_fas.append(g)
    return list_of_species_lists, list_of_idlists, list_of_fas
# this will correlated IDs based on Genus_species or full id name across datasets (all of which need to be opened in class Fasta
def correlate_ids(list_of_id_lists):
    numlists = len(list_of_id_lists)
    output_list = []
    used = []
    for item in list_of_id_lists:
        # list1
        for ids in item:
            # id1 of list1
            name = ids
            if name in used:
                pass
            else:
                used.append(name)
                id_index_list = []
                id_index_list.append(ids)
                # check the index of that id in each list and append to "id_index_list"
                # if the id is not in that list, should append "NA"
                for eachlist in list_of_id_lists:
                    try:
                        index = eachlist.index(ids)
                    except ValueError:
                        index = "NA"
                    id_index_list.append(index)
                # add the result of scanning that id to overall output list.
                output_list.append(id_index_list)
    # output list looks like:
    # outputlist = [ ["Cat",1,2,3,4] , ["Dog", 2,1,13,14] ]
    print("Correlated sequences")
    print(output_list)
    return output_list

#does the actual concat and prints the output file.
def Concatenate(indexed_ids_list, fasta_class_list, output_name):

    # indexed_ids_list in form [ ["Cat", 1,2,3],["dog",2,1,2] ]
    # this part will create a new, concatenated .fasta file
    # requires indexed_ids_list, fasta_class_list, function that returns
    # number of sites. self.number_of_sites
    lenlist_final = []
    with open(output_name, "w") as new:
        for item in indexed_ids_list:
            #lenlist will collect lengths of each seq to verify at the end that each seq is the same length as each other (nothing was messed up)
            lenlist = []
            #write the first id
            new.write(">" + item[0].strip()+"\n")
            fas_num = 0
            allseq = ""
            #ensure that for each fas in fasta class list, all sequences are of the same length.
            for fas in fasta_class_list:
                fas_num += 1
                # fas_num keeps track of what number fasta we are on, which
                # correlates to the index of the index in indexed_ids_list
                search_index = item[fas_num]
                # search_index will be "2" if this is the first fasta (fas_num = 1) and item = ["cat", "2", "3", "23"]
                #represents "look at the second sequence of fasta (fas)"
                # if search_index is NA, generate a str of "-" that is n long
                if search_index == "NA":
                    ndash = fas.number_of_sites()
                    retreived_seq = ""
                    for i in range(int(ndash)):
                        retreived_seq = retreived_seq + ("-")
                else:
                    retreived_seq = fas.seqs[search_index]
                    # retreived_seq wil be something like "the 22nd sequence in
                    # object Fas's sequence list... " or "BLAHSEQUENCEDATA"
                #be sure to remove any \n that may have been saved up in there.
                retreived_seq = re.sub("\n", "", retreived_seq)
                #save how long the retreived seq is for future verification that all are same length
                lenlist.append(len(retreived_seq))
                #count is how many characters have been written to new concat file
                count = 0
                #retreived seq is like "Cat_gene_1"
                #allseq is like "Cat_gene_1"+"Cat_gene_2"+"Cat_gene_3"
                allseq = allseq + retreived_seq
            #add length total for verification purposes
            lenlist_final.append(len(allseq))
            newseq = ""
            #print in 80 character lines so it looks nice / fasta standard
            for letter in allseq:
                if count > 79:
                    count = 0
                    newseq = newseq + ("\n")
                newseq = newseq + letter
                count += 1
            #write the sequence to your new file.
            new.write(newseq.strip()+"\n")
        #verification bit
        for length in lenlist_final:
            if length == lenlist_final[0]:
                pass
            else:
                print("ERROR your concat sequences are not all of the same length something has gone horribly wrong! aborting.")
                raise SystemExit
    return output_name


#necessary class
class Fasta:
    def __init__(self, name="whatever"):
        # all ids should be stripped and have ">" removed for reasons.
        # for now, sequences do not have any stripping applied
        self.name = name
        self.ids = []
        self.original_ids = []
        self.original_seqs = []
        self.seqs = []
        self.species_list = []

    def gen_original_lists(self, fastaname):
        try:
            with open(fastaname) as fastafile:
                for line in fastafile:
                    if "\n" == line:
                        pass
                    if ">" in line:
                        # write the previous AA seq
                        try:
                            AAseq = AAseq.strip()
                            self.seqs.append(AAseq)
                            self.original_seqs.append(AAseq)
                        except:
                            pass
                            # initialize a new AAseq
                        AAseq = ""
                        # format the seqID
                        newline = line.strip()
                        newline = newline.strip(">")
                        # write the seqID
                        self.ids.append(newline)
                        self.original_ids.append(newline)
                    else:
                        AAseq = AAseq + line
                AAseq=AAseq.strip()
                # catch the last AAseq pass
                self.seqs.append(AAseq)
                self.original_seqs.append(AAseq)
            print("Initial sequence and ID lists created for "+self.name+". Contains " +
                  str(len(self.ids)) + " sequences")
        except UnboundLocalError:
            print("probably this file :" + fastaname +
                  " has nothing in it. skipping.")
            pass
        except IOError:
            print("no file named: " + fastaname +
                  " exists... creating a blank file")
            with open(fastaname, "w") as new:
                pass
            print("hopefully you intended that!")
    def number_seqs(self):
        a = len(self.ids)
        return a
    def number_of_sites(self, num = 0):
        testseq = self.original_seqs[num]
        testseq = re.sub("\n", "", testseq)
        #print(testseq)
        #print (len(self.original_seqs[num]))
        #print (len(testseq))
        return len(testseq)
    def gen_species_lists(self):
        for item in self.ids:
                            # item will be "Nostoc_punctiforme_PCC_73102|gi#|186468349" or "Blah|Rank|Nostoc_punctiforme_PCC_73102|gi#|186468349"
                            # for now, ignores anything that isn't Genus_species.
                            # for example, ignores strain, ., things with an extra
                            # word, etc.
            taxon = re.sub("([^_]*)([A-Z][a-z]*_[a-z]*)(.*)", "\\2", item)
            if "#" in taxon:
                print ("TAXON error in gen_species_lists():" + taxon)
            self.species_list.append(taxon)
        return self.species_list

# parser

if __name__ == "__main__":

    print("Running in terminal")  
    import sys
    import argparse
    import os
    import re
    parser = argparse.ArgumentParser(description="All")
    parser.add_argument("directory", nargs='?', default=os.getcwd(), type=str, help="type name of directory to run in (where .nex resides)")
    parser.add_argument("-f", "--fasta", action = "store", default = False, help="give each pre-aligned fasta file to concatenate within quotes, seperated by spaces. eg \"first.fasta second.fasta third.fasta\" ")
    parser.add_argument("-o", "--output", action = "store", default = "Concat.fasta", help="provide a name for your output")
    parser.add_argument("-s", "--species_mode", action = "store_true", default = False, help="use species_name, not entire seqID")

    args = parser.parse_args()
    #change dir if given
    print("hopefully your seqIDs (or species_names with -s) are identical across each file and have no repeats.")
    try:
        os.chdir(args.directory)
    except:
        print ("didn't change dir")
    #run the thing
    if args.fasta is False:
        print("Please give a list of .fasta files!")
    fastas = args.fasta.split()
    if len(fastas) is 1:
        print("Need more than one file to concatenate???")

    if args.species_mode is True:
        master(fastas, args.output, True)
    else:
        master(fastas, args.output)
    print("Well, it ought to be done now!")
