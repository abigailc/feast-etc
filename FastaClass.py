
#use as a module. import FastaClass
#this is a reusable definition of Fasta that allows on-the-fly editing of sequence ids and sequences.
#add shorten in? to comply with blast output.
#last edit abigailc@Artemis Sept 8 2016


class Fasta:
	def __init__(self, name):
		#all ids should be stripped and have ">" removed for reasons.
		#for now, sequences do not have any stripping applied
		self.name = name
		self.ids = []
		self.original_ids = []
		self.original_seqs = []
		self.seqs = []
	def gen_original_lists(self, fastaname):
		with open(fastaname) as fastafile:
			for line in fastafile:
				if ">" in line:
					#write the previous AA seq
					try:
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
					self.ids.append(newline)
					self.original_ids.append(newline)
				else:
					AAseq = AAseq+line
			#catch the last AAseq pass
			self.seqs.append(AAseq)
			self.original_seqs.append(AAseq)
		print("Initial sequence and ID lists created. Contains "+str(len(self.ids))+" sequences")
	def load_info_swap(self, info_file_in):
		#reads a file of form
		#   originalID
		#   changedID
		#and generates self.ids from that file.
		kid = "no"
		vid = "no"
		CTdict = {}
		with open (infofile) as old:
			for line in old:
				#first pass: gets key (original ID)
				#second pass: gets value (new ID)
				#if we have no info, get key
				if kid == "no":
					key = line.strip()
					kid = "yes"
					continue
				elif kid == "yes":
					#if we have key and value, record.
					if vid == "yes":
						CTdict[key]=value
						vid = "no"
						kid = "no"
						continue
					#if we have key but no value, get value.
					if vid == "no":
						value = line.strip()
						vid = "yes"
			#catch the final pass
			CTdict[key]=value
		for item in self.original_ids:
			index = self.original_ids.index(item)
			newid = CTdict[item]
			self.ids[index] = newestid
		#done
		#troubleshooting: do not preform this operation after any that change self.ids. this op must be done first, or in a seperate command.

	def gen_new_fasta(self, new_fasta_name):
		#this should print the changed seqids and changed AA sequences to file.
		newfasta = new_fasta_name
		# print(len(self.original_ids))
		# print(len(self.ids))
		# print(len(self.original_seqs))
		# print(len(self.seqs))
		with open (newfasta, "w") as new:
			for i in range(len(self.original_ids)):
				new.write(">"+self.ids[i])
				# print(i)		#
				#unclear if this needs a "\n" after it... check.#TODO
				new.write(self.seqs[i])
		print("Finished, your new fasta file is located at "+newfasta)
		#done

	def swap_in_newick(self, old_newick_name, new_file_name):
		#this replaces the tip names in a newick file. sometimes works on nexus files too, but I havent extensively tested it.
		newick = old_newick_name
		newnewick = new_file_name
		with open (newick) as old:
			with open (newnewick, "w") as new:
				for line in old:
					for item in self.original_ids:
						index = self.original_ids.index(item)
						line = line.replace(item, self.ids[index])
						new.write(line)
		print("finished, tip-replaced-newick file at: "+newnewick)
		#done

	def swap_in_nexus(self):
		print ("You didn't implement this yet. try using newick replace, it might work")
		pass
		#something
		#to-do, try nexus replace in the meantime, it should work
	def gen_info(self, info_file_name):
		#writes a file of form
		#   originalID
		#   changedID
		with open(info_file_name, "w") as inf:
			listlength = len(self.original_ids)
			if listlength != len(self.ids):
				print ("List lengths do not match! FATAL ERROR")
				print (self.original_ids)
				print (self.ids)
				raiseSystemExit
			for i in range(listlength):
				inf.write(self.original_ids[i])
				inf.write(self.ids[i]+"\n")
		print("Info file was generated. Named "+info_file_name)
		#done
