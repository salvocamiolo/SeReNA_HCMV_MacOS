# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'genomeEntry.ui'
#
# Created by: PyQt5 UI code generator 5.9.2
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5.QtCore import QObject, QThread, pyqtSignal
import os
from Bio import SeqIO

from Bio.Seq import Seq
from Bio import pairwise2
import time
import os
from os import listdir
from os.path import isfile, join
import time





# Step 1: Create a worker class
class Worker(QObject):
	finished = pyqtSignal()
	progress = pyqtSignal(str)

	def __init__(self,ii,g,sampleID, sampleType, sampleDate, sampleCompartment):
		super(QObject, self).__init__()
		self.installationDirectory = ii
		self.genomeFileName = g
		self.sampleID = sampleID
		self.sampleType = sampleType
		self.sampleDate = sampleDate
		self.sampleCompartment = sampleCompartment


	def reverseComplement(self,sequence):
		outSequence = ""
		accepted = ["A","C","T","G","a","c","t","g","N","n"]
		for a in range(len(sequence)-1,-1,-1):
			if sequence[a]=="A" or sequence[a]=="a":
				outSequence+="T"
			elif sequence[a]=="T" or sequence[a]=="t":
				outSequence+="A"
			elif sequence[a]=="G" or sequence[a]=="g":
				outSequence+="C"
			elif sequence[a]=="C" or sequence[a]=="c":
				outSequence+="G"
			else:
				outSequence+="N"
			

		return outSequence

	def run(self):
		"""Long-running task."""
		genomeFileNameList = []
		genomeFileNameMetadata = []
		if ".bulk" in self.genomeFileName:
			infile = open(self.genomeFileName)
			while True:
				line = infile.readline().rstrip()
				if not line:
					break
				fields =  line.split("\t")
				genomeFileNameList.append(fields[0])
				genomeFileNameMetadata.append((fields[1],fields[2],fields[3],fields[4]))
			infile.close()
		
		else:
			genomeFileNameList.append(self.genomeFileName)
			genomeFileNameMetadata.append((self.sampleID,self.sampleType,self.sampleDate,self.sampleCompartment))


		for c in range(len(genomeFileNameList)):

			prot2map = []
			tempFolder = genomeFileNameMetadata[c][0]
			os.system("mkdir "+tempFolder)

			prefix = genomeFileNameMetadata[c][1]+"-"+genomeFileNameMetadata[c][0]+"-"+genomeFileNameMetadata[c][2]+"-"+genomeFileNameMetadata[c][3]
			cdsFasta = open(tempFolder+"/"+prefix+"_cds.fasta","w")
			proteinFasta = open(tempFolder+"/"+prefix+"_proteins.fasta","w")
			logFile = open(tempFolder+"/"+prefix+"_annotationLog.txt","w")



			for seq_record in SeqIO.parse(genomeFileNameList[c],"fasta"):
				genomeSeq = str(seq_record.seq)
				assemblyName = str(seq_record.id)
		
			onlyfiles = [f for f in listdir(self.installationDirectory+"data/proteinDB/") if isfile(join(self.installationDirectory+"data/proteinDB/", f))]
			for f in onlyfiles:
				if f[-13:] == "_models.fasta":
					prot2map.append(f)

			for f in prot2map:
			
				numCodonRefines = 0
				#Find the best model for the gene *********************************************************
				locus = f.replace("_models.fasta","")
				self.progress.emit("Searching CDS for "+locus+"....")
				protSeqs = SeqIO.to_dict(SeqIO.parse(self.installationDirectory+"data/proteinDB/"+f,"fasta"))
				#print "Choosing best match for protein",f
							

				os.system(self.installationDirectory+"src/conda/bin/makeblastdb -dbtype nucl -in "+genomeFileNameList[c])
				os.system(self.installationDirectory+"src/conda/bin/tblastn -query "+self.installationDirectory +
							"data/proteinDB/"+f+" -db "+genomeFileNameList[c]+" -outfmt 6 -max_intron_length 350000 | sort -k 12rn,12rn  > "+tempFolder+"/preOutputBlast.txt")
				time.sleep(1)
				os.system("head -20 "+tempFolder+"/preOutputBlast.txt > "+tempFolder+"/preOutputBlast2.txt")
				time.sleep(1)
				os.system("sort -k 3rn,3rn -k12rn,12rn "+tempFolder+"/preOutputBlast2.txt > "+tempFolder+"/outputBlast.txt")
				os.system("rm -f "+tempFolder+"/preOutputBlast* ")
				blastFile = open(tempFolder+"/outputBlast.txt")
				bestCoverage = 0
				blastLine = blastFile.readline().rstrip()
				blastField = blastLine.split("\t")
				alignmentLength = int(blastField[7]) - int(blastField[6])
				bestCoverage = alignmentLength
				bestProt = blastField[0]

				for hit in range(5):
					if float(blastField[2]) > 99.8:
						break
					blastLine = blastFile.readline().rstrip()
					if not blastLine:
						break
					blastField = blastLine.split("\t")
					alignmentLength = int(blastField[7]) - int(blastField[6])
					if alignmentLength > bestCoverage:
						bestCoverage = alignmentLength
						bestProt = blastField[0]


				sequence = str(protSeqs[bestProt].seq)
				tempFasta = open(tempFolder+"/tempFasta.fasta","w")
				tempFasta.write(">"+locus+"\n"+str(protSeqs[bestProt].seq)+"\n")
				tempFasta.close()

				#Run Exonerate on the best model *********************************************************


				os.system(self.installationDirectory+"src/conda/bin/exonerate --model protein2genome "+tempFolder+"/tempFasta.fasta "+genomeFileNameList[c]+" --showtargetgff -s 0 -n 1 --forcegtag --minintron 35 --maxintron 10000  > "+tempFolder+"/outputExonerate")



				#Check Exonerate output *****************************************************************
				#Check the proteins gave a match in the target genome
				exResult = open(tempFolder+"/outputExonerate")
				line = exResult.readline().rstrip()
				while not "Query range:" in line:
					line = exResult.readline().rstrip()
					if line is None:
						print("WARNING no protein found for ", locus)
						logFile.write("Missing gene: locus "+locus+" did not provide any alignment")
				exResult.close()


				#Reconstruct Exons  **********************************************************************
				exResult = open(tempFolder+"/outputExonerate")
				while not line == "# --- START OF GFF DUMP ---":
					line = exResult.readline().rstrip()
					if line is None:
						#print "WARNING no protein found for ", locus
						logFile.write("Missing gene: locus "+locus+" did not provide any alignment")
				for a in range(10):
					line = exResult.readline().rstrip()

				gene = {}
				exon = {}
				warnings = []
				#Collect the exonerate output
				while not line == "# --- END OF GFF DUMP ---":
					line = exResult.readline().rstrip()
					if line is None:
						#print "WARNING no protein found for ", locus
						logFile.write("Missing gene: locus "+locus+" did not provide any alignment")
					fields = line.split("\t")
					if not line == "# --- END OF GFF DUMP ---":
						if fields[2]=="gene":
							if not locus in gene:
								gene[locus] = (fields[3],fields[4],fields[6])
							else:
								if (int(fields[4]) - int(fields[3])) > (int(int(gene[locus][1])) - int(int(gene[locus][0]))):
									gene[locus] = (fields[3],fields[4],fields[6])

						if fields[2] == "exon":
							if not locus in exon:
								exon[locus] = []
							exon[locus].append((fields[3],fields[4],fields[6]))
							#if "frameshifts" in fields[8]:
							#	warnings.append("Frameshifts in exon "+str(fields[3])+" "+str(fields[4])+" "+fields[8])

				newList = sorted(exon[locus], key=lambda x: x[1])
				exon[locus] = newList


				#Reconstruct CDS  ***********************************************************
				cdsSeq = ""
				if exon[locus][0][2]=="+":   #************* Positive strand
					for item in exon[locus]:
						cdsSeq+=genomeSeq[int(item[0])-1:int(item[1])]
					cdsSeq += genomeSeq[int(item[1]):int(item[1])+3]

				else: # *********************** Negative Strand
					for a in range(len(exon[locus])-1,-1,-1):
						cdsSeq+= self.reverseComplement(genomeSeq[ int(exon[locus][a][0])-1: int(exon[locus][a][1]) ])
					cdsSeq += self.reverseComplement(genomeSeq[int(exon[locus][a][0])-4 :int(exon[locus][a][0])-1])


				notes = ""

				# Check CDS integrity *********************************************************
				foundStartCodon = True
				foundStopCodon = True
				if not  cdsSeq[:3]=="ATG" or not (cdsSeq[-3:]=="TGA" or cdsSeq[-3:]=="TAA" or cdsSeq[-3:]=="TAG" ):

					notes = "\n"+locus+"\n"
					notes += "Either the start or the stop codon was not found. Searching nearby....\n"




					# Look for ATG at the beginning of the sequence or closely ********************
					if exon[locus][0][2]=="+":   #************* Positive strand
						if not cdsSeq[:3]=="ATG":
							foundStartCodon = False
							print("- Looking for start codon upstream....\n")
							for a in range(len(sequence)-int(len(cdsSeq)/3)+30):
								newStart = genomeSeq[int(exon[locus][0][0])-1-a*3-3:int(exon[locus][0][0])-1-a*3]
								if newStart == "ATG":
									notes += "Valid start codon found upstream!\n"
									exon[locus][0]=(int(exon[locus][0][0])-a*3-3,exon[locus][0][1],exon[locus][0][2])
									gene[locus] = (int(exon[locus][0][0])-a*3-3,int(gene[locus][1]),gene[locus][2])
									cdsSeq = ""
									for item in exon[locus]:
										cdsSeq+=genomeSeq[int(item[0])-1:int(item[1])]
									foundStartCodon = True
									notes += "- Start codon refined  "+str(a)+" codons upstream\n"
									numCodonRefines = a
									break
								if newStart == "TGA" or newStart == "TAA" or newStart=="TAG":
									foundStartCodon = False
									notes += "- Found a stop codon while searching for start codon upstream!\nStart codon could not be found upstream\n"
									break
							#If the new start codon was not found in the region upstream then the downstream region is searched
							if foundStartCodon == False:
								notes += "- Looking for start codon downstream....\n"
								for a in range(len(sequence)-int(len(cdsSeq)/3)+30):
									newStart = genomeSeq[int(exon[locus][0][0])-1+a*3:int(exon[locus][0][0])-1+a*3+3]
									if newStart == "ATG":
										notes += "Valid start codon found downstream!\n"
										exon[locus][0]=(int(exon[locus][0][0])+a*3,exon[locus][0][1],exon[locus][0][2])
										gene[locus] = (int(exon[locus][0][0])+a*3,int(gene[locus][1]),gene[locus][2])
										cdsSeq = ""
										for item in exon[locus]:
											cdsSeq+=genomeSeq[int(item[0])-1:int(item[1])]
										foundStartCodon = True
										notes += "- Start codon refined  "+str(a)+" codons downstream\n"
										numCodonRefines = a
										break
									if newStart == "TGA" or newStart == "TAA" or newStart=="TAG":
										foundStartCodon = False
										notes += "- Found a stop codon while searching for start codon downstream!\nStart codon could not be found downstream\n"
										break
								if foundStartCodon == False:
									notes += "Start codon could not be found at this stage\n"

					else: # *********************** Negative Strand
						if not cdsSeq[:3]=="ATG" or not (cdsSeq[:3] =="TTG" and locus=="RL6"): #RL6 start with alternative start codon
							foundStartCodon = False
							notes += "- Looking for start codon upstream....\n"
							for a in range(len(sequence)-int(int(len(cdsSeq)/3))+30):
								#print "New start codons"
								newStart = self.reverseComplement(genomeSeq[int(exon[locus][-1][1])+a*3:int(exon[locus][-1][1])+a*3+3])
								#print newStart
								if newStart == "ATG":
									notes += "Valid start codon found upstream!\n"
									exon[locus][-1]=(int(exon[locus][-1][0]),int(exon[locus][-1][1])+a*3+3,exon[locus][-1][2])
									gene[locus] = (int(gene[locus][0]),int(exon[locus][-1][1])+a*3+3,gene[locus][2])
									cdsSeq = ""
									for a1 in range(len(exon[locus])-1,-1,-1):
										cdsSeq+=self.reverseComplement(genomeSeq[ int(exon[locus][a1][0])-1: int(exon[locus][a1][1]) ])
									foundStartCodon = True
									notes += "- Start codon refined  "+str(a)+" codons upstream\n"
									numCodonRefines = a
									break
								if newStart == "TGA" or newStart == "TAA" or newStart=="TAG":
									foundStartCodon = False
									notes += "- Found a stop codon while searching for start codon upstream!\nStart codon could not be found upstream\n"
									break
						#If the new start codon was not found in the region upstream then the downstream region is searched
							if foundStartCodon == False:
								notes += "- Looking for start codon downstream....\n"

								for a in range(len(sequence)-int(len(cdsSeq)/3)+30):
									#print "New start codons"
									newStart = self.reverseComplement(genomeSeq[int(exon[locus][-1][1])-a*3-3:int(exon[locus][-1][1])-a*3])

									#print newStart
									if newStart == "ATG":
										notes += "Valid start codon found downstream!\n"
										exon[locus][-1]=(int(exon[locus][-1][0]),int(exon[locus][-1][1])-a*3,exon[locus][-1][2])
										gene[locus] = (int(gene[locus][0]),int(exon[locus][-1][1])-a*3,gene[locus][2])
										cdsSeq = ""
										for a1 in range(len(exon[locus])-1,-1,-1):
											cdsSeq+=self.reverseComplement(genomeSeq[ int(exon[locus][a1][0])-1: int(exon[locus][a1][1]) ])
										foundStartCodon = True
										notes += "- Start codon refined  "+str(a)+" codons downstream\n"
										numCodonRefines = a
										break
								if newStart == "TGA" or newStart == "TAA" or newStart=="TAG":
									foundStartCodon = False
									notes += "- Found a stop codon while searching for start codon downstream! \nStart codon could not be found downstream\n"
									break

								if foundStartCodon == False:
									notes += "Start codon could not be found at this stage\n"




					# Look for Stop codon at the end of the sequence or closely ********************
					if exon[locus][0][2]=="+":   #************* Positive strand
						#print "We start from", cdsSeq[-3:]
						if not cdsSeq[-3:]=="TGA" and not cdsSeq[-3:]=="TAA" and not cdsSeq[-3:]=="TAG":
							notes += "- Looking for stop codon downstream....\n"
							foundStopCodon = False
							for a in range(len(sequence)-int(len(cdsSeq)/3)+30):
								newStop = genomeSeq[int(exon[locus][-1][1])+a*3:int(exon[locus][-1][1])+a*3+3]
								if newStop == "TAA" or newStop=="TGA" or newStop=="TAG":
									notes += "Valid stop codon found downstream\n"
									exon[locus][0]=(int(exon[locus][0][0]),int(exon[locus][-1][1])+a*3,exon[locus][0][2])
									gene[locus] = (int(gene[locus][0]), int(exon[locus][-1][1])+a*3, gene[locus][2])
									cdsSeq = ""
									for item in exon[locus]:
										cdsSeq+=genomeSeq[int(item[0])-1:int(item[1])]
									foundStopCodon = True
									notes += "- Stop codon refined " + \
										str(a)+" codon downstream\n"
									break


					else: # *********************** Negative Strand
						if not cdsSeq[-3:]=="TGA" and not cdsSeq[-3:]=="TAA" and not cdsSeq[-3:]=="TAG":
							notes += "- Looking for stop codon downstream....\n"
							for a in range(len(sequence)-int(len(cdsSeq)/3)+30):
								#print "New Stop codons"
								newStop = self.reverseComplement(genomeSeq[int(exon[locus][0][0])-1-a*3-3:int(exon[locus][0][0])-1-a*3])
								#print a,newStop
								if newStop == "TAA" or newStop=="TGA" or newStop=="TAG":
									notes += "Valid stop codon found downstream\n"
									exon[locus][0]=(int(exon[locus][0][0])-a*3,exon[locus][0][1],exon[locus][0][2])
									gene[locus] =(int(exon[locus][0][0])-a*3, int(gene[locus][1]),gene[locus][2])
									cdsSeq = ""
									for b in range(len(exon[locus])-1,-1,-1):
										cdsSeq+=self.reverseComplement(genomeSeq[ int(exon[locus][b][0])-1: int(exon[locus][b][1]) ])
									foundStopCodon = True
									notes += "- Stop codon refined " + \
										str(a)+" codon(s) downstream\n"
									break



				if foundStartCodon == True and foundStopCodon == True and numCodonRefines <5 and abs(int(len(cdsSeq)/3) - len(sequence)) <10: # Write gff and cds file ********************
					#warnFile.write("A valid ORF for gene "+locus+" after prediction refinement\n")
					if  notes == "":
						notes = "A valid ORF has been found for gene "+locus+"!\n"
					else:
						notes += "A valid ORF has been found for gene "+locus+"!\n"

					#  ******************* Check CDS integrity
					cdsGood = True
					plausablePrediction = True
					for a in range(0,len(cdsSeq)-3,+3):
						if cdsSeq[a:a+3] == "TAA" or cdsSeq[a:a+3] == "TGA" or cdsSeq[a:a+3] == "TAG":
							# Check if the shorter sequence is compatible with one of the models
							newProtLen = float(a)/3.0
							plausablePrediction = False
							newProt = Seq(cdsSeq[:a+3]).translate()
							for protein in protSeqs:
								score = pairwise2.align.localxx(newProt, protSeqs[protein].seq,score_only=True)
								if score / float(len(protSeqs[protein].seq)) >= 0.8:
									plausablePrediction = True
									if exon[locus][0][2]=="+":  #Check it if the strand is positive
										newmRNALength = 0
										newExonSet = {}
										if not locus in newExonSet:
											newExonSet[locus] = []
										for item in exon[locus]:
											if int(item[1])-int(item[0]) + newmRNALength > a+3:
												newExonSet[locus].append((int(item[0]),int(item[0]) + a - newmRNALength -1 ,item[2]))
												exon[locus] = newExonSet[locus]
												gene[locus] = (int(gene[locus][0]),int(item[0]) + a - newmRNALength -1 ,gene[locus][2])
												break
											else:
												newmRNALength += int(item[1]) - int(item[0])
												newExonSet[locus].append((int(item[0]),int(item[1]),item[2]))
										cdsSeq = cdsSeq[:a+3]
									else:  #Check it if the strand is negative
										newmRNALength = 0
										newExonSet = {}
										if not locus in newExonSet:
											newExonSet[locus] = []
										for a in range(len(exon[locus])-1,-1,-1):
											if int(exon[locus][a][1])-int(exon[locus][a][0]) + newmRNALength > a+3:
												newExonSet[locus].append((int(exon[locus][a][1]) - a  + newmRNALength, int(exon[locus][a][1]),exon[locus][a][2]))
												#print "Previous exon locus",exon[locus]
												exon[locus] = newExonSet[locus]
												#print "after exon locus",exon[locus]
												#print gene[locus]
												#gene[locus] = (int(exon[locus][a][1]) - a  + newmRNALength, int(gene[locus][1]), gene[locus][2])
												break

											else:
												newmRNALength = int(exon[locus][a][1]) - int(exon[locus][a][0])
												newExonSet[locus].append((int(exon[locus][a][0]),int(exon[locus][a][1]),exon[locus][a][2]))
										cdsSeq = cdsSeq[:a+3]
										#print exon[locus][0][1]
										exon[locus][0] = (exon[locus][0][0], exon[locus][0][1]-6,exon[locus][0][2])
										#print exon[locus][0][1]
									break
								else:
									plausablePrediction = False

							if plausablePrediction == False:
								break


					if plausablePrediction == False:
						cdsGood = False

						notes += "- The coding sequence is interrupted by a stop codon\n"


					if not len(cdsSeq)%3 == 0:
						cdsGood = False
						notes += "- The coding sequence length is not multiple of 3\n"




					if cdsGood == True:  # CDS passed quality check
						cdsFasta.write(">"+locus+" +\n"+cdsSeq+"\n")
						proteinFasta.write(">"+locus+" +\n"+str((Seq(cdsSeq)).translate())+"\n")
					else:
						cdsFasta.write(">"+locus+" +_pseudo\n"+cdsSeq+"\n")
						proteinFasta.write(">"+locus+" +_pseudo\n"+str((Seq(cdsSeq)).translate())+"\n")
					logFile.write(notes+"\n")




				else:
				#if foundStartCodon == False or foundStopCodon == False or abs(int(len(cdsSeq)/3) - len(sequence))>10:
					notes += "The prediction was not successfull! Now attempting a refinement....\n"
					

					
					
					print("Annotation needs refinement")
					exResult.close()
					#  ***************************************************************************
					#  ************************* Annotation refinement ***************************
					#  ***************************************************************************


					os.system(self.installationDirectory+"src/conda/bin/exonerate --model protein2genome "+tempFolder+"/tempFasta.fasta "+genomeFileNameList[c]+" --showtargetgff -s 0 -n 1 --refine full --forcegtag --minintron 35 --maxintron 10000 > "+tempFolder+"/outputExonerate")

					#Check Exonerate output *****************************************************************
					#Check the proteins gave a match in the target genome
					exResult = open(tempFolder+"/outputExonerate")
					line = exResult.readline().rstrip()
					while not "Query range:" in line:
						line = exResult.readline().rstrip()
						if line is None:
							print("WARNING no protein found for ", locus)
							warnings.append(
								"Missing gene: locus "+locus+" did not provide any alignment")
					exResult.close()

					#Reconstruct Exons  **********************************************************************
					exResult = open(tempFolder+"/outputExonerate")
					while not line == "# --- START OF GFF DUMP ---":
						line = exResult.readline().rstrip()
						if line is None:
							#print "WARNING no protein found for ", locus
							logFile.write(
								"Missing gene: locus "+locus+" did not provide any alignment")
					for a in range(10):
						line = exResult.readline().rstrip()

					gene = {}
					exon = {}
					warnings = []
					#Collect the exonerate output
					while not line == "# --- END OF GFF DUMP ---":
						line = exResult.readline().rstrip()
						if line is None:
							#print "WARNING no protein found for ", locus
							logFile.write(
								"Missing gene: locus "+locus+" did not provide any alignment")
						fields = line.split("\t")
						if not line == "# --- END OF GFF DUMP ---":
							if fields[2] == "gene":
								if not locus in gene:
									gene[locus] = (
										fields[3], fields[4], fields[6])
								else:
									if (int(fields[4]) - int(fields[3])) > (int(int(gene[locus][1])) - int(int(gene[locus][0]))):
										gene[locus] = (
											fields[3], fields[4], fields[6])

							if fields[2] == "exon":
								if not locus in exon:
									exon[locus] = []
								exon[locus].append(
									(fields[3], fields[4], fields[6]))
								
					newList = sorted(exon[locus], key=lambda x: x[1])
					exon[locus] = newList

					#Reconstruct CDS  ***********************************************************
					cdsSeq = ""
					if exon[locus][0][2] == "+":  # ************* Positive strand
						for item in exon[locus]:
							cdsSeq += genomeSeq[int(item[0]) -
												1:int(item[1])]
						cdsSeq += genomeSeq[int(item[1]):int(item[1])+3]

					else:  # *********************** Negative Strand
						for a in range(len(exon[locus])-1, -1, -1):
							cdsSeq += self.reverseComplement(
								genomeSeq[int(exon[locus][a][0])-1: int(exon[locus][a][1])])
						cdsSeq += self.reverseComplement(
							genomeSeq[int(exon[locus][a][0])-4:int(exon[locus][a][0])-1])



					# Check CDS integrity *********************************************************
					foundStartCodon = True
					foundStopCodon = True
					if not cdsSeq[:3] == "ATG" or not (cdsSeq[-3:] == "TGA" or cdsSeq[-3:] == "TAA" or cdsSeq[-3:] == "TAG"):


						notes += "Either the start or the stop codon was not found. Searching nearby....\n"

						# Look for ATG at the beginning of the sequence or closely ********************
						if exon[locus][0][2] == "+":  # ************* Positive strand
							if not cdsSeq[:3] == "ATG":
								foundStartCodon = False
								print("- Looking for start codon upstream....\n")
								for a in range(len(sequence)-int(len(cdsSeq)/3)+30):
									newStart = genomeSeq[int(
										exon[locus][0][0])-1-a*3-3:int(exon[locus][0][0])-1-a*3]
									if newStart == "ATG":
										notes += "Valid start codon found upstream!\n"
										exon[locus][0] = (
											int(exon[locus][0][0])-a*3-3, exon[locus][0][1], exon[locus][0][2])
										gene[locus] = (
											int(exon[locus][0][0])-a*3-3, int(gene[locus][1]), gene[locus][2])
										cdsSeq = ""
										for item in exon[locus]:
											cdsSeq += genomeSeq[int(
												item[0])-1:int(item[1])]
										foundStartCodon = True
										notes += "- Start codon refined  " + \
											str(a)+" codons upstream\n"
										numCodonRefines = a
										break
									if newStart == "TGA" or newStart == "TAA" or newStart == "TAG":
										foundStartCodon = False
										notes += "- Found a stop codon while searching for start codon upstream!\nStart codon could not be found upstream\n"
										break
								#If the new start codon was not found in the region upstream then the downstream region is searched
								if foundStartCodon == False:
									notes += "- Looking for start codon downstream....\n"
									for a in range(len(sequence)-int(len(cdsSeq)/3)+30):
										newStart = genomeSeq[int(
											exon[locus][0][0])-1+a*3+3:int(exon[locus][0][0])-1+a*3]
										if newStart == "ATG":
											notes += "Valid start codon found downstream!\n"
											exon[locus][0] = (
												int(exon[locus][0][0])+a*3+3, exon[locus][0][1], exon[locus][0][2])
											gene[locus] = (
												int(exon[locus][0][0])+a*3+3, int(gene[locus][1]), gene[locus][2])
											cdsSeq = ""
											for item in exon[locus]:
												cdsSeq += genomeSeq[int(
													item[0])-1:int(item[1])]
											foundStartCodon = True
											notes += "- Start codon refined  " + \
												str(a) + \
												" codons upstream\n"
											numCodonRefines = a
											break
										if newStart == "TGA" or newStart == "TAA" or newStart == "TAG":
											foundStartCodon = False
											notes += "- Found a stop codon while searching for start codon downstream!\nStart codon could not be found downstream\n"
											break

						else:  # *********************** Negative Strand
							# RL6 start with alternative start codon
							if not cdsSeq[:3] == "ATG" or not (cdsSeq[:3] == "TTG" and locus == "RL6"):
								foundStartCodon = False
								notes += "- Looking for start codon upstream....\n"
								for a in range(len(sequence)-int(len(cdsSeq)/3)+30):
									#print "New start codons"
									newStart = self.reverseComplement(
										genomeSeq[int(exon[locus][-1][1])+a*3:int(exon[locus][-1][1])+a*3+3])
									#print newStart
									if newStart == "ATG":
										notes += "Valid start codon found upstrean!\n"
										exon[locus][-1] = (int(exon[locus][-1][0]), int(
											exon[locus][-1][1])+a*3+3, exon[locus][-1][2])
										gene[locus] = (int(gene[locus][0]), int(
											exon[locus][-1][1])+a*3+3, gene[locus][2])
										cdsSeq = ""
										for a1 in range(len(exon[locus])-1, -1, -1):
											cdsSeq += self.reverseComplement(
												genomeSeq[int(exon[locus][a1][0])-1: int(exon[locus][a1][1])])
										foundStartCodon = True
										notes += "- Start codon refined  " + \
											str(a)+" codons upstream\n"
										numCodonRefines = a
										break
									if newStart == "TGA" or newStart == "TAA" or newStart == "TAG":
										foundStartCodon = False
										notes += "- Found a stop codon while searching for start codon upstream!\nStart codon could not be found upstream\n"
										break
							#If the new start codon was not found in the region upstream then the downstream region is searched
								if foundStartCodon == False:
									notes += "- Looking for start codon downstream....\n"
									for a in range(len(sequence)-int(len(cdsSeq)/3)+30):
										#print "New start codons"
										newStart = self.reverseComplement(genomeSeq[int(exon[locus][-1][1])-a*3-3:int(exon[locus][-1][1])-a*3])
										#print newStart
										if newStart == "ATG":
											notes += "Valid start codon found dowstream!\n"
											exon[locus][-1] = (int(exon[locus][-1][0]), int(exon[locus][-1][1])-a*3, exon[locus][-1][2])
											gene[locus] = (int(gene[locus][0]), int(
												exon[locus][-1][1])-a*3-3, gene[locus][2])
											cdsSeq = ""
											for a1 in range(len(exon[locus])-1, -1, -1):
												cdsSeq += self.reverseComplement(
													genomeSeq[int(exon[locus][a1][0])-1: int(exon[locus][a1][1])])
											foundStartCodon = True
											notes += "- Start codon refined  " + \
												str(a) + \
												" codons dowstream\n"
											numCodonRefines = a
											break
										if newStart == "TGA" or newStart == "TAA" or newStart == "TAG":
											foundStartCodon = False
											notes += "- Found a stop codon while searching for start codon downstream! \nStart codon could not be found downstream\n"
											break
								if foundStartCodon == False:
									notes += "Start codon could not be found at this stage\n"


						# Look for Stop codon at the end of the sequence or closely ********************
						if exon[locus][0][2] == "+":  # ************* Positive strand
							if not cdsSeq[-3:] == "TGA" and not cdsSeq[-3:] == "TAA" and not cdsSeq[-3:] == "TAG":
								notes += "- Looking for stop codon downstream....\n"
								foundStopCodon = False
								for a in range(len(sequence)-int(len(cdsSeq)/3)+30):
									newStop = genomeSeq[int(
										exon[locus][-1][1])+a*3:int(exon[locus][-1][1])+a*3+3]
									if newStop == "TAA" or newStop == "TGA" or newStop == "TAG":
										notes += "Valid stop codon found downstream!\n"
										exon[locus][0] = (int(exon[locus][0][0]), int(
											exon[locus][-1][1])+a*3, exon[locus][0][2])
										gene[locus] = (int(gene[locus][0]), int(
											exon[locus][-1][1])+a*3, gene[locus][2])
										cdsSeq = ""
										for item in exon[locus]:
											cdsSeq += genomeSeq[int(
												item[0])-1:int(item[1])]
										foundStopCodon = True
										notes += "- Stop codon refined " + \
											str(a) + \
											" codon downstream\n"
										break

						else:  # *********************** Negative Strand
							if not cdsSeq[-3:] == "TGA" and not cdsSeq[-3:] == "TAA" and not cdsSeq[-3:] == "TAG":
								foundStopCodon = False
								notes += "- Looking for stop codon downstream....\n"
								for a in range(len(sequence)-int(len(cdsSeq)/3)+30):
									#print "New Stop codons"
									newStop = self.reverseComplement(
										genomeSeq[int(exon[locus][0][0])-1-a*3-3:int(exon[locus][0][0])-1-a*3])

									if newStop == "TAA" or newStop == "TGA" or newStop == "TAG":
										notes += "Valid stop codon found downstream!\n"

										exon[locus][0] = (
											int(exon[locus][0][0])-a*3, exon[locus][0][1], exon[locus][0][2])
										gene[locus] = (
											int(exon[locus][0][0])-a*3, int(gene[locus][1]), gene[locus][2])

										cdsSeq = ""
										for b in range(len(exon[locus])-1, -1, -1):
											cdsSeq += self.reverseComplement(
												genomeSeq[int(exon[locus][b][0])-1: int(exon[locus][b][1])])
										foundStopCodon = True
										notes += "- Stop codon refined " + \
											str(a) + \
											" codon downstream\n"
										break

						logFile.write(notes)#+"\n")

					# Write gff and cds file ********************
					if foundStartCodon == True and foundStopCodon == True:

						logFile.write("A valid ORF for gene "+locus+" after prediction refinement!\n")
						if notes == "":
							notes = "A valid ORF for gene "+locus+" after prediction refinement!\n"
						else:
							notes += "A valid ORF has been found for gene "+locus+" after prediction refinement!\n"

						#  ******************* Check CDS integrity
						cdsGood = True
						plausablePrediction = True

						for a in range(0, len(cdsSeq)-3, +3):
							plausablePrediction = True
							if cdsSeq[a:a+3] == "TAA" or cdsSeq[a:a+3] == "TGA" or cdsSeq[a:a+3] == "TAG":
								# Check if the shorter sequence is compatible with one of the models
								newProtLen = float(a)/3.0
								plausablePrediction = False
								newProt = Seq(cdsSeq[:a+3]).translate()
								for protein in protSeqs:
									score = pairwise2.align.localxx(newProt, protSeqs[protein].seq,score_only=True)
									if score / float(len(protSeqs[protein].seq)) >= 0.8:
										plausablePrediction = True
										# Check it if the strand is positive
										if exon[locus][0][2] == "+":
											newmRNALength = 0
											newExonSet = {}
											if not locus in newExonSet:
												newExonSet[locus] = []
											for item in exon[locus]:
												if int(item[1])-int(item[0]) + newmRNALength > a+3:
													newExonSet[locus].append(
														(int(item[0]), int(item[0]) + a - newmRNALength - 1, item[2]))
													exon[locus] = newExonSet[locus]
													gene[locus] = (int(gene[locus][0]), int(
														item[0]) + a - newmRNALength - 1, gene[locus][2])
													break
												else:
													newmRNALength += int(
														item[1]) - int(item[0])
													newExonSet[locus].append(
														(int(item[0]), int(item[1]), item[2]))
											cdsSeq = cdsSeq[:a+3]
										else:  # Check it if the strand is negative
											newmRNALength = 0
											newExonSet = {}
											if not locus in newExonSet:
												newExonSet[locus] = []
											for a in range(len(exon[locus])-1, -1, -1):
												if int(exon[locus][a][1])-int(exon[locus][a][0]) + newmRNALength > a+3:
													newExonSet[locus].append(
														(int(exon[locus][a][1]) - a + newmRNALength, int(exon[locus][a][1]), exon[locus][a][2]))
													#print "Previous exon locus",exon[locus]
													exon[locus] = newExonSet[locus]
													#print "after exon locus",exon[locus]
													#print gene[locus]
													#gene[locus] = (int(exon[locus][a][1]) - a  + newmRNALength, int(gene[locus][1]), gene[locus][2])
													break

												else:
													newmRNALength = int(
														exon[locus][a][1]) - int(exon[locus][a][0])
													newExonSet[locus].append(
														(int(exon[locus][a][0]), int(exon[locus][a][1]), exon[locus][a][2]))
											cdsSeq = cdsSeq[:a+3]
											#print exon[locus][0][1]
											exon[locus][0] = (
												exon[locus][0][0], exon[locus][0][1]-6, exon[locus][0][2])
											#print exon[locus][0][1]
										break
									else:
										plausablePrediction = False

								if plausablePrediction == False:
									break




						if plausablePrediction==False:

							cdsGood = False

							notes += "- The coding sequence is interrupted by a stop codon\n"
							#break

						if not len(cdsSeq) % 3 == 0:
							cdsGood = False
							notes += "- The coding sequence length is not multiple of 3\n"
							

						if cdsGood == True:  # CDS passed quality check
							cdsFasta.write(">"+locus+" +\n"+cdsSeq+"\n")
							proteinFasta.write(">"+locus+" +\n"+str((Seq(cdsSeq)).translate())+"\n")

						else: # CDS DID NOT passed quality check
							cdsFasta.write(">"+locus+" +_pseudo\n"+cdsSeq+"\n")
							proteinFasta.write(">"+locus+" +_pseudo\n"+str((Seq(cdsSeq)).translate())+"\n")
							logFile.write(notes+"\n")
							


					else:
						logFile.write("An incomplete ORF for gene "+locus+" after prediction refinement!\n")
						if notes == "":
							notes = "An incomplete ORF for gene "+locus+" after prediction refinement!\n"
						else:
							notes += "An imcoplete ORF has been found for gene "+locus+" after prediction refinement!\n"

						#  ******************* Check CDS integrity
						cdsGood = True
						plausablePrediction = True
						for a in range(0, len(cdsSeq)-3, +3):
							if cdsSeq[a:a+3] == "TAA" or cdsSeq[a:a+3] == "TGA" or cdsSeq[a:a+3] == "TAG":
								# Check if the shorter sequence is compatible with one of the models
								newProtLen = float(a)/3.0
								plausablePrediction = False
								newProt = Seq(cdsSeq[:a+3]).translate()
								for protein in protSeqs:
									score = pairwise2.align.localxx(newProt, protSeqs[protein].seq,score_only=True)
									if score / float(len(protSeqs[protein].seq)) >= 0.8:
										plausablePrediction = True
										# Check it if the strand is positive
										if exon[locus][0][2] == "+":
											newmRNALength = 0
											newExonSet = {}
											if not locus in newExonSet:
												newExonSet[locus] = []
											for item in exon[locus]:
												if int(item[1])-int(item[0]) + newmRNALength > a+3:
													newExonSet[locus].append(
														(int(item[0]), int(item[0]) + a - newmRNALength - 1, item[2]))
													exon[locus] = newExonSet[locus]
													gene[locus] = (int(gene[locus][0]), int(
														item[0]) + a - newmRNALength - 1, gene[locus][2])
													break
												else:
													newmRNALength += int(
														item[1]) - int(item[0])
													newExonSet[locus].append(
														(int(item[0]), int(item[1]), item[2]))
											cdsSeq = cdsSeq[:a+3]
										else:  # Check it if the strand is negative
											newmRNALength = 0
											newExonSet = {}
											if not locus in newExonSet:
												newExonSet[locus] = []
											for a in range(len(exon[locus])-1, -1, -1):
												if int(exon[locus][a][1])-int(exon[locus][a][0]) + newmRNALength > a+3:
													newExonSet[locus].append(
														(int(exon[locus][a][1]) - a + newmRNALength, int(exon[locus][a][1]), exon[locus][a][2]))
													#print "Previous exon locus",exon[locus]
													exon[locus] = newExonSet[locus]
													#print "after exon locus",exon[locus]
													#print gene[locus]
													#gene[locus] = (int(exon[locus][a][1]) - a  + newmRNALength, int(gene[locus][1]), gene[locus][2])
													break

												else:
													newmRNALength = int(
														exon[locus][a][1]) - int(exon[locus][a][0])
													newExonSet[locus].append(
														(int(exon[locus][a][0]), int(exon[locus][a][1]), exon[locus][a][2]))
											cdsSeq = cdsSeq[:a+3]
											#print exon[locus][0][1]
											exon[locus][0] = (
												exon[locus][0][0], exon[locus][0][1]-6, exon[locus][0][2])
											#print exon[locus][0][1]
										break
						if plausablePrediction==False:
							cdsGood = False
							
							notes += "- The coding sequence is interrupted by a stop codon\n"


						if not len(cdsSeq) % 3 == 0:
							cdsGood = False
							notes += "- The coding sequence length is not multiple of 3\n"
							

						if cdsGood == True:  # CDS passed quality check
							

							cdsFasta.write(">"+locus+" +\n"+cdsSeq+"\n")
							proteinFasta.write(">"+locus+" +\n"+str((Seq(cdsSeq)).translate())+"\n")

						else: # CDS DID NOT passed quality check
							logFile.write(notes+"\n")
							cdsFasta.write(">"+locus+" +_pseudo\n"+cdsSeq+"\n")
							proteinFasta.write(">"+locus+" +_pseudo\n"+str((Seq(cdsSeq)).translate())+"\n")

							













				os.system("rm -f "+tempFolder+"/tempFasta.fasta "+tempFolder+"/output*")
			os.system("mv "+tempFolder+" "+self.installationDirectory+"/DB/")







class Ui_Form(object):
	def setupUi(self, Form, installationDirectory):
		Form.setObjectName("Form")
		Form.resize(837, 312)
		self.label = QtWidgets.QLabel(Form)
		self.label.setGeometry(QtCore.QRect(10, 70, 141, 16))
		self.label.setObjectName("label")
		self.label_2 = QtWidgets.QLabel(Form)
		self.label_2.setGeometry(QtCore.QRect(10, 10, 131, 16))
		self.label_2.setObjectName("label_2")
		self.label_3 = QtWidgets.QLabel(Form)
		self.label_3.setGeometry(QtCore.QRect(10, 130, 151, 16))
		self.label_3.setObjectName("label_3")
		self.sampleID_lineEdit = QtWidgets.QLineEdit(Form)
		self.sampleID_lineEdit.setGeometry(QtCore.QRect(10, 30, 151, 21))
		self.sampleID_lineEdit.setObjectName("sampleID_lineEdit")
		self.patientCode_lineEdit = QtWidgets.QLineEdit(Form)
		self.patientCode_lineEdit.setGeometry(QtCore.QRect(10, 90, 151, 21))
		self.patientCode_lineEdit.setObjectName("patientCode_lineEdit")
		self.sampleDate_lineEdit = QtWidgets.QLineEdit(Form)
		self.sampleDate_lineEdit.setGeometry(QtCore.QRect(10, 150, 151, 21))
		self.sampleDate_lineEdit.setObjectName("sampleDate_lineEdit")
		self.label_4 = QtWidgets.QLabel(Form)
		self.label_4.setGeometry(QtCore.QRect(10, 190, 151, 16))
		self.label_4.setObjectName("label_4")
		self.compartment_lineEdit = QtWidgets.QLineEdit(Form)
		self.compartment_lineEdit.setGeometry(QtCore.QRect(10, 210, 151, 21))
		self.compartment_lineEdit.setObjectName("compartment_lineEdit")
		self.label_5 = QtWidgets.QLabel(Form)
		self.label_5.setGeometry(QtCore.QRect(10, 250, 151, 16))
		self.label_5.setObjectName("label_5")
		self.fileName_lineEdit = QtWidgets.QLineEdit(Form)
		self.fileName_lineEdit.setGeometry(QtCore.QRect(10, 270, 501, 21))
		self.fileName_lineEdit.setObjectName("fileName_lineEdit")
		self.openFile_button = QtWidgets.QPushButton(Form)
		self.openFile_button.setGeometry(QtCore.QRect(520, 265, 131, 32))
		self.openFile_button.setObjectName("openFile_button")
		self.textEdit = QtWidgets.QTextEdit(Form)
		self.textEdit.setGeometry(QtCore.QRect(340, 30, 471, 221))
		self.textEdit.setObjectName("textEdit")
		self.label_6 = QtWidgets.QLabel(Form)
		self.label_6.setGeometry(QtCore.QRect(340, 10, 131, 16))
		self.label_6.setObjectName("label_6")
		self.pushButton = QtWidgets.QPushButton(Form)
		self.pushButton.setGeometry(QtCore.QRect(650, 265, 171, 32))
		self.pushButton.setObjectName("pushButton")

		self.retranslateUi(Form)
		QtCore.QMetaObject.connectSlotsByName(Form)

		self.openFile_button.clicked.connect(self.openGenomeFile)
		self.pushButton.clicked.connect(self.runAnnoation)


	def openGenomeFile(self):
		filename = QtWidgets.QFileDialog.getOpenFileName(None,"Select genome fasta file","./","*.fasta *.fa *.bulk")
		self.fileName_lineEdit.setText(str(filename[0]))





	def runAnnoation(self):

			if len(self.sampleID_lineEdit.text()) == 0 or len(self.patientCode_lineEdit.text()) == 0 or len(self.sampleDate_lineEdit.text()) == 0 or len(self.compartment_lineEdit.text()) == 0:
				msg = QtWidgets.QMessageBox()
				msg.setText("Sample type/number/date/compartment are mandatory fields")
				msg.exec() 
				return

			if len(self.fileName_lineEdit.text()) == 0:
				msg = QtWidgets.QMessageBox()
				msg.setText("A genomic sequence in fasta format must be provided")
				msg.exec() 
				return


		 	# Step 2: Create a QThread object
			self.thread = QThread()
			# Step 3: Create a worker object
			self.worker = Worker(installationDirectory,str(self.fileName_lineEdit.text()),str(self.sampleID_lineEdit.text()),self.patientCode_lineEdit.text(),self.compartment_lineEdit.text(),self.sampleDate_lineEdit.text())
			# Step 4: Move worker to the thread
			self.worker.moveToThread(self.thread)
			# Step 5: Connect signals and slots
			self.thread.started.connect(self.worker.run)
			self.worker.finished.connect(self.thread.quit)
			#self.worker.finished.connect(self.worker.deleteLater)
			#self.thread.finished.connect(self.thread.deleteLater)
			self.worker.progress.connect(self.reportProgress)
			# Step 6: Start the thread
			self.thread.start()



		
				
	def reportProgress(self,outputLog):
		self.textEdit.append(str(outputLog))
		self.textEdit.repaint()


	

	def retranslateUi(self, Form):
		_translate = QtCore.QCoreApplication.translate
		Form.setWindowTitle(_translate("Form", "Entry from genome sequence"))
		self.label.setText(_translate("Form", "Sample name (required)"))
		self.label_2.setText(_translate("Form", "Sample type (required)"))
		self.label_3.setText(_translate("Form", "Sample date (MMDDYY)"))
		self.label_4.setText(_translate("Form", "Compartment"))
		self.label_5.setText(_translate("Form", "File name"))
		self.openFile_button.setText(_translate("Form", "Open file"))
		self.label_6.setText(_translate("Form", "Log"))
		self.pushButton.setText(_translate("Form", "Retrieve sequences"))


import sys
if __name__ == "__main__":
	app = QtWidgets.QApplication(sys.argv)
	Form = QtWidgets.QWidget()
	installationDirectory = sys.argv[1]
	ui = Ui_Form()
	ui.setupUi(Form,installationDirectory)
	Form.show()
	sys.exit(app.exec_())