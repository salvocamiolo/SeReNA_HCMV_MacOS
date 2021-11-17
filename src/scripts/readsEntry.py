# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'readsEntry.ui'
#
# Created by: PyQt5 UI code generator 5.9.2
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5.QtCore import QObject, QThread, pyqtSignal
import os
from Bio import SeqIO
from Bio.Seq import Seq
import numpy as np
import time

class Worker(QObject):
	finished = pyqtSignal()
	progress = pyqtSignal(str)

	def __init__(self,ii, read1,read2,numThreads,sampleID, sampleType, sampleDate, sampleCompartment):
		super(QObject, self).__init__()
		self.installationDirectory = ii
		self.read1 = read1
		self.read2 = read2
		self.numThreads = numThreads
		self.sampleID = sampleID
		self.sampleType = sampleType
		self.sampleDate = sampleDate
		self.sampleCompartment = sampleCompartment

	def run(self):
		installationDirectory = self.installationDirectory
		read1 = self.read1
		read2 = self.read2
		numThreads = self.numThreads
		pseudoGenes = []

		prefix = self.sampleType+"-"+self.sampleID+"-"+self.sampleDate+"-"+self.sampleCompartment

		tempFolder = prefix+"_temp"
		os.system("mkdir "+tempFolder)
		os.system(installationDirectory+"/DB/"+self.sampleID)
		logFile = open(tempFolder+"/"+prefix+"_annotationLog.txt","w")
		cdsFile = open(tempFolder+"/"+prefix+"_cds.fasta","w")
		proteinsFile = open(tempFolder+"/"+prefix+"_proteins.fasta","w")

		for geneFile in os.listdir(installationDirectory+"/data/elongedCDS"):
			if geneFile.endswith("_elongedCDS.fasta") == True:
				gene2Analyse = geneFile.replace("_elongedCDS.fasta","")
				self.progress.emit("Log_Analyzing gene "+gene2Analyse)
				self.progress.emit("Log_Finding gene specific reads.... ")
				os.system(installationDirectory+"/src/conda/bowtie2 -1 "+read1+" -2 "+ read2 + " -p "+numThreads+" -x "+installationDirectory+"/data/elongedCDS/"+geneFile+" -S "+tempFolder+"/alignment.sam ")
				self.progress.emit("Log_Parsing.... ")
				os.system(installationDirectory+"src/conda/bin/samtools view -bS -F4 -h "+tempFolder+"/alignment.sam > "+tempFolder+"/alignment.bam")
				os.system(installationDirectory+"src/conda/bin/samtools sort -o "+tempFolder+"/alignment_sorted.bam "+tempFolder+"/alignment.bam")
				os.system(installationDirectory+"src/conda/bin/samtools bam2fq "+tempFolder+"/alignment_sorted.bam > "+tempFolder+"/reads.fastq")
				self.progress.emit("Log_Assembling.... ")
				os.system("rm -rf "+tempFolder+"/outputSpades")
				os.system(installationDirectory+"src/conda/bin/spades.py -s  "+tempFolder+"/reads.fastq --careful --cov-cutoff auto -o  "+tempFolder+"/outputSpades --trusted-contigs "+installationDirectory+"/data/elongedCDS/"+geneFile+" > "+tempFolder+"/null")

				if os.path.isfile(tempFolder+"/outputSpades/scaffolds.fasta") == True:
					avGeneLen = []
					for seq_record in SeqIO.parse(installationDirectory+"/data/elongedCDS/"+geneFile,"fasta"):
						avGeneLen.append(len(str(seq_record.seq))-600)
					avGeneLen = np.mean(avGeneLen)

					numSeq = 0
					for seq_record in SeqIO.parse(tempFolder+"/outputSpades/scaffolds.fasta","fasta"):
						assembledSequence = str(seq_record.seq)
						numSeq+=1
						if numSeq>0:
							break
					
					outfile = open(tempFolder+"/assembledSequence.fasta","w")
					outfile.write(">assembled\n"+assembledSequence+"\n")
					outfile.close()
					print("A sequence was assembled covering the",float(len(assembledSequence))/avGeneLen,"of the expected")



					#Annotate CDS
					os.system(installationDirectory+"/src/conda/bin/makeblastdb -dbtype nucl -in "+tempFolder+"/assembledSequence.fasta")
					protName = (geneFile.replace("_elongedCDS.fasta",""))+"_models.fasta"

					protSeqs = {}
					for seq_record in SeqIO.parse(installationDirectory+"/data/proteinDB/"+protName,"fasta"):
						if not str(seq_record.id) in protSeqs:
							protSeqs[str(seq_record.id)] = str(seq_record.seq)
					
					strand = "+"
					os.system(installationDirectory+"/src/conda/bin/tblastn -query "+installationDirectory+"/data/proteinDB/"+protName+" -db "+tempFolder+"/assembledSequence.fasta -outfmt 6 -max_intron_length 350000 | sort -k 12rn,12rn  > "+tempFolder+"/preOutputBlast.txt")
					time.sleep(1)
					os.system("head -20 "+tempFolder+"/preOutputBlast.txt > "+tempFolder+"/preOutputBlast2.txt")
					time.sleep(1)
					os.system("sort -k 3rn,3rn -k12rn,12rn "+tempFolder+"/preOutputBlast2.txt > "+tempFolder+"/outputBlast.txt")
					os.system("rm "+tempFolder+"/preOutputBlast*")
					blastFile = open(tempFolder+"/outputBlast.txt")
					bestCoverage = 0
					blastLine = blastFile.readline().rstrip()
					blastField = blastLine.split("\t")
					if len(blastField)>1:
						alignmentLength = int(blastField[7]) - int(blastField[6])
						bestCoverage = alignmentLength
						bestProt = blastField[0]
						if int(blastField[8]) > int(blastField[9]):
							strand = "-"


						
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
								print(int(blastField[8]),int(blastField[9]))
								if int(blastField[8]) > int(blastField[9]):
									strand = "-"


						print("Best prot",bestProt)

						outfile = open(tempFolder+"/bestProt.fasta","w")
						outfile.write(">bestProt\n"+protSeqs[bestProt]+"\n")
						outfile.close()

						if strand=="-":
							os.system("rm -rf "+tempFolder+"/assembledSequence.fasta*")
							outfile = open(tempFolder+"/assembledSequence.fasta","w")
							assembledSequence = str( (Seq(assembledSequence)).reverse_complement() )
							outfile.write(">assembled\n"+assembledSequence+"\n")
							outfile.close()


						os.system(installationDirectory+"src/conda/bin/exonerate --model protein2genome "+tempFolder+"/bestProt.fasta "+tempFolder+"/assembledSequence.fasta --showtargetgff -s 0 -n 1 --forcegtag --minintron 35 --maxintron 10000  > "+tempFolder+"/outputExonerate")
						os.system(" awk '$3==\"exon\"' "+tempFolder+"/outputExonerate > "+tempFolder+"/exons")
						infile = open(tempFolder+"/exons")
						cdsSeq = ""
						numExon = 0
						while True:
							line = infile.readline().rstrip()
							if not line:
								break
							fields = line.split("\t")
							print(fields[3],fields[4],fields[5],fields[6])
							print(assembledSequence[int(fields[3])-1:int(fields[4])])
							print("\n")
							cdsSeq+= assembledSequence[int(fields[3])-1:int(fields[4])]
						cdsSeq+=assembledSequence[int(fields[4]):int(fields[4])+3]
						infile.close()
						print(cdsSeq)
						proteinSequence = str((Seq(cdsSeq)).translate())
						if float(len(cdsSeq))/float(avGeneLen) >= 0.8:
							if not "*" in proteinSequence[:-1] and ( (cdsSeq[:3] == "ATG" and not gene2Analyse == "UL30A") or gene2Analyse=="UL30A") and (cdsSeq[-3:]=="TGA" or cdsSeq[-3:]=="TAG" or cdsSeq[-3:]=="TAA"): #there is not stop codon within the CDS
								cdsFile.write(">"+gene2Analyse+"\n"+cdsSeq+"\n")
								proteinsFile.write(">"+gene2Analyse+"\n"+proteinSequence+"\n")
								self.progress.emit("Fou_"+gene2Analyse)
							else:
								cdsFile.write(">"+gene2Analyse+" Pseudo\n"+cdsSeq+"\n")
								proteinsFile.write(">"+gene2Analyse+" Pseudo\n"+proteinSequence+"\n")
								pseudoGenes.append(gene2Analyse)
								print("The coding sequence is interrupted by a stop codon!")
								self.progress.emit("Pse_"+gene2Analyse)
						else:
							print(gene2Analyse,"has not been found!")
							self.progress.emit("Not_"+gene2Analyse)
					else:
						print(gene2Analyse,"has not been found!")
						self.progress.emit("Not_"+gene2Analyse)

				else:
					print(gene2Analyse,"has not been found!")
					self.progress.emit("Not_"+gene2Analyse)

		
		os.system("mkdir "+installationDirectory+"/DB/"+self.sampleID)
		os.system("mv "+tempFolder+"/*_annotationLog.txt "+installationDirectory+"/DB/"+self.sampleID+"/")
		os.system("mv "+tempFolder+"/*_cds.fasta "+installationDirectory+"/DB/"+self.sampleID+"/")
		os.system("mv "+tempFolder+"/*_proteins.fasta "+installationDirectory+"/DB/"+self.sampleID+"/")
		os.system("rm -rf "+tempFolder)

				




class Ui_Form(object):
	def setupUi(self, Form, installationDirectory):
		Form.setObjectName("Form")
		Form.resize(775, 484)
		self.read1_lineEdit = QtWidgets.QLineEdit(Form)
		self.read1_lineEdit.setGeometry(QtCore.QRect(10, 106, 231, 21))
		self.read1_lineEdit.setObjectName("read1_lineEdit")
		self.read2_lineEdit = QtWidgets.QLineEdit(Form)
		self.read2_lineEdit.setGeometry(QtCore.QRect(10, 186, 231, 21))
		self.read2_lineEdit.setObjectName("read2_lineEdit")
		self.label = QtWidgets.QLabel(Form)
		self.label.setGeometry(QtCore.QRect(10, 86, 60, 16))
		self.label.setObjectName("label")
		self.label_2 = QtWidgets.QLabel(Form)
		self.label_2.setGeometry(QtCore.QRect(10, 166, 60, 16))
		self.label_2.setObjectName("label_2")
		self.read1_pushButton = QtWidgets.QPushButton(Form)
		self.read1_pushButton.setGeometry(QtCore.QRect(247, 101, 113, 32))
		self.read1_pushButton.setObjectName("read1_pushButton")
		self.read2_pushButton = QtWidgets.QPushButton(Form)
		self.read2_pushButton.setGeometry(QtCore.QRect(250, 182, 113, 32))
		self.read2_pushButton.setObjectName("read2_pushButton")
		self.label_3 = QtWidgets.QLabel(Form)
		self.label_3.setGeometry(QtCore.QRect(10, 246, 121, 16))
		self.label_3.setObjectName("label_3")
		self.numThreads_lineEdit = QtWidgets.QLineEdit(Form)
		self.numThreads_lineEdit.setGeometry(QtCore.QRect(10, 266, 113, 21))
		self.numThreads_lineEdit.setObjectName("numThreads_lineEdit")
		self.foundGenes_textEdit = QtWidgets.QTextEdit(Form)
		self.foundGenes_textEdit.setGeometry(QtCore.QRect(400, 106, 104, 181))
		self.foundGenes_textEdit.setObjectName("foundGenes_textEdit")
		self.label_4 = QtWidgets.QLabel(Form)
		self.label_4.setGeometry(QtCore.QRect(400, 86, 60, 16))
		self.label_4.setObjectName("label_4")
		self.notFoundGenes_textEdit = QtWidgets.QTextEdit(Form)
		self.notFoundGenes_textEdit.setGeometry(QtCore.QRect(520, 106, 104, 181))
		self.notFoundGenes_textEdit.setObjectName("notFoundGenes_textEdit")
		self.label_5 = QtWidgets.QLabel(Form)
		self.label_5.setGeometry(QtCore.QRect(520, 86, 71, 16))
		self.label_5.setObjectName("label_5")
		self.pseudoGenes_textEdit = QtWidgets.QTextEdit(Form)
		self.pseudoGenes_textEdit.setGeometry(QtCore.QRect(640, 106, 104, 181))
		self.pseudoGenes_textEdit.setObjectName("pseudoGenes_textEdit")
		self.label_6 = QtWidgets.QLabel(Form)
		self.label_6.setGeometry(QtCore.QRect(640, 86, 60, 16))
		self.label_6.setObjectName("label_6")
		self.label_7 = QtWidgets.QLabel(Form)
		self.label_7.setGeometry(QtCore.QRect(10, 326, 60, 16))
		self.label_7.setObjectName("label_7")
		self.logArea_textEdit = QtWidgets.QTextEdit(Form)
		self.logArea_textEdit.setGeometry(QtCore.QRect(10, 346, 621, 121))
		self.logArea_textEdit.setObjectName("logArea_textEdit")
		self.run_pushButton = QtWidgets.QPushButton(Form)
		self.run_pushButton.setGeometry(QtCore.QRect(640, 440, 113, 32))
		self.run_pushButton.setObjectName("run_pushButton")
		self.label_8 = QtWidgets.QLabel(Form)
		self.label_8.setGeometry(QtCore.QRect(10, 10, 91, 16))
		self.label_8.setObjectName("label_8")
		self.label_9 = QtWidgets.QLabel(Form)
		self.label_9.setGeometry(QtCore.QRect(170, 10, 91, 16))
		self.label_9.setObjectName("label_9")
		self.label_10 = QtWidgets.QLabel(Form)
		self.label_10.setGeometry(QtCore.QRect(320, 10, 91, 16))
		self.label_10.setObjectName("label_10")
		self.label_11 = QtWidgets.QLabel(Form)
		self.label_11.setGeometry(QtCore.QRect(470, 10, 91, 16))
		self.label_11.setObjectName("label_11")
		self.sampleName_lineEdit = QtWidgets.QLineEdit(Form)
		self.sampleName_lineEdit.setGeometry(QtCore.QRect(10, 30, 113, 21))
		self.sampleName_lineEdit.setObjectName("sampleName_lineEdit")
		self.sampleType_lineEdit = QtWidgets.QLineEdit(Form)
		self.sampleType_lineEdit.setGeometry(QtCore.QRect(170, 30, 113, 21))
		self.sampleType_lineEdit.setObjectName("sampleType_lineEdit")
		self.sampleDate_lineEdit = QtWidgets.QLineEdit(Form)
		self.sampleDate_lineEdit.setGeometry(QtCore.QRect(320, 30, 113, 21))
		self.sampleDate_lineEdit.setObjectName("sampleDate_lineEdit")
		self.compartment_lineEdit = QtWidgets.QLineEdit(Form)
		self.compartment_lineEdit.setGeometry(QtCore.QRect(470, 30, 113, 21))
		self.compartment_lineEdit.setObjectName("compartment_lineEdit")

		self.retranslateUi(Form)
		QtCore.QMetaObject.connectSlotsByName(Form)

		self.read1_pushButton.clicked.connect(self.openRead1)
		self.read2_pushButton.clicked.connect(self.openRead2)
		self.run_pushButton.clicked.connect(self.runAssembler)
		

	def openRead1(self):
		filename = QtWidgets.QFileDialog.getOpenFileName(None,"Select file for read 1","./","*.fastq *.f1")
		self.read1_lineEdit.setText(str(filename[0]))

	def openRead2(self):
		filename = QtWidgets.QFileDialog.getOpenFileName(None,"Select file for read 2","./","*.fastq *.f1")
		self.read2_lineEdit.setText(str(filename[0]))

	def runAssembler(self):
		
		if len(self.sampleName_lineEdit.text()) == 0 or len(self.sampleType_lineEdit.text()) == 0 or len(self.sampleDate_lineEdit.text()) == 0 or len(self.compartment_lineEdit.text()) == 0:
			msg = QtWidgets.QMessageBox()
			msg.setText("Sample type/number/date/compartment are mandatory fields")
			msg.exec() 
			return
		
		if len(self.read1_lineEdit.text()) == 0:
			msg = QtWidgets.QMessageBox()
			msg.setText("No file has been selected for read 1")
			msg.exec() 
			return

		if len(self.read2_lineEdit.text()) == 0:
			msg = QtWidgets.QMessageBox()
			msg.setText("No file has been selected for read 2")
			msg.exec() 
			return

		if len(self.numThreads_lineEdit.text()) == 0 or (self.numThreads_lineEdit.text()).isnumeric() == False:
			msg = QtWidgets.QMessageBox()
			msg.setText("Please select a numeric value for the number of threads field")
			msg.exec() 
			return



		# Step 2: Create a QThread object
		self.thread = QThread()
		# Step 3: Create a worker object
		self.worker = Worker(installationDirectory,str(self.read1_lineEdit.text()),str(self.read2_lineEdit.text()),str(self.numThreads_lineEdit.text()), str(self.sampleName_lineEdit.text()),str(self.sampleType_lineEdit.text()),self.sampleDate_lineEdit.text(),self.compartment_lineEdit.text())
		# Step 4: Move worker to the thread
		self.worker.moveToThread(self.thread)
		# Step 5: Connect signals and slots
		self.thread.started.connect(self.worker.run)
		self.worker.finished.connect(self.thread.quit)
		self.worker.finished.connect(self.worker.deleteLater)
		self.thread.finished.connect(self.thread.deleteLater)
		self.worker.progress.connect(self.reportProgress)
		# Step 6: Start the thread
		self.thread.start()

	def reportProgress(self,outputLog):
		if outputLog[:4]=="Log_":
			self.logArea_textEdit.append(outputLog[4:])
			self.logArea_textEdit.repaint()
		if outputLog[:4] == "Fou_":
			self.foundGenes_textEdit.append(outputLog[4:])
			self.foundGenes_textEdit.repaint()
		if outputLog[:4] == "Not_":
			self.notFoundGenes_textEdit.append(outputLog[4:])
			self.notFoundGenes_textEdit.repaint()
		if outputLog[:4] == "Pse_":
			self.pseudoGenes_textEdit.append(outputLog[4:])
			self.pseudoGenes_textEdit.repaint()


	def retranslateUi(self, Form):
		_translate = QtCore.QCoreApplication.translate
		self.numThreads_lineEdit.setText("1")
		Form.setWindowTitle(_translate("Form", "Entry from reads"))
		self.label.setText(_translate("Form", "Read 1"))
		self.label_2.setText(_translate("Form", "Read 2"))
		self.read1_pushButton.setText(_translate("Form", "Open"))
		self.read2_pushButton.setText(_translate("Form", "Open"))
		self.label_3.setText(_translate("Form", "Number of threads"))
		self.label_4.setText(_translate("Form", "Found"))
		self.label_5.setText(_translate("Form", "Not found"))
		self.label_6.setText(_translate("Form", "Pseudo"))
		self.label_7.setText(_translate("Form", "Log"))
		self.run_pushButton.setText(_translate("Form", "Run"))
		self.label_8.setText(_translate("Form", "Sample name"))
		self.label_9.setText(_translate("Form", "Sample type"))
		self.label_10.setText(_translate("Form", "Sample date"))
		self.label_11.setText(_translate("Form", "Compartment"))

import sys
if __name__ =="__main__":
	App = QtWidgets.QApplication(sys.argv)
	Form = QtWidgets.QWidget()
	ui = Ui_Form()
	installationDirectory = sys.argv[1]
	ui.setupUi(Form, installationDirectory)
	Form.show()
	sys.exit(App.exec_())


