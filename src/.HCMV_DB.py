# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'mainWindow.ui'
#
# Created by: PyQt5 UI code generator 5.9.2
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5.QtWidgets import QFileDialog
import os
from Bio import SeqIO


class Ui_MainWindow(object):
	def setupUi(self, MainWindow):
		MainWindow.setObjectName("MainWindow")
		MainWindow.resize(760, 600)
		self.centralwidget = QtWidgets.QWidget(MainWindow)
		self.centralwidget.setObjectName("centralwidget")
		self.patientTypeList = QtWidgets.QListWidget(self.centralwidget)
		self.patientTypeList.setGeometry(QtCore.QRect(20, 30, 121, 111))
		self.patientTypeList.setObjectName("patientTypeList")
		self.label = QtWidgets.QLabel(self.centralwidget)
		self.label.setGeometry(QtCore.QRect(20, 10, 101, 16))
		self.label.setObjectName("label")
		self.patientNumberList = QtWidgets.QListWidget(self.centralwidget)
		self.patientNumberList.setGeometry(QtCore.QRect(150, 30, 121, 111))
		self.patientNumberList.setObjectName("patientNumberList")
		self.label_2 = QtWidgets.QLabel(self.centralwidget)
		self.label_2.setGeometry(QtCore.QRect(150, 10, 101, 16))
		self.label_2.setObjectName("label_2")
		self.label_3 = QtWidgets.QLabel(self.centralwidget)
		self.label_3.setGeometry(QtCore.QRect(20, 150, 101, 16))
		self.label_3.setObjectName("label_3")
		self.sampleDateList = QtWidgets.QListWidget(self.centralwidget)
		self.sampleDateList.setGeometry(QtCore.QRect(20, 170, 121, 121))
		self.sampleDateList.setObjectName("sampleDateList")
		self.compartmentList = QtWidgets.QListWidget(self.centralwidget)
		self.compartmentList.setGeometry(QtCore.QRect(150, 170, 121, 121))
		self.compartmentList.setObjectName("compartmentList")
		self.label_4 = QtWidgets.QLabel(self.centralwidget)
		self.label_4.setGeometry(QtCore.QRect(150, 150, 101, 16))
		self.label_4.setObjectName("label_4")
		self.label_5 = QtWidgets.QLabel(self.centralwidget)
		self.label_5.setGeometry(QtCore.QRect(360, 10, 101, 16))
		self.label_5.setObjectName("label_5")
		self.geneList = QtWidgets.QListWidget(self.centralwidget)
		self.geneList.setGeometry(QtCore.QRect(360, 30, 101, 261))
		self.geneList.setObjectName("geneList")
		self.label_7 = QtWidgets.QLabel(self.centralwidget)
		self.label_7.setGeometry(QtCore.QRect(570, 10, 141, 16))
		self.label_7.setObjectName("label_7")
		self.numSeqText = QtWidgets.QTextEdit(self.centralwidget)
		self.numSeqText.setGeometry(QtCore.QRect(570, 30, 161, 21))
		self.numSeqText.setObjectName("numSeqText")
		self.addCDSButton = QtWidgets.QPushButton(self.centralwidget)
		self.addCDSButton.setGeometry(QtCore.QRect(570, 100, 171, 32))
		self.addCDSButton.setObjectName("addCDSButton")
		self.addProteinsButton = QtWidgets.QPushButton(self.centralwidget)
		self.addProteinsButton.setGeometry(QtCore.QRect(570, 140, 171, 32))
		self.addProteinsButton.setObjectName("addProteinsButton")
		self.treeButton = QtWidgets.QPushButton(self.centralwidget)
		self.treeButton.setGeometry(QtCore.QRect(570, 220, 171, 32))
		self.treeButton.setObjectName("treeButton")
		self.alignButton = QtWidgets.QPushButton(self.centralwidget)
		self.alignButton.setGeometry(QtCore.QRect(570, 180, 171, 32))
		self.alignButton.setObjectName("alignButton")
		self.clearButton = QtWidgets.QPushButton(self.centralwidget)
		self.clearButton.setGeometry(QtCore.QRect(570, 260, 171, 32))
		self.clearButton.setObjectName("clearButton")
		self.label_6 = QtWidgets.QLabel(self.centralwidget)
		self.label_6.setGeometry(QtCore.QRect(20, 310, 141, 16))
		self.label_6.setObjectName("label_6")
		self.selectedSequencesArea = QtWidgets.QPlainTextEdit(self.centralwidget)
		self.selectedSequencesArea.setGeometry(QtCore.QRect(20, 330, 721, 211))
		self.selectedSequencesArea.setObjectName("selectedSequencesArea")
		self.checkBox = QtWidgets.QCheckBox(self.centralwidget)
		self.checkBox.setGeometry(QtCore.QRect(570, 60, 121, 20))
		self.checkBox.setObjectName("checkBox")
		MainWindow.setCentralWidget(self.centralwidget)
		self.menubar = QtWidgets.QMenuBar(MainWindow)
		self.menubar.setGeometry(QtCore.QRect(0, 0, 760, 24))
		self.menubar.setObjectName("menubar")
		self.menuNew_entry = QtWidgets.QMenu(self.menubar)
		self.menuNew_entry.setObjectName("menuNew_entry")
		self.menuNew_entry_2 = QtWidgets.QMenu(self.menuNew_entry)
		self.menuNew_entry_2.setObjectName("menuNew_entry_2")
		self.menuNew_sequence = QtWidgets.QMenu(self.menuNew_entry)
		self.menuNew_sequence.setObjectName("menuNew_sequence")
		#self.menuSync = QtWidgets.QMenu(self.menubar)
		#self.menuSync.setObjectName("menuSync")
		self.menuFile = QtWidgets.QMenu(self.menubar)
		self.menuFile.setObjectName("menuFile")
		self.menuEdit = QtWidgets.QMenu(self.menubar)
		self.menuEdit.setObjectName("menuEdit")
		MainWindow.setMenuBar(self.menubar)
		self.statusbar = QtWidgets.QStatusBar(MainWindow)
		self.statusbar.setObjectName("statusbar")
		MainWindow.setStatusBar(self.statusbar)
		self.actionFrom_genome = QtWidgets.QAction(MainWindow)
		self.actionFrom_genome.setObjectName("actionFrom_genome")
		self.actionFrom_reads = QtWidgets.QAction(MainWindow)
		self.actionFrom_reads.setObjectName("actionFrom_reads")
		self.actionFrom_NCBI = QtWidgets.QAction(MainWindow)
		self.actionFrom_NCBI.setObjectName("actionFrom_NCBI")
		self.actionSave_alignment = QtWidgets.QAction(MainWindow)
		self.actionSave_alignment.setEnabled(False)
		self.actionSave_alignment.setObjectName("actionSave_alignment")

		self.actionUpdateDB = QtWidgets.QAction(MainWindow)
		self.actionUpdateDB.setEnabled(False)
		self.actionUpdateDB.setObjectName("actionUpdateDB")


		self.actionSave_tree = QtWidgets.QAction(MainWindow)
		self.actionSave_tree.setEnabled(False)
		self.actionSave_tree.setObjectName("actionSave_tree")
		self.actionSave_sequences = QtWidgets.QAction(MainWindow)
		self.actionSave_sequences.setEnabled(False)
		self.actionSave_sequences.setObjectName("actionSave_sequences")
		self.actionQuit = QtWidgets.QAction(MainWindow)
		self.actionQuit.setObjectName("actionQuit")
		self.actionModify_entry = QtWidgets.QAction(MainWindow)
		self.actionModify_entry.setObjectName("actionDelet_entry")
		self.actionDelete_entry = QtWidgets.QAction(MainWindow)
		self.actionDelete_entry.setObjectName("actionDelete_entry")

		
		self.actionAnnotated_from_reads_fastq = QtWidgets.QAction(MainWindow)
		self.actionAnnotated_from_reads_fastq.setObjectName("actionAnnotated_from_reads_fastq")
		self.actionFrom_genome_fasta = QtWidgets.QAction(MainWindow)
		self.actionFrom_genome_fasta.setObjectName("actionFrom_genome_fasta")
		self.actionFrom_NCBI_accessio_number = QtWidgets.QAction(MainWindow)
		self.actionFrom_NCBI_accessio_number.setObjectName("actionFrom_NCBI_accessio_number")
		self.actionBlank_entry = QtWidgets.QAction(MainWindow)
		self.actionBlank_entry.setObjectName("actionBlank_entry")
		self.actionTo_existing_entry = QtWidgets.QAction(MainWindow)
		self.actionTo_existing_entry.setObjectName("actionTo_existing_entry")
		self.actionProtein = QtWidgets.QAction(MainWindow)
		self.actionProtein.setObjectName("actionProtein")
		self.menuNew_entry_2.addAction(self.actionAnnotated_from_reads_fastq)
		self.menuNew_entry_2.addAction(self.actionFrom_genome_fasta)
		self.menuNew_entry_2.addAction(self.actionFrom_NCBI_accessio_number)
		self.menuNew_entry_2.addAction(self.actionBlank_entry)
		self.menuNew_sequence.addAction(self.actionTo_existing_entry)
		self.menuNew_sequence.addAction(self.actionProtein)
		self.menuNew_entry.addAction(self.menuNew_entry_2.menuAction())
		self.menuNew_entry.addAction(self.menuNew_sequence.menuAction())
		self.menuFile.addAction(self.actionUpdateDB)
		self.menuFile.addAction(self.actionSave_sequences)
		self.menuFile.addAction(self.actionSave_alignment)
		self.menuFile.addAction(self.actionSave_tree)
		self.menuFile.addAction(self.actionQuit)
		self.menuEdit.addAction(self.actionModify_entry)
		self.menuEdit.addAction(self.actionDelete_entry)
		self.menubar.addAction(self.menuFile.menuAction())
		self.menubar.addAction(self.menuEdit.menuAction())
		self.menubar.addAction(self.menuNew_entry.menuAction())
		#self.menubar.addAction(self.menuSync.menuAction())

		self.retranslateUi(MainWindow)
		QtCore.QMetaObject.connectSlotsByName(MainWindow)

		self.geneList.addItems(["IRS1","RL1","RL10","RL11","RL12","RL13","RL5A","RL6","RL8A","RL9A","TRS1","UL1","UL10","UL100","UL102","UL103","UL104","UL105","UL11","UL111A","UL112","UL114","UL115","UL116","UL117","UL119","UL120","UL121","UL122","UL123","UL124","UL128","UL13","UL130","UL131A","UL132","UL133","UL135","UL136","UL138","UL139","UL14","UL140","UL141","UL142","UL144","UL145","UL146","UL147","UL147A","UL148","UL148A","UL148B","UL148C","UL148D","UL150","UL150A","UL15A","UL16","UL17","UL18","UL19","UL2","UL20","UL21A","UL22A","UL23","UL24","UL25","UL26","UL27","UL29","UL30","UL30A","UL31","UL32","UL33","UL34","UL35","UL36","UL37","UL38","UL4","UL40","UL41A","UL42","UL43","UL44","UL45","UL46","UL47","UL48","UL48A","UL49","UL5","UL50","UL51","UL52","UL53","UL54","UL55","UL56","UL57","UL6","UL69","UL7","UL70","UL71","UL72","UL73","UL74","UL74A","UL75","UL76","UL77","UL78","UL79","UL8","UL80","UL80.5","UL82","UL83","UL84","UL85","UL86","UL87","UL88","UL89","UL9","UL91","UL92","UL93","UL94","UL95","UL96","UL97","UL98","UL99","US1","US10","US11","US12","US13","US14","US15","US16","US17","US18","US19","US2","US20","US21","US22","US23","US24","US26","US27","US28","US29","US3","US30","US31","US32","US33A","US34","US34A","US6","US7","US8","US9"])
		self.geneList.setSelectionMode(QtWidgets.QListWidget.MultiSelection)
		self.compartmentList.setSelectionMode(QtWidgets.QListWidget.MultiSelection)
		self.sampleDateList.setSelectionMode(QtWidgets.QListWidget.MultiSelection)
		self.patientTypeList.setSelectionMode(QtWidgets.QListWidget.MultiSelection)
		self.patientNumberList.setSelectionMode(QtWidgets.QListWidget.MultiSelection)
		self.menubar.setNativeMenuBar(False)
		self.actionQuit.triggered.connect(self.exitProgram)
		self.actionFrom_genome_fasta.triggered.connect(self.selectGenomeEntry)
		self.actionAnnotated_from_reads_fastq.triggered.connect(self.selectReadsEntry)
		self.actionBlank_entry.triggered.connect(self.createBlankEntry)
		self.actionModify_entry.triggered.connect(self.modifyEntry)
		
		self.actionDelete_entry.triggered.connect(self.deleteEntry)
		self.actionFrom_NCBI_accessio_number.triggered.connect(self.addEntryFromNCBI)
		self.actionSave_sequences.triggered.connect(self.saveSequencesToFile)
		self.actionSave_alignment.triggered.connect(self.saveAlignmentToFile)
		self.actionUpdateDB.triggered.connect(self.updateDB)
		self.actionUpdateDB.setEnabled(True)
		self.actionSave_tree.triggered.connect(self.saveTreeToFile)
		self.actionTo_existing_entry.triggered.connect(self.addNucleotideSequence)
		self.actionProtein.triggered.connect(self.addProteinSequence)
		self.patientTypeList.itemClicked.connect(self.addSampleIDs)
		self.patientNumberList.itemClicked.connect(self.addSampleDetails)
		self.sampleDateList.itemClicked.connect(self.collectSequences)
		self.compartmentList.itemClicked.connect(self.collectSequences)
		self.geneList.itemClicked.connect(self.collectSequences)
		self.addCDSButton.clicked.connect(self.addCDSsequences)
		self.addProteinsButton.clicked.connect(self.addProteinSequences)
		self.alignButton.clicked.connect(self.alignSequences)
		self.clearButton.clicked.connect(self.clearLogArea)
		self.treeButton.clicked.connect(self.createTree)
		self.cdsSeqs = {}
		self.proteinSeqs = {}

	


		#Check and report the sample types
		typesSet = set()
		for sampleID in os.listdir(installationDirectory+"/DB/"):
			if not sampleID[0] == ".":
				for filename in os.listdir(installationDirectory+"/DB/"+sampleID+"/"):
					if "_cds.fasta" in filename:
						typesSet.add((filename.split("-"))[0]) 
		typesList = list(typesSet)
		typesList = sorted(typesList)
		for a in typesList:
			self.patientTypeList.addItem(a)


	def updateDB(self):
		self.patientTypeList.clear()
		self.patientNumberList.clear()
		self.sampleDateList.clear()
		self.compartmentList.clear()
		typesSet = set()
		for sampleID in os.listdir(installationDirectory+"/DB/"):
			if not sampleID[0] == ".":
				for filename in os.listdir(installationDirectory+"/DB/"+sampleID+"/"):
					if "_cds.fasta" in filename:
						typesSet.add((filename.split("-"))[0]) 
		typesList = list(typesSet)
		typesList = sorted(typesList)
		for a in typesList:
			self.patientTypeList.addItem(a)

	def addEntryFromNCBI(self):
		os.system(installationDirectory+"/src/conda/bin/python "+installationDirectory+"/src/scripts/getNCBI.py "+installationDirectory+" &")




	def createTree(self):
		from Bio import Phylo
		import matplotlib

		if len(self.selectedSequencesArea.toPlainText()) < 10 or (self.selectedSequencesArea.toPlainText()).count(">") < 2:
			msg = QtWidgets.QMessageBox()
			msg.setText("At least two sequences must be available in the preview window" )
			msg.exec() 
			return
			
		outfile = open(".temp140875.fasta","w")
		outfile.write(self.selectedSequencesArea.toPlainText())
		outfile.close()

		outfile = open(".temp1408752.fasta","w")
		for seq_record in SeqIO.parse(".temp140875.fasta","fasta"):
			seqID = str(seq_record.id)
			seqSeq = str(seq_record.seq).replace("-","")
			outfile.write(">"+seqID+"\n"+seqSeq+"\n")
		outfile.close()
		os.system(installationDirectory+"/src/conda/bin/clustalo -i .temp1408752.fasta > .aln140875.fasta")
		os.system(installationDirectory+"/src/conda/bin/fasttree -nt .aln140875.fasta > .tree")
		tree = Phylo.read(".tree", "newick")
		tree.ladderize()  # Flip branches so deeper clades are displayed at top
		Phylo.draw(tree)
		infile = open(".tree")
		self.selectedSequencesArea.clear()
		while True:
			line = infile.readline().rstrip()
			if not line:
				break
			self.selectedSequencesArea.insertPlainText(line)
		infile.close()

		self.actionSave_sequences.setEnabled(False)
		self.actionSave_alignment.setEnabled(False)
		self.actionSave_tree.setEnabled(True)
		os.system("rm -rf .tree .aln140875.fasta .temp1408752.fasta .temp140875.fasta")
	
	def saveSequencesToFile(self):
		options = QFileDialog.Options()
		options |= QFileDialog.DontUseNativeDialog
		fileName, _ = QFileDialog.getSaveFileName(None,"Save sequences","","FastaFiles (*.fasta, *.fa)", options=options)
		if fileName.endswith(".fasta") == True or fileName.endswith(".fa"):
			pass
		else:
			fileName = fileName+".fasta"
		outfile = open(fileName,"w")
		outfile.write(self.selectedSequencesArea.toPlainText())
		outfile.close()

	def saveAlignmentToFile(self):
		options = QFileDialog.Options()
		options |= QFileDialog.DontUseNativeDialog
		fileName, _ = QFileDialog.getSaveFileName(None,"Save alignment","","Alignment file (*.aln)", options=options)
		if fileName.endswith(".aln") == True:
			pass
		else:
			fileName = fileName+".aln"
		outfile = open(fileName,"w")
		outfile.write(self.selectedSequencesArea.toPlainText())
		outfile.close()

	def saveTreeToFile(self):
		options = QFileDialog.Options()
		options |= QFileDialog.DontUseNativeDialog
		fileName, _ = QFileDialog.getSaveFileName(None,"Save tree","","Tree file (*.nwk)", options=options)
		if fileName.endswith(".nwk") == True:
			pass
		else:
			fileName = fileName+".nwk"
		outfile = open(fileName,"w")
		outfile.write(self.selectedSequencesArea.toPlainText())
		outfile.close()


	def clearLogArea(self):
		self.selectedSequencesArea.clear()
		self.actionSave_sequences.setEnabled(False)
		self.actionSave_alignment.setEnabled(False)
		self.actionSave_tree.setEnabled(False)



	def alignSequences(self):

		if len(self.selectedSequencesArea.toPlainText()) < 10 or (self.selectedSequencesArea.toPlainText()).count(">") < 2:
			msg = QtWidgets.QMessageBox()
			msg.setText("At least two sequences must be available in the preview window" )
			msg.exec() 
			return

		
		tempFile = open("temp140875.fasta","w")
		tempFile.write(self.selectedSequencesArea.toPlainText())
		tempFile.close()

		

		self.selectedSequencesArea.clear()
		self.selectedSequencesArea.insertPlainText("Aligning sequences with clustalx....")
		os.system(installationDirectory+"/src/conda/bin/clustalo -i temp140875.fasta > aln140875.fasta")
		os.system("rm temp140875.fasta")
		self.selectedSequencesArea.clear()
		tempFile = open("aln140875.fasta")
		while True:
			line = tempFile.readline()
			if not line:
				break
			self.selectedSequencesArea.insertPlainText(line)
		tempFile.close()
		
		os.system("rm aln140875.fasta")
		self.actionSave_sequences.setEnabled(False)
		self.actionSave_alignment.setEnabled(True)
		
		


	def addCDSsequences(self):
		if len(self.cdsSeqs) == 0:
			msg = QtWidgets.QMessageBox()
			msg.setText("No gene selected.\nPlease select at least one sample type/number/date/compartment/gene." )
			msg.exec() 
			return
		else:
			for sequence in self.cdsSeqs:
				self.selectedSequencesArea.insertPlainText(">"+sequence+"\n"+self.cdsSeqs[sequence]+"\n")
				self.actionSave_sequences.setEnabled(True)

	def addProteinSequences(self):
		if len(self.proteinSeqs) == 0:
			msg = QtWidgets.QMessageBox()
			msg.setText("No gene selected.\nPlease select at least one sample type/number/date/compartment/gene." )
			msg.exec() 
			return
		else:
			for sequence in self.proteinSeqs:
				self.selectedSequencesArea.insertPlainText(">"+sequence+"\n"+self.proteinSeqs[sequence]+"\n")
				self.actionSave_sequences.setEnabled(True)

	def collectSequences(self):
		selectedPatientIDs = []
		[selectedPatientIDs.append(x.data()) for x in self.patientNumberList.selectedIndexes()]

		selectedPatientTypes = []
		[selectedPatientTypes.append(x.data()) for x in self.patientTypeList.selectedIndexes()]	

		selectedDate = set()
		[selectedDate.add(x.data()) for x in self.sampleDateList.selectedIndexes()]	

		selectedCompartments = set()
		[selectedCompartments.add(x.data()) for x in self.compartmentList.selectedIndexes()]

		selectedGenes = []
		[selectedGenes.append(x.data()) for x in self.geneList.selectedIndexes()]	

		availableCDSSequences = {}
		availableProteinsSequences = {}
		if len(selectedPatientTypes)>0 and len(selectedPatientIDs)>0 and len(selectedDate)>0 and len(selectedCompartments)>0 and len(selectedGenes)>0:
			for t in selectedPatientTypes:
				for i in selectedPatientIDs:
					for d in selectedDate:
						for c in selectedCompartments:
							for gene in selectedGenes:
								prefix = t+"-"+i+"-"+d+"-"+c

								filename = installationDirectory+"/DB/"+i+"/"+t+"-"+i+"-"+d+"-"+c+"_cds.fasta"
								if os.path.isfile(filename) == True:
									for seq_record in SeqIO.parse(filename,"fasta"):
										if str(seq_record.id) == gene:
											if not gene+"_"+prefix in availableCDSSequences:
												if self.checkBox.isChecked() == True:
													if not "pseudo" in str(seq_record.description) and not "Pseudo" in str(seq_record.description):
														availableCDSSequences[gene+"_"+prefix] = str(seq_record.seq)
												else:
													availableCDSSequences[gene+"_"+prefix] = str(seq_record.seq)


								filename = installationDirectory+"/DB/"+i+"/"+t+"-"+i+"-"+d+"-"+c+"_proteins.fasta"
								if os.path.isfile(filename) == True:
									for seq_record in SeqIO.parse(filename,"fasta"):
										if str(seq_record.id) == gene:
											if not gene+"_"+prefix in availableProteinsSequences:
												availableProteinsSequences[gene+"_"+prefix] = str(seq_record.seq)

		self.numSeqText.setText(str(len(availableCDSSequences)))
		self.proteinSeqs = availableProteinsSequences
		self.cdsSeqs = availableCDSSequences



	def addSampleDetails(self):
		selectedPatientIDs = []
		[selectedPatientIDs.append(x.data()) for x in self.patientNumberList.selectedIndexes()]

		details2add = set()
		self.sampleDateList.clear()
		self.compartmentList.clear()
		self.geneList.clear()
		self.numSeqText.setText("")
		self.cdsSeqs = {}
		self.proteinSeqs = {}
		self.geneList.addItems(["IRS1","RL1","RL10","RL11","RL12","RL13","RL5A","RL6","RL8A","RL9A","TRS1","UL1","UL10","UL100","UL102","UL103","UL104","UL105","UL11","UL111A","UL112","UL114","UL115","UL116","UL117","UL119","UL120","UL121","UL122","UL123","UL124","UL128","UL13","UL130","UL131A","UL132","UL133","UL135","UL136","UL138","UL139","UL14","UL140","UL141","UL142","UL144","UL145","UL146","UL147","UL147A","UL148","UL148A","UL148B","UL148C","UL148D","UL150","UL150A","UL15A","UL16","UL17","UL18","UL19","UL2","UL20","UL21A","UL22A","UL23","UL24","UL25","UL26","UL27","UL29","UL30","UL30A","UL31","UL32","UL33","UL34","UL35","UL36","UL37","UL38","UL4","UL40","UL41A","UL42","UL43","UL44","UL45","UL46","UL47","UL48","UL48A","UL49","UL5","UL50","UL51","UL52","UL53","UL54","UL55","UL56","UL57","UL6","UL69","UL7","UL70","UL71","UL72","UL73","UL74","UL74A","UL75","UL76","UL77","UL78","UL79","UL8","UL80","UL80.5","UL82","UL83","UL84","UL85","UL86","UL87","UL88","UL89","UL9","UL91","UL92","UL93","UL94","UL95","UL96","UL97","UL98","UL99","US1","US10","US11","US12","US13","US14","US15","US16","US17","US18","US19","US2","US20","US21","US22","US23","US24","US26","US27","US28","US29","US3","US30","US31","US32","US33A","US34","US34A","US6","US7","US8","US9"])
		for sampleID in os.listdir(installationDirectory+"/DB/"):
			if not sampleID[0] == ".":
				for filename in os.listdir(installationDirectory+"/DB/"+sampleID+"/"):
					if "_cds.fasta" in filename:
						if (filename.split("-"))[1] in selectedPatientIDs:
							details2add.add(( (filename.split("-"))[2],((filename.split("-"))[3]).replace("_cds.fasta","")  ))
		
		for detail in details2add:
			self.sampleDateList.addItem(detail[0])
			self.compartmentList.addItem(detail[1])



	def addSampleIDs(self):
		selectedPatientTypes = []
		
		[selectedPatientTypes.append(x.data()) for x in self.patientTypeList.selectedIndexes()]
		#print(selectedPatientTypes)

		sample2add = set()
		self.patientNumberList.clear()
		self.sampleDateList.clear()
		self.compartmentList.clear()
		self.geneList.clear()
		self.numSeqText.setText("")
		self.cdsSeqs = {}
		self.proteinSeqs = {}
		self.geneList.addItems(["IRS1","RL1","RL10","RL11","RL12","RL13","RL5A","RL6","RL8A","RL9A","TRS1","UL1","UL10","UL100","UL102","UL103","UL104","UL105","UL11","UL111A","UL112","UL114","UL115","UL116","UL117","UL119","UL120","UL121","UL122","UL123","UL124","UL128","UL13","UL130","UL131A","UL132","UL133","UL135","UL136","UL138","UL139","UL14","UL140","UL141","UL142","UL144","UL145","UL146","UL147","UL147A","UL148","UL148A","UL148B","UL148C","UL148D","UL150","UL150A","UL15A","UL16","UL17","UL18","UL19","UL2","UL20","UL21A","UL22A","UL23","UL24","UL25","UL26","UL27","UL29","UL30","UL30A","UL31","UL32","UL33","UL34","UL35","UL36","UL37","UL38","UL4","UL40","UL41A","UL42","UL43","UL44","UL45","UL46","UL47","UL48","UL48A","UL49","UL5","UL50","UL51","UL52","UL53","UL54","UL55","UL56","UL57","UL6","UL69","UL7","UL70","UL71","UL72","UL73","UL74","UL74A","UL75","UL76","UL77","UL78","UL79","UL8","UL80","UL80.5","UL82","UL83","UL84","UL85","UL86","UL87","UL88","UL89","UL9","UL91","UL92","UL93","UL94","UL95","UL96","UL97","UL98","UL99","US1","US10","US11","US12","US13","US14","US15","US16","US17","US18","US19","US2","US20","US21","US22","US23","US24","US26","US27","US28","US29","US3","US30","US31","US32","US33A","US34","US34A","US6","US7","US8","US9"])
		for sampleID in os.listdir(installationDirectory+"/DB/"):
			if not sampleID[0] == ".":
				for filename in os.listdir(installationDirectory+"/DB/"+sampleID+"/"):
					if "_cds.fasta" in filename:
						if (filename.split("-"))[0] in selectedPatientTypes:
							sample2add.add((filename.split("-"))[1])

		self.patientNumberList.addItems(sorted(list(sample2add)))





	def exitProgram(self):
		exit()


	def addNucleotideSequence(self):
		os.system(installationDirectory+"/src/conda/bin/python "+installationDirectory+"/src/scripts/insertNucleotideSeq.py "+installationDirectory+" &")

	def addProteinSequence(self):
		os.system(installationDirectory+"/src/conda/bin/python "+installationDirectory+"/src/scripts/insertProteinSeq.py "+installationDirectory+" &")



	def selectGenomeEntry(self):

		os.system(installationDirectory+"/src/conda/bin/python "+installationDirectory+"/src/scripts/genomeEntry.py "+installationDirectory+" &") 

	def selectReadsEntry(self):
		os.system(installationDirectory+"/src/conda/bin/python "+installationDirectory+"/src/scripts/readsEntry.py "+installationDirectory+" &")


	def createBlankEntry(self):
		os.system(installationDirectory+"/src/conda/bin/python "+installationDirectory+"/src/scripts/blankEntry.py "+installationDirectory+" &")


	def modifyEntry(self):
		os.system(installationDirectory+"/src/conda/bin/python "+installationDirectory+"/src/scripts/modifyEntry.py "+installationDirectory+" &")
	
	def deleteEntry(self):
		os.system(installationDirectory+"/src/conda/bin/python "+installationDirectory+"/src/scripts/removeEntry.py "+installationDirectory+" &")
	
	def retranslateUi(self, MainWindow):
		_translate = QtCore.QCoreApplication.translate
		MainWindow.setWindowTitle(_translate("MainWindow", "SeReNA"))
		self.label.setText(_translate("MainWindow", "Sample type"))
		self.label_2.setText(_translate("MainWindow", "Sample number"))
		self.label_3.setText(_translate("MainWindow", "Sample date"))
		self.label_4.setText(_translate("MainWindow", "Compartment"))
		self.label_5.setText(_translate("MainWindow", "Gene"))
		self.label_7.setText(_translate("MainWindow", "Found sequences"))
		self.numSeqText.setHtml(_translate("MainWindow", "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.0//EN\" \"http://www.w3.org/TR/REC-html40/strict.dtd\">\n"
"<html><head><meta name=\"qrichtext\" content=\"1\" /><style type=\"text/css\">\n"
"p, li { white-space: pre-wrap; }\n"
"</style></head><body style=\" font-family:\'.AppleSystemUIFont\'; font-size:13pt; font-weight:400; font-style:normal;\">\n"
"<p align=\"right\" style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\">0</p></body></html>"))
		self.addCDSButton.setText(_translate("MainWindow", "Add CDS"))
		self.addProteinsButton.setText(_translate("MainWindow", "Add proteins"))
		self.treeButton.setText(_translate("MainWindow", "Tree"))
		self.alignButton.setText(_translate("MainWindow", "Align"))
		self.clearButton.setText(_translate("MainWindow", "Clear preview"))
		self.label_6.setText(_translate("MainWindow", "Preview"))
		self.checkBox.setText(_translate("MainWindow", "Exclude pseudo"))
		self.menuNew_entry.setTitle(_translate("MainWindow", "Insert"))
		self.menuNew_entry_2.setTitle(_translate("MainWindow", "New entry"))
		self.menuNew_sequence.setTitle(_translate("MainWindow", "New sequence"))
		#self.menuSync.setTitle(_translate("MainWindow", "Update DB"))
		self.menuFile.setTitle(_translate("MainWindow", "File"))
		self.menuEdit.setTitle(_translate("MainWindow", "Edit"))
		self.actionFrom_genome.setText(_translate("MainWindow", "From genome (fasta)"))
		self.actionFrom_reads.setText(_translate("MainWindow", "From reads (fastq)"))
		self.actionFrom_NCBI.setText(_translate("MainWindow", "From NCBI accession no."))
		self.actionSave_alignment.setText(_translate("MainWindow", "Save alignment (fasta)"))
		self.actionUpdateDB.setText(_translate("MainWindow", "Update DB"))
		self.actionSave_tree.setText(_translate("MainWindow", "Save tree"))
		self.actionSave_sequences.setText(_translate("MainWindow", "Save sequences (fasta)"))
		self.actionQuit.setText(_translate("MainWindow", "Quit"))
		self.actionModify_entry.setText(_translate("MainWindow", "Modify entry"))
		self.actionDelete_entry.setText(_translate("MainWindow","Delete entry"))
		self.actionAnnotated_from_reads_fastq.setText(_translate("MainWindow", "From reads (fastq)"))
		self.actionFrom_genome_fasta.setText(_translate("MainWindow", "From genome (fasta)"))
		self.actionFrom_NCBI_accessio_number.setText(_translate("MainWindow", "From NCBI accessio number "))
		self.actionBlank_entry.setText(_translate("MainWindow", "Blank entry"))
		self.actionTo_existing_entry.setText(_translate("MainWindow", "Nucleotide"))
		self.actionProtein.setText(_translate("MainWindow", "Protein"))
		

import sys 
if __name__ == "__main__":
	app = QtWidgets.QApplication(sys.argv)
	win = QtWidgets.QMainWindow()
	ui = Ui_MainWindow()
	ui.setupUi(win)

	win.show()
	sys.exit(app.exec())