# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'insertNucleotideSeq.ui'
#
# Created by: PyQt5 UI code generator 5.9.2
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets
import os
from Bio import SeqIO

class Ui_Form(object):
	def setupUi(self, Form, installationDirectory):
		Form.setObjectName("Form")
		Form.resize(635, 361)
		self.label_2 = QtWidgets.QLabel(Form)
		self.label_2.setGeometry(QtCore.QRect(10, 20, 71, 16))
		self.label_2.setObjectName("label_2")
		self.sampleID_lineEdit = QtWidgets.QLineEdit(Form)
		self.sampleID_lineEdit.setGeometry(QtCore.QRect(10, 40, 113, 21))
		self.sampleID_lineEdit.setObjectName("sampleID_lineEdit")
		self.label_5 = QtWidgets.QLabel(Form)
		self.label_5.setGeometry(QtCore.QRect(510, 20, 101, 16))
		self.label_5.setObjectName("label_5")
		self.sampleType_lineEdit = QtWidgets.QLineEdit(Form)
		self.sampleType_lineEdit.setGeometry(QtCore.QRect(180, 40, 113, 21))
		self.sampleType_lineEdit.setObjectName("sampleType_lineEdit")
		self.label_4 = QtWidgets.QLabel(Form)
		self.label_4.setGeometry(QtCore.QRect(350, 20, 81, 16))
		self.label_4.setObjectName("label_4")
		self.label_3 = QtWidgets.QLabel(Form)
		self.label_3.setGeometry(QtCore.QRect(180, 20, 101, 16))
		self.label_3.setObjectName("label_3")
		self.sampleDate_lineEdit = QtWidgets.QLineEdit(Form)
		self.sampleDate_lineEdit.setGeometry(QtCore.QRect(350, 40, 113, 21))
		self.sampleDate_lineEdit.setObjectName("sampleDate_lineEdit")
		self.compartment_lineEdit = QtWidgets.QLineEdit(Form)
		self.compartment_lineEdit.setGeometry(QtCore.QRect(510, 40, 113, 21))
		self.compartment_lineEdit.setObjectName("compartment_lineEdit")
		self.label_6 = QtWidgets.QLabel(Form)
		self.label_6.setGeometry(QtCore.QRect(10, 140, 151, 16))
		self.label_6.setObjectName("label_6")
		self.textEdit = QtWidgets.QTextEdit(Form)
		self.textEdit.setGeometry(QtCore.QRect(10, 160, 611, 151))
		self.textEdit.setObjectName("textEdit")
		self.pushButton = QtWidgets.QPushButton(Form)
		self.pushButton.setGeometry(QtCore.QRect(470, 320, 151, 32))
		self.pushButton.setObjectName("pushButton")
		self.label_7 = QtWidgets.QLabel(Form)
		self.label_7.setGeometry(QtCore.QRect(10, 80, 71, 16))
		self.label_7.setObjectName("label_7")
		self.comboBox = QtWidgets.QComboBox(Form)
		self.comboBox.setGeometry(QtCore.QRect(10, 100, 111, 26))
		self.comboBox.setObjectName("comboBox")
		self.comboBox.addItem("")

		self.retranslateUi(Form)
		QtCore.QMetaObject.connectSlotsByName(Form)


		self.pushButton.clicked.connect(self.addSequence)
		self.comboBox.addItems(["IRS1","RL1","RL10","RL11","RL12","RL13","RL5A","RL6","RL8A","RL9A","TRS1","UL1","UL10","UL100","UL102","UL103","UL104","UL105","UL11","UL111A","UL112","UL114","UL115","UL116","UL117","UL119","UL120","UL121","UL122","UL123","UL124","UL128","UL13","UL130","UL131A","UL132","UL133","UL135","UL136","UL138","UL139","UL14","UL140","UL141","UL142","UL144","UL145","UL146","UL147","UL147A","UL148","UL148A","UL148B","UL148C","UL148D","UL150","UL150A","UL15A","UL16","UL17","UL18","UL19","UL2","UL20","UL21A","UL22A","UL23","UL24","UL25","UL26","UL27","UL29","UL30","UL30A","UL31","UL32","UL33","UL34","UL35","UL36","UL37","UL38","UL4","UL40","UL41A","UL42","UL43","UL44","UL45","UL46","UL47","UL48","UL48A","UL49","UL5","UL50","UL51","UL52","UL53","UL54","UL55","UL56","UL57","UL6","UL69","UL7","UL70","UL71","UL72","UL73","UL74","UL74A","UL75","UL76","UL77","UL78","UL79","UL8","UL80","UL80.5","UL82","UL83","UL84","UL85","UL86","UL87","UL88","UL89","UL9","UL91","UL92","UL93","UL94","UL95","UL96","UL97","UL98","UL99","US1","US10","US11","US12","US13","US14","US15","US16","US17","US18","US19","US2","US20","US21","US22","US23","US24","US26","US27","US28","US29","US3","US30","US31","US32","US33A","US34","US34A","US6","US7","US8","US9"])

	def addSequence(self):
		prefix = self.sampleType_lineEdit.text()+"-"+self.sampleID_lineEdit.text()+"-"+self.sampleDate_lineEdit.text()+"-"+self.compartment_lineEdit.text()
		
		if os.path.isdir(installationDirectory+"/DB/"+self.sampleID_lineEdit.text()) == False:
			msg = QtWidgets.QMessageBox()
			msg.setText("The sample ID is not present in the database")
			msg.exec() 
			return
		if os.path.isfile(installationDirectory+"/DB/"+self.sampleID_lineEdit.text()+"/"+prefix+"_cds.fasta") == False:
			msg = QtWidgets.QMessageBox()
			msg.setText("Corrupted entry! The CDS file is not present in the specified sample")
			msg.exec() 
			return

		if str(self.comboBox.currentText()) == "Select":
			msg = QtWidgets.QMessageBox()
			msg.setText("Please select a gene name")
			msg.exec() 
			return

		if str(self.textEdit.toPlainText()) == "":
			msg = QtWidgets.QMessageBox()
			msg.setText("Sequence missing!")
			msg.exec() 
			return

		#check sequence
		sequence = (str(self.textEdit.toPlainText())).replace("\n","")
		admittedNucleotides = ['A','C','G','T','N']
		position = 0
		for a in sequence.upper():
			position+=1
			if not a in admittedNucleotides:
				msg = QtWidgets.QMessageBox()
				msg.setText("Found irregular nucleotide "+a+" at position "+str(position) )
				msg.exec() 
				return

		#check if the sequence is present in the database
		geneToInsert = str(self.comboBox.currentText())	
		geneExists = False
		foundSequence = ""
		sequences = {}
		for seq_record in SeqIO.parse(installationDirectory+"/DB/"+self.sampleID_lineEdit.text()+"/"+prefix+"_cds.fasta","fasta"):
			if str(seq_record.id) == geneToInsert:
				geneExists = True
				foundSequence = str(seq_record.seq)
				foundDescription = str(seq_record.description)
			if not str(seq_record.id) in sequences:
				sequences[str(seq_record.id)] = seq_record
		if geneExists == True:
			"""qm = QtWidgets.QMessageBox
			qm.question(self,'', "There is already a sequence for this gene. Do you want to overwrite it?", qm.Yes | qm.No)"""
			msgBox = QtWidgets.QMessageBox()
			msgBox.setText("A sequence for the gene is already present in DB")
			msgBox.setInformativeText("Do you want to overwrite it?")
			msgBox.setStandardButtons(QtWidgets.QMessageBox.Yes | QtWidgets.QMessageBox.No)
			msgBox.setDefaultButton(QtWidgets.QMessageBox.No)
			msgBox.setDetailedText("The following sequence for gene "+geneToInsert+" was found:\n"+foundSequence)
			ret = msgBox.exec_()

			if ret == QtWidgets.QMessageBox.No:
				return
			else:
				outfile = open(installationDirectory+"/DB/"+self.sampleID_lineEdit.text()+"/"+prefix+"_cds.fasta","w")
				for gene in sequences:
					if not gene == geneToInsert:
						SeqIO.write(sequences[gene],outfile,"fasta")
					else:
						outfile.write(">"+foundDescription+"\n"+(self.textEdit.toPlainText()).replace("\n","")+"\n")
				msg = QtWidgets.QMessageBox()
				msg.setText("Sequence has been added to DB")
				msg.exec() 
				self.textEdit.clear()
				self.sampleID_lineEdit.setText("")
				self.sampleType_lineEdit.setText("")
				self.sampleDate_lineEdit.setText("")
				self.compartment_lineEdit.setText("")
				self.comboBox.setCurrentIndex(0)
		else:
			outfile = open(installationDirectory+"/DB/"+self.sampleID_lineEdit.text()+"/"+prefix+"_cds.fasta","w")
			for gene in sequences:
				SeqIO.write(sequences[gene],outfile,"fasta")
			
			outfile.write(">"+geneToInsert+"\n"+(self.textEdit.toPlainText()).replace("\n","")+"\n")
			msg = QtWidgets.QMessageBox()
			msg.setText("Sequence has been added to DB")
			msg.exec() 
			self.textEdit.clear()
			self.sampleID_lineEdit.setText("")
			self.sampleType_lineEdit.setText("")
			self.sampleDate_lineEdit.setText("")
			self.compartment_lineEdit.setText("")
			self.comboBox.setCurrentIndex(0)
		outfile.close()






	def retranslateUi(self, Form):
		_translate = QtCore.QCoreApplication.translate
		Form.setWindowTitle(_translate("Form", "Insert nucleotide sequence"))
		self.label_2.setText(_translate("Form", "Sample ID"))
		self.label_5.setText(_translate("Form", "Compartment"))
		self.label_4.setText(_translate("Form", "Sample date"))
		self.label_3.setText(_translate("Form", "Sample type"))
		self.label_6.setText(_translate("Form", "Nucleotide sequence"))
		self.pushButton.setText(_translate("Form", "Add sequence to DB"))
		self.label_7.setText(_translate("Form", "Gene"))
		self.comboBox.setItemText(0, _translate("Form", "Select"))


import sys
if __name__ == "__main__":
	App = QtWidgets.QApplication(sys.argv)
	Form = QtWidgets.QWidget()
	ui = Ui_Form()
	installationDirectory = sys.argv[1]
	ui.setupUi(Form, installationDirectory)
	Form.show()
	sys.exit(App.exec())