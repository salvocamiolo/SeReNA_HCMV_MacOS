# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'getNCBI.ui'
#
# Created by: PyQt5 UI code generator 5.9.2
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets
import os
from Bio import SeqIO
from Bio.Seq import Seq


class Ui_Form(object):
	def setupUi(self, Form, installationDirectory):
		Form.setObjectName("Form")
		Form.resize(337, 251)
		self.accessionNumber_lineEdit = QtWidgets.QLineEdit(Form)
		self.accessionNumber_lineEdit.setGeometry(QtCore.QRect(10, 30, 171, 21))
		self.accessionNumber_lineEdit.setObjectName("accessionNumber_lineEdit")
		self.label = QtWidgets.QLabel(Form)
		self.label.setGeometry(QtCore.QRect(10, 10, 181, 16))
		self.label.setObjectName("label")
		self.sampleDate_lineEdit = QtWidgets.QLineEdit(Form)
		self.sampleDate_lineEdit.setGeometry(QtCore.QRect(10, 150, 171, 21))
		self.sampleDate_lineEdit.setObjectName("sampleDate_lineEdit")
		self.label_2 = QtWidgets.QLabel(Form)
		self.label_2.setGeometry(QtCore.QRect(10, 130, 181, 16))
		self.label_2.setObjectName("label_2")
		self.compartment_lineEdit = QtWidgets.QLineEdit(Form)
		self.compartment_lineEdit.setGeometry(QtCore.QRect(10, 210, 171, 21))
		self.compartment_lineEdit.setObjectName("compartment_lineEdit")
		self.label_3 = QtWidgets.QLabel(Form)
		self.label_3.setGeometry(QtCore.QRect(10, 190, 181, 16))
		self.label_3.setObjectName("label_3")
		self.sampleType_lineEdit = QtWidgets.QLineEdit(Form)
		self.sampleType_lineEdit.setGeometry(QtCore.QRect(10, 90, 171, 21))
		self.sampleType_lineEdit.setObjectName("sampleType_lineEdit")
		self.label_4 = QtWidgets.QLabel(Form)
		self.label_4.setGeometry(QtCore.QRect(10, 70, 181, 16))
		self.label_4.setObjectName("label_4")
		self.pushButton = QtWidgets.QPushButton(Form)
		self.pushButton.setGeometry(QtCore.QRect(210, 206, 113, 32))
		self.pushButton.setObjectName("pushButton")

		self.retranslateUi(Form)
		QtCore.QMetaObject.connectSlotsByName(Form)

		self.pushButton.clicked.connect(self.getEntry)

	def getEntry(self):
		if self.accessionNumber_lineEdit.text() == "":
			msg = QtWidgets.QMessageBox()
			msg.setText("Please add an accession number...")
			msg.exec() 
			return

		if len(self.sampleDate_lineEdit.text()) == 0 or len(self.sampleType_lineEdit.text())==0 or len(self.compartment_lineEdit.text()) == 0:
			msg = QtWidgets.QMessageBox()
			msg.setText("Sample type/date/compartment are mandatory fields")
			msg.exec() 
			return
		
		os.system(installationDirectory+"/src/conda/bin/esearch -db nucleotide -query "+self.accessionNumber_lineEdit.text()+" | "+installationDirectory+"/src/conda/bin/efetch -format fasta_cds_na > .temp140875.fasta")
		prefix = self.sampleType_lineEdit.text()+"-"+self.accessionNumber_lineEdit.text()+"-"+self.sampleDate_lineEdit.text()+self.compartment_lineEdit.text()

		sequences = {}
		geneNames = ["IRS1","RL1","RL10","RL11","RL12","RL13","RL5A","RL6","RL8A","RL9A","TRS1","UL1","UL10","UL100","UL102","UL103","UL104","UL105","UL11","UL111A","UL112","UL114","UL115","UL116","UL117","UL119","UL120","UL121","UL122","UL123","UL124","UL128","UL13","UL130","UL131A","UL132","UL133","UL135","UL136","UL138","UL139","UL14","UL140","UL141","UL142","UL144","UL145","UL146","UL147","UL147A","UL148","UL148A","UL148B","UL148C","UL148D","UL150","UL150A","UL15A","UL16","UL17","UL18","UL19","UL2","UL20","UL21A","UL22A","UL23","UL24","UL25","UL26","UL27","UL29","UL30","UL30A","UL31","UL32","UL33","UL34","UL35","UL36","UL37","UL38","UL4","UL40","UL41A","UL42","UL43","UL44","UL45","UL46","UL47","UL48","UL48A","UL49","UL5","UL50","UL51","UL52","UL53","UL54","UL55","UL56","UL57","UL6","UL69","UL7","UL70","UL71","UL72","UL73","UL74","UL74A","UL75","UL76","UL77","UL78","UL79","UL8","UL80","UL80.5","UL82","UL83","UL84","UL85","UL86","UL87","UL88","UL89","UL9","UL91","UL92","UL93","UL94","UL95","UL96","UL97","UL98","UL99","US1","US10","US11","US12","US13","US14","US15","US16","US17","US18","US19","US2","US20","US21","US22","US23","US24","US26","US27","US28","US29","US3","US30","US31","US32","US33A","US34","US34A","US6","US7","US8","US9"]
		foundGenes = 0
		outfile = open(".tempCDS_140875.fasta","w")
		outfile2 = open(".tempProteins_140875.fasta","w")
		for seq_record in SeqIO.parse(".temp140875.fasta","fasta"):
			geneDescription = str(seq_record.description)
			geneID = ((geneDescription.split("[gene=")[1]).split("]"))[0]
			geneSeq = str(seq_record.seq)
			protSeq = str(Seq(geneSeq).translate())
			if geneID in geneNames:
				foundGenes+=1
				outfile.write(">"+geneID+"\n"+geneSeq+"\n")
			if "*" in protSeq[:-1]:
				outfile2.write(">"+geneID+" _pseudo\n"+protSeq+"\n")
			else:
				outfile2.write(">"+geneID+"\n"+protSeq+"\n")
			

		outfile.close()
		outfile2.close()
		if (os.stat(".tempCDS_140875.fasta").st_size == 0) == False:
			os.system("mkdir "+installationDirectory+"/DB/"+self.accessionNumber_lineEdit.text())
			os.system("mv .tempCDS_140875.fasta "+installationDirectory+"/DB/"+self.accessionNumber_lineEdit.text()+"/"+prefix+"_cds.fasta")
			os.system("mv .tempProteins_140875.fasta "+installationDirectory+"/DB/"+self.accessionNumber_lineEdit.text()+"/"+prefix+"_proteins.fasta")
		msg = QtWidgets.QMessageBox()

		msg.setText("Entry successfully created!\n"+str(foundGenes)+" CDS/proteins have been retrieved.")
		msg.exec() 

		


	def retranslateUi(self, Form):
		_translate = QtCore.QCoreApplication.translate
		Form.setWindowTitle(_translate("Form", "Import from NCBI"))
		self.label.setText(_translate("Form", "GenBank accession number"))
		self.sampleDate_lineEdit.setText(_translate("Form", "XXXXXX"))
		self.label_2.setText(_translate("Form", "Sample date"))
		self.compartment_lineEdit.setText(_translate("Form", "X"))
		self.label_3.setText(_translate("Form", "Sample compartment"))
		self.sampleType_lineEdit.setText(_translate("Form", "NCBI"))
		self.label_4.setText(_translate("Form", "Sample type"))
		self.pushButton.setText(_translate("Form", "Get record"))

import sys
if __name__ == "__main__":
	import sys
	app = QtWidgets.QApplication(sys.argv)
	Form = QtWidgets.QWidget()
	installationDirectory = sys.argv[1]
	ui = Ui_Form()
	ui.setupUi(Form,installationDirectory)
	Form.show()
	sys.exit(app.exec_())