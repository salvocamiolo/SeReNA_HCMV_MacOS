# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'blankEntry.ui'
#
# Created by: PyQt5 UI code generator 5.9.2
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets
import os

class Ui_Form(object):
	def setupUi(self, Form, installationDirectory):
		Form.setObjectName("Form")
		Form.resize(304, 189)
		self.compartment_lineEdit = QtWidgets.QLineEdit(Form)
		self.compartment_lineEdit.setGeometry(QtCore.QRect(170, 100, 113, 21))
		self.compartment_lineEdit.setObjectName("compartment_lineEdit")
		self.sampleDate_lineEdit = QtWidgets.QLineEdit(Form)
		self.sampleDate_lineEdit.setGeometry(QtCore.QRect(20, 100, 113, 21))
		self.sampleDate_lineEdit.setObjectName("sampleDate_lineEdit")
		self.label_11 = QtWidgets.QLabel(Form)
		self.label_11.setGeometry(QtCore.QRect(170, 80, 91, 16))
		self.label_11.setObjectName("label_11")
		self.label_10 = QtWidgets.QLabel(Form)
		self.label_10.setGeometry(QtCore.QRect(20, 80, 91, 16))
		self.label_10.setObjectName("label_10")
		self.sampleName_lineEdit = QtWidgets.QLineEdit(Form)
		self.sampleName_lineEdit.setGeometry(QtCore.QRect(20, 40, 113, 21))
		self.sampleName_lineEdit.setObjectName("sampleName_lineEdit")
		self.sampleType_lineEdit = QtWidgets.QLineEdit(Form)
		self.sampleType_lineEdit.setGeometry(QtCore.QRect(170, 40, 113, 21))
		self.sampleType_lineEdit.setObjectName("sampleType_lineEdit")
		self.label_9 = QtWidgets.QLabel(Form)
		self.label_9.setGeometry(QtCore.QRect(170, 20, 91, 16))
		self.label_9.setObjectName("label_9")
		self.label_8 = QtWidgets.QLabel(Form)
		self.label_8.setGeometry(QtCore.QRect(20, 20, 91, 16))
		self.label_8.setObjectName("label_8")
		self.pushButton = QtWidgets.QPushButton(Form)
		self.pushButton.setGeometry(QtCore.QRect(70, 150, 181, 32))
		self.pushButton.setObjectName("pushButton")

		self.retranslateUi(Form)
		QtCore.QMetaObject.connectSlotsByName(Form)

		self.pushButton.clicked.connect(self.createNewEntry)

	def createNewEntry(self):
		#to implement the quality check of the entry
		if self.sampleType_lineEdit.text() == "":
			msg = QtWidgets.QMessageBox()
			msg.setText("A sample type must be specified!")
			msg.exec() 
			return

		if self.sampleName_lineEdit.text() == "":
			msg = QtWidgets.QMessageBox()
			msg.setText("A sample ID must be specified!")
			msg.exec() 
			return

		if self.sampleDate_lineEdit.text() == "":
			msg = QtWidgets.QMessageBox()
			msg.setText("A valid date must be specified!")
			msg.exec() 
			return

		if self.compartment_lineEdit.text() == "":
			msg = QtWidgets.QMessageBox()
			msg.setText("The compartment must be specified!")
			msg.exec() 
			return

		prefix = self.sampleType_lineEdit.text()+"-"+self.sampleName_lineEdit.text()+"-"+self.sampleDate_lineEdit.text()+"-"+self.compartment_lineEdit.text()
		outputFolder = installationDirectory+"/DB/"+self.sampleName_lineEdit.text()
		if os.path.isdir(outputFolder):
			msg = QtWidgets.QMessageBox()
			msg.setText("The sample ID is alread present in the database! Please choose another name")
			msg.exec() 
			return
		else:
			os.system("mkdir "+outputFolder)
			os.system("touch "+outputFolder+"/"+prefix+"_cds.fasta")
			os.system("touch "+outputFolder+"/"+prefix+"_proteins.fasta")
			msg = QtWidgets.QMessageBox()
			msg.setText("New query ID successfully created!")
			msg.exec() 
			self.sampleType_lineEdit.setText("")
			self.sampleName_lineEdit.setText("")
			self.sampleDate_lineEdit.setText("")
			self.compartment_lineEdit.setText("")
			return

			




	def retranslateUi(self, Form):
		_translate = QtCore.QCoreApplication.translate
		Form.setWindowTitle(_translate("Form", "New blank entry"))
		self.label_11.setText(_translate("Form", "Compartment"))
		self.label_10.setText(_translate("Form", "Sample date"))
		self.label_9.setText(_translate("Form", "Sample type"))
		self.label_8.setText(_translate("Form", "Sample name"))
		self.pushButton.setText(_translate("Form", "Create new blank entry"))


import sys
if __name__ == "__main__":
	App = QtWidgets.QApplication(sys.argv)
	Form = QtWidgets.QWidget()
	ui = Ui_Form()
	installationDirectory = sys.argv[1]
	ui.setupUi(Form, installationDirectory)
	Form.show()
	sys.exit(App.exec())