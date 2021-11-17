# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'modifyEntry.ui'
#
# Created by: PyQt5 UI code generator 5.9.2
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets
import os

class Ui_Form(object):
	def setupUi(self, Form, installationDirectory):
		Form.setObjectName("Form")
		Form.resize(379, 523)
		self.updateButton = QtWidgets.QPushButton(Form)
		self.updateButton.setGeometry(QtCore.QRect(230, 480, 131, 32))
		self.updateButton.setObjectName("updateButton")
		self.groupBox = QtWidgets.QGroupBox(Form)
		self.groupBox.setGeometry(QtCore.QRect(20, 20, 341, 201))
		self.groupBox.setObjectName("groupBox")
		self.compartment_lineEdit = QtWidgets.QLineEdit(self.groupBox)
		self.compartment_lineEdit.setGeometry(QtCore.QRect(180, 140, 113, 21))
		self.compartment_lineEdit.setObjectName("compartment_lineEdit")
		self.label_2 = QtWidgets.QLabel(self.groupBox)
		self.label_2.setGeometry(QtCore.QRect(20, 50, 71, 16))
		self.label_2.setObjectName("label_2")
		self.sampleID_lineEdit = QtWidgets.QLineEdit(self.groupBox)
		self.sampleID_lineEdit.setGeometry(QtCore.QRect(20, 70, 113, 21))
		self.sampleID_lineEdit.setObjectName("sampleID_lineEdit")
		self.label_5 = QtWidgets.QLabel(self.groupBox)
		self.label_5.setGeometry(QtCore.QRect(180, 120, 101, 16))
		self.label_5.setObjectName("label_5")
		self.sampleType_lineEdit = QtWidgets.QLineEdit(self.groupBox)
		self.sampleType_lineEdit.setGeometry(QtCore.QRect(180, 70, 113, 21))
		self.sampleType_lineEdit.setObjectName("sampleType_lineEdit")
		self.label_4 = QtWidgets.QLabel(self.groupBox)
		self.label_4.setGeometry(QtCore.QRect(20, 120, 81, 16))
		self.label_4.setObjectName("label_4")
		self.label_3 = QtWidgets.QLabel(self.groupBox)
		self.label_3.setGeometry(QtCore.QRect(180, 50, 101, 16))
		self.label_3.setObjectName("label_3")
		self.sampleDate_lineEdit = QtWidgets.QLineEdit(self.groupBox)
		self.sampleDate_lineEdit.setGeometry(QtCore.QRect(20, 140, 113, 21))
		self.sampleDate_lineEdit.setObjectName("sampleDate_lineEdit")
		self.groupBox_2 = QtWidgets.QGroupBox(Form)
		self.groupBox_2.setGeometry(QtCore.QRect(20, 250, 341, 201))
		self.groupBox_2.setObjectName("groupBox_2")
		self.compartment_lineEdit_2 = QtWidgets.QLineEdit(self.groupBox_2)
		self.compartment_lineEdit_2.setGeometry(QtCore.QRect(180, 140, 113, 21))
		self.compartment_lineEdit_2.setObjectName("compartment_lineEdit_2")
		self.label_6 = QtWidgets.QLabel(self.groupBox_2)
		self.label_6.setGeometry(QtCore.QRect(20, 50, 71, 16))
		self.label_6.setObjectName("label_6")
		self.sampleID_lineEdit_2 = QtWidgets.QLineEdit(self.groupBox_2)
		self.sampleID_lineEdit_2.setGeometry(QtCore.QRect(20, 70, 113, 21))
		self.sampleID_lineEdit_2.setObjectName("sampleID_lineEdit_2")
		self.label_7 = QtWidgets.QLabel(self.groupBox_2)
		self.label_7.setGeometry(QtCore.QRect(180, 120, 101, 16))
		self.label_7.setObjectName("label_7")
		self.sampleType_lineEdit_2 = QtWidgets.QLineEdit(self.groupBox_2)
		self.sampleType_lineEdit_2.setGeometry(QtCore.QRect(180, 70, 113, 21))
		self.sampleType_lineEdit_2.setObjectName("sampleType_lineEdit_2")
		self.label_8 = QtWidgets.QLabel(self.groupBox_2)
		self.label_8.setGeometry(QtCore.QRect(20, 120, 81, 16))
		self.label_8.setObjectName("label_8")
		self.label_9 = QtWidgets.QLabel(self.groupBox_2)
		self.label_9.setGeometry(QtCore.QRect(180, 50, 101, 16))
		self.label_9.setObjectName("label_9")
		self.sampleDate_lineEdit_2 = QtWidgets.QLineEdit(self.groupBox_2)
		self.sampleDate_lineEdit_2.setGeometry(QtCore.QRect(20, 140, 113, 21))
		self.sampleDate_lineEdit_2.setObjectName("sampleDate_lineEdit_2")

		self.retranslateUi(Form)
		QtCore.QMetaObject.connectSlotsByName(Form)


		self.updateButton.clicked.connect(self.updateRecord)

	def updateRecord(self):

		if len(self.sampleID_lineEdit.text()) == 0 or len(self.compartment_lineEdit.text()) == 0 or len(self.sampleType_lineEdit.text()) == 0 or len(self.sampleDate_lineEdit.text()) == 0 or len(self.sampleID_lineEdit_2.text()) == 0 or len(self.compartment_lineEdit_2.text()) == 0 or len(self.sampleType_lineEdit_2.text()) == 0 or len(self.sampleDate_lineEdit_2.text()) == 0:
			msg = QtWidgets.QMessageBox()
			msg.setText("Please fill all the fields to modify an entry in the DB")
			msg.exec() 
			return

		oldRecordPrefix = installationDirectory+"/DB/"+self.sampleID_lineEdit.text()+"/"+self.sampleType_lineEdit.text()+"-"+self.sampleID_lineEdit.text()+"-"+self.sampleDate_lineEdit.text()+"-"+self.compartment_lineEdit.text()
		oldRecordFolder = installationDirectory+"/DB/"+self.sampleID_lineEdit.text()

		if os.path.isfile(oldRecordPrefix+"_cds.fasta")==True or os.path.isfile(oldRecordPrefix+"_proteins.fasta")==True:
			newRecordFolder = installationDirectory+"/DB/"+self.sampleID_lineEdit_2.text()
			newRecordPrefix = installationDirectory+"/DB/"+self.sampleID_lineEdit_2.text()+"/"+self.sampleType_lineEdit_2.text()+"-"+self.sampleID_lineEdit_2.text()+"-"+self.sampleDate_lineEdit_2.text()+"-"+self.compartment_lineEdit_2.text()
			if not oldRecordFolder == newRecordFolder:
				os.system("mv "+oldRecordFolder+" "+newRecordFolder)

			if os.path.isfile(oldRecordPrefix+"_cds.fasta")==True:
				os.system("mv "+oldRecordPrefix+"_cds.fasta "+newRecordPrefix+"_cds.fasta")
			
			if os.path.isfile(oldRecordPrefix+"_proteins.fasta")==True:
				os.system("mv "+oldRecordPrefix+"_proteins.fasta "+newRecordPrefix+"_proteins.fasta")

			msg = QtWidgets.QMessageBox()
			msg.setText("The record was successfully updated!")
			msg.exec()
			return
		
		else:
			msg = QtWidgets.QMessageBox()
			msg.setText("The record you entered does not seem to exist!")
			msg.exec() 
			return




	def retranslateUi(self, Form):
		_translate = QtCore.QCoreApplication.translate
		Form.setWindowTitle(_translate("Form", "Modify entry"))
		self.updateButton.setText(_translate("Form", "Update record"))
		self.groupBox.setTitle(_translate("Form", "Old metadata"))
		self.label_2.setText(_translate("Form", "Sample ID"))
		self.label_5.setText(_translate("Form", "Compartment"))
		self.label_4.setText(_translate("Form", "Sample date"))
		self.label_3.setText(_translate("Form", "Sample type"))
		self.groupBox_2.setTitle(_translate("Form", "New metadata"))
		self.label_6.setText(_translate("Form", "Sample ID"))
		self.label_7.setText(_translate("Form", "Compartment"))
		self.label_8.setText(_translate("Form", "Sample date"))
		self.label_9.setText(_translate("Form", "Sample type"))

import sys
if __name__ == "__main__":
	App = QtWidgets.QApplication(sys.argv)
	Form = QtWidgets.QWidget()
	ui = Ui_Form()
	installationDirectory = sys.argv[1]
	ui.setupUi(Form, installationDirectory)
	Form.show()
	sys.exit(App.exec())