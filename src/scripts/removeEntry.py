# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'removeEntry.ui'
#
# Created by: PyQt5 UI code generator 5.9.2
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets
import os

class Ui_Form(object):
	def setupUi(self, Form, installationDirectory):
		Form.setObjectName("Form")
		Form.resize(365, 314)
		self.label = QtWidgets.QLabel(Form)
		self.label.setGeometry(QtCore.QRect(10, 20, 181, 16))
		self.label.setObjectName("label")
		self.remove_pushButton = QtWidgets.QPushButton(Form)
		self.remove_pushButton.setGeometry(QtCore.QRect(180, 30, 181, 32))
		self.remove_pushButton.setObjectName("remove_pushButton")
		self.pushButton = QtWidgets.QPushButton(Form)
		self.pushButton.setGeometry(QtCore.QRect(180, 60, 181, 32))
		self.pushButton.setObjectName("pushButton")
		self.listWidget = QtWidgets.QListWidget(Form)
		self.listWidget.setGeometry(QtCore.QRect(10, 40, 161, 261))
		self.listWidget.setObjectName("listWidget")
		self.listWidget.setSelectionMode(QtWidgets.QListWidget.MultiSelection)
		self.retranslateUi(Form)
		QtCore.QMetaObject.connectSlotsByName(Form)

		patientList = []
		for folderName in os.listdir(installationDirectory+"/DB/"):
			if not folderName[0] == ".":
				patientList.append(folderName)
		patientList = sorted(patientList)
		for patient in patientList: 
			self.listWidget.addItem(patient)

		self.pushButton.clicked.connect(self.exitWindow)
		self.remove_pushButton.clicked.connect(self.removeEntry)

	def removeEntry(self):
		entry2remove = []
		[entry2remove.append(x.data()) for x in self.listWidget.selectedIndexes()]	
		
		removedString = ""
		for patient in entry2remove:
			removedString+= patient+"\n"
		msgBox = QtWidgets.QMessageBox()
		msgBox.setText("The following entries will be permanently deleted from the DB:\n\n"+removedString+"\nDo you wish to continue?")
		msgBox.setStandardButtons(QtWidgets.QMessageBox.Yes | QtWidgets.QMessageBox.No)
		msgBox.setDefaultButton(QtWidgets.QMessageBox.No)
		ret = msgBox.exec_()

		if ret == QtWidgets.QMessageBox.No:
			return
		else:
			for patient in entry2remove:
				os.system("rm -rf "+installationDirectory+"/DB/"+patient)
				msgBox = QtWidgets.QMessageBox()
				msgBox.setText("Entries deleted!")
				ret = msgBox.exec_()
			
		self.listWidget.clear()
		patientList = []
		for folderName in os.listdir(installationDirectory+"/DB/"):
			if not folderName[0] == ".":
				patientList.append(folderName)
		patientList = sorted(patientList)
		for patient in patientList: 
			self.listWidget.addItem(patient)



	
	def exitWindow(self):
		exit(1)


	def retranslateUi(self, Form):
		_translate = QtCore.QCoreApplication.translate
		Form.setWindowTitle(_translate("Form", "Delete entry"))
		self.label.setText(_translate("Form", "Select entries to remove"))
		self.remove_pushButton.setText(_translate("Form", "Remove entry from DB"))
		self.pushButton.setText(_translate("Form", "Cancel"))

import sys
if __name__ =="__main__":
	App = QtWidgets.QApplication(sys.argv)
	Form = QtWidgets.QWidget()
	ui = Ui_Form()
	installationDirectory = sys.argv[1]
	ui.setupUi(Form, installationDirectory)
	Form.show()
	sys.exit(App.exec_())