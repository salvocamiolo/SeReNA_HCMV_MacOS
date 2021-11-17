#!/bin/bash

cd -- "$(dirname -- "$BASH_SOURCE")"

installationPath=$(pwd)
echo "#!"$installationPath"/src/conda/bin/python" > SeReNA
echo "installationDirectory = \""$installationPath"/\"" >> SeReNA
echo " " >> SeReNA
cat SeReNA ./src/.HCMV_DB.py >temp ; mv temp SeReNA
chmod +x SeReNA

if  test -f "./src/conda/bin/conda"; then
	echo "Miniconda3 already installed"
else

	echo "Installing miniconda 3"
	cd src
	if  test -f "installation.log"; then
	rm installation.log
	fi
	touch installation.log
	bash Miniconda3-latest-MacOSX-x86_64.sh -b -p ./conda >> installation.log
	cd ../
	if  test -f "./src/conda/bin/conda"; then
	echo "Miniconda3 successfully installed"
	./src/conda/bin/conda config --set notify_outdated_conda false
	fi
fi


echo "Installing pyqt. Please wait...."
./src/conda/bin/conda install -c anaconda -y pyqt=5.9.2 >> installation.log
echo "Checking pyqt installation...."
./src/conda/bin/conda list > ./src/condaList
if grep -Fq pyqt ./src/condaList; then
	echo "pyqt was successfully installed"
else
	echo "pyqt was not installed. Please check file installation.log for details"
fi

echo "Installing exonerate. Please wait...."
./src/conda/bin/conda install -c bioconda -y exonerate=2.4.0 >> installation.log
echo "Checking exonerate installation...."
./src/conda/bin/conda list > ./src/condaList
if grep -Fq exonerate ./src/condaList; then
	echo "exonerate was successfully installed"
else
	echo "exonerate was not installed. Please check file installation.log for details"
fi

echo "Installing biopython. Please wait...."
	./src/conda/bin/conda install -c bioconda -y  biopython>> installation.log
	echo "Checking biopython installation...."
	./src/conda/bin/conda list > ./src/condaList
	if grep -Fq biopython ./src/condaList; then
		echo "biopython was successfully installed"
	else
		echo "biopython was not installed. Please check file installation.log for details"
	fi


	echo "Installing blast. Please wait...."
	./src/conda/bin/conda install -c bioconda -y  blast=2.9.0 >> installation.log
	echo "Checking blast installation...."
	./src/conda/bin/conda list > ./src/condaList
	if grep -Fq blast ./src/condaList; then
		echo "blast was successfully installed"
	else
		echo "blast was not installed. Please check file installation.log for details"
	fi


	echo "Installing samtools. Please wait...."
	./src/conda/bin/conda install -c bioconda -y samtools >> installation.log
	echo "Checking samtools installation...."
	./src/conda/bin/conda list > ./src/condaList
	if grep -Fq samtools ./src/condaList; then
		echo "samtools was successfully installed"
	else
		echo "samtools was not installed. Please check file installation.log for details"
	fi
	cp ./data/samtools ./src/conda/bin/

	echo "Installing spades. Please wait...."
	./src/conda/bin/conda install -c bioconda -y spades >> installation.log
	echo "Checking spades installation...."
	./src/conda/bin/conda list > ./src/condaList
	if grep -Fq spades ./src/condaList; then
		echo "spades was successfully installed"
	else
		echo "spades was not installed. Please check file installation.log for details"
	fi

	echo "Installing bowtie2. Please wait...."
	./src/conda/bin/conda install -c bioconda/label/broken -y  bowtie2 >> installation.log
	echo "Checking bowtie2 installation...."
	./src/conda/bin/conda list > ./src/condaList
	if grep -Fq -w bowtie2 ./src/condaList; then
		echo "bowtie2 was successfully installed"
	else
		echo "bowtie2 was not installed. Please check file installation.log for details"
	fi


	echo "Installing exonerate. Please wait...."
	./src/conda/bin/conda install -c bioconda -y exonerate >> installation.log
	echo "Checking exonerate installation...."
	./src/conda/bin/conda list > ./src/condaList
	if grep -Fq exonerate ./src/condaList; then
		echo "exonerate was successfully installed"
	else
		echo "exonerate was not installed. Please check file installation.log for details"
	fi


	echo "Installing clustal omega. Please wait...."
	./src/conda/bin/conda install -c bioconda -y clustalo >> installation.log
	echo "Checking clustalo installation...."
	./src/conda/bin/conda list > ./src/condaList
	if grep -Fq clustalo ./src/condaList; then
		echo "clustalo was successfully installed"
	else
		echo "clustalo was not installed. Please check file installation.log for details"
	fi


	echo "Installing fasttree. Please wait...."
	./src/conda/bin/conda install -c bioconda -y fasttree >> installation.log
	echo "Checking fasttree installation...."
	./src/conda/bin/conda list > ./src/condaList
	if grep -Fq fasttree ./src/condaList; then
		echo "fasttree was successfully installed"
	else
		echo "fasttree was not installed. Please check file installation.log for details"
	fi



	echo "Installing matplotlib. Please wait...."
	./src/conda/bin/conda install -c conda-forge matplotlib -y  >> installation.log
	echo "Checking matplotlib installation...."
	./src/conda/bin/conda list > ./src/condaList
	if grep -Fq matplotlib ./src/condaList; then
		echo "matplotlib was successfully installed"
	else
		echo "matplotlib was not installed. Please check file installation.log for details"
	fi

	echo "Installing entrez-direct. Please wait...."
	./src/conda/bin/conda install -c bioconda entrez-direct -y  >> installation.log
	echo "Checking entrez-direct installation...."
	./src/conda/bin/conda list > ./src/condaList
	if grep -Fq "entrez-direct" ./src/condaList; then
		echo "entrez-direct was successfully installed"
	else
		echo "entrez-direct was not installed. Please check file installation.log for details"
	fi

	./src/conda/bin/conda install -c bioconda java-jdk -y  >> installation.log
	

	#echo "Installing pillow. Please wait...."
	#./src/conda/bin/conda install -c anaconda -y pillow >> installation.log
	#echo "Checking pillow installation...."
	#./src/conda/bin/conda list > condaList
	#if grep -Fq pillow condaList; then
	#	echo "Pillow was successfully installed"
	#else
	#	echo "Pillow was not installed. Please check file installation.log for details"
	#fi

	