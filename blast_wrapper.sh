echo "start"
if [[ $1 == "" ]]
then
	echo "########################################"
	echo "arguments"
	echo "\$1 folder with fastas - query"
	echo "\$2 path to database or folder with dbs - target"
	echo "\$3 threshold evalue"
	echo "\$4 method: tblastx, blastn"
	echo "\$5 cpu: 4"
	echo "\$6 cluster: y / n"
	echo "\$7 list (for multiple dbs)"
	echo "########################################"
else 
	if [[ $6 == y ]]
	then
		echo "loading modules"
		module load ncbi-blast
	else 
		echo "no modules loaded"
	fi
	echo "erasing the content of the blast output file..."
	rm *.blast
	if [ -z "$7" ]
	  then
	    echo "single db blast option selected"
	    echo "blasting queries:"
		declare -i CT2=0
		declare -i TOTAL=$(ls $1/*.fas | cat | wc -l) #number of ahe files
		for f in $(ls $1/*.fas)
		do
			COUNT=0
			zo=$(( CT2*100 / TOTAL ))
			echo -ne "                                                    \r"
			echo -ne $zo"%\treading... \r"
			while read LINE
			do
				if [[ $LINE =~ ">" ]] && [[ $COUNT == 0 ]]
				then
					echo \>$f > fasextr.blast
					COUNT=1
					continue
				elif [[ $COUNT == 1 ]]
				then
					if [[ $LINE =~ ">" ]]
					then
						break
		 			else
						echo $LINE >> fasextr.blast
			 		fi
			 	fi
			done < $f
			name=$(echo $f | rev | cut -d'/' -f-1 | rev)
		 	CT2=$CT2+1
		 	echo -ne "                                                                          \r"
		 	echo -ne $zo"% blast $name against $2...\r"
		 	touch output.blast
		 	$4 -db $2 -query fasextr.blast -out output.blast -outfmt 6 -num_threads $5 -evalue $3
		 	blname=$(echo $2 | rev | cut -d'/' -f-1 | rev )
		 	#cat output.blast | sort -k1,1 -k2,2 -k11,11g -k12,12nr | sort -u -k2,2 --merge | sort -k1,1 -k11,11g -k12,12nr >> $blname.blast
		 	cat output.blast | sort -k1,1 -k2,2 -k11,11g -k12,12nr >> $blname.blast
		done
	  else
	  	echo "multiple db blast option selected, number of dbs: $(cat $7 | wc -l)"
	  	echo "blasting queries:"
		declare -i CT2=0
		declare -i TOTAL=$(ls $1/*.fas | cat | wc -l)*$(cat $7 | wc -l) #number of files
		zo=$(( CT2*100 / TOTAL ))
		for f in $(ls $1/*.fas)
		do
			COUNT=0
			echo -ne "                                                    \r"
			echo -ne $zo"%\treading... \r"
			while read LINE
			do
				if [[ $LINE =~ ">" ]] && [[ $COUNT == 0 ]]
				then
					echo \>$f > fasextr.blast
					COUNT=1
					continue
				elif [[ $COUNT == 1 ]]
				then
					if [[ $LINE =~ ">" ]]
					then
						break
		 			else
						echo $LINE >> fasextr.blast
			 		fi
			 	fi
			done < $f
			name=$(echo $f | rev | cut -d'/' -f-1 | rev)
			while read fdb
			do
				CT2=$CT2+1
				echo -ne "                                                                          \r"
		 		echo -ne $zo"% blast $name against $fdb...\r"
		 		zo=$(( CT2*100 / TOTAL ))
		 		touch output.blast
		 		$4 -db $2/$fdb -query fasextr.blast -out output.blast -outfmt 6 -num_threads $5 -evalue $3
		 		#cat output.blast | sort -k1,1 -k2,2 -k11,11g -k12,12nr | sort -u -k2,2 --merge | sort -k1,1 -k11,11g -k12,12nr >> $fdb.blast
		 		cat output.blast | sort -k1,1 -k2,2 -k11,11g -k12,12nr >> $fdb.blast
		 	done < $7
		done
	fi
	echo ""
	rm output.blast fasextr.blast
fi
echo "done"
