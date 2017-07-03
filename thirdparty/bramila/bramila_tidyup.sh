#!/bin/bash
toolboxfolder='/m/nbe/scratch/braindata/shared/toolboxes/bramila/bramila/';
if [ "x$1" == "x" ]; then
	echo "Input parameter missing. The usage of the script is:"
	echo "./bramila_tidyup.sh <name_of_folder_with_preprocessing_results>";
	exit -1;
fi
if [ "x$2" == "x" ]; then
	echo "You have not specified a second parameter. I assume you want to do a dry-run, no files will be deleted. To actually delete files use:"
	echo "./bramila_tidyup.sh <name_of_folder_with_preprocessing_results> delete";
fi
actualrun=0;
prefix="This is a dry run. The script would have deleted "
if [ "x$2" == "xdelete" ]; then
	actualrun=1;
	prefix="Deleting ";
fi




for n in $(cat $toolboxfolder/bramila_files.txt|grep -v '#');do
		
	file=$(echo $1/$n); 
	if [ ! -w $file ]; then
		echo "Requested resource $file does not exist or is not writable by this user";
		continue;
	else
		echo "$prefix $file"
		if [ $actualrun -eq 1 ]; then
			if [ ${n: -1} == "/" ]; then
		 
				echo "rm -rfv $file"
				rm -rfv $file
			else
				echo "rm -fv $file"
				rm -fv $file
			fi
		fi
	fi
done

