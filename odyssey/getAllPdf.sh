find . -name '*.pdf' -print > file_list.txt
tar -czvf resultPdf.tar.gz -T file_list.txt
