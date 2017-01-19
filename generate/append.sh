
cat hits_train.csv | sed "s/\[//g"  | sed "s/\]//g" >> hits_train_merged.csv
