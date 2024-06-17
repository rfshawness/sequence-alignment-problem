# changed argument specifications
for file in "../datapoints/$1"*.txt
do
   python3 basic_3.py "$file" "$2".txt
done