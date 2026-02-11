
for a in *.raw 
do
 ThermoRawFileParser.sh -i=$a  -o=./ -f=2
done 
