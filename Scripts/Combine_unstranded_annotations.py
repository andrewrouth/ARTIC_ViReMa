import argparse
parser = argparse.ArgumentParser()
parser.add_argument("InBED", help= "Input BED File")
parser.add_argument("OutBED", help= "Root of output files")
parser.add_argument("-BED12", action='store_true', help= "Select ViReMa output was in BED12")
parser.add_argument("-Stranded", action='store_true', help= "Select if coverage from ViReMa output was in BED12 files and as '-Stranded'")
args = parser.parse_args()

if args.Stranded:
    Stranded = True
else:
    Stranded = False
if args.BED12:
    BED12 = True
else:
    BED12 = False

class RecEvent(object):
    def __init__(self, line, Name):
        if BED12: ###Just realized this is actually BED10.....Need to revisit ViReMa
            [self.Ref, self.Start, self.Stop, 
             self.Type, self.Count, self.Dir, 
             self.CSL, self.CSR, self.Seq1, self.Seq2] = line
        else:
            [self.Ref, self.Start, self.Stop, 
             self.Type, self.Count, self.Dir] = line
        self.Name = Name
        self.Count = int(self.Count)
        self.CSL = int(self.CSL)
        self.CSR = int(self.CSR)
    def __str__(self):
        return self.Name

Dict = {}      
Out = open(str(args.OutBED),'w')
with open(str(args.InBED),'r') as In:
    Header = In.readline()
    Out.write(Header)
    Data = In.readline().rstrip().split()
    while Data:
        if Data[5] == '+':
            Name = Data[0] + ":" + Data[1] +  "_to_" + Data[2]
        else:
            Name = Data[0] + ":" + Data[2] +  "_to_" + Data[1]
        if Name in Dict:
            Dict[Name].Count += int(Data[4])
            if BED12:
                if Stranded:    #If -Stranded, then Coverage seperated reported for each stranded
                    Dict[Name].CSL += int(Data[6])
                    Dict[Name].CSR += int(Data[7])
                else:           #If not -Stranded, then Coverage will be same for both and already combined
                    pass
            else:
                pass
        else:
            Dict[Name] = RecEvent(Data, Name)
        Data = In.readline().rstrip().split()
        
for i in Dict:
    if BED12:
        line = '\t'.join([Dict[i].Ref, Dict[i].Start, Dict[i].Stop, 
                          Dict[i].Type, str(Dict[i].Count), '+/-', 
                          str(Dict[i].CSL), str(Dict[i].CSR), Dict[i].Seq1, Dict[i].Seq2])
    else:
        line = '\t'.join([Dict[i].Ref, Dict[i].Start, Dict[i].Stop, 
                          Dict[i].Type, str(Dict[i].Count), '+/-',])
    Out.write(line + '\n')
Out.close()