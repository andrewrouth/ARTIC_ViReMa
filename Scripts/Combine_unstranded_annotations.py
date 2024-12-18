import argparse
parser = argparse.ArgumentParser()
parser.add_argument("InBED", help= "Input BED File")
parser.add_argument("OutBED", help= "Root of output files")
parser.add_argument("-CountCombine",  action='store_true', help= "")
args = parser.parse_args()

if args.CountCombine:
    CountCombine = True
else:
    CountCombine = False

class RecEvent(object):
    def __init__(self, line, Name):
        [self.Ref, self.Start, self.Stop, 
         self.Type, self.Count, self.Dir, 
         self.CSL, self.CSR, self.Seq1, self.Seq2] = line
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
        #print(Data)
        if Data[0] == 'track':
            print("New Header in middle of file: ", Data)
        else:
            if Data[5] == '-':
                Data[2],Data[1] = Data[1],Data[2]
                #Name = Data[0] + ":" + Data[1] +  "_to_" + Data[2]
            else:
                #Name = Data[0] + ":" + Data[2] +  "_to_" + Data[1]
                pass
            Name = Data[0] + ":" + Data[1] +  "_to_" + Data[2]
            if Name in Dict:
                Dict[Name].Count += int(Data[4])
                if CountCombine:
                    Dict[Name].CSL += int(Data[6])
                    Dict[Name].CSR += int(Data[7])
                else:
                    pass
            else:
                Dict[Name] = RecEvent(Data, Name)
        Data = In.readline().rstrip().split()
        
for i in Dict:
    line = '\t'.join([Dict[i].Ref, Dict[i].Start, Dict[i].Stop, 
                      Dict[i].Type, str(Dict[i].Count), '+/-', 
                      str(Dict[i].CSL), str(Dict[i].CSR), Dict[i].Seq1, Dict[i].Seq2])
    Out.write(line + '\n')
Out.close()
