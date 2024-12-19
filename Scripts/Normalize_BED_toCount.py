import argparse

parser = argparse.ArgumentParser()
parser.add_argument("Root", help= "Root of BED Files")
parser.add_argument("--Mode", help="Modes: F, B, R, L: F = Frequency, B = Both Cutting sites, R = Right, L = Left")
parser.add_argument("--MicroInDel_Length", help="Minimun Coverage at cutting site")
parser.add_argument("--MinCov", help="Minimun Coverage at cutting site")
parser.add_argument("--MinCount", help="Minimum number of reads of Recombination event")
args = parser.parse_args()

## Handle options
if args.Mode:
    Mode = str(args.Mode)
else:
    Mode = 'B'
if args.MinCov:
    MinCov = int(args.MinCov)
else:
    MinCov = 1
if args.MinCount:
    MinCount = int(args.MinCount)
else:
    MinCount = 1
if args.MicroInDel_Length:
    MicroInDel_Length = int(args.MicroInDel_Length)
else:
    MicroInDel_Length = 0

class RecEvent(object):
    def __init__(self, line):
        [self.Ref, self.Start, self.Stop, self.Type, self.Count, self.Dir, self.CSL, self.CSR, self.Seq1, self.Seq2] = line
        self.Name = self.Start +  "_to_" + self.Stop + "_#_" + self.Count
        self.Count = int(self.Count)
        self.Gap = max(int(self.Stop), int(self.Start)) - min(int(self.Stop), int(self.Start)) - 1
        self.CSL = int(self.CSL)
        self.CSR = int(self.CSR)
        self.CSB = (self.CSR + self.CSL)/2
        try:
            self.BCount = self.Count/self.CSB
            self.LCount = self.Count/self.CSL
            self.RCount = self.Count/self.CSR
        except:
            self.Type = 'NoCov'
    def __str__(self):
        return self.Name

Events = {}
Refs = set()
with open(str(args.Root + '_ViReMa_comb_Recombination_Results_noDir_WA1coords.bed'), "r") as In:
    Header = In.readline()
    Data = In.readline()
    while Data:
        Data = Data.split()
        Name = Data[1] +  "_to_" + Data[2] + "_#_" + Data[4]
        Refs.add(Data[0])
        Lib = Data[0] + ":" + Data[5]
        if Lib in Events:
            Events[Lib][Name] = RecEvent(Data)
        else:
            Events[Lib] = {}
            Events[Lib][Name] = RecEvent(Data)
        Data = In.readline()

Temp = []
for Lib in Events:
    for Name in Events[Lib]:
        if Events[Lib][Name].Type != 'NoCov':
            if Events[Lib][Name].Gap <= MicroInDel_Length:
                #if int(Events[Lib][Name].Start) + 1 < int(Events[Lib][Name].Stop):
                Events[Lib][Name].Type = 'u' + Events[Lib][Name].Type
            else:
                Events[Lib][Name].Type = Events[Lib][Name].Type
            #else:  ##DVG
                #if int(Events[Lib][Name].Start) < int(Events[Lib][Name].Stop):
                    #Events[Lib][Name].Type = Events[Lib][Name].Type
                #else:
                    #Events[Lib][Name].Type = 'Insertion'
            Temp.append([Events[Lib][Name].Ref, Events[Lib][Name].Start, Events[Lib][Name].Stop, 
                         Events[Lib][Name].Type, str(Events[Lib][Name].BCount), Events[Lib][Name].Dir])
        else:
            pass

Temp.sort(key=lambda a:float(a[4]), reverse=True)
with open(str(args.Root + "_normalised.bed"), "w") as Out:
    Out.write(Header)
    for i in Temp:
        Out.write('\t'.join(i) + '\n')
