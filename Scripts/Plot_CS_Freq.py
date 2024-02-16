import argparse

parser = argparse.ArgumentParser()
parser.add_argument("Root", help= "Root of BED Files")
parser.add_argument("BED", help= "BED File")
parser.add_argument("FASTA", help= "Reference FASTA File")
parser.add_argument("--Mode", help="Modes: F, B, R, L: F = Frequency, B = Both Cutting sites, R = Right, L = Left")
parser.add_argument("-CoVData", action='store_true', help="Select if partitioning out sgmRNAs in CoV data")
parser.add_argument("--sgmRNA_Threshold", help="Set threshold nt for sgmRNA search. Default = 80")
parser.add_argument("--sgmRNA_Fuzz", help="Set threshold nt for sgmRNA search +/- from canoncial acceptor site. Default = 5")
parser.add_argument("-Ends", action='store_true', help="Select if partitioning out end fusions (e.g. in CoV data)")
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
if args.CoVData:
    CoVData = True
else:
    CoVData = False
if args.sgmRNA_Threshold:
    sgmRNA_Threshold = int(args.sgmRNA_Threshold)
else:
    sgmRNA_Threshold = 80
if args.sgmRNA_Fuzz:
    sgmRNA_Fuzz = int(args.sgmRNA_Fuzz)
else:
    sgmRNA_Fuzz = 5
if args.Ends:
    Ends = True
else:
    Ends = False
if args.MicroInDel_Length:
    MicroInDel_Length = int(args.MicroInDel_Length)
else:
    MicroInDel_Length = 5

Canonical_sgmRNAs = {'21557':  'Spike',
                    '25386':  'ORF3a',
                    '26238':  'E',
                    '26474':  'M',
                    '27042':  'ORF6',
                    '27389':  'ORF7a',
                    '27761':  'ORF7b',
                    '27889':  'ORF8',
                    '28261':  'N',
                    '28285':  'N*'}

RefSeqs = {}
with open(str(args.FASTA),'r') as FastaIn:
    data = FastaIn.readline().rstrip()
    while data:
        if data[0] == '>':
            Name = data[1:].split()[0]
            RefSeqs[Name] = ''   
            Seq = ''
        else:
            RefSeqs[Name] += data
        data = FastaIn.readline().rstrip()

def Rev_Comp(Seq):
        Seq = Seq.upper()
        basecomplement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
        letters = list(Seq)
        letters = [basecomplement[base] for base in letters]
        return ''.join(letters)[::-1]

class RecEvent(object):
    def __init__(self, line):
        [self.Ref, self.Start, self.Stop, self.Type, self.Count, self.Dir, self.CSL, self.CSR, self.Seq1, self.Seq2] = line
        self.Name = self.Start +  "_to_" + self.Stop + "_#_" + self.Count
        if self.Dir == '-':
            self.Start, self.Stop = self.Stop, self.Start
        else:
            pass
        self.Count = int(self.Count)
        self.Gap = max(int(self.Stop), int(self.Start)) - min(int(self.Stop), int(self.Start)) - 1
        self.CSL = int(self.CSL)
        self.CSR = int(self.CSR)
        self.CSB = (self.CSR + self.CSL)/2.
        try:
            self.BCount = self.Count/self.CSB
            self.LCount = self.Count/self.CSL
            self.RCount = self.Count/self.CSR
        except:
            self.Type = 'NoCov'
        if int(self.Start) < sgmRNA_Threshold:
            #print(self.Name)
            ##sgmRNA:
            #if self.Start == '70' and self.Stop in Canonical_sgmRNAs:
            sgmRNARange = range(int(self.Stop) - sgmRNA_Fuzz, int(self.Stop) + sgmRNA_Fuzz)
            #if any([str(i) in Canonical_sgmRNAs for i in sgmRNARange]):
            self.sgmRNA = 'NonCanonical'
            for i in sgmRNARange:
                if str(i) in Canonical_sgmRNAs:
                    self.sgmRNA = Canonical_sgmRNAs[str(i)]
                    break ##forces to find first matched (i.e. 5'-most) sgmRNA
                else: 
                    pass
            print(self.sgmRNA, self.Name)
        else:
            self.sgmRNA = False
    def __str__(self):
        return self.Name

Events = {}
Refs = set()
with open(str(args.BED), "r") as In:
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

with open(str(args.Root + '.coverage-stats.txt'), "r") as In:
    CovHeader = In.readline()
    Data = In.readline()
    [ref, startpos, endpos, numreads, covbases, coverage, meandepth, meanbaseq, meanmapq] = Data.split()
        
sgmRNA_Count = 0
EndFus_Count = 0
Deletion_Count = 0
Insertion_Count = 0
uIns_Count = 0
uDel_Count = 0

sgmRNA_Counts = {'Spike': 0,
                'ORF3a': 0,
                'E': 0,
                'M': 0,
                'ORF6': 0,
                'ORF7a': 0,
                'ORF7b': 0,
                'ORF8': 0,
                'N': 0,
                'N*': 0,
                'NonCanonical': 0}

for Lib in Events:
    for Name in Events[Lib]:
        if Events[Lib][Name].Type == 'NoCov':
            pass
        elif int(Events[Lib][Name].Start) < sgmRNA_Threshold and CoVData:
            Events[Lib][Name].Type = 'sgmRNA'
            sgmRNA_Count += Events[Lib][Name].Count
            sgmRNA_Counts[Events[Lib][Name].sgmRNA] += Events[Lib][Name].Count
        elif int(Events[Lib][Name].Stop) < 25 and Ends:
            Events[Lib][Name].Type = 'EndFus'
            EndFus_Count += Events[Lib][Name].Count
        elif Events[Lib][Name].Gap <= MicroInDel_Length:
            # if Events[Lib][Name].Dir == "+":
                if int(Events[Lib][Name].Start) + 1 < int(Events[Lib][Name].Stop):
                    Events[Lib][Name].Type = 'uDel'
                    uDel_Count += Events[Lib][Name].Count
                else:
                    Events[Lib][Name].Type = 'uIns'
                    uIns_Count += Events[Lib][Name].Count
            # elif Events[Lib][Name].Dir == "-":
            #     if int(Events[Lib][Name].Start) < int(Events[Lib][Name].Stop):
            #         Events[Lib][Name].Type = 'uIns'
            #         uIns_Count += Events[Lib][Name].Count
            #     else:
            #         Events[Lib][Name].Type = 'uDel'
            #         uDel_Count += Events[Lib][Name].Count
        else:  ##DVG
            # if Events[Lib][Name].Dir == "+":
                if int(Events[Lib][Name].Start) < int(Events[Lib][Name].Stop):
                    Events[Lib][Name].Type = 'Deletion'
                    Deletion_Count += Events[Lib][Name].Count
                else:
                    Events[Lib][Name].Type = 'Insertion'
                    Insertion_Count += Events[Lib][Name].Count
            # elif Events[Lib][Name].Dir == "-":
            #     if int(Events[Lib][Name].Start) < int(Events[Lib][Name].Stop):
            #         Events[Lib][Name].Type = 'Insertion'
            #         Insertion_Count += Events[Lib][Name].Count
            #     else:
            #         Events[Lib][Name].Type = 'Deletion'
            #         Deletion_Count += Events[Lib][Name].Count

sgmRNA_Jfreq = (sgmRNA_Count/int(numreads)) * 10000
EndFus_Jfreq = (EndFus_Count/int(numreads)) * 10000
Deletion_Jfreq = (Deletion_Count/int(numreads)) * 10000
Insertion_Jfreq = (Insertion_Count/int(numreads)) * 10000
uIns_Jfreq = (uIns_Count/int(numreads)) * 10000
uDel_Jfreq = (uDel_Count/int(numreads)) * 10000

sgmRNA_JFreqs = {}
for i in sgmRNA_Counts:
    sgmRNA_JFreqs[i] = (sgmRNA_Counts[i]/int(numreads)) * 10000

with open(str(args.Root + "_report.txt"), "w") as Out:
    #Out.write(Header)
    Out.write(''.join(['sgmRNA_Jfreq:\t', str(sgmRNA_Jfreq), 
                        '\nEndFus_Jfreq:\t', str(EndFus_Jfreq), 
                        '\nDeletion_Jfreq:\t', str(Deletion_Jfreq), 
                        '\nInsertio_Jfreq:\t', str(Insertion_Jfreq), 
                        '\nuIns_Jfreq:\t', str(uIns_Jfreq), 
                        '\nuDel_Jfreq:\t', str(uDel_Jfreq)]))
    for i in sgmRNA_JFreqs:
        Out.write('\n' + i + ':\t' + str(sgmRNA_JFreqs[i]))

print('sgmRNA_Jfreq:\t', sgmRNA_Jfreq, 
        '\nEndFus_Jfreq:\t', EndFus_Jfreq, 
        '\nDeletion_Jfreq:\t', Deletion_Jfreq, 
        '\nInsertio_Jfreq:\t', Insertion_Jfreq, 
        '\nuIns_Jfreq:\t', uIns_Jfreq, 
        '\nuDel_Jfreq:\t', uDel_Jfreq)
for i in sgmRNA_JFreqs:
    print(i + ':\t' +  str(sgmRNA_JFreqs[i]))

##Make Normalized Files
Temp = []
for Lib in Events:
    for Name in Events[Lib]:
        if Events[Lib][Name].Count >= MinCount:
            if Events[Lib][Name].Type in ['sgmRNA','EndFus'] and Events[Lib][Name].CSL >= MinCov:
                Temp.append([Events[Lib][Name].Ref, Events[Lib][Name].Start, Events[Lib][Name].Stop, 
                             Events[Lib][Name].Type, str(Events[Lib][Name].LCount), Events[Lib][Name].Dir])
            elif Events[Lib][Name].Type in ['Deletion', 'Insertion', 'uDel', 'uIns'] and Events[Lib][Name].CSB >= MinCov:
                Temp.append([Events[Lib][Name].Ref, Events[Lib][Name].Start, Events[Lib][Name].Stop, 
                             Events[Lib][Name].Type, str(Events[Lib][Name].BCount), Events[Lib][Name].Dir])
            else:
                print("Excluded NoCov: ", Name)


Temp.sort(key=lambda a:float(a[4]), reverse=True)
with open(str(args.Root + "_normalised.bed"), "w") as Out:
    Out.write(Header)
    for i in Temp:
        Out.write('\t'.join(i) + '\n')

