#!/usr/bin/env python
import os,sys
import numpy as np
import argparse as arg

def printnumblocks(List,NER,fmt):
  N     = len(List)          # Total number of elements in List
  NER   = int(NER)           # Number of element per row
  NFR   = N/NER              # Number of full row
  NLR   = N-NFR*NER          # Number of element in the last row
  if ( N >= NER):
    Fmt = ((fmt)*NER+"\n")*(N/NER)+((fmt)*NLR+"\n")
  else:
    Fmt = ((fmt)*N+"\n")
  print (Fmt % tuple(List))
  pass

def stringconverter(In):

  Out = []
  blocks = In.split(',')
  for block in blocks:
    blk = block.split('-')
    if len(blk) > 1:
      sel =  map(int,blk)
      Out += range(sel[0],sel[1]+1,1)
    else:
      Out += [int(block)]
  return Out

#
# Atom selection to perform alignment
#
def atomselect(Selection):
  try: Sel = stringconverter(Selection)
  except:
    print ("\n ERROR! \n Specify a valid selecion string (e.g. 1-15,12,13)\n")
    sys.exit()
  print ("\n Selection will be performed on the following %d atoms:" % len(Sel))
  printnumblocks(Sel,15,"%3d")
  Sel = (np.array(stringconverter(Selection))-1).tolist()
  return Sel

if __name__ == '__main__':

  parser = arg.ArgumentParser(description="Calculate the dCT index for a listed subset of atoms")
  parser.add_argument('-inp','-i',metavar='InpFile',help='Input File containing the coordinates and the charges',required=True)
  parser.add_argument('-selatom',help='Atom selection',default='all')
  parser.add_argument('-out','-o',help='Output file with coordinates + baricenters',default=None)
  args = parser.parse_args()

  InpFile    = args.inp
  InpSel   = args.selatom
  OutFile    = args.out

  if OutFile is None: OutFile = "PRJ_%s" % InpFile
  print ("\n The Input File is  : %s" % InpFile) 
  print (" The Output File is : %s" % OutFile)
# Take the information from the input file 
  AtomSel = atomselect(InpSel)
  coord=[]
  AtomicNumber=[]
  EsCharge=[]
  GsCharge=[]
  with open ("%s" % (InpFile), 'r') as f:
    for line in f:
      try:
        AN=line.split()[0] 
        x=float(line.split()[1])
        y=float(line.split()[2])
        z=float(line.split()[3])
        esC=float(line.split()[4])
        gsC=float(line.split()[5])
        coord.append([x,y,z])
        AtomicNumber.append(AN)
        EsCharge.append(esC)
        GsCharge.append(gsC)
      except:
        break
#Transform the vector in np.array
  coord=np.array(coord)    
  EsCharge=np.array(EsCharge)    
  GsCharge=np.array(GsCharge)
  DiffCharges=EsCharge-GsCharge
#Calculate the Different Charge Contribution
  gsCharge=0.
  esCharge=0.
  diffCharge=0.
  for i,item in enumerate(AtomSel):
    #print i,item,GsCharge[item-1]
    gsCharge=gsCharge+GsCharge[item]   
    esCharge=esCharge+EsCharge[item]   
    diffCharge=diffCharge+DiffCharges[item]
  print (" The GS total charge for the fragment is   : %12.6f" % gsCharge)
  print (" The Es total charge for the fragment is   : %12.6f" % esCharge)
  print (" The Diff total charge for the fragment is : %12.6f \n" % diffCharge)
#Calculate the Baricenters
  BarPos=np.zeros(3)
  BarNeg=np.zeros(3)
  posCharge=0.
  negCharge=0.
  for i,item in enumerate(AtomSel):
    if DiffCharges[item] > 0:
      posCharge=posCharge+DiffCharges[item]
      BarPos[:]=BarPos[:]+coord[item,:]*DiffCharges[item]
    else:
      negCharge=negCharge+DiffCharges[item]
      BarNeg[:]=BarNeg[:]+coord[item,:]*DiffCharges[item]
  BarPos=BarPos/posCharge
  BarNeg=BarNeg/negCharge
  dCT=np.linalg.norm(BarPos-BarNeg)
  print ("\n The dCT value is : %12.6f angstrom \n" % dCT) 
#Write the coordinates+baricenters
  out=open(OutFile, 'w')
  out.write("%d \n" % (len(DiffCharges)+2))
  out.write("\n")
  for i in range(len(DiffCharges)):
    out.write("%-3s  %8.4f  %8.4f  %8.4f  \n" % (AtomicNumber[i],coord[i][0],coord[i][1],coord[i][2]))
  out.write("2    %8.4f  %8.4f  %8.4f  \n" % (BarPos[0],BarPos[1],BarPos[2]))
  out.write("2    %8.4f  %8.4f  %8.4f  \n" % (BarNeg[0],BarNeg[1],BarNeg[2]))
