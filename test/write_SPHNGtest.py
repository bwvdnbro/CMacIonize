import struct
import numpy as np

def write_block(file, data):
  sizestruct = struct.pack('i', len(data))
  file.write(sizestruct)
  file.write(data)
  file.write(sizestruct)

file = open("SPHNGtest.dat", "wb")

# write header
# most info here is ignored by SPHNGSnapshotDensityFunction, except the units
write_block(file, struct.pack('c'*7, *"ignored"))
write_block(file, struct.pack('c'*2, *"FT"))
write_block(file, struct.pack('i', 44))
write_block(file, struct.pack('c'*4, *"tags"))
npart = np.zeros(44, dtype = 'i')
npart[0] = 100
npart[6] = 1
write_block(file, struct.pack('i'*44, *npart))
write_block(file, struct.pack('c'*7, *"ignored"))
write_block(file, struct.pack('c'*7, *"ignored"))
write_block(file, struct.pack('c'*7, *"ignored"))
write_block(file, struct.pack('i', 1))
write_block(file, struct.pack('c'*4, *"tags"))
write_block(file, struct.pack('i', 100))
write_block(file, struct.pack('i', 30))
tags = ['gt              ','dtmax           ','gamma           ',
        'rhozero         ','RK2             ','escap           ',
        'tkin            ','tgrav           ','tterm           ',
        'anglostx        ','anglosty        ','anglostz        ',
        'specang         ','ptmassin        ','tmag            ',
        'Bextx           ','Bexty           ','Bextz           ',
        'hzero           ','uzero_n2        ','hmass           ',
        'gapfac          ','                ','sdprof          ',
        'rorbit_orig     ','min_rplan       ','max_rplan       ',
        'planetesimalmass','coremass_orig   ','coremass        ']
tagstr = ""
for tag in tags:
  tagstr += struct.pack('c'*16, *tag)
write_block(file, tagstr)
header = np.zeros(30, dtype = 'd')
write_block(file, struct.pack('d'*30, *header))
write_block(file, struct.pack('c'*7, *"ignored"))
write_block(file, struct.pack('i', 4))
units = [1., 1., 1., 1.] # cgs units
write_block(file, struct.pack('d'*4, *units))
write_block(file, struct.pack('i', 2))

# write data block
