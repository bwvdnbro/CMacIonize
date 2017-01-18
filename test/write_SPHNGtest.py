import struct
import numpy as np

def write_block(file, data):
  sizestruct = struct.pack('i', len(data))
  file.write(sizestruct)
  file.write(data)
  file.write(sizestruct)

# generate random data arrays
x = np.random.rand(100)
y = np.random.rand(100)
z = np.random.rand(100)
m = np.random.rand(100)
h = np.random.rand(100)
vx = np.random.rand(100)
vy = np.random.rand(100)
vz = np.random.rand(100)
u = np.random.rand(100)
rho = np.random.rand(100)

file = open("SPHNGtest.dat", "wb")

# write header
# most info here is ignored by SPHNGSnapshotDensityFunction, except the units
write_block(file, struct.pack('c'*7, *"ignored"))
write_block(file, struct.pack('c'*2, *"FT"))
write_block(file, struct.pack('i', 44))
tags = ['nparttot        ','n1              ','n2              ',
        'nreassigntot    ','naccretetot     ','nkilltot        ',
        'nblocks         ','iyr             ','idum            ',
        'iv              ','iv              ','iv              ',
        'iv              ','iv              ','iv              ',
        'iv              ','iv              ','iv              ',
        'iv              ','iv              ','iv              ',
        'iv              ','iv              ','iv              ',
        'iv              ','iv              ','iv              ',
        'iv              ','iv              ','iv              ',
        'iv              ','iv              ','iv              ',
        'iv              ','iv              ','iv              ',
        'iv              ','iv              ','iv              ',
        'iv              ','iv              ','iplanetesimals  ',
        'irotpot         ','idragscheme     ']
tagstr = ""
for tag in tags:
  tagstr += struct.pack('c'*16, *tag)
write_block(file, tagstr)
npart = np.zeros(44, dtype = 'i')
npart[0] = 100
npart[6] = 1
write_block(file, struct.pack('i'*44, *npart))
write_block(file, struct.pack('c'*7, *"ignored"))
write_block(file, struct.pack('c'*7, *"ignored"))
write_block(file, struct.pack('c'*7, *"ignored"))
write_block(file, struct.pack('i', 1))
write_block(file, struct.pack('c'*4, *"tags"))
write_block(file, struct.pack('L', 100))
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
write_block(file, struct.pack('i', 1))
write_block(file, struct.pack('i', 4))
tags = ['udist           ','umass           ','utime           ',
        'umagfd          ']
tagstr = ""
for tag in tags:
  tagstr += struct.pack('c'*16, *tag)
write_block(file, tagstr)
units = [1., 1., 1., 1.] # cgs units
write_block(file, struct.pack('d'*4, *units))
write_block(file, struct.pack('i', 2))

# write data block
numbers = [2, 1, 0, 0, 1, 9, 4, 0]
write_block(file, struct.pack('LIIIIIIII', 100, *numbers))
numbers = [1, 0, 0, 0, 0, 9, 0, 0]
write_block(file, struct.pack('LIIIIIIII', 1, *numbers))
write_block(file, struct.pack('c'*4, *"tags"))
numsteps = np.zeros(100, dtype = 'I')
write_block(file, struct.pack('I'*100, *numsteps))
write_block(file, struct.pack('c'*4, *"tags"))
write_block(file, struct.pack('c'*7, *"ignored"))
write_block(file, struct.pack('c'*16, *"iphase          "))
# make one random sink
iphase = np.zeros(100, dtype = 'c')
iphase[42] = -1
write_block(file, iphase)
write_block(file, struct.pack('c'*16, *"iunique         "))
iunique = np.zeros(100, dtype = 'i')
write_block(file, struct.pack('i'*100, *iunique))
write_block(file, struct.pack('c'*16, *"x               "))
write_block(file, struct.pack('d'*100, *x))
write_block(file, struct.pack('c'*16, *"y               "))
write_block(file, struct.pack('d'*100, *y))
write_block(file, struct.pack('c'*16, *"z               "))
write_block(file, struct.pack('d'*100, *z))
write_block(file, struct.pack('c'*16, *"m               "))
write_block(file, struct.pack('d'*100, *m))
write_block(file, struct.pack('c'*16, *"h               "))
write_block(file, struct.pack('d'*100, *h))
write_block(file, struct.pack('c'*16, *"vx              "))
write_block(file, struct.pack('d'*100, *vx))
write_block(file, struct.pack('c'*16, *"vy              "))
write_block(file, struct.pack('d'*100, *vy))
write_block(file, struct.pack('c'*16, *"vz              "))
write_block(file, struct.pack('d'*100, *vz))
write_block(file, struct.pack('c'*16, *"u               "))
write_block(file, struct.pack('d'*100, *u))
write_block(file, struct.pack('c'*16, *"rho             "))
write_block(file, struct.pack('d'*100, *rho))
# igrad related blocks (number depends on the value of the last entry in
# numbers)
write_block(file, struct.pack('c'*7, *"ignored"))
write_block(file, struct.pack('c'*7, *"ignored"))
write_block(file, struct.pack('c'*7, *"ignored"))
write_block(file, struct.pack('c'*7, *"ignored"))
write_block(file, struct.pack('c'*7, *"ignored"))
write_block(file, struct.pack('c'*7, *"ignored"))
# sink particle data
for i in range(10):
  write_block(file, struct.pack('c'*4, *"tags"))
  write_block(file, struct.pack('c'*7, *"ignored"))

# we don't write data after this point, since we don't need to read it...

# untagged file
file = open("SPHNGtest_notags.dat", "wb")

# write header
# most info here is ignored by SPHNGSnapshotDensityFunction, except the units
write_block(file, struct.pack('c'*7, *"ignored"))
write_block(file, struct.pack('c'*2, *"FX"))
write_block(file, struct.pack('i', 44))
npart = np.zeros(44, dtype = 'i')
npart[0] = 100
npart[6] = 1
write_block(file, struct.pack('i'*44, *npart))
write_block(file, struct.pack('c'*7, *"ignored"))
write_block(file, struct.pack('c'*7, *"ignored"))
write_block(file, struct.pack('c'*7, *"ignored"))
write_block(file, struct.pack('i', 1))
write_block(file, struct.pack('L', 100))
write_block(file, struct.pack('i', 30))
header = np.zeros(30, dtype = 'd')
write_block(file, struct.pack('d'*30, *header))
write_block(file, struct.pack('i', 1))
write_block(file, struct.pack('i', 4))
units = [1., 1., 1., 1.] # cgs units
write_block(file, struct.pack('d'*4, *units))
write_block(file, struct.pack('i', 2))

# write data block
numbers = [2, 1, 0, 0, 1, 9, 4, 0]
write_block(file, struct.pack('LIIIIIIII', 100, *numbers))
numbers = [1, 0, 0, 0, 0, 9, 0, 0]
write_block(file, struct.pack('LIIIIIIII', 1, *numbers))
numsteps = np.zeros(100, dtype = 'I')
write_block(file, struct.pack('I'*100, *numsteps))
write_block(file, struct.pack('c'*7, *"ignored"))
# make one random sink
iphase = np.zeros(100, dtype = 'c')
iphase[42] = -1
write_block(file, iphase)
iunique = np.zeros(100, dtype = 'i')
write_block(file, struct.pack('i'*100, *iunique))
write_block(file, struct.pack('d'*100, *x))
write_block(file, struct.pack('d'*100, *y))
write_block(file, struct.pack('d'*100, *z))
write_block(file, struct.pack('d'*100, *m))
write_block(file, struct.pack('d'*100, *h))
write_block(file, struct.pack('d'*100, *vx))
write_block(file, struct.pack('d'*100, *vy))
write_block(file, struct.pack('d'*100, *vz))
write_block(file, struct.pack('d'*100, *u))
write_block(file, struct.pack('d'*100, *rho))
write_block(file, struct.pack('c'*7, *"ignored"))
write_block(file, struct.pack('c'*7, *"ignored"))
write_block(file, struct.pack('c'*7, *"ignored"))
write_block(file, struct.pack('c'*7, *"ignored"))
# sink particle data
for i in range(10):
  write_block(file, struct.pack('c'*7, *"ignored"))

# we don't write data after this point, since we don't need to read it...

# write reference file
reffile = open("SPHNG_data.txt", 'w')
for i in range(100):
  if not i == 42:
    reffile.write("{x:.14e}\t{y:.14e}\t{z:.14e}\t{m:.14e}\t{h:.14e}\n".format(
      x = x[i], y = y[i], z = z[i], m = m[i], h = h[i]))
