import struct
import numpy as np


def write_block(file, data):
    sizestruct = struct.pack("i", len(data))
    file.write(sizestruct)
    file.write(data)
    file.write(sizestruct)


# generate random data arrays
x = np.random.rand(100)
y = np.random.rand(100)
z = np.random.rand(100)
h = np.random.rand(100).astype("f")

file = open("Phantomtest.dat", "wb")

# write header
# most info here is ignored by PhantomSnapshotDensityFunction, except the units
write_block(file, struct.pack("7s", b"ignored"))
write_block(file, struct.pack("2s", b"FT"))

# write ints block
write_block(file, struct.pack("i", 2))
tags = [b"nblocks         ", b"ints block      "]
tagstr = b""
for tag in tags:
    tagstr += struct.pack("16s", tag)
write_block(file, tagstr)
vals = np.zeros(2, dtype="i")
vals[0] = 1
vals[1] = 42
write_block(file, struct.pack("i" * 2, *vals))

# write int8s block
write_block(file, struct.pack("i", 2))
tags = [b"ignored         ", b"int8s block     "]
tagstr = b""
for tag in tags:
    tagstr += struct.pack("16s", tag)
write_block(file, tagstr)
vals = np.zeros(2, dtype="b")
vals[0] = 1
vals[1] = 2
write_block(file, struct.pack("b" * 2, *vals))

# write int16s block
write_block(file, struct.pack("i", 2))
tags = [b"ignored         ", b"int16s block   "]
tagstr = b""
for tag in tags:
    tagstr += struct.pack("16s", tag)
write_block(file, tagstr)
vals = np.zeros(2, dtype="h")
vals[0] = 1
vals[1] = 2
write_block(file, struct.pack("h" * 2, *vals))

# write int32s block
write_block(file, struct.pack("i", 2))
tags = [b"ignored         ", b"int32s block   "]
tagstr = b""
for tag in tags:
    tagstr += struct.pack("16s", tag)
write_block(file, tagstr)
vals = np.zeros(2, dtype="i")
vals[0] = 1
vals[1] = 2
write_block(file, struct.pack("i" * 2, *vals))

# write int64s block
write_block(file, struct.pack("i", 2))
tags = [b"npartoftype     ", b"int64s block   "]
tagstr = b""
for tag in tags:
    tagstr += struct.pack("16s", tag)
write_block(file, tagstr)
vals = np.zeros(2, dtype="l")
vals[0] = 100
vals[1] = 2
write_block(file, struct.pack("l" * 2, *vals))

# write reals block
write_block(file, struct.pack("i", 2))
tags = [b"massoftype      ", b"tag not used    "]
tagstr = b""
for tag in tags:
    tagstr += struct.pack("16s", tag)
write_block(file, tagstr)
vals = np.zeros(2, dtype="d")
vals[0] = 0.01
vals[1] = 42.0
write_block(file, struct.pack("d" * 2, *vals))

# write real4s block
write_block(file, struct.pack("i", 2))
tags = [b"massoftype      ", b"tag not used    "]
tagstr = b""
for tag in tags:
    tagstr += struct.pack("16s", tag)
write_block(file, tagstr)
vals = np.zeros(2, dtype="f")
vals[0] = 0.01
vals[1] = 42.0
write_block(file, struct.pack("f" * 2, *vals))

# write real8s block
write_block(file, struct.pack("i", 2))
tags = [b"umass           ", b"udist           "]
tagstr = b""
for tag in tags:
    tagstr += struct.pack("16s", tag)
write_block(file, tagstr)
vals = np.zeros(2, dtype="d")
vals[0] = 1.0
vals[1] = 1.0
write_block(file, struct.pack("d" * 2, *vals))

write_block(file, struct.pack("i", 2))

# write data block
numbers = [0, 0, 0, 0, 0, 3, 1, 0]
write_block(file, struct.pack("liiiiiiii", 100, *numbers))
numbers = [0, 0, 0, 0, 0, 0, 0, 0]
write_block(file, struct.pack("liiiiiiii", 0, *numbers))

write_block(file, struct.pack("16s", b"x               "))
write_block(file, struct.pack("d" * 100, *x))
write_block(file, struct.pack("16s", b"y               "))
write_block(file, struct.pack("d" * 100, *y))
write_block(file, struct.pack("16s", b"z               "))
write_block(file, struct.pack("d" * 100, *z))
write_block(file, struct.pack("16s", b"h               "))
write_block(file, struct.pack("f" * 100, *h))

# we don't write data after this point, since we don't need to read it...

# write reference file
reffile = open("Phantom_data.txt", "w")
for i in range(100):
    reffile.write(
        "{x:.14e}\t{y:.14e}\t{z:.14e}\t{h:.14e}\n".format(
            x=x[i], y=y[i], z=z[i], h=h[i]
        )
    )
