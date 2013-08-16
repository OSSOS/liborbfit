import os
import math
import ephem

__author__ = 'jjk'

class orbfit(object):

    def __init__(self):
        self.orbfit = ctypes.CDLL('orbfit/liborbfit.so')

import ctypes
import tempfile
mpc_lines="""     HL7j2    C2013 04 03.62926 17 12 01.16 +04 13 33.3          24.1 R      568
     HL7j2    C2013 04 04.58296 17 11 59.80 +04 14 05.5          24.0 R      568
     HL7j2    C2013 05 03.52252 17 10 38.28 +04 28 00.9          23.4 R      568
     HL7j2    C2013 05 08.56725 17 10 17.39 +04 29 47.8          23.4 R      568"""


orbfit = ctypes.CDLL('orbfit/liborbfit.so')

# call predict with agbfile, jdate, obscode
orbfit.predict.restype = ctypes.POINTER(ctypes.c_double * 5)
orbfit.predict.argtypes = [ ctypes.c_char_p, ctypes.c_float, ctypes.c_int ]

# class fitradec with mpcfile, abgfile, resfile
orbfit.fitradec.restype = ctypes.c_int
orbfit.fitradec.argtypes = [ ctypes.c_char_p, ctypes.c_char_p, ctypes.c_char_p ]


mpc_file = tempfile.NamedTemporaryFile()
abg_file = tempfile.NamedTemporaryFile()
res_file = tempfile.NamedTemporaryFile()

for line in mpc_lines:
    mpc_file.file.write(line)

mpc_file.file.flush()
mpc_file.seek(0)

result = orbfit.fitradec(ctypes.c_char_p(mpc_file.name), ctypes.c_char_p(abg_file.name),
                         ctypes.c_char_p(res_file.name))

print result

res_file.file.seek(0)
for line in res_file.file.readlines():
    print line.strip()

jd = ctypes.c_float(2456519.500000)
obscode = ctypes.c_int(568)
abg_filename = ctypes.c_char_p(abg_file.name)

# result = ctypes.c_double * 5
result = orbfit.predict(abg_filename, jd, obscode)


ra = ephem.hours(math.radians(result.contents[0]))
dec = ephem.degrees(math.radians(result.contents[1]))
print ra, dec


