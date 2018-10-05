#!/usr/bin/python

import base64
import struct

data = [36,1,
					3, 
					5 ,
					7 ,
					9 ,
					11 ,
					13, 15 ,
					17 ]
fmt = "<" + len(data)*"i"
binary_data = struct.pack(fmt,*data)
print("data: {}, fmt: {}, binary_data: {}".format(data, fmt, binary_data))
print base64.b64encode(binary_data)


decoded = base64.b64decode("AAAQQgAAgD8AAEBA")
print struct.unpack("fff",decoded)
