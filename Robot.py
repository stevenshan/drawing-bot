# Sends instructions to robot via serial port

import serial as serial_module
from time import sleep

input = raw_input("Filename: ")

try:
	with open(input, "r") as fo:
		file = fo.read().split("\n")
		fo.close()
except IOError:
	raise ValueError("File could not be opened")

serial = serial_module.Serial("/dev/cu.HB-01-DevB",9600,timeout = 2)

print("Connected to: " + serial.portstr)

#this will store the line
line = []


for i,x in enumerate(file):
	if x != "":
		serial.write(x + "#")
		read = serial.read()
		print x
		if i != len(file) - 1:
			while read.find("%") == -1:
	# 			if read != "":
	# 				print read
				read = serial.read()

# serial.close()