# Sends instructions to robot via serial port

import serial as serial_module
from time import sleep

# Input filename for file with coordinate set
input = raw_input("Filename: ")

# Open coordinate file
try:
    with open(input, "r") as fo:
        file = fo.read().split("\n")
        fo.close()
except IOError:
    raise ValueError("File could not be opened")

# Serial interface with robot
serial = serial_module.Serial("/dev/cu.HB-01-DevB", 9600, timeout=2)

print("Connected to: " + serial.portstr)

# this will store the line
line = []

# Send all coordinates in file 
for i, x in enumerate(file):
    if x != "":
        serial.write(x + "#")
        read = serial.read()
        print x
        if i != len(file) - 1:
            while read.find("%") == -1:
                #if read != "":
                #print read
                read = serial.read()

# serial.close()
