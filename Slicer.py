import Run_Slicer
import os
import Unit_Converter
import re

def clear():
	os.system("clear")
	
clear()

filename = raw_input("Filename of .svg file: ")

file_preg = re.compile("^.*[\.](svg|SVG)$")
unit_number_preg = re.compile("([0-9]+[\.]?[0-9]*)[^A-Za-z0-9]*([A-Za-z]*)")

if not file_preg.match(filename):
	raise ValueError("Error: Must be .svg file")

slicer_obj = Run_Slicer.Slicer(filename)

input = ""
while input != "yes" and input != "y" and input != "no" and input != "n":
	input = raw_input("\nDo you want to optimize the travel paths? [y/n]: ").lower()

if input == "y" or input == "yes":
	slicer_obj.shorten_toolpaths()

dimensions = [0, 0, 0]
dimension_labels = ["Width", "Height", "Margins"]
dimension_defaults = [216, 280, False]
dimensions_okay = False

print "\nDimensions of standard printer paper: 8.5 x 11 in or 216 x 280 mm"

while dimensions_okay is False:
	print "\n"
	for x in range(0, 3):
		flag = 0
		while flag >= 0:
			if flag != 0:
				print "Invalid dimension (include units maybe?)"
			input = raw_input(dimension_labels[x] + " of paper (press enter to use default): ")
			if input == "":
				dimensions[x] = dimension_defaults[x]
				flag = -1
			else:
				f = True
				try:
					temp = list(unit_number_preg.findall(input)[0])
					dimensions[x] = Unit_Converter.convert(float(temp[0]), temp[1])
				except:
					flag += 1
					f = False
				if f is True:
					flag = -1
	t = "Are the dimensions of " + str(dimensions[0]) + "mm x " + str(dimensions[1]) + "mm with " + ("default" if dimensions[2] is False else str(dimensions[2]) + "mm") + " margins right? [y/n]:"
	temp = raw_input(t).lower()
	if temp == "y" or temp == "yes":
		dimensions_okay = True
		
input = ""
while input != "yes" and input != "y" and input != "no" and input != "n":
	input = raw_input("\nLock aspect ratio scaling? [y/n]: ").lower()

lock_ratio = True if (input == "yes" or input == "y") else False

slicer_obj.set_paper(dimensions[0], dimensions[1], dimensions[2], lock_ratio)

print "\n///////////////////////////////////"
print   "/            Options              /"
print   "///////////////////////////////////"
print   "/ 1.) Continue                    /"
print   "/ 2.) Paper preview               /"
print   "/ 3.) Plot coordinates            /"
print   "/ 4.) Export coordinates          /"
print   "/ 5.) Export scaled coordinates   /"
print   "///////////////////////////////////\n"

print "Choose one of the above options to preview before exporting:"

input = ""
while input != "1":
	input = raw_input("\nOption: ")
	if input == "2":
		slicer_obj.plot_paper()
	elif input == "3":
		slicer_obj.plot_coordinates()
	elif input == "4" or input == "5":
		i = raw_input("File to export to (default is slicer_coordinates.csv): ")
		if i == "":
			if input == "4":
				slicer_obj.export_coordinates()
			else:
				slicer_obj.export_paper()
		else:
			if input == "4":
				slicer_obj.export_coordinates(i + ".csv")
			else:
				slicer_obj.export_paper(i + ".csv")
	elif input != "1":
		print "Invalid option"
		
print "\n"
input = raw_input("Filename, without .acode extension, to export to (default is \"autopen\"): ")
if input == "":
	slicer_obj.export()
else:
	slicer_obj.export(input + ".acode")
