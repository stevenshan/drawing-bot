import re
import os.path
import Vector_to_Coordinate as Slicer_Analyzer
import csv
import matplotlib.pyplot as plt
import math
import Unit_Converter

# Class that contains methods for converting several sets of coordinates to instructions for robot to follow
# Similar in function to a slicer for a 3D printer or CNC machine


class Slicer():

    VALID_TAGS = [
        "svg",
        "circle",
        "ellipse",
        "line",
        "path",
        "polygon",
        "polyline",
        "rect"]
    VALID_UNITS = ["em", "ex", "px", "pt", "pc", "cm", "mm", "in"]
    WIDTH = 216
    HEIGHT = 280
    DEFAULT_MARGIN = 25
    preg_units = False
    FLOAT_DECIMAL_PNTS = 5

    """
	_file_name
	_raw_file
	_vectors_arr
	_coordinate_grid
	_svg_tag
	_width
	_height
	_scale_x
	_scale_y
	_viewbox
	_paper_coordinates
	_paper_dimensions
	"""

    def __init__(self, filename):
        self._file_name = filename
        if not self.preg_units:
            self.preg_units = "|".join(self.VALID_UNITS)
        self._raw_file = ""
        self._vectors_arr = []
        self._coordinate_grid = []
        self._paper_coordinates = False
        self._svt_tag = ""
        self.read(self._file_name)
        self.clean_raw()
        self.analyze_obj()
        self.bound_coordinates()

    def init_analyzer(self, res=False, err=False):
        self.analyzer_root = Slicer_Analyzer.analysis_methods(
            self, res, err, self.FLOAT_DECIMAL_PNTS)
        self.analyzer = {
            "circle": self.analyzer_root.circle,
            "ellipse": self.analyzer_root.ellipse,
            "line": self.analyzer_root.line,
            "path": self.analyzer_root.path,
            "polygon": self.analyzer_root.polygon,
            "polyline": self.analyzer_root.polyline,
            "rect": self.analyzer_root.rect
        }

    def print_paper(self, opt=False):
        return self.print_coordinates(opt, self._paper_coordinates)

    def print_coordinates(self, opt=False, grid=False):
        if grid is False:
            grid = self._coordinate_grid
        if not opt:
            print grid
            return False
        elif opt == "ARRAY":
            temp = []
            for g in grid:
                for c in g:
                    temp.append(c)
            print temp
        elif opt == "RETURN":
            temp = []
            for g in grid:
                for c in g:
                    temp.append(c)
            return temp

    def export(self, file="autopen.acode"):
        if self._paper_coordinates is False:
            raise ValueError("Error: Cannot export, need to set paper size")
        else:
            fl = open(file, 'w')
            fl.truncate()
            current_position = [0.0, 0.0]
            pen_position = 0
            direction = 0
            min_y = 0
            for g in range(0, len(self._paper_coordinates)):
                group = self._paper_coordinates[g]
                for c in range(0, len(group)):
                    coordinate = group[c]
                    y_coord = -1 * self._paper_coordinates[g][c][1]
                    self._paper_coordinates[g][c][1] *= -1
                    if y_coord < min_y:
                        min_y = y_coord
            min_y *= -1
            for g in range(0, len(self._paper_coordinates)):
                group = self._paper_coordinates[g]
                for c in range(0, len(group)):
                    coordinate = group[c]
                    coordinate[1] += min_y
                    dist = self.dist(
                        coordinate[0],
                        coordinate[1],
                        current_position[0],
                        current_position[1])
                    deg = False
                    if dist != 0:
                        if c == 0:
                            fl.write("P 0\n")
                            pen_position = 0
                        if float(coordinate[0] - current_position[0]) == 0:
                            #deg = math.pi / 2
			    deg = 0
			    if (coordinate[1] < current_position[1]):
				deg = math.pi
                        else:
                            slope = float(
                                coordinate[1] - current_position[1]) / float(coordinate[0] - current_position[0])
                            deg = math.pi/2 - abs(math.atan(slope))
			    if (slope <= 0 and coordinate[0] < current_position[0]):
				deg = -deg
                    if deg is not False:
                        temp = direction
                        direction = deg
                        deg -= temp	
			if abs(deg) > math.pi:
			    deg %= math.pi 
                        fl.write("M " +
                                 (str("%." +
                                      str(self.FLOAT_DECIMAL_PNTS) +
                                      "f") %
                                     deg) +
                                 " " +
                                 (str("%." +
                                      str(self.FLOAT_DECIMAL_PNTS) +
                                      "f") %
                                     dist) +
                                 "\n")
			#print([int(x) for x in current_position], " to ", [int(x) for x in coordinate])
                        current_position = coordinate
                    if pen_position == 0:
                        fl.write("P 1\n")
                        pen_position = 1
# 			fl.write("P 0")
            fl.close()

    def export_paper(self, file=False):
        if file is False:
            self.export_coordinates(self.print_paper("RETURN"))
        else:
            self.export_coordinates(self.print_paper("RETURN"), file)

    def export_coordinates(self, coord=False, file="slicer_coordinates.csv"):
        ar = []
        if coord is False:
            ar = self.print_coordinates("RETURN")
        else:
            ar = coord
        fl = open(file, 'w')

        writer = csv.writer(fl)
        writer.writerow(['x', 'y'])  # if needed
        for values in ar:
            writer.writerow(values)

        fl.close()

    def plot_paper(self):
        self.plot_coordinates(True, self._paper_dimensions)

    def plot_coordinates(self, paper=False, ranges=False):
        ar = []
        if paper is False:
            ar = self.print_coordinates("RETURN")
        else:
            ar = self.print_paper("RETURN")
        x = []
        y = []
        for i in ar:
            x.append(i[0])
            y.append(-i[1])
        plt.plot(x, y, "ro")
        if ranges is not False:
            plt.xlim(ranges[0])
            plt.ylim([-ranges[1][1], -ranges[1][0]])
        plt.gca().set_aspect('equal', adjustable='box')
        plt.show()

    def read(self, filename):
        if os.path.exists(filename):
            try:
                with open(filename, "r") as fo:
                    self._raw_file = fo.read()
                    fo.close()
            except IOError:
                raise ValueError("File could not be opened")
        else:
            raise ValueError("SVG file not found: ", filename)

    def clean_raw(self):
        self._raw_file = ' '.join(self._raw_file.replace("\n", "").split())
        self._vectors_arr = re.findall(
            r'(<[^a-zA-Z/]?(' + "|".join(self.VALID_TAGS) + ')[^>]*>)', self._raw_file)
        if(len(self._vectors_arr) == 0):
            raise ValueError("No valid SVG features found")

    @staticmethod
    def convert_to_mm(num, unit):
        return Unit_Converter.convert(num, unit)

    def init_svg(self, obj):
        arr_args = re.findall(
            r'([a-zA-Z]+)[ ]*=[ ]*[\'\"]([^(\"|\')]*)[\'\"]', obj[0])
        args = {}
        for x in arr_args:
            args[x[0].lower()] = x[1]
        self._svt_tag = obj[0]
        w = 0
        h = 0
        # Dimensions not specified | Response: use default dimensions _WIDTH and
        # _HEIGHT
        if "width" not in args or "height" not in args:
            print "Warning: SVG dimensions not found. Proceeding with default dimensions of ", self.WIDTH, "mm x ", self.HEIGHT, "mm"
            w = str(self.WIDTH) + "mm"
            h = str(self.HEIGHT) + "mm"
        else:  # Otherwise proceed with given dimensions
            w = args["width"].lower()
            h = args["height"].lower()

        f = 0
        if not re.match(
            "^([0-9]+[\.]?[0-9]*)(" + self.preg_units + ")$",
                w):  # if width attribute is invalid
            p = re.match("^([0-9]+[\.]?[0-9]*)", w)
            if p:  # if units are invalid
                w = p.group(0) + "xx"
                f = 1
            else:
                print "Warning: Invalid SVG dimensions. Proceeding with default dimensions of ", self.WIDTH, "mm x ", self.HEIGHT, "mm"
                w = str(self.WIDTH) + "mm"
                h = str(self.HEIGHT) + "mm"
                f = 2

        if f != 2:
            if not re.match(
                "^([0-9]+[\.]?[0-9]*)(" + self.preg_units + ")$",
                    h):  # if height attribute is invalid
                p = re.match("^([0-9]+[\.]?[0-9]*)", h)
                if p:
                    h = p.group(0) + "xx"
                    if f == 0:
                        w = re.match("^([0-9]+[\.]?[0-9]*)", w).group(0) + "xx"
                    f = 1
                else:
                    print "Warning: Invalid SVG dimensions. Proceeding with default dimensions of ", self.WIDTH, "mm x ", self.HEIGHT, "mm"
                    w = str(self.WIDTH) + "mm"
                    h = str(self.HEIGHT) + "mm"
                    f = 2
            elif f == 1:
                h = re.match("^([0-9]+[\.]?[0-9]*)", h).group(0) + "xx"

        preg = re.compile("([0-9]+[\.]?[0-9]*)(" + self.preg_units + "|xx)")
        w = list(preg.findall(w)[0])
        h = list(preg.findall(h)[0])
        w[0] = float(w[0])
        h[0] = float(h[0])

        if f == 1:
            print "Warning: SVG dimension units invalid. Removing units and scaling to default margin of 1\" on ", self.WIDTH, "mm x ", self.HEIGHT, "mm page"
            if h[0] / w[0] > (self.HEIGHT - 2 * self.DEFAULT_MARGIN) / \
                    (self.WIDTH - 2 * self.DEFAULT_MARGIN):  # Scale based on height
                temp = h[0]
                h[0] = self.HEIGHT - 2 * self.DEFAULT_MARGIN
                w[0] *= h[0] / temp
            else:  # Scale based on width
                temp = w[0]
                w[0] = self.WIDTH - 2 * self.DEFAULT_MARGIN
                h[0] *= w[0] / temp
            print "         Assigned dimension of ", w[0], "mm x ", h[0], "mm"
        if f == 0:  # Make sure everything is in millimeters
            if w[1] != "mm":
                print "Warning: SVG width dimensions converted from " + w[1] + " to mm"
                w[0] = Slicer.convert_to_mm(w[0], w[1])
                w[1] = "mm"
            if h[1] != "mm":
                print "Warning: SVG height dimensions converted from " + h[1] + " to mm"
                h[0] = Slicer.convert_to_mm(h[0], h[1])
                h[1] = "mm"

        if "viewbox" in args:
            self._viewbox = re.findall("[-]?[0-9]+[\.]*[0-9]*", args["viewbox"])
            if(len(self._viewbox) != 4):
                print "Warning: Ignoring viewbox attribute"
                self._viewbox = [0, 0, w[0], h[0]]
            else:
                for i in range(len(self._viewbox)):
                    self._viewbox[i] = float(self._viewbox[i])
        else:
            print "Warning: Viewbox not specified in SVG"
            self._viewbox = [0, 0, w[0], h[0]]

        self._scale_x = w[0] / self._viewbox[2]
        self._scale_y = h[0] / self._viewbox[3]
        self._width = w[0]
        self._height = h[0]
        self._viewbox[2] += self._viewbox[0]
        self._viewbox[3] += self._viewbox[1]

        self.init_analyzer(2 / (self._scale_x + self._scale_y))

    def analyze_obj(self):
        j = 0
        while self._vectors_arr[j][1] != "svg":
            j += 1
        self.init_svg(self._vectors_arr[j])
        i = 0
        for obj in self._vectors_arr:
            if i != j:
                if obj[1] == "svg":
                    raise ValueError("Error: Multiple SVGs found")
                arr_args = re.findall(
                    r'([a-zA-Z]+)[ ]*=[ ]*[\'\"]([^(\"|\')]*)[\'\"]', obj[0])
                args = {}
                for x in arr_args:
                    args[x[0].lower()] = x[1]
                if len(args) > 0:
                    self._coordinate_grid.extend(self.analyzer[obj[1]](args))
            i += 1

    def coordinate_check_out(self, coord):
        if coord[0] < self._viewbox[0] or coord[1] < self._viewbox[1] or coord[0] > self._viewbox[2] or coord[1] > self._viewbox[3]:
            return True
        return False

    def locate_intersection_bound(self, c1, c2):
        if c2[0] == c1[0]:  # slope is undefined
            if c1[0] > self._viewbox[0] and c1[0] < self._viewbox[2]:
                y_min = min([c1[1], c2[1]])
                y_max = max([c1[1], c2[1]])
                intersections = []
                if y_min < self._viewbox[1] and y_max > self._viewbox[1]:
                    intersections.append([c1[0], self._viewbox[1]])
                if y_min < self._viewbox[3] and y_max > self._viewbox[3]:
                    intersections.append([c1[0], self._viewbox[3]])
                if len(intersections) != 0:
                    return intersections
            return False
        else:
            intersections = []
            # y = mx + b
            m = float(c2[1] - c1[1]) / float(c2[0] - c1[0])
            b = -m * c1[0] + c1[1]

            y_min = min([c1[1], c2[1]])
            y_max = max([c1[1], c2[1]])

            # check with left boundary
            y = m * self._viewbox[0] + b
            if y > self._viewbox[1] and y < self._viewbox[3] and y_min < y and y_max > y:
                intersections.append([self._viewbox[0], y])

            # check with right boundary
            y = m * self._viewbox[2] + b
            if y > self._viewbox[1] and y < self._viewbox[3] and y_min < y and y_max > y:
                intersections.append([self._viewbox[2], y])

            if m != 0:
                # check with top boundary
                x = (self._viewbox[1] * 1.0 - b) / m
                x_min = min([c1[0], c2[0]])
                x_max = max([c1[0], c2[0]])
                if x > self._viewbox[0] and x < self._viewbox[2] and x_min < x and x_max > x:
                    intersections.append([x, self._viewbox[1]])

                # check with bottom boundary
                x = (self._viewbox[3] * 1.0 - b) / m
                if x > self._viewbox[0] and x < self._viewbox[2] and x_min < x and x_max > x:
                    intersections.append([x, self._viewbox[3]])

            if len(intersections) == 0:
                return False
            return intersections

    @staticmethod
    def dist(x1, y1, x2, y2):
        return math.sqrt(math.pow(x2 - x1, 2) + math.pow(y2 - y1, 2))

    def bound_coordinates(self):
        index = -1
        _len = len(self._coordinate_grid)
        while index + 1 < _len:
            index += 1
            x = self._coordinate_grid[index]
            current_pnt_out = False
            if len(x) > 1:
                y = 0
                cond = False
                while y < cond or cond is False:
                    b = [self.coordinate_check_out(
                        x[y]), self.coordinate_check_out(x[y + 1])]
                    intersections = self.locate_intersection_bound(
                        x[y], x[y + 1])
                    intersections_disp = []
                    if intersections:
                        for w in intersections:
                            intersections_disp.append(
                                self.dist(x[y][0], x[y][1], w[0], w[1]))
                    if b[0] and b[1]:
                        if not intersections:
                            if y + 2 < len(x):
                                self._coordinate_grid.insert(
                                    index + 1, x[y + 1:])
                                if y != 0:
                                    self._coordinate_grid[index] = x[0:y]
                                else:
                                    del self._coordinate_grid[index]
                                    index -= 1
                                x = x[0:y]
                            else:
                                self._coordinate_grid[index] = self._coordinate_grid[index][0:y]
                                x = self._coordinate_grid[index]
                        elif len(intersections) != 2:
                            raise ValueError("Error: Invalid intersection (0)")
                        else:
                            min_index = intersections_disp.index(
                                min(intersections_disp))
                            x[y] = intersections[min_index]
                            x.insert(y + 1,
                                     intersections[(0 if min_index == 1 else 1)])
                            self._coordinate_grid[index] = x[0:y + 1]
                            self._coordinate_grid.insert(index + 1, x[y + 1:])
                            x = x[0:y + 1]
                    elif b[1] and not b[0]:  # Goes in to out
                        if intersections == False or len(intersections) != 1:
                            raise ValueError("Error: Invalid intersection (1)")
                        else:
                            if y + 1 >= len(x) - 1:
                                x[y + 1] = intersections[0]
                                self._coordinate_grid[index] = x
                            else:
                                x.insert(y + 1, intersections[0])
                                self._coordinate_grid[index] = x[0:y + 2]
                                self._coordinate_grid.insert(
                                    index + 1, x[y + 2:])
                                x = x[0:y + 2]
                    elif b[0] and not b[1]:  # Goes out to in
                        if intersections == False or len(intersections) != 1:
                            raise ValueError("Error: Invalid intersection (2)")
                        else:
                            x[y] = intersections[0]
                            self._coordinate_grid[index] = x

                    cond = len(x) - 1
                    y += 1
            else:
                del self._coordinate_grid[index]
                index -= 1
            _len = len(self._coordinate_grid)

        x = 0
        _len = len(self._coordinate_grid)
        while x < len(self._coordinate_grid):
            if len(self._coordinate_grid[x]) == 0:
                del self._coordinate_grid[x]
                x -= 1
            _len = len(self._coordinate_grid)
            x += 1

    def shorten_toolpaths(self):
        index = 0
        min_dist = False
        for x in range(0, len(self._coordinate_grid)):
            temp = min([
                self.dist(0, 0, self._coordinate_grid[x][0][0], self._coordinate_grid[x][0][1]),
                self.dist(0, 0, self._coordinate_grid[x][-1][0], self._coordinate_grid[x][-1][1])
            ])

            if min_dist == False or temp < min_dist:
                index = x
                min_dist = temp

        for x in range(0, len(self._coordinate_grid) - 1):
            index = 0
            min_dist = False
            for y in range(x + 1, len(self._coordinate_grid)):
                temp = min([self.dist(self._coordinate_grid[x][0][0],
                                      self._coordinate_grid[x][0][1],
                                      self._coordinate_grid[y][0][0],
                                      self._coordinate_grid[y][0][1]),
                            self.dist(self._coordinate_grid[x][0][0],
                                      self._coordinate_grid[x][0][1],
                                      self._coordinate_grid[y][-1][0],
                                      self._coordinate_grid[y][-1][1]),
                            self.dist(self._coordinate_grid[x][-1][0],
                                      self._coordinate_grid[x][-1][1],
                                      self._coordinate_grid[y][0][0],
                                      self._coordinate_grid[y][0][1]),
                            self.dist(self._coordinate_grid[x][-1][0],
                                      self._coordinate_grid[x][-1][1],
                                      self._coordinate_grid[y][-1][0],
                                      self._coordinate_grid[y][-1][1]),
                            ])
                if min_dist == False or temp < min_dist:
                    index = x
                    min_dist = temp
            temp = self._coordinate_grid[x + 1]
            self._coordinate_grid[x + 1] = self._coordinate_grid[index]
            self._coordinate_grid[index] = temp

    def set_paper(
            self,
            width=False,
            height=False,
            margins=False,
            lock_ratio=False):
        if width is False:
            width = self.WIDTH
        if height is False:
            height = self.HEIGHT
        if margins is False:
            margins = self.DEFAULT_MARGIN

        print_width = width - margins * 2
        print_height = height - margins * 2

        scale_x = print_width / self._width
        scale_y = print_height / self._height

        x_disp = 0
        y_disp = 0

        if lock_ratio is True:
            scl = min([scale_x, scale_y])
            scale_x = scl
            scale_y = scl
            x_disp = (print_width - scale_x * self._width) / 2.0
            y_disp = (print_height - scale_y * self._height) / 2.0

        if scale_x < 1 or scale_y < 1:
            print "Warning: Image was scaled down to fit paper"

        self._paper_coordinates = []
        self._paper_dimensions = [[0, width], [0, height]]

        for x in range(0, len(self._coordinate_grid)):
            self._paper_coordinates.append([])
            for c in range(0, len(self._coordinate_grid[x])):
                self._paper_coordinates[-1].append([])
                self._paper_coordinates[x][c] = [
                    (self._coordinate_grid[x][c][0] -
                     self._viewbox[0]) *
                    self._scale_x *
                    scale_x +
                    margins +
                    x_disp,
                    (self._coordinate_grid[x][c][1] -
                     self._viewbox[1]) *
                    self._scale_y *
                    scale_y +
                    margins +
                    y_disp]
