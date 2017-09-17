import re
import math

"""
Sources for approximation of arc to bezier curve:
	1.) https://www.spaceroots.org/documents/ellipse/elliptical-arc.pdf
	2.) https://www.joecridge.me/content/pdf/bezier-arcs.pdf
	
Source code for converting arc to bezier curves:
	1.) https://github.com/fontello/svgpath
"""

class analysis_methods():

	def __init__ (self,parent,res = False,err = False,decimal = 5):
		self._parent = parent
		if res == False:
			self.RESOLUTION = 1
		else:
			self.RESOLUTION = res
		if err == False:
			self.RESOLUTION_ERR = 0.01
		else:
			self.RESOLUTION_ERR = err
		self.FLOAT_DECIMAL_PNTS = decimal

	VALID_ARGUMENTS = [
		["cx", "cy", "r"], #circle
		["cx", "cy", "rx", "ry"], #ellipse
		["x1", "y1", "x2", "y2"], #line
		["d"], #path
		[], #polygon | same as polyline
		["points"], #polyline
		["x", "y", "width", "height"], #rect
	]
		
	@staticmethod
	def check_args (const, args):
		for i in analysis_methods.VALID_ARGUMENTS[const]:
			if i not in args:
				raise ValueError("Invalid SVG element arguments (0)")
			elif args[i] == "":
				raise ValueError("Empty SVG element argument")
		
	@staticmethod
	def float_args (const,args):
		for i in analysis_methods.VALID_ARGUMENTS[const]:
			args[i] = float(args[i])
			if i not in args:
				raise ValueError("Invalid SVG element arguments (1)")
			elif args[i] == "":
				raise ValueError("Empty SVG element argument")
			else:
				args[i] = float(args[i])
				
	@staticmethod
	def dist (x1, y1, x2, y2):
		return math.sqrt(math.pow(x2-x1, 2) + math.pow(y2-y1,2))
	
	def circle (self,args):
		const = 0
		analysis_methods.float_args(const, args)
		iv = float(self.RESOLUTION) / args["r"]
		deg = 0
		point_array = [[]]
		while deg < 2 * math.pi:
			point_array[0].append([math.cos(deg) * args["r"] + args["cx"], math.sin(deg) * args["r"] + args["cy"]])
			deg += iv
		return point_array
		
	@staticmethod
	def ellipse_coord (a,b,x):
		return math.sqrt((math.pow(float(a)*b,2) - math.pow(float(b)*x,2))/math.pow(a, 2))
		
	def approx_ellipse(self,n1,n2,a,b,xi,yi):
		xx = xi-(n2+n1)/2.0
		yx = self.ellipse_coord(a,b,xx)
		dist = self.dist(xi, yi, xx, yx)
		if math.fabs(dist - self.RESOLUTION) / self.RESOLUTION <= self.RESOLUTION_ERR: # condition
			return [xx, yx]
		else:
			if dist > self.RESOLUTION:
				n2 = (n2+n1)/2.0
			else:
				n1 = (n2+n1)/2.0
			return self.approx_ellipse(n1, n2, a, b, xi, yi)
	
	def ellipse (self,args):
		const = 1
		analysis_methods.float_args(const, args)
		xi = args["rx"] / 2.0
		point_arr = [[xi, 0]]
		while point_arr[-1][0] > 0:
			pnt = self.approx_ellipse(0, self.RESOLUTION, args["rx"] / 2.0, args["ry"] / 2.0, point_arr[-1][0], point_arr[-1][1])
			point_arr.append(pnt)
		point_arr[-1] = [0, args["ry"]/2] # Quadrant 1 of ellipse
			
		l = len(point_arr) - 1
		i = l
		while i >= 0: # Flip quadrant 1 over y-axis
			if i != l:
				point_arr.append([point_arr[i][0] * -1, point_arr[i][1]])
			i -= 1
			
		l = len(point_arr) - 1
		i = l
		while i >= 0: # Flip quadrant 1 over x-axis
			if i != l and i != 0:
				point_arr.append([point_arr[i][0] + args["cx"], point_arr[i][1] * -1 + args["cy"]])
			point_arr[i][0] += args["cx"]
			point_arr[i][1] += args["cy"]
			i -= 1
		
		return [point_arr]
		
	@staticmethod
	def line (args):
		const = 2
		analysis_methods.check_args(const, args)
		return [[
			[float(args["x1"]), float(args["y1"])],
			[float(args["x2"]), float(args["y2"])],
		]]
		
	@staticmethod
	def path_break_args (arg):
		ret = re.findall("([-]?[0-9]+[\.]?[0-9]*)[^-0-9]*", arg)
		i = 0
		for x in ret:
			ret[i] = float(ret[i])
			i += 1
		return ret
		
	@staticmethod
	def calc_cubic_bezier (t,arg,w):
		Ax = ((1 - t) * w[0] ) + (t * arg[0])
		Ay = ((1 - t) * w[1] ) + (t * arg[1])
		Bx = ((1 - t) * arg[0]) + (t * arg[2])
		By = ((1 - t) * arg[1]) + (t * arg[3])
		Cx = ((1 - t) * arg[2]) + (t * arg[4])
		Cy = ((1 - t) * arg[3]) + (t * arg[5])
		Dx = ((1 - t) * Ax) + (t * Bx)
		Dy = ((1 - t) * Ay) + (t * By)
		Ex = ((1 - t) * Bx) + (t * Cx)
		Ey = ((1 - t) * By) + (t * Cy)
		return [((1 - t) * Dx) + (t * Ex), ((1 - t) * Dy) + (t * Ey)]
		
	@staticmethod
	def calc_quadratic_bezier (t,arg,w):
		x = (w[0] - 2 * arg[0] + arg[2])*math.pow(t,2) + (2 * arg[0] - 2 * w[0]) * t + w[0]
		y = (w[1] - 2 * arg[1] + arg[3])*math.pow(t,2) + (2 * arg[1] - 2 * w[1]) * t + w[1]
		return [x,y]
		
	def approx_bezier (self,type,t1,t2,arg,x,y,xi,yi,res = False):
		if res == False:
			res = self.RESOLUTION_ERR
		coord = getattr(self, "calc_" + type + "_bezier")((t2+t1)/2.0, arg, [xi,yi])
		dist = self.dist(x,y,coord[0],coord[1])
		if t2 == 1 and dist < self.RESOLUTION:
			return [1, [arg[-2], arg[-1]]]
		if math.fabs(dist - self.RESOLUTION) / self.RESOLUTION <= res: # condition
			return [(t2+t1)/2.0,coord]
		else:
			if dist > self.RESOLUTION:
				t2 = (t2+t1)/2.0
			else:
				t1 = (t2+t1)/2.0
			return self.approx_bezier(type, t1, t2, arg, x, y, xi, yi, res)
			
	def generate_bezier (self, w, arg, first_coordinate):
		#Algorithm Backup: Approximate all points
		t = 0
		dist = 0
		last_coordinate = [w[-1][0], w[-1][1]]
		while t == 0 or dist > self.RESOLUTION:
			approximation = self.approx_bezier("cubic", t, 1, arg, last_coordinate[0], last_coordinate[1], first_coordinate[0], first_coordinate[1])
			t = approximation[0]
			last_coordinate = approximation[1]
			w.append(approximation[1])
			dist = self.dist(arg[4], arg[5], approximation[1][0], approximation[1][1])
		w.append([arg[4], arg[5]])
						
#		Algorithm 2: Use approximation from first interval to extrapolate
# 		approximation = self.approx_bezier("cubic", 0, 1, arg, w[-1][0], w[-1][1], w[-1][0], w[-1][1])[0]
# 		dist = 0
# 		i = 1
# 		while i == 1 or dist > self.RESOLUTION:
# 			new_coordinate = self.calc_cubic_bezier(i * approximation, arg, first_coordinate)
# 			dist = self.dist(arg[4], arg[5], new_coordinate[0], new_coordinate[1])
# 			w.append(new_coordinate)
# 			i += 1
# 		w.append([arg[4], arg[5]])
		
	@staticmethod
	def unit_vector_angle(ux, uy, vx, vy):
		sign = -1 if (ux * vy - uy * vx < 0) else 1
		dot = ux * vx + uy * vy
		if dot > 1.0:
			dot = 1.0
		if dot < -1.0:
			dot = -1.0
		return sign * math.acos(dot)

	@staticmethod
	def get_arc_center(x1, y1, x2, y2, fa, fs, rx, ry, sin_phi, cos_phi):
		TAU = math.pi * 2.0
		x1p =  cos_phi*(x1-x2)/2.0 + sin_phi*(y1-y2)/2.0
		y1p = -sin_phi*(x1-x2)/2.0 + cos_phi*(y1-y2)/2.0

		rx_sq = rx * rx
		ry_sq = ry * ry
		x1p_sq = x1p * x1p
		y1p_sq = y1p * y1p

		radicant = (rx_sq * ry_sq) - (rx_sq * y1p_sq) - (ry_sq * x1p_sq)
	
		if (radicant < 0):
			radicant = 0;

		radicant /=   (rx_sq * y1p_sq) + (ry_sq * x1p_sq);
		radicant = math.sqrt(radicant) * (-1 if fa == fs else 1)

		cxp = radicant *  rx/ry * y1p
		cyp = radicant * -ry/rx * x1p


		cx = cos_phi*cxp - sin_phi*cyp + (x1+x2)/2
		cy = sin_phi*cxp + cos_phi*cyp + (y1+y2)/2

		v1x =  (x1p - cxp) / rx
		v1y =  (y1p - cyp) / ry
		v2x = (-x1p - cxp) / rx
		v2y = (-y1p - cyp) / ry

		theta1 = analysis_methods.unit_vector_angle(1, 0, v1x, v1y)
		delta_theta = analysis_methods.unit_vector_angle(v1x, v1y, v2x, v2y)
		
		if (fs == False and delta_theta > 0):
			delta_theta -= TAU

		if (fs == True and delta_theta < 0):
			delta_theta += TAU
			
		return [cx, cy, theta1, delta_theta]

	@staticmethod
	def approximate_unit_arc(theta1, delta_theta):
		alpha = 4.0/3.0 * math.tan(delta_theta/4.0)
		x1 = math.cos(theta1)
		y1 = math.sin(theta1)
		x2 = math.cos(theta1 + delta_theta)
		y2 = math.sin(theta1 + delta_theta)
		
		return [x1, y1, x1 - y1*alpha, y1 + x1*alpha, x2 + y2*alpha, y2 - x2*alpha, x2, y2]

	@staticmethod
	def arc_to_bezier (x1, y1, x2, y2, fa, fs, rx, ry, phi):
		TAU = math.pi * 2.0

		sin_phi = math.sin(phi * TAU / 360)
		cos_phi = math.cos(phi * TAU / 360)

		x1p =  cos_phi*(x1-x2)/2 + sin_phi*(y1-y2)/2
		y1p = -sin_phi*(x1-x2)/2 + cos_phi*(y1-y2)/2

		if (x1p == 0 and y1p == 0) or (rx == 0 or ry == 0):
			return [];

		rx = math.fabs(rx)
		ry = math.fabs(ry)

		y = (x1p * x1p) / (rx * rx) + (y1p * y1p) / (ry * ry);
		if (y > 1):
			rx *= math.sqrt(y)
			ry *= math.sqrt(y)

		cc = analysis_methods.get_arc_center(x1, y1, x2, y2, fa, fs, rx, ry, sin_phi, cos_phi)
		
		result = []
		theta1 = cc[2]
		delta_theta = cc[3]
	
		segments = max([int(math.ceil(math.fabs(delta_theta) / (TAU / 4))), 1])
		delta_theta /= segments
		for i in range(0, segments):
			result.append(analysis_methods.approximate_unit_arc(theta1, delta_theta))
			theta1 += delta_theta
											
		return_arr = []
		for curve in result:
			for i in range(0, len(curve), 2):
				x = curve[i]
				y = curve[i + 1]

				x *= rx
				y *= ry

				xp = cos_phi*x - sin_phi*y
				yp = sin_phi*x + cos_phi*y

				curve[i] = xp + cc[0]
				curve[i + 1] = yp + cc[1]
			return_arr.append(curve)
		return return_arr
		
	def convert_rel (self, rel, arg, w):
		if rel:
			rel = False
			i = 0
			for c in arg:
				if i % 2 == 0: # x-axis
					arg[i] += w[-1][0]
				else: # y-axis
					arg[i] += w[-1][1]
				i += 1
		
	def path (self,args):
		const = 3
		analysis_methods.check_args(const, args)
		_args = re.findall("(([MLHVCSQTAZmlhvcsqtaz])[-0-9, \.]*)", args["d"])
		array = []
		w = False
		last_command = False
		index = -1
		skip_flag = False
		for x in _args:
			index += 1
			x = list(x)
			rel = False
			if not x[1].isupper():
				x[1] = x[1].upper()
				rel = True
			
			arg = analysis_methods.path_break_args(x[0][1:])
			if x[1] == "M":
				if w == False:
					w = []
				elif len(w) > 0:
					array.append(w)
					w = []
				if len(arg) % 2 != 0:
					raise ValueError("Invalid 'M' value in path SVG element")
				else:
					if len(arg) == 2:
						if not rel or len(array) == 0 or len(array[-1]) == 0:
							w.append([arg[0], arg[1]])
						else:
							w.append([arg[0] + array[-1][-1][0], arg[1] + array[-1][-1][1]])
					else:
						cmd = "m" if rel else "M"
						_args.insert(index+1,[cmd + " "  + " ".join([str('{:.' + str(self.FLOAT_DECIMAL_PNTS) + 'f}').format(n) for n in arg[0:2]]),cmd])
						cmd = "l" if rel else "L"
						for u in range(2, len(arg), 2):
							_args.insert(index+1+u/2,[cmd + " " + " ".join([str('{:.' + str(self.FLOAT_DECIMAL_PNTS) + 'f}').format(n) for n in arg[u:u+2]]),cmd])
			elif x[1] == "Z" and w != False:
				w_len = len(w)
				if w_len >= 3:
					w.append(w[0])
				array.append(w)
				w = []
			else:
				if w == False:
					w = []
					w.append([0.0, 0.0])
				elif len(w) == 0:
					w.append(arr[-1][-1])
					
				if x[1] == "H":
					if len(arg) != 1:
						raise ValueError("Invalid 'H' value in path SVG element")
					else:
						if not rel:
							w.append([arg[0], w[-1][1]])
						else:
							w.append([w[-1][0] + arg[0], w[-1][1]])
				elif x[1] == "V":
					if len(arg) != 1:
						raise ValueError("Invalid 'V' value in path SVG element")
					else:
						if not rel:
							w.append([w[-1][0], arg[0]])
						else:
							w.append([w[-1][0], w[-1][1] + arg[0]])
				elif x[1] == "A":
					if len(arg) != 7:
						raise ValueError("Invalid 'A' value in path SVG element")
					else:
						if rel:
							arg[5] += w[-1][0]
							arg[6] += w[-1][1]
							
						arg[3] = True if arg[3] == 1 else False
						arg[4] = True if arg[4] == 1 else False
						
						bezier_approximations = self.arc_to_bezier(w[-1][0], w[-1][1], arg[5], arg[6], arg[3], arg[4], arg[0], arg[1], arg[2])
						cmd = "M " + (str("%." + str(self.FLOAT_DECIMAL_PNTS) + "f") % w[-1][0]) + " " + (str("%." + str(self.FLOAT_DECIMAL_PNTS) + "f") % w[-1][1]) + " "
						
						for b in bezier_approximations:
							b = [str("%." + str(self.FLOAT_DECIMAL_PNTS) + "f") % number for number in b]						
							cmd += "C " + " ".join(b[2:]) + " "
						w.extend(self.path({"d": cmd})[0][1:])	
				else:
					if x[1] == "C" or x[1] == "S":
						if ((len(arg) < 6 or len(arg) % 6 != 0) and x[1] == "C") or ((len(arg) < 4 or len(arg) % 4 != 0) and x[1] == "S"):
							raise ValueError("Invalid 'C/S' value in path SVG element")
						else:
							if x[1] == "S":
								if len(arg) > 4:
									skip_flag = True
									temp = "S"
									temp = temp if not rel else temp.lower()
									for u in range(0, len(arg), 4):
										_args.insert(index+1+u/4,[temp + " " + " ".join([str('{:.' + str(self.FLOAT_DECIMAL_PNTS) + 'f}').format(n) for n in arg[u:u+4]]),temp])
								else:
									self.convert_rel(rel, arg, w)
									if last_command[0] != "C":
										x[1] = "Q"
									else:
										x1 = w[-1][0] * 2 - last_command[1][2]
										y1 = w[-1][1] * 2 - last_command[1][3]
										arg = [x1, y1] + arg
										x[1] = "C"
						
							if x[1] == "C":
								if len(arg) > 6:
									skip_flag = True
									temp = "C"
									temp = temp if not rel else temp.lower()
									for u in range(0, len(arg), 6):
										_args.insert(index+1+u/6,[temp + " " + " ".join([str('{:.' + str(self.FLOAT_DECIMAL_PNTS) + 'f}').format(n) for n in arg[u:u+6]]),temp])
								else:
									self.convert_rel(rel, arg, w)
									first_coordinate = [w[-1][0], w[-1][1]]
									self.generate_bezier(w, arg, first_coordinate)
				
					if x[1] == "Q" or x[1] == "T":
						if ((len(arg) < 4 or len(arg) % 4 != 0) and x[1] == "Q") or ((len(arg) < 2 or len(arg) % 2 != 0) and x[1] == "T"):
							raise ValueError("Invalid 'Q' value in path SVG element")
						else:
					
							if x[1] == "T":
								if len(arg) > 2:
									skip_flag = True
									temp = "T"
									temp = temp if not rel else temp.lower()
									for u in range(0, len(arg), 2):
										_args.insert(index+1+u/2,[temp + " " + " ".join([str('{:.' + str(self.FLOAT_DECIMAL_PNTS) + 'f}').format(n) for n in arg[u:u+2]]),temp])
								else:
									self.convert_rel(rel, arg, w)
									if last_command[0] != "Q":
										x[1] = "L"
									else:
										x1 = w[-1][0] * 2 - last_command[1][0]
										y1 = w[-1][1] * 2 - last_command[1][1]
										arg = [x1, y1] + arg
										x[1] = "Q"
							
							if x[1] == "Q":
								if len(arg) > 4:
									skip_flag = True
									temp = "Q"
									temp = temp if not rel else temp.lower()
									for u in range(0, len(arg), 4):
										_args.insert(index+1+u/4,[temp + " " + " ".join([str('{:.' + str(self.FLOAT_DECIMAL_PNTS) + 'f}').format(n) for n in arg[u:u+4]]),temp])
								else:
									self.convert_rel(rel, arg, w)
									first_coordinate = [w[-1][0], w[-1][1]]

									# Algorithm Backup: Approximate all points
									t = 0
									dist = 0
									last_coordinate = [w[-1][0], w[-1][1]]
									while t == 0 or dist > self.RESOLUTION:
										approximation = self.approx_bezier("quadratic", t, 1, arg, last_coordinate[0], last_coordinate[1], first_coordinate[0], first_coordinate[1])
										t = approximation[0]
										last_coordinate = approximation[1]
										w.append(approximation[1])
										dist = self.dist(arg[2], arg[3], approximation[1][0], approximation[1][1])
									w.append([arg[2], arg[3]])
						
									# Algorithm 2: Use approximation from first interval to extrapolate
	# 								approximation = self.approx_bezier("quadratic", 0, 1, arg, w[-1][0], w[-1][1], w[-1][0], w[-1][1])[0]
	# 								dist = 0
	# 								i = 1
	# 								while i == 1 or dist > self.RESOLUTION:
	# 									new_coordinate = self.calc_quadratic_bezier(i * approximation, arg, first_coordinate)
	# 									dist = self.dist(arg[2], arg[3], new_coordinate[0], new_coordinate[1])
	# 									w.append(new_coordinate)
	# 									i += 1
	# 								w.append([arg[2], arg[3]])
						
					if x[1] == "L":
						if len(arg) > 2 and len(arg) % 2  == 0:
							cmd = "l" if rel else "L"
							_args.insert(index+1,[cmd + " " + " ".join([str('{:.' + str(self.FLOAT_DECIMAL_PNTS) + 'f}').format(n) for n in arg[0:2]]),cmd])
							cmd = "l" if rel else "L"
							for u in range(2, len(arg), 2):
								_args.insert(index+1+u/2,[cmd + " " + " ".join([str('{:.' + str(self.FLOAT_DECIMAL_PNTS) + 'f}').format(n) for n in arg[u:u+2]]),cmd])
						elif len(arg) != 2:
							raise ValueError("Invalid 'L' value in path SVG element")
						else:
							self.convert_rel(rel, arg, w)
							w.append([arg[0], arg[1]])	
			if not skip_flag:
				last_command = [x[1], arg]	
		if w != False and len(w) > 0:
			array.append(w)
		return array
		
	@staticmethod	
	def polygon (args):
		const = 4
 		point_arr = analysis_methods.polyline(args)
 		point_arr.append(point_arr[0])
 		return [point_arr]

	@staticmethod
	def polyline (args):
		const = 5
		analysis_methods.check_args(const, args)
		_args = args["points"]
		preg = re.compile("^([ ]*([-]?[0-9]+[\.]?[0-9]*,[-]?[0-9]+[\.]?[0-9]*)[ ]*)*$")
		if not preg.search(_args):
			raise ValueError("Invalid polyline SVG element argument value")
		point_arr = re.findall("([0-9]+[\.]?[0-9]*,[0-9]+[\.]?[0-9]*)", _args)
		i = 0
		for x in point_arr:
			j = x.find(",")
			point_arr[i] = [float(x[0:j]), float(x[j+1:])]
			i += 1
 		return [point_arr]
	
	@staticmethod
	def rect (args):
		const = 6
		analysis_methods.check_args(const, args)
		_args = [args["x"], args["y"], args["width"], args["height"]]
		preg = re.compile("^[-]?[0-9]+[\.]?[0-9]*$")
		i = 0
		for x in _args:
			if(not preg.match(x)):
				raise ValueError("Invalid rect SVG element argument value")
			else:
				_args[i] = float(x)
			i += 1
		return [[
 			[_args[0], _args[1]],
 			[_args[0] + _args[2], _args[1]],
 			[_args[0] + _args[2], _args[1] + _args[3]],
			[_args[0], _args[1] + _args[3]],
 			[_args[0], _args[1]]			
 		]];