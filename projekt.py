import tkinter
from tkinter import messagebox
from tkinter import *
import math


def multiply(m, array):

	for i in array:
		i = int(i)

	for i in range(0, len(array), 2):
		result_x = array[i] * m[0][0] + array[i + 1] * m[0][1]
		result_y = array[i] * m[1][0] + array[i + 1] * m[1][1]
		array[i] = result_x
		array[i + 1] = result_y
	return array

class PolyObject:
	def __init__(self, canvas, x, y, tag):
		self.canvas = canvas
		self.x = x
		self.y = y
		self.tag = tag
		self.is_active = 0

	def draw(self, coords, color):
		""" Coords - budem posielat suradnice bodov na vykreslenie poly-lineu
		"""
		print(type(coords))
		if self.is_active:
			self.id = self.canvas.create_polygon(*coords, fill=color, outline='red', width=3, tag=self.tag)
		else:
			self.id = self.canvas.create_polygon(*coords, fill=color, outline='black', width=1, tag=self.tag)

	def set_active(self, state):
		""" Zmena outline farby
		"""
		if state == 1:
			self.canvas.itemconfig(self.id, outline='red', width=3)
			self.is_active = 1
		else:
			self.canvas.itemconfig(self.id, outline='black', width=1)
			self.is_active = 0

	def set_coords(self, coords):
		self.coords = coords

	def get_centroid(self):
		self.coords.append(self.coords[0])
		self.coords.append(self.coords[1])
		A = 0
		for i in range(0, len(self.coords) - 2, 2):
			A += self.coords[i]*self.coords[i+3] - self.coords[i+2]*self.coords[i+1]
		A = A*(1/2)

		Cx = Cy = 0
		for i in range(0, len(self.coords) - 2, 2):
			Cx += (self.coords[i]+self.coords[i+2])*(self.coords[i]*self.coords[i+3] - self.coords[i+2]*self.coords[i+1])
			Cy += (self.coords[i+1]+self.coords[i+3])*(self.coords[i]*self.coords[i+3] - self.coords[i+2]*self.coords[i+1])

		self.coords.pop()
		self.coords.pop()
		x = Cx/(6*A)
		y = Cy/(6*A)

		return [x, y]

	def rotate(self):
		matrix = [[math.cos(math.radians(15)), math.sin(math.radians(15))],
				  [-math.sin(math.radians(15)), math.cos(math.radians(15))]]

		center = self.get_centroid()

		## translate polygon to the center of coordinate system
		for i in range(0, len(self.coords), 2):
			self.coords[i] -= center[0]
			self.coords[i+1] -= center[1]

		## rotate polygon
		self.coords = multiply(matrix, self.coords)

		## translate rotated polygon back to the original position
		for i in range(0, len(self.coords), 2):
			self.coords[i] += center[0]
			self.coords[i+1] += center[1]



	def move(self, xx, yy):
		# matrix = [[1, 0, 0],
		# 		  [0, 1, 0],
		# 		  [xx - self.x, yy - self.y, 1]]
		# self.x = xx - self.x
		# self.y = yy - self.y
		# self.coords = multiply(matrix, self.coords)
		self.canvas.move(self.id, xx, yy)


class Playground(Tk):
	def __init__(self):
		super().__init__()
		main = Frame(self, bg='red', width=1000)
		main.pack()

		self.labelframe = LabelFrame(main, text="Rotation", width=200)
		self.labelframe.pack(side=LEFT, fill=BOTH)

		rotate = Button(self.labelframe, text = "Rotate", width = 13, command = self.rotate_object)
		rotate.pack(side = TOP) 

		self.playground = Canvas(main, width=800, height=600, highlightthickness=0)
		self.playground.pack(side=LEFT, fill=BOTH)
		self.update()

		self.polygons_array = []
		self.obj_id = None
		self.load_from_file('objects.txt')
		self.playground.bind("<ButtonPress-1>", self.click)
		self.playground.bind("<B1-Motion>", self.drag)
		self.playground.bind("<ButtonRelease-1>", self.drop)

	def rotate_object(self):

		if self.p1.is_active:
			# print(self.playground.coords(self.p1.id))
			self.obj_id = self.p1.id
			self.p1.rotate()
		elif self.p2.is_active:
			# print(self.playground.coords(self.p2.id))
			self.obj_id = self.p2.id
			# self.p2.set_coords(rotated_points)
			self.p2.rotate()
		else:
			messagebox.showwarning("Warning","No object is selected!!!")
			return

		self.redraw_canvas()
		
	def load_from_file(self, file):
		""" Naciatnie suradnic zo suboru
		"""
		coords = []
		f = open(file, 'r')

		for line in f:
			coords.append([line.strip(' ')])
		f.close()
		self.create_polygons(coords)

	def create_polygons(self, polygons):
		""" Z nacitanych suradnic vytvor objekty polygonov a vykresli ich
		"""
		coords_1 = polygons[0]
		coords_2 = polygons[1]
		self.p1 = PolyObject(self.playground, 100, 100, 'object1')
		self.p1.draw(coords_1, 'white')
		# print('Suradnice bodov objektu p1: ', self.playground.coords(self.p1.id))
		self.p2 = PolyObject(self.playground, 400, 400, 'object2')
		self.p2.draw(coords_2, 'navy')
		self.polygons_array.append(self.p1)
		self.polygons_array.append(self.p2)
		# print(self.polygons_array)

		self.p1.set_coords(self.playground.coords(self.p1.id))
		self.p2.set_coords(self.playground.coords(self.p2.id))

	def redraw_canvas(self):
		self.playground.delete('all')
		self.p1.draw(self.p1.coords, 'white')
		self.p2.draw(self.p2.coords, 'navy')


	def click(self, event):
		self.obj_id = event.widget.find_closest(event.x, event.y)[0]
		if self.playground.gettags(event.widget.find_closest(event.x, event.y)):
			self.obj_tag = self.playground.gettags(event.widget.find_closest(event.x, event.y))[0]
		if 'object' in self.obj_tag:
			self.playground.lift(self.obj_id)
		if 'object1' in self.obj_tag:
			self.p1.set_active(1)
			self.p2.set_active(2)	
		elif 'object2' in self.obj_tag:
			self.p2.set_active(1)
			self.p1.set_active(2)

		# self.initial_coords = self.playground.coords(self.obj_id)
		self.ex, self.ey = event.x, event.y

	def drop(self, event):
		if self.p1.is_active:
			self.p1.set_coords(self.playground.coords(self.p1.id))
		elif self.p2.is_active:
			self.p2.set_coords(self.playground.coords(self.p2.id))
		self.p1.get_centroid()
		self.p2.get_centroid()

	def drag(self, event):
		# if self.obj_tag == ('object' + str(self.obj_id)):
			self.playground.move(self.obj_id, event.x - self.ex, event.y - self.ey)
			self.ex, self.ey = event.x, event.y

	def v_clip(self):
		pass

	def clip_vertex(self):
		pass

	def clip_edge(self):
		pass


if __name__ == '__main__':
	mw = Playground()
	mw.mainloop()
