import tkinter
from tkinter import messagebox
from tkinter import *
import math


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

	def move(self, xx, yy):
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
		self.load_from_file('objects.txt')
		self.playground.bind("<ButtonPress-1>", self.click)
		self.playground.bind("<B1-Motion>", self.drag)

	def multiply(self, m, array):
		'''funkcia na nasobenie bodu a matice, vracia pole prepocitanych bodov
		'''
		for i in array:
			i = int(i)

		for i in range(0, len(array), 2):
			result_x = array[i]*m[0][0] + array[i+1]*m[0][1]
			result_y = array[i]*m[1][0] + array[i+1]*m[1][1]
			array[i] = result_x
			array[i+1] = result_y
		return array   

	def rotate_object(self):
		## pocitanie rotacie po osi pomocou matice, pocitanie sa realizuje na povodnych suradniciach
		## vykreslenie na plochu
		matrix = [[math.cos(math.radians(15)),math.sin(math.radians(15))],
				  [-math.sin(math.radians(15)),math.cos(math.radians(15))]]
		
		if self.p1.is_active:
			print(self.playground.coords(self.p1.id))
			rotated_points = self.multiply(matrix, self.playground.coords(self.p1.id))
			self.p1.set_coords(rotated_points)
		elif self.p2.is_active:
			print(self.playground.coords(self.p2.id))
			rotated_points = self.multiply(matrix, self.playground.coords(self.p2.id))
			self.p2.set_coords(rotated_points)
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
		print('Suradnice bodov objektu p1: ', self.playground.coords(self.p1.id))
		self.p2 = PolyObject(self.playground, 300, 300, 'object2')
		self.p2.draw(coords_2, 'navy')
		self.polygons_array.append(self.p1)
		self.polygons_array.append(self.p2)
		print(self.polygons_array)

		self.p1.set_coords(self.playground.coords(self.p1.id))
		self.p2.set_coords(self.playground.coords(self.p2.id))

	def redraw_canvas(self):
		self.playground.delete('all')
		self.p1.draw(self.p1.coords, 'white')
		self.p2.draw(self.p2.coords, 'navy')


	def click(self, event):
		self.obj_id = event.widget.find_closest(event.x, event.y)[0]
		print(self.obj_id)
		if self.playground.gettags(event.widget.find_closest(event.x, event.y)):
			self.obj_tag = self.playground.gettags(event.widget.find_closest(event.x, event.y))[0]
		if 'object' in self.obj_tag:
			self.playground.lift(self.obj_id)
		if 'object1' in self.obj_tag:
			self.p1.set_active(1)
			self.p2.set_active(2)
		if 'object2' in self.obj_tag:
			self.p2.set_active(1)
			self.p1.set_active(2)

		self.initial_coords = self.playground.coords(self.obj_id)
		self.ex, self.ey = event.x, event.y

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
