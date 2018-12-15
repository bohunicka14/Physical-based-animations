from tkinter import *


class PolyObject:
    def __init__(self, canvas, x, y, tag):
        self.canvas = canvas
        self.x = x
        self.y = y
        self.tag = tag

    def draw(self, coords, color):
        """ Coords - budem posielat suradnice bodov na vykreslenie poly-lineu
        """
        print(type(coords))
        self.id = self.canvas.create_polygon(*coords, fill=color, outline='black', tag=self.tag)

    def is_active(self, state):
        """ Zmena outline farby
        """
        if state == 1:
            self.canvas.itemconfig(self.id, outline='red', width=3)
        else:
            self.canvas.itemconfig(self.id, outline='black', width=1)

    def move(self, xx, yy):
        self.canvas.move(self.id, xx, yy)


class Playground(Tk):
    def __init__(self):
        super().__init__()
        parent = Frame(self, bg='red')
        parent.pack(side=BOTTOM, fill=BOTH, expand=TRUE)
        self.playground = Canvas(parent, width=800, height=600, highlightthickness=0)
        self.playground.pack()
        self.update()

        self.polygons_array = []
        self.load_from_file('objects.txt')
        self.playground.bind("<ButtonPress-1>", self.click)
        self.playground.bind("<B1-Motion>", self.drag)

    def load_from_file(self, file):
        """ Naciatnie suradnic zo suboru
        """
        coords = []
        f = open(file, 'r')

        for line in f:
            coords.append([line.strip(' ')])
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

    def click(self, event):
        self.obj_id = event.widget.find_closest(event.x, event.y)[0]
        # print(self.obj_id)
        if self.playground.gettags(event.widget.find_closest(event.x, event.y)):
            self.obj_tag = self.playground.gettags(event.widget.find_closest(event.x, event.y))[0]
        if 'object' in self.obj_tag:
            self.playground.lift(self.obj_id)
        if 'object1' in self.obj_tag:
            self.p1.is_active(1)
            self.p2.is_active(2)
        if 'object2' in self.obj_tag:
            self.p2.is_active(1)
            self.p1.is_active(2)

        self.initial_coords = self.playground.coords(self.obj_id)
        self.ex, self.ey = event.x, event.y

    def drag(self, event):
        if self.obj_tag == ('object' + str(self.obj_id)):
            self.playground.move(self.obj_tag, event.x - self.ex, event.y - self.ey)
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
