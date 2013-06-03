import numpy as np
import PySide.QtGui as qtg
import PySide.QtCore as qtc
import Image as Image
import landscape_generator as lg

paddock_color = (20, 200, 20)

class LSMainWindow(qtg.QMainWindow):
	def __init__(self):
		super(LSMainWindow, self).__init__()
		
		self._suitability = 1
		
		self._landscape = None
		self._landscapePMI = None
		
		#Set up the view for the image display
		self._scene = qtg.QGraphicsScene()
		self._sceneView = qtg.QGraphicsView(self._scene, self)
		self._sceneView.setEnabled(False)
		self.setCentralWidget(self._sceneView)
		
		#Initialise variables for selection using the mouse
		self._lastPressLocation = None
		self._currentMouseLocation = None
		self._mousePressed = False
		self._selectionArea = None

		#Add the open file menu option
		self.addFileMenu()
		self.addLandscapeMenu()

		self.loadImage()
		
		self.setWindowTitle("Landscaper")
		self.show()

	
	def addFileMenu(self):
		openImage = qtg.QAction(qtg.QIcon("open.png"), "Open Image", self)
		openImage.setShortcut("Ctrl+O")
		openImage.triggered.connect(self.loadImage)

		saveImage = qtg.QAction(qtg.QIcon("save.png"), "Save Image", self)
		saveImage.setShortcut("Ctrl+S")
		saveImage.triggered.connect(self.saveImage)

		saveArray = qtg.QAction(qtg.QIcon("save.png"), "Save Array", self)
		saveArray.setShortcut("Ctrl+Shift+S")
		saveArray.triggered.connect(self.saveLandscapeArray)

		self.addMenuWithActions("&File", [openImage, saveImage, saveArray])
	
	
	def addLandscapeMenu(self):
		detect = qtg.QAction(qtg.QIcon("open.png"), "Detect Regions", self)
		detect.triggered.connect(self.detectRegions)
		
		suitability = qtg.QAction(qtg.QIcon("open.png"), "Set Suitability", self)
		suitability.triggered.connect(self.setSuitability)
		
		self.addMenuWithActions("&Landscape", [detect, suitability])
	

	def addMenuWithActions(self, menuName, actions):
		menubar = self.menuBar()
		menu = menubar.addMenu(menuName)
		for a in actions:
			menu.addAction(a)


	def showFileDialog(self, title, save=False):
		if save:
			filename, _ = qtg.QFileDialog.getSaveFileName(self, title, ".")
		else:
			filename, _ = qtg.QFileDialog.getOpenFileName(self, title, ".")
		
		return filename


	def loadImage(self):
		#Get the filename
		filename = self.showFileDialog("Open image file")
		
		#Load the new image and add it to the scene
		self.setLandscapeImage(qtg.QPixmap(filename))
	
		
	def saveImage(self):
		filename = self.showFileDialog("Save image file", True)
		self._landscapePMI.pixmap().save(filename, quality=100)
	
	
	def saveLandscapeArray(self):	
		filename = self.showFileDialog("Save image file", True)
		lg.write_landscape(self._landscape, filename)


	def detectRegions(self):
		if self._landscapePMI != None and self._landscape != None:
			#Have the user select a color to use
			#c1 = qtg.QColorDialog.getColor().toTuple()[0:3]
			
			#Convert the landscpae pix map to an image
			image = self._landscapePMI.pixmap().toImage()
			width, height = (image.width(), image.height())
		
			ri, rj = self._lastRightPressLocation
			c1 = qtg.QColor.fromRgb(image.pixel(int(ri), int(rj))).toTuple()[0:3]
			
			#Go through each pixel in the image
			for i in range(width):
				for j in range(height):
					#Get the pixel color value as an rgb tuple
					c2 = qtg.QColor.fromRgb(image.pixel(i, j)).toTuple()[0:3]
					
					#If the pixel color matches the selected color within a set threshold
					if all([v1 - 5 < v2 < v1 + 5 for v1, v2 in zip(c1, c2)]):
						#Set the landscape suitability
						self._landscape[i, j] = self._suitability
			print "found %d pixels of color " % np.sum(self._landscape == self._suitability) + str(c1)
	
	def setSuitability(self):
		self._suitability, _ = qtg.QInputDialog.getDouble(self,
			"set landscape suitablity", "suitability", value=self._suitability, minValue=0., maxValue=1., step=0.001)
	
	
	def setLandscapeImage(self, pixmap):
		#If a previous image exists, remove it from the scene
		if self._landscapePMI != None:
			self._scene.removeItem(self._landscapePMI)
	
		#Add the pixmap to scene and adjust the size of the view
		self._landscapePMI = self._sceneView.scene().addPixmap(pixmap)
		self._sceneView.setSceneRect(self._landscapePMI.pixmap().rect())

		#Setup the a corresponding numpy array to hold the suitability data
		size = pixmap.size()
		self._landscape = np.zeros((size.width(), size.height()))

		
	def mousePressEvent(self, event):
		if event.button() == qtc.Qt.LeftButton:
			self._lastPressLocation = (event.x(), event.y())
		elif event.button() == qtc.Qt.RightButton:
			self._lastRightPressLocation = (event.x(), event.y())
			

	def mouseReleaseEvent(self, event):
		if event.button() == qtc.Qt.LeftButton:
			self._currentMouseLocation = (event.x(), event.y())
			self.cropToSelectionArea()

	
	def mouseMoveEvent(self, event):
		self._currentMouseLocation = (event.x(), event.y())
		self.drawSelectionArea()

			
	def drawSelectionArea(self):
		if self._selectionArea != None:
			self._sceneView.scene().removeItem(self._selectionArea)
		
		if (self._lastPressLocation != None
						and self._currentMouseLocation != None):
			#Calculate the selection area, which must be square
			px, py = self._lastPressLocation
			cx, cy = self._currentMouseLocation
			width = max([cx - px, cy-py])
			self._selectionArea = self._sceneView.scene().addRect(px, py, width, width)
			self._sceneView.setSceneRect(self._sceneView.sceneRect())


	def cropToSelectionArea(self):
		self._sceneView.scene().removeItem(self._selectionArea)

		#Generate the new pixmap and add it to the scene
		rect = self._selectionArea.rect().toRect()
		self.setLandscapeImage(self._landscapePMI.pixmap().copy(rect))

		
		
"""
	def OnLeftClick(self, event):
		pos = event.GetPosition()
		self.suitable_areas.append(pos)
		
		memdc = wx.MemoryDC()
		memdc.SelectObject(self.BufferBmp)
		self.DrawImagePanel(memdc)
		self.Refresh()
	
	
	def OnRightClick(self, event):
		img = self.BufferBmp.ConvertToImage()
		ls = imageToLandscape(wxImageToPIL(img))
		lg.write_landscape(ls, "/Users/davidsavage/Desktop/Merredin.csv")
	
	
	def OnPaint(self, event):
		dc = wx.PaintDC(self)
		dc.BeginDrawing()
		if self.BufferBmp:
			dc.DrawBitmap(self.BufferBmp, 0 ,0, True)
		
		dc.EndDrawing()

	
	def DrawImagePanel(self, dc):
		#try:
		dc.BeginDrawing()
		dc.Clear()
		
		dc.DrawBitmap(self.background_img, 0, 0, True)
		
		#Draw the suitable areas
		dc.SetBrush(wx.Brush(paddock_color, wx.SOLID))
		dc.SetPen(wx.Pen(paddock_color, 1, wx.SOLID))

		for i, j in self.suitable_areas:
			dc.DrawRectangle(i - 3, j - 3, 6, 6)
		
		dc.EndDrawing()
		return True
		#except:
		#	return False

	
	def setBackgroundImage(self, img):
		w, h = img.GetSize()
		img = img.Resize((min(w, h), min(w, h)), (0, 0))
		w, h = img.GetSize()
				
		#self.suitable_areas.append((w/2, h/2))
		self.suitable_areas.append((600, 100))		
		self.background_img = img.ConvertToBitmap()				
		self.SetVirtualSize((w, h))
		self.SetScrollRate(20,20)
		
		self.BufferBmp = wx.EmptyBitmap(w, h)
		memdc = wx.MemoryDC()
		memdc.SelectObject(self.BufferBmp)
		self.DrawImagePanel(memdc)
	
	
	def loadImage(self):
		img_file = "/Users/david/Desktop/Mullewa.png"	
		img = wx.Image(img_file)
		self.setBackgroundImage(img)
		


def imageToLandscape(img):
	nx, ny = img.size
	nx = min(nx, ny)
	ny = min(nx, ny)
	
	ls = np.zeros((nx, ny))
	
	for index in np.ndindex(img.size):
		if index[0] < nx and index[1] < ny:
			p = img.getpixel(index)
			
			if p == paddock_color:
				ls[index] = 0.5
			else:
				ls[index] = 0
	
	return ls


def wxImageToPIL(wximg):
	w, h = wximg.GetSize()
	data = wximg.GetData()
	
	redImg = Image.new("L", (w, h))
	redImg.fromstring(data[0::3])
	greenImg = Image.new("L", (w, h))
	greenImg.fromstring(data[1::3])	
	blueImg = Image.new("L", (w, h))
	blueImg.fromstring(data[2::3])
	
	pilImg = Image.merge("RGB", (redImg, greenImg, blueImg))
	return pilImg
"""
if __name__ == "__main__":
	#Create the wx app instance
	app = qtg.QApplication([])
	window = LSMainWindow()
	
	app.exec_()
