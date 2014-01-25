# Convert from PIL Image and Numpy array to PyPNG format
# By Tim Sheerman-Chase, 2013
# This code is released as CC0 and public domain code. Do what you will!

import png
from PIL import Image
import numpy as np

class NumpyArrayToPyPngAdapter:
	def __init__(self, img):
		self.img = img

	def __len__(self):
		return len(self.img)

	def __getitem__(self, key):
		row = self.img[key]
		out = []
		for p in row:
			if hasattr(p, '__iter__'):
				# Multi-channel image
				out.extend(p)
			else:
				# Single channel image
				out.append(p)
		return out

class PilImageToPyPngAdapter:
	def __init__(self, im):
		self.im = im
		self.iml = im.load()

	def __len__(self):
		return self.im.size[1]

	def __getitem__(self, row):
		out = []
		for col in range(self.im.size[0]):
			px = self.iml[col, row]
			if hasattr(px, '__iter__'):
				# Multi-channel image
				out.extend(px)
			else:
				# Single channel image
				out.append(px)
		return out

def PyPngToPilImage(width, height, rows, args):
	m = 'RGB'
	planes = int(args['planes'])
	if planes == 1: m = 'L'
	if planes == 4: m = 'RGBA'
	im = Image.new(m, (width, height))
	iml = im.load()
	for rowNum, row in enumerate(rows):
		for colNum in range(width):
			if planes > 1:
				iml[colNum, rowNum] = tuple(row[colNum * planes:(colNum + 1) * planes])
			else:
				iml[colNum, rowNum] = row[colNum]
	return im

def savePNG(fname, ima):
	
	fi = open(fname, "wb")
	wri = png.Writer(height=ima.shape[0], width=ima.shape[1], greyscale=(ima.ndim == 2))
	wri.write(fi, NumpyArrayToPyPngAdapter(ima))
	fi.flush()
	fi.close()


if __name__ == "__main__":
	im = Image.open("2.jpg")
	# im = im.convert("RGB") # Single channel "L" should also work
	ima = np.array(im)
	savePNG("test.png", ima)
# 	fi = open("test.png","wb")
# 	wri = png.Writer(size=im.size, greyscale=(len(im.mode)==1))
# 	wri.write(fi, NumpyArrayToPyPngAdapter(ima))
# 	#wri.write(fi, PilImageToPyPngAdapter(im))
# 	fi.flush()
# 	fi.close()
	
# 	read = png.Reader(file=open("test.png","rb"))
# 	readRet = read.read()
# 	print readRet
# 	pilim = PyPngToPilImage(*readRet)
# 	pilim.show()


