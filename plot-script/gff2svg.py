#!/usr/bin/env python
#Yuzhen Ye, May 22, 2018
#see example usage in gff2svg-test.php
import sys
import os

typecol = {}
def svg_header(height=200,width=640,divid="holder"):
	s = '<div id="' + divid + '">\n'
   	s += '<svg height="%d" width="%d" version="1.1" xmlns="http://www.w3.org/2000/svg" style="overflow: hidden; position: relative; ">' % (height, width)
	s += ' <defs> <marker id="mCircle" markerWidth="40" markerHeight="40" refx="20" refy="20"> <circle cx="20" cy="20" r="12" style="stroke: #000000; fill:#ffffff;"/> </marker>'
	s += '    <marker id="mArrow" markerWidth="13" markerHeight="13" refx="2" refy="6" orient="auto"> <path d="M2,2 L2,11 L10,6 L2,2" style="fill: #000000;" />'
	s += '    </marker>'
	s += '   </defs>'
	return s

def svg_end():
	s = 'Sorry, your browser does not support inline SVG.\n'
	s += '</svg>'
	return s

def text_str(scale=1,x0=30,y0=0,x=0,y=0,text="test", anchor="middle", col="black", display="block", fontsize=15):
	s = '   <text text-anchor="%s" x="%d" y="%d" fill="%s" font-size="%d" style="display:%s;">%s</text>' % (anchor, int(x0+x*scale),y0+y,col,fontsize,display,text)
	return s

def line_str(x0=0,y0=0,seqlen=400,scale=1,col="black"):
	s = '   <path d="M%d,%d L%d,%d" style="stroke: %s; stroke-width: 2px; fill: none;"></path>' % (x0,y0,x0+seqlen*scale,y0, col)
	return s

def get_col(what="gene", type="unk"):
	global typecol
	if what[:10] == 'antiRepeat':
		return "red"	
	elif what[:7] == 'repeat:':
		return "white"

	#for idx in range(len(typedes)): #same type, different subtypes, use the same color
	#	if typedes[idx] in type: 
	#		return typecol[idx]
	if type in typecol:
		return typecol[type]
	else: 
		return "white"

def poly_str(ele=1,scale=1,x0=30,y0=0,yh=10,dy=20,arrow=10,beg=100,end=200,offset=0,dir='+',what="gene",type="type-II",strokewidth=1,strokecol="black",fillcol="white",des=[],tag="",divid="div",loc="down", showtxt=""):
	bx = beg - offset
	ex = end - offset
	fillcol = get_col(what, type)
	if arrow > (ex-bx)*scale*0.5:
		arrow = (ex-bx)*scale*0.5
			
	s = ''
	s += '   <polygon id="%s_ele_%d" style="stroke:%s; stroke-width:%d; fill:%s; " locus_tag="%s" points="' % (divid, ele, strokecol, strokewidth, fillcol, des)
	if what == "contig":
		s += "%d,%d %d,%d %d,%d %d,%d" % (bx*scale+x0,y0-yh,bx*scale+x0,y0+yh,ex*scale+x0,y0+yh,ex*scale+x0,y0-yh)
	elif what[:7] == "repeat:":
		#warning: ex, bx artificially changed for cosmetiic purpose
		if (ex - bx) * scale < 10:
			bx = bx - 5.0 / scale	
			ex = ex + 5.0 / scale
		s += "%d,%d %d,%d %d,%d %d,%d %d,%d %d,%d" % (bx*scale+x0,y0,bx*scale+x0+arrow,y0+yh,ex*scale+x0-arrow,y0+yh,ex*scale+x0,y0,ex*scale+x0-arrow,y0-yh,bx*scale+x0+arrow,y0-yh)
	elif what[:10] == "antiRepeat":
		x = (bx + 0.5*(ex-bx))*scale+x0
		s += "%d,%d %d,%d %d,%d" % (x,y0,x-4,y0-yh,x+4,y0-yh)
	else: #gene
		if dir == '+':
			s += "%d,%d %d,%d %d,%d %d,%d %d,%d" % (bx*scale+x0,y0-yh,bx*scale+x0,y0+yh,ex*scale+x0-arrow,y0+yh,ex*scale+x0,y0,ex*scale+x0-arrow,y0-yh)
		else:
			s += "%d,%d %d,%d %d,%d %d,%d %d,%d" % (bx*scale+x0,y0,bx*scale+x0+arrow,y0+yh,ex*scale+x0,y0+yh,ex*scale+x0,y0-yh,bx*scale+x0+arrow,y0-yh)
	s += '"></polygon>\n'
	text="%s %s\nLocation %d-%d bp" % (what, type, beg, end)
	s1 = "" #hidden elements to be drawn at the end
	if (what[:7] == 'repeat:') or (what[:10] == 'antiRepeat'):
		textheight = 2 * 20
	else:
		textheight = 3 * 25
	textheight += len(des) * 25
	if loc == 'up':
		yhigh, ylow, tloc = y0 - yh - textheight, y0 - yh - 3, y0 + yh - textheight 
	else:
		yhigh, ylow, tloc = y0 + yh + 3, y0 + yh + textheight, y0 + 3*yh
	if loc == 'up':
		textheight = -textheight
	s1 += '   <polygon visibility="hidden" style="stroke:%s; stroke-width:%d; fill:%s; " points="%d,%d %d,%d %d,%d %d,%d">\n' % ("white", strokewidth, "#CCCCCC", bx*scale+x0-10,ylow,bx*scale+x0-10,yhigh,bx*scale+x0+300,yhigh,bx*scale+x0+300,ylow)
	s1 += '     <set attributeName="visibility" from="hidden" to="visible" begin="%s_ele_%d.mouseover" end="%s_ele_%d.mouseout" />\n' % (divid, ele, divid, ele)
	s1 += '   </polygon>'
	s1 += '   <text id="des_%d" text-anchor="%s" x="%d" y="%d" fill="%s" font-size="16" visibility="hidden">\n' % (ele, "start", x0+bx*scale,tloc,"black")
	if what[:7] == "repeat:":
		s1 += '     <tspan x="%d">%s</tspan>\n' % (x0+bx*scale, "CRISPR array") 
		s1 += '     <tspan x="%d" dy="%d">Location: %d-%d bp</tspan>\n' % (x0+bx*scale, dy, beg, end) 
	elif what[:10] == "antiRepeat":
		s1 += '     <tspan x="%d">%s</tspan>\n' % (x0+bx*scale, what) 
		s1 += '     <tspan x="%d" dy="%d">Location: %d-%d bp</tspan>\n' % (x0+bx*scale, dy, beg, end) 
	else:
		s1 += '     <tspan x="%d">Gene: %s</tspan>\n' % (x0+bx*scale, what) 
		s1 += '     <tspan x="%d" dy="%d">Type: %s</tspan>\n' % (x0+bx*scale, dy, type) 
		s1 += '     <tspan x="%d" dy="%d">Location: %d-%d bp</tspan>\n' % (x0+bx*scale, dy, beg, end) 
	for ades in des:
		s1 += '     <tspan x="%d" dy="%d">%s</tspan>\n' % (x0+bx*scale, dy, ades)

	s1 += '     <set attributeName="visibility" from="hidden" to="visible" begin="%s_ele_%d.mouseover" end="%s_ele_%d.mouseout" />' % (divid, ele, divid, ele)
	s1 += '   </text>'
	s1c = s1

	if showtxt:
		s += '   <text id="%s_ele_%d" text-anchor="%s" x="%d" y="%d" fill="%s" font-size="%d" style="display:%s;">%s</text>' % (divid + "t", ele, "middle", (bx + (ex-bx)*0.5)*scale+x0,y0+5,"blue",13,"block",showtxt)
		s1c = s1.replace(divid, divid + "t")
		s1 += s1c

	return s, s1

def oval_str(x0=0, y0=0, scale=1, x=100, y=100, text="", rx=9, ry=9,fillcol="yellow", ele=1, divid="div",beg=0, end=0, up=True):
	if rx < 2:
		rx = 2
	s = '   <ellipse id="%s_ele_%d" cx="%d" cy="%d" rx="%d" ry="%d", stroke="black" stroke-width="1" fill="%s"></ellipse>' % (divid, ele, x0+x*scale, y0+y, rx, ry, fillcol)
	if text:
		if up:
			s += "\n" + text_str(x0=x0, y0=y0, x=x, y=y-1.2*ry, scale=scale,text=text, anchor="middle",col="black")
		else:
			s += "\n" + text_str(x0=x0, y0=y0, x=x, y=y+2.4*ry, scale=scale,text=text, anchor="middle",col="black")
		yhigh, ylow = y0+y+ry, y0+y+4*ry
		xleft, xright = int(x0+x*scale)-120, int(x0+x*scale) + 120
		if xleft < 0:
			xright -= xleft
			xleft = 0
		s += '   <polygon visibility="hidden" style="stroke:%s; stroke-width:%d; fill:%s; " points="%d,%d %d,%d %d,%d %d,%d">\n' % ("white", 1, "#CCCCCC", xleft,yhigh,xleft,ylow,xright,ylow,xright,yhigh)
		s += '     <set attributeName="visibility" from="hidden" to="visible" begin="%s_ele_%d.mouseover" end="%s_ele_%d.mouseout" />\n' % (divid, ele, divid, ele)
		s += '   </polygon>'
		s += '   <text id="des_%d" text-anchor="%s" x="%d" y="%d" fill="%s" font-size="16" visibility="hidden">%s\n' % (ele, "middle", int(x0+x*scale),y0+y+3*ry,"black","Locus " + text + ": " + str(beg) + "-" + str(end) + " bp")
		s += '      <set attributeName="visibility" from="hidden" to="visible" begin="%s_ele_%d.mouseover" end="%s_ele_%d.mouseout" />' % (divid, ele, divid, ele)
		s += '   </text>'
	return s

def main():
	if len(sys.argv) < 5:
		sys.exit("Need input: gfffile id divid colorfile")

	gfffile = sys.argv[1]
	seqid = sys.argv[2]
	divid = sys.argv[3]
	colorfile = sys.argv[4]

	global typecol;
	inf = open(colorfile, "r")
	for aline in inf:
		subs = aline.strip().split()
		typecol[subs[0]] = subs[2]
	inf.close()

	#get array/cas annonations from gff file
	inf = open(gfffile, "r")
	features = []
	for aline in inf:
		if aline[0] == '#':
			continue
		else:
			subs = aline.strip().split("\t")	
			if len(subs) != 9:
				sys.exit("wrong gff input")
			if subs[0] != seqid:
				continue
			if subs[2] == 'region':
				seqlen = int(subs[4]) - int(subs[3]) + 1
				continue
			notes = subs[-1].split(";")
			des, detail = "unk", "unk"
			name = ""
			for anote in notes:
				if "des=" in anote:
					des = anote[4:]
				if "note=" in anote:
					detail = anote[5:]
				if "Name=" in anote:
					name = anote[5:]
			if name != "":
				des = name + ":" + des
			features.append([subs[0], int(subs[3]), int(subs[4]), subs[6], des, detail])
	inf.close()

	#define blocks
	ifarray, ifcas = False, False
	lastbeg, lastend = 0, 0
	blockidx = []
	for idx in range(len(features)):
		beg, end = features[idx][1:3]	
		#print "idx", idx, "beg", beg, "end", end
		if idx and beg - lastend >= 10000:
			blockidx.append(idx-1)
		lastbeg, lastend = beg, end
	blockidx.append(len(features) - 1)
	#print "blockidx", blockidx
	blocksize = []
	for idx in range(0, len(blockidx)):
		if idx:
			bf = blockidx[idx - 1] + 1
		else:
			bf = 0
		ef = blockidx[idx]
		blocksize.append(features[ef][2] - features[bf][1])
	maxblocksize = max(blocksize)

	#print svg
	if (len(sys.argv) > 5) and (sys.argv[5] == "small"):
		width0 = 600
	else:
		width0 = 800
	x0, y0, width, height = 80, 25, width0, 70 + 70 * len(blockidx)
	print svg_header(height=height+y0, width=width+x0+300, divid=divid)

	#print overview
	print text_str(x0=x0, y0=y0, text="Overview", anchor="start")
	y0 += 25
	scale = width * 1.0 / seqlen
	##overview line
	bf = 0
	print line_str(x0=x0,y0=y0,seqlen=seqlen,scale=scale)
	##blocks on the overview line
	ele = 0
	up = True
	for idx in range(len(blockidx)):
		if idx:
			bf = blockidx[idx - 1] + 1
		ef = blockidx[idx] 
		if idx and (features[bf][1] - features[bf-1][2]) < seqlen * 0.02: #blocks close to each other
			up = not up
		r = 0.5 * (features[ef][2]-features[bf][1])
		ele += 1
		#print "block", idx, "begin-feature", bf, "end-feature", ef, "r", r, "scale", scale
		#raw_input("type enter to continue..")
		print oval_str(ele=ele, x0=x0,y0=y0,x=features[bf][1] + r,y=0,rx=r*scale,ry=9,scale=scale,text=str(idx+1), divid=divid, beg=features[bf][1], end=features[ef][2],up=up)

	#zoom in on individual blocks
	bf = 0
	scale = width * 1.0 / maxblocksize
	#print "scale for blocks", scale, "blocksize", maxblocksize, "seqlen", seqlen
	sidx = 1
	s_hidden = ""
	for idx in range(len(blockidx)):
		y0 += 50
		if idx:
			bf = blockidx[idx - 1] + 1
		##line for the block
		print text_str(x0=x0, y0=y0, text="Zoom in on locus " + str(idx + 1), anchor="start")
		y0 += 20
		print line_str(x0=x0,y0=y0,seqlen=blocksize[idx]*1.01,scale=scale)
		##arrows for the genes
		if idx == len(blockidx) - 1: 
			loc = "up"
		else:
			loc = "down"
		for f in range(bf, blockidx[idx] + 1):
			des = features[f][4]
			offset = features[bf][1]
			ele += 1
			if des[:7] == 'repeat:': 
				subs = features[f][-1].split("-r")
				copy = int(subs[-1][5:])
				deslist = ["Num of repeats: " + str(copy),  "Spacers: s" + str(sidx) + "-s" + str(sidx + copy - 1 - 1)]  
				sidx += copy - 1
				s, s1 = poly_str(ele=ele, x0=x0,y0=y0,what=features[f][4],type=features[f][-1],beg=features[f][1],end=features[f][2],offset=offset,scale=scale,des=deslist,divid=divid,loc=loc,showtxt="x" + str(copy))
				s_hidden += s1
				print s
			elif des[:10] == 'antiRepeat':
				subs = features[f][-1].split("-")
				deslist = ["Strand: " + features[f][3], "Match length: " + subs[0] + " bp",  "Mismatch: " + subs[1] + " bp"]
				s, s1 = poly_str(ele=ele, x0=x0,y0=y0,what=features[f][4][:10],type=features[f][-1],beg=features[f][1],end=features[f][2],offset=offset,scale=scale,des=deslist,divid=divid,fillcol="red",loc=loc)
				s_hidden += s1
				print s
			else:
				s, s1 = poly_str(ele=ele, x0=x0,y0=y0,dir=features[f][3], what=features[f][4],type=features[f][-1],beg=features[f][1],end=features[f][2],offset=offset,scale=scale,des=[],divid=divid,loc=loc)
				s_hidden += s1
				print s

	print s_hidden
	print svg_end()
	print "</div>"

main()
