from ctypes import *

class Complex_c(Structure):
	_fields_=[
		('x',c_double),
		('y',c_double)
		]
	def __init__(self,n=0+0j):
		self.x=c_double(n.real);
		self.y=c_double(n.imag);
def toComplex_c(vals):
	cs=(Complex_c*len(vals))();
	for i in range(len(vals)):
		cs[i]=Complex_c(vals[i])
		#arr[i].x=(vals[i].real);
		#arr[i].y=(vals[i].imag);
	return cs;
def toComplex(cs):
	vals=[]
	for v in cs:
		vals.append(complex(v.x,v.y));
	return vals;
import os

libdir=os.path.split(os.path.realpath(__file__))[0]; 
cfftlib=CDLL(libdir+"/cFFT.dylib")
cfftlib.FFT.argtypes=(POINTER(Complex_c),POINTER(Complex_c),c_int)
cfftlib.IFFT.argtypes=(POINTER(Complex_c),POINTER(Complex_c),c_int)

def complexs(len):
	return (Complex_c*len)();
def fft(f,power):
	F=(Complex_c*(1<<power))();
	cfftlib.FFT(f,F,power);
	return F;
def ifft(F,power):
	f=(Complex_c*(1<<power))();
	cfftlib.IFFT(F,f,power);
	return f;
if __name__=='__main__':
	f=complexs(1024);
	fft(f,10);
	import time;
	for i in range(1024*1024):
		f=complexs(1024)
