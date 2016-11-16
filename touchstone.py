#!/usr/bin/env python

import numpy as np
from numpy import complex128 as complex

import math
import cmath
import time
import re

def readmodel(sNpfilename, **kwargs):
    """
    Read a touchstone file and output a ScatteringModel.
    
    See the help for ScatteringModel for the kwargs.
    """
    f, S = read(sNpfilename)
    kwargs['name'] = sNpfilename
    return ScatteringModel(f,S,**kwargs)

class ScatteringModel(object):
    """
    The scattering model is a class to represent an N-port S-parameter.
    
    Contains a numpy.ndarray with dimensions F x N x N, where F is the
    number of frequency points. 
    
    Options:
    oe_ordering = True (default)
    mixedmode = False (default)
    """
    def __init__(self, farray = np.zeros([1]), Sarray = np.zeros([1,1,1]), **kwargs):
        self.farray = farray # frequency array
        self.Sarray = Sarray # ScatteringMatrix array
        
        #self.attributes = {'oe_ordering': True,
        #                   'mixedmode': False}
        self.oe_ordering = True
        self.mixedmode = False
        
        # Modify default parameters.  TODO: Allow the file input to do this in readmodel.
        for ikey, ival in kwargs.items():        
            #self.attributes[ikey] = ival
            self.__setattr__(ikey, ival)            

        # When translating the S-parameter data from an emout file from an AFS
        # sweep, the frequencies are out of order.  So, reordering is required
        # to make sense of th data.
        
    def __copy__(self):
        M = ScatteringModel()
        for ikey, ival in self.__dict__.items():
            M.__setattr__(ikey,ival)
        return M

    def getfrequency(self, slice_sel=None):
        if slice_sel is None:
            return self.farray
        else:
            print( slice_sel )
            return self.farray[slice_sel]

    def __getitem__(self, slice_sel):
        """Return the ScatteringMatrix for the index slice specified."""
        return self.Sarray[slice_sel]
    
    def __len__(self):
        """Return the number of frequencies."""
        return len(self.farray)

    def getTrace(self, port_neg, port_pos):
        if self.mixedmode:
            def label2portindex(s):                
                N = self.getPorts()
                mode = s[0].lower()
                mixedmodeport = int(s[1:])
                ind = mixedmodeport-1
                if mode == 'c':
                    ind += N/2 # common mode ports are listed after differential                    
                return ind
                    
            ind0 = label2portindex(port_neg)
            ind1 = label2portindex(port_pos)
            return self.Sarray[:,ind0,ind1]
        else:
            # Expect an integer or a string with regular (not 0-based index) for
            # the input or V+ (port_pos) and the output or V- port (port_neg).            
            ind0 = int(port_neg) - 1
            ind1 = int(port_pos) - 1
            return self.Sarray[:,ind0,ind1]

    def concat(self, S2):
        """
        Concatenate the data of this model and the S2 scattering model.
        
        The only data to really connect is the frequency and S array.  At this point I have done nothing to interleave the frequency bands or deal with any overlap in frequency converage.
        """
        if not isinstance(S2, ScatteringModel):
            raise TypeError
        
        f = np.hstack((self.farray, S2.farray))
        
        if self.mixedmode and not S2.mixedmode:
            S2.toMixedMode()
            S = np.vstack((self.Sarray, S2.Sarray))
        elif not self.mixedmode and S2.mixedmode:
            S2.toSingleEnded()
            S = np.vstack((self.Sarray, S2.Sarray))
        else:
            S = np.vstack((self.Sarray, S2.Sarray))        

        M = self.__copy__()
        M.farray = f
        M.Sarray = S
        return M        
    
    def getPorts(self):
        return np.size(self.Sarray, 1) # Assume that it is square.

    def reorderPorts(self, reordering_l):
        """
        Renumber the ports according to the provided list.  

        reordering_l - lists the old port numbers in the positions of the new
        port numbers. I.e., reordering_l[0] = 5 would have old Port 5 be the new
        Port 1.
        """
        self.Sarray = PortReorder(reordering_l, self.Sarray)

    
    def toMixedMode(self):
        """
        Change to mixed mode if self.mixedmode = False. 
        Do nothing if it is True.
        """
        if not self.mixedmode:
            mode = {True: 1, False: 0}[self.oe_ordering]
            self.Sarray = MixedMode(self.Sarray, mode)
            self.mixedmode = True
            
        return self       
        
    def toSingleEnded(self):
        """
        Change to single-ended mode if self.mixedmode = True. 
        Do nothing if it is False.
        """
        if self.mixedmode:
            mode = {True: 1, False: 0}[self.oe_ordering]
            self.Sarray = MixedMode(self.Sarray, mode, convert='To Single-Ended Mode')
            self.mixedmode = False
            
        return self
            
    def Renormalize(self, Zo_vector = None):
        """
        If Renormalize is run without a Zo_vector given, the S-parameters will
        be reset back to the original 50 ohm (or 100 ohm / 25 ohm for mixed
        modes).
        """
        if Zo_vector and hasattr(self, 'Soriginal'):
            self.Sarray = self.Soriginal[:]            
        else:
            Zo_base = [50,]*self.getPorts()     
            self.Soriginal = self.Sarray[:]
            snew = renormalize(self.Sarray, Zo_base, Zo_vector)
            self.Sarray = snew
        return self.Sarray
        
    def TransmissionMatrix(self):
        """
        Output the T-parameter equivalent of the scattering matrix
        if the port numbers are even.  Assumes an odd and even
        port ordering.
        """
        pass
    
    def T(self):
        """
        Short form of TransmissionMatrix.
        """
        return self.TransmissionMatrix()
    
    def getZ(self, zo=None):
        """
        Translate the model to a Z matrix representation.
        
        The function accepts Zo as an optional argument.  This will override the
        models previous Zo if assigned.  If neither the argument Zo or the 
        model Zo is specified, 50 ohms is assumed.
        """
        
        # Check if zo was specified and if it already exists as an attribute.
        if zo is None:
            if hasattr(self, 'Zo'):
                zo = self.Zo
            else:
                zo = 50.             
        
        mattype = self.Sarray.dtype
        Z = np.zeros(self.Sarray.shape, dtype=mattype)
        for k in range(Z.shape[0]):
            s = self.Sarray[k,:,:]
            u = np.eye(s.shape[0], dtype=mattype )
            Z[k] = zo * np.dot( np.linalg.inv(u - s), (u + s) )           
        
        return Z

    def getY(self, zo=None):
        """
        Translate the model to a Y matrix representation.
        
        The function accepts Zo as an optional argument.  This will override the
        models previous Zo if assigned.  If neither the argument Zo or the 
        model Zo is specified, 50 ohms is assumed.
        """
        
        # Check if zo was specified and if it already exists as an attribute.
        if zo is None:
            if hasattr(self, 'Zo'):
                zo = self.Zo
            else:
                zo = 50.             
        yo = 1./zo
        
        mattype = self.Sarray.dtype
        Y = np.zeros(self.Sarray.shape, dtype=mattype)
        for k in range(Y.shape[0]):
            s = self.Sarray[k,:,:]
            u = np.eye(s.shape[0], dtype=mattype )
            Y[k] = yo * np.dot( np.linalg.inv(u + s), (u - s) )           
        
        return Y
    
    def writeS(self, filename):
        """
        The filename argument assumes the base will be extracted from the model
        itself.
        """
        N = self.getPorts()
        fn_out = filename + '.s' + str(N) + 'p'        
        if hasattr(self, 'Zo'):
            Zo = self.Zo
        else:
            Zo = 50.
        write(fn_out, self.farray, self.Sarray, 'HZ', 'S', 'RI', Zo)        
    
    def writeZ(self, filename):
        """
        The filename argument assumes the base will be extracted from the model
        itself.
        """
        N = self.getPorts()
        fn_out = filename + '.z' + str(N) + 'p'       
        if hasattr(self, 'Zo'):
            Zo = self.Zo
        else:
            Zo = 50.
        write(fn_out, self.farray, self.getZ(Zo), 'HZ', 'Z', 'RI', Zo)

    def writecsv(self, filename=None):
        """
        Output the S-parameter data to a tabulated csv form. (only dB magnitude at
        present)
        """
        N = self.getPorts()
        M = len(self) # Returns the number of frequency points.
    
        Smag_sl = ['|S[%d,%d]|' % (n1+1,n2+1) for n1 in range(N) for n2 in range(N)]
        header = 'Frequency,' + ','.join(Smag_sl)

        Sout = header + '\n'
        for m in range(M):
            Sout += '%.6g,' % self.farray[m]
            for n1 in range(N):
                if n1 < N-1:
                    TERMCHAR = ','
                else:
                    TERMCHAR = '\n'
                v = 20.*np.log10(np.abs(self[m,n1,:]))
                Sout += ','.join(['%.6g' % vi for vi in v]) + TERMCHAR
        
        if filename is not None:
            fout = open(filename, 'w')
            fout.write(Sout)
            fout.close()
        else:
            return Sout

    def interpolate(self, frequency):
        return self.linear_interpolate(frequency)
    
    def linear_interpolate(self, frequency):
        # Use searchsorted to find the correct frequency.
        pass
    


###########################################################################

def write(fn_out, freqArray, paramArray, frequnit, datatype, dataform, Refimp):

    ## Math Space-Saving Functions #########
    def magangleEntry(S):
        x = np.abs(S)
        y = np.angle(S) * 180. / math.pi
        return x, y

    def dbEntry(S):
        x = 20. * np.log10(np.abs(S))
        y = np.angle(S) * 180. / math.pi        
        return x, y

    def realimagEntry(S):
        return S.real, S.imag

    def GetEntryPair(S, dataFormat):
        result = {
            'db': dbEntry,
            'ma': magangleEntry,
            'ri': realimagEntry
            }[dataFormat.lower()](S)
        return result
    ########################################
    
    if fn_out is None or fn_out == "":
        print( "No filename given." )
        return None
    filename = fn_out.strip()

    # Determine the number of ports.
    extstMatch = re.search('\.[syz][1-9][0-9]*p$', filename)

    n = 0 # number of ports
    if extstMatch:
        n = int(filename[(extstMatch.start()+2):(extstMatch.end()-1)])
        #print( "Touchstone with ", n, " ports." )
    else:
        print( "Port number could not be determined from the file extension." )
        return None

    # Check the port number above and make sure it does exceed the
    # amount given in the paramArray.
    if n > np.size(paramArray[0], 0) or n > np.size(paramArray[0], 1):
        print( "Only ", np.size(paramArray[0], 0), " ports detected in paramArray." )
        return None

    # Start a string for the file output.
    outStr = '! %d-port %s-parameter file \n' % (n, datatype)
    outStr += '! Created by cramsens.touchstone.write(), a Python routine: %s \n' % time.ctime()

    # Add the format line
    outStr += '# %s %s %s R %g \n' % (frequnit.upper(), datatype.upper(), dataform.upper(), Refimp)
    outStr += '! \n! \n'

    # Figure out the line dimensions for the output.
    
    if n <= 4:
        EntriesPerLine = { 1 : 1, 2 : 4, 3 : 3, 4 : 4 }[n] # ports : entries per line
    else:
        EntriesPerLine = 4

    FreqNum = np.size(freqArray, 0)
    for I in range(0, FreqNum):
        outStr += "%.12g " % freqArray[I]

        for J in range(0,n):
            elemCount = 0
            for K in range(0,n):
                if elemCount >= EntriesPerLine:
                    outStr += '\n'
                    elemCount = 0

                x, y = GetEntryPair(paramArray[I][J][K], dataform)
                outStr += "%.12g %.12g " % (x, y)
                
                elemCount += 1
            if n == 2: # 2 ports is a special case.
                pass
            else:
                outStr += '\n'
        outStr += '\n'
    outStr += '\n'
    #print( outStr )

    # Write the file.
    fid = open(filename, 'w')
    fid.write(outStr)
    fid.close()
    
    return True

###########################################################################

def read(fn_in):

    ## Math Space-Saving Functions #########
    def magangleEntry(mag, angle):
        return mag*cmath.exp(1j*angle*np.pi/180.)

    def dbEntry(dbmag, angle):
        return 10.**(dbmag/20.)*cmath.exp(1j*angle*np.pi/180.)

    def realimagEntry(x,y):
        return (x + 1j * y)

    def EntryMod(x, y, dataFormat):
        result = {
            'db': dbEntry,
            'ma': magangleEntry,
            'ri': realimagEntry
            }[dataFormat](x,y)
        return result
    ########################################
    
    if fn_in is None or fn_in == "":
        print( "No filename given." )
        return None
    filename = fn_in.strip()

    # Determine the number of ports.
    extstMatch = re.search('\.[szy][1-9][0-9]*p$', filename)

    n = 0 # number of ports
    if extstMatch:
        n = int(filename[(extstMatch.start()+2):(extstMatch.end()-1)])
        #print( "Touchstone with ", n, " ports." )
    else:
        print( "Port number could not be determined from the file extension." )
        return None

    # Open the file and split into lines.
    fid = open(filename, 'r')
    fileStrList = [s.strip().lower() for s in fid.read().splitlines() if s.strip()]
    fid.close()

    commentList = [s for s in fileStrList if re.search('^[ \t]*!', s)]
    headerList = [s for s in fileStrList if re.search('^[ \t]*#', s)]
    dataList = [s for s in fileStrList if (not re.search('^[ \t]*[!#]', s))]

    # Parse the header line.
    iline = headerList[0]
    if re.search('[ \t]+hz', iline):
        freqMult = 1.
    elif re.search('[ \t]+khz', iline):
        freqMult = 1.e3
    elif re.search('[ \t]+mhz', iline):
        freqMult = 1.e6
    elif re.search('[ \t]+ghz', iline):
        freqMult = 1.e9
        
    if re.search('[ \t]+db', iline):
        dataform = 'db';
    elif re.search('[ \t]+ma', iline):
        dataform = 'ma';
    elif re.search('[ \t]+ri', iline):
        dataform = 'ri';

    if re.search('[ \t]+s[ \t]+', iline):
        datatype = 's';
    elif re.search('[ \t]+z[ \t]+', iline):
        datatype = 'z';
    elif re.search('[ \t]+y[ \t]+', iline):
        datatype = 'y';
    elif re.search('[ \t]+h[ \t]+', iline):
        datatype = 'h';
    elif re.search('[ \t]+g[ \t]+', iline):
        datatype = 'g';

    # Determine the reference impedance.
    mRefImp = re.search('[ \t]+r[ \t]+', iline)
    if mRefImp:
        strRefImp = iline[mRefImp.end():].strip()
        if (len(strRefImp) > 0) and re.search('[0-9e]*',strRefImp):
            Refimp = float(strRefImp)
        else:
            print( "Error determining the reference impedance from line..." )
            print( iline )
            return None

    # Put all of the entries into a single flattened list.
    # Remove comments at the end of data lines.
    entryList = []
    for iline in dataList:
        m1 = re.search('!.*', iline)
        if m1:
            ilist = iline[:m1.start()].strip().split()
        else:
            ilist = iline.strip().split()

        for iarg in ilist:
            entryList.append(iarg)

    # Fill out the frequency and parameter lists
    Nentries = len(entryList)
    NumFreq = int(Nentries / (2*(n**2) + 1))

    freqArray = np.zeros(NumFreq)
    paramArray= np.zeros((NumFreq,n,n), dtype=complex)
    for I in range(0,NumFreq):
        istart = (2*(n**2) + 1) * I
        iend = (2*(n**2) + 1) * (I + 1)
        #print( I, freqArray[I] )
        #print( istart, entryList[istart], len(entryList) )
        #print( iend, entryList[iend] )
        
        freqArray[I] = float(entryList[istart])*freqMult
        
        paramdata = entryList[(istart+1):iend]
        indexCnt = 0
        for J in range(0,n):
            for K in range(0,n):
                paramArray[I,J,K] = EntryMod(float(paramdata[2*indexCnt]),
                                               float(paramdata[2*indexCnt+1]),
                                               dataform)
                #print( indexCnt, 2*indexCnt, 2*indexCnt + 1 )
                indexCnt += 1

    return freqArray, paramArray

###########################################################################

def ExtractParam(paramArray, ii, jj):
    """
    Extract a particular parameter index from the 3D array
    and output it as a 1D array.  Keep in mind that the indices start
    at 0 rather than 1. So, S11 is actually at ii = 0, jj = 0.
    """
    #FreqNum = np.size(paramArray, 0)
    #singleparam = []
    #for mat1 in paramArray:
    #    singleparam.append(mat1[ii][jj])
    #return np.array(singleparam)
    return [np.array(s)[ii][jj] for s in paramArray]

###########################################################################

def CombineParams(paramGroupArray):
    """
    Take a array/matrix of 1D arrays of the params, and put it into
    an array of 2D matrices/arrays.
    """

    FreqNum = np.size(paramGroupArray[0][0],0)
    m = np.size(paramGroupArray, 0)
    n = np.size(paramGroupArray, 1)

    outparams = np.zeros((FreqNum, m, n))
    for I in range(0,FreqNum):
        for J in range(0,m):
            for K in range(0,n):
                outparams[I][J][K] = paramGroupArray[J][K][I]
    return outparams

################################################################################

def MixedMode(S, mode = 1, convert='To Mixed Mode'):
    """
    Convert regular S-parameters to mixed-mode parameters such that.
    mode = 0 ==> p1 -> p3 and p2 -> p4
    mode = 1 ==> p1 -> p2 and p3 -> p4 (oe ordering)

    Smixed = |Sd1d1 Sd1d2 Sd1c1 Sd1c2|
             |Sd2d1 Sd2d2 Sd2c1 Sd2c2|
             |Sc1d1 Sc1d2 Sc1c1 Sc1c2|
             |Sc2d1 Sc2d2 Sc2c1 Sc2c2|

    """
    Nports = np.size(S, 1)
    
    a = None
    av = []
    if mode == 0:
        for I in range(int(Nports/2)):
            aav = [0,]*Nports
            aav[2*I] = 1.0
            aav[2*I+1] = -1.0
            av.append(aav)
        for I in range(int(Nports/2)):
            aav = [0,]*Nports
            aav[2*I] = 1.0
            aav[2*I+1] = 1.0
            av.append(aav)
    elif mode == 1:
        for I in range(int(Nports/4)):
            aav = [0,]*Nports
            aav[4*I] = 1.0
            aav[4*I+2] = -1.0            
            av.append(aav)
            aav = [0,]*Nports
            aav[4*I+1] = 1.0
            aav[4*I+3] = -1.0
            av.append(aav)
        for I in range(int(Nports/4)):
            aav = [0,]*Nports
            aav[4*I] = 1.0
            aav[4*I+2] = 1.0            
            av.append(aav)
            aav = [0,]*Nports
            aav[4*I+1] = 1.0
            aav[4*I+3] = 1.0
            av.append(aav)
            
    a = np.matrix(av)/np.sqrt(2.0)
    #print( a )
    
    #print( np.size(S,0), np.size(S,1), np.size(S,2) )
    #print( np.size(Smixed,0), np.size(Smixed,1), np.size(Smixed,2) )
    Smatrix = [np.matrix(s) for s in S]
    if convert == 'To Mixed Mode':
        Smixed = np.array([a*smat*np.transpose(a) for smat in Smatrix])        
        return Smixed
    elif convert == 'To Single-Ended Mode':
        Sse = np.array([np.transpose(a)*smat*a for smat in Smatrix])
        return Sse


################################################################################

def OddEvenPortReorderMatrix(N):
    """
    Create a matrix that reorders the ports of an S-parameter vector such
    that the odds are listed first, and the evens are listed second. So,
    [a1,a2,a3,a4] is reordered to [a1,a3,a2,a4].  This is a useful transformation
    for prepping the S-parameter matrix to convert it to T-parameters and back.
    
    Part of the assumption is that the odd ports are network inputs and the even
    ports are network outputs.  If the ordering is [1,N/2] as inputs and [N/2+1,N]
    as outputs, this reordering is not required before converting to T-parameters.
    """  
    
    if N % 2 == 1:
        print( "N (%d) must be an even number for this reordering matrix." )
        return None
    
    return np.hstack( [ np.eye(N)[:,::2], np.eye(N)[:,1::2] ] )

################################################################################

def CombineNPortTouchstones(Ports, PortFileMap):
    """
    ts filename, p1, p2, p3, p4
    """
    
    partialSparams = []
    for filemap in PortFileMap:
        print( filemap )
        thisMap = filemap[:]
        model = readmodel(thisMap.pop(0))
        
        # After popping the name off, there should just be a map from VNA port 
        # to New port number.
        partialSparams.append({'model': model,\
                               'portmap': thisMap})
    
    # Arbitrarily take the first model and assume they have the same frequency arrays.
    thisModel = partialSparams[0]['model']
    fnew = thisModel.farray
    Nfreq = len(fnew)
    
    # Create a new model with the number of ports specified above.
    snew = np.zeros((Nfreq, Ports, Ports), dtype=np.complex)
    for param in partialSparams:
        model = param['model']
        portmap = param['portmap']
        
        for m in range(len(portmap)):
            for n in range(len(portmap)):
                mnew = int(portmap[m])-1
                nnew = int(portmap[n])-1
                print( "mvna %d -> mdut %d\tnvna %d -> ndut %d" % (m, mnew, n, nnew) )
                # I'm not getting my imaginary numbers below. 
                snew[:,mnew,nnew] = model.Sarray[:,m,n]
    
    write('test.s%dp' % Ports, fnew, snew, 'hz', 's', 'ri', 50)
    #write(fn_out, freqArray, paramArray, frequnit, datatype, dataform, Refimp)
    
    return fnew, snew
        
################################################################################

def PortReorder(reordering_l, S):
    N = np.size(S,1)
    lam = np.zeros((len(reordering_l),N))
    
    for ii, old_port in enumerate(reordering_l):
        lam[ii, old_port-1] = 1.

    Sreordered = np.array([np.dot(np.dot(lam, s), lam.T) for s in S])
    
    return Sreordered

################################################################################
dB = lambda y : 20.0*np.log10(np.abs(y))

def ReturnLoss2x2(model_l, plot_args_l, **kwargs):
    plot_mixedmode = kwargs.get('mixedmode', True)
    
    if plot_mixedmode:
        for mi in model_l:
            mi.toMixedMode() # set to mixed mode if not alerady.
        port_l = ['d1', 'd2', 'c1', 'c2'] # Assume 4 port configurations for this routine.
    else:
        for mi in model_l:
            mi.toSingleEnded()
        port_l = range(1,5)
    
    trace_l = zip(port_l,port_l) # Return loss observes the same port as the input.
    return Plot2x2(model_l, plot_args_l, trace_l, **kwargs)

# def InsertionLoss2x1(model_l, plot_args_l, **kwargs):
#     plot_mixedmode = kwargs.get('mixedmode', True)
#     
#     if plot_mixedmode:
#         for mi in model_l:
#             mi.toMixedMode() # set to mixed mode if not alerady.
#         port_l = [('d2','d1'), ('c2','c1')] # Assume 4 port configurations for this routine.
#     else:
#         for mi in model_l:
#             mi.toSingleEnded()
#         port_l = [(2,1), (4,3)]
#     
#     # Return loss observes the same port as the input.
#     return Plot2x1(model_l, plot_args_l, port_l, **kwargs)

def TransModeConversion2x2(model_l, plot_args_l, **kwargs):
    mixedmode_previously = []
    for mi in model_l:
        mixedmode_previously.append(mi.mixedmode)
        mi.toMixedMode()
    
    trace_l = [('c2','d1'), ('d2','c1'),('c1','d2'),('d1','c2')]
    return Plot2x2(model_l, plot_args_l, trace_l, **kwargs)

def ReflModeConversion2x2(model_l, plot_args_l, **kwargs):
    mixedmode_previously = []
    for mi in model_l:
        mixedmode_previously.append(mi.mixedmode)
        mi.toMixedMode()
    
    trace_l = [('c1','d1'), ('d1','c1'),('c2','d2'),('d2','c2')]
    return Plot2x2(model_l, plot_args_l, trace_l, **kwargs)

def Plot2x2(model_l, plot_args_l, trace_l, **kwargs):
    from matplotlib.pyplot import figure    
    from cramsens.plot_tools import remove_xtick_labels, set_xlim, set_ylim, savefigure
    
    fig = figure()
    fig.subplots_adjust(left=0.09, bottom=0.08, right=0.91, top=0.92, 
                        hspace=0.08, wspace=0.07)
    
    if 'title' in kwargs: fig.suptitle(kwargs['title'], fontsize=16)
          
    ax1 = fig.add_subplot(2,2,1)
    remove_xtick_labels(ax1)
    ax2 = fig.add_subplot(2,2,2,sharex=ax1)
    ax2.yaxis.tick_right()
    ax2.yaxis.set_label_position("right")
    remove_xtick_labels(ax2)
    ax3 = fig.add_subplot(2,2,3,sharex=ax1)
    ax4 = fig.add_subplot(2,2,4,sharex=ax1)
    ax4.yaxis.tick_right()
    ax4.yaxis.set_label_position("right")
    ax_l = [ax1, ax2, ax3, ax4]
    
    
    
    xlabels = ['',]*2 + ['Frequency (GHz)',]*2
    ylabels = []    
    for pout, pi in trace_l:
        ylabels.append(r'$|S_{%s%s}|$ (dB)' % (str(pout).upper(),str(pi).upper()))
        
    for ii, axi in enumerate(ax_l):
        axi.grid(True)
        axi.set_xlabel(xlabels[ii])
        axi.set_ylabel(ylabels[ii])
    
    xfcn = kwargs.get('xfcn', lambda x : x/1e9)
    yfcn = kwargs.get('yfcn', dB)
    for mi, pai in zip(model_l, plot_args_l):
        for axi, porti in zip(ax_l,trace_l):
            axi.plot(xfcn(mi.farray), yfcn(mi.getTrace(porti[0],porti[1])), **pai)
    
    if 'ylim' in kwargs:
        for axi in ax_l:
            set_ylim(axi, kwargs['ylim'])
        
    if 'ylim_l' in kwargs:
        for axi, ylim_i in zip(ax_l, kwargs['ylim_l']): 
            set_ylim(axi, ylim_i)
            
    if 'xlim' in kwargs:
        set_xlim(ax1, kwargs['xlim'])

    ax1.legend(loc='best')

    if 'filename' in kwargs: savefigure(fig, kwargs['filename'])
                
    fig.show()
    return fig

        
def zin_from_s21(model, port1=1, port2=2, zo=50.0):
    """Based on the expression in Istvan Novak's '99 DesignCon presentation,
    this routine uses the expression on page 30 of his presentation.  Lp = 0
    is assumed, since it is--for practical purposes--calibrated out of the 
    measurement."""
    
    s21 = model.getTrace(port1, port2)
    zin = s21 * (zo/2.0) / (1.0 - s21)
    
    return zin

def ztrans_from_s21(model, port1=1, port2=2, zo=50.0):
    """similar"""
    s21 = model.getTrace(port1, port2)

###############################################################################
def renormalize(S, Zprev, Znew):    
    inv = np.linalg.inv
    
    Snew = np.zeros(S.shape, dtype=np.complex128)
    
    Zp = np.array(Zprev)
    Zn = np.array(Znew)
    R = np.matrix(np.diag(((Zn - Zp)/(Zn + Zp))))    
    A = np.matrix(np.diag((np.sqrt(Zn/Zp)/(Zn+Zp))))
    I = np.matrix(np.eye(Snew.shape[1]))
    
    for ii in range(Snew.shape[0]):
        s0 = S[ii,:,:]
        s1 = inv(A) * (s0-R) * inv((I - (R*s0))) * A
        Snew[ii,:,:] = s1
        
    return Snew
    
    
    
    
    
    
    
    