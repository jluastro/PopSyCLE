#! /usr/bin/env python
#-----------------------------------------------------------------------------
#
#EBF (Efficient Binary Format) Software Library and Utilities
#Copyright (c) 2012 Sanjib Sharma
#All rights reserved.

# Copyright Notice and License Terms for 
# EBF (Efficient Binary Format) Software Library and Utilities
# -----------------------------------------------------------------------------

# EBF (Efficient Binary Format) Software Library and Utilities
# Copyright (c) 2012 Sanjib Sharma
# All rights reserved.

# Redistribution and use in source and binary forms, with or without 
# modification, are permitted for any purpose (including commercial purposes) 
# provided that the following conditions are met:

# 1. Redistributions of source code must retain the above copyright notice, 
#    this list of conditions, and the following disclaimer.

# 2. Redistributions in binary form must reproduce the above copyright notice, 
#    this list of conditions, and the following disclaimer in the documentation 
#    and/or materials provided with the distribution.

# 3. In addition, redistributions of modified forms of the source or binary 
#    code must carry prominent notices stating that the original code was 
#    changed and the date of the change.

# 4. All publications or advertising materials mentioning features or use of 
#    this software are asked, but not required, to acknowledge 
#    and credit the contributors.

# 5. Neither the name of the copyright holders nor the names of its contributors 
#    may be used to endorse or promote products derived from this software 
#    without specific prior written permission.


# DISCLAIMER: 
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
# ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
#-----------------------------------------------------------------------------
#
# This file is part of EBF.  The full EBF copyright notice, including       
# terms governing use, modification, and redistribution, is contained in    
# the files COPYING and COPYRIGHT,  which can be found at the root        
# of the source code distribution tree.     
#--------------------------------------------------------------------------    

# Copied https://github.com/segasai/ebfpy on 2019-04-24

"""
A module to read and write data in ebf format. 

 .. moduleauthor: Sanjib Sharma <bugsanjib at gmail com>

EBF is a binary format for storing data. It is designed to   
read and write data, easily and efficiently. 

- Store multiple data items in one file, each having a unique tag name

  + tagnames follow the convention of unix style pathname e.g. /x or /mydata/x
  + this allows hierarchical storage of data

- Automatic type and endian conversion  
- Support for mutiple programming languages

  + data can easily read in C, C++, Fortran, Java, IDL and Matlab
  + facilitates easy distribution of data 

- Comprehensive numpy support

  + data is read back as numpy arrays
  + almost any numpy array can be written
  + Nested numpy structures are also supported

- Read and write directly a recursive dictionary of numpy arrays

To install  
::

$pip install ebfpy           OR
$pip install ebfpy --user    OR

Alternatively
::

$tar -zxvf ebfpy_x.x.x.tar.gz
$cd ebfpy_x.x.x
$python setup.py install --user                            OR 
$python setup.py install --user --install-scripts=mypath   OR
$python setup.py install  --install-scripts=mypath 


The --install_scripts option if specified 
determines the installation location of the command line script ebftkpy, 
the ebf module is always installed in a standard location. 
It is better to set this manually (to something like '/usr/local/bin' 
or somewhere in home) because the standard script installation location might 
not be in your search path. With *--user* option generally the scripts are 
installed in *~/.local/bin/*.

To run the test suite just do (from within folder ebfpy_x.x.x)
::

$./ebf.py 

Example:

Write specific numpy arrays.

>>> import ebf
>>> import numpy
>>> x = numpy.random.rand(2,5)
>>> y = numpy.random.rand(2,5)
>>> ebf.write('check.ebf', '/x', x, "w")
>>> ebf.write('check.ebf', '/y', y, "a")

Write in a different path within an ebf file .

>>> ebf.write('check.ebf', '/mypath/x', x, "a")
>>> ebf.write('check.ebf', '/mypath/y', y, "a")

Read back the written arrays

>>> x1 = ebf.read('check.ebf', '/x')
>>> y1 = ebf.read('check.ebf', '/mypath/y')

Read all items in an ebf path as a dictionary
such that data["x"] is same as x1
such that data["y"] is same as y1

>>> data = ebf.read('check.ebf', '/mypath/')

Check the contents of the file.

>>> ebf.info('check.ebf')
check.ebf 2460 bytes
------------------------------------------------------------------
name                           dtype    endian  unit       dim       
------------------------------------------------------------------
/.ebf/info                     int64    little             [5]       
/.ebf/htable                   int8     little             [1256]    
/x                             float64  little             [2 5]     
/y                             float64  little             [2 5]     
/mypath/x                      float64  little             [2 5]     
/mypath/y                      float64  little             [2 5]     

    
Split a structure and write individual data items in 
path "/mypath/" in an ebf file.

>>> dth = numpy.dtype([('data_u1', 'u1', (2, 5)), ('data_u2', 'u2', (2, 5))])
>>> data = numpy.zeros(1, dtype = dth)
>>> ebf.write('check.ebf', '/mypath/', data, "w")
>>> data1 = ebf.read('check.ebf', '/mypath/')
>>> ebf.info('check.ebf') 
check.ebf 1906 bytes
------------------------------------------------------------------
name                           dtype    endian  unit       dim       
------------------------------------------------------------------
/.ebf/info                     int64    little             [5]       
/.ebf/htable                   int8     little             [1256]    
/mypath/data_u1                uint8    little             [2 5]     
/mypath/data_u2                uint16   little             [2 5] 


Write a nested structure and read it back.

>>> dth = numpy.dtype([('data_u1', 'u1', (2, 5)), ('data_u2', 'u2', (2, 5))])
>>> dth1 = numpy.dtype([('data_u1', 'u1', (2, 5)), ('point1', dth, (1, ))])
>>> data = numpy.zeros(10, dtype = dth1)
>>> ebf.write("check.ebf", "/data", data, "w")
>>> data1 = ebf.read("check.ebf", "/data")
>>> ebf.info("check.ebf")
check.ebf 2247 bytes
------------------------------------------------------------------
name                           dtype    endian  unit       dim       
------------------------------------------------------------------
/.ebf/info                     int64    little             [5]       
/.ebf/htable                   int8     little             [1256]    
/data                          struct   little             [10]      
structure definition:
ver-1 
struct {
uint8 data_u1 2  2  5  ;
struct {
uint8 data_u1 2  2  5  ;
uint16 data_u2 2  2  5  ;
} point1 1 1   ; 
} anonymous 1 1   ; 

Write a string and read it back as string.
Note, return type is numpy.ndarray, hence have to use tostring() 
method to convert it back to string.

>>> x = "abcdefghijkl"
>>> ebf.write("check.ebf", "/mystr", numpy.array(x), "w")
>>> y = ebf.read("check.ebf", "/mystr").tostring()

Write a list of string and read it back as numpy.ndarray of type numpy.string

>>> x = ["abc", "abcdef"]
>>> ebf.write("check.ebf", "/mystr", numpy.array(x), "w")
>>> y = ebf.read("check.ebf", "/mystr")
>>> print y[0] == "abc",y[1] == "abcdef"
True True


Write with units and read it back.

>>> data = numpy.zeros(1, dtype = "int32")
>>> ebf.write('check.ebf', '/data', data, "w",dataunit="100 m/s")
>>> print, ebf.unit('check.ebf', '/data')


Check if a data item is present.

>>> ebf.containsKey('check.ebf', '/data')


"""

# 0.0.9  with open replaced with try and finally in update_ind for pyver<2.5
# 0.0.8 ebf.update_ind()  added
# 0.0.7 EbfHeader modified to take into account writing structures with non native endian format 
# In 0.0.4 write mode 'u' and 'e' added
# update mode works to replace existing data of same size and shape
# expand mode works to extend data size, first dim only. Provided other dims are same and data is 
# of same enidianness
# In 0.0.3 recon can be 1 and 2 
# path without data but additional dirs had problem in recon. This is solved 
# _getHeader replaced with getHeader
# update_ind has option ind=None
# write() with mode u and e, has been modified now tests for dtype match and allows scalar to be passed
# think about precision preserving csv,ssv
# added following in __LoadMap to prevent false entries 
#        except RuntimeError:
#            del _EbfMap.ltable[filename]
#            raise RuntimeError('Cannot Read File:'+filename)    
# 2014Oct27 added in def keys() ability to read structure dtype.names
# 2014ONov13 added dictnpstruct and npstruct2dict and outpath option in copy
# 2014ONov14 improved the cat module so that there is no loss of precision when printing 
# April 2015 all None comparison changed to is not and is 
# Overflow Runtime warning suppressed in ebflthash using np.seterr


import numpy
import sys
import time
import os
from itertools import groupby
from operator import itemgetter


__version__ = "0.0.14"


class _EbfUtils(object):
        
    @staticmethod
    def createPathNode(name):
        """
        Create a path node for path Tree
        """
        node={}
        node['files']=[]
        node['dirs']={}
        node['name']=name
        return node
    
    @staticmethod
    def addPathToTree(node,path_list):
        """
        add a list of paths to Tree through root node
        """
        for path1 in path_list:
            if path1.count('/') == 0:
                node['files'].append(path1)
            else:
                x=path1.split('/',1)
                if not (x[0]+'/' in node['dirs']):
                    node['dirs'][x[0]+'/']=_EbfUtils.createPathNode(node['name']+x[0]+'/')
                _EbfUtils.addPathToTree(node['dirs'][x[0]+'/'],[x[1]])    
                
    @staticmethod
    def printPathTree(node):
        """
        print the tree
        """
        if node['name'] != "":
            print('ls ',node['name'])
            print(node['files'])
            print(list(node['dirs'].keys()))
            print() 
            for x in list(node['dirs'].keys()):
                _EbfUtils.printPathTree(node['dirs'][x])
    
    @staticmethod
    def searchPathTree(node,dataname):
        if node['name'] == dataname :
            return node
        else:
            for key in list(node['dirs'].keys()):
                nodef=_EbfUtils.searchPathTree(node['dirs'][key],dataname)
                if nodef['name'] == dataname:
                    return nodef        
        return node

    @staticmethod
    def getKeysRecursive(node):
        """
        print the tree
        """
        keys=[]
        for key in node['files']:
            keys.append(node['name']+key)
        for x in list(node['dirs'].keys()):
            keys=keys+_EbfUtils.getKeysRecursive(node['dirs'][x])
        return keys

    @staticmethod
    def get_byteorder(data):
        if sys.byteorder == 'little':
            sorder='<'
        else:
            sorder='>'
        if data.dtype.names is not None:
            x=[]
            for key in data.dtype.names:
                if data[key].dtype.byteorder == '<':
                    x.append('<')
                elif data[key].dtype.byteorder == '>':
                    x.append('>')
                elif data[key].dtype.byteorder == '=':
                    x.append(sorder)
            x=numpy.array(x)
            if x.size>0:
                if numpy.where(x==x[0])[0].size != x.size:
                    raise RuntimeError("EBF: error in all fields are not of same byte order")
                return x[0]
            else:
                return '|'
        else:
            if data.dtype.byteorder == '=':
                return sorder
            else:
                return data.dtype.byteorder



class _EbfMap(object):
    """
    A class that is used to get location of data objects in an ebf file
    """
    ltable={}
    @staticmethod
    def printout():
        print(_EbfMap.ltable)
        
    @staticmethod
    def __loadMap(filename):
        """
        load the location table
        """
        if filename in _EbfMap.ltable:
            del _EbfMap.ltable[filename]
        _EbfMap.ltable[filename] = {}
            
            
        
        try:
            if os.path.isfile(filename)==False:
                raise RuntimeError('File not found:'+filename)    
                
            fp1 = open(filename, 'rb')
            fp1.seek(0, 2)
            filesize = fp1.tell()
            fp1.seek(0, 0)
    
            while fp1.tell() < filesize:
                mypos = fp1.tell()
                header = _EbfHeader()
                header.read(fp1)
                _EbfMap.ltable[filename][header.name.lower()] = mypos
                fp1.seek(header.capacity(), 1)
    
            if fp1.tell() != filesize:
                fp1.close()
                del _EbfMap.ltable[filename]
                raise RuntimeError('EBFCorrupt')    
            else:
                fp1.close()
            key_list=list(_EbfMap.ltable[filename].keys())
            _EbfMap.ltable[filename]['pathtree']=_EbfUtils.createPathNode('')
            _EbfUtils.addPathToTree(_EbfMap.ltable[filename]['pathtree'],key_list)    
            _EbfMap.ltable[filename]['checksum'] = _EbfMap.getCheckSum(filename)
        except RuntimeError:
            del _EbfMap.ltable[filename]
            raise RuntimeError('Cannot Read File:'+filename)    
            

    @staticmethod
    def keys(filename):
        """
        check if a data object exists in a file
        loads file into map if it does not exist
        """
        if (filename in _EbfMap.ltable) == 0 :
            _EbfMap.__loadMap(filename)
        if _EbfMap.ltable[filename]['checksum'] != _EbfMap.getCheckSum(filename):
            _EbfMap.__loadMap(filename)
                
        keys1=list(_EbfMap.ltable[filename].keys())
        if (keys1.count('checksum') > 0):   
            keys1.remove('checksum')
        if (keys1.count('pathtree') > 0):   
            keys1.remove('pathtree')
        return keys1
    


    @staticmethod
    def get(filename, dataname,option=1):
        """
        get the location of the data object
        """
        if (filename in _EbfMap.ltable) == 0 :
            _EbfMap.__loadMap(filename)
        elif option == 1:    
            if _EbfMap.ltable[filename]['checksum'] != _EbfMap.getCheckSum(filename):
                _EbfMap.__loadMap(filename)
        #1/0
        if dataname in _EbfMap.ltable[filename]:
            return _EbfMap.ltable[filename][dataname]
        else:
            return -1

    @staticmethod
    def getCheckSum(filename):
        fp1 = open(filename, 'rb')
        checksum=numpy.int64(0)
        data = fp1.read(1)
        if (data != ""):
            fp1.seek(0,0)
            header = _EbfHeader()
            header.read(fp1)
            if ((header.datatype == 3)&(header.name == '/.ebf/info')):
                checksum = numpy.fromstring(fp1.read(8), dtype = 'int64')            
                if header.flagswap == 1:
                    checksum = checksum.byteswap(True)
                checksum=checksum[0]
        fp1.close()                    
        return checksum
    
    @staticmethod
    def clear():
        """
        remove a file from map
        """
        _EbfMap.ltable={}
    

#    @staticmethod
#    def __put(filename,dataname,location):
#        fp1 = open(filename, 'rb+')
#        checksum=numpy.int64(0)
#        data = fp1.read(1)
#        if (data != ""):
#            fp1.seek(0,0)
#            header = _EbfHeader()
#            header.read(fp1)
#            datapos=fp1.tell()
#            if ((header.datatype == 3)&(header.name == '/.ebf/hinfo')):
#                checksum = numpy.fromstring(fp1.read(header.elements()*header.datasize), dtype = 'int64')            
#                if header.flagswap == 1:
#                    checksum = checksum.byteswap(True)
#                mystr='('+dataname+', '+str(location)+') '
#                checksum[0]=numpy.int64(_EbfTable.ebfckhash(mystr,checksum[0]))
#                if header.flagswap == 1:
#                    checksum = checksum.byteswap(True)
#                fp1.seek(datapos,0)
#                fp1.write(checksum.tostring('C'))            
#        fp1.close()                    
#        
#        _EbfMap.ltable[filename][dataname.lower()] = location
#        _EbfMap.ltable[filename]['checksum'] = checksum[0]
#        _EbfUtils.addPathToTree(_EbfMap.ltable[filename]['pathtree'],[dataname.lower()])



class _TypeManager(object):
    """
    A class to convert data type strings to ebf data type integer codes and vice versa
    """
    typedicts = {}
    typedictl = {}
    typelists = ['', 'S1', 'i4', 'i8', 'f4', 'f8', 'i2', '', 'V', 'i1', 'u1', 'u2', 'u4', 'u8']
    typelistl = ['', 'char', 'int32', 'int64', 'float32', 'float64', 'int16', '', 'struct', 'int8', 'uint8', 'uint16', 'uint32', 'uint64']
    for i in range(len(typelistl)):
        if(typelists[i] != ''):
            typedicts[typelists[i]] = i               
            typedictl[typelistl[i]] = i
    typedicts['b1']=9
    @staticmethod        
    def stoi (typename):
        """
        python type string to ebf type code    
        """
        if typename in _TypeManager.typedicts:
            return _TypeManager.typedicts[typename]                               
        if typename in _TypeManager.typedictl:
            return _TypeManager.typedictl[typename]
        else:
            print('datatype=',typename)
            raise RuntimeError("Ebf error: unrecognized data type ")

    @staticmethod        
    def containsKey (typename):
        """
        check if if type string is a valid supported type by ebf    
        """
        if typename in _TypeManager.typedicts:
            return True                               
        if typename in _TypeManager.typedictl:
            return True
        return False
        
    @staticmethod        
    def itos_s (i):
        """
        ebf type code to python type string short form   
        """
        if (i < 1)|(i == 7): 
            raise RuntimeError("Ebf error: unrecognized data type index "+str(i))
        return _TypeManager.typelists[i]                               
    
    @staticmethod        
    def itos_l (i):
        """
        ebf type code to python type string long form   
        """
        if (i < 1)|(i == 7) :
            raise RuntimeError("Ebf error: unrecognized data type index "+str(i))
        return _TypeManager.typelistl[i]                               
    
    @staticmethod        
    def stos_l (typename):
        """
        python type string short to long form   
        """
        return _TypeManager.typelistl[_TypeManager.stoi(typename)]
        
    @staticmethod        
    def stos_s (typename):
        """
        python type string long to short form   
        """
        return _TypeManager.typelists[_TypeManager.stoi(typename)]





class _EbfHeader:
    """
    a class to read and write ebf data object headers
    """
    def get_dtype(self):
        if self.datatype == 8:
            dth1 = numpy.dtype(sdef2descr(self.sdef)[0])
        else:
            dth1 = numpy.dtype(_TypeManager.itos_s(self.datatype))
            if self.datatype == 7:
                dth1 = numpy.dtype('S'+str(self.datasize))
            if self.datatype == 1:
                dth1 = numpy.dtype('S'+str(self.dim[-1]))
        return dth1
        
    def capacity(self):
        """
        size in bytes of data item in an ebf file
        """
        return self.capacity_
    
    def elements(self):
        """
        size in bytes of data item in an ebf file
        """
        return numpy.prod(self.dim)

    def getshape(self):
        """
        shape of data item in an ebf file
        """
        return list(self.dim)
    

        
    def read(self, fp1):
        """
        read the header from file
        """
        sig0 = numpy.fromstring(fp1.read(3), dtype = 'S3')[0]
        fp1.seek(-3, 1)
        sig1 = numpy.fromstring(fp1.read(8), dtype = 'S8')[0]
        version = numpy.fromstring(fp1.read(4), dtype = 'int8')
        fp1.seek(-12, 1)
        if sig0 == 'EBF' :                    
            self.__read10(fp1)
        elif (version[0] == 1):
            self.__read11(fp1)
        else:
            raise RuntimeError('EBF unrecognized header format')
    


    def __read10(self, fp1):
        """
        read the header from file
        """
        sig = numpy.fromstring(fp1.read(6), dtype = 'int8')
        self.name = numpy.fromstring(fp1.read(100), dtype = 'S100')[0]
        self.name=self.name.replace('\x00',' ')
        self.name=self.name.strip().lower().decode('ascii')
        
        unused = numpy.fromstring(fp1.read(2), dtype = 'int8')
        unused = numpy.fromstring(fp1.read(36), dtype = 'S1')
        self.endiantest = numpy.array(numpy.fromstring(fp1.read(4), dtype = 'int32')[0])
        self.datatype = numpy.array(numpy.fromstring(fp1.read(4), dtype = 'int32')[0])
        self.datasize = numpy.array(numpy.fromstring(fp1.read(4), dtype = 'int32')[0])
        rank = numpy.array(numpy.fromstring(fp1.read(4), dtype = 'int32')[0])
        unused = numpy.fromstring(fp1.read(32), dtype = 'i4')
        self.dim = numpy.fromstring(fp1.read(64), dtype = 'i8')
        self.flagswap = 0
        self.headersize = numpy.array(256, dtype = "int32")
#        self.unitsize = numpy.array(0,dtype = "int32")
#        self.namesize = numpy.array(len(self.name), dtype = "int32")


#        if sig.tostring() != 'EBF>>>':
#            raise RuntimeError('EBF file signature error ')
        
        if self.endiantest != 256:
            self.flagswap = 1
            self.endiantest.byteswap(True)
            if self.endiantest != 256:
                fp1.close()
                raise RuntimeError('EBF unrecognized header')
            self.datatype.byteswap(True)
            self.datasize.byteswap(True)
            rank.byteswap(True)
            self.dim.byteswap(True)

        self.dim.resize((rank, )) 
        self.dataunit = ""
        self.sdef = ""
        self.capacity_=numpy.array(self.elements()*self.datasize,dtype='int64')

        if self.datatype > 13:
            fp1.close()
            raise RuntimeError('EBF Type Code unrecognized')


    def __read11(self, fp1):
        """
        read the header from file
        """
        self.flagswap = 0
        self.headerpos = fp1.tell()
        sig = numpy.fromstring(fp1.read(8), dtype = 'int8')
        sig1 =  numpy.array((-118, 69, 66, 70, -82, 43, -81, 10), dtype = "int8")
        if (sum(sig == sig1)!=8):
            fp1.close()
            raise RuntimeError('Ebf signature does not match')
                                 
        self.version = numpy.fromstring(fp1.read(4), dtype = 'int8')        
        temp = numpy.fromstring(fp1.read(4*8), dtype = 'int32')
        self.flags = numpy.fromstring(fp1.read(4), dtype = 'int8')
        self.capacity_ = numpy.fromstring(fp1.read(8), dtype = 'int64')[0]
        
        if temp[0] != 1684234849:
            temp.byteswap(True)
            self.capacity_=self.capacity_.byteswap()
            self.flagswap = 1
            if temp[0] != 1684234849:
                fp1.close()
                raise RuntimeError('EBF unrecognized header')
        
        self.endiantest = numpy.array(temp[0])
        self.headersize = numpy.array(temp[1])
        namesize = numpy.array(temp[2])
        self.datatype = numpy.array(temp[3])
        self.datasize = numpy.array(temp[4])
        rank = numpy.array(temp[5])
        unitsize = numpy.array(temp[6])
        sdefsize = numpy.array(temp[7])                

        if rank > 0:
            self.dim = numpy.fromstring(fp1.read(8*rank), dtype = 'int64')
            if self.flagswap == 1:
                self.dim.byteswap(True)
        else:
            self.dim=numpy.array((),'int64')
#            self.dim=numpy.zeros(0,'int64')
            
        self.name = numpy.fromstring(fp1.read(namesize), dtype = 'S1').tostring()
        self.name=self.name.lower().decode('ascii')
        self.dataunit = ""
        if unitsize > 0:
            self.dataunit = numpy.fromstring(fp1.read(unitsize), dtype = 'S1').tostring().decode('ascii')
            
        self.sdef = ""
        if sdefsize > 0:
            self.sdef = numpy.fromstring(fp1.read(sdefsize), dtype = 'S1').tostring().decode('ascii')


        self.datapos = self.headerpos+self.headersize    
        
        if self.datapos < fp1.tell():
            fp1.close()
            raise RuntimeError('EBF file error reading header')
        else:
            fp1.seek(self.datapos)     

        if self.datatype > 13:
            fp1.close()
            raise RuntimeError('EBF Type Code unrecognized')
        

    def rename(self,name):
        extrasize=self.headersize-(56 + len(name) + len(self.dataunit) + len(self.sdef) +  8 * self.dim.size)
        extrasize =extrasize- (len(name)-len(self.name))
        if extrasize < 2:
            raise RuntimeError('EBF: Not enough extra space in header to rename')
        self.name=name.lower()        
        


    def write(self, fp1):
        """
        write the header to file
        """
        self.headerpos = fp1.tell()
        self.__write11(fp1)
        self.datapos = fp1.tell()
        if self.headersize != (self.datapos-self.headerpos):
            fp1.close()
            raise RuntimeError('EBF error while writing header mismatch of size')
            

    def __write10(self, fp1):
        """
        write the header to file
        """
        if len(self.name) > 100: 
            fp1.close()
            raise RuntimeError('Data object name too large')
        if self.dim.size > 8: 
            fp1.close()
            raise RuntimeError('Data object rank too large')
        rank = numpy.array(self.dim.size,dtype='int32')
        if rank == 0:
            rank=numpy.array(1,dtype='int32')
        
        if self.flagswap == 1:
            rank.byteswap(True)            
        sig = "EBF>>>"    
        fp1.write(sig)
        fp1.write(self.name)
        if len(self.name) < 100: 
            unused = numpy.zeros(100-len(self.name),  dtype = "int8")
            fp1.write(unused.tostring('C'))    
        version = numpy.array([1, 1], dtype = "int8")
        fp1.write(version.tostring('C'))
        unused = numpy.zeros(36,  dtype = "int8")
        fp1.write(unused.tostring('C'))
        fp1.write(self.endiantest.tostring('C'))
        fp1.write(self.datatype.tostring('C'))
        fp1.write(self.datasize.tostring('C'))
        fp1.write(rank.tostring('C'))
        unused = numpy.zeros(32, dtype = "int8")
        fp1.write(unused.tostring('C'))
        if self.dim.size == 0:
            dim = numpy.array(1, dtype = "int64")
            fp1.write(dim.tostring('C'))
            unused = numpy.zeros(8*7, dtype = "int8")
            fp1.write(unused.tostring('C'))            
        else:
            fp1.write(self.dim.tostring('C'))
            if self.dim.size < 8:
                unused = numpy.zeros(8*(8-self.dim.size), dtype = "int8")
                fp1.write(unused.tostring('C'))

    def __write11(self, fp1):
        """
        write the header to file- general form for arrays
        """
#        pos = fp1.tell()
        rank = numpy.array(self.dim.size, dtype = "int32")
        namesize = numpy.array(len(self.name), dtype = "int32")
        unitsize = numpy.array(len(self.dataunit), dtype = "int32")
        sdefsize = numpy.array(len(self.sdef), dtype = "int32")

        sig =  numpy.array((-118, 69, 66, 70, -82, 43, -81, 10), dtype = "int8")
        version = numpy.array([1, 1, 0, 0], dtype = "int8")
        
        if self.flagswap == 1:
            temp=sig.tostring('C')+version.tostring('C')+self.endiantest.byteswap().tostring('C')
            temp+=self.headersize.byteswap().tostring('C')+namesize.byteswap().tostring('C')+self.datatype.byteswap().tostring('C')
            temp+=self.datasize.byteswap().tostring('C')+rank.byteswap().tostring('C')+unitsize.byteswap().tostring('C')+sdefsize.byteswap().tostring('C')
            temp+=self.flags.tostring('C')+self.capacity_.byteswap().tostring('C')+self.dim.byteswap().tostring('C')
            temp+=self.name.encode('ascii')+self.dataunit.encode('ascii')+self.sdef.encode('ascii')
        else:
            temp=sig.tostring('C')+version.tostring('C')+self.endiantest.tostring('C')
            temp+=self.headersize.tostring('C')+namesize.tostring('C')+self.datatype.tostring('C')
            temp+=self.datasize.tostring('C')+rank.tostring('C')+unitsize.tostring('C')+sdefsize.tostring('C')
            temp+=self.flags.tostring('C')+self.capacity_.tostring('C')+self.dim.tostring('C')
            temp+=self.name.encode('ascii')+self.dataunit.encode('ascii')+self.sdef.encode('ascii')
        
        
        extra = self.headersize-len(temp)
    
        if extra<0:
            fp1.close()
            raise RuntimeError('EBF error while writing header mismatch of size')
        else:
            if extra>0:
                temp1 = numpy.zeros(extra, dtype = "int8")+60
                temp+=temp1.tostring()                
            fp1.write(temp)            

        
 

    def create(self, tagname, data, dataunit, sdef):
        """
        create the header from numpy array
        """
                    
        self.name = tagname.strip().lower()
        self.dim = numpy.array(data.shape, dtype = "int64")
        '''Treat Scalars as rank 1 and dim[0] = 1'''
#        if data.ndim == 0:
#            self.dim = numpy.array([1], dtype = "int64")
#        else:
#            self.dim = numpy.array(data.shape, dtype = "int64")
            
        self.datasize = numpy.array(data.itemsize, dtype = "int32")
        self.flags = numpy.array([0, 0, 0, 0], dtype = "int8")

        if data.dtype.char == 'V':
            '''Voids'''
            self.datatype = numpy.array(_TypeManager.stoi('V'), dtype = "int32")
        elif data.dtype.char == 'S':
            
            '''Strings make S1 and adjust rank, dim, datasize'''
            self.datatype = numpy.array(_TypeManager.stoi('S1'), dtype = "int32")
            self.datasize = numpy.array(1, dtype = "int32")
            
            self.dim=list(data.shape)
            self.dim.append(data.itemsize)
            self.dim=numpy.array(self.dim,dtype='int64')
            
#            if self.dim.size == 0:
#                self.dim = numpy.array([1], dtype = "int64")                
#            if data.itemsize > 1:
#                if(self.dim[self.dim.size-1] == 1):
#                    self.dim[self.dim.size-1] = data.itemsize
#                else:
#                    rank = self.dim.size+1
#                    self.dim.resize((rank, ))
#                    self.dim[rank-1] = data.itemsize
            ''' type 7 not implemented '''
    #            self.datatype = numpy.array(typedict['S'], dtype="int32")
    #            self.datasize = numpy.array(data.itemsize, dtype="int32")
            
        else:
            '''Rest'''
            self.datatype = numpy.array(_TypeManager.stoi(data.dtype.str[1:3]), dtype = "int32")
#            if typedict.has_key(data.dtype.str[1:3]):
#                self.datatype = numpy.array(typedict[data.dtype.str[1:3]], dtype="int32")
#            else:
#                raise RuntimeError('datatype ' + data.dtype.str + ' is not supported')

#        if (typelist.count(data.dtype.str[1:3]) == 0) & (data.dtype.char != 'V') & (data.dtype.char != 'S'):
#            raise RuntimeError("Data type" + data.dtype.char + " not supported")

        self.flagswap = 0
        if data.dtype.char == 'V':
            border=_EbfUtils.get_byteorder(data)
#            if sys.byteorder == 'little':
#                borderd='<'
#            else:
#                borderd='>'
#
#            border_a=[]                
#            for key in data.dtype.names:
#                if data[key].dtype.byteorder == '>':
#                    border_a.append('>')
#                elif data[key].dtype.byteorder == '<':
#                    border_a.append('>')
#                elif data[key].dtype.byteorder == '=':
#                    border_a.append(borderd)
#                    
#            border_a=numpy.array(border_a)
#            if numpy.where(border_a==border_a[0])[0].size != border_a.size:
#                raise RuntimeError("EBF: error in _EbfHeader.create all fields are not of same byte order")
#            border=border_a[0]
                
            if (sys.byteorder == 'little' and (border == '>' )) or (sys.byteorder == 'big' and (border == '<')):
                self.flagswap = 1        
        else:
            if (sys.byteorder == 'little' and (data.dtype.byteorder == '>' )) or (sys.byteorder == 'big' and (data.dtype.byteorder == '<')):
                self.flagswap = 1

        self.dataunit = dataunit    
        self.endiantest = numpy.array(numpy.int32(1684234849))
        self.sdef = sdef    
#        self.fields = fields

#        if (self.rank == 1) & (self.dim[0] == 1) & (self.datatype != 8):    
#            self.extra = numpy.array([45, 45, 45, 45, 60, -79, 62, -78], dtype = "int8").tostring()
#            self.unitsize = numpy.array(len(self.dataunit), dtype = "int32")
#            self.version = numpy.array([1, 2, 0, 0], dtype = "int8")
#            self.headersize = numpy.array((40 + len(self.name) + len(self.extra)), dtype = "int32")
#        else:    
#            temp = numpy.zeros(64, dtype = "int8")
#            temp[56:64] = numpy.array([45, 45, 45, 45, 60, -79, 62, -78], dtype = "int8")
#            self.extra = temp.tostring()
#            self.version = numpy.array([1, 1, 0, 0], dtype = "int8")
#            self.headersize = numpy.array((44 + len(self.name) + len(self.dataunit) + len(self.extra) + 8 * self.rank), dtype = "int32")

        if (numpy.prod(self.dim) == 1) & (self.datatype != 8):    
            extrasize = 16
        else:
            extrasize = 64
        self.headersize = numpy.array((56 + len(self.name) + len(self.dataunit) + len(self.sdef) + extrasize + 8 * self.dim.size), dtype = "int32")
        
        self.capacity_=numpy.array(self.elements()*self.datasize,dtype='int64')
            



def clearEbfMap():
    """
    Clears cached information about all files. This 
    could be used to conserve memory after a lot of different files have been 
    read. 

    >>> ebf.clearEbfMap()


    """    
    _EbfMap.ltable={}


    

class _EbfTable(object):
    """
    A class that is used to get location of data objects in an ebf file

    """
    itype=numpy.dtype([('keyloc','int64'),('keysize','int64'),('value','int64'),('next','int64'),('tnext','int64')])
    htype=numpy.dtype([('version','int8',(8,)),('endiantest','int64'),('headersize','int64'),('datatypesize','int64'),('hash_algo','int64'),('current','int64'),('htcapacity','int64'),('itemsize','int64'),('itempos','int64'),('itemcapacity','int64'),('keypos','int64'),('keycapacity','int64'),('keyposcur','int64')])
    def __init__(self,filename,mode):
        self.filename=filename
        self.mode=mode
        self.__setup(self.mode)        
        
    def close(self):
        if(self.fp1.closed == False):           
            self.fp1.close()

    def __read_hvalue(self,i):
        self.fp1.seek(self.data[2]+self.header['headersize']+self.header['datatypesize']*i,0)        
        temp=numpy.fromstring(self.fp1.read(8),dtype='int64')[0]
        if self.flagswap == 1:
            return temp.byteswap()
        else:
            return temp
        
    def __write_hvalue(self,i,hvalue):
        self.fp1.seek(self.data[2]+self.header['headersize']+self.header['datatypesize']*i,0)        
        if self.flagswap == 1:
            self.fp1.write(numpy.int64(hvalue).byteswap().tostring('C'))
        else:                                        
            self.fp1.write(numpy.int64(hvalue).tostring('C'))            
        
    def __read_node(self,i):
        self.fp1.seek(self.data[2]+self.header['itempos']+i*self.header['itemsize'],0)
        node=numpy.fromstring(self.fp1.read(self.header['itemsize']),dtype=_EbfTable.itype)
        if self.flagswap == 1:
            node=node.byteswap(True)
        return node[0]
    
    def __write_node(self,i,item):
        self.fp1.seek(self.data[2]+self.header['itempos']+i*self.header['itemsize'],0)
        if self.flagswap == 1:
            self.fp1.write(numpy.array(item).byteswap().tostring('C'))
        else:
            self.fp1.write(item.tostring('C'))
    
    def __read_key(self,item):
        self.fp1.seek(self.data[2]+self.header['keypos']+item['keyloc'],0)
        return self.fp1.read(item['keysize']).decode('ascii')
    
    def __write_key(self,item,key):
        self.fp1.seek(self.data[2]+self.header['keypos']+item['keyloc'],0)
        self.fp1.write(key.encode('ascii'))
        
    def __read_header(self):
        if (self.data[1]>0)and(self.data[2]>0)and(self.data[3]==1):
            """ read header"""
            self.fp1.seek(self.data[2],0)
            self.header=numpy.fromstring(self.fp1.read(_EbfTable.htype.itemsize),dtype=_EbfTable.htype)
            if self.header['endiantest'] != 1684234849:
                self.header=self.header.byteswap(True)
                self.flagswap=1
            self.header=self.header[0]                
            if(self.header['version'][0] != 1)or(self.header['endiantest'] != 1684234849):
                self.ecode=11
        else:
                self.ecode=12                
    
    def __write_header(self):
        self.fp1.seek(self.data[2],0)
        if self.flagswap == 1:
            self.fp1.write(numpy.array(self.header).byteswap().tostring('C'))
        else:
            self.fp1.write(self.header.tostring('C'))
    
    def getKeyValsHT(self):
        """
        get the location of the data object
        """        
        
        """ HT 120KOPS, old HT 23 KOPS, IT 9 KOPS"""
        keys=[]
        values=[]        
        self.fp1.seek(self.data[2]+self.header['headersize'],0);        
        hvalue=numpy.fromstring(self.fp1.read(8*self.header['htcapacity']),dtype='int64')
        if self.flagswap == 1:
            hvalue=hvalue.byteswap(True)
        self.fp1.seek(self.data[2]+self.header['itempos'],0)
        nodearr=numpy.fromstring(self.fp1.read(self.header['itemsize']*self.header['itemcapacity']),dtype=_EbfTable.itype)
        if self.flagswap == 1:
            nodearr=nodearr.byteswap(True)
        self.fp1.seek(self.data[2]+self.header['keypos'],0)
        keyarr=self.fp1.read(self.header['keycapacity'])                       
        if self.ecode==0:
            for i in numpy.arange(0,self.header['htcapacity']):
                loc1=hvalue[i]
                if loc1 != 0:
                    item1=nodearr[loc1]
                    keys.append(keyarr[item1['keyloc']:item1['keyloc']+item1['keysize']])
                    values.append(item1['value'])
                    while item1['next'] != -1:
                        item1=nodearr[item1['next']]
                        keys.append(keyarr[item1['keyloc']:item1['keyloc']+item1['keysize']])
                        values.append(item1['value'])                    
        return keys,values


    
    def __expand(self,keyno,keysize):
        """ 
        expand an existing hash table 
        deletes old and transfers its contents to the new expanded table.
        While expanding extra item is added (trash previous htbale) hence it is important to 
        check the correct capacity  within expand function.
        """
        if(self.ecode != 0):
            self.close()
            raise RuntimeError('EBF: error, cannot expand htable')        
        factor=self.header['keycapacity']/self.header['itemcapacity']    
        capacity=self.header['itemcapacity']
        keys1,values1=self.getKeyValsHT()
        loc=self.data[1]
        self.close()

        """ get name for trashing """
        i=0
        while _EbfTable.get(self.filename,'/.tr/.ebf/htable.t'+str(i)) >= 0:
            i=i+1
                               

        """ rename """
        fp1 = open(self.filename, 'rb+')
        fp1.seek(loc,0)        
        ebfh=_EbfHeader()
        ebfh.read(fp1)
        if(ebfh.name != '/.ebf/htable'):
            fp1.close()
            raise RuntimeError('EBF error: htable not found')                
        ebfh.rename('/.tr/.ebf/htable.t'+str(i))
        fp1.seek(loc,0)
        ebfh.write(fp1)
        fp1.close()
        keys1.append('/.tr/.ebf/htable.t'+str(i))
        values1.append(loc)
        
        temp=0
        for key1 in keys1:
            temp=temp+len(key1)                                
        if capacity <= 0:
            capacity=16
        while (keyno+len(keys1)+1) > capacity:
            capacity=capacity*2
        while (keysize+temp+1) > (capacity*factor):
            capacity=capacity*2            
        
        
        """ create """
        _EbfTable.__create(self.filename,capacity,self.flagswap)
        
        """ add key values """
        self.__setup("rb+")
        if self.ecode != 0:
            self.close()
            raise RuntimeError('EBF error: __expand,  htable is closed')
        else:
            values1[keys1.index('/.ebf/htable')]=self.data[1]            
            for key1,value1 in zip(keys1,values1):    
                self.__add(key1,value1)                
                
        self.close()                       
        self.__setup(self.mode)    
        
                        
    def __setup(self,mode):
        self.ecode=0
        self.hecode=0
        self.flagswap=0
        self.ebfh=_EbfHeader()
        self.data=numpy.zeros(4,dtype='int64')
        self.dpos=-1
        self.fp1 = open(self.filename, mode)        
        if (self.fp1.closed == False):
            self.fp1.seek(0,0)
            self.ebfh.read(self.fp1)            
            self.dpos=self.fp1.tell()
            if (self.ebfh.name == '/.ebf/info') & (self.ebfh.elements() >= 4) & (self.ebfh.datatype ==3):
                """ get htable location """
                self.data=numpy.fromstring(self.fp1.read(self.ebfh.capacity()),dtype=_TypeManager.itos_s(self.ebfh.datatype))
                if self.ebfh.flagswap == 1:
                    self.data=self.data.byteswap(True)
                self.__read_header()                                    
            else:
                self.flagswap=0    
                self.ecode=10
                self.hecode=10
        else:
            self.ecode=5
            self.hecode=5
            
#            return ebfh, dpos,numpy.zeros(3,dtype='int64'),numpy.zeros(1,dtype=__EbfTable.htype)[0],flagswap
    def __getfromfp(self,key):
        """
        get the location of the data object
        """
        location=-1
        if self.ecode == 0:                
            keyhash=_EbfTable.ebflthash(key,self.header['htcapacity'])
            loc1=self.__read_hvalue(keyhash)   
            if loc1 != 0:
                item1=self.__read_node(loc1)
                while (self.__read_key(item1) != key) and (item1['next'] != -1):
                    item1=self.__read_node(item1['next'])
                if self.__read_key(item1) == key:
                    location=item1['value']                   
        else:        
            self.fp1.seek(0,2)
            filesize=self.fp1.tell()
            self.fp1.seek(0,0)
            ebfh=_EbfHeader()
            while self.fp1.tell() < filesize:
                location1=self.fp1.tell()
                ebfh.read(self.fp1)
                if ebfh.name == key:
                    location=location1
                    break
                else:
                    self.fp1.seek(ebfh.capacity(),1)
                    
        return location
            
    def __add(self,key,value):
            self.header['current']=self.header['current']+1
            """ check for space """
            """ +1 needed as index is 1 based"""
            if (self.header['current']+1) > self.header['itemcapacity']:
                self.close()
                raise RuntimeError("not enough space for more keys")
            """ +1 not needed as index is zero based, but is just kept for simplicity"""
            if (self.header['keyposcur']+len(key)+1) > self.header['keycapacity']:
                self.close()
                raise RuntimeError("not enough space for more keys")            
            
            """ create item """
            item=numpy.zeros(1,dtype=_EbfTable.itype)[0]            
            item['keyloc']=self.header['keyposcur']
            item['keysize']=len(key)
            item['value']=value
            item['next']=-1
            item['tnext']=-1
            self.header['keyposcur']=self.header['keyposcur']+len(key)            
            
            """ write current number of items"""
            self.__write_header()
            self.__write_node(self.header['current'],item)
            self.__write_key(item,key)                        
            
            """ write to hash table or item updated with pointer"""
            keyhash=_EbfTable.ebflthash(key,self.header['htcapacity'])
            loc1=self.__read_hvalue(keyhash)
            if loc1 != 0:
                item1=self.__read_node(loc1)
                while item1['next'] != -1:
                    loc1=item1['next']
                    item1=self.__read_node(loc1)
                item1['next']=self.header['current']
                self.__write_node(loc1,item1)
            else:
                self.__write_hvalue(keyhash, self.header['current'])
                    
    @staticmethod
    def ebflthash(mystr,capacity):
# does not work with older numpy versions <1.9
#        numpy.warnings.simplefilter("ignore",RuntimeWarning)
#        numpy.warnings.simplefilter("default",RuntimeWarning)
        old_settings=numpy.seterr(over='ignore')
        ehash=numpy.uint64(5381)
        y=numpy.uint64(numpy.fromstring(mystr,dtype='int8'))    
        for i in y:
            ehash=ehash*numpy.uint64(33)+i
        numpy.seterr(**old_settings)
        return numpy.int64(ehash%numpy.uint64(capacity))    

    @staticmethod
    def ebflthash1(mystr,capacity):
        ehash=numpy.uint64(5381)
        y=numpy.uint64(numpy.fromstring(mystr,dtype='int8'))    
        for i in y:
            ehash=ehash*numpy.uint64(33)+i
        return numpy.int64(ehash%numpy.uint64(capacity))    
        

    @staticmethod
    def ebfckhash(mystr,hash1):
        old_settings=numpy.seterr(over='ignore')
        if (hash1 == 0):
            ehash=numpy.int64(5381)
        else:
            ehash=numpy.int64(hash1)
        y=numpy.int64(numpy.fromstring(mystr,dtype='int8'))    
        for i in y:
            ehash=ehash*numpy.int64(33)+i
        numpy.seterr(**old_settings)
        return numpy.int64(ehash)    
    
    @staticmethod
    def remove(filename, key):
        """
        get the location of the data object
        """        
        fileht=_EbfTable(filename,'rb+')
        if(fileht.fp1.closed == True):
            raise RuntimeError("Ebf error: unable ot open file- "+filename)
        ecode1=0
        
        if (fileht.ecode == 0):                             
            keyhash=_EbfTable.ebflthash(key,fileht.header['htcapacity'])
            loc1=fileht.__read_hvalue(keyhash)
            ecode1=1        
            if loc1 != 0:            
                item1=fileht.__read_node(loc1)
                if fileht.__read_key(item1) == key:                             
                    if item1['next'] == -1:
                        temp=numpy.int64(0)
                    else:    
                        temp=numpy.int64(item1['next'])
                    fileht.__write_hvalue(keyhash, temp)
                    ecode1=0
                else:
                    locp=-1
                    locc=loc1
                    while (fileht.__read_key(item1) != key) and (item1['next'] != -1):
                        locp=locc        
                        locc=item1['next']        
                        item1=fileht.__read_node(item1['next'])                            
                    if fileht.__read_key(item1) == key:
                        itemp=fileht.__read_node(locp)
                        itemp['next']=item1['next']
                        fileht.__write_node(locp, itemp)
                        ecode1=0            
        
        fileht.close()
        if (fileht.ecode != 0):                             
            raise RuntimeError("EBF: error in __remove(), ecode!=0")
            
        if(ecode1 != 0):
            raise RuntimeError(['EBF error: item to remove not found- '+filename+':'+key]);

        
    @staticmethod
    def __create(filename,capacity,option):
        """ create a new hash table with a given capacity """
        """ position is updated """                                       
        header=numpy.zeros(1,dtype=_EbfTable.htype)[0]
        header['version'] = numpy.zeros(8,dtype='int8')
        header['version'][0]=1
        header['endiantest']=1684234849
        header['headersize'] = _EbfTable.htype.itemsize
        header['datatypesize'] = 8
        header['hash_algo'] = 1
        header['current'] = 0
        header['htcapacity'] = capacity+capacity/2
        header['itemsize'] = _EbfTable.itype.itemsize
        header['itempos']=header['headersize']+header['datatypesize']*header['htcapacity']        
        header['itemcapacity'] = capacity
        header['keypos']=header['headersize']+header['datatypesize']*header['htcapacity']+header['itemcapacity']*header['itemsize']        
        header['keycapacity'] = 20*header['itemcapacity']
        header['keyposcur'] = 0
        
        table = numpy.zeros(header['htcapacity'],dtype='int64')
        items = numpy.zeros(header['itemcapacity'],dtype=numpy.dtype(_EbfTable.itype))
        
        fp1 = open(filename, 'rb+')
        if(fp1.closed == True):
            raise RuntimeError('EBF error: cannot open file- '+filename);
        
        ebfh=_EbfHeader()
        ebfh.read(fp1)
        if (ebfh.name != '/.ebf/info') or (ebfh.elements() < 4) or (ebfh.datatype !=3):        
            fp1.close()
            raise RuntimeError('EBF: error, /.ebf/info not found')        

        dpos_ck=fp1.tell()            
        data=numpy.fromstring(fp1.read(ebfh.capacity()),dtype=_TypeManager.itos_s(ebfh.datatype))
        if ebfh.flagswap == 1:
            data=data.byteswap(True)
        
        
        """Write htable """
        ebfh2=_EbfHeader()
        data1=numpy.zeros(header['keypos']+header['keycapacity'],dtype='int8')
        ebfh2.create('/.ebf/htable',data1,'','')
        fp1.seek(0,2)
        location=fp1.tell()
        if(option == 1):
            ebfh2.flagswap=1
        ebfh2.write(fp1)
        
        offset=fp1.tell()
        if(option == 1):
            fp1.write(numpy.array(header).byteswap().tostring('C'))
            fp1.write(table.byteswap().tostring('C'))
            fp1.write(numpy.array(items).byteswap().tostring('C'))
            x=numpy.zeros(header['keycapacity'],dtype='int8')
            fp1.write(x.byteswap().tostring('C'))        
        else:
            fp1.write(header.tostring('C'))
            fp1.write(table.tostring('C'))
            fp1.write(items.tostring('C'))
            x=numpy.zeros(header['keycapacity'],dtype='int8')
            fp1.write(x.tostring('C'))        
        offset1=fp1.tell()
            
        data[1]=location        
        data[2]=offset
        data[3]=1
        fp1.seek(dpos_ck)        
        if ebfh.flagswap == 1:
            fp1.write(data.byteswap().tostring('C'))
        else:
            fp1.write(data.tostring('C'))
                        
        fp1.close()
        
        if (offset1-offset) != ebfh2.capacity():
            fp1 = open(filename, 'rb+')
            temp=offset1-offset
            data1=numpy.zeros(1,temp,dtype='int8')
            ebfh2.setHeader('/.ebf/htable',data1,'','')
            fp1.seek(location,'bof')
            ebfh2.write(fp1)
            fp1.close()
            raise RuntimeError('Ebf error: something wrong in create()')   


                
    @staticmethod
    def init(filename):
        """create /.ebf/info and /.ebf/htable in an empty file"""
        """then add them to hash table and update checksum"""
        fp1 = open(filename, 'wb')
        if (fp1.closed == True):
            raise RuntimeError('Ebf error: cannot open file- '+filename)
        else:                                            
            ebfh = _EbfHeader()
            data=numpy.zeros(5,dtype='int64')
            ebfh.create('/.ebf/info', data, '', '')
            ebfh.write(fp1)
            fp1.write(data.tostring('C'))
            fp1.close()
            _EbfTable.__create(filename,16,0)
            keys,values=_EbfTable.getKeyValsIT(filename)
            _EbfTable.put(filename,keys,values)

    @staticmethod
    def init_swap(filename):
        """create /.ebf/info and /.ebf/htable in an empty file"""
        """then add them to hash table and update checksum"""
        fp1 = open(filename, 'wb')
        if (fp1.closed == True):
            raise RuntimeError('Ebf error: cannot open file- '+filename)
        else:                                            
            ebfh = _EbfHeader()
            data=numpy.zeros(5,dtype='int64')
            ebfh.create('/.ebf/info', data, '', '')
            ebfh.flagswap=1
            ebfh.write(fp1)
            fp1.write(data.byteswap().tostring('C'))
            fp1.close()
            _EbfTable.__create(filename,16,1)
            keys,values=_EbfTable.getKeyValsIT(filename)
            _EbfTable.put(filename,keys,values)
        
    @staticmethod
    def put(filename,key,value):
                    
        if type(key) is str:
            key=[key]
            value=[value]
        i=0    
        for key1 in key:
            key[i]=key[i].strip().lower()        
            i=i+1
        
        fileht=_EbfTable(filename,"rb+")                        
        """ check if item present """                        
        if (fileht.fp1.closed == True):
            raise RuntimeError('Ebf error: cannot open file- '+filename)
                
        if (fileht.ecode == 0):                                                                                                         
            """ check capapcity """
            for key1 in key:
                temp=fileht.__getfromfp(key1)                                        
                if temp >= 0:
                    print(filename,":",key1," already present hence exiting put")
                    fileht.ecode = 10
               
            if (fileht.ecode == 0):                                                                                                         
                temp=0
                for key1 in key:
                    temp=temp+len(key1)                                
                capacity=fileht.header['itemcapacity']
                if capacity <= 0:
                    capacity=16
                while fileht.header['current']+len(key)+1 > capacity:
                    capacity=capacity*2
                while (fileht.header['keyposcur']+temp+1) > capacity*(fileht.header['keycapacity']/fileht.header['itemcapacity'])  :
                    capacity=capacity*2            
                """ expand capacity if needed """
                if capacity != fileht.header['itemcapacity']:
                    fileht.__expand(len(key),temp)
                                    
                """ add key values """
                for key1,value1 in zip(key,value):    
                    fileht.__add(key1,value1)
                    
        if (fileht.hecode==0):
            """ compute cksum """                        
            cksum=fileht.data[0]    
            for key1,value1 in zip(key,value):
                mystr='('+key1+', '+str(value1)+')'
                cksum=numpy.int64(_EbfTable.ebfckhash(mystr,cksum))                
            """ update checksum """
            fileht.fp1.seek(fileht.dpos,0)
            fileht.data[0]=cksum
            if fileht.ebfh.flagswap == 1:
                fileht.fp1.write(fileht.data[0].byteswap().tostring('C'))
            else:
                fileht.fp1.write(fileht.data[0].tostring('C'))
                
                
        fileht.close()
        return 1
                

    @staticmethod
    def get(filename,key):
        """
        get the location of the data object
        """
        key=key.lower()        
        fileht=_EbfTable(filename,"rb")
        if (fileht.fp1.closed == True):
            raise RuntimeError('Ebf error: cannot open file- '+filename)
        location=fileht.__getfromfp(key)                                                            

        fileht.close()
        return location
            
        
            
        
    @staticmethod        
    def display(filename):
        """
        get the location of the data object
        """
        fileht=_EbfTable(filename,'rb')
        j=0
        l=0
        if fileht.ecode == 0:
            for i in numpy.arange(0,fileht.header['htcapacity']):
                loc1=fileht.__read_hvalue(i)
                if loc1 == 0:
                    j=j+1
                else:                        
                    k=1
                    item1=fileht.__read_node(loc1)
                    while item1['next'] != -1:
                        item1=fileht.__read_node(item1['next'])
                        k=k+1

                    if k > 1:
                        l=l+1
            print('total',j*1.0/(fileht.header['htcapacity']),l*1.0/(fileht.header['itemcapacity']))
        else:
            print('Ebf: error in dispaly(), ecode!=0')
                
        fileht.close()
        
        
    @staticmethod        
    def display_htab(filename):
        """
        get the location of the data object
        """
        fileht=_EbfTable(filename,'rb')
        print("{0:<24s}{1:s}".format("filename:",filename))
        print("\t {0:>24s}{1:<}".format("ecode=",fileht.ecode))
        print("\t {0:>24s}{1:<}".format("hecode=",fileht.hecode))
        print("\t {0:>24s}{1:<}".format("info_flagswap=",fileht.ebfh.flagswap))
        print("\t {0:>24s}{1:<}".format("flagswap=",fileht.flagswap))
        print("\t {0:>24s}{1:<24}{2:<24}{3:<24}".format("/.ebf/info=",fileht.data[0],fileht.data[1],fileht.data[2]))
        if fileht.ecode == 0:
            print('header:')
            print("\t {0:>24s}{1:<}".format("count=",fileht.header['current']))
            print("\t {0:>24s}{1:<}".format("htcapapcity=",fileht.header['htcapacity']))
            print("\t {0:>24s}{1:<}".format("itemcapapcity=",fileht.header['itemcapacity']))
            print("\t {0:>24s}{1:<}".format("keycapapcity=",fileht.header['keycapacity']))
            print("\t {0:>24s}{1:<}".format("keyposcur=",fileht.header['keyposcur']))
        fileht.close()
        
        [keys,values]=_EbfTable.getKeyVals(filename)
        print('Key Value List:  nsize=',len(keys))
        for key1,value1 in zip(keys,values):
            print('{0:>24s}{1:4s}{2:<}'.format(key1,' -> ',value1))
    
    @staticmethod
    def getKeyValsIT(filename):
        ebfh=_EbfHeader()
        fp1 = open(filename, 'rb')
        if (fp1.closed == True):
            raise RuntimeError('Ebf error: cannot open file- '+filename)
        fp1.seek(0,2)
        filesize=fp1.tell()
        keys=[]
        values=[]
        fp1.seek(0,0)
        while fp1.tell() < filesize:
            location=fp1.tell()
            ebfh.read(fp1)
            keys.append(ebfh.name)
            values.append(location)
            fp1.seek(ebfh.capacity(),1)
        fp1.close()
        return keys,values
    
    @staticmethod
    def getKeyVals(filename):
        fileht=_EbfTable(filename,'rb')
        if fileht.ecode == 0:
            keys,values=fileht.getKeyValsHT()
            fileht.close()
        else:
            fileht.close()
            keys,values=_EbfTable.getKeyValsIT(filename)
            
        return keys,values
                               

def keys(filename,dataname):
    keys1=[]
    if dataname[-1] == '/':
        keys=_EbfTable.getKeyVals(filename)[0]
        for i in range(len(keys)):
            if keys[i].startswith(dataname)and(keys[i].startswith('/.ebf/') == False)and(keys[i].startswith('/.tr/') == False):
                keys1.append(keys[i].split(dataname,1)[1]) 
    else:
        header=getHeader(filename,dataname)
        if header.datatype == 8:
            keys1 = numpy.dtype(sdef2descr(header.sdef)[0]).names           
    return  keys1

def rename(filename,oldkey,newkey):
    """ 
    Rename a data item in an ebf file

    Args:
        filename: string

        oldkey: a string, the name of key to rename

        newkey: a string, the new name. If new key is blank '', then  a \
        name of the form '/.tr'+oldkey+'.X' is created. Here X is a an \ 
        integer greater than equal to zero, which is incremented each \ 
        time the item with same name is deleted. 

    Example:

    >>> ebf.rename('check.ebf','/x1','/x2')
    >>> ebf.rename('check.ebf','/x1','')    

    """
    
    oldkey=oldkey.strip().lower()
    newkey=newkey.strip().lower()
    loc=_EbfTable.get(filename,oldkey)                               
    if (newkey != oldkey)and(loc >= 0):            
        if newkey == '':
            i=0
            while _EbfTable.get(filename,'/.tr'+oldkey+'.'+str(i)) != -1:
                i=i+1
                if(i > 1000000):
                    raise RuntimeError('EBF: error, too many deleted items')
            newkey='/.tr'+oldkey+'.'+str(i)

        loc1=_EbfTable.get(filename,newkey)                               
        if (loc < 0):
            raise RuntimeError('EBF error: data item/key not found')
        
        if (loc1 > 0):
            raise RuntimeError('EBF error: a key with given name already exists')
                        
        """ rename """
        fp1 = open(filename, 'rb+')
        fp1.seek(loc,0)        
        ebfh=_EbfHeader()
        ebfh.read(fp1)
        ebfh.rename(newkey)
        fp1.seek(loc,0)
        ebfh.write(fp1)
        fp1.close()
        _EbfTable.remove(filename,oldkey)
        _EbfTable.put(filename,newkey,loc)



#----------------------------------------------------------------------------------------------

def unit(filename, dataname):
    """
    Get physical units of the data type if supplied in file or else empty string

    Args:
        filename(str):

        dataname(str):

   Returns:
        str. 

   Example:

   >>> ebf.unit("check.ebf","/x")
    

    """
    location = _EbfTable.get(filename, dataname)
    if location < 0:
        raise RuntimeError("Ebf error: Data object "+dataname+" not found")
    if location >= 0:
        fp1 = open(filename, 'rb')
        header = _EbfHeader()
        fp1.seek(location, 0)
        header.read(fp1)
        fp1.close()            
        if header.datatype == 8:
            return sdef2descr(header.sdef)[1]
        else:
            return header.dataunit
    else:
        return ""


def getHeader(filename, dataname):
    """
    Get header of the data item

    Args:
        filename(str):

        dataname(str):

   Returns:
        str. 

   Example:

   >>> ebf.getHeader("check.ebf","/x")
    

    """
    location = _EbfTable.get(filename, dataname)
    if location < 0:
        raise RuntimeError("Ebf error: Data object "+dataname+" not found")
    if location >= 0:
        fp1 = open(filename, 'rb')
        header = _EbfHeader()
        fp1.seek(location, 0)
        header.read(fp1)
        fp1.close()            
        return header
    else:
        return ""
    
    

def sdef2descr(sdef):
    temp=sdef.split('\n',1)[0].strip()
    if temp == 'ver-1':
        return [__sdef2descrv1(sdef),'']
    elif temp == 'ver-2':
        return __sdef2descrv2(sdef)
    else:
        raise RuntimeError('Ebf unrecognized sdef version'+temp)

#def __descr2sdef(descr):
#    return __descr2sdefv1(descr)
def descr2sdef(descr,units=''):
#    if units == '':
#        units='NULL'
    return __descr2sdefv2(descr,units)
#    return __descr2sdefv1(descr)

#def __descr2sizev3(descr):
#    size=len(descr)
#    for temp in descr:
#        if type(temp[1]) == type([]):
#            size=size+__descr2sizev3(temp[1])
#    return size



def __descr2sdefv2(descr,units,ic=0):
    status=0
    if ic == 0:
#        taglist=['ver-2 ','<sdef>','anonymous 8 1 1 '+str(len(descr))]
        taglist=['ver-2 ','<sdef>','anonymous,8,1,1,'+str(len(descr))]
        status=1
    else:
        taglist=[]
        
    units1=[]
    if units == '':
        for temp in descr:
            units1.append('')
    else:
        units1=units
        
    for i,temp in enumerate(descr):
        shape=''
        if len(temp) > 2:
            shape=temp[2]
            if type(shape) == tuple :
                shape=list(shape)
            if type(shape) != list :
                shape=[shape]
        else:
            shape=[]
        if type(temp[1]) == type([]):
#            tag=temp[0]+' 8 '+str(len(shape))+' '+' '.join(map(str, shape))+' '+str(len(temp[1]))+' '+str(units1[i])
            tag=temp[0]+','+','.join(map(str, [8,len(shape)]+shape+[len(temp[1])]))
#            tags=tag.split()
#            tag=' '.join(tags)
            taglist.append(tag)
            taglist=taglist+__descr2sdefv2(temp[1],units1[i],ic=1)
        else:
            if temp[1][1] == 'S':
                shape.append(int(temp[1][2:]))
                type1=1
            else:
                type1=_TypeManager.stoi(temp[1][1:])
            if len(units1[i]) >0:
                units2=','+units1[i]
            else:   
                units2=''
            tag=temp[0]+','+','.join(map(str, [type1,len(shape)]+shape+[0]))+units2
#            tag=temp[0]+' '+str(type1)+' '+str(len(shape))+' '+' '.join(map(str, shape))+' 0 '+str(units1[i]
#            tags=tag.split()
#            tag=' '.join(tags)
            taglist.append(tag)
        
    if status == 1:
        taglist.append('</sdef>')
        return '\n'.join(taglist)
    else:
        return taglist


def __sdef2descrv2(sdef,begin=0,nsize=1):
    status=0        
    if begin == 0:
        sdef=sdef.split('\n')
        begin=0
        while sdef[begin].split()[0] != '<sdef>':
            begin=begin+1
        words=[x.strip() for x in sdef[begin+1].split(',')]
        rank=int(words[2])
        nsize=int(words[3+rank])
        begin=begin+2
        status=1
    dt=[]
    units=[]
    i=0
    while i < nsize:
        words=[x.strip() for x in sdef[begin].split(',')]
#        words=sdef[begin].split()
        rank=int(words[2])
        shape=[int(temp) for temp in words[3:3+rank]]
        n_fields=int(words[3+rank])
        if len(words) > 4+rank:
#            units1=words[4+rank]
            units1=','.join(words[4+rank:len(words)])
        else:
            units1='NULL'
        begin=begin+1
        if int(words[1]) == 8:
            [dtv,units1,begin]=__sdef2descrv2(sdef,begin,n_fields)
            dt.append((words[0],dtv,tuple(shape)))
            units.append(units1)
        elif int(words[1]) == 1:
            dt.append((words[0],'S'+str(shape[-1]),tuple(shape[0:-1])))            
            units.append(units1)
        else:
            dt.append((words[0],_TypeManager.itos_l(int(words[1])),tuple(shape)))
            units.append(units1)
        i=i+1
        
    if status ==1:
        if sdef[begin].split()[0] != '</sdef>':
            raise RuntimeError('Problem reading sdef')
#        return [dt[0][1],units[0]]
        return [dt,units]
    else:
        return [dt,units,begin]

#
#def __descr2sdefv3(descr,units,ic=0):
#    status=0
#    if ic == 0:
#        taglist=['ver-3 ','anonymous 8 1 1 '+str(len(descr))+' 2 ']
#        ic=2
#        status=1
#    else:
#        taglist=[]
#        
#    units1=[]
#    if units == 'NULL':
#        for temp in descr:
#            units1.append('NULL')
#    else:
#        units1=units
#        
#    ic=ic+len(descr)
#    mylist=[]
#    for i,temp in enumerate(descr):
#        shape=''
#        if len(temp) > 2:
#            shape=temp[2]
#            if type(shape) == tuple :
#                shape=list(shape)
#            if type(shape) != list :
#                shape=[shape]
#        else:
#            shape=[]
#        if type(temp[1]) == type([]):
#            tag=temp[0]+' 8 '+str(len(shape))+' '+' '.join(map(str, shape))+' '+str(len(temp[1]))+' '+str(ic)+' '+str(units1[i])
#            mylist.append([temp[1],ic,units1[i]])
#            ic=ic+__descr2sizev3(temp[1])
#        else:
#            if temp[1][1] == 'S':
#                shape.append(int(temp[1][2:]))
#                type1=1
#            else:
#                type1=_TypeManager.stoi(temp[1][1:])
#            tag=temp[0]+' '+str(type1)+' '+str(len(shape))+' '+' '.join(map(str, shape))+' 0 '+str(ic)+' '+str(units1[i])
#        taglist.append(tag)
#        
#    for temp in mylist:
#        taglist=taglist+__descr2sdefv3(temp[0],temp[2],ic=temp[1])
#    if status == 1:
#        return '\n'.join(taglist)
#    else:
#        return taglist
#        
#def __sdef2descrv3(sdef,begin=1,nsize=1):        
#    if begin == 1:
#        sdef=sdef.split('\n')
#    dt=[]
#    units=[]
#    for tag in sdef[begin:begin+nsize]:
#        words=tag.split()
#        rank=int(words[2])
#        shape=[int(temp) for temp in words[3:3+rank]]
#        n_fields=int(words[3+rank])
#        fields=int(words[3+rank+1])
#        if len(words) > 5+rank:
#            units1=words[5+rank]
#        else:
#            units1='NULL'
#        if int(words[1]) == 8:
#            [dtv,units1]=__sdef2descrv3(sdef,fields,n_fields)
#            dt.append((words[0],dtv,tuple(shape)))
#            units.append(units1)
#        elif int(words[1]) == 1:
#            dt.append((words[0],'S'+str(shape[-1]),tuple(shape[0:-1])))            
#            units.append(units1)
#        else:
#            dt.append((words[0],_TypeManager.itos_l(int(words[1])),tuple(shape)))
#            units.append(units1)
#    if begin ==1:
#        return [dt[0][1],units[0]]
#    else:
#        return [dt,units]
    
def __descr2sdefv1(descr,name='anonymous',dshape=()):
    if name == 'anonymous':
        mystr='ver-1 \n'+'struct {\n'
        dshape=(1,)
    else:    
        mystr = 'struct {\n'
    for tag1 in descr:
        if(type(tag1[1]) == type([])):
            if len(tag1) >2 :
                tagcur = __descr2sdefv1(tag1[1],tag1[0],tag1[2])
            else:
                tagcur = __descr2sdefv1(tag1[1],tag1[0])
            #+' '+tag1[0]+' '+str(len(tag1[2]))+' '+str(tag1[2])+' ;\n'
        else:
            
            if (tag1[1][1] == 'S'):
                datatype = 'char'
            else:
                datatype = _TypeManager.stos_l(str(tag1[1][1:]))
            
            if (datatype == 'char') & (int(tag1[1][2:]) > 1):
                if(len(tag1) > 2):
                    tagcur = datatype+' '+str(tag1[0])+' '+str(len(tag1[2])+1)+' '+str(tag1[2])+tag1[1][2:]+' ;\n'
                else:
                    tagcur = datatype+' '+str(tag1[0])+' 1 '+tag1[1][2:]+' ;\n'
            else:
                if(len(tag1) > 2):
                    tagcur = datatype+' '+str(tag1[0])+' '+str(len(tag1[2]))+' '+str(tag1[2])+' ;\n'
                else:
                    tagcur = datatype+' '+str(tag1[0])+' 0  ;\n'
                        
        tagcur = tagcur.replace('(', ' ')
        tagcur = tagcur.replace(')', ' ')
        tagcur = tagcur.replace(',', ' ')
        mystr += tagcur
        
        tagcur=str(dshape)
        tagcur = tagcur.replace('(', ' ')
        tagcur = tagcur.replace(')', ' ')
        tagcur = tagcur.replace(',', ' ')

    mystr += '} '+ name +' '+str(len(dshape))+ tagcur+' ; \n'            
    return mystr

def __sdef2descrv1(wl, ic=0):
    if ic == 0:
        wl=wl.split()
    dth = []
    i=ic
    while wl[i] != '{':
        i=i+1
        
    i=i+1
    while i < len(wl):
        if wl[i] == '}':
            break
        datatype=_TypeManager.stos_s(wl[i])
        if datatype == 'V':
            datatype=__sdef2descrv1(wl,i)
            i=i+1
            l=1
            while l != 0:
                i=i+1
                if wl[i] == '{':
                    l=l+1
                if wl[i] == '}':
                    l=l-1                
                    
        name=wl[i+1]
        rank=int(wl[i+2])
        dims=[]
        for j in range(0,rank):
            dims.append(int(wl[i+3+j]))
        i=i+3+rank
        if (datatype == 'S1'):
            if dims[rank-1] > 1:
                temp = dims.pop(rank-1)
                datatype = 'S'+str(temp)
                rank=rank-1
        if rank > 0:        
            dth.append((name,datatype, tuple(dims)))
        else:
            dth.append((name,datatype))
                
        if wl[i] == ';':
            i=i+1
        else:
            raise RuntimeError("EBF: missing ; in sdef")
        
            
    return numpy.dtype(dth)

def containsKey(filename,dataname):
    """    
    Check if a data item is present in an ebf  file.     
    
    Args:

        filename : a string specifying filename

        dataname : name of data item
    
    Returns:
        1 if an item is present else 0
    
    Example:

    >>> ebf.containsKey('check.ebf','/x')
    
    
    """
    if(_EbfTable.get(filename, dataname) < 0):
        return 0
    else:
        return 1

                    
def read(filename, path = '/' ,recon=0,ckon=1,begin=0,end=None):
    """
    Read data from an ebf file

    Args:

         filename(str) :

         path(str)     : tagname of data to be read from the ebf file \
         or a path to the data items within the file. If ending with +\
         then all arrays in the same path having same size as the 
         specfied array are read. Useful to load tables where 
         individual columns are written separately.  

         recon(integer): Should be 1 if one wants to load data \
         objects recursively. Should be 0 if one wants to load \
         data objects only under current path.Defualt is 0.

         ckon  : option that determines if checksum is to be compared with \
         checksum on file. Default is to compare, but if there is \
         little possibility of file being externally modified then it can 
         be set to 0.  

    Returns:
         numpy.ndarray or a dictionary of numpy.ndarray. If multiple \
         data items are to be read as a dictionary, the path must end \
         with '/' in the later case.    

    """

    path = path.strip()
    mydict = {}
    x = ''
    rows=-1
    if path.endswith('+'):
        path=path.rstrip('+')
        if path.endswith('/'):
            raise RuntimeError('EBF Error: in read()- tagname ending with /+')
        temp=getHeader(filename,path).getshape()
        if len(temp) == 0:
            temp=numpy.array([1],dtype='int64')
        rows=temp[0]
        temp=path.rpartition('/')                        
        path=temp[0]+temp[1]
        
        


    """ Needed here to make sure the map is loaded """    
    if (path.endswith('/') == 0):
        location=_EbfTable.get(filename, path.lower())
        if location >= 0:
            fp1 = open(filename, 'rb')
            fp1.seek(location, 0)        
    
            header = _EbfHeader()
            header.read(fp1)
            if header.datatype == 8:
#                dth1 = numpy.dtype(__sdef2descr(header.sdef.split()))
                dth1 = numpy.dtype(sdef2descr(header.sdef)[0])
            else:
                dth1 = _TypeManager.itos_s(header.datatype)
            
            """for character array change the dim to convert the last dim to string """            
            if (header.dim.size > 0):
                if (end is None):
                    end1=header.dim[0]
                else:
                    end1=int(end)
                if (end1 < 0):
                    end1=header.dim[0]+end1
                if end1 > (header.dim[0]):
                    end1=header.dim[0]
                    
                begin1=int(begin)
                if begin1 > (header.dim[0]-1):
                    begin1=header.dim[0]-1
                # note for header.dim[0] this can be negative so do as below
                if begin1 < 0:
                    begin1=0
                if end1 < begin1:
                    print('ebf Warning, begin>end')
                    end1=begin1
                if begin1 > 0:
                    fp1.seek(begin1*header.datasize*header.elements()//header.dim[0],1) 
                if (end1-begin1) != header.dim[0]:
                    header.dim[0]=end1-begin1
                    
            block_size = header.elements()*header.datasize            
            
            if header.datatype == 7:
                dth1 = dth1+str(header.datasize)
            if header.datatype == 1:
                dth1 = 'S'+str(header.dim[-1])
                header.datasize=header.datasize*header.dim[-1]
                if header.dim.size == 1:
                    header.dim[0] = 1
                else:                
                    header.dim = header.dim[0:len(header.dim)-1]
                
                
            x = numpy.fromstring(fp1.read(block_size), dtype = dth1)
            if header.flagswap == 1:
                x = x.byteswap(True)
            x = x.reshape(header.getshape())
            fp1.close()
            return x
    if (path.endswith('/') == 1):
        location=_EbfMap.get(filename, path.lower(),ckon)
        node=_EbfUtils.searchPathTree(_EbfMap.ltable[filename]['pathtree'],path.lower())
        if (node['name'] == path.lower()):
            if len(node['files']) > 0:
                for key in node['files']:
                    if rows > -1:
                        temp=getHeader(filename,node['name']+key).getshape()
                        if len(temp) == 0:
                            temp=numpy.array([1],dtype='int64')
                        if temp[0] == rows:
                            mydict[key] = read(filename,node['name']+key,recon,0,begin,end)                            
                    else:
                        mydict[key] = read(filename,node['name']+key,recon,0)
            if (recon > 0)&(len(node['dirs']) > 0):
                if recon == 2:
                    for key in list(node['dirs'].keys()):
                        mydict[key.strip('/')] = read(filename,node['name']+key,recon,0)                                    
                if recon == 1:
                    for key in list(node['dirs'].keys()):
                        if (key.startswith('.ebf/') == False)and(key.startswith('.tr/') == False):
                            mydict[key.strip('/')] = read(filename,node['name']+key,recon,0)                                    
        
        
    if len(mydict) == 0:
        print(filename+":"+path)
        raise RuntimeError("Ebf error from read(): requested object not found- ")

    return mydict

    
def write(filename, tagname, data, mode, dataunit = ""):
    """
    Write data to a file

    Args:
        filename(str):

        tagname(str) : the name of data to be written to the ebf file or \
        a path ending with '/' if multiple items are to be written

        data(numpy.ndarray)    : data to be to be written

        mode(str)    : writing mode, "w" to write a fresh file or "a" \
        to append an existing file

    Kwargs:
        dataunit(str): units of data default is a blank string


    """
    if (mode == 'w')|(mode == 'wb'):
        _EbfTable.init(filename)
        mode1='ab'
        mode='ab'
    elif (mode == 'a')|(mode == 'ab'):
        mode1='ab'
        mode='ab'
    elif (mode == 'u'):
        mode1='rb+'    
    elif (mode == 'e'):
        mode1='rb+'    
    else:
        raise RuntimeError("mode must be 'w', 'a' 'u' or 'e' ")
        
         
        

    """Due to dict has to check numpy.void or else could have tested data.dtype.char """
    if (type(data) is not numpy.ndarray)&(type(data) is not numpy.void)&(type(data) is not dict):
        data=numpy.array(data)
    if mode == 'e':
        if data.ndim == 0:
            data=data.reshape(1)
        
#        raise RuntimeError('Data to be written must be of type nummpy.ndarray or numpy.void or dict')
    tagname = tagname.strip()
    if tagname.endswith('/'):
        if type(data) is dict:
            
            mykeys=list(data.keys())
            if type(dataunit) == str:
                dataunitl=[]
                for name in mykeys:
                    dataunitl.append(dataunit)
            else:
                dataunitl=dataunit
            if type(dataunitl) != list:
                raise RuntimeError('Ebf Error: dataunit must be a list')
            if len(dataunitl) != len(mykeys):
                raise RuntimeError('Ebf Error: length of dataunit list must match length of data dict')            
            i=0    
            
            for name in list(data.keys()):
                if (type(data[name]) is dict)|(type(data[name]) is numpy.void):
                    write(filename, tagname+name+'/', data[name], mode, dataunitl[i])
                else:
                    write(filename, tagname+name, data[name], mode, dataunitl[i])
                i=i+1
                
        elif (data.dtype.char == 'V'):
            if(data.size >= 1):
                data1 = data
                if data1.size == 1: 
                    data1 = numpy.squeeze(data1)

                mykeys=data1.dtype.names    
                if type(dataunit) == str:
                    dataunitl=[]
                    for name in mykeys:
                        dataunitl.append(dataunit)
                else:
                    dataunitl=dataunit
                if type(dataunitl) != list:
                    raise RuntimeError('Ebf Error: dataunit must be a list')
                if len(dataunitl) != len(mykeys):
                    raise RuntimeError('Ebf Error: length of dataunit list must match length of data dict')                                        
                i=0   
                
                for name in data1.dtype.names:
                    if data1[name].dtype.char == 'V':
                        write(filename, tagname+name+'/', data1[name], mode, dataunitl[i])
                    else:
                        write(filename, tagname+name, data1[name], mode, dataunitl[i])
                    i=i+1
            else:
                raise RuntimeError('size of ndarray must be at least one')

        else: 
            print(type(data), tagname)
            raise RuntimeError('with path ending with /, data.dtype.char should be V i.e., structure. Here '+data.dtype.char)
                   
                                         
    else:
        location=_EbfTable.get(filename, tagname.lower())
        header = _EbfHeader()
        if data.dtype.char == 'V':        
#            header.create(tagname, data, dataunit, "ver-1 \n"+__descr2sdef(data.dtype.descr,dshape=data.shape))
            header.create(tagname, data,'',descr2sdef(data.dtype.descr, dataunit))
        else:    
            header.create(tagname, data, dataunit, "")

        fp1  =  open(filename, mode1)
        if mode == 'u' :
            if location >= 0 :
                fp1.seek(location, 0)            
                header1 = _EbfHeader()
                header1.read(fp1)
                fp1.seek(location, 0)
                if (header1.get_dtype() != header.get_dtype()):
                    fp1.close()
                    raise RuntimeError('Data to be updated not present '+tagname)                                                     
                            
                if (header1.datatype != header.datatype)|(header1.datasize != header.datasize)|(header.elements() != header1.elements()):
                    fp1.close()
                    raise RuntimeError('Data to be updated not present '+tagname)                                                     
                header=header1
                if header.flagswap == 1:
                    header.flagswap=0
            else:
                fp1.close()
                raise RuntimeError('Data to be updated not present '+tagname)
            header.write(fp1)        
        elif mode == 'e' :
            if location >= 0 :
                fp1.seek(location, 0)            
                header1 = _EbfHeader()
                header1.read(fp1)
                fp1.seek(location, 0)            
                if (header1.get_dtype() != header.get_dtype()):
                    fp1.close()
                    raise RuntimeError('Data to be updated not present '+tagname)                                                     
                if (header1.datatype != header.datatype)|(header1.datasize != header.datasize)|(header1.dim.size != header.dim.size):
                    fp1.close()
                    
                    raise RuntimeError('Data to be updated not present or of mismatch size'+tagname)                                                     
                if header1.flagswap == 1:
                    fp1.close()
                    raise RuntimeError('Data is of different endian format '+tagname)                                                     
                fp1.seek(0, 2)
                locend=fp1.tell()            
                fp1.seek(location, 0)            
                if location+header1.headersize+header1.capacity() != locend:
                    fp1.close()
                    raise RuntimeError('Cannot update as not last item '+tagname)
                if header1.dim.size > 1:
                    if numpy.prod(header1.dim[1:]) != numpy.prod(header.dim[1:]) :
                        fp1.close()
                        raise RuntimeError('Cannot update as rank do not match '+tagname)
                dataend=location+header1.headersize+header1.elements()*header1.datasize
                header1.dim[0]=header.dim[0]+header1.dim[0]
                if header1.capacity_<(header1.elements()*header1.datasize):
                    header1.capacity_=(header1.elements()*header1.datasize)
                header=header1
            else:
                fp1.close()
                raise RuntimeError('Data to be updated not present '+tagname)                            
            header.write(fp1)        
            fp1.seek(dataend, 0)            
        else:
            if location >= 0 :
                fp1.close()
                raise RuntimeError('Data to be written already present '+tagname)
            location=fp1.tell()
            header.write(fp1)        
            
        fp1.write(data.tostring('C'))
        fp1.close()
        if (mode1 == 'ab'):
            _EbfTable.put(filename,tagname,location)

def join(files,path,outfile,outpath,mode):            
    data0=read(files[0],'/')
    if mode == 'w':
        initialize(outfile)
        mode='a'
    for key in list(data0.keys()):
        if containsKey(outfile,outpath+key) == 0:
            dataunit=unit(files[0],path+key)
            efile=EbfFile(outfile,outpath+key,mode,dataunit)
            for file1 in files:
                data=read(file1,path+key)
                efile.write(data)
            efile.close()
        else:
            print(("item="+outpath+key+" already present. Hence, skipping"))



def dict2npstruct(data,basekey=None,keylist=None):    
    """
    Convert a python dict containing numpy arrays to numpy struct

    Args:
        data        :

        basekey(str): Only those items in dict whose size match that of data[bsekey] will 
        be used.

        keylist(str): list of keys to beused when constructing npstruct

    """
    if keylist is None:
        keylist=list(data.keys())
    if basekey is None:
        nsize=data[keylist[0]].size
    else:
        nsize=data[basekey].size        

    dt=[]
    for key in keylist:
        if data[key].size == nsize:
            dt.append((key,data[key].dtype))

    data1=None
    if len(dt)>0:
        data1=numpy.zeros(nsize,dtype=dt)
        for key in data1.dtype.names:
            data1[key.lower()]=data[key]
    return data1



def npstruct2dict(data):
    """
    Convert an array of numpy struct to a python dict of numpy arrays

    Args:
        data        :
    """
    data1={}
    for x in data.dtype.names:
        data1[x.lower()]=data[x]
    return data1



#def islast(filename,tagname):
#    location=_EbfTable.get(filename, tagname.lower())
#    fp1  =  open(filename, mode1)
#    fp1.seek(location, 0)            
#    header1 = _EbfHeader()
#    header1.read(fp1)
#    fp1.seek(location, 0)            
#    fp1.seek(0, 2)
#    locend=fp1.tell()            
#    fp1.close()
#    if location+header1.headersize+header1.capacity() == locend:
#        return True
#    else:
#        return False


def read_ind(filename,tagname,ind):
    """
    read data from specified locations in a file 

    Args:
        filename(str):

        tagname(str) : the name of data to be read 

        ind(str)    : list or array of indices to be read

    """
    ind=numpy.array(ind,dtype='int64')        
    if ind.ndim==0:
        ind=numpy.array([ind],dtype='int64')
        efile=EbfFile(filename,tagname,'r',cache=min(1000,1))
#        data=efile.read_ind(numpy.array([ind]))[0]
        data=efile.read_ind(ind)[0]
    else:
        if ind.ndim >1:
            raise RuntimeError('ind must be 1 dimensional or scalar')
        efile=EbfFile(filename,tagname,'r',cache=min(1000,len(ind)))
        data=efile.read_ind(numpy.array(ind))        
    efile.close()
    return data

def update_ind(filename,dataname,data,ind=None):
    """
    Update existing data array in a file at user given index positions. 

    Args:
        filename(str):

        dataname(str) : the name of data to be upated 

        data    : data to be updated
        
        ind    : indices of the array on file that 
        needs to be  updated.

    """
        

    if sys.byteorder == 'little':
        sorder='<'
    else:
        sorder='>'

    location=_EbfTable.get(filename, dataname.lower())
    if location >= 0:
#        with open(filename,'rb+') as fp1:
        fp1=open(filename,'rb+')
        try:
            fp1.seek(location, 0)            
            header = _EbfHeader()
            header.read(fp1)
            datalocation=fp1.tell()
        
            if header.datatype == 8:
                dth1 = numpy.dtype(sdef2descr(header.sdef)[0])
#                if type(data) != numpy.void:
#                    data=numpy.array([data],dtype=data.dtype)
#                if type(data) != numpy.ndarray:
#                    raise RuntimeError('EbfError: data must be numpy array of void')
                if (dth1.names != data.dtype.names):
                    raise RuntimeError('EbfError: data.dtype.names do not match info on file')
            else:
                dth1 = numpy.dtype(_TypeManager.itos_s(header.datatype))
                if header.datatype == 7:
                    dth1 = numpy.dtype('S'+str(header.datasize))
                if header.datatype == 1:
                    dth1 = numpy.dtype('S'+str(header.dim[-1]))
                    header.datasize=header.datasize*header.dim[-1]
                    if header.dim.size == 1:
                        header.dim[0] = 1
                    else:                
                        header.dim = header.dim[0:len(header.dim)-1]
                
            shape = header.getshape()
            rest=1
            nsize=header.elements()
            if (len(shape)>1)and(nsize>0):
                rest=header.elements()/shape[0]
                nsize=shape[0]
            
                
            if header.flagswap==1:
                if sorder=='<':
                    sorder='>'
                elif sorder == '>': 
                    sorder='<'
                                    
            data=numpy.array(data,dtype=dth1)
            if data.ndim ==0:
                data=numpy.array([data],dtype=data.dtype)
                    
            if ind is None:
                ind=numpy.arange(shape[0])
                allset=True
            else:
                allset=False
                ind=numpy.array(ind,dtype='int64')
                if ind.ndim ==0:
                    ind=numpy.array([ind],dtype='int64')
                if ind.ndim >1:
                    raise RuntimeError('ind must be 1 dimensional or scalar')
                    
            if numpy.max(ind) >= nsize:
                raise RuntimeError('EbfError: index supplied is out of bound with data on file')
            if ind.size*rest != data.size:
                print(ind.size,data.size,rest)
                raise RuntimeError('EbfError: size of data not equal to size of ind')
            if rest != data[0].size:
                raise RuntimeError('EbfError: shape of data does not match')
            
            dorder=_EbfUtils.get_byteorder(data)
            if (sorder!=dorder):
                data=data.byteswap()

            fp1.seek(datalocation, 0)            
            if allset:
                fp1.write(data.tostring('C'))
            else:
                icur=0    
                inda=numpy.argsort(ind)
                for i in inda:
                    if ind[i] != icur: 
                        fp1.seek(datalocation+ind[i]*header.datasize*data[0].size, 0)
                        icur=ind[i]
                    # to handle strings [i:i+1] needed instead of [i]
                    fp1.write(data[i:i+1].tostring('C'))
                    icur+=1
                
        finally:
            fp1.close()
    else:
        raise RuntimeError('EbfError: data not found in file')



def iterate(filename,tagname,cache):
    """
    An iterator to read in data, part by part of a given size.
    Useful for reading big arrays which are difficult to fit in RAM.

    Args:
        filename(str):

        tagname(str) : the name of data to be read.
        Multiple items of same size can be read by appending a + sign 
        
        

        cache(int)    : no of data items to read at a time
        
    Example:
    
    >>> temp=0.0
    >>> for x in ebf.iterate('check.ebf','/x',1000):
    >>>     temp=temp+np.sum(x)   
    
    To read all items whose size match with size of "/x"
    

    >>> temp=0.0
    >>> for data in ebf.iterate('check.ebf','/x+',1000):
    >>>     temp=temp+np.sum(data['/x'])   
    

    """
    header=getHeader(filename,tagname.rstrip('+'))
    begin=0
    end=cache
    while begin < header.dim[0]:
        data=read(filename,tagname,begin=begin,end=end)
        yield data
        begin=end
        end =end+cache
        if end > header.dim[0]:
            end=header.dim[0]


class EbfFile():
    def __init__(self, filename,path,mode,dataunit='',cache=100):
        self.filename=filename
        self.path=path
        self.fp=None
        self.cache=cache
        self.begin=-1
        self.end=0
        self.mode=mode
        self.defined=False
        self.dataunit=dataunit
        if self.mode == 'w':
            _EbfTable.init(filename)
            self.mode='a'
        if (self.mode != 'a')&(self.mode != 'r'):
            print('mode=',self.mode)
            raise RuntimeError("mode must be 'r' , 'w' or 'a' ")

            
        if (self.path.endswith('/') == 0):
            self.location=_EbfTable.get(self.filename, self.path.lower())
            if self.mode == 'r':
                self._read_init()
            else: 
                self._write_init()
                
    def _read_init(self):
        if self.location >= 0:
            self.fp = open(self.filename, 'rb')
            self.fp.seek(self.location, 0)            
            self.header = _EbfHeader()
            self.header.read(self.fp)
            self.datalocation = self.fp.tell()
            if self.header.datatype == 8:
#                self.dtype = numpy.dtype(__sdef2descr(self.header.sdef.split()))
                [dt,units]=sdef2descr(self.header.sdef)
                self.units=units
                self.dtype = numpy.dtype(dt)
            else:
                self.dtype = _TypeManager.itos_s(self.header.datatype)
                self.units=self.header.dataunit
        
            """for character array change the dim to convert the last dim to string """
            if self.header.datatype == 7:
                self.dtype = self.dtype + str(self.header.datasize)
            if self.header.datatype == 1:
                self.dtype = 'S' + str(self.header.dim[-1])
                self.header.datasize = self.header.dim[-1]
                self.header.datatype = 7
                if self.header.dim.size == 1:
                    self.header.dim[0] = 1
                else:                
                    self.header.dim = self.header.dim[0:len(self.header.dim) - 1]
                    
            self.shape = list(self.header.getshape())
            self.rank = self.header.dim.size                        
#            self.rest = self.header.elements() / self.shape[0]
#            self.elements = self.shape[0]
            self.rest = 1
            self.elements = self.header.elements()
            if self.rank>1:
                self.rest = self.header.elements() / self.shape[0]
                self.elements = self.shape[0]
                
            self.datasize = self.header.datasize * self.rest
        
    def _write_init(self):
        if self.location < 0:
            self.fp = open(self.filename, 'rb+')
            self.fp.seek(0, 2)
            self.location=self.fp.tell()            
            
        
    def read(self,i,nsize=1):
        if hasattr(i,'__len__'):
            raise RuntimeError('must be scalar')
        if ((i+nsize)>self.end) or (i<self.begin):
            self.begin=i
            if self.cache < nsize:
                self.end=i+nsize
            else:
                self.end=i+self.cache
            if self.end> self.elements:
                self.end=self.elements
            if self.begin>= self.end:
                self.begin=self.end-1
            self.fp.seek(self.datalocation+self.begin*self.datasize, 0)
            self.x = numpy.fromstring(self.fp.read((self.end-self.begin)*self.datasize), dtype = self.dtype)
            if self.header.flagswap == 1:
                self.x = self.x.byteswap(True)

            if self.rank > 1:
                self.shape[0]=(self.end-self.begin)
                self.x=self.x.reshape(self.shape)

        if (i+nsize)<=self.elements:
            if nsize>1:
                return self.x[i-self.begin:i-self.begin+nsize].copy()
            else:
                return self.x[i-self.begin].copy()
        else:
            return None

    def read_ind(self, ind):
        # This method looks for groups of contiguous indices in 'ind' and
        # loads those blocks of memory with a single copy command.
        # This method is ~50% faster than looping over each index and running:
        #   data[i] = self.read(ind[i])
        if numpy.max(ind) < self.elements:
            # Find the ascending order of indices
            ind1 = numpy.argsort(ind)
            # Create an empty data array for copying the data into
            data = numpy.zeros(len(ind),dtype=self.dtype)
            # Begin at the start of the data array
            begin_data = 0
            # groupy with this lambda returns a list of generators, with 
            # each generator containing a contiguous 
            # block of indices from 'ind'
            for k, g in groupby(enumerate(ind[ind1]),
                                lambda ix : ix[0] - ix[1]):
                # converts the generator into a list 
                ind_grp = list(map(itemgetter(1), g))
                # first index of the contiguous block
                begin = ind_grp[0]
                # number of elements in the contiguous block
                nsize = len(ind_grp)
                # read from memory all the elements in the contiguous block
                d = self.read(begin, nsize=nsize)
                # final the final index into data of the contiguous block
                end_data = begin_data + nsize
                # place the extracted data into the data array
                data[begin_data:end_data] = d
                # Increment the place in the data array
                # by the size of the contiguous block
                begin_data += nsize
            return data[ind1]
        else:
            return None
            
    def write(self,data):
        if len(data) > 0:    
            if self.fp != None:
                if self.defined == False:
                    self.defined=True
                    self.header = _EbfHeader()
                    self.datatype = data.dtype
                    if data.dtype.char == 'V':
                        self.header.create(self.path, data,'',descr2sdef(data.dtype.descr,self.dataunit))
                    else:    
                        self.header.create(self.path, data, self.dataunit, "")
                    self.header.write(self.fp)
                if data.dtype != self.datatype:
                    try:
                        temp=numpy.array(data,dtype=self.datatype)
                        self.fp.write(temp.tostring('C'))
                    except:
                        self.close()
                        raise RuntimeError("EbfFile.write() error: Cannot convert types")
                else:
                    self.fp.write(data.tostring('C'))
            else:
                raise RuntimeError("EbfFile.write() error: file is closed")
            
            

    def close(self):
        if self.fp != None:
            if ((self.mode == 'w')|(self.mode == 'a'))and(self.fp.tell()>self.location):
                temp=self.header.elements()/self.header.dim[0]
                datawritten=self.fp.tell()-(self.location+self.header.headersize)
                bufsize=datawritten%(temp*self.header.datasize)
                if bufsize > 0:     
                    x=numpy.zeros(bufsize,dtype='int8')
                    self.fp.write(x.tostring('C'))
                    datawritten=datawritten+bufsize
                self.header.dim[0]=datawritten/(temp*self.header.datasize)
                if datawritten == 0:
                    self.header.dim=numpy.zeros(1, dtype = "int64")                            
                if datawritten > self.header.capacity_:
                    self.header.capacity_=datawritten
                self.fp.seek(self.location,0)
                self.header.write(self.fp)        

                self.fp.close()
                self.fp=None            
                _EbfTable.put(self.filename,self.path,self.location)
            else:
                self.fp.close()
                self.fp=None            
            self.mode=None
            self.filename=None
            self.location=None
            self.header=None
            self.path=None
            self.defined=False
            
    def __del__(self):
        self.close()

def initialize(filename):
    """
    Initialize a file for writing with mode='w'.
    After this one can use mode='a' to write rest of the items.
        
    Args:
         filename(str):
         
    Example:

    >>> ebf.initialize('check.ebf')
    >>> ebf.write('check.ebf','/x',[0,1,2],'a')
    >>> ebf.write('check.ebf','/y',[0,1,2],'a')
    is same as
    >>> ebf.write('check.ebf','/x',[0,1,2],'w')
    >>> ebf.write('check.ebf','/y',[0,1,2],'a')
    
    """
    
    _EbfTable.init(filename)

def info(filename,option=0):
    """
    Get summary of the contents of a file

    Args:
         filename(str):
         
    Kwargs:     

    Example:

    >>> ebf.info('check.ebf')

    """

    fp1 = open(filename, 'rb')
    fp1.seek(0,2)
    filesize = fp1.tell()
    fp1.seek(0,0)
    print(filename, filesize, 'bytes ') 
    print('------------------------------------------------------------------')
    print("{0:30s} {1:8s} {2:7s} {3:10s} {4:10s}".format('name', 'dtype', 'endian', 'unit', 'dim'))
    print('------------------------------------------------------------------')
    header = _EbfHeader()
    while fp1.tell() < filesize:
        header.read(fp1)
        en = sys.byteorder
        if header.flagswap == 1:
            if en == 'little':
                en = 'big' 
            else:
                en = 'little' 
    
        print("{0:30s} {1:8s} {2:7s} {3:10s} {4:10s}".format(header.name, _TypeManager.itos_l(header.datatype), en, header.dataunit, str(header.dim)))
        fp1.seek(header.capacity(), 1)
            
        if (option == 1) and (header.datatype == 8):
            print("structure definition:")    
            print(header.sdef)    
        
    if fp1.tell() != filesize:
        raise RuntimeError('EBFCorrupt')    
    else:
        fp1.close()
        
def check(filename):
    """
    check if the file is not corrupted

    Args:
         filename(str):
         
    Kwargs:     

    Example:

    >>> ebf.check('check.ebf')

    """

    fp1 = open(filename, 'rb')
    fp1.seek(0,2)
    filesize = fp1.tell()
    fp1.seek(0,0)
    header = _EbfHeader()
    ecode=0
    while fp1.tell() < filesize:
        location=fp1.tell()
        header.read(fp1)
        if(location != _EbfTable.get(filename, header.name)):
            ecode=1
            break
        fp1.seek(header.capacity(), 1)
        if(header.datasize*header.elements() > header.capacity()):
            ecode=1
            break
                    
    if fp1.tell() != filesize:
        ecode=2
    fp1.close()
    return ecode;

    

        
def stat(filename, tagname,recon=0):
    """
    Get statistics of a data item

    Args:
         filename(str):
         
         tagname(str):
         
    Kwargs:     

    Example:

    >>> ebf.stat('check.ebf','/x /y ')

    """
        
    tagname=tagname.lower()
    keysin=tagname.split()
    for key in keysin:
        if key.endswith('/'):
            location=_EbfMap.get(filename, key,1)
            break
    keys=[]    
    for key in keysin:
        if key.endswith('/'):
            if recon == 0:
                nodef=_EbfUtils.searchPathTree(_EbfMap.ltable[filename]['pathtree'],key)
                for key1 in nodef['files']:
                    keys.append(nodef['name']+key1)
            else:
                keys=keys+_EbfUtils.getKeysRecursive(_EbfUtils.searchPathTree(_EbfMap.ltable[filename]['pathtree'],key))
            
        elif containsKey(filename,key) == 1:
            keys.append(key)
        else:            
            raise RuntimeError('EBF Error: in stat(), key not present in input file')
            
                    
    print("{0:15s} {1:>10s} {2:>12s} {3:>12s} {4:>12s} {5:>12s}".format("name","items", "min", "max", "mean", "stddev"))
    for dataname in keys:
        data=read(filename,dataname)
        if data.dtype.type != numpy.string_:
            data=numpy.float64(read(filename,dataname))
            print("{0:15s} {1:10d} {2:12.4f} {3:12.4f} {4:12.4f} {5:12.4f}".format(dataname, data.size, (numpy.min(data)), (numpy.max(data)), numpy.mean(data), numpy.std(data)))
    
def cat(filename, tagname,delimiter=' ',tableon=0):
    """
    print data items in ascii format

    Args:
         filename(str):
         
         tagname(str):
         
    Kwargs:
         delimiter(str) - ' ' or ', ' for csv

    Example:

    >>> ebf.cat('check.ebf','/x /y',', ')
    >>> ebf.cat('check.ebf','/x+',', ')
    >>> ebf.cat('check.ebf','/x+',', ',1)

    """
    
    
    """ check for / and initialize _EbfMap for directory walk"""
    keys=tagname.lower().strip().split()
    for key in keys:
        if key.endswith('/'):
            location=_EbfMap.get(filename, key,1)
            break
        
    numpy.set_printoptions(threshold='nan',precision=17)        
    if tableon == 1:                           
        data={}
        i=0
        skeys=[]
        for key in keys:
            datat=read(filename,key,0,0)
            if type(datat) == dict:
                for key1 in list(datat.keys()):
                    if datat[key1].ndim == 2:
                        for j in numpy.arange(0,datat[key1].shape[1]):
                            data[key1+"_"+str(i)]=datat[key1][:,j]
                            i=i+1                
                            skeys.append(key1+"_"+str(i))
                    elif datat[key1].ndim == 1:
                        data[key1]=datat[key1]
                        i=i+1                
                        skeys.append(key1)
                    elif datat[key1].ndim == 0:
                        data[key1]=numpy.array([datat[key1]])
                        i=i+1                
                        skeys.append(key1)
                    else:
                        raise RuntimeError('EBF Error: cannot print array with ndim >2')                    
            else:
                if datat.ndim == 0:
                    data[str.rpartition(key,'/')[2]]=numpy.array([datat])
                else:
                    data[str.rpartition(key,'/')[2]]=datat
                skeys.append(str.rpartition(key,'/')[2])
                i=i+1
        
        for key in list(data.keys()):
            if data[key].dtype.kind == 'V':
                for key2 in data[key].dtype.names:
                    data[key2]=data[key][key2]
                    skeys.append(key2)                                    
                    i=i+1
                del data[key]
                skeys.remove(key)
                i=i-1
        if len(data) != i:
            raise RuntimeError('Error in  ebf.cat(), duplicate input keys')

        width=0
        formatstring=[]
        formatstringh=[]
        for key in skeys:
            if (data[key].dtype.kind=='f') or (data[key].dtype.kind=='c'):
                formatstring.append("{0:>25.17}")
                formatstringh.append("{0:>25}")
            elif (data[key].dtype.kind=='S'):
                formatstring.append("{0:>"+str(data[key].dtype.itemsize)+"}")            
                formatstringh.append("{0:>"+str(data[key].dtype.itemsize)+"}")            
            else:
                formatstring.append("{0:>25}")
                formatstringh.append("{0:>25}")

        print(delimiter.join(formatstringh[j].format(key) for j,key in enumerate(skeys)))
        elements=min([data[key].size if data[key].ndim<2 else data[key].shape[0] for key in skeys])
        for i in numpy.arange(0,elements):
            print(delimiter.join(formatstring[j].format(data[key][i]) for j,key in enumerate(skeys)))

        
#        width=0
#        for key in data.keys():
#            if width < len("{0:<}".format(data[key][0])):
#                width=len("{0:<}".format(data[key][0]))                
#        formatstring="{0:>"+str(width+4)+"}"                
##        skeys=sorted(data.keys())        
#        elements=min([data[key].size if data[key].ndim<2 else data[key].shape[0] for key in skeys])
#        print delimiter.join(formatstring.format(key) for key in skeys)
#        for i in numpy.arange(0,elements):
#            print delimiter.join(formatstring.format(data[key][i]) for key in skeys)
        
#        elements=-1
#        for key in data.keys():
#            if data[key].ndim == 0:
#                elements1=data[key].size
#            else:
#                elements1=data[key].shape[0]
#            if elements == -1:
#                elements=elements1                
#            elif elements > elements1:
#                elements=elements1                
#
#        for i in numpy.arange(0,elements):
#            if data[skeys[0]].ndim == 0:
#                print delimiter.join(formatstring.format(data[key]) for key in skeys)
#            else:
#                print delimiter.join(formatstring.format(data[key][i]) for key in skeys)
        
                
    else:
        data={}
        i=0
        for key in keys:
            datat=read(filename,key,0,0)
            if type(datat) == dict:
                for key1 in list(datat.keys()):
                    width=8
                    while width < len(key1):
                        width=width*2
                    formatstring="{0:<"+str(width)+"}"
                    temp=formatstring.format(str.rpartition(key1,'/')[2])+'= '
                    if(datat[key1].dtype.type==numpy.string_):
                        datat1=datat[key1].tostring()
                    else:
                        datat1=numpy.array_str(numpy.squeeze(datat[key1]))                    
                    if len(temp)+len(datat1) > 64:
                        temp=temp+'\n'
                    print(temp+datat1)
            else:
                width=8
                while width < len(key):
                    width=width*2
                formatstring="{0:<"+str(width)+"}"
                temp=formatstring.format(str.rpartition(key,'/')[2])+'= '
                if(datat.dtype.type==numpy.string_):
                    datat=datat.tostring()
                    datalen=1
                else:
                    datat=numpy.squeeze(datat)
                    datalen=datat.size
                if len(temp)+datalen > 64:
                    if len(keys) > 1:
                        print(temp)
                    print("[")
                    for d in datat:
                        print(d)
                    print("]")
                                        
                else:
                    if len(keys) > 1:
                        print(temp,datat)
                    else:
                        print(datat)
    
    numpy.set_printoptions()    
    
    



def swapEndian(filename):
    """
    Swaps the endianess of the file. Little to Big or Big to Little

    Args:
         filename(str):

    Example:

    >>> ebf.swapEndian("check.ebf")

    """
    
    filename_out = filename.rpartition('.ebf')[0]+'_swap.ebf'

    fp1 = open(filename, 'rb')
    fp1.seek(0,2)
    filesize=fp1.tell()
    fp1.seek(0,0)
    header1=_EbfHeader()
    header2=_EbfHeader()    
    header1.read(fp1)
    fp1.seek(0,0)
    if header1.flagswap == 0:
        flagswap=1
        _EbfTable.init_swap(filename_out)
    else:
        flagswap=0
        _EbfTable.init(filename_out)
                
    fout = open(filename_out, 'rb+')
    fout.seek(0,2)
    keys=[]
    values=[]
    while fp1.tell() < filesize:
        header1.read(fp1)
        loc=fp1.tell()
        if header1.datatype == 8:
#            dth1 = numpy.dtype(__sdef2descr(header1.sdef.split()))
            dth1 = numpy.dtype(sdef2descr(header1.sdef)[0])
        else:
            dth1 = _TypeManager.itos_s(header1.datatype)
        if header1.datatype == 1:
            dth1 = 'S'+str(header1.dim[-1])
            
        dblock_size=header1.elements()*header1.datasize    
        
        if ((header1.name.startswith('/.ebf/')==False) and (header1.name.startswith('/.tr/')==False)):                  
            data = numpy.fromstring(fp1.read(dblock_size), dtype = dth1)
            if (header1.flagswap == 1) and (flagswap==0):
                data = data.byteswap(True)
            if (header1.flagswap == 0) and (flagswap==1):
                data = data.byteswap(True)
                
            if header1.datatype == 1:
                data = data.reshape(header1.getshape()[0:-1])
            else:
                data = data.reshape(header1.getshape())
            header2.create(header1.name,data,header1.dataunit,header1.sdef)

            if flagswap == 1:
                header2.flagswap = 1            
                
            keys.append(header2.name)
            values.append(fout.tell())
            header2.write(fout)            
            fout.write(data.tostring('C'))
                        
        fp1.seek(loc+header1.capacity(), 0)                
        
    filesize1=fp1.tell()
    fp1.close()
    fout.close()
    if filesize1 != filesize:
        raise RuntimeError('EBFCorrupt')
    if len(keys) > 0:
        _EbfTable.put(filename_out, keys, values)    
        
def copy(filename1,filename2,mode='a',tagnames='',outpath=None):
    """
    copy data items from one file to another

    Args:
         filename1(str):
         
         filename2(str):
         
         mode(str)     : 'w' or 'a'   
         
         tagnames(str) : if blank then copies all items or else one can \         
         supply space separated list of data items as a single string
         
         outpath(str): Path ending with '/' into which to copy items

    Example:

    >>> ebf.copy("check1.ebf",'check2.ebf','/x /y','w')
    >>> ebf.copy("check1.ebf",'check2.ebf','/x')
    >>> ebf.copy("check1.ebf",'check2.ebf')

    """
    if tagnames == '':
        keys=_EbfTable.getKeyVals(filename1)[0]
    else:
        keys=tagnames.split()
        for key in keys:
            if containsKey(filename1,key) == 0:
                raise RuntimeError('EBF Error: in copy(), key not present in input file')
            
    keyst=keys
    keys=[]
    for key in keyst:
        if (key.startswith('/.ebf/') == False)and(key.startswith('/.tr/') == False):
            keys.append(key)
            
            
    if os.path.isfile(filename2) == False:
        mode='w'
                
    if mode == 'a':
        for key in keys:
            if containsKey(filename2,key) == 1:
                raise RuntimeError('EBF Error: in copy(), key already present in output file')
    elif mode != 'w':
        raise RuntimeError('EBF Error: in copy(), mode must be w or a')
    
    if mode == 'w':
        initialize(filename2)
        mode='a'        

    if outpath is not None:
        if outpath[-1]!='/':    
            raise RuntimeError('EBF Error: in copy(), outpath must end in /')

        
    for key in keys:
        data=read(filename1,key)
        if outpath is None:    
            write(filename2,key,data,mode)
        else:
            write(filename2,outpath,data,mode)
            

    


def _checkSpeed():
    print('Test speed read and write-->')
    
        
    start=time.time()    
    nsize=1000
    data1 = numpy.zeros(2, dtype = 'int32')+1
    print('\n item read/write speed: in Kilo operations per second KOPS:')
    print('size of each item',data1.size)
    print('Number of items',nsize)
    
    write('check.ebf', '/x0', data1, 'w')    
    for i in numpy.arange(1, nsize):
        write('check.ebf', '/x'+str(i), data1, 'a')
    print('\t Writing speed=', nsize*1e-3/(time.time()-start), ' Kops (' , (time.time()-start),' s)')
#    info('test12.ebf')

    start=time.time()    
#    tot=1
    for i in numpy.arange(0, nsize):
        y = read('check.ebf', '/x'+str(i))
#        print y[0]
#        tot = y[0]+1
    print('\t Reading speed=', nsize*1e-3/(time.time()-start), ' Kops (',(time.time()-start),' s)')

    
    print('\n get key list:')
    start = time.time()    
    keys=_EbfTable.getKeyVals('check.ebf')[0]
    print('\t Reading speed Key vals HT=', nsize*1e-3/(time.time()-start), ' Kops (' ,(time.time()-start),'s) ',' keys=',len(keys))
    
    start = time.time()    
    keys=_EbfTable.getKeyValsIT('check.ebf')[0]
    print('\t Reading speed Key vals IT =', nsize*1e-3/(time.time()-start), ' Kops (' ,(time.time()-start),' s) ',' keys=',len(keys))



    print('\n data read/write speed: in MB/s:')
    nsize = [10, 100, 1000, 10000, 100000, 1000000, 10000000, 100000000]
    nsize = [10000000]
    for j in nsize:
        print(j, 10)
        data1 = numpy.zeros(j, dtype = 'int32')
        start = time.time()    
        for i in numpy.arange(1, 10):
            write('check.ebf', '/x1', data1, 'w')
        print('\t Writing speed=', data1.size*data1.itemsize*1e-6*10/(time.time()-start), ' MB/s (', (time.time()-start),' s)')


        start = time.time()    
        for i in numpy.arange(1, 10):
            y = read('check.ebf', '/x1')
        print('\t Reading speed=', data1.size*data1.itemsize*1e-6*10/(time.time()-start), ' MB/s (', (time.time()-start),' s)')


def diff(filename1,filename2):
    """
    Perform diff operation on two files. Ignores data items starting with "/." which are 
    for internal use. If file contents are same it does not print anything.

    Args:
         filename1(str):
         
         filename2(str):

    Example:

    >>> ebf.diff("check1.ebf","check2.ebf")

    """
    keys1=_EbfTable.getKeyVals(filename1)[0]
    keys2=_EbfTable.getKeyVals(filename2)[0]
    
    temp=keys1
    keys1=[]
    for key in temp:
        if (key.startswith('/.ebf/') == False)and(key.startswith('/.tr/') == False):
            keys1.append(key)
    temp=keys2
    keys2=[]
    for key in temp:
        if (key.startswith('/.ebf/') == False)and(key.startswith('/.tr/') == False):
            keys2.append(key)
            
    if len(keys1) != len(keys2):
        print('files differ: unequal number of data itmes, ', len(keys1),' and ',len(keys2))
            
    count_differ=0
    count_match=0
    for key in keys1:
        if keys2.count(key) == 1:
            data1=read(filename1,key).tostring() 
            data2=read(filename2,key).tostring()
            if data1 != data2:
                print('data item->',key,' differs')
                count_differ=count_differ+1
            else:
                count_match=count_match+1
        else:
            print('data item->',key,' not present in second file')
    
    if count_match != len(keys1):
        print(len(keys1)-count_match,' data items differ out of',len(keys1), 'items in first file')
        
    
    


def _usage():
#        print 'To run test suite'
#        print '     ebftkpy -test'
#        print 'To get summary of file'
#        print '     ebftkpy -list filename'
#        print '     ebftkpy filename'
#        print 'To print a data item'
#        print '     ebftkpy -cat filename tagname'            
#        print 'To get statistics of data in file'
#        print '     ebftkpy -stat filename'
#        print 'To swap endianess of file'
#        print '     ebftkpy -swap filename'
#        print 'To check speed of input output'
#        print '     ebftkpy -speed filename'        
        print("NAME:")
        print('\t >>EBF<<  (Efficient and Easy to use Binary File Format)')
        print("\t ebftkpy 0.0.1 - a toolkit for  EBF  files")
        print("\t Copyright (c) 2012 Sanjib Sharma ")
        print("USAGE:")
        print("\t ebftkpy\t -list filename")
        print("\t ebftkpy\t  filename  (same as -list)")
        print("\t ebftkpy\t -cat filename \"TagName1 TagName2 ..\"")
        print("\t ebftkpy\t -csv filename \"TagName1 TagName2 ..\"")
        print("\t ebftkpy\t -ssv filename \"TagName1 TagName2 ..\"")
        print("\t ebftkpy\t -stat filename \"TagName1 TagName2 ..\"")
        print("\t ebftkpy\t -swap filename")
        print("\t ebftkpy\t -copy src_file dest_file")
        print("\t ebftkpy\t -copy src_file dest_file TagName")
        print("\t ebftkpy\t -diff  filename1 filename2")
        print("\t ebftkpy\t -rename  filename1 tagname_old tagname_new")
        print("\t ebftkpy\t -remove  filename1 tagname")
        print("\t ebftkpy\t -htab filename")
        print("DESCRIPTION:")
        print("\t -list    ","view headers/TagNames of data in file ")    
        print("\t -cat     ","print data in ascii format")
        print("\t          ","e.g., for \"TagName1\" a record of rank 2 with")
        print("\t          ","dimensions N and 3 will print a Nx3 table,")
        print("\t          ","for \"TagName2\" a record of rank 1 with dimension N")
        print("\t          ","will print a column of size N")
        print("\t          ","multiple tags can be specified as space separated ")
        print("\t          ","strings as \"TagName1 TagName2\"  ")
        print("\t          ","but the condition is that the number of elements in")
        print("\t          ","each record should be same. This will print a Nx4 table")
        print("\t -csv     ","print data in csv tabular format, syntax same as cat")
        print("\t -ssv     ","print data in csv tabular format, but delimitier as space")
        print("\t -stat    ","print min max mean stddev of specified data tags")
        print("\t -swap    ","swap the endianness of a file, output file has") 
        print("\t          ","suffix _swap.ebf")
        print("\t -copy    ","copy contents of one file to another or only a tag")
        print("\t -diff    ","difference of two data items in two ebf files")
        print("\t -rename  ","rename a data item")
        print("\t -remove  ","remove a data item. It is renamed with prefix /.tr/ ")
        print("\t          ","which can be restored using rename if needed")
        print("\t -htab    ","get information about internal hashtable")
        print("CONTACT:")
        print("http://ebfformat.sourceforge.net")
    


    
import unittest
class _ebf_test(unittest.TestCase):
    
    def setUp(self):
        self.seq = list(range(10))
        
#    def test_expand(self):
#        x1=numpy.zeros(10)        
#        x2=numpy.ones(10)        
        
        
    def test_ebfht(self):        
        print('Testing ebftable get, put and remove-->')        
        _EbfTable.init('check.txt')
        nsize=100
        
        start = time.time()    
        for i in numpy.arange(0, nsize):
            _EbfTable.put('check.txt','/x'+str(i),i*10)
        print('Writing  ', nsize*1e-3/(time.time()-start), ' Kops') 
                
        start = time.time()    
        x=numpy.zeros(nsize)
        for i in numpy.arange(0, nsize):
            x[i]=_EbfTable.get('check.txt','/x'+str(i))
        print('Reading  ', nsize*1e-3/(time.time()-start), ' Kops') 
        
        status1=1
        for i in numpy.arange(0, nsize):
            if x[i] != i*10:
                status1=0
                print(i,x[i])
                
        self.assertEqual(status1, 1)
        
        
        start = time.time()    
        x=numpy.zeros(nsize)
        for i in numpy.arange(0, nsize/2):
            _EbfTable.remove('check.txt','/x'+str(i))
        print('Removing ', nsize*1e-3/(time.time()-start), ' Kops') 
        
        for i in numpy.arange(0, nsize):
            x[i]=_EbfTable.get('check.txt','/x'+str(i))
            
        status2=1
        for i in numpy.arange(0, nsize/2):
            if x[i] != -1:
                status2=0
        for i in numpy.arange(nsize/2, nsize):
            if x[i] != i*10:
                status2=0
        
        
        self.assertEqual(status2, 1)
        
#        start = time.time()    
#        keys=_EbfTable.getKeyValsIT('check.txt')[0]
#        print 'Reading Key vals raw', 1e3*1e-3/(time.time()-start), ' Kops' ,(time.time()-start),len(keys)
    
#        start = time.time()    
#        fileht=_EbfTable('check.txt','rb')
#        keys=_EbfTable.getKeyVals('check.txt')[0]
#        fileht.close()
#        print 'Reading Key vals HT', 1e3*1e-3/(time.time()-start), ' Kops' ,(time.time()-start),len(keys)
        
    
    
        
        
    def test_header256(self):
        print("Testing header256-->")
        ebfdir='data/'
        data=read(ebfdir+'header256.ebf','/')
        x=data["xvar"]
        y=data["yvar"]
        z=numpy.arange(0,10)
        self.assertEqual(x.size,z.size)
        self.assertEqual(y.size,z.size)
        self.assertEqual(numpy.sum(x==z),x.size)
        self.assertEqual(numpy.sum((y-10)==z),x.size)
        
        
        
    def testtable(self):
        """ Check rename """
        print("Testing ebftable-->")
        x = numpy.arange(-65636, -65636-128,-1, dtype = "int64")
        write("check_table.ebf", "/x1", x, "w")      
        write("check_table.ebf", "/x2", x, "a")      
        write("check_table.ebf", "/x3", x, "a")      
        rename("check_table.ebf","/x1","/x5")
        rename("check_table.ebf","/x2",'')
        y1=read('check_table.ebf','/x5')
        y2=read('check_table.ebf','/.tr/x2.0')
        self.assertEqual(numpy.sum(x ==y1),x.size)
        self.assertEqual(numpy.sum(x ==y2),x.size)
        
        swapEndian('check_table.ebf')
        write("check_table_swap.ebf", "/x6", x, "a")      
        write("check_table_swap.ebf", "/x7", x, "a")      
        y1=read('check_table_swap.ebf','/x5')
        y2=read('check_table_swap.ebf','/x3')
        y3=read('check_table_swap.ebf','/x6')
        y4=read('check_table_swap.ebf','/x7')
#        info('check_table_swap.ebf')
        self.assertEqual(numpy.sum(x ==y1),x.size)
        self.assertEqual(numpy.sum(x ==y2),x.size)
        self.assertEqual(numpy.sum(x ==y3),x.size)
        self.assertEqual(numpy.sum(x ==y4),x.size)
        
        
    def teststring(self):
        """ Check string read write """
        print("Testing string read/write-->")
        x = "ebcdefgh"
        write("check.ebf", "/mystr", numpy.array(x), "w")      
        y = read("check.ebf", "/mystr").tostring()
        self.assertEqual(x, y)
        x = numpy.array(['aa','ba','ca','da'])
        write("check.ebf", "/mystr", x, "w")      
        y = read("check.ebf", "/mystr")
        self.assertEqual(numpy.all(x==y),True)
        
        x = numpy.array(['a','b','c','d'])
        write("check.ebf", "/mystr", x, "w")      
        y = read("check.ebf", "/mystr")
        self.assertEqual(numpy.all(x==y),True)
        
        
    def testdataunit(self):
        """ Check data units read write"""
        print("Testing data unit-->")
        write('check.ebf', '/data', numpy.zeros(1, dtype = "int32"), "w", dataunit = "100 m/s")
        self.assertEqual("100 m/s", unit('check.ebf', '/data'))
        
    def testexceptions(self):
        """ Check overwrite protection, write in between and then check """
        print("Testing exceptions-->")
        x = numpy.zeros(10, dtype = 'int32')
        write('check.ebf', '/x', x, "w")
        write('check1.ebf', '/x', x, "w")
        self.assertRaises(RuntimeError, write, 'check.ebf', '/x', x, "a")
        self.assertRaises(RuntimeError, read, 'check1.ebf', '/x1')
        self.assertRaises(IOError, read, 'check123.ebf', '/x1')
        
    def testchecksum(self):
        print("Testing checksum-->")
        nsize=10
        x1=numpy.zeros(nsize,dtype='float32')
        x2=numpy.zeros(nsize,dtype='float64')
        write("check.ebf","/x1",x1,"w")
        write("check.ebf","/x2",x2,"a")
        write("check.ebf","/single/x1",x1[0:1],"a")
        write("check.ebf","/single/x2",x2[0:1],"a")
        checksum=read("check.ebf",'/.ebf/info')
#        info("check.ebf")
        print(checksum);
        print("(EBF, 0)    hash=",_EbfTable.ebfckhash("(EBF, 0) ",0))
        print("(EBF, 1000) hash=",_EbfTable.ebfckhash("(EBF, 1000)",1000))

    def testmultiple(self):
        """ Check overwrite protection, write in between and then check """
        print("Testing mutiple read write-->")
        x = numpy.zeros(10, dtype = 'int32')
        write('check1.ebf', '/x1', x, "w")
        write('check1.ebf', '/x2', x, "a")
        write('check2.ebf', '/x1', x, "w")
        y = read('check2.ebf', '/x1')
        self.assertEqual(y.size, x.size)
        write('check2.ebf', '/x2', x, "a")
        y = read('check2.ebf', '/x1')
        self.assertEqual(y.size, x.size)
        write('check1.ebf', '/y1', x, "w")
        y = read('check2.ebf', '/x1')
        self.assertEqual(y.size, x.size)
        write('check1.ebf', '/y2', x, "a")        
        y = read('check2.ebf', '/x1')
        self.assertEqual(y.size, x.size)
        write('check2.ebf', '/x1', x, "w")
        y = read('check2.ebf', '/x1')
        self.assertEqual(y.size, x.size)
        write('check2.ebf', '/x2', x, "a")
        y1 = read('check2.ebf', '/x1')
        x1 = read('check1.ebf', '/y1')
        self.assertEqual(numpy.sum(x == x1), x.size)
        self.assertEqual(numpy.sum(x == y1), x.size)
        self.assertRaises(RuntimeError, read, 'check1.ebf', '/x1')
        
    def write_master(self):
        # make sure the shuffled sequence does not lose any elements
        print("write master test file-->")
        data = {}
        keys = ["x1", "x2", "x3", "x4", "x5", "x6", "x9", "x10", "x11", "x12", "x13"]
        x=numpy.arange(0,128,dtype='int8')
        data["x1"]  =  numpy.array(x)
#        data["x1"]  =  numpy.fromstring(x.tostring(),dtype='S1')
        data["x9"] = numpy.arange(-128, 0, dtype = 'int8')
        data["x6"] = numpy.arange(-256, -256-128,-1, dtype = 'int16')
        data["x2"] = numpy.arange(-65636, -65636-128,-1, dtype = "int32")
        data["x3"] = numpy.arange(-4294967296, -4294967296-128,-1, dtype = "int64")
        data["x4"] = numpy.array(numpy.linspace(1.23e20, 128.23e20, 128), dtype = "float32")
        data["x5"] = numpy.array(numpy.linspace(1.23456789e200, 128.23456789e200, 128), dtype = "float64")
        data["x10"] = numpy.arange(128, 128+128, dtype = 'uint8')
        data["x11"] = numpy.arange(256, 256 + 128, dtype = 'uint16')
        data["x12"] = numpy.arange(65636, 65636 + 128, dtype = 'uint32')
        data["x13"] = numpy.arange(4294967296, 4294967296 + 128, dtype = 'uint64')
#        ebfdir='/home/sharma/sw/share/ebf/'
        ebfdir='data/'
        
        write(ebfdir+'master_test1.ebf', '/', data, 'w')    
        for key in keys:
            if key != 'x1':
                newsize = data[key].size/4
                data[key] = data[key].reshape(4, newsize)
                
        
#        write('/home/sharma/ebf_demo/test1.ebf', '/', data, 'w',"100 km/s")    
#        write('/home/sharma/ebf_demo/test1.ebf', '/2d/', data, 'a',"kpc")    


        dt1 = []        
        dt1.append(("x2", data["x2"].dtype, data["x2"].shape))
        dt1.append(("x3", data["x3"].dtype, data["x3"].shape))         
        
        data2 = numpy.zeros(2, dtype = dt1)
        data2[0]["x2"][:][:] = data["x2"][:][:]
        data2[0]["x3"][:][:] = data["x3"][:][:]
        data2[1]["x2"][:][:] = data["x2"][:][:]
        data2[1]["x3"][:][:] = data["x3"][:][:]
        
        dt2 = []        
        dt2.append(("x2", data["x2"].dtype, data["x2"].shape))
        dt2.append(("x3", data["x3"].dtype, data["x3"].shape))         
        dt2.append(('point1', dt1, (2, )))
        
        data3 = numpy.zeros(1, dtype = dt2)[0]
        data3['x2'][:][:] = data["x2"][:][:]
        data3['x3'][:][:] = data["x3"][:][:]
        data3['point1'][0] = data2[0].copy()
        data3['point1'][1] = data2[0].copy()            
        
        dt2 = []        
        dt2.append(("x2", data["x2"].dtype, data["x2"].shape))
        dt2.append(("x3", data["x3"].dtype, data["x3"].shape))         
        dt2.append(('point1', dt1, (1, )))
        
        dt3 = []        
        dt3.append(("x2", data["x2"].dtype, data["x2"].shape))
        dt3.append(("x3", data["x3"].dtype, data["x3"].shape))         
        dt3.append(('point2', dt2, (1, )))

        
        data4 = numpy.zeros(1, dtype = dt2)[0]
        data4['x2'][:][:] = data["x2"][:][:]
        data4['x3'][:][:] = data["x3"][:][:]
        data4['point1'][0] = data2[0].copy()
        
        
        data5 = numpy.zeros(1, dtype = dt3)[0]
        data5['x2'][:][:] = data["x2"][:][:]
        data5['x3'][:][:] = data["x3"][:][:]
        data5['point2'][0] = data4.copy()
            
        write(ebfdir+'master_test1.ebf', '/dir1/data_struct', data2, 'a')    
        write(ebfdir+'master_test1.ebf', '/dir1/data_struct_rec', data3, 'a')    
        write(ebfdir+'master_test1.ebf', '/dir1/data_struct_rec2', data5, 'a')    
        swapEndian(ebfdir+'master_test1.ebf')
        self.assertEqual(1, 1)
        
        
            
            

    def test_read_masterfile(self):
        # make sure the shuffled sequence does not lose any elements
        print("Testing read masterfile-->")
        data = {}
        keys = ["x1", "x2", "x3", "x4", "x5", "x6", "x9", "x10", "x11", "x12", "x13"]
        x=numpy.arange(0,128,dtype='int8')
        data["x1"]  =  numpy.array(x)
#        data["x1"]  =  numpy.fromstring(x.tostring(),dtype='S1')
        data["x9"] = numpy.arange(-128, 0, dtype = 'int8')
        data["x6"] = numpy.arange(-256, -256-128,-1, dtype = 'int16')
        data["x2"] = numpy.arange(-65636, -65636-128,-1, dtype = "int32")
        data["x3"] = numpy.arange(-4294967296, -4294967296-128,-1, dtype = "int64")
        data["x4"] = numpy.array(numpy.linspace(1.23e20, 128.23e20, 128), dtype = "float32")
        data["x5"] = numpy.array(numpy.linspace(1.23456789e200, 128.23456789e200, 128), dtype = "float64")
        data["x10"] = numpy.arange(128, 128+128, dtype = 'uint8')
        data["x11"] = numpy.arange(256, 256 + 128, dtype = 'uint16')
        data["x12"] = numpy.arange(65636, 65636 + 128, dtype = 'uint32')
        data["x13"] = numpy.arange(4294967296, 4294967296 + 128, dtype = 'uint64')
        
#        ebfdir='/home/sharma/sw/share/ebf/'
        ebfdir='data/'
        filename1=ebfdir+'master_test1.ebf'
        filename2=ebfdir+'master_test1_swap.ebf'
        
        data1={}
        for key in keys:
            data1[key] = data[key].copy()

        for key in keys:
            if key != 'x1':
                newsize = data[key].size//4
                data[key] = data[key].reshape(4, newsize)        


        dt1 = []        
        dt1.append(("x2", data["x2"].dtype, data["x2"].shape))
        dt1.append(("x3", data["x3"].dtype, data["x3"].shape))         
        
        data2 = numpy.zeros(2, dtype = dt1)
        data2[0]["x2"][:][:] = data["x2"][:][:]
        data2[0]["x3"][:][:] = data["x3"][:][:]
        data2[1]["x2"][:][:] = data["x2"][:][:]
        data2[1]["x3"][:][:] = data["x3"][:][:]
        


        
                     
        datar1 = read(filename1, "/")
        for key in list(data1.keys()):
            self.assertEqual(numpy.sum(datar1[key] == data1[key]), data1[key].size)

        
        datar2 = read(filename1, "/dir1/data_struct_rec")[0]
        datar3 = read(filename1, "/dir1/data_struct_rec2")[0]        
        datar4 = read(filename1, "/dir1/data_struct")[1]
        for key in ['x2','x3']:
            self.assertEqual(numpy.sum(datar4[key] == data[key]), data[key].size)
            self.assertEqual(numpy.sum(datar2[key] == data[key]), data[key].size)
            self.assertEqual(numpy.sum(datar2['point1'][1][key] == data[key]), data[key].size)
            self.assertEqual(datar2[key].shape, data[key].shape)
            self.assertEqual(datar2['point1'][1][key].shape, data[key].shape)
            self.assertEqual(datar3['point2'][0][key].shape, data[key].shape)
            self.assertEqual(datar4[key].shape, data[key].shape)
            self.assertEqual(numpy.sum(datar3['point2']['point1'][key] == data[key]), data[key].size)

            
        datar1 = read(filename2, "/")
        for key in list(data1.keys()):
            self.assertEqual(numpy.sum(datar1[key] == data1[key]), data1[key].size)
        datar2 = read(filename2, "/dir1/data_struct_rec")[0]
        datar3 = read(filename2, "/dir1/data_struct_rec2")[0]        
        datar4 = read(filename1, "/dir1/data_struct")[1]
        for key in ['x2','x3']:
            self.assertEqual(numpy.sum(datar4[key] == data[key]), data[key].size)
            self.assertEqual(numpy.sum(datar2[key] == data[key]), data[key].size)
            self.assertEqual(numpy.sum(datar2['point1'][1][key] == data[key]), data[key].size)
            self.assertEqual(datar2[key].shape, data[key].shape)
            self.assertEqual(datar2['point1'][1][key].shape, data[key].shape)
            self.assertEqual(datar3['point2'][0][key].shape, data[key].shape)
            self.assertEqual(datar4[key].shape, data[key].shape)
            self.assertEqual(numpy.sum(datar3['point2']['point1'][key] == data[key]), data[key].size)
        
        

    def testcorrectness1(self):
        # make sure the shuffled sequence does not lose any elements
        print("Testing correctness ver-1 of read write-->")
        data = {}
        keys = ["x1", "x2", "x3", "x4", "x5", "x6", "x9", "x10", "x11", "x12", "x13"]
        data["x1"]  =  numpy.array(['EBF','EBFGH'])
        data["x2"] = numpy.arange(-65636, -65636 - 128, -1, dtype = "int32")
        data["x3"] = numpy.arange(-4294967296, -4294967296 - 128, -1, dtype = "int64")
        data["x4"] = numpy.array(numpy.linspace(3.1459 * 1e-30, (3.14159 + 127) * 1e-30, 128), dtype = "float32")
        data["x5"] = numpy.array(numpy.linspace(3.1459 * 1e-300, (3.14159 + 127) * 1e-300, 128), dtype = "float64")
        data["x6"] = numpy.arange(-256, -256 - 128, -1, dtype = 'int16')
        data["x9"] = numpy.arange(-126, 2, dtype = 'int8')
        data["x10"] = numpy.arange(0, 128, dtype = 'uint8')
        data["x11"] = numpy.arange(256, 256 + 128, dtype = 'uint16')
        data["x12"] = numpy.arange(65636, 65636 + 128, dtype = 'uint32')
        data["x13"] = numpy.arange(4294967296, 4294967296 + 128, dtype = 'uint64')
                
        for key in keys:
            if key != "x1":
                newsize = data[key].size//4
                data[key] = data[key].reshape(4, newsize)
                        
                
        dt = []        
        for key in keys:
#            if key !=  "x1":
            dt.append((key, data[key].dtype, data[key].shape))
         
        data2 = numpy.zeros(1, dtype = dt)

        for key in keys:
#            if key != "x1":
            data2[key][:][:] = data[key][:][:]
        
        dth = list(dt)
        dth.append(('point1', dt, (2, )))
        data22 = numpy.zeros(1, dtype = dth)[0]
        data22['point1'][0] = data2[0].copy()
        data22['point1'][1] = data2[0].copy()
        for key in keys:
            data22[key][:][:] = data[key][:][:]

        dth = list(dt)
        dth.append(('point1', dt, (1, )))
        data23 = numpy.zeros(1, dtype = dth)[0]
        data23['point1'][0] = data2[0].copy()
        for key in keys:
            data23[key][:][:] = data[key][:][:]
            
            
        
        fout = open("check.ebf", "wb")
        fout.close()
        write("check.ebf", "/", data, "w")
        write("check.ebf", "/struct1/data2", data2, "a")
        write("check.ebf", "/dir3/", data2, "a")
        data2 = data2[0]
        write("check.ebf", "/dir4/", data2, "a")
        write("check.ebf", "/struct2/data22", data22, "a")
        write("check.ebf", "/dir5/", data23, "a")
        
        for key in keys:
            write("check.ebf", "/dir1/"+key, data[key], "a")
            write("check.ebf", "/dir2/"+key, data2[key], "a")
#        info('check.ebf')
             
        data1 = read("check.ebf", "/")
        data3 = read("check.ebf", "/struct1/data2")[0]
        data33 = read("check.ebf", "/struct2/data22")
        self.assertEqual(len(data), len(data1))        
        for key in keys:
            x1 = read("check.ebf", "/"+key)            
            x2 = read("check.ebf", "/dir1/"+key)            
            x3 = read("check.ebf", "/dir2/"+key)            
            x4 = read("check.ebf", "/dir3/"+key)            
            x5 = read("check.ebf", "/dir4/"+key)            
            x5 = read("check.ebf", "/dir5/"+key)
#            print 'here',key,x1.dtype,data[key].dtype            
            self.assertEqual(numpy.sum(x1 == data[key]), data[key].size)
            self.assertEqual(numpy.sum(x2 == data[key]), data[key].size)
            self.assertEqual(numpy.sum(x3 == data[key]), data[key].size)
            self.assertEqual(numpy.sum(x4 == data[key]), data[key].size)
            self.assertEqual(numpy.sum(x5 == data[key]), data[key].size)
            self.assertEqual(x1.shape, data[key].shape)
            self.assertEqual(x2.shape, data[key].shape)
            self.assertEqual(x3.shape, data[key].shape)
            self.assertEqual(x4.shape, data[key].shape)
            self.assertEqual(numpy.sum(data1[key] == data[key]), data[key].size)
            self.assertEqual(numpy.sum(data3[key] == data[key]), data[key].size)
            self.assertEqual(numpy.sum(data33['point1'][1][key] == data[key]), data[key].size)
            self.assertEqual(data1[key].shape, data[key].shape)
            self.assertEqual(data3[key].shape, data[key].shape)
            self.assertEqual(data33['point1'][1][key].shape, data[key].shape)
        swapEndian("check.ebf")    
#        info('check.ebf')
        
        data1 = read("check_swap.ebf", "/")
        data3 = read("check_swap.ebf", "/struct1/data2")[0]
        data33 = read("check_swap.ebf", "/struct2/data22")
        self.assertEqual(len(data), len(data1))        
#        info('check.ebf')
#        info('check_swap.ebf')
        for key in keys:
#            print key
            x1 = read("check_swap.ebf", "/"+key)            
            x2 = read("check_swap.ebf", "/dir1/"+key)            
            x3 = read("check_swap.ebf", "/dir2/"+key)            
            x4 = read("check_swap.ebf", "/dir3/"+key)            
            x5 = read("check_swap.ebf", "/dir4/"+key)            
            x5 = read("check_swap.ebf", "/dir5/"+key)            
            self.assertEqual(numpy.sum(x1 == data[key]), data[key].size)
            self.assertEqual(numpy.sum(x2 == data[key]), data[key].size)
            self.assertEqual(numpy.sum(x3 == data[key]), data[key].size)
            self.assertEqual(numpy.sum(x4 == data[key]), data[key].size)
            self.assertEqual(numpy.sum(x5 == data[key]), data[key].size)
            self.assertEqual(x1.shape, data[key].shape)
            self.assertEqual(x2.shape, data[key].shape)
            self.assertEqual(x3.shape, data[key].shape)
            self.assertEqual(x4.shape, data[key].shape)
            self.assertEqual(numpy.sum(data1[key] == data[key]), data[key].size)
            self.assertEqual(numpy.sum(data3[key] == data[key]), data[key].size)
            self.assertEqual(numpy.sum(data33['point1'][1][key] == data[key]), data[key].size)
            self.assertEqual(data1[key].shape, data[key].shape)
            self.assertEqual(data3[key].shape, data[key].shape)
            self.assertEqual(data33['point1'][1][key].shape, data[key].shape)


    def test_correctness2(self):
        # make sure the shuffled sequence does not lose any elements
        print("Testing correctness ver-2 of read write-->")
        data = {}
        keys = ["x1", "x2", "x3", "x4", "x5", "x6", "x9", "x10", "x11", "x12", "x13"]
        x=numpy.arange(0,128,dtype='int8')
        data["x1"]  =  numpy.array(x)
#        data["x1"]  =  numpy.fromstring(x.tostring(),dtype='S1')
        data["x9"] = numpy.arange(-128, 0, dtype = 'int8')
        data["x6"] = numpy.arange(-256, -256-128,-1, dtype = 'int16')
        data["x2"] = numpy.arange(-65636, -65636-128,-1, dtype = "int32")
        data["x3"] = numpy.arange(-4294967296, -4294967296-128,-1, dtype = "int64")
        data["x4"] = numpy.array(numpy.linspace(1.23e20, 128.23e20, 128), dtype = "float32")
        data["x5"] = numpy.array(numpy.linspace(1.23456789e200, 128.23456789e200, 128), dtype = "float64")
        data["x10"] = numpy.arange(128, 128+128, dtype = 'uint8')
        data["x11"] = numpy.arange(256, 256 + 128, dtype = 'uint16')
        data["x12"] = numpy.arange(65636, 65636 + 128, dtype = 'uint32')
        data["x13"] = numpy.arange(4294967296, 4294967296 + 128, dtype = 'uint64')
        for key in keys:
            if key != 'x1':
                newsize = data[key].size//4
                data[key] = data[key].reshape(4, newsize)
                
        
#        write('/home/sharma/ebf_demo/test1.ebf', '/', data, 'w',"100 km/s")    
#        write('/home/sharma/ebf_demo/test1.ebf', '/2d/', data, 'a',"kpc")    


        dt1 = []        
        dt1.append(("x2", data["x2"].dtype, data["x2"].shape))
        dt1.append(("x3", data["x3"].dtype, data["x3"].shape))         
        
        data2 = numpy.zeros(2, dtype = dt1)
        data2[0]["x2"][:][:] = data["x2"][:][:]
        data2[0]["x3"][:][:] = data["x3"][:][:]
        data2[1]["x2"][:][:] = data["x2"][:][:]
        data2[1]["x3"][:][:] = data["x3"][:][:]
        
        dt2 = []        
        dt2.append(("x2", data["x2"].dtype, data["x2"].shape))
        dt2.append(("x3", data["x3"].dtype, data["x3"].shape))         
        dt2.append(('point1', dt1, (2, )))
        data3 = numpy.zeros(1, dtype = dt2)[0]
        data3['x2'][:][:] = data["x2"][:][:]
        data3['x3'][:][:] = data["x3"][:][:]
        data3['point1'][0] = data2[0].copy()
        data3['point1'][1] = data2[0].copy()            
        
        dt2 = []        
        dt2.append(("x2", data["x2"].dtype, data["x2"].shape))
        dt2.append(("x3", data["x3"].dtype, data["x3"].shape))         
        dt2.append(('point1', dt1, (1, )))
        
        dt3 = []        
        dt3.append(("x2", data["x2"].dtype, data["x2"].shape))
        dt3.append(("x3", data["x3"].dtype, data["x3"].shape))         
        dt3.append(('point2', dt2, (1, )))

        
        data4 = numpy.zeros(1, dtype = dt2)[0]
        data4['x2'][:][:] = data["x2"][:][:]
        data4['x3'][:][:] = data["x3"][:][:]
        data4['point1'][0] = data2[0].copy()
        
        
        data5 = numpy.zeros(1, dtype = dt3)[0]
        data5['x2'][:][:] = data["x2"][:][:]
        data5['x3'][:][:] = data["x3"][:][:]
        data5['point2'][0] = data4[0].copy()


        
        
        
        fout = open("check.ebf", "wb")
        fout.close()
        write("check.ebf", "/", data, "w")
        write("check.ebf", "/dir1/struct", data2, "a")
        write("check.ebf", "/dir1/struct_single", data2[0], "a")
        write("check.ebf", "/dir1/struct_rec", data3, "a")
        write("check.ebf", "/dir1/struct_rec2", data5, "a")
        write("check.ebf", "/struct_split/", data2[0], "a")
        write("check.ebf", "/dir1/struct_rec_split/", data4, "a")
        
#        info('check.ebf')
             
        datar1 = read("check.ebf", "/")
        datar2 = read("check.ebf", "/struct_split/")
        temp = read("check.ebf", "/dir1/struct")
        self.assertEqual(temp.size, data2.size)
        datar3=temp[0]
        datar4=temp[1]
        datar5 = read("check.ebf", "/dir1/struct_rec_split/",recon=1)
        datar6=datar5['point1']
        datar7 = read("check.ebf", "/dir1/struct_single")
        self.assertEqual(datar7.size, 1)
        datar8 = read("check.ebf", "/dir1/struct_rec")
        self.assertEqual(datar8.size, 1)
        datar9 = datar8['point1'][1]
#        datar5=datar5[0]['point1'][1]
        datar10 = read("check.ebf", "/dir1/struct_rec2")
        datar11=datar10['point2'][0]
        
        for key in data2[0].dtype.names:
            x1 = read("check.ebf", "/struct_split/"+key)            
            self.assertEqual(numpy.sum(x1 == data[key]), data[key].size)
            self.assertEqual(numpy.sum(datar2[key] == data[key]), data[key].size)
            self.assertEqual(numpy.sum(datar3[key] == data[key]), data[key].size)
            self.assertEqual(numpy.sum(datar4[key] == data[key]), data[key].size)
            self.assertEqual(numpy.sum(datar5[key] == data[key]), data[key].size)
            self.assertEqual(numpy.sum(datar6[key] == data[key]), data[key].size)
            self.assertEqual(numpy.sum(datar7[key] == data[key]), data[key].size)
            self.assertEqual(numpy.sum(datar8[key] == data[key]), data[key].size)
            self.assertEqual(numpy.sum(datar9[key] == data[key]), data[key].size)
            self.assertEqual(x1.shape, data[key].shape)
            self.assertEqual(datar2[key].shape, data[key].shape)
            self.assertEqual(datar3[key].shape, data[key].shape)
            self.assertEqual(datar4[key].shape, data[key].shape)
            self.assertEqual(datar5[key].shape, data[key].shape)
            self.assertEqual(datar6[key].shape, data[key].shape)
            self.assertEqual(datar7[key].shape, data[key].shape)
            self.assertEqual(datar8[key].shape, data[key].shape)
            self.assertEqual(datar9[key].shape, data[key].shape)
            self.assertEqual(datar11[key].shape, data[key].shape)

        for key in list(data.keys()):
            x1 = read("check.ebf", "/"+key)
#            print key, data[key].shape, x1.shape,x1.dtype            
            self.assertEqual(numpy.sum(x1 == data[key]), data[key].size)
            self.assertEqual(numpy.sum(datar1[key] == data[key]), data[key].size)
            self.assertEqual(datar1[key].shape, data[key].shape)
            self.assertEqual(x1.shape, data[key].shape)


    def test_correctness3(self):
        # make sure the shuffled sequence does not lose any elements
        print("Testing correctness ver-3 of read write-->")
        data = {}
        data["x2"] = numpy.arange(-65636, -65636-128,-1, dtype = "int32")
        data["x3"] = numpy.arange(-4294967296, -4294967296-128,-1, dtype = "int64")
        data["x4"] = numpy.array(numpy.linspace(1.23e20, 128.23e20, 128), dtype = "float32")
        data["x5"] = numpy.array(numpy.linspace(1.23456789e200, 128.23456789e200, 128), dtype = "float64")
        
        write("check.ebf", "/", data, "w")
        write("check.ebf", "/rec/rec1/", data, "a")
        
        datar = read("check.ebf", "/",recon=2)
        
        for key in list(data.keys()):
            self.assertEqual(numpy.sum(datar[key] == data[key]), data[key].size)
            self.assertEqual(numpy.sum(datar['rec']['rec1'][key] == data[key]), data[key].size)
        
    def test_readpart(self):
        print("Testing partial read -->")
        x=numpy.linspace(0,99,100)
        write("check.ebf", "/x", x, "w")
        y1 = read("check.ebf", "/x",begin=0)
        y2 = read("check.ebf", "/x",begin=0,end=-10)
        y3 = read("check.ebf", "/x",begin=0,end=20)
        y4 = read("check.ebf", "/x",begin=20,end=30)
        y5 = read("check.ebf", "/x",begin=10,end=-10)
        self.assertTrue(numpy.all(y1==x))
        self.assertTrue(numpy.all(y2==x[0:-10]))
        self.assertTrue(numpy.all(y3==x[0:20]))
        self.assertTrue(numpy.all(y4==x[20:30]))
        self.assertTrue(numpy.all(y5==x[10:-10]))
        
        x=numpy.ones((10,5,2))
        write("check.ebf", "/x", x, "w")
        y1 = read("check.ebf", "/x",begin=0,end=-1)
        self.assertTrue(numpy.all(y1==x[0:-1,:,:]))
        
    def test_iterate(self):
        print("Testing iterate -->")
        x=numpy.linspace(0,99,100)
        write("check.ebf", "/x", x, "w")
        
        cache=12
        begin=0
        end=cache
        for data in iterate('check.ebf','/x',cache):
            self.assertTrue(numpy.all(data==x[begin:end]))
            begin=end
            end=end+cache
            if end>100:
                end=100
                
        x=numpy.linspace(0,99,100)
        y=numpy.linspace(100,199,100)
        z=numpy.linspace(200,199,100)
        write('check.ebf','/x',x,'w')
        write('check.ebf','/y',y,'a')
        write('check.ebf','/z',z,'a')
        cache=12
        begin=0
        end=cache
        for data in iterate('check.ebf','/x+',cache):
            self.assertTrue(numpy.all(data["x"]==x[begin:end]))
            self.assertTrue(numpy.all(data["y"]==y[begin:end]))
            self.assertTrue(numpy.all(data["z"]==z[begin:end]))
            begin=end
            end=end+cache
            if end>100:
                end=100
                      
    def test_ebffile(self):
        print("Testing ebffile -->")
        dt=[('x','float64'),('y','float64')]
        data=numpy.zeros(100,dtype=dt)
        data['x']=numpy.linspace(0,99,100)
        write("check.ebf", "/data", data, "w")
        efile=EbfFile('check.ebf','/data','r')
        data1=efile.read(0,100)
        self.assertTrue(numpy.all(data1['x']==data['x']))
        self.assertTrue(numpy.all(data1['y']==data['y']))
        
    def test_structunits(self):
        print("Testing ebffile -->")
        dt=[('x','float64',(2,)),('y','float64',(5,))]
        data=numpy.zeros(100,dtype=dt)
        write("check.ebf", "/data", data, "w",dataunit='')
#        info("check.ebf",1)
        
        dt=[('x','float64'),('y','float64')]
        data=numpy.zeros(100,dtype=dt)
        data['x']=numpy.linspace(0,99,100)
        write("check.ebf", "/data", data, "w",dataunit='')
#        info("check.ebf",1)
        units1=['u=km/s,l=\\alpha','u=m/s,l=\\gamma']
        write("check.ebf", "/data", data, "w",dataunit=units1)
#        info("check.ebf",1)
        units2=unit("check.ebf", "/data")
        print('printing ',units2)
        self.assertTrue(numpy.all(units1==units2))
        
    def test_read_ind(self):
        print("Testing read_ind -->")
        dt=[('x','float64'),('y','float64')]
        nsize=100
        data=numpy.zeros(nsize,dtype=dt)
        data['x']=numpy.linspace(0,nsize-1,nsize)
        data['y']=numpy.linspace(0,nsize-1,nsize)
        write('check.ebf','/data',data,'w')
        ind=numpy.array([5,9,2,10,5])
        data1=read_ind('check.ebf','/data',ind)
        self.assertTrue(numpy.all(data['x'][ind]==data1['x']))
        self.assertTrue(numpy.all(data['y'][ind]==data1['y']))
        data1=read_ind('check.ebf','/data',1)
        print(type(data1),hasattr(data1,'__len__'))
        self.assertTrue(numpy.all(type(data1)==numpy.void))
        
        write('check.ebf','/data',numpy.array([]),'w')
        data1=read_ind('check.ebf','/data',[0])
        self.assertTrue(numpy.all(data1==None))
        
        write('check.ebf','/data',10,'w')
        data1=read_ind('check.ebf','/data',0)
        self.assertTrue(data1==10)
        

    def test_update(self):
        dt=[('x','float64'),('y','float64')]
        nsize=100        
        data=numpy.zeros(nsize,dtype=dt)
        write('check.ebf','/data',data,'w')
        data['x']=numpy.linspace(0,nsize-1,nsize)
        data['y']=numpy.linspace(0,nsize-1,nsize)
        write('check.ebf','/data',data,'u')
        datar=read('check.ebf','/data')
        self.assertTrue(numpy.all(data==datar))
        data=numpy.concatenate([data,data])
        try:
            write('check.ebf','/data',data,'u')
            self.assertTrue(False)
        except RuntimeError:
            self.assertTrue(True)
        dt=[('x','float64'),('y','int64')]
        nsize=100        
        data=numpy.zeros(nsize,dtype=dt)
        try:
            write('check.ebf','/data',data,'u')
            self.assertTrue(False)
        except RuntimeError:
            self.assertTrue(True)
            
            
#        datar=read('check.ebf','/data')
#        self.assertTrue(numpy.all(data==datar))
                

    def test_extend(self):
        dt=[('x','float64'),('y','float64')]
        nsize=100        
        data=numpy.zeros(nsize,dtype=dt)
        write('check.ebf','/data',data,'w')
        data['x']=numpy.linspace(0,nsize-1,nsize)
        data['y']=numpy.linspace(0,nsize-1,nsize)
        write('check.ebf','/data',data[0],'e')
        write('check.ebf','/data',data,'e')
        datar=read('check.ebf','/data')
        self.assertTrue(numpy.all(datar.size==(2*nsize+1)))
        self.assertTrue(numpy.all(datar[nsize+1:2*nsize+1]==data))
        self.assertTrue(numpy.all(datar[nsize]==data[0]))
        
    def test_update_ind(self):
        print("Testing update_ind -->")
        dt=[('x','float64'),('y','float64'),('z','S10')]
        nsize=100

        data=numpy.zeros(nsize,dtype=dt)
        write('check.ebf','/data',data,'w')        
        ind=numpy.arange(nsize)
        data['x'][ind]=numpy.linspace(0,nsize-1,nsize)[ind]
        data['y'][ind]=numpy.linspace(0,nsize-1,nsize)[ind]
        update_ind('check.ebf','/data',data)
        datar=read('check.ebf','/data')
        self.assertTrue(numpy.all(data==datar))

        
        data=numpy.zeros(nsize,dtype=dt)
        write('check.ebf','/data',data,'w')        
        ind=numpy.arange(1)+20
        data['x'][ind]=numpy.linspace(0,nsize-1,nsize)[ind]
        data['y'][ind]=numpy.linspace(0,nsize-1,nsize)[ind]
        update_ind('check.ebf','/data',data[ind[0]],ind[0])
        datar=read('check.ebf','/data')
        self.assertTrue(numpy.all(data==datar))
        
        
        data=numpy.zeros(nsize,dtype=dt)
        write('check.ebf','/data',data,'w')        
        ind=numpy.arange(50)+20
        data['x'][ind]=numpy.linspace(0,nsize-1,nsize)[ind]
        data['y'][ind]=numpy.linspace(0,nsize-1,nsize)[ind]
        update_ind('check.ebf','/data',data[ind],ind)
        datar=read('check.ebf','/data')
        self.assertTrue(numpy.all(data==datar))
        
        write('check.ebf','/data',data,'w')        
        dt1=[('x','int32'),('y','int32'),('z','S10')]
        data1=numpy.zeros(nsize,dtype=dt1)
        data1['x'][ind]=numpy.arange(nsize,dtype='int32')[ind]
        data1['y'][ind]=numpy.arange(nsize,dtype='int32')[ind]
        update_ind('check.ebf','/data',data1[ind],ind)
        datar=read('check.ebf','/data')
        self.assertTrue(numpy.all(data==datar))
        
        x=numpy.zeros(nsize,dtype='S10')
        write('check.ebf','/x',x,'w')        
        x1=numpy.arange(nsize)
        update_ind('check.ebf','/x',x1,numpy.arange(nsize))
        datar=read('check.ebf','/x')
        self.assertTrue(numpy.all(numpy.float64(x1)==numpy.float64(datar)))
        
        


        x=numpy.zeros(nsize,dtype='float64')
        write('check.ebf','/x',x,'w')        
        x1=numpy.array(numpy.arange(nsize),dtype='S10')
#        x1=numpy.array(x1,dtype='float64')
        update_ind('check.ebf','/x',x1,numpy.arange(nsize))
        datar=read('check.ebf','/x')
        self.assertTrue(numpy.all(numpy.float64(x1)==numpy.float64(datar)))
        
        write('check.ebf','/x',10,'w')        
        update_ind('check.ebf','/x',20,0)
        x=read('check.ebf','/x')
        self.assertTrue(x==20)

    def test_dict2npstruct(self):
        print("Testing dict2npstruct and npstruct2dict -->")
        dt=[('x','float64'),('y','float64'),('z','S10')]
        nsize=100
        data1=numpy.zeros(nsize,dtype=dt)
        data1['x']=numpy.arange(nsize,dtype='int32')
        data1['y']=numpy.arange(nsize,dtype='int32')+10
        data1['y']=numpy.arange(nsize,dtype='int32')+20
        data2=npstruct2dict(data1)
        for key in data1.dtype.names:
            self.assertTrue(numpy.all(data1[key]==data2[key]))
        
        data2['extra']=numpy.arange(10)
        data3=dict2npstruct(data2,basekey='x')
        for key in data1.dtype.names:
            self.assertTrue(numpy.all(data1[key]==data3[key]))
        self.assertTrue(len(data3.dtype.names)==3)
        
        data3=dict2npstruct(data1,keylist=['x','y','z'])
        for key in data1.dtype.names:
            self.assertTrue(numpy.all(data1[key]==data3[key]))
        
        data3=dict2npstruct(data2,keylist=['y','z'])
        for key in data3.dtype.names:
            self.assertTrue(numpy.all(data1[key]==data3[key]))
        self.assertTrue(len(data3.dtype.names)==2)
        self.assertTrue('y' in data3.dtype.names)
        self.assertTrue('z' in data3.dtype.names)

    def test_copy(self):
        print("Testing copy -->")
        dt=[('x','float64'),('y','float64'),('z','S10')]
        nsize=100
        data1=numpy.zeros(nsize,dtype=dt)
        data1['x']=numpy.arange(nsize,dtype='int32')
        data1['y']=numpy.arange(nsize,dtype='int32')
        write('check.ebf','/data',data1,'w')
        copy('check.ebf','check1.ebf','w','/data','/')
        data2=read('check1.ebf','/')
        self.assertTrue(type(data2)==dict)
        for key in data1.dtype.names:
            self.assertTrue(numpy.all(data1[key]==data2[key]))
            
#    def test_cat(self):
#        print "Testing cat -->"
#        dt=[('x','float64'),('y','float64'),('z','S30')]
#        nsize=10
#        data1=numpy.zeros(nsize,dtype=dt)
#        data1['x']=numpy.arange(nsize,dtype='int32')
#        data1['y']=numpy.arange(nsize,dtype='int32')
#        data1['z']='absdfeffffffffffffffffffllllmmmmm'
#        data1['x'][-1]=-1.234567890123456789e+101
#        write('check.ebf','/data',data1,'w')
#        write('check.ebf','/x1',data1['x'],'a')
#        write('check.ebf','/y1',data1['y'],'a')
#        write('check.ebf','/z1',data1['z'],'a')
#        cat('check.ebf','/z1 /x1 /y1 /data',' ',1)

    
#def __test_scalar():
#    cat('//home/sharma/sw/share/ebf/scalar.ebf','/x1+',' ',1)
#    ebfdir='data/'
#    data=read(ebfdir+'scalar.ebf','/')
#    x1=read(ebfdir+'scalar.ebf','/x1')
#    x2=read(ebfdir+'scalar.ebf','/x2')
#    print data
#    print x1.shape
#    print x2.shape
#    
#    x=numpy.zeros(1,dtype=[('x','int32'),('y','int32',(1,)),('z','int32')])
#    print type(x[0])
#    write('check.ebf','/',x[0],'w')
#    y=numpy.zeros(0,dtype='int32')
#    write('check.ebf','/xy1',y,'a')
#    y=numpy.int64(1)
#    
#    z=[1,2,3]
#    mystr='This is'
#    write('check.ebf','/y1',y,'a')
#    write('check.ebf','/z1',z,'a')
#    write('check.ebf','/mystr',mystr,'a')
#    info('check.ebf')
#    x=read('check.ebf','/y1')
#    print x == 1
#    print type(z)
    
#def __check():
#    file1='/work1/sharma/Projects/Stellarhalo/data/halo02.ebf' 
#    cat(file1,'/log')
#    raise RuntimeError('')

if __name__  ==  '__main__':
#    __test_scalar()
#    _checkSpeed()
    unittest.main()
#    __check()

    if len(sys.argv) == 1:
        _usage()
    elif len(sys.argv) == 2:
        if sys.argv[1] == '-speed':
            _checkSpeed()
        elif sys.argv[1] == '-help':
            _usage()
        else:
            info(sys.argv[1])
    else:
        if sys.argv[1] == '-list':
            info(sys.argv[2],1)
        elif sys.argv[1] == '-stat':
            stat(sys.argv[2],sys.argv[3])
        elif sys.argv[1] == '-print':
            cat(sys.argv[2],sys.argv[3],' ',0)
        elif sys.argv[1] == '-cat':
            cat(sys.argv[2],sys.argv[3],' ',0)
        elif sys.argv[1] == '-ssv':
            cat(sys.argv[2],sys.argv[3],' ',1)
        elif sys.argv[1] == '-csv':
            cat(sys.argv[2],sys.argv[3],', ',1)                                
        elif sys.argv[1] == '-swap':
            swapEndian(sys.argv[2])
        elif sys.argv[1] == '-diff':
            diff(sys.argv[2],sys.argv[3])
        elif sys.argv[1] == '-htab':
            _EbfTable.display_htab(sys.argv[2])
        elif sys.argv[1] == '-copy':
            if len(sys.argv) == 4: 
                copy(sys.argv[2],sys.argv[3],'a')
            elif len(sys.argv) == 5:    
                copy(sys.argv[2],sys.argv[3],'a',sys.argv[4])
            elif len(sys.argv) == 6:    
                copy(sys.argv[2],sys.argv[3],'a',sys.argv[4],sys.argv[5])
            else:
                _usage()
        elif sys.argv[1] == '-rename':
            if len(sys.argv) == 5: 
                rename(sys.argv[2],sys.argv[3],sys.argv[4])            
            else:
                _usage()                
        elif sys.argv[1] == '-remove':
            if len(sys.argv) == 4: 
                rename(sys.argv[2],sys.argv[3],'')
            else:
                _usage()                       
        elif sys.argv[1] == '-join':
            if len(sys.argv) == 5: 
                from glob import glob
                if '*' in sys.argv[2]:
                    filelist=glob(sys.argv[2])
                    filelist.sort()
                    print(filelist)
                    join(filelist,'/',sys.argv[3],'/',sys.argv[4])
                else:
                    join(sys.argv[2],'/',sys.argv[3],'/',sys.argv[4])
            else:
                _usage()
        else:
            _usage()
        


    
    
